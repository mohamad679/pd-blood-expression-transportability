"""Recheck one frozen GEO scope's ComBat matrices and retain numerical stderr.

This diagnostic repeats ComBat fits and transforms only. It does not refit M1–M5
or use PPMI. The registered full sensitivity outputs are not modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
import warnings
from contextlib import redirect_stderr
from pathlib import Path

import numpy as np

import geo_pipeline as gp
import run_geo_combat_sensitivity as combat


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/geo_combat_numeric_audit_20261020"


def finite_matrix(label: str, original: np.ndarray, adjusted: np.ndarray,
                  corrected_mask: np.ndarray) -> dict:
    if original.shape != adjusted.shape or corrected_mask.shape != (len(original),):
        raise ValueError(f"Shape mismatch at {label}")
    nan_count = int(np.isnan(adjusted).sum())
    inf_count = int(np.isinf(adjusted).sum())
    if nan_count or inf_count:
        raise ValueError(f"Nonfinite adjusted values at {label}: NaN={nan_count}, Inf={inf_count}")
    unchanged_small = bool(np.array_equal(original[~corrected_mask],
                                          adjusted[~corrected_mask]))
    if not unchanged_small:
        raise ValueError(f"Small/unseen-batch samples changed at {label}")
    return {
        "partition": label,
        "sample_count": int(len(original)),
        "gene_count": int(original.shape[1]),
        "corrected_count": int(corrected_mask.sum()),
        "uncorrected_count": int((~corrected_mask).sum()),
        "uncorrected_samples_exactly_unchanged": unchanged_small,
        "adjusted_nan_count": nan_count,
        "adjusted_inf_count": inf_count,
        "adjusted_min": float(np.min(adjusted)),
        "adjusted_max": float(np.max(adjusted)),
    }


def check_estimates(label: str, estimates: dict) -> None:
    for key, value in estimates.items():
        values = np.asarray(value)
        if values.dtype.kind in "iuf" and not np.isfinite(values).all():
            raise ValueError(f"Nonfinite ComBat estimate at {label}: {key}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", default="GSE99039/outer/0")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to reuse nonempty diagnostic directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    data, input_hashes = gp.load_geo()
    splits = gp.split_manifest(data, input_hashes)
    payload = json.dumps(splits, indent=2, ensure_ascii=False, allow_nan=False).encode()
    split_hash = hashlib.sha256(payload).hexdigest()
    if split_hash != combat.EXPECTED_SPLIT_SHA256:
        raise ValueError(f"Frozen split hash mismatch: {split_hash}")
    split = next((item for item in splits["splits"] if item["scope"] == args.scope), None)
    if split is None or args.scope == "pooled/final":
        raise ValueError(f"Expected a frozen GEO evaluation scope: {args.scope}")
    batch_lookup, batch_hashes = combat.load_batch_lookup()
    sample_lookup = {str(sample): i for i, sample in enumerate(data.sample)}
    train, test, inner = combat.scope_data(split, data, sample_lookup)
    full_ids = list(split["train"]) + list(split["test"])
    full_data = data.take(np.asarray([sample_lookup[s] for s in full_ids], dtype=int))

    records = []
    total = 2 * (len(inner) + 1) + 1
    started = time.monotonic()
    stderr_path = output / "combat_stderr.log"
    with stderr_path.open("w", encoding="utf-8") as stderr_stream, redirect_stderr(stderr_stream):
        with warnings.catch_warnings():
            warnings.simplefilter("always", RuntimeWarning)
            for arm in ("inductive_no_disease", "inductive_preserve_disease"):
                preserve = arm == "inductive_preserve_disease"
                precheck = combat.precheck_train_only(train, inner, batch_lookup, preserve)
                blocked = [row for row in precheck if row["not_evaluable_reason"]]
                if blocked:
                    raise ValueError(f"ComBat arm blocked at {args.scope}: {blocked}")
                arm_dir = output / arm
                arm_dir.mkdir()
                partitions = [(f"inner/{i}", train.take(train_ix),
                               train.take(validation_ix))
                              for i, (train_ix, validation_ix) in enumerate(inner)]
                partitions.append(("outer_refit", train, test))
                for name, fit_data, heldout in partitions:
                    label = f"{arm}/{name}"
                    fit_batch = combat.get_batches(fit_data, batch_lookup)
                    disease = fit_data.y if preserve else None
                    state = combat.fit_partial_combat(
                        fit_data.X, fit_batch, disease, arm_dir / "neurocombat.log"
                    )
                    check_estimates(label, state["estimates"])
                    adjusted_fit, fit_mask = combat.apply_partial_combat(
                        fit_data.X, fit_batch, state, arm_dir / "neurocombat.log"
                    )
                    records.append({"arm": arm, "stage": name, "role": "train",
                                    **finite_matrix(label + "/train", fit_data.X,
                                                    adjusted_fit, fit_mask)})
                    heldout_batch = combat.get_batches(heldout, batch_lookup)
                    adjusted_heldout, heldout_mask = combat.apply_partial_combat(
                        heldout.X, heldout_batch, state, arm_dir / "neurocombat.log"
                    )
                    records.append({"arm": arm, "stage": name, "role": "heldout",
                                    **finite_matrix(label + "/heldout", heldout.X,
                                                    adjusted_heldout, heldout_mask)})
                    print(f"ComBat numeric audit {len(records)//2}/{total}: {label}", flush=True)

            arm = "transductive_no_disease"
            arm_dir = output / arm
            arm_dir.mkdir()
            full_batch = combat.get_batches(full_data, batch_lookup)
            state = combat.fit_partial_combat(
                full_data.X, full_batch, disease=None, log_path=arm_dir / "neurocombat.log"
            )
            check_estimates(arm, state["estimates"])
            adjusted, mask = combat.apply_partial_combat(
                full_data.X, full_batch, state, arm_dir / "neurocombat.log"
            )
            records.append({"arm": arm, "stage": "full_scope", "role": "all_expression",
                            **finite_matrix(arm + "/full_scope", full_data.X, adjusted, mask)})
            print(f"ComBat numeric audit {total}/{total}: {arm}", flush=True)

    stderr_text = stderr_path.read_text()
    report = {
        "status": "matrix_finiteness_passed_with_runtime_warnings" if "RuntimeWarning:" in stderr_text else "passed",
        "purpose": "repeat ComBat fits/transforms for numerical diagnosis only; no model refit",
        "scope": args.scope,
        "split_manifest_sha256": split_hash,
        "input_sha256": input_hashes,
        "batch_metadata_sha256": batch_hashes,
        "neuroCombat_version": combat.PINNED_NEUROCOMBAT,
        "minimum_training_batch_n": combat.MIN_TRAIN_BATCH_N,
        "partition_fit_count": total,
        "checked_matrix_count": len(records),
        "all_adjusted_matrices_finite": True,
        "stderr_file": stderr_path.name,
        "runtime_warning_count": stderr_text.count("RuntimeWarning:"),
        "divide_by_zero_warning_count": stderr_text.count("RuntimeWarning: divide by zero"),
        "records": records,
        "elapsed_minutes": (time.monotonic() - started) / 60,
        "PPMI_accessed": False,
    }
    combat.write_json(output / "SUMMARY.json", report)
    print(f"GEO_COMBAT_NUMERIC_AUDIT_MATRICES_FINITE: {output}", flush=True)


if __name__ == "__main__":
    main()
