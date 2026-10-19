"""One-scope leakage/runtime pilot for the registered GEO matching sensitivity."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo
from pipeline.audit_geo_age_sex_matching import match_ids
from pipeline.prediction_warning_capture import instrument_predict

SCOPE = "GSE99039/outer/1"
OUT = ROOT / "private_data" / "geo_age_sex_matching_pilot_20261019_GSE99039_outer_1"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
AUDIT = ROOT / "data/manifests/input_audit.json"
DECISION = ROOT / "docs/GEO_AGE_SEX_MATCHING_DECISION_20261019_EN.md"
MATCH_COUNTS = ROOT / "results/geo_age_sex_matching_audit_20261019.json"
WARNING_LOG = OUT / "prediction_warnings.jsonl"
MIN_CLASS_N = 5

_ORIGINAL_CHECKED_FIT = geo.checked_fit
_ORIGINAL_PREP = geo.Prep
_PREP_FIT_SAMPLE_HASHES = []


class AuditedPrep(_ORIGINAL_PREP):
    """Record exact sample membership used to fit every preprocessing object."""

    def __init__(self, train, max_genes=1000):
        super().__init__(train, max_genes)
        _PREP_FIT_SAMPLE_HASHES.append(geo.array_hash(train.sample))


def instrumented_checked_fit(*args, **kwargs):
    fit, records = _ORIGINAL_CHECKED_FIT(*args, **kwargs)
    context = dict(args[5] if len(args) > 5 else kwargs["context"])
    context.setdefault("model", args[2] if len(args) > 2 else kwargs["model"])
    context.setdefault("lambda", args[3] if len(args) > 3 else kwargs["lam"])
    context.setdefault("rho", args[4] if len(args) > 4 else kwargs["rho"])
    return instrument_predict(fit, context, WARNING_LOG), records


def _indices_for_pairs(data, pair_rows):
    paired_ids = {sample for pair in pair_rows
                  for sample in (pair["case_sample"], pair["control_sample"])}
    return np.asarray([i for i, sample in enumerate(data.sample)
                       if str(sample) in paired_ids], dtype=int)


def _split_hash(splits):
    return hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite prior matching pilot: {OUT}")
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    count_audit = json.loads(MATCH_COUNTS.read_text(encoding="utf-8"))
    if count_audit["status"] != "metadata_count_audit_complete":
        raise ValueError("Metadata matching audit is not complete")
    data, input_hashes = geo.load_geo()
    if input_hashes != audit["input_sha256"]:
        raise ValueError("Prepared GEO inputs differ from frozen input audit")
    splits = geo.split_manifest(data, input_hashes)
    split_hash = _split_hash(splits)
    if split_hash != release["split_manifest_sha256"]:
        raise ValueError("Reconstructed splits differ from the frozen GEO release")
    found = [entry for entry in splits["splits"] if entry["scope"] == SCOPE]
    if len(found) != 1:
        raise ValueError("Frozen pilot scope is not unique")
    split = found[0]
    count_row = next(row for row in count_audit["scopes"] if row["scope"] == SCOPE)
    if not count_row["tuning_scope_evaluable"]:
        raise ValueError("Selected pilot scope fails the preregistered size threshold")

    train_match = match_ids(data, split["train"])
    train_ix = _indices_for_pairs(data, train_match["pairs"])
    train = data.take(train_ix)
    train_sample_set = set(map(str, train.sample))
    test_match = match_ids(data, split["test"])
    if not test_match["passes_minimum_class_n"]:
        raise ValueError("Outer test fails the preregistered size threshold")

    local = {str(sample): i for i, sample in enumerate(train.sample)}
    inner, expected_fit_hashes, inner_counts = [], [], []
    for fold, frozen in enumerate(split["inner"]):
        restricted_train = [sample for sample in frozen["train"]
                            if sample in train_sample_set]
        restricted_validation = [sample for sample in frozen["validation"]
                                 if sample in train_sample_set]
        a_pairs = match_ids(data, restricted_train)
        b_pairs = match_ids(data, restricted_validation)
        if not b_pairs["passes_minimum_class_n"]:
            raise ValueError(f"Inner validation {fold} fails the preregistered threshold")
        ai = _indices_for_pairs(train, a_pairs["pairs"])
        bi = _indices_for_pairs(train, b_pairs["pairs"])
        a, b = train.take(ai), train.take(bi)
        if set(a.donor.tolist()) & set(b.donor.tolist()):
            raise AssertionError(f"Donor leakage in matched inner fold {fold}")
        a_ids, b_ids = set(map(str, a.sample)), set(map(str, b.sample))
        if a_ids & b_ids or not (a_ids | b_ids) <= train_sample_set:
            raise AssertionError(f"Matched membership failure in inner fold {fold}")
        train_local = np.asarray([local[str(sample)] for sample in a.sample], dtype=int)
        valid_local = np.asarray([local[str(sample)] for sample in b.sample], dtype=int)
        inner.append((train_local, valid_local))
        expected_fit_hashes.append(geo.array_hash(a.sample))
        inner_counts.append({
            "fold": fold, "restricted_train_before_rematching": len(restricted_train),
            "train_pairs": a_pairs["matched_pairs"],
            "train_samples_after_rematching": a_pairs["matched_n"],
            "validation_restricted_before_rematching": len(restricted_validation),
            "validation_pairs": b_pairs["matched_pairs"],
            "validation_samples_after_rematching": b_pairs["matched_n"],
            "validation_age_gap_years": b_pairs["age_gap_years"],
            "validation_minimum_pass": b_pairs["passes_minimum_class_n"],
        })
    expected_fit_hashes.append(geo.array_hash(train.sample))

    OUT.mkdir(parents=True)
    started, cpu_started = time.monotonic(), time.process_time()
    geo.checked_fit = instrumented_checked_fit
    geo.Prep = AuditedPrep
    try:
        prep, fits, decisions = geo.train_suite(train, inner, OUT, SCOPE)
    finally:
        geo.checked_fit = _ORIGINAL_CHECKED_FIT
        geo.Prep = _ORIGINAL_PREP
    prep_hashes_match = _PREP_FIT_SAMPLE_HASHES == expected_fit_hashes
    if not prep_hashes_match:
        raise AssertionError("Preprocessing fit partitions differ from matched inner-train/S")
    warning_count = (sum(1 for line in WARNING_LOG.read_text(encoding="utf-8").splitlines()
                         if line.strip()) if WARNING_LOG.exists() else 0)
    selected = {
        name: decision.get("selected", decision.get("hyperparameters"))
        for name, decision in decisions.items()
    }
    result = {
        "status": "matching_pilot_complete",
        "sensitivity": "within_source_recorded_sex_age_matched_GEO",
        "scope": SCOPE, "seconds": time.monotonic() - started,
        "process_CPU_seconds": time.process_time() - cpu_started,
        "outer_train_matched_pairs": train_match["matched_pairs"],
        "outer_train_samples": len(train.y),
        "outer_test_matched_pairs": test_match["matched_pairs"],
        "outer_test_scored": False,
        "inner_fold_count": len(inner), "inner_partitions": inner_counts,
        "preprocessing_fit_sample_hashes": _PREP_FIT_SAMPLE_HASHES,
        "expected_inner_train_plus_final_fit_hashes": expected_fit_hashes,
        "preprocessing_fit_membership_matches": prep_hashes_match,
        "selected_hyperparameters": selected,
        "prediction_warning_count": warning_count,
        "prediction_warning_log": str(WARNING_LOG),
        "published_pipeline_sha256": geo.digest(Path(geo.__file__)),
        "published_split_manifest_sha256": split_hash,
        "input_sha256": input_hashes,
        "decision_sha256": geo.digest(DECISION),
        "PPMI_loaded": False, "PPMI_scores_computed": False,
        "outer_holdout_performance_metrics_computed": False,
        "interpretation": "One-scope correctness/runtime pilot only; inner CV scores tune models and are not a sensitivity effect estimate.",
        "python": sys.version, "platform": platform.platform(),
    }
    geo.write_json(OUT / "RESULTS.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
