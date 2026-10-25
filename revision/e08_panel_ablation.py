#!/usr/bin/env python3
"""E8: M4 leave-one-panel-out fits; frozen GEO stage before PPMI scoring."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import warnings

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from score_ppmi_after_lock import load_locked_input
from postlock_common import (private_output, require_execution_authorization,
                             rng_for, sha256, synthetic_ranks, write_aggregate)


def paired_auc_interval(y: np.ndarray, donors: np.ndarray, full: np.ndarray,
                        ablated: np.ndarray, stream: int, draws: int = 5000) -> list[float]:
    """Resample donors by diagnosis, using identical indices for paired scores."""
    groups = {}
    for i, donor in enumerate(map(str, donors)):
        groups.setdefault(donor, []).append(i)
    by_class = {label: [] for label in (0, 1)}
    for ix in groups.values():
        labels = np.unique(y[ix])
        if len(labels) != 1:
            raise ValueError("Donor has inconsistent diagnosis in E8 bootstrap")
        by_class[int(labels[0])].append(np.asarray(ix, int))
    if not all(by_class.values()):
        raise ValueError("Both classes required for E8 bootstrap")
    rng = rng_for(8, stream)
    values = np.empty(draws, float)
    for repetition in range(draws):
        sampled = []
        for label in (0, 1):
            donor_groups = by_class[label]
            sampled.extend(donor_groups[i] for i in rng.integers(
                0, len(donor_groups), len(donor_groups)))
        ix = np.concatenate(sampled)
        values[repetition] = (roc_auc_score(y[ix], ablated[ix])
                              - roc_auc_score(y[ix], full[ix]))
    return np.quantile(values, [.025, .975]).astype(float).tolist()


class PanelOmissionPrep(geo.Prep):
    def __init__(self, train: geo.Data, omitted_panel: int | None):
        self.omitted_panel = omitted_panel
        super().__init__(train)

    def baseline_raw(self, data: geo.Data) -> np.ndarray:
        full = super().baseline_raw(data)
        if self.omitted_panel is None:
            return full
        if not 0 <= self.omitted_panel < 6:
            raise ValueError("Omitted panel index must be 0..5")
        return np.delete(full, 4 + self.omitted_panel, axis=1)


def fit_ablation(train: geo.Data, omitted_panel: int | None, scope: str,
                 out: Path):
    out.mkdir(parents=True, exist_ok=False)
    prep = PanelOmissionPrep(train, omitted_panel)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        fit, _ = geo.checked_fit(prep, train, "M4", 0., 0.,
                                 dict(scope=scope, stage="final_refit", panel=omitted_panel),
                                 out / "fits.jsonl")
    prep.save(out / "preprocessing.npz")
    np.savez_compressed(out / "M4_model.npz", intercept=fit.intercept_, coef=fit.coef_)
    return prep, fit, {"fit_warning_count": len(captured),
                       "baseline_features_retained": [int(i) for i in prep.keep],
                       "standardized_coefficients": np.asarray(fit.coef_, float).tolist()}


def score(prep: PanelOmissionPrep, fit, data: geo.Data) -> np.ndarray:
    x, _ = prep.transform(data, "M4")
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        pred = np.asarray(fit.predict(x), float)
    if pred.shape != data.y.shape or not np.isfinite(pred).all():
        raise ValueError("Nonfinite M4 ablation predictions")
    return pred


def loaded_score(out: Path, data: geo.Data, panel: int | None) -> np.ndarray:
    with np.load(out / "preprocessing.npz", allow_pickle=False) as z:
        prep = PanelOmissionPrep.__new__(PanelOmissionPrep)
        prep.omitted_panel = panel
        for key in ("age_mean", "age_sd", "base_mean", "base_sd", "keep",
                    "selected", "gene_mean", "gene_sd", "genes", "fit_samples_hash"):
            value = z[key]
            setattr(prep, key, value.item() if value.ndim == 0 else value.copy())
    if not np.array_equal(prep.genes, data.genes):
        raise ValueError("Ablation gene universe differs from scoring data")
    with np.load(out / "M4_model.npz", allow_pickle=False) as z:
        intercept, coef = float(z["intercept"]), np.asarray(z["coef"], float)
    x, _ = prep.transform(data, "M4")
    if x.shape[1] != len(coef):
        raise ValueError("Ablation design and model dimensions differ")
    from scipy.special import expit
    p = expit(intercept + x @ coef)
    if not np.isfinite(p).all():
        raise ValueError("Nonfinite frozen ablation score")
    return p


def settings() -> list[tuple[str, int | None]]:
    return [("full", None)] + [(f"omit_{i}_{geo.PANELS[i].replace(' ', '_')}", i)
                               for i in range(6)]


def geo_stage(output: Path) -> dict:
    data, source_hashes = geo.load_geo()
    result = {"analysis": "E8", "stage": "GEO_complete_before_PPMI",
              "script_sha256": sha256(Path(__file__)),
              "geo_input_sha256": source_hashes,
              "transfer": {}, "pooled": {}}
    sources = sorted(set(map(str, data.source)))
    if len(sources) != 2:
        raise ValueError("Exactly two GEO sources expected")
    reference_predictions = {}
    for label, panel in settings():
        pooled_prep, pooled_fit, pooled_info = fit_ablation(
            data, panel, f"E8/pooled/{label}", output / "pooled" / label)
        result["pooled"][label] = pooled_info
        for source_ix, source in enumerate(sources):
            train = data.take(np.flatnonzero(data.source == source))
            test = data.take(np.flatnonzero(data.source != source))
            prep, fit, info = fit_ablation(train, panel,
                f"E8/{source}/transfer/{label}", output / source / label)
            key = f"{source}_to_{str(test.source[0])}"
            prediction = score(prep, fit, test)
            row = {"AUROC": float(roc_auc_score(test.y, prediction)), **info}
            if label == "full":
                reference_predictions[key] = prediction
                row["delta_from_full_M4_CI95_percentile"] = [0., 0.]
            else:
                row["delta_from_full_M4_CI95_percentile"] = paired_auc_interval(
                    test.y, test.donor, reference_predictions[key], prediction,
                    stream=100 + source_ix * 10 + int(panel))
            result["transfer"].setdefault(key, {})[label] = row
        print(f"E8 GEO {label} complete ({len(result['pooled'])}/7 settings)", flush=True)
    for direction, rows in result["transfer"].items():
        reference = rows["full"]["AUROC"]
        for row in rows.values():
            row["delta_from_full_M4"] = row["AUROC"] - reference
    write_aggregate(output / "E8_GEO_AGGREGATE.json", result)
    return result


def ppmi_stage(output: Path, bundle: Path, binding: Path) -> dict:
    geo_report = json.loads((output / "E8_GEO_AGGREGATE.json").read_text(encoding="utf-8"))
    if (geo_report.get("stage") != "GEO_complete_before_PPMI"
            or geo_report.get("script_sha256") != sha256(Path(__file__))):
        raise ValueError("E8 GEO stage or its code hash changed")
    binding_obj = json.loads(binding.read_text(encoding="utf-8"))
    ranks, age, sex, y, donors, genes = load_locked_input(bundle, binding_obj)
    manifest = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    data = geo.Data(ranks, age, sex, y, np.repeat("PPMI", len(y)), donors,
                    donors.copy(), genes, markers)
    result = {"analysis": "E8", "stage": "PPMI_scored_after_frozen_GEO",
              "exploratory": True, "pooled_fit_only": True, "models": {}}
    reference_prediction = None
    for label, panel in settings():
        p = loaded_score(output / "pooled" / label, data, panel)
        row = {"AUROC": float(roc_auc_score(y, p))}
        if label == "full":
            reference_prediction = p
            row["delta_from_full_M4_CI95_percentile"] = [0., 0.]
        else:
            row["delta_from_full_M4_CI95_percentile"] = paired_auc_interval(
                y, donors, reference_prediction, p, stream=200 + int(panel))
        result["models"][label] = row
    reference = result["models"]["full"]["AUROC"]
    for row in result["models"].values():
        row["delta_from_full_M4"] = row["AUROC"] - reference
    write_aggregate(output / "E8_PPMI_AGGREGATE.json", result)
    print("E8_COMPLETE: PPMI aggregate output written under ~/PPMI_private")
    return result


def synthetic_self_test() -> None:
    rng = np.random.default_rng(33)
    n, g = 48, 40
    ranks = synthetic_ranks(n, g, 33)
    y = np.array([0, 1] * (n // 2))
    data = geo.Data(ranks, rng.uniform(50, 75, n), rng.integers(0, 2, n), y,
                    np.repeat("SYNTHETIC", n), np.array([f"s{i}" for i in range(n)]),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"G{i:04d}" for i in range(g)]),
                    [np.arange(i*3, i*3+3) for i in range(6)])
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        for panel in (None, 2):
            folder = Path(temporary) / ("full" if panel is None else "omit")
            prep, fit, _ = fit_ablation(data, panel, "E8/synthetic", folder)
            direct = score(prep, fit, data)
            restored = loaded_score(folder, data, panel)
            assert np.allclose(direct, restored, atol=1e-10, rtol=0)
            assert len(prep.base_mean) == (10 if panel is None else 9)
        same_ci = paired_auc_interval(y, data.donor, direct, direct, 999, draws=20)
        assert same_ci == [0., 0.]
    print("E8_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--stage", choices=("geo", "ppmi"))
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if args.stage is None or args.authorization_record is None or args.output_dir is None:
        parser.error("Real E8 requires stage, authorization and output directory")
    require_execution_authorization(args.authorization_record, analysis="E8")
    os.umask(0o077)
    output = private_output(args.output_dir)
    if args.stage == "geo":
        if output.exists():
            raise ValueError("Use a fresh E8 GEO stage directory")
        output.mkdir(mode=0o700, parents=True)
        geo_stage(output)
    else:
        if args.bundle is None or args.binding is None:
            parser.error("PPMI stage requires bundle and binding")
        ppmi_stage(output, args.bundle, args.binding)


if __name__ == "__main__":
    main()
