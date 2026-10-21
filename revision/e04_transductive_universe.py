#!/usr/bin/env python3
"""E4: unlabeled cross-platform rank-universe selection, staged GEO then PPMI."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
from scipy.special import expit
from sklearn.metrics import average_precision_score, log_loss, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from score_ppmi_after_lock import load_locked_input
from e03_ppmi_nested_cv import fit_pair, score_model
from postlock_common import (private_output, require_execution_authorization,
                             sha256, synthetic_ranks, write_aggregate)


def least_shifted(geo_ranks: np.ndarray, ppmi_ranks: np.ndarray,
                  genes: np.ndarray) -> np.ndarray:
    if (geo_ranks.ndim != 2 or ppmi_ranks.ndim != 2
            or geo_ranks.shape[1] != ppmi_ranks.shape[1]
            or geo_ranks.shape[1] != len(genes)):
        raise ValueError("Aligned common-universe matrices required")
    shift = np.abs(geo_ranks.mean(axis=0) - ppmi_ranks.mean(axis=0))
    order = np.lexsort((np.asarray(genes).astype(str), shift))
    return np.sort(order[:math.ceil(len(genes) * .5)])


class EligiblePrep(geo.Prep):
    def __init__(self, train: geo.Data, eligible: np.ndarray, max_genes: int = 1000):
        self.eligible = np.asarray(eligible, int)
        super().__init__(train, max_genes=max_genes)
        variability = []
        for source in np.unique(train.source):
            x = train.X[train.source == source]
            variability.append(np.median(np.abs(x - np.median(x, axis=0)), axis=0))
        mad = np.min(variability, axis=0)
        eligible = self.eligible[mad[self.eligible] > 0]
        order = np.lexsort((train.genes[eligible], -mad[eligible]))
        self.selected = eligible[order[:max_genes]]
        if len(self.selected) == 0:
            raise ValueError("No eligible positive-MAD gene")
        selected = train.X[:, self.selected]
        self.gene_mean, self.gene_sd = selected.mean(axis=0), selected.std(axis=0)
        if np.any(self.gene_sd <= 0):
            raise ValueError("Selected E4 gene has zero training SD")


def save_fit(out: Path, prep, m4, m5, decision: dict) -> None:
    prep.save(out / "preprocessing.npz")
    for name, fit in (("M4", m4), ("M5", m4 if m5 is None else m5)):
        np.savez_compressed(out / f"{name}_model.npz",
                            intercept=fit.intercept_, coef=fit.coef_)
    (out / "decision.json").write_text(json.dumps(decision, sort_keys=True), encoding="utf-8")


def restore_prediction(out: Path, data: geo.Data, name: str) -> np.ndarray:
    decision = json.loads((out / "decision.json").read_text(encoding="utf-8"))
    actual = "M4" if name == "M5" and decision["M5_reuses_M4"] else name
    with np.load(out / "preprocessing.npz", allow_pickle=False) as z:
        prep = EligiblePrep.__new__(EligiblePrep)
        for key in ("age_mean", "age_sd", "base_mean", "base_sd", "keep",
                    "selected", "gene_mean", "gene_sd", "genes", "fit_samples_hash"):
            value = z[key]
            setattr(prep, key, value.item() if value.ndim == 0 else value.copy())
    if not np.array_equal(prep.genes, data.genes):
        raise ValueError("Frozen E4 gene universe differs from input")
    x, _ = prep.transform(data, actual)
    with np.load(out / f"{actual}_model.npz", allow_pickle=False) as z:
        intercept, coef = float(z["intercept"]), np.asarray(z["coef"], float)
    if len(coef) != x.shape[1]:
        raise ValueError("Frozen E4 design/model mismatch")
    p = expit(intercept + x @ coef)
    if not np.isfinite(p).all():
        raise ValueError("Nonfinite E4 prediction")
    return p


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    return {"AUROC": float(roc_auc_score(y, p)),
            "AP": float(average_precision_score(y, p)),
            "log_loss": float(log_loss(y, p, labels=[0, 1]))}


def ppmi_data(bundle: Path, binding: Path) -> geo.Data:
    obj = json.loads(binding.read_text(encoding="utf-8"))
    ranks, age, sex, y, donors, genes = load_locked_input(bundle, obj)
    manifest = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    return geo.Data(ranks, age, sex, y, np.repeat("PPMI", len(y)), donors,
                    donors.copy(), genes, markers)


def geo_stage(output: Path, bundle: Path, binding: Path) -> None:
    geo_data, hashes = geo.load_geo()
    target = ppmi_data(bundle, binding)
    if not np.array_equal(geo_data.genes, target.genes):
        raise ValueError("GEO/PPMI gene universe differs")
    eligible = least_shifted(geo_data.X, target.X, geo_data.genes)
    result = {"analysis": "E4", "stage": "GEO_complete_before_PPMI_scoring",
              "label_free_universe_selection": True,
              "transductive_nonprospective": True,
              "script_sha256": sha256(Path(__file__)),
              "fit_helper_sha256": sha256(Path(__file__).with_name("e03_ppmi_nested_cv.py")),
              "input_bundle_sha256": sha256(bundle), "geo_input_sha256": hashes,
              "common_genes": len(geo_data.genes), "retained_genes": len(eligible),
              "eligible_gene_ids_sha256": geo.array_hash(geo_data.genes[eligible]),
              "transfer": {}}
    np.save(output / "eligible_indices.npy", eligible, allow_pickle=False)
    for source in sorted(set(map(str, geo_data.source))):
        train = geo_data.take(np.flatnonzero(geo_data.source == source))
        test = geo_data.take(np.flatnonzero(geo_data.source != source))
        inner = geo.donor_folds(train, 5, f"E4/{source}/transfer/inner")
        out = output / "transfer" / source
        factory = lambda data, max_genes: EligiblePrep(data, eligible, max_genes)
        prep, m4, m5, decision = fit_pair(train, inner, out,
            f"E4/{source}/transfer", prep_factory=factory)
        save_fit(out, prep, m4, m5, decision)
        p4 = score_model(prep, m4, test, "M4")
        p5 = p4 if m5 is None else score_model(prep, m5, test, "M5")
        result["transfer"][f"{source}_to_{str(test.source[0])}"] = {
            "M4": metrics(test.y, p4), "M5": metrics(test.y, p5),
            "M5_minus_M4_AUROC": float(roc_auc_score(test.y, p5)-roc_auc_score(test.y, p4)),
            "selection": decision}
        print(f"E4 GEO transfer {source} complete", flush=True)
    inner = geo.donor_folds(geo_data, 5, "E4/pooled/final/inner")
    out = output / "pooled"
    factory = lambda data, max_genes: EligiblePrep(data, eligible, max_genes)
    prep, m4, m5, decision = fit_pair(geo_data, inner, out,
        "E4/pooled/final", prep_factory=factory)
    save_fit(out, prep, m4, m5, decision)
    result["pooled_selection"] = decision
    write_aggregate(output / "E4_GEO_AGGREGATE.json", result)
    print("E4 GEO stage frozen; PPMI scoring is a separate command", flush=True)


def ppmi_stage(output: Path, bundle: Path, binding: Path) -> None:
    report = json.loads((output / "E4_GEO_AGGREGATE.json").read_text(encoding="utf-8"))
    if (report.get("stage") != "GEO_complete_before_PPMI_scoring"
            or report.get("script_sha256") != sha256(Path(__file__))
            or report.get("fit_helper_sha256") != sha256(Path(__file__).with_name("e03_ppmi_nested_cv.py"))
            or report.get("input_bundle_sha256") != sha256(bundle)):
        raise ValueError("E4 frozen GEO stage, code, or bundle changed")
    target = ppmi_data(bundle, binding)
    p4 = restore_prediction(output / "pooled", target, "M4")
    p5 = restore_prediction(output / "pooled", target, "M5")
    result = {"analysis": "E4", "stage": "PPMI_after_frozen_GEO",
              "transductive_nonprospective": True, "exploratory": True,
              "M4": metrics(target.y, p4), "M5": metrics(target.y, p5),
              "M5_minus_M4_AUROC": float(roc_auc_score(target.y, p5)-roc_auc_score(target.y, p4))}
    write_aggregate(output / "E4_PPMI_AGGREGATE.json", result)
    print("E4_COMPLETE: private aggregate PPMI output written")


def synthetic_self_test() -> None:
    rng = np.random.default_rng(41)
    n, g = 48, 40
    ranks = synthetic_ranks(n, g, 41)
    target = synthetic_ranks(16, g, 42)
    genes = np.array([f"G{i:04d}" for i in range(g)])
    eligible = least_shifted(ranks, target, genes)
    assert len(eligible) == 20
    data = geo.Data(ranks, rng.uniform(50, 75, n), rng.integers(0, 2, n),
                    np.array([0, 1] * (n//2)), np.repeat("SYNTHETIC", n),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"s{i}" for i in range(n)]), genes,
                    [np.arange(i*3, i*3+3) for i in range(6)])
    from e03_ppmi_nested_cv import donor_folds
    inner = list(donor_folds(data, 2, 700))
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        out = Path(temporary) / "fit"
        prep, m4, m5, decision = fit_pair(data, inner, out, "E4/synthetic",
            max_genes=8, grid=(.1,), rhos=(0.,), expansion_rounds=0,
            prep_factory=lambda data, max_genes: EligiblePrep(data, eligible, max_genes))
        save_fit(out, prep, m4, m5, decision)
        direct = score_model(prep, m4, data, "M4")
        restored = restore_prediction(out, data, "M4")
        assert np.allclose(direct, restored, atol=1e-10, rtol=0)
        assert set(prep.selected).issubset(set(eligible))
    print("E4_SYNTHETIC_TEST_PASSED")


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
    if any(value is None for value in (args.stage, args.authorization_record,
                                       args.output_dir, args.bundle, args.binding)):
        parser.error("Real E4 requires stage, authorization, bundle, binding and output")
    require_execution_authorization(args.authorization_record, analysis="E4")
    os.umask(0o077)
    output = private_output(args.output_dir)
    if args.stage == "geo":
        if output.exists():
            raise ValueError("Use a fresh E4 GEO output directory")
        output.mkdir(mode=0o700, parents=True)
        geo_stage(output, args.bundle, args.binding)
    else:
        ppmi_stage(output, args.bundle, args.binding)


if __name__ == "__main__":
    main()
