#!/usr/bin/env python3
"""E1: same-fit baseline/gene logit decomposition, aggregate output only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from frozen_geo_scorer import FrozenGEOScorer
from score_ppmi_after_lock import load_locked_input
from postlock_common import private_output, require_execution_authorization, summary, write_aggregate


def decompose(design: np.ndarray, intercept: float, coefficients: np.ndarray,
              n_baseline: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    design = np.asarray(design, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)
    if (design.ndim != 2 or design.shape[1] != len(coefficients)
            or not 0 <= n_baseline <= len(coefficients)
            or not np.isfinite(design).all() or not np.isfinite(coefficients).all()):
        raise ValueError("Invalid frozen M5 design or coefficients")
    baseline = design[:, :n_baseline] @ coefficients[:n_baseline]
    gene = design[:, n_baseline:] @ coefficients[n_baseline:]
    logit = intercept + baseline + gene
    if not np.allclose(logit, intercept + design @ coefficients,
                       rtol=0, atol=1e-8):
        raise AssertionError("M5 component reconstruction failed")
    return baseline, gene, logit


def component_arrays(scorer: FrozenGEOScorer, ranks: np.ndarray,
                     age: np.ndarray, sex: np.ndarray, genes: np.ndarray):
    scorer.validate_input(ranks, age, sex, genes)
    reused = scorer.decisions.get("M5", {}).get("selected", {}).get("reuses_baseline")
    if reused:
        if reused != "M4":
            raise ValueError("Unexpected M5 baseline reuse")
        x = scorer.design(ranks, age, sex, "M4")
        intercept, coef = scorer.models["M4"]
        baseline, _ignored, logit = decompose(x, intercept, coef, len(coef))
        gene = np.zeros(len(ranks), dtype=float)
    else:
        x = scorer.design(ranks, age, sex, "M5")
        intercept, coef = scorer.models["M5"]
        baseline, gene, logit = decompose(x, intercept, coef, len(scorer.keep))
    if not np.allclose(expit(logit), scorer.predict_models(ranks, age, sex, genes)["M5"],
                       rtol=0, atol=1e-8):
        raise AssertionError("Decomposition differs from frozen M5 predictions")
    return baseline, gene


def absolute_smd(a: np.ndarray, b: np.ndarray) -> float | None:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float(abs(a.mean() - b.mean()) / pooled) if pooled > 0 else None


def validate_scopes(scopes: list, frozen: dict) -> None:
    expected = {entry["scope"]: entry for entry in frozen["splits"]
                if entry["scope"] != "pooled/final"}
    if not isinstance(scopes, list) or len(scopes) != len(expected):
        raise ValueError("E1 requires every frozen outer and transfer scope exactly once")
    names = [entry["scope"] for entry in scopes]
    if len(set(names)) != len(names) or set(names) != set(expected):
        raise ValueError("E1 scope names differ from the frozen split manifest")
    for item in scopes:
        locked = expected[item["scope"]]
        if (item["train_samples"] != locked["train"]
                or item["heldout_samples"] != locked["test"]):
            raise ValueError("E1 samples/order differ from frozen GEO splits")


def verify_fit_samples(model_dir: Path, train_samples: np.ndarray) -> None:
    with np.load(model_dir / "preprocessing.npz", allow_pickle=False) as z:
        actual = str(z["fit_samples_hash"].item())
    if actual != geo.array_hash(train_samples):
        raise ValueError("E1 model preprocessing was fit on different GEO samples")


def scope_summary(scorer, geo_data, train_ix, heldout_ix,
                  ppmi_ranks, ppmi_age, ppmi_sex, ppmi_y,
                  within_cohort_outer: bool):
    train = geo_data.take(train_ix)
    heldout = geo_data.take(heldout_ix)
    if set(train.donor) & set(heldout.donor):
        raise ValueError("GEO scope has donor overlap")
    if not np.array_equal(geo_data.genes, scorer.genes):
        raise ValueError("GEO model gene order differs")
    gt_b, gt_g = component_arrays(scorer, train.X, train.age, train.sex, train.genes)
    gh_b, gh_g = component_arrays(scorer, heldout.X, heldout.age, heldout.sex, heldout.genes)
    pp_b, pp_g = component_arrays(scorer, ppmi_ranks, ppmi_age, ppmi_sex, scorer.genes)
    gt, gh = train.y == 0, heldout.y == 0
    pc, pd = ppmi_y == 0, ppmi_y == 1
    ratio = None
    if within_cohort_outer and gt.sum() > 1 and gh.sum() > 1:
        denom = gt_g[gt].std(ddof=1)
        ratio = float(gh_g[gh].std(ddof=1) / denom) if denom > 0 else None
    return {
        "GEO_train_controls": {"baseline": summary(gt_b[gt]), "gene": summary(gt_g[gt])},
        "GEO_heldout_controls": {"baseline": summary(gh_b[gh]), "gene": summary(gh_g[gh])},
        "PPMI_controls_same_fit": {"baseline": summary(pp_b[pc]), "gene": summary(pp_g[pc])},
        "PPMI_PD_same_fit": {"baseline": summary(pp_b[pd]), "gene": summary(pp_g[pd])},
        "absolute_gene_SMD_heldout_GEO_vs_PPMI_controls": absolute_smd(gh_g[gh], pp_g[pc]),
        "gene_SD_ratio_heldout_to_train": ratio,
        "material_gene_shift_rule": bool(absolute_smd(gh_g[gh], pp_g[pc]) >= .5)
            if absolute_smd(gh_g[gh], pp_g[pc]) is not None else None,
        "scale_instability_rule": bool(ratio >= 2) if ratio is not None else None,
    }


def synthetic_self_test() -> None:
    x = np.array([[1., 2., 3.], [2., -1., 4.], [0., 0., 1.]])
    coef = np.array([.5, -.25, 2.])
    b, g, eta = decompose(x, .3, coef, 2)
    assert np.allclose(eta, .3 + x @ coef)
    assert np.allclose(g, x[:, 2] * 2.)
    assert absolute_smd(np.array([0., 1., 2.]), np.array([1., 2., 3.])) > 0
    assert summary(b)["n"] == 3
    locked = {"splits": [{"scope": "A/outer/0", "train": ["a"], "test": ["b"]},
                         {"scope": "pooled/final", "train": ["a", "b"], "test": []}]}
    validate_scopes([{"scope": "A/outer/0", "train_samples": ["a"],
                      "heldout_samples": ["b"], "model_dir": "unused"}], locked)
    try:
        validate_scopes([{"scope": "A/outer/0", "train_samples": ["b"],
                          "heldout_samples": ["a"]}], locked)
    except ValueError:
        pass
    else:
        raise AssertionError("Wrong frozen split was accepted")
    print("E1_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--scope-manifest", type=Path,
                        help="JSON list: scope, model_dir, train_samples, heldout_samples")
    parser.add_argument("--locked-splits", type=Path,
                        default=ROOT / "data/manifests/splits.json")
    parser.add_argument("--final-dir", type=Path)
    parser.add_argument("--marker-manifest", type=Path)
    parser.add_argument("--sex-artifact", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    fields = (args.authorization_record, args.bundle, args.binding,
              args.scope_manifest, args.final_dir, args.marker_manifest,
              args.sex_artifact, args.output)
    if any(value is None for value in fields):
        parser.error("Real E1 execution requires all input paths and authorization")
    require_execution_authorization(args.authorization_record, analysis="E1")
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    pp_x, pp_age, pp_sex, pp_y, _donors, pp_genes = load_locked_input(args.bundle, binding)
    geo_data, _ = geo.load_geo()
    if not np.array_equal(geo_data.genes, pp_genes):
        raise ValueError("GEO/PPMI common gene order differs")
    scopes = json.loads(args.scope_manifest.read_text(encoding="utf-8"))
    frozen = json.loads(args.locked_splits.read_text(encoding="utf-8"))
    validate_scopes(scopes, frozen)
    sample_lookup = {str(s): i for i, s in enumerate(geo_data.sample)}
    result = {"analysis": "E1", "status": "complete_exploratory", "scopes": {}}
    for item in scopes:
        name = str(item["scope"])
        if name in result["scopes"]:
            raise ValueError("Duplicate E1 scope")
        ti = np.array([sample_lookup[str(s)] for s in item["train_samples"]], int)
        hi = np.array([sample_lookup[str(s)] for s in item["heldout_samples"]], int)
        verify_fit_samples(Path(item["model_dir"]), geo_data.sample[ti])
        scorer = FrozenGEOScorer(Path(item["model_dir"]), args.marker_manifest, args.sex_artifact)
        result["scopes"][name] = scope_summary(scorer, geo_data, ti, hi,
            pp_x, pp_age, pp_sex, pp_y, within_cohort_outer="/outer/" in name)
    verify_fit_samples(args.final_dir, geo_data.sample)
    final = FrozenGEOScorer(args.final_dir, args.marker_manifest, args.sex_artifact)
    pb, pg = component_arrays(final, pp_x, pp_age, pp_sex, pp_genes)
    gb, gg = component_arrays(final, geo_data.X, geo_data.age, geo_data.sex, geo_data.genes)
    result["final_GEO_model_PPMI"] = {
        "controls": {"baseline": summary(pb[pp_y == 0]), "gene": summary(pg[pp_y == 0])},
        "PD": {"baseline": summary(pb[pp_y == 1]), "gene": summary(pg[pp_y == 1])}}
    result["final_GEO_model_GEO_training_controls"] = {
        "baseline": summary(gb[geo_data.y == 0]),
        "gene": summary(gg[geo_data.y == 0])}
    write_aggregate(output, result)
    print("E1_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
