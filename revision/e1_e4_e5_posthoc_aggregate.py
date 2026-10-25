#!/usr/bin/env python3
"""Aggregate-only post hoc E1 component and E4/E5 paired intervals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from frozen_geo_scorer import FrozenGEOScorer
from e01_logit_decomposition import component_arrays
from e04_transductive_universe import ppmi_data, restore_prediction
from postlock_common import private_output, write_aggregate

SEED = [20261014, 41, 0]
DRAW_COUNT = 5000


def donor_bootstrap_indices(y: np.ndarray, donors: np.ndarray):
    groups: dict[str, list[int]] = {}
    for i, donor in enumerate(map(str, donors)):
        groups.setdefault(donor, []).append(i)
    by_class: dict[int, list[np.ndarray]] = {0: [], 1: []}
    for indices in groups.values():
        labels = np.unique(y[indices])
        if len(labels) != 1 or int(labels[0]) not in by_class:
            raise ValueError("Inconsistent donor diagnosis")
        by_class[int(labels[0])].append(np.asarray(indices, dtype=int))
    if not by_class[0] or not by_class[1]:
        raise ValueError("Both diagnosis classes are required")
    rng = np.random.default_rng(np.random.SeedSequence(SEED))
    for _ in range(DRAW_COUNT):
        sampled = []
        for label in (0, 1):
            donor_groups = by_class[label]
            sampled.extend(donor_groups[j] for j in rng.integers(
                0, len(donor_groups), size=len(donor_groups)))
        yield np.concatenate(sampled)


def absolute_smd_from_moments(a: dict, b: dict) -> float:
    denom = np.sqrt((float(a["sample_sd"]) ** 2 + float(b["sample_sd"]) ** 2) / 2)
    if denom <= 0:
        raise ValueError("Pooled SD is zero")
    return float(abs(float(a["mean"]) - float(b["mean"]))) / float(denom)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--e1-aggregate", required=True, type=Path)
    parser.add_argument("--geo-final", required=True, type=Path)
    parser.add_argument("--sex-artifact", required=True, type=Path)
    parser.add_argument("--e4-dir", required=True, type=Path)
    parser.add_argument("--e5-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = private_output(args.output)
    if output.exists():
        raise ValueError("Refusing to overwrite a prior aggregate")

    target = ppmi_data(args.bundle, args.binding)
    e1 = json.loads(args.e1_aggregate.read_text(encoding="utf-8"))
    if e1.get("analysis") != "E1" or e1.get("status") != "complete_exploratory":
        raise ValueError("Expected completed E1 aggregate")
    pooled_smd = absolute_smd_from_moments(
        e1["final_GEO_model_GEO_training_controls"]["gene"],
        e1["final_GEO_model_PPMI"]["controls"]["gene"])
    scorer = FrozenGEOScorer(args.geo_final,
        ROOT / "protocol/v3_marker_manifest.json", args.sex_artifact)
    _baseline, gene = component_arrays(scorer, target.X, target.age,
                                       target.sex, target.genes)
    e1_point = float(roc_auc_score(target.y, gene))

    configurations = {"E4": args.e4_dir / "pooled"}
    for name in ("deviance_only", "panel_exclusion_only", "combined"):
        configurations[f"E5/{name}"] = args.e5_dir / "variants" / name / "pooled"
    prior_e4 = json.loads((args.e4_dir / "E4_PPMI_AGGREGATE.json").read_text())
    prior_e5 = json.loads((args.e5_dir / "E5_PPMI_AGGREGATE.json").read_text())
    predictions = {}
    results = {}
    for name, folder in configurations.items():
        p4 = restore_prediction(folder, target, "M4")
        p5 = restore_prediction(folder, target, "M5")
        a4, a5 = float(roc_auc_score(target.y, p4)), float(roc_auc_score(target.y, p5))
        reported = prior_e4 if name == "E4" else prior_e5["variants"][name.split("/")[1]]
        if not (abs(a4 - float(reported["M4"]["AUROC"])) <= 1e-12
                and abs(a5 - float(reported["M5"]["AUROC"])) <= 1e-12):
            raise AssertionError("Frozen model rescore differs from prior aggregate")
        predictions[name] = (p4, p5)
        results[name] = {"M4_AUROC": a4, "M5_AUROC": a5,
                         "M5_minus_M4_AUROC": a5 - a4,
                         "point_parity_with_prior_aggregate": True}

    draws = {name: np.empty(DRAW_COUNT, float) for name in predictions}
    gene_draws = np.empty(DRAW_COUNT, float)
    for repetition, ix in enumerate(donor_bootstrap_indices(target.y, target.donor)):
        y = target.y[ix]
        gene_draws[repetition] = roc_auc_score(y, gene[ix])
        for name, (p4, p5) in predictions.items():
            draws[name][repetition] = (roc_auc_score(y, p5[ix])
                                       - roc_auc_score(y, p4[ix]))
    for name in results:
        results[name]["M5_minus_M4_CI95_percentile"] = (
            np.quantile(draws[name], [.025, .975]).astype(float).tolist())
    write_aggregate(output, {
        "analysis": "E1_E4_E5_posthoc_aggregate_diagnostics",
        "status": "complete_descriptive",
        "requested_after_original_results_seen": True,
        "bootstrap_draws": DRAW_COUNT,
        "seed_sequence": SEED,
        "E1_final_pooled": {
            "absolute_gene_SMD_GEO_training_controls_vs_PPMI_controls": pooled_smd,
            "signed_gene_component_AUROC_PPMI": e1_point,
            "signed_gene_component_AUROC_CI95_percentile": (
                np.quantile(gene_draws, [.025, .975]).astype(float).tolist()),
        },
        "E4_E5": results,
    })
    print("E1_E4_E5_POSTHOC_COMPLETE: aggregate only")


if __name__ == "__main__":
    main()
