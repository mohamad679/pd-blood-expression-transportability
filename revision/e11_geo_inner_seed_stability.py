#!/usr/bin/env python3
"""Post hoc, GEO-only audit of M5 tuning across shared inner-fold seeds."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from e03_ppmi_nested_cv import fit_pair, score_model
from e04_transductive_universe import EligiblePrep
from e05_weighted_deviance_panel_exclusion import panel_excluded_indices

CONDITIONS = (
    ("auroc_full", "AUROC", False),
    ("deviance_full", "weighted_deviance", False),
    ("auroc_panel_excluded", "AUROC", True),
)


def quartiles(values: list[float]) -> dict:
    if not values:
        return {"median": None, "range": None, "IQR": None}
    a = np.asarray(values, float)
    return {"median": float(np.median(a)),
            "range": [float(np.min(a)), float(np.max(a))],
            "IQR": np.quantile(a, [.25, .75]).astype(float).tolist()}


def summarize(rows: list[dict]) -> dict:
    summary = {}
    for source in sorted({row["source"] for row in rows}):
        by_source = [row for row in rows if row["source"] == source]
        summary[source] = {"conditions": {}, "paired_differences": {}}
        for condition, _, _ in CONDITIONS:
            cells = [row["conditions"][condition] for row in by_source]
            successes = [cell for cell in cells if cell["status"] == "complete"]
            failures = len(cells) - len(successes)
            categories = ["no_gene" if cell["M5_reuses_M4"]
                          else format(cell["selected_lambda"], ".12g")
                          for cell in successes]
            counts = Counter(categories)
            no_gene = counts["no_gene"]
            deltas = [cell["delta_AUROC"] for cell in successes]
            summary[source]["conditions"][condition] = {
                "attempts": len(cells), "successes": len(successes),
                "failures": failures,
                "failure_types": dict(Counter(cell["error_type"] for cell in cells
                                               if cell["status"] == "failed")),
                "lambda_category_counts": dict(sorted(counts.items())),
                "rho_counts": dict(sorted(Counter(
                    "no_gene" if cell["M5_reuses_M4"]
                    else format(cell["selected_rho"], ".12g")
                    for cell in successes).items())),
                "no_gene_count": no_gene,
                "no_gene_fraction_of_attempts": no_gene / len(cells),
                "no_gene_stable_by_predeclared_rule": failures <= 2 and no_gene >= 16,
                "lambda_stable_by_predeclared_rule": (
                    failures <= 2 and max(counts.values(), default=0) >= 16),
                "delta_AUROC": quartiles(deltas),
                "delta_at_least_plus_0_05_count": sum(d >= .05 for d in deltas),
                "positive_direction_rule_met": (
                    failures <= 2 and sum(d >= .05 for d in deltas) >= 16
                    and bool(deltas) and np.median(deltas) >= .05),
            }
        for left, right in (("deviance_full", "auroc_full"),
                            ("auroc_panel_excluded", "auroc_full"),
                            ("deviance_full", "auroc_panel_excluded")):
            paired = [row["conditions"][left]["delta_AUROC"]
                      - row["conditions"][right]["delta_AUROC"]
                      for row in by_source
                      if row["conditions"][left]["status"] == "complete"
                      and row["conditions"][right]["status"] == "complete"]
            summary[source]["paired_differences"][f"{left}_minus_{right}"] = {
                "paired_seeds": len(paired), **quartiles(paired)}
    summary["both_directions_stable_positive_increment"] = {
        condition: all(summary[source]["conditions"][condition]
                       ["positive_direction_rule_met"]
                       for source in summary if source != "both_directions_stable_positive_increment")
        for condition, _, _ in CONDITIONS}
    return summary


def run(output: Path, seeds: int = 20) -> None:
    if seeds != 20:
        raise ValueError("The hashed addendum fixes exactly 20 seeds")
    if output.exists():
        raise ValueError("Use a fresh GEO-only output directory")
    data, input_hashes = geo.load_geo()
    excluded_eligible = panel_excluded_indices(data)
    output.mkdir(parents=True)
    started = time.monotonic()
    rows = []
    sources = sorted(set(map(str, data.source)))
    for source in sources:
        train = data.take(np.flatnonzero(data.source == source))
        target = data.take(np.flatnonzero(data.source != source))
        for seed_index in range(seeds):
            scope = f"E11/{source}/seed_{seed_index}/inner"
            # Construct once, then reuse the same folds in all three conditions.
            inner = list(geo.donor_folds(train, 5, scope))
            row = {"source": source, "target_source": str(target.source[0]),
                   "seed_index": seed_index, "inner_scope": scope,
                   "conditions": {}}
            for condition, objective, exclude_panels in CONDITIONS:
                factory = (lambda part, max_genes: EligiblePrep(
                    part, excluded_eligible, max_genes)) if exclude_panels else None
                folder = output / source / f"seed_{seed_index:02d}" / condition
                fit_scope = f"E11/{source}/seed_{seed_index}/{condition}"
                try:
                    prep, m4, m5, decision = fit_pair(
                        train, inner, folder, fit_scope, prep_factory=factory,
                        objective=objective)
                    p4 = score_model(prep, m4, target, "M4")
                    p5 = p4 if m5 is None else score_model(prep, m5, target, "M5")
                    selected = decision["selected"]
                    row["conditions"][condition] = {
                        "status": "complete",
                        "selected_lambda": None if selected["baseline"]
                                           else float(selected["lambda"]),
                        "selected_rho": None if selected["baseline"]
                                        else float(selected["rho"]),
                        "M5_reuses_M4": bool(decision["M5_reuses_M4"]),
                        "persistent_boundary": bool(decision["persistent_boundary"]),
                        "M4_AUROC": float(roc_auc_score(target.y, p4)),
                        "M5_AUROC": float(roc_auc_score(target.y, p5)),
                        "delta_AUROC": float(roc_auc_score(target.y, p5)
                                             - roc_auc_score(target.y, p4)),
                    }
                except Exception as error:
                    row["conditions"][condition] = {
                        "status": "failed", "error_type": type(error).__name__}
            folder = output / source / f"seed_{seed_index:02d}"
            (folder / "aggregate.json").write_text(
                json.dumps(row, indent=2, sort_keys=True, allow_nan=False),
                encoding="utf-8")
            rows.append(row)
            completed = len(rows)
            total = len(sources) * seeds
            elapsed = time.monotonic() - started
            eta = elapsed / completed * (total - completed)
            print(f"E11 GEO-only {completed}/{total} seed-directions | "
                  f"elapsed {elapsed/3600:.1f}h | ETA ~{eta/3600:.1f}h",
                  flush=True)
    result = {"analysis": "E11_GEO_inner_seed_stability", "post_hoc": True,
              "GEO_only": True, "PPMI_scored": False, "seeds_per_direction": seeds,
              "shared_inner_folds_across_conditions": True,
              "geo_input_sha256": input_hashes, "summary": summarize(rows),
              "per_seed_aggregate": rows}
    (output / "E11_GEO_AGGREGATE.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8")
    print(f"E11_GEO_ONLY_COMPLETE: {output}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
