#!/usr/bin/env python3
"""Aggregate stored E6 replicate metrics without any model fit or PPMI rescore."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np


def distribution(values: list[float], no_gene: list[bool]) -> dict:
    a = np.asarray(values, float)
    if len(a) != 200 or not np.isfinite(a).all() or len(no_gene) != len(a):
        raise ValueError("E6 must have 200 finite paired replicate differences")
    zeros = a == 0.0
    nogene = np.asarray(no_gene, bool)
    return {
        "mean": float(np.mean(a)),
        "median": float(np.median(a)),
        "IQR": np.quantile(a, [.25, .75]).astype(float).tolist(),
        "range": [float(np.min(a)), float(np.max(a))],
        "count_delta_below_zero": int(np.sum(a < 0)),
        "fraction_delta_below_zero": float(np.mean(a < 0)),
        "count_delta_exactly_zero": int(np.sum(zeros)),
        "fraction_delta_exactly_zero": float(np.mean(zeros)),
        "count_no_gene": int(np.sum(nogene)),
        "fraction_no_gene": float(np.mean(nogene)),
        "zero_delta_and_no_gene_count": int(np.sum(zeros & nogene)),
        "zero_delta_without_no_gene_count": int(np.sum(zeros & ~nogene)),
        "no_gene_without_zero_delta_count": int(np.sum(~zeros & nogene)),
        "empirical_percentile_of_minus_0_069_fraction_at_or_below":
            float(np.mean(a <= -0.069)),
        "empirical_percentile_of_minus_0_069_percent_at_or_below":
            float(100 * np.mean(a <= -0.069)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Use a fresh output path; do not overwrite earlier results")
    final = json.loads((args.input / "E6_AGGREGATE.json").read_text())
    if final.get("status") != "complete_exploratory" or final.get("attempts") != 200 \
            or final.get("successful_attempts") != 200 or final.get("failed_attempts") != 0:
        raise ValueError("E6 final aggregate and error ledger must be complete")
    paths = sorted(args.input.glob("attempt_*/attempt_aggregate.json"))
    rows = [json.loads(path.read_text()) for path in paths]
    if len(rows) != 200 or sorted(row["repetition"] for row in rows) != list(range(200)) \
            or any(row["status"] != "complete" for row in rows):
        raise ValueError("E6 attempt ledger is incomplete")
    rows.sort(key=lambda row: row["repetition"])
    no_gene = [bool(row["M5_selection"]["M5_reuses_M4"]) for row in rows]
    result = {
        "analysis": "E6_stored_replicate_posthoc_aggregate",
        "posthoc_after_E6_complete": True,
        "fit_or_rescore_performed": False,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "attempts": len(rows),
        "threshold_minus_0_069_is_rounded_primary_delta": True,
        "percentile_definition": "fraction of replicate deltas <= -0.069",
        "training_only": distribution(
            [row["M5_minus_M4_AUROC_training_only"] for row in rows], no_gene),
        "training_plus_PPMI": distribution(
            [row["M5_minus_M4_AUROC_training_plus_PPMI"] for row in rows], no_gene),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True),
                           encoding="utf-8")
    print(f"E6_POSTHOC_AGGREGATE_COMPLETE: {args.output}")


if __name__ == "__main__":
    main()
