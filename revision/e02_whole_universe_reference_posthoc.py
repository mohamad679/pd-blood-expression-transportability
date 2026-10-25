#!/usr/bin/env python3
"""Aggregate, label-free reference for E2's selected-gene rank shift.

This post-primary, post-E2 diagnostic was requested after E2 was inspected. It
does not change the timestamped amendment or the original E2 output.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from score_ppmi_after_lock import load_locked_input
from postlock_common import private_output, write_aggregate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    ppmi, _age, _sex, _labels, _donors, genes = load_locked_input(args.bundle, binding)
    development, _ = geo.load_geo()
    if not np.array_equal(genes, development.genes):
        raise ValueError("Frozen GEO and PPMI common-gene order differs")
    if not np.isfinite(ppmi).all() or not np.isfinite(development.X).all():
        raise ValueError("Nonfinite rank input")

    groups = {str(source): development.X[development.source == source]
              for source in np.unique(development.source)}
    groups["GEO_pooled"] = development.X
    mean_ppmi = ppmi.mean(axis=0)
    rows = {}
    for source, matrix in groups.items():
        displacement = np.abs(mean_ppmi - matrix.mean(axis=0))
        rows[source] = {
            "n_GEO": int(len(matrix)),
            "n_PPMI": int(len(ppmi)),
            "n_genes": int(len(genes)),
            "absolute_delta_median": float(np.median(displacement)),
            "absolute_delta_q25": float(np.quantile(displacement, .25)),
            "absolute_delta_q75": float(np.quantile(displacement, .75)),
        }
    write_aggregate(output, {
        "analysis": "E2_whole_universe_reference_posthoc",
        "status": "complete_descriptive",
        "label_free": True,
        "requested_after_E2_result_seen": True,
        "sources": rows,
    })
    print("E2_WHOLE_UNIVERSE_REFERENCE_COMPLETE: aggregate only")


if __name__ == "__main__":
    main()
