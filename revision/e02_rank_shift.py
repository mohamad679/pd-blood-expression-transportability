#!/usr/bin/env python3
"""E2: label-free rank-shift diagnostic; synthetic test only in phase 2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from score_ppmi_after_lock import load_locked_input
from postlock_common import (private_output, require_execution_authorization,
                             summary, synthetic_ranks, write_aggregate)


def quantiles(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    return {"median": float(np.median(values)),
            "q25": float(np.quantile(values, .25)),
            "q75": float(np.quantile(values, .75)),
            "q05": float(np.quantile(values, .05)),
            "q95": float(np.quantile(values, .95)),
            "absolute_median": float(np.median(np.abs(values))),
            "absolute_max": float(np.max(np.abs(values)))}


def rank_shift(geo_ranks: np.ndarray, geo_source: np.ndarray,
               ppmi_ranks: np.ndarray, genes: np.ndarray,
               selected_nonzero_ids: set[str]) -> dict:
    geo_ranks = np.asarray(geo_ranks, float)
    ppmi_ranks = np.asarray(ppmi_ranks, float)
    genes = np.asarray(genes)
    if (geo_ranks.ndim != 2 or ppmi_ranks.ndim != 2
            or geo_ranks.shape[1] != ppmi_ranks.shape[1]
            or geo_ranks.shape[1] != len(genes)
            or geo_source.shape != (len(geo_ranks),)
            or not np.isfinite(geo_ranks).all() or not np.isfinite(ppmi_ranks).all()):
        raise ValueError("Aligned finite rank matrices and source vector are required")
    selected = np.array([str(g) in selected_nonzero_ids for g in genes])
    if selected.sum() != len(selected_nonzero_ids):
        raise ValueError("Nonzero M5 genes not fully contained in common universe")
    if selected.sum() == 0:
        return {"analysis": "E2", "status": "not_evaluable_M5_gene_block_zero",
                "label_free": True, "common_genes": int(len(genes)),
                "nonzero_M5_genes": 0}
    mean_ppmi = ppmi_ranks.mean(axis=0)
    output = {"analysis": "E2", "status": "complete_exploratory",
              "label_free": True, "common_genes": int(len(genes)),
              "nonzero_M5_genes": int(selected.sum()), "by_source": {}}
    groups = {str(s): geo_ranks[geo_source == s] for s in np.unique(geo_source)}
    groups["GEO_pooled"] = geo_ranks
    for name, matrix in groups.items():
        geo_mean = matrix.mean(axis=0)
        rho = float(spearmanr(geo_mean, mean_ppmi).statistic)
        if not np.isfinite(rho):
            raise ValueError("Undefined whole-universe rank correlation")
        shift = mean_ppmi[selected] - geo_mean[selected]
        output["by_source"][name] = {"n_GEO": int(len(matrix)),
                                      "n_PPMI": int(len(ppmi_ranks)),
                                      "whole_universe_spearman": rho,
                                      "nonzero_gene_delta_rank": quantiles(shift)}
    return output


def nonzero_m5_gene_ids(final_dir: Path) -> set[str]:
    with np.load(final_dir / "preprocessing.npz", allow_pickle=False) as z:
        genes = z["genes"]
        selected = z["selected"]
        n_baseline = len(z["keep"])
    with np.load(final_dir / "M5_model.npz", allow_pickle=False) as z:
        coef = z["coef"]
    if len(coef) != n_baseline + len(selected):
        raise ValueError("Frozen M5 coefficient layout is inconsistent")
    return {str(g) for g in genes[selected[np.flatnonzero(coef[n_baseline:] != 0.)]]}


def synthetic_self_test() -> None:
    genes = np.array([f"G{i:03d}" for i in range(30)])
    geo_ranks = synthetic_ranks(20, 30, 11)
    ppmi_ranks = synthetic_ranks(12, 30, 12)
    source = np.array(["A"] * 10 + ["B"] * 10)
    result = rank_shift(geo_ranks, source, ppmi_ranks, genes,
                        {"G001", "G004", "G013"})
    assert result["nonzero_M5_genes"] == 3
    assert set(result["by_source"]) == {"A", "B", "GEO_pooled"}
    assert all(-1 <= row["whole_universe_spearman"] <= 1
               for row in result["by_source"].values())
    print("E2_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--geo-final", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(x is None for x in (args.authorization_record, args.bundle, args.binding,
                               args.geo_final, args.output)):
        parser.error("Real run requires authorization, bundle, binding, GEO final fit and output")
    require_execution_authorization(args.authorization_record, analysis="E2")
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    ppmi_ranks, _age, _sex, _labels, _donors, ppmi_genes = load_locked_input(args.bundle, binding)
    geo_data, _ = geo.load_geo()
    if not np.array_equal(geo_data.genes, ppmi_genes):
        raise ValueError("GEO/PPMI common gene order differs")
    result = rank_shift(geo_data.X, geo_data.source, ppmi_ranks, ppmi_genes,
                        nonzero_m5_gene_ids(args.geo_final))
    write_aggregate(output, result)
    print("E2_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
