#!/usr/bin/env python3
"""E9: locked-prediction sporadic de novo subgroup, conditional on source audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from score_ppmi_after_lock import load_locked_input
from postlock_common import (percentile_interval, private_output,
                             require_execution_authorization, stratified_indices,
                             write_aggregate)


def subgroup_result(y: np.ndarray, p4: np.ndarray, p5: np.ndarray,
                    sporadic_de_novo: np.ndarray, draws: int = 5000) -> dict:
    y, p4, p5, mask = (np.asarray(y), np.asarray(p4, float),
                       np.asarray(p5, float), np.asarray(sporadic_de_novo, bool))
    if (y.ndim != 1 or any(a.shape != y.shape for a in (p4, p5, mask))
            or not np.isin(y, [0, 1]).all()
            or not np.isfinite(p4).all() or not np.isfinite(p5).all()):
        raise ValueError("Aligned finite locked predictions required")
    sy, s4, s5 = y[mask], p4[mask], p5[mask]
    counts = {"Control": int((sy == 0).sum()), "PD": int((sy == 1).sum())}
    result = {"analysis": "E9", "counts": counts, "draws": draws,
              "fixed_models": True, "exploratory": True}
    if min(counts.values()) < 10:
        result["status"] = "not_evaluable_fewer_than_10_per_class"
        return result
    d = stratified_indices(sy, draws, 9)
    deltas = np.array([roc_auc_score(sy[ix], s5[ix]) - roc_auc_score(sy[ix], s4[ix])
                       for ix in d], float)
    result.update(status="complete_exploratory",
                  M4_AUROC=float(roc_auc_score(sy, s4)),
                  M5_AUROC=float(roc_auc_score(sy, s5)),
                  M5_minus_M4=float(roc_auc_score(sy, s5)-roc_auc_score(sy, s4)),
                  M5_minus_M4_CI95_percentile=percentile_interval(deltas))
    return result


def read_status(path: Path, donors: np.ndarray) -> np.ndarray:
    values = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"PATNO", "sporadic_de_novo"}.issubset(reader.fieldnames or []):
            raise ValueError("Status file needs PATNO and sporadic_de_novo columns")
        for row in reader:
            donor = str(row["PATNO"]).strip()
            value = str(row["sporadic_de_novo"]).strip().lower()
            if donor in values or value not in {"true", "false"}:
                raise ValueError("Duplicate donor or uncertain subgroup status")
            values[donor] = value == "true"
    if not set(map(str, donors)).issubset(values):
        raise ValueError("Subgroup status is incomplete")
    return np.array([values[str(d)] for d in donors], bool)


def synthetic_self_test() -> None:
    rng = np.random.default_rng(99)
    y = np.array([0, 1] * 25)
    p4 = rng.uniform(.1, .9, len(y))
    p5 = rng.uniform(.1, .9, len(y))
    mask = np.ones(len(y), bool)
    result = subgroup_result(y, p4, p5, mask, draws=30)
    assert result["status"] == "complete_exploratory"
    assert len(result["M5_minus_M4_CI95_percentile"]) == 2
    mask[:40] = False
    assert subgroup_result(y, p4, p5, mask, draws=10)["status"].startswith("not_evaluable")
    print("E9_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--eligibility-audit", type=Path)
    parser.add_argument("--status-csv", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--locked-predictions", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(v is None for v in (args.authorization_record, args.eligibility_audit,
                                args.status_csv, args.bundle, args.binding,
                                args.locked_predictions, args.output)):
        parser.error("Real E9 requires authorization and every named source")
    require_execution_authorization(args.authorization_record, analysis="E9")
    output = private_output(args.output)
    audit = json.loads(args.eligibility_audit.read_text(encoding="utf-8"))
    if audit.get("genetic_PD_included") is not True or audit.get("genetic_status_reproducible") is not True:
        write_aggregate(output, {"analysis": "E9", "status": "not_evaluable_genetic_status"})
        print("E9_NOT_EVALUABLE: source-data criterion was not met")
        return
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    _ranks, _age, _sex, y, donors, _genes = load_locked_input(args.bundle, binding)
    with np.load(args.locked_predictions, allow_pickle=False) as z:
        py, p4, p5 = z["y"], z["M4"], z["M5"]
    if not np.array_equal(y, py):
        raise ValueError("Locked prediction row order differs from bound input")
    mask = read_status(args.status_csv, donors)
    write_aggregate(output, subgroup_result(y, p4, p5, mask))
    print("E9_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
