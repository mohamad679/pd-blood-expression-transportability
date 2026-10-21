#!/usr/bin/env python3
"""E10: same-month CBC concordance for all locked PPMI PD and control donors."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from audit_ppmi_cbc_coverage import canonical_patno
from audit_ppmi_collection_linkage import parse_month_year
from run_ppmi_postlock_secondary import read_collection_months
from score_ppmi_after_lock import load_locked_input
from postlock_common import (private_output, require_execution_authorization,
                             rng_for, synthetic_ranks, write_aggregate)

# The PPMI hematology code map is fixed before E10 execution. Lymphocyte
# percentage is nonspecific for T/B; RBC/Hb cannot validate reticulocytes.
COMPARISONS = (
    ("T_cells_vs_lymphocyte_percent", 0, "HMT16"),
    ("B_lineage_vs_lymphocyte_percent", 1, "HMT16"),
    ("Monocytic_vs_monocytes_absolute", 2, "HMT10"),
    ("Monocytic_vs_monocytes_percent", 2, "HMT17"),
    ("Neutrophils_vs_neutrophils_absolute", 3, "HMT8"),
    ("Neutrophils_vs_neutrophils_percent", 3, "HMT15"),
    ("Erythroid_vs_RBC", 4, "HMT3"),
    ("Erythroid_vs_hemoglobin", 4, "HMT40"),
    ("Platelet_vs_platelet_count", 5, "HMT13"),
)


def read_cbc(cbc: Path, donors: np.ndarray, collection_months: dict):
    wanted = set(map(str, donors))
    codes = {item[2] for item in COMPARISONS}
    rows = defaultdict(list)
    with cbc.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"PATNO", "EVENT_ID", "LTSTCODE", "LCOLLDT", "COLLTM",
                    "LSIRES", "LSIUNIT"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CBC export lacks predeclared columns")
        for row in reader:
            donor = canonical_patno(row.get("PATNO"))
            code = (row.get("LTSTCODE") or "").strip().upper()
            if (donor not in wanted or code not in codes
                    or collection_months.get(donor) is None
                    or parse_month_year(row.get("LCOLLDT")) != collection_months[donor]):
                continue
            event = (row.get("EVENT_ID") or "").strip().upper()
            tier = 0 if event == "BL" else 1 if event == "SC" else 2
            try:
                value = float((row.get("LSIRES") or "").replace(",", "").strip())
                if not np.isfinite(value):
                    value = None
            except ValueError:
                value = None
            rows[donor, code].append((tier, (row.get("LCOLLDT") or "").strip(),
                                      (row.get("COLLTM") or "").strip(), value,
                                      (row.get("LSIUNIT") or "").strip()))
    return rows


def select_assay(rows: list[tuple]) -> tuple[float | None, str | None, str]:
    if not rows:
        return None, None, "missing"
    top = min(row[0] for row in rows)
    highest = [row for row in rows if row[0] == top]
    if len({(r[1], r[2]) for r in highest}) != 1 or len({(r[3], r[4]) for r in highest}) != 1:
        return None, None, "conflicting_highest_tier"
    value, unit = highest[0][3], highest[0][4]
    if value is None or not unit:
        return None, None, "non_numeric_or_blank_unit"
    return value, unit, "usable"


def compare(ranks: np.ndarray, y: np.ndarray, donors: np.ndarray,
            marker_indices: list[np.ndarray], collection_months: dict,
            cbc_rows: dict, draws: int = 5000) -> dict:
    output = {"analysis": "E10", "status": "complete_descriptive",
              "same_month_any_visit": True, "n_candidates": len(donors),
              "comparison_count": len(COMPARISONS), "comparisons": {}}
    for comparison_ix, (name, panel, code) in enumerate(COMPARISONS):
        sample_ix, values, units = [], [], []
        excluded = defaultdict(int)
        for i, donor in enumerate(map(str, donors)):
            if collection_months.get(donor) is None:
                excluded["unknown_collection_month"] += 1
                continue
            value, unit, state = select_assay(cbc_rows.get((donor, code), []))
            if state == "usable":
                sample_ix.append(i)
                values.append(value)
                units.append(unit)
            else:
                excluded[state] += 1
        unit_set = sorted(set(units))
        row = {"panel": geo.PANELS[panel], "assay_code": code,
               "n_before_unit_check": len(sample_ix), "excluded": dict(excluded),
               "SI_units": unit_set}
        if len(unit_set) != 1:
            row["status"] = "not_evaluable_missing_or_mixed_SI_units"
            output["comparisons"][name] = row
            continue
        ix = np.asarray(sample_ix, int)
        value = np.asarray(values, float)
        panel_score = ranks[ix][:, marker_indices[panel]].mean(axis=1)
        row["n"] = len(ix)
        row["n_PD"] = int(np.sum(y[ix] == 1))
        row["n_control"] = int(np.sum(y[ix] == 0))
        if len(ix) < 3 or len(np.unique(value)) < 2 or len(np.unique(panel_score)) < 2:
            row["status"] = "not_evaluable_too_few_or_constant"
            output["comparisons"][name] = row
            continue
        rho = float(spearmanr(panel_score, value).statistic)
        rng = rng_for(10, comparison_ix)
        boot = []
        for _ in range(draws):
            pick = rng.integers(0, len(ix), len(ix))
            if len(np.unique(value[pick])) > 1 and len(np.unique(panel_score[pick])) > 1:
                boot.append(float(spearmanr(panel_score[pick], value[pick]).statistic))
        row.update(status="complete_descriptive", Spearman=rho,
                   bootstrap_valid_draws=len(boot),
                   CI95_percentile=(np.quantile(boot, [.025, .975]).tolist()
                                    if len(boot) >= .9*draws else None))
        output["comparisons"][name] = row
    return output


def synthetic_self_test() -> None:
    import tempfile
    ranks = synthetic_ranks(40, 30, 71)
    donors = np.array([str(i) for i in range(40)])
    y = np.array([0, 1] * 20)
    markers = [np.arange(i*3, i*3+3) for i in range(6)]
    months = {donor: (2026, 10) for donor in donors}
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "cbc.csv"
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["PATNO", "EVENT_ID", "LTSTCODE",
                                                       "LCOLLDT", "COLLTM", "LSIRES", "LSIUNIT"])
            writer.writeheader()
            for donor in donors:
                for _name, _panel, code in COMPARISONS:
                    writer.writerow({"PATNO": donor, "EVENT_ID": "BL", "LTSTCODE": code,
                                     "LCOLLDT": "10/2026", "COLLTM": "08:00",
                                     "LSIRES": float(donor)+1, "LSIUNIT": "unit"})
        rows = read_cbc(path, donors, months)
        result = compare(ranks, y, donors, markers, months, rows, draws=10)
        assert result["comparison_count"] == 9
        assert all(row["n"] == 40 for row in result["comparisons"].values())
    print("E10_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--collection-csv", type=Path)
    parser.add_argument("--cbc-csv", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(v is None for v in (args.authorization_record, args.bundle, args.binding,
                                args.collection_csv, args.cbc_csv, args.output)):
        parser.error("Real E10 requires authorization and all locked inputs")
    require_execution_authorization(args.authorization_record, analysis="E10")
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    ranks, _age, _sex, y, donors, genes = load_locked_input(args.bundle, binding)
    manifest = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    months = read_collection_months(args.collection_csv, donors)
    rows = read_cbc(args.cbc_csv, donors, months)
    write_aggregate(output, compare(ranks, y, donors, markers, months, rows))
    print("E10_COMPLETE: private aggregate output written")


if __name__ == "__main__":
    main()
