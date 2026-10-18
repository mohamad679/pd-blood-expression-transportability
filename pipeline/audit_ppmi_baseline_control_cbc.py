#!/usr/bin/env python3
"""Aggregate-only same-month CBC availability in the frozen 161 baseline controls.

The candidate definition reproduces the prior aggregate clinical join. PPMI
files are read locally; output contains counts and code-value frequencies only.
No expression, sample-level outcomes, scores, correlations, or participant-level
rows are read or written.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

try:
    from .audit_ppmi_cbc_coverage import MARKERS, canonical_patno, has_numeric_result, sample_pair
    from .audit_ppmi_collection_linkage import parse_month_year
except ImportError:  # pragma: no cover - direct script execution
    from audit_ppmi_cbc_coverage import MARKERS, canonical_patno, has_numeric_result, sample_pair
    from audit_ppmi_collection_linkage import parse_month_year


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def eligible_baseline_controls(metadata_path: Path, age_path: Path, status_path: Path) -> set[tuple[str, str]]:
    metadata = rows(metadata_path)
    ages = rows(age_path)
    statuses = rows(status_path)

    status_index: dict[str, dict[str, str]] = {}
    for record in statuses:
        patno = canonical_patno(record.get("PATNO"))
        if not patno or patno in status_index:
            raise ValueError("Participant status keys are missing or nonunique")
        status_index[patno] = record

    age_values: dict[tuple[str, str], set[float]] = {}
    for record in ages:
        patno = canonical_patno(record.get("PATNO"))
        event = (record.get("EVENT_ID") or "").strip()
        if not patno or not event:
            continue
        value = (record.get("AGE_AT_VISIT") or "").strip()
        try:
            age = float(value)
        except ValueError:
            continue
        if 0 < age < 120:
            age_values.setdefault((patno, event), set()).add(age)

    eligible: set[tuple[str, str]] = set()
    for record in metadata:
        if not (
            record.get("Disease Status") == "Healthy Control"
            and record.get("Clinical Event") == "BL"
            and record.get("Month") == "M0"
            and (record.get("Long RNA-seq") or "").strip().casefold() == "yes"
            and record.get("QC") == "HighQC"
        ):
            continue
        patno = canonical_patno(record.get("PATNO"))
        if not patno:
            continue
        status = status_index.get(patno)
        if not status or status.get("COHORT_DEFINITION") != "Healthy Control":
            continue
        if len(age_values.get((patno, "BL"), set())) != 1:
            continue
        if record.get("Sex") not in {"Male", "Female"}:
            continue
        if record.get("Diagnosis") != "Control" or record.get("Case Control") != "Control":
            continue
        pair = sample_pair(record)
        if pair:
            eligible.add(pair)
    if len(eligible) != 161:
        raise ValueError(f"Expected the previously audited 161 baseline controls, found {len(eligible)}")
    return eligible


def read_collection(path: Path, eligible: set[tuple[str, str]]) -> tuple[dict, Counter]:
    required = {"PATNO", "EVENT_ID", "BLDCOLL", "BLDCLLDT"}
    dates: dict[tuple[str, str], tuple[int, int] | None] = {}
    codes: Counter[str] = Counter()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Whole-blood collection export lacks required columns")
        for record in reader:
            patno = canonical_patno(record.get("PATNO"))
            event = (record.get("EVENT_ID") or "").strip().upper()
            pair = (patno, event)
            if patno and event and pair in eligible:
                if pair in dates:
                    raise ValueError("Duplicate collection records at an eligible baseline visit")
                dates[pair] = parse_month_year(record.get("BLDCLLDT"))
                codes[(record.get("BLDCOLL") or "").strip() or "blank"] += 1
    return dates, codes


def read_cbc(path: Path, eligible: set[tuple[str, str]], collection_months: dict) -> dict:
    required = {"PATNO", "EVENT_ID", "LTSTCODE", "LSIRES", "LUSRES", "LCOLLDT"}
    any_row_pairs: set[tuple[str, str]] = set()
    marker_pairs: set[tuple[str, str]] = set()
    same_month_by_code: dict[str, set[tuple[str, str]]] = {code: set() for code in MARKERS}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CBC export lacks required columns")
        for record in reader:
            patno = canonical_patno(record.get("PATNO"))
            event = (record.get("EVENT_ID") or "").strip().upper()
            pair = (patno, event)
            if not patno or not event or pair not in eligible:
                continue
            any_row_pairs.add(pair)
            code = (record.get("LTSTCODE") or "").strip().upper()
            if code not in MARKERS or not has_numeric_result(record):
                continue
            marker_pairs.add(pair)
            month = parse_month_year(record.get("LCOLLDT"))
            if month and month == collection_months.get(pair):
                same_month_by_code[code].add(pair)

    same_month_any = set().union(*same_month_by_code.values())
    return {
        "eligible_controls_with_any_CBC_row_at_exact_visit": len(any_row_pairs),
        "eligible_controls_with_numeric_core_marker_at_exact_visit": len(marker_pairs),
        "eligible_controls_with_numeric_core_marker_in_same_reported_month": len(same_month_any),
        "same_month_numeric_marker_pair_counts": {
            MARKERS[code]: len(same_month_by_code[code]) for code in sorted(MARKERS)
        },
        "same_month_both_absolute_neutrophil_and_lymphocyte": len(
            same_month_by_code["HMT8"] & same_month_by_code["HMT9"]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-metadata", type=Path, required=True)
    parser.add_argument("--age-at-visit", type=Path, required=True)
    parser.add_argument("--participant-status", type=Path, required=True)
    parser.add_argument("--whole-blood-archive", type=Path, required=True)
    parser.add_argument("--cbc-current", type=Path, required=True)
    parser.add_argument("--cbc-archived", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    eligible = eligible_baseline_controls(args.rna_metadata, args.age_at_visit, args.participant_status)
    collection_dates, code_counts = read_collection(args.whole_blood_archive, eligible)
    collection_months = {pair: month for pair, month in collection_dates.items() if month}
    result = {
        "purpose": "Aggregate same-month CBC availability among the previously audited baseline healthy-control set",
        "eligibility_definition": "2021 RNA metadata: healthy control, BL/M0, Long RNA-seq Yes, HighQC; exact 2026 participant-status cohort; one valid BL age; sex recorded; portal Diagnosis and Case Control both Control",
        "matching": "Exact normalized PATNO plus visit/EVENT_ID; no fuzzy or adjacent-visit matching",
        "privacy": "Aggregate counts and code-value frequencies only; no IDs, dates, expression, outcomes, scores or correlations written",
        "eligible_baseline_healthy_controls": len(eligible),
        "whole_blood_collection_archive": {
            "eligible_controls_with_exact_visit_record": len(collection_dates),
            "eligible_controls_with_month_year_collection_date": len(collection_months),
            "BLDCOLL_code_counts_uninterpreted": dict(sorted(code_counts.items())),
        },
        "cbc_current": read_cbc(args.cbc_current, eligible, collection_months),
        "cbc_archived": read_cbc(args.cbc_archived, eligible, collection_months),
        "temporal_rule": "Exact participant+visit and same reported calendar month only; adjacent months excluded. Day-level same-day or nearest-record selection is not possible at month/year precision.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Aggregate-only baseline-control CBC audit written to {args.output}")


if __name__ == "__main__":
    main()
