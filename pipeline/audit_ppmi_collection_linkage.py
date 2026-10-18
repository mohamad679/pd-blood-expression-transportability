#!/usr/bin/env python3
"""Aggregate-only exact visit crosswalk for PPMI whole-blood collection and CBC.

Restricted PPMI files are read locally. The output contains aggregate counts
only: no participant identifiers, event-level rows, dates, or clinical values.
This script does not read disease labels or gene expression and computes no
marker score or association.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

try:  # Support both direct script execution and package-style test imports.
    from .audit_ppmi_cbc_coverage import MARKERS, canonical_patno, has_numeric_result, read_samples
except ImportError:  # pragma: no cover - exercised by the command-line entry point
    from audit_ppmi_cbc_coverage import MARKERS, canonical_patno, has_numeric_result, read_samples


def parse_month_year(value: str | None) -> tuple[int, int] | None:
    """Parse the month/year precision present in the supplied PPMI exports."""
    value = (value or "").strip()
    match = re.fullmatch(r"(0?[1-9]|1[0-2])/(\d{4})", value)
    if not match:
        return None
    return int(match.group(2)), int(match.group(1))


def read_collection(path: Path, rna_pairs: set[tuple[str, str]]) -> dict:
    required = {"PATNO", "EVENT_ID", "BLDCOLL", "BLDCLLDT", "INFODT"}
    row_count = 0
    valid_rows = 0
    status_counts: Counter[str] = Counter()
    date_precision: Counter[str] = Counter()
    records: dict[tuple[str, str], tuple[str, tuple[int, int] | None]] = {}
    duplicate_pairs = 0
    overlap_status: Counter[str] = Counter()
    overlap_dates: Counter[str] = Counter()

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Whole-blood collection file missing required columns: {sorted(required)}")
        for row in reader:
            row_count += 1
            patno = canonical_patno(row.get("PATNO"))
            event = (row.get("EVENT_ID") or "").strip().upper()
            if not patno or not event:
                continue
            valid_rows += 1
            pair = (patno, event)
            status = (row.get("BLDCOLL") or "").strip() or "blank"
            collection_month = parse_month_year(row.get("BLDCLLDT"))
            date_precision["month_year" if collection_month else "blank" if not (row.get("BLDCLLDT") or "").strip() else "other"] += 1
            status_counts[status] += 1
            if pair in records:
                duplicate_pairs += 1
            else:
                records[pair] = (status, collection_month)
            if pair in rna_pairs:
                overlap_status[status] += 1
                overlap_dates["month_year" if collection_month else "blank_or_other"] += 1

    return {
        "rows": row_count,
        "valid_rows": valid_rows,
        "unique_participant_visit_pairs": len(records),
        "duplicate_pair_rows_after_first": duplicate_pairs,
        "status_code_counts": dict(sorted(status_counts.items())),
        "collection_date_precision_counts": dict(sorted(date_precision.items())),
        "rna_visit_overlap_pairs": set(records) & rna_pairs,
        "collection_month_by_pair": {pair: value[1] for pair, value in records.items()},
        "overlap_status_code_counts": dict(sorted(overlap_status.items())),
        "overlap_collection_date_precision_counts": dict(sorted(overlap_dates.items())),
    }


def read_cbc_month_coverage(path: Path, overlap_pairs: set[tuple[str, str]], collection_months: dict) -> dict:
    required = {"PATNO", "EVENT_ID", "LTSTCODE", "LSIRES", "LUSRES", "LCOLLDT"}
    pairs_with_marker: set[tuple[str, str]] = set()
    same_month_pairs: set[tuple[str, str]] = set()
    date_precision: Counter[str] = Counter()
    rows = 0

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"CBC file missing required columns: {sorted(required)}")
        for row in reader:
            rows += 1
            patno = canonical_patno(row.get("PATNO"))
            event = (row.get("EVENT_ID") or "").strip().upper()
            code = (row.get("LTSTCODE") or "").strip().upper()
            if not patno or not event or code not in MARKERS or not has_numeric_result(row):
                continue
            pair = (patno, event)
            if pair not in overlap_pairs:
                continue
            month = parse_month_year(row.get("LCOLLDT"))
            date_precision["month_year" if month else "blank_or_other"] += 1
            if month is None:
                continue
            pairs_with_marker.add(pair)
            if month == collection_months.get(pair):
                same_month_pairs.add(pair)

    return {
        "rows": rows,
        "numeric_core_marker_rows_for_overlapping_visits_date_precision": dict(sorted(date_precision.items())),
        "rna_visit_pairs_with_numeric_core_marker": len(pairs_with_marker),
        "rna_visit_pairs_with_marker_in_same_reported_month": len(same_month_pairs),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-metadata", type=Path, required=True)
    parser.add_argument("--whole-blood-archive", type=Path, required=True)
    parser.add_argument("--cbc-current", type=Path, required=True)
    parser.add_argument("--cbc-archived", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples, sample_counts = read_samples(args.rna_metadata)
    rna_pairs = set(samples)
    collection = read_collection(args.whole_blood_archive, rna_pairs)
    overlap_pairs = collection["rna_visit_overlap_pairs"]
    collection_months = collection["collection_month_by_pair"]
    result = {
        "purpose": "Aggregate exact PPMI RNA visit to whole-blood collection and CBC date-precision audit only",
        "matching": "Exact normalized PATNO + EVENT_ID/visit; no fuzzy matching",
        "privacy": "No individual identifiers, event rows, dates, labels, gene expression, scores, or associations are written",
        "rna_sample_metadata": {
            "long_rna_rows": sample_counts.get("long_rna_rows", 0),
            "valid_rows": len(samples),
            "unique_patno_visit_pairs": len(rna_pairs),
            "unmatched_or_unparseable_rows": sample_counts.get("unmatched_or_unparseable_sample_rows", 0),
        },
        "whole_blood_collection_archive": {
            key: value
            for key, value in collection.items()
            if key not in {"rna_visit_overlap_pairs", "collection_month_by_pair"}
        },
        "cbc_current": read_cbc_month_coverage(args.cbc_current, overlap_pairs, collection_months),
        "cbc_archived": read_cbc_month_coverage(args.cbc_archived, overlap_pairs, collection_months),
        "interpretation_limit": (
            "The supplied collection and CBC date fields are month/year precision. If both fields represent collection dates "
            "as their names imply, same-month pairs are within 30 days; exact day differences and nearest-record selection "
            "across adjacent months cannot be determined. "
            "No date-based eligibility rule is declared passed by this audit."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Aggregate-only collection crosswalk audit written to {args.output}")


if __name__ == "__main__":
    main()
