#!/usr/bin/env python3
"""Count CBC visit/date coverage for frozen baseline controls, across all visits.

Only PATNO, EVENT_ID and LCOLLDT are retained from CBC exports. Measurement
values, gene expression and scores are not inspected or used. Output is
aggregate metadata only. A same-month row is availability, not proof of an
available numeric result for any particular marker.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

try:
    from .audit_ppmi_baseline_control_cbc import eligible_baseline_controls
    from .audit_ppmi_cbc_coverage import canonical_patno
    from .audit_ppmi_collection_linkage import parse_month_year
except ImportError:
    from audit_ppmi_baseline_control_cbc import eligible_baseline_controls
    from audit_ppmi_cbc_coverage import canonical_patno
    from audit_ppmi_collection_linkage import parse_month_year


def read_baseline_collection(path: Path, eligible: set[tuple[str, str]]) -> tuple[dict, Counter]:
    """Use baseline dates; optionally count BLDCOLL without interpreting it."""
    dates = {}
    codes: Counter[str] = Counter()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"PATNO", "EVENT_ID", "BLDCLLDT"}.issubset(reader.fieldnames or []):
            raise ValueError("Collection export lacks participant/visit/date fields")
        for row in reader:
            pair = (canonical_patno(row.get("PATNO")), (row.get("EVENT_ID") or "").strip().upper())
            if pair not in eligible:
                continue
            if pair in dates:
                raise ValueError("Ambiguous baseline collection rows")
            dates[pair] = parse_month_year(row.get("BLDCLLDT"))
            if "BLDCOLL" in row:
                codes[(row["BLDCOLL"] or "").strip() or "blank"] += 1
    return dates, codes


def audit_table(path: Path, baseline_months: dict[str, tuple[int, int] | None]) -> dict:
    """Project metadata columns and count people, without using CBC results."""
    row_counts: Counter[str] = Counter()
    participants: dict[str, set[str]] = defaultdict(set)
    same_month_rows: Counter[str] = Counter()
    same_month_people: dict[str, set[str]] = defaultdict(set)
    month_precision: Counter[str] = Counter()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        required = ("PATNO", "EVENT_ID", "LCOLLDT")
        if not set(required).issubset(header):
            raise ValueError("CBC export lacks participant/visit/date fields")
        indices = [header.index(name) for name in required]
        for fields in reader:
            if not fields:
                continue
            if len(fields) != len(header):
                raise ValueError("CBC export contains a malformed CSV row")
            participant_raw, event_raw, date_raw = (fields[i] for i in indices)
            patno = canonical_patno(participant_raw)
            if patno not in baseline_months:
                continue
            event = event_raw.strip().upper() or "blank"
            month = parse_month_year(date_raw)
            row_counts[event] += 1
            participants[event].add(patno)
            month_precision["month_year" if month else "blank_or_other"] += 1
            if month is not None and month == baseline_months[patno]:
                same_month_rows[event] += 1
                same_month_people[event].add(patno)

    all_people = set().union(*participants.values()) if participants else set()
    all_same = set().union(*same_month_people.values()) if same_month_people else set()
    same_bl = same_month_people.get("BL", set())
    same_sc = same_month_people.get("SC", set())
    return {
        "eligible_participants_with_any_row_any_visit": len(all_people),
        "eligible_participants_with_any_row_same_BL_collection_month": len(all_same),
        "same_month_BL_priority_participants": len(same_bl),
        "same_month_SC_priority_participants_without_BL": len(same_sc - same_bl),
        "same_month_other_visit_only_participants": len(all_same - same_bl - same_sc),
        "date_precision_row_counts": dict(sorted(month_precision.items())),
        "by_EVENT_ID": {
            event: {
                "CBC_table_rows": row_counts[event],
                "unique_eligible_participants": len(participants[event]),
                "same_BL_collection_month_rows": same_month_rows[event],
                "same_BL_collection_month_unique_participants": len(same_month_people[event]),
            }
            for event in sorted(row_counts)
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rna-metadata", "age-at-visit", "participant-status", "whole-blood-archive", "cbc-current", "cbc-archived", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    eligible = eligible_baseline_controls(args.rna_metadata, args.age_at_visit, args.participant_status)
    if {event for _, event in eligible} != {"BL"}:
        raise ValueError("Frozen control set does not consist solely of baseline visits")
    if len({patno for patno, _ in eligible}) != len(eligible):
        raise ValueError("Frozen control set has repeated participants")
    collection_dates, code_counts = read_baseline_collection(args.whole_blood_archive, eligible)
    baseline_months = {patno: collection_dates.get((patno, "BL")) for patno, _ in eligible}
    result = {
        "purpose": "CBC metadata availability across every visit of the frozen baseline healthy controls",
        "eligible_baseline_healthy_controls": len(eligible),
        "controls_with_BL_collection_month": sum(month is not None for month in baseline_months.values()),
        "CBC_columns_retained": ["PATNO", "EVENT_ID", "LCOLLDT"],
        "measurement_values_used": False,
        "expression_or_scores_used": False,
        "BLDCOLL_codes_uninterpreted_not_used_for_eligibility": dict(sorted(code_counts.items())),
        "current": audit_table(args.cbc_current, baseline_months),
        "archived": audit_table(args.cbc_archived, baseline_months),
        "interpretation": "Availability only. Per-marker numeric n and record ambiguities remain unevaluated; no association or marker selection performed.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
