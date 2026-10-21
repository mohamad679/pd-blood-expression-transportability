"""Read locked PPMI candidate metadata and emit aggregate descriptions only.

No participant identifier, candidate row, expression value, or prediction is printed
or written by this script. The output path must remain under ~/PPMI_private.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import Counter
from pathlib import Path
from zipfile import ZipFile


PRIVATE = (Path.home() / "PPMI_private").resolve()


def numeric(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def clean_category(value):
    if value is None or str(value).strip() == "":
        return "<missing>"
    return str(value).strip()


def stats(values):
    values = [value for value in values if value is not None]
    return {"n_available": len(values), "mean": statistics.mean(values) if values else None,
            "sd": statistics.stdev(values) if len(values) > 1 else None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--collection", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(PRIVATE):
        raise ValueError("Aggregate output must stay under ~/PPMI_private")
    by_group = {"PD": [], "Control": []}
    with ZipFile(PRIVATE / "PPMI_WINDOWS_PRIVATE_BUNDLE.zip") as archive:
        with archive.open("candidate_metadata.jsonl") as stream:
            for line in stream:
                row = json.loads(line)
                group = row.get("historical_label")
                if group not in by_group:
                    raise ValueError("Unexpected historical diagnosis category")
                by_group[group].append(row)
    result = {"source": "locked candidate_metadata.jsonl member of PPMI_WINDOWS_PRIVATE_BUNDLE.zip",
              "scope": "selected 472 locked candidates; not a pre-QC denominator",
              "groups": {}}
    for group, rows in by_group.items():
        portal = [r.get("portal_metadata") or {} for r in rows]
        result["groups"][group] = {
            "n": len(rows),
            "age_at_BL_visit": stats([numeric(r.get("age_at_BL_visit")) for r in rows]),
            "recorded_sex": dict(Counter(clean_category(r.get("sex")) for r in rows)),
            "portal_study_arm": dict(Counter(clean_category(p.get("Study Arm")) for p in portal)),
            "portal_genetic_status": dict(Counter(clean_category(p.get("Genetic Status")) for p in portal)),
            "portal_diagnosis": dict(Counter(clean_category(p.get("Diagnosis")) for p in portal)),
            "portal_UPDRS3_score": stats([numeric(p.get("UPDRS3 score")) for p in portal]),
            "portal_age_at_diagnosis": stats([numeric(p.get("Age at diagnosis")) for p in portal]),
            "not_available_in_source": ["Hoehn and Yahr stage", "LEDD at draw", "exact collection day"],
        }
    if sum(x["n"] for x in result["groups"].values()) != 472:
        raise ValueError("Locked candidate count mismatch")
    if args.collection:
        candidate_keys = set()
        for rows in by_group.values():
            for row in rows:
                pid = str(row["PATNO"] or "").strip()
                if pid.upper().startswith("PP-"):
                    pid = pid[3:]
                candidate_keys.add((pid, str(row["EVENT_ID"]).strip().upper()))
        collection_rows = {}
        with args.collection.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not {"PATNO", "EVENT_ID", "BLDCLLDT"}.issubset(reader.fieldnames or []):
                raise ValueError("Collection export lacks required columns")
            for row in reader:
                pid = str(row.get("PATNO") or "").strip()
                if pid.upper().startswith("PP-"):
                    pid = pid[3:]
                key = (pid, str(row.get("EVENT_ID") or "").strip().upper())
                if key in candidate_keys:
                    collection_rows.setdefault(key, []).append(str(row.get("BLDCLLDT") or "").strip())
        reason_counts = Counter()
        for key in candidate_keys:
            values = collection_rows.get(key, [])
            if not values:
                reason_counts["no_matching_baseline_collection_row"] += 1
            elif len(values) > 1:
                reason_counts["multiple_matching_rows"] += 1
            elif not values[0]:
                reason_counts["one_row_blank_month"] += 1
            elif not re.fullmatch(r"(0?[1-9]|1[0-2])/(\d{4})", values[0]):
                reason_counts["one_row_nonblank_unparseable_month"] += 1
            else:
                reason_counts["one_row_parseable_month"] += 1
        result["collection_month_reasons"] = dict(reason_counts)
        result["collection_source_basename"] = args.collection.name
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    output.chmod(0o600)
    print("PPMI_AGGREGATE_METADATA_COMPLETE: 472 candidates; output remains private")


if __name__ == "__main__":
    main()
