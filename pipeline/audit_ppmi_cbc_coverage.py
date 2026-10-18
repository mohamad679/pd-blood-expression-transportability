#!/usr/bin/env python3
"""Aggregate-only audit of PPMI long-RNA sample coverage by CBC visit.

The script reads restricted input files locally and writes aggregate counts only.
It does not inspect case/control labels, expression values, or calculate marker
scores or associations. PATNO and visit are matched exactly; no fuzzy matching.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


MARKERS = {
    "HMT2": "hematocrit",
    "HMT3": "RBC",
    "HMT7": "WBC",
    "HMT8": "neutrophils_absolute",
    "HMT9": "lymphocytes_absolute",
    "HMT13": "platelets",
    "HMT15": "neutrophils_percent",
    "HMT16": "lymphocytes_percent",
    "HMT40": "hemoglobin",
}


def canonical_patno(value: str | None) -> str | None:
    """Normalize only plain numeric PATNO or the verified PP-<digits> form."""
    value = (value or "").strip()
    match = re.fullmatch(r"(?:PP-)?([0-9]+)", value, flags=re.IGNORECASE)
    return match.group(1) if match else None


def sample_pair(row: dict[str, str]) -> tuple[str, str] | None:
    patno = canonical_patno(row.get("PATNO"))
    patno_visit = (row.get("PATNO Visit") or "").strip()
    match = re.fullmatch(r"([0-9]+)_([A-Za-z][A-Za-z0-9]*)", patno_visit)
    if not patno or not match or patno != match.group(1):
        return None
    return patno, match.group(2).upper()


def has_numeric_result(row: dict[str, str]) -> bool:
    for column in ("LSIRES", "LUSRES"):
        value = (row.get(column) or "").strip()
        if value:
            try:
                float(value.replace(",", ""))
                return True
            except ValueError:
                continue
    return False


def read_samples(path: Path) -> tuple[list[tuple[str, str]], dict[str, int]]:
    samples: list[tuple[str, str]] = []
    counts = Counter()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"PATNO", "PATNO Visit", "Long RNA-seq"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"RNA metadata missing required columns: {sorted(required)}")
        for row in reader:
            if (row.get("Long RNA-seq") or "").strip().casefold() != "yes":
                continue
            counts["long_rna_rows"] += 1
            pair = sample_pair(row)
            if pair is None:
                counts["unmatched_or_unparseable_sample_rows"] += 1
            else:
                samples.append(pair)
    return samples, dict(counts)


def read_cbc(path: Path) -> dict:
    pairs: set[tuple[str, str]] = set()
    numeric_pairs: set[tuple[str, str]] = set()
    marker_pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    marker_records: dict[tuple[str, ...], Counter] = defaultdict(Counter)
    marker_units: dict[str, set[tuple[str, str]]] = defaultdict(set)
    row_counts = Counter()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "PATNO", "EVENT_ID", "LTSTCODE", "LSIRES", "LUSRES", "LCOLLDT", "COLLTM"
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"CBC file missing required columns: {sorted(required)}")
        row_counts["data_rows"] = 0
        for row in reader:
            row_counts["data_rows"] += 1
            patno = canonical_patno(row.get("PATNO"))
            event = (row.get("EVENT_ID") or "").strip().upper()
            if not patno or not event:
                row_counts["unmatched_or_unparseable_cbc_rows"] += 1
                continue
            pair = (patno, event)
            pairs.add(pair)
            if has_numeric_result(row):
                numeric_pairs.add(pair)
                code = (row.get("LTSTCODE") or "").strip().upper()
                if code in MARKERS:
                    marker_pairs[code].add(pair)
                    record_key = (
                        patno,
                        event,
                        code,
                        (row.get("LCOLLDT") or "").strip(),
                        (row.get("COLLTM") or "").strip(),
                    )
                    signature = (
                        (row.get("LSIRES") or "").strip(),
                        (row.get("LSIUNIT") or "").strip(),
                        (row.get("LUSRES") or "").strip(),
                        (row.get("LUSUNIT") or "").strip(),
                    )
                    marker_records[record_key][signature] += 1
                    marker_units[code].add((signature[1], signature[3]))
    return {
        "pairs": pairs,
        "numeric_pairs": numeric_pairs,
        "marker_pairs": marker_pairs,
        "marker_records": marker_records,
        "marker_units": marker_units,
        "data_rows": row_counts["data_rows"],
        "unmatched_rows": row_counts["unmatched_or_unparseable_cbc_rows"],
    }


def compare_marker_versions(current: dict, archived: dict) -> dict:
    current_records = current["marker_records"]
    archived_records = archived["marker_records"]
    shared_keys = current_records.keys() & archived_records.keys()
    identical = 0
    different = 0
    ambiguous = 0
    per_code = defaultdict(lambda: Counter())
    for key in shared_keys:
        code = key[2]
        left = current_records[key]
        right = archived_records[key]
        if sum(left.values()) != 1 or sum(right.values()) != 1:
            ambiguous += 1
            per_code[code]["shared_keys_with_multiple_rows"] += 1
        if left == right:
            identical += 1
            per_code[code]["identical_value_and_unit_keys"] += 1
        else:
            different += 1
            per_code[code]["different_value_or_unit_keys"] += 1
    return {
        "scope": "numeric rows for the prespecified hematology marker test codes; matched by exact PATNO, EVENT_ID, LTSTCODE, collection date and time",
        "current_unique_marker_measurement_keys": len(current_records),
        "archived_unique_marker_measurement_keys": len(archived_records),
        "shared_marker_measurement_keys": len(shared_keys),
        "shared_keys_with_identical_value_and_units": identical,
        "shared_keys_with_different_value_or_units": different,
        "shared_keys_ambiguous_due_to_duplicate_rows": ambiguous,
        "marker_code_details": {
            code: {
                "current_si_us_unit_pairs": sorted([list(x) for x in current["marker_units"].get(code, set())]),
                "archived_si_us_unit_pairs": sorted([list(x) for x in archived["marker_units"].get(code, set())]),
                **dict(per_code.get(code, {})),
            }
            for code in sorted(MARKERS)
        },
    }


def sample_coverage(samples: list[tuple[str, str]], cbc: dict) -> dict[str, int]:
    sample_pairs = set(samples)
    matched = sample_pairs & cbc["pairs"]
    numeric = sample_pairs & cbc["numeric_pairs"]
    result = {
        "sample_rows_with_valid_patno_visit": len(samples),
        "unique_sample_patno_visit_pairs": len(sample_pairs),
        "unique_sample_participants": len({patno for patno, _ in sample_pairs}),
        "pairs_with_any_cbc_record": len(matched),
        "pairs_with_any_numeric_cbc_result": len(numeric),
    }
    result["sample_rows_with_any_cbc_record"] = sum(pair in cbc["pairs"] for pair in samples)
    result["sample_rows_with_any_numeric_cbc_result"] = sum(pair in cbc["numeric_pairs"] for pair in samples)
    result["marker_pair_coverage"] = {
        MARKERS[code]: len(sample_pairs & pairs)
        for code, pairs in sorted(cbc["marker_pairs"].items())
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-metadata", type=Path, required=True)
    parser.add_argument("--cbc-current", type=Path, required=True)
    parser.add_argument("--cbc-archived", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples, sample_counts = read_samples(args.rna_metadata)
    current = read_cbc(args.cbc_current)
    archived = read_cbc(args.cbc_archived)
    sample_pairs = set(samples)
    current_sample_pairs = sample_pairs & current["pairs"]
    archived_sample_pairs = sample_pairs & archived["pairs"]
    combined = {
        "pairs": current["pairs"] | archived["pairs"],
        "numeric_pairs": current["numeric_pairs"] | archived["numeric_pairs"],
        "marker_pairs": defaultdict(set),
    }
    for source in (current, archived):
        for code, pairs in source["marker_pairs"].items():
            combined["marker_pairs"][code].update(pairs)

    result = {
        "purpose": "CBC-to-long-RNA metadata coverage only; no outcome/score/association analysis",
        "matching": "exact canonical PATNO + exact visit code from PATNO Visit; no fuzzy matching",
        "temporal_resolution": "visit-level only; same-day alignment is not assessed by this audit",
        "sample_metadata": sample_counts,
        "current_cbc": {
            "data_rows": current["data_rows"],
            "unique_participant_visit_pairs": len(current["pairs"]),
            "unmatched_or_unparseable_rows": current["unmatched_rows"],
            "coverage": sample_coverage(samples, current),
        },
        "archived_cbc": {
            "data_rows": archived["data_rows"],
            "unique_participant_visit_pairs": len(archived["pairs"]),
            "unmatched_or_unparseable_rows": archived["unmatched_rows"],
            "coverage": sample_coverage(samples, archived),
        },
        "version_overlap_at_sample_visit_level": {
            "matched_in_both_versions": len(current_sample_pairs & archived_sample_pairs),
            "matched_only_in_current_version": len(current_sample_pairs - archived_sample_pairs),
            "matched_only_in_archived_version": len(archived_sample_pairs - current_sample_pairs),
        },
        "marker_version_concordance": compare_marker_versions(current, archived),
        "union_for_coverage_audit_only": {
            "unique_participant_visit_pairs": len(combined["pairs"]),
            "coverage": sample_coverage(samples, combined),
        },
        "marker_test_codes": MARKERS,
        "disclaimer": "This is not a disease-transfer or external-validation analysis; only visit-level measurement availability is summarized.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Aggregate-only CBC coverage audit written to {args.output}")


if __name__ == "__main__":
    main()
