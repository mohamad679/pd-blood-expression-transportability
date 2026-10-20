"""Aggregate-only collection timing audit for frozen PPMI candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


BUNDLE_SHA256 = "1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def patient(value: object) -> str:
    result = str(value or "").strip()
    return result[3:] if result.upper().startswith("PP-") else result


def month(value: object) -> tuple[int, int] | None:
    match = re.fullmatch(r"(0?[1-9]|1[0-2])/(\d{4})", str(value or "").strip())
    return (int(match.group(2)), int(match.group(1))) if match else None


def read_candidates(bundle: Path) -> dict[tuple[str, str], str]:
    with zipfile.ZipFile(bundle) as archive:
        lines = archive.read("candidate_metadata.jsonl").decode("utf-8").splitlines()
    candidates = {}
    for line in lines:
        row = json.loads(line)
        key = (patient(row["PATNO"]), str(row["EVENT_ID"]).strip().upper())
        group = row["historical_label"]
        if not key[0] or key[1] != "BL" or group not in {"PD", "Control"} or key in candidates:
            raise ValueError("Frozen candidate labels, visit, or uniqueness failed validation")
        candidates[key] = group
    if len(candidates) != 472 or Counter(candidates.values()) != {"PD": 311, "Control": 161}:
        raise ValueError("Frozen candidate bundle does not match its aggregate manifest")
    return candidates


def read_collection(path: Path) -> tuple[dict[tuple[str, str], list[tuple[int, int] | None]], int]:
    records: dict[tuple[str, str], list[tuple[int, int] | None]] = defaultdict(list)
    count = 0
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"PATNO", "EVENT_ID", "BLDCLLDT"}.issubset(reader.fieldnames or []):
            raise ValueError("Collection export lacks a required column")
        for row in reader:
            count += 1
            pid = patient(row.get("PATNO"))
            event = str(row.get("EVENT_ID") or "").strip().upper()
            if pid and event:
                records[(pid, event)].append(month(row.get("BLDCLLDT")))
    return records, count


def audit(candidates: dict[tuple[str, str], str], records: dict[tuple[str, str], list[tuple[int, int] | None]]) -> dict:
    groups: dict[str, Counter[str]] = {"PD": Counter(), "Control": Counter()}
    for key, group in candidates.items():
        hits = records.get(key, [])
        counts = groups[group]
        counts["candidates"] += 1
        if not hits:
            counts["no_collection_row"] += 1
        elif len(hits) == 1:
            counts["one_row_parseable_month"] += hits[0] is not None
            counts["one_row_missing_or_unparseable_month"] += hits[0] is None
        else:
            counts["multiple_rows"] += 1
            if all(value is not None for value in hits) and len(set(hits)) == 1:
                counts["multiple_rows_same_parseable_month"] += 1
            else:
                counts["multiple_rows_conflicting_or_missing_month"] += 1
    return {"status": "metadata_only_collection_timing_complete_medication_pending",
            "join": "exact normalized PATNO plus baseline EVENT_ID",
            "date_precision": "MM/YYYY; no exact day inference",
            "groups": {group: dict(sorted(counts.items())) for group, counts in groups.items()}}


def self_test() -> None:
    candidates = {("1", "BL"): "PD", ("2", "BL"): "PD", ("3", "BL"): "Control"}
    records = {("1", "BL"): [(2020, 1)], ("2", "BL"): [(2020, 1), (2020, 2)],
               ("3", "BL"): [None]}
    result = audit(candidates, records)
    assert result["groups"]["PD"]["one_row_parseable_month"] == 1
    assert result["groups"]["PD"]["multiple_rows_conflicting_or_missing_month"] == 1
    assert result["groups"]["Control"]["one_row_missing_or_unparseable_month"] == 1
    assert month("02/2020") == (2020, 2) and month("02/20") is None
    print("SYNTHETIC_SELF_TEST_PASSED")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--collection", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.bundle or not args.collection or not args.output:
        parser.error("--bundle, --collection, and --output are required")
    start = time.monotonic()
    try:
        print("Stage 1/3: validating frozen candidate bundle", flush=True)
        if digest(args.bundle) != BUNDLE_SHA256:
            raise ValueError("Candidate bundle hash differs from frozen input")
        candidates = read_candidates(args.bundle)
        print("Stage 2/3: counting collection linkage", flush=True)
        records, collection_rows = read_collection(args.collection)
        result = audit(candidates, records)
        result["input_sha256"] = {"bundle": BUNDLE_SHA256, "collection": digest(args.collection)}
        result["collection_total_rows"] = collection_rows
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Stage 3/3: aggregate report saved; elapsed {time.monotonic()-start:.1f}s", flush=True)
        print("PPMI_CANDIDATE_COLLECTION_TIMING_AUDIT_COMPLETE", flush=True)
        return 0
    except Exception as error:
        print(f"PPMI_CANDIDATE_COLLECTION_TIMING_AUDIT_FAILED: {type(error).__name__}; no private row details printed", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
