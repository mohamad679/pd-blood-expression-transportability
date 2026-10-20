"""Aggregate-only descriptive LEDD medication timing audit for PPMI candidates."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from audit_ppmi_candidate_collection_timing import (
    BUNDLE_SHA256,
    digest,
    month,
    patient,
    read_candidates,
    read_collection,
)


LEDD_SHA256 = "cdab712e45d2a57d4e81b27fa414b40f228a9d7216edf31d55891beb69b4e10b"
COLLECTION_SHA256 = "e609f9adb835e3654a5de2d5c02169549e095cec1d36a1cc7365b7d1027080bc"


def read_ledd(path: Path) -> tuple[dict[str, list[tuple[bool, tuple[int, int] | None, tuple[int, int] | None, bool]]], int]:
    records = defaultdict(list)
    row_count = 0
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"PATNO", "LEDTRT", "STARTDT", "STOPDT"}.issubset(reader.fieldnames or []):
            raise ValueError("LEDD export lacks a required column")
        for row in reader:
            row_count += 1
            pid = patient(row.get("PATNO"))
            if not pid:
                continue
            stop_text = (row.get("STOPDT") or "").strip()
            records[pid].append((bool((row.get("LEDTRT") or "").strip()),
                                 month(row.get("STARTDT")), month(stop_text), not stop_text))
    return records, row_count


def classify_record(
    treatment_named: bool,
    start: tuple[int, int] | None,
    stop: tuple[int, int] | None,
    stop_blank: bool,
    draw: tuple[int, int],
) -> str:
    if not treatment_named or start is None or (stop is None and not stop_blank):
        return "indeterminate_invalid_or_missing_boundary"
    if stop is not None and stop < start:
        return "indeterminate_inconsistent_interval"
    if start > draw or (stop is not None and stop < draw):
        return "not_ongoing_by_dates"
    if start == draw or stop == draw:
        return "indeterminate_same_month_boundary"
    if start < draw and (stop_blank or (stop is not None and stop > draw)):
        return "recorded_ongoing_use"
    raise AssertionError("Date classification is incomplete")


def audit(candidates: dict[tuple[str, str], str],
          collection: dict[tuple[str, str], list[tuple[int, int] | None]],
          ledd: dict[str, list[tuple[bool, tuple[int, int] | None, tuple[int, int] | None, bool]]]) -> dict:
    groups = {"PD": Counter(), "Control": Counter()}
    record_outcomes = Counter()
    for key, group in candidates.items():
        out = groups[group]
        out["candidates"] += 1
        dates = collection.get(key, [])
        if len(dates) != 1 or dates[0] is None:
            out["collection_month_unknown"] += 1
            continue
        meds = ledd.get(key[0], [])
        out["collection_month_known"] += 1
        if not meds:
            out["no_ledd_log_rows"] += 1
            continue
        states = [classify_record(*record, dates[0]) for record in meds]
        record_outcomes.update(states)
        if "recorded_ongoing_use" in states:
            out["recorded_ongoing_ledd_use"] += 1
        elif any(state.startswith("indeterminate") for state in states):
            out["ledd_timing_indeterminate"] += 1
        else:
            out["no_recorded_ongoing_ledd_use"] += 1
    for group, counts in groups.items():
        assert counts["collection_month_known"] + counts["collection_month_unknown"] == counts["candidates"]
        assert (counts["no_ledd_log_rows"] + counts["recorded_ongoing_ledd_use"] +
                counts["ledd_timing_indeterminate"] + counts["no_recorded_ongoing_ledd_use"] ==
                counts["collection_month_known"])
    return {"status": "descriptive_ledd_timing_complete_no_exclusions",
            "scope": "recorded dopaminergic medication in LEDD log only; not all medications",
            "candidate_exclusions_due_to_medication": 0,
            "groups": {group: dict(sorted(counts.items())) for group, counts in groups.items()},
            "medication_record_timing_counts_with_known_collection_month": dict(sorted(record_outcomes.items())),
            "interpretation_limit": "No LEDD row does not establish untreated status; MM/YYYY dates cannot resolve same-month boundaries"}


def self_test() -> None:
    draw = (2020, 6)
    assert classify_record(True, (2020, 5), None, True, draw) == "recorded_ongoing_use"
    assert classify_record(True, (2020, 5), (2020, 7), False, draw) == "recorded_ongoing_use"
    assert classify_record(True, (2020, 6), None, True, draw) == "indeterminate_same_month_boundary"
    assert classify_record(True, (2020, 5), (2020, 6), False, draw) == "indeterminate_same_month_boundary"
    assert classify_record(True, (2020, 7), None, True, draw) == "not_ongoing_by_dates"
    assert classify_record(True, (2020, 5), (2020, 4), False, draw) == "indeterminate_inconsistent_interval"
    candidates = {("1", "BL"): "PD", ("2", "BL"): "PD", ("3", "BL"): "Control"}
    collection = {("1", "BL"): [draw], ("2", "BL"): [None], ("3", "BL"): [draw]}
    ledd = {"1": [(True, (2020, 5), None, True)]}
    result = audit(candidates, collection, ledd)
    assert result["groups"]["PD"]["recorded_ongoing_ledd_use"] == 1
    assert result["groups"]["PD"]["collection_month_unknown"] == 1
    assert result["groups"]["Control"]["no_ledd_log_rows"] == 1
    print("SYNTHETIC_SELF_TEST_PASSED")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--collection", type=Path)
    parser.add_argument("--ledd", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not all((args.bundle, args.collection, args.ledd, args.output)):
        parser.error("--bundle, --collection, --ledd, and --output are required")
    started = time.monotonic()
    try:
        print("Stage 1/4: validating frozen input hashes", flush=True)
        if (digest(args.bundle) != BUNDLE_SHA256 or digest(args.collection) != COLLECTION_SHA256 or
                digest(args.ledd) != LEDD_SHA256):
            raise ValueError("Frozen input hash mismatch")
        print("Stage 2/4: reading candidate and collection metadata", flush=True)
        candidates = read_candidates(args.bundle)
        collection, collection_rows = read_collection(args.collection)
        print("Stage 3/4: classifying LEDD date intervals", flush=True)
        ledd, ledd_rows = read_ledd(args.ledd)
        result = audit(candidates, collection, ledd)
        result["input_sha256"] = {"bundle": BUNDLE_SHA256, "collection": COLLECTION_SHA256, "ledd": LEDD_SHA256}
        result["input_row_counts"] = {"collection": collection_rows, "ledd": ledd_rows}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Stage 4/4: aggregate report saved; elapsed {time.monotonic()-started:.1f}s", flush=True)
        print("PPMI_LEDD_AT_COLLECTION_AUDIT_COMPLETE", flush=True)
        return 0
    except Exception as error:
        print(f"PPMI_LEDD_AT_COLLECTION_AUDIT_FAILED: {type(error).__name__}; no private row details printed", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
