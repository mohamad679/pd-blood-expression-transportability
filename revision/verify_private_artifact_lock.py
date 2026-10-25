#!/usr/bin/env python3
"""Hash-only pre-execution comparison against the previous lock inventory.

No participant rows are parsed or printed. ZIP members are streamed through
SHA-256, and the output contains aggregate counts and artifact names only on
failure. Two post-lock prediction files have local baseline hashes, not public
reference hashes, and are explicitly counted separately.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def digest_member(path: Path, member: str) -> str:
    h = hashlib.sha256()
    with zipfile.ZipFile(path) as archive, archive.open(member) as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def contains_string(value: object, expected: str) -> bool:
    if isinstance(value, dict):
        return any(contains_string(child, expected) for child in value.values())
    if isinstance(value, list):
        return any(contains_string(child, expected) for child in value)
    return value == expected


def run(inventory: Path, output: Path) -> dict:
    if output.exists():
        raise ValueError("Use a fresh lock-audit output file")
    with inventory.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 39:
        raise ValueError("Expected the 39-item expanded Phase-1 lock inventory")
    failures = []
    reference_count = reported_count = local_count = 0
    for row in rows:
        name = row["artifact"]
        path = Path(row["source_path"])
        expected = row["expected_sha256"]
        try:
            if row["status"] == "REPORTED_MATCH":
                reported_count += 1
                value = json.loads(path.read_text(encoding="utf-8"))
                match = bool(expected) and contains_string(value, expected)
            else:
                actual = (digest_member(path, row["zip_member"])
                          if row["zip_member"] else digest_file(path))
                if expected:
                    reference_count += 1
                    match = actual == expected
                else:
                    local_count += 1
                    match = actual == row["actual_sha256"]
            if not match:
                failures.append(name)
        except (OSError, KeyError, ValueError, zipfile.BadZipFile):
            failures.append(name)
    receipt = {
        "kind": "pre_postlock_private_artifact_hash_audit_v1",
        "result": "PASS" if not failures and (reference_count, reported_count, local_count) == (34, 3, 2) else "FAIL",
        "reference_hash_match": not failures and reference_count == 34,
        "verified_file_count": len(rows) - len(failures),
        "reference_sha256_match_count": reference_count,
        "previously_reported_hash_count": reported_count,
        "local_baseline_only_count": local_count,
        "limitation": "Three hashes come from the prior aggregate report; two post-lock prediction files lack a public reference hash and are compared only with the previous local inventory.",
        "inventory_sha256": digest_file(inventory),
        "mismatches": failures,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.inventory, args.output)
    print("PRIVATE_LOCK_AUDIT:", result["result"],
          "reference", result["reference_sha256_match_count"],
          "reported", result["previously_reported_hash_count"],
          "local_baseline_only", result["local_baseline_only_count"],
          "failures", len(result["mismatches"]))
    if result["result"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
