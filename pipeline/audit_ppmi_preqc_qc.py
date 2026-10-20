"""Mac-local, metadata-only PPMI pre-QC audit. Aggregate output only."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


def clean(value: object) -> str:
    return " ".join(str(value or "").strip().lower().replace("_", " ").split())


def patient(value: object) -> str:
    result = str(value or "").strip()
    return result[3:] if result.upper().startswith("PP-") else result


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        header = stream.readline()
        delimiter = "\t" if header.count("\t") > header.count(",") else ","
        stream.seek(0)
        reader = csv.DictReader(stream, delimiter=delimiter)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Input table has missing or duplicate column names")
        return list(reader)


def require_columns(rows: list[dict[str, str]], columns: set[str], name: str) -> None:
    if not rows:
        raise ValueError(f"{name} has no rows")
    missing = columns - set(rows[0])
    if missing:
        raise ValueError(f"{name} lacks required columns: {', '.join(sorted(missing))}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def portal_group(row: dict[str, str]) -> str | None:
    triple = tuple(clean(row.get(k)) for k in ("Disease Status", "Diagnosis", "Case Control"))
    if triple == ("idiopathic pd", "pd", "case"):
        return "PD"
    if triple == ("healthy control", "control", "control"):
        return "Control"
    return None


def age_valid(values: list[str]) -> bool:
    ages = set()
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number) and 0 < number < 120:
            ages.add(number)
    return len(ages) == 1


def candidate_keys(bundle: Path) -> set[tuple[str, str]]:
    with zipfile.ZipFile(bundle) as archive:
        lines = archive.read("candidate_metadata.jsonl").decode("utf-8").splitlines()
    keys = set()
    for line in lines:
        row = json.loads(line)
        key = (patient(row["PATNO"]), str(row["release_metadata"]["Sample"]).strip())
        if not all(key) or key in keys:
            raise ValueError("Candidate bundle has empty or duplicate participant/sample keys")
        keys.add(key)
    if len(keys) != 472:
        raise ValueError("Candidate bundle size differs from the frozen provisional set")
    return keys


def audit(
    release: list[dict[str, str]],
    portal: list[dict[str, str]],
    age: list[dict[str, str]],
    candidates: set[tuple[str, str]],
) -> dict:
    require_columns(release, {"PATNO", "Sample", "HudAlphaID", "CLINICAL_EVENT", "DIAGNOSIS", "QCflagIR3", "GENDER"}, "Release metadata")
    require_columns(portal, {"PATNO", "HudAlphaSampleName", "Clinical Event", "Disease Status", "Diagnosis", "Case Control", "Long RNA-seq", "QC", "Sex"}, "Portal metadata")
    require_columns(age, {"PATNO", "EVENT_ID", "AGE_AT_VISIT"}, "Age metadata")
    pidx: dict[str, set[int]] = defaultdict(set)
    for index, row in enumerate(portal):
        name = (row.get("HudAlphaSampleName") or "").strip()
        if name:
            pidx[name].add(index)
    aidx: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in age:
        pid = patient(row.get("PATNO"))
        event = (row.get("EVENT_ID") or "").strip().upper()
        if pid and event:
            aidx[(pid, event)].append(row.get("AGE_AT_VISIT", ""))

    stages: dict[str, Counter[str]] = {group: Counter() for group in ("PD", "Control")}
    qc: dict[str, Counter[tuple[str, str]]] = {group: Counter() for group in stages}
    broad_qc: dict[str, Counter[str]] = {group: Counter() for group in stages}
    participants: dict[str, set[str]] = {group: set() for group in stages}
    eligible_participants: dict[str, set[str]] = {group: set() for group in stages}
    eligible: set[tuple[str, str]] = set()
    eligible_keys_seen: Counter[tuple[str, str]] = Counter()
    for row in release:
        if (row.get("CLINICAL_EVENT") or "").strip().upper() != "BL":
            continue
        diagnosis = clean(row.get("DIAGNOSIS"))
        if diagnosis not in {"pd", "control"}:
            continue
        group = "PD" if diagnosis == "pd" else "Control"
        count = stages[group]
        count["broad_release_rows"] += 1
        broad_flag = clean(row.get("QCflagIR3"))
        broad_qc[group][broad_flag if broad_flag in {"pass", "fail"} else "other_or_missing"] += 1
        pid = patient(row.get("PATNO"))
        if not pid or clean(pid) in {"null", "nan", "pool"}:
            count["missing_or_pool_id"] += 1
            continue
        names = {(row.get(k) or "").strip() for k in ("HudAlphaID", "Sample")}
        hits = set().union(*(pidx.get(name, set()) for name in names if name))
        if len(hits) != 1:
            count["no_portal_match" if not hits else "ambiguous_portal_match"] += 1
            continue
        prow = portal[next(iter(hits))]
        if patient(prow.get("PATNO")) != pid:
            count["portal_id_discordant"] += 1
            continue
        if (prow.get("Clinical Event") or "").strip().upper() != "BL":
            count["portal_not_baseline"] += 1
            continue
        if clean(prow.get("Long RNA-seq")) != "yes":
            count["portal_not_long_rna"] += 1
            continue
        if portal_group(prow) != group:
            count["historical_group_discordant_or_other"] += 1
            continue
        count["strict_linked_pre_qc_rows"] += 1
        participants[group].add(pid)
        rq = clean(row.get("QCflagIR3"))
        release_qc = rq if rq in {"pass", "fail"} else "other_or_missing"
        pq = clean(prow.get("QC"))
        portal_qc = "highqc" if pq == "highqc" else ("missing" if not pq else "other")
        qc[group][(release_qc, portal_qc)] += 1
        if release_qc != "pass" or portal_qc != "highqc":
            continue
        count["both_qc_pass_rows"] += 1
        if not age_valid(aidx.get((pid, "BL"), [])):
            count["post_qc_age_invalid_or_discordant"] += 1
            continue
        sex = clean(row.get("GENDER"))
        if sex not in {"male", "female"} or clean(prow.get("Sex")) != sex:
            count["post_qc_sex_invalid_or_discordant"] += 1
            continue
        count["post_qc_metadata_eligible_rows"] += 1
        eligible_participants[group].add(pid)
        key = (pid, (row.get("Sample") or "").strip())
        if not key[1]:
            count["post_qc_missing_release_sample"] += 1
        else:
            eligible.add(key)
            eligible_keys_seen[key] += 1

    result = {"status": "metadata_only_pre_second_lock", "unit": "release_sample_row",
              "definitions": "strict linked baseline Long RNA denominator before either QC flag; historical labels agree across release and portal"}
    result["groups"] = {}
    for group in ("PD", "Control"):
        denom = stages[group]["strict_linked_pre_qc_rows"]
        q = qc[group]
        release_fail = sum(n for (r, _), n in q.items() if r == "fail")
        portal_not_highqc = sum(n for (_, p), n in q.items() if p == "other")
        any_exclusion = sum(n for (r, p), n in q.items() if r != "pass" or p != "highqc")
        unknown = sum(n for (r, p), n in q.items() if r == "other_or_missing" or p == "missing")
        result["groups"][group] = {
            "stages": dict(sorted(stages[group].items())),
            "strict_denominator_unique_participants": len(participants[group]),
            "post_qc_metadata_eligible_unique_participants": len(eligible_participants[group]),
            "broad_release_qc_counts": dict(sorted(broad_qc[group].items())),
            "qc_cross_tab": [
                {"release_qc": r, "portal_qc": p, "rows": n}
                for (r, p), n in sorted(q.items())
            ],
            "rates": {
                "explicit_release_fail": {"numerator": release_fail, "denominator": denom, "fraction": release_fail / denom if denom else None},
                "portal_other_not_highqc": {"numerator": portal_not_highqc, "denominator": denom, "fraction": portal_not_highqc / denom if denom else None},
                "any_qc_exclusion": {"numerator": any_exclusion, "denominator": denom, "fraction": any_exclusion / denom if denom else None},
                "unknown_qc_state": {"numerator": unknown, "denominator": denom, "fraction": unknown / denom if denom else None},
            },
        }
    repeated = sum(n - 1 for n in eligible_keys_seen.values() if n > 1)
    missing_candidates = len(candidates - eligible)
    result["reconciliation"] = {
        "candidate_bundle_unique_keys": len(candidates),
        "candidate_keys_missing_from_metadata_eligible_pool": missing_candidates,
        "metadata_eligible_keys_not_in_candidate_bundle": len(eligible - candidates),
        "duplicate_metadata_eligible_keys": repeated,
        "candidate_subset_check_passed": missing_candidates == 0 and repeated == 0,
    }
    if sum(x["stages"].get("strict_linked_pre_qc_rows", 0) for x in result["groups"].values()) == 0:
        raise ValueError("No strict linked pre-QC rows; check input release and portal versions")
    return result


def self_test() -> None:
    def r(pid: str, sample: str, group: str, qc: str) -> dict[str, str]:
        return {"PATNO": pid, "Sample": sample, "HudAlphaID": sample, "CLINICAL_EVENT": "BL", "DIAGNOSIS": group, "QCflagIR3": qc, "GENDER": "Male"}
    def p(pid: str, sample: str, group: str, qc: str) -> dict[str, str]:
        return {"PATNO": pid, "HudAlphaSampleName": sample, "Clinical Event": "BL", "Disease Status": "Idiopathic PD" if group == "PD" else "Healthy Control", "Diagnosis": group, "Case Control": "Case" if group == "PD" else "Control", "Long RNA-seq": "Yes", "QC": qc, "Sex": "Male"}
    release = [r("1", "a", "PD", "Pass"), r("2", "b", "PD", "Fail"), r("3", "c", "Control", "Pass")]
    portal = [p("1", "a", "PD", "HighQC"), p("2", "b", "PD", "LowQC"), p("3", "c", "Control", "HighQC")]
    ages = [{"PATNO": str(i), "EVENT_ID": "BL", "AGE_AT_VISIT": "60"} for i in (1, 2, 3)]
    out = audit(release, portal, ages, {("1", "a"), ("3", "c")})
    assert out["groups"]["PD"]["stages"]["strict_linked_pre_qc_rows"] == 2
    assert out["groups"]["PD"]["rates"]["explicit_release_fail"]["fraction"] == 0.5
    assert out["reconciliation"]["candidate_subset_check_passed"]
    portal.append(p("1", "a", "PD", "HighQC"))
    out2 = audit(release, portal, ages, {("3", "c")})
    assert out2["groups"]["PD"]["stages"]["ambiguous_portal_match"] == 1
    print("SYNTHETIC_SELF_TEST_PASSED")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--release", type=Path)
    parser.add_argument("--portal", type=Path)
    parser.add_argument("--age", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if any(getattr(args, key) is None for key in ("release", "portal", "age", "bundle", "output")):
        parser.error("--release, --portal, --age, --bundle and --output are required")
    start = time.monotonic()
    try:
        print("Stage 1/4: reading complete metadata tables", flush=True)
        release, portal, age = (read_table(path) for path in (args.release, args.portal, args.age))
        print("Stage 2/4: reading candidate metadata only", flush=True)
        candidates = candidate_keys(args.bundle)
        print("Stage 3/4: counting pre-QC and post-QC groups", flush=True)
        result = audit(release, portal, age, candidates)
        result["input_sha256"] = {key: sha256(getattr(args, key)) for key in ("release", "portal", "age", "bundle")}
        result["input_row_counts"] = {"release": len(release), "portal": len(portal), "age": len(age)}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Stage 4/4: aggregate report written; elapsed {time.monotonic() - start:.1f}s", flush=True)
        print("PPMI_PREQC_QC_AUDIT_COMPLETE", flush=True)
        return 0
    except Exception as error:
        # Deliberately never print row values, identifiers, or traceback.
        print(f"PPMI_PREQC_QC_AUDIT_FAILED: {type(error).__name__}; no private row details printed", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
