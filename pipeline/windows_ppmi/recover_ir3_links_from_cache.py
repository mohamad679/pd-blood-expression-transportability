"""Recover exact sample-prefix links from the already-created private IR3 cache.

Runs locally in seconds. No 42 GB archive access. Outputs aggregate counts only.
Keep IR3_PRIVATE_CACHE.zip on the approved personal computer.
"""
import argparse
import csv
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from audit_ir3_all_once import from_bytes, release_summary, norm_id


def release_index(meta_blob):
    text = meta_blob.decode("utf-8-sig")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    headers = next(reader)
    target_cols = [(i, h) for i, h in enumerate(headers)
                   if h in {"HudAlphaID", "Sample", "Specimen Bar Code"}]
    idx = defaultdict(set)
    for row_i, row in enumerate(reader):
        for col_i, _ in target_cols:
            if col_i < len(row):
                key = row[col_i].strip()
                if key:
                    idx[key].add(row_i)
    return idx, {h: sum(name == h for _, name in target_cols) for h in set(h for _, h in target_cols)}


def read_cached_csv(z, name):
    return list(csv.DictReader(io.StringIO(z.read("inputs/"+name).decode("utf-8-sig"))))


def prefix_key(value):
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    return base.split(".", 1)[0].strip()


def candidate_keys(value):
    """Full basename plus every complete dot-delimited prefix; no fuzzy IDs."""
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    keys = {base}
    keys.update(base[:i] for i, c in enumerate(base) if c == ".")
    return {key for key in keys if key}


def composite_patno_event(value):
    """PPMI guide p41: PPMI-Phase2-IR1.3174.V08...."""
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    found = re.search(r"PPMI-Phase[12]-IR[123]\.(\d+)\.(BL|V\d{2})\.", base, re.I)
    if not found:
        return None
    return norm_id(found.group(1)), found.group(2).upper()


def sample_candidates(value, exact_index, visit_index):
    exact = set().union(*(exact_index.get(key, set()) for key in candidate_keys(value)))
    visit = composite_patno_event(value)
    by_visit = visit_index.get(visit, set()) if visit else set()
    if exact and by_visit and exact != by_visit:
        return set(), "identifier_vs_visit_conflict"
    if exact:
        return exact, "exact_dot_prefix"
    if by_visit:
        return by_visit, "guide_composite_PATNO_visit"
    return set(), "unmatched"


def run(cache_path):
    with zipfile.ZipFile(cache_path) as z:
        index = json.loads(z.read("PRIVATE_index.json"))
        meta_blob = z.read("inputs/metaDataIR3_original")
        release = from_bytes(meta_blob)
        release_idx, sources = release_index(meta_blob)
        visit_idx = defaultdict(set)
        for i, row in enumerate(release):
            pid = norm_id(row.get("PATNO"))
            event = (row.get("CLINICAL_EVENT") or "").strip().upper()
            if pid and event:
                visit_idx[(pid, event)].add(i)
        portal = read_cached_csv(z, "meta_data.11192021.csv")
        ages = read_cached_csv(z, "Age_at_visit_13Oct2026.csv")
        statuses = read_cached_csv(z, "Participant_Status_13Oct2026.csv")
        release_audit, candidates = release_summary(release, portal, ages, statuses)
        missing_vectors = sum(1 for rec in index if rec.get("cache_entry") not in z.NameToInfo)
    counts = Counter()
    matches_by_row = Counter()
    for rec in index:
        header_matches, header_method = sample_candidates(rec["count_header"], release_idx, visit_idx)
        file_matches, file_method = sample_candidates(rec["member_name"], release_idx, visit_idx)
        if header_method == "identifier_vs_visit_conflict" or file_method == "identifier_vs_visit_conflict":
            counts["identifier_vs_visit_conflict"] += 1
            continue
        if header_matches and file_matches and header_matches != file_matches:
            counts["header_vs_filename_conflict"] += 1
            continue
        matches = header_matches or file_matches
        if len(matches) == 1:
            row_i = next(iter(matches))
            matches_by_row[row_i] += 1
            counts["uniquely_matched_count_files"] += 1
            if header_matches and file_matches:
                counts["two_sources_agree"] += 1
            elif header_matches:
                counts["header_only_" + header_method] += 1
            else:
                counts["filename_only_" + file_method] += 1
        elif len(matches) > 1:
            counts["ambiguous_metadata_match"] += 1
        else:
            counts["no_documented_format_match"] += 1
    linked = Counter()
    linked_patnos = defaultdict(set)
    for row_i, (label, pid) in candidates.items():
        n = matches_by_row[row_i]
        if n == 1:
            linked[label] += 1
            linked_patnos[label].add(pid)
        elif n > 1:
            linked["primary_candidate_multiple_count_tables"] += 1
        else:
            linked["primary_candidate_no_count_table"] += 1
    return {
        "audit_version": "cache_documented_ID_recovery_2.0",
        "scope": "already-created local private cache only; no raw archive access",
        "privacy": "aggregate counts only; no identifiers, filenames, count values or row data",
        "cache_entries": len(index), "cache_vectors_missing": missing_vectors,
        "release_rows": len(release), "candidate_rows_before_expression_link": len(candidates),
        "release_identifier_columns_indexed": sources,
        "matching_rule": "Exact dot-delimited prefixes against release identifiers; fallback to PPMI guide p41 composite PATNO+visit format. Header and filename links must agree whenever both resolve.",
        "link_results": dict(sorted(counts.items())),
        "release_rows_with_exactly_one_count_table": sum(v == 1 for v in matches_by_row.values()),
        "release_rows_with_multiple_count_tables": sum(v > 1 for v in matches_by_row.values()),
        "primary_candidate_counts_after_expression_link": dict(sorted(linked.items())),
        "primary_candidate_distinct_PATNO_after_expression_link": {k: len(v) for k, v in sorted(linked_patnos.items())},
        "notes": [
            "Release and portal metadata linked by original sample name, PATNO, baseline visit, sex and historical labels inside the local cache.",
            "A failed prefix match leaves eligibility unverified; never equate a filename pattern with participant identity without this exact metadata link.",
            "No classifier, AUC or external performance calculation was done."
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--output", type=Path, default=Path("IR3_LINK_RECOVERY_AUDIT.json"))
    args = parser.parse_args()
    try:
        result = run(args.cache)
        args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    except (OSError, ValueError, KeyError, UnicodeError, zipfile.BadZipFile, csv.Error):
        print("Cache audit failed. Check that the complete private cache and original script are beside this script. No private values printed.")
        return 1
    print("Finished. Send only IR3_LINK_RECOVERY_AUDIT.json. Keep the private cache local.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
