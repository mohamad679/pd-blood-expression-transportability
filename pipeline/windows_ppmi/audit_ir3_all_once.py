"""Single-pass local IR3 audit and private count cache. Python 3.10+ stdlib only.

Place this script, common_genes_gencode_v29.tsv, and three PPMI portal CSVs
beside each other. Send ONLY IR3_COMPLETE_AUDIT.json to the analyst.
Keep IR3_PRIVATE_CACHE.zip local: it contains individual data and counts.
"""
import argparse
import array
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import tarfile
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

PORTAL = "meta_data.11192021.csv"
AGE = "Age_at_visit_13Oct2026.csv"
STATUS = "Participant_Status_13Oct2026.csv"
GENES = "common_genes_gencode_v29.tsv"
META_NAMES = {"metadatair3.tsv", "metadatair3.csv", "metadatair3.txt"}
KNOWN_SUFFIXES = (
    ".Aligned.sortedByCoord.out.bam", ".Aligned.sortedByCoord.out", ".bam",
    ".featureCounts.txt", ".featurecounts.txt", ".counts.txt",
    ".counts", ".count", ".txt", ".tsv", ".csv", ".gz",
)


def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def from_bytes(blob):
    text = blob.decode("utf-8-sig")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    return list(csv.DictReader(io.StringIO(text), delimiter=delimiter))


def norm_id(s):
    s = (s or "").strip()
    if s.upper().startswith("PP-"):
        s = s[3:]
    return s


def cleaned(s):
    return " ".join((s or "").strip().lower().replace("_", " ").split())


def valid_age(s):
    try:
        age = float(s)
    except (ValueError, TypeError):
        return False
    return math.isfinite(age) and 0 < age < 120


def normalized_gene(s):
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", s)


def format_shape(s):
    # Never emit input tokens. Keep only separator structure.
    base = s.rsplit("/", 1)[-1]
    mask = "".join("T" if c.isalnum() else (c if c in "._-/" else "?") for c in base)
    return re.sub(r"T+", "TOKEN", mask)[:120]


def key_variants(value):
    base = value.replace("\\", "/").rsplit("/", 1)[-1]
    variants = {value, base}
    current = base
    for _ in range(5):
        suffix = next((x for x in KNOWN_SUFFIXES if current.lower().endswith(x.lower())), None)
        if suffix is None:
            break
        current = current[:-len(suffix)]
        variants.add(current)
    return {x for x in variants if x}


def index_release(rows):
    index = defaultdict(set)
    for i, row in enumerate(rows):
        for col in ("HudAlphaID", "Sample"):
            value = (row.get(col) or "").strip()
            if value:
                index[value].add(i)
    return index


def match_release(count_header, member_name, release_index):
    matches = set()
    for value in (count_header, member_name):
        for component in [value, *value.replace("\\", "/").split("/")]:
            for key in key_variants(component):
                matches.update(release_index.get(key, ()))
    return matches


def make_portal_index(rows):
    index = defaultdict(set)
    for i, row in enumerate(rows):
        name = (row.get("HudAlphaSampleName") or "").strip()
        if name:
            index[name].add(i)
    return index


def classify_portal(row):
    status = cleaned(row.get("Disease Status"))
    diagnosis = cleaned(row.get("Diagnosis"))
    case = cleaned(row.get("Case Control"))
    if (status, diagnosis, case) == ("idiopathic pd", "pd", "case"):
        return "idiopathic_PD_three_labels_agree"
    if (status, diagnosis, case) == ("healthy control", "control", "control"):
        return "healthy_control_three_labels_agree"
    return "other_or_discordant"


def age_index(rows):
    out = defaultdict(list)
    for row in rows:
        key = (norm_id(row.get("PATNO")), (row.get("EVENT_ID") or "").strip().upper())
        if key[0] and key[1]:
            out[key].append(row.get("AGE_AT_VISIT", ""))
    return out


def status_index(rows):
    out = defaultdict(list)
    for row in rows:
        pid = norm_id(row.get("PATNO"))
        if pid:
            out[pid].append(row)
    return out


def release_summary(release, portal, ages, statuses):
    pidx = make_portal_index(portal)
    aidx = age_index(ages)
    sidx = status_index(statuses)
    counters = Counter()
    candidate_patnos = {"idiopathic_PD_three_labels_agree": set(), "healthy_control_three_labels_agree": set()}
    eligible_release_rows = {}
    eligible_cov = defaultdict(Counter)
    plate_labels = defaultdict(set)
    for release_i, row in enumerate(release):
        qc = cleaned(row.get("QCflagIR3"))
        event = (row.get("CLINICAL_EVENT") or "").strip().upper()
        diagnosis = cleaned(row.get("DIAGNOSIS"))
        diagnosis = diagnosis if diagnosis in {"pd", "control", "prodromal", "swedd"} else "other_or_unrecognised"
        counters[("release_QC", qc if qc in {"pass", "fail"} else "other")] += 1
        counters[("release_event", event if event in {"BL", "V02", "V04", "V06", "V08"} else "other")] += 1
        counters[("release_diagnosis", diagnosis)] += 1
        if qc != "pass" or event != "BL":
            continue
        counters[("BL_Pass_diagnosis", diagnosis)] += 1
        pid = norm_id(row.get("PATNO"))
        if not pid or pid.lower() in {"null", "nan", "pool"}:
            counters[("candidate", "missing_or_pool_PATNO")] += 1
            continue
        names = {(row.get(c) or "").strip() for c in ("HudAlphaID", "Sample")}
        hit = set().union(*(pidx.get(name, set()) for name in names if name))
        if not hit:
            counters[("candidate", "no_portal_sample_match")] += 1
            continue
        if len(hit) != 1:
            counters[("candidate", "ambiguous_portal_sample_match")] += 1
            continue
        prow = portal[next(iter(hit))]
        if norm_id(prow.get("PATNO")) != pid:
            counters[("candidate", "PATNO_discordant")] += 1
            continue
        if (prow.get("Clinical Event") or "").strip().upper() != "BL":
            counters[("candidate", "portal_visit_not_BL")] += 1
            continue
        if cleaned(prow.get("Long RNA-seq")) != "yes" or cleaned(prow.get("QC")) != "highqc":
            counters[("candidate", "portal_not_LongRNA_HighQC")] += 1
            continue
        if cleaned(row.get("GENDER")) not in {"male", "female"}:
            counters[("candidate", "release_sex_unrecognised")] += 1
            continue
        psex = cleaned(prow.get("Sex"))
        if psex not in {"male", "female"}:
            counters[("candidate", "portal_sex_unrecognised")] += 1
            continue
        if psex in {"male", "female"} and psex != cleaned(row.get("GENDER")):
            counters[("candidate", "sex_discordant")] += 1
            continue
        label = classify_portal(prow)
        if label == "other_or_discordant":
            counters[("candidate", "portal_labels_not_primary")] += 1
            continue
        if (label.startswith("idiopathic") and diagnosis != "pd") or (label.startswith("healthy") and diagnosis != "control"):
            counters[("candidate", "release_diagnosis_discordant")] += 1
            continue
        avalues = aidx.get((pid, "BL"), [])
        valid_values = {float(x) for x in avalues if valid_age(x)}
        if len(valid_values) != 1:
            counters[("candidate", "BL_age_missing_invalid_or_discordant")] += 1
            continue
        counters[("candidate", label)] += 1
        candidate_patnos[label].add(pid)
        eligible_release_rows[release_i] = (label, pid)
        eligible_cov[label]["female" if psex == "female" else "male"] += 1
        eligible_cov[label]["age_0_to_49" if next(iter(valid_values)) < 50 else
                            ("age_50_to_69" if next(iter(valid_values)) < 70 else "age_70_plus")] += 1
        for field in ("Neutrophil Score", "Neutrophils (%)", "Lymphocytes (%)",
                      "Basophils (%)", "Eosinophils (%)", "RIN Value", "Plate", "Box"):
            if (prow.get(field) or "").strip():
                eligible_cov[label]["nonempty_" + field] += 1
        plate = (row.get("Plate") or "").strip()
        if plate:
            plate_labels[plate].add(label)
        counters[("candidate_current_status", "found" if pid in sidx else "not_found")] += 1
    counts = [{"dimension": dim, "category": cat, "rows": n} for (dim, cat), n in sorted(counters.items())]
    return {"release_rows": len(release), "portal_rows": len(portal), "age_rows": len(ages),
            "participant_status_rows": len(statuses), "aggregate_counts": counts,
            "distinct_PATNO_in_primary_candidates": {k: len(v) for k, v in candidate_patnos.items()},
            "primary_candidate_covariates": {k: dict(sorted(v.items())) for k, v in sorted(eligible_cov.items())},
            "primary_candidate_plates": {"distinct": len(plate_labels),
                "both_case_and_control": sum(len(v) > 1 for v in plate_labels.values())},
            "note": "Candidate label uses historical portal agreement, not current 2026 cohort. No future status exclusion. Final count requires linked expression and donor review."}, eligible_release_rows


def count_file(handle, member_name, reference, release_index):
    # A tar stream is not seekable, so decode one line at a time.
    reader = csv.reader((line.decode("utf-8-sig") for line in handle), delimiter="\t")
    header = next((row for row in reader if row and not row[0].startswith("#")), None)
    if not header or [h.lower() for h in header[:6]] != ["geneid", "chr", "start", "end", "strand", "length"] or len(header) != 7:
        return None
    matches = match_release(header[6], member_name, release_index)
    values = array.array("d", [math.nan]) * len(reference)
    gene_index = {g: i for i, g in enumerate(reference)}
    stats = Counter()
    seen = set()
    lengths_match = lengths_checked = 0
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        stats["gene_rows"] += 1
        if len(row) != 7:
            stats["wrong_width"] += 1
            continue
        gene = normalized_gene(row[0])
        if gene in seen:
            stats["duplicate_gene_ids"] += 1
        seen.add(gene)
        try:
            count = float(row[6]); length = float(row[5])
        except ValueError:
            stats["non_numeric_count_or_length"] += 1
            continue
        if not math.isfinite(count) or count < 0:
            stats["invalid_count"] += 1
            continue
        if not math.isfinite(length) or length <= 0:
            stats["invalid_length"] += 1
            continue
        if not count.is_integer():
            stats["fractional_count"] += 1
        pos = gene_index.get(gene)
        if pos is not None:
            values[pos] = count
            lengths_checked += 1
            lengths_match += int(length == reference[gene])
    present = sum(math.isfinite(x) for x in values)
    info = {"stats": dict(stats), "common_genes_present": present,
            "common_genes_missing": len(reference)-present,
            "common_length_equal_to_reference": lengths_match,
            "common_length_checked": lengths_checked,
            "release_metadata_match_count": len(matches),
            "header_format": format_shape(header[6]),
            "filename_format": format_shape(member_name.replace("\\", "/"))}
    return info, values, (next(iter(matches)) if len(matches) == 1 else None), header[6]


def audit(archive, directory, report_path, cache_path):
    t0 = time.monotonic()
    with (directory / GENES).open("r", encoding="utf-8-sig", newline="") as ref_handle:
        ref_rows = list(csv.DictReader(ref_handle, delimiter="\t"))
    reference = {r["gene_id"]: int(r["exon_union_length_bp"]) for r in ref_rows}
    if len(reference) != 16364 or len(ref_rows) != len(reference) or min(reference.values()) <= 0:
        raise ValueError("Invalid reference")
    portal = read_csv(directory / PORTAL)
    ages = read_csv(directory / AGE)
    statuses = read_csv(directory / STATUS)
    report = {"audit_version": "all_once_1.0", "archive_bytes": archive.stat().st_size,
              "scope": "complete tar.gz stream, all featureCounts tables, release and three portal CSVs",
              "privacy": "Aggregate output only; no participant IDs, sample IDs, original paths or expression values.",
              "archive_end_reached": False, "archive_members_seen": 0,
              "count_files_seen": 0, "quant_sf_files_seen": 0,
              "count_tables_verified": 0, "metadata_found": False}
    counts = Counter()
    rows_with_unique_link = set()
    matches_by_release = Counter()
    sizes = []
    coverage = []
    length_ratios = []
    first_formats = []
    release = None
    release_index = None
    eligible_release_rows = {}
    private_index = []
    with zipfile.ZipFile(cache_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=3, allowZip64=True) as cache:
        cache.write(directory / GENES, "inputs/" + GENES)
        for filename in (PORTAL, AGE, STATUS):
            cache.write(directory / filename, "inputs/" + filename)
        with tarfile.open(archive, "r|gz") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                report["archive_members_seen"] += 1
                name = member.name.replace("\\", "/")
                base = name.rsplit("/", 1)[-1].lower()
                if base in META_NAMES:
                    if release is not None:
                        counts["multiple_release_metadata_files"] += 1
                        continue
                    if member.size > 64 * 1024 * 1024:
                        raise ValueError("Release metadata over 64 MiB limit")
                    blob = tar.extractfile(member).read()
                    release = from_bytes(blob)
                    release_index = index_release(release)
                    report["metadata_found"] = True
                    report["release_audit"], eligible_release_rows = release_summary(release, portal, ages, statuses)
                    cache.writestr("inputs/metaDataIR3_original", blob)
                    continue
                if base.endswith(".sf"):
                    report["quant_sf_files_seen"] += 1
                is_counts = "featurecount" in base or any(x.lower() in {"counts", "featurecounts"} for x in name.split("/"))
                if not is_counts:
                    continue
                report["count_files_seen"] += 1
                if release_index is None:
                    counts["count_files_before_release_metadata"] += 1
                    continue
                handle = tar.extractfile(member)
                if base.endswith(".gz"):
                    import gzip
                    handle = gzip.GzipFile(fileobj=handle)
                result = count_file(handle, name, reference, release_index)
                if result is None:
                    counts["candidate_file_not_single_sample_featureCounts"] += 1
                    continue
                info, vector, rowindex, header = result
                i = report["count_tables_verified"]
                report["count_tables_verified"] += 1
                for k, v in info["stats"].items():
                    counts[k] += v
                counts["metadata_match_count_" + str(info["release_metadata_match_count"] if info["release_metadata_match_count"] < 3 else "3plus")] += 1
                coverage.append(info["common_genes_present"])
                sizes.append(info["stats"].get("gene_rows", 0))
                if info["common_length_checked"]:
                    length_ratios.append(info["common_length_equal_to_reference"] / info["common_length_checked"])
                if len(first_formats) < 8:
                    first_formats.append({"header_format": info["header_format"],
                                          "filename_format": info["filename_format"],
                                          "metadata_matches": info["release_metadata_match_count"]})
                if rowindex is not None:
                    rows_with_unique_link.add(rowindex)
                    matches_by_release[rowindex] += 1
                # Cache is private: it contains individual gene counts and original ID mapping.
                entry = f"counts/{i:06d}.f64"
                cache.writestr(entry, vector.tobytes())
                private_index.append({"cache_entry": entry, "release_row_index": rowindex,
                                      "member_name": name, "count_header": header})
                if report["count_tables_verified"] % 250 == 0:
                    print("Checked", report["count_tables_verified"], "count tables...", flush=True)
            else:
                report["archive_end_reached"] = True
        cache.writestr("PRIVATE_index.json", json.dumps(private_index))
        cache.writestr("PRIVATE_cache_note.txt", "CONTAINS INDIVIDUAL PPMI DATA. Keep on the approved local machine. Never send as an aggregate report.\n")
    report["numeric_and_structure_totals"] = dict(sorted(counts.items()))
    report["count_table_gene_rows"] = {"min": min(sizes) if sizes else None, "max": max(sizes) if sizes else None}
    report["common_genes_per_count_table"] = {"reference_genes": len(reference),
        "min": min(coverage) if coverage else None, "max": max(coverage) if coverage else None,
        "tables_missing_any_reference_gene": sum(n < len(reference) for n in coverage)}
    report["length_agreement_fraction_per_table"] = {"min": min(length_ratios) if length_ratios else None,
        "max": max(length_ratios) if length_ratios else None}
    report["matched_release_rows"] = len(rows_with_unique_link)
    report["release_rows_with_multiple_count_tables"] = sum(n > 1 for n in matches_by_release.values())
    linked_primary = Counter()
    linked_patnos = defaultdict(set)
    for release_i, (label, pid) in eligible_release_rows.items():
        if matches_by_release[release_i] == 1:
            linked_primary[label] += 1
            linked_patnos[label].add(pid)
        elif matches_by_release[release_i] > 1:
            linked_primary["excluded_multiple_count_tables"] += 1
        else:
            linked_primary["excluded_no_unique_count_table"] += 1
    report["primary_candidate_rows_with_one_uniquely_linked_count_table"] = dict(sorted(linked_primary.items()))
    report["primary_candidate_distinct_PATNO_with_count_table"] = {k: len(v) for k, v in sorted(linked_patnos.items())}
    report["primary_candidates_with_multiple_rows_per_PATNO"] = {
        k: linked_primary[k] - len(v) for k, v in sorted(linked_patnos.items())}
    report["first_eight_redacted_formats"] = first_formats
    report["elapsed_seconds"] = round(time.monotonic()-t0, 1)
    report["private_cache_written"] = cache_path.exists()
    report["private_cache_bytes"] = cache_path.stat().st_size if cache_path.exists() else 0
    report["limitations"] = [
        "This is structural and metadata audit, not a classifier or external validation.",
        "Automatic sample matching is conservative; unmatched count files remain cached for local follow-up without re-reading the 42 GB archive.",
        "Historical portal labels and current status were not conflated. Release metadata has no age column; BL ages are joined locally by exact PATNO and event.",
        "Private cache contains participant-level information and must not be transferred as an aggregate result."
    ]
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--input-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "IR3_COMPLETE_AUDIT.json"
    cache_path = args.output_dir / "IR3_PRIVATE_CACHE.zip"
    temporary_cache = args.output_dir / "IR3_PRIVATE_CACHE.inprogress.zip"
    if cache_path.exists():
        print("Private cache already exists. Move or rename it before retrying to avoid overwriting individual data.")
        return 2
    try:
        report = audit(args.archive, args.input_dir, report_path, temporary_cache)
        temporary_cache.replace(cache_path)
    except (OSError, ValueError, KeyError, EOFError, tarfile.TarError, csv.Error, UnicodeError, zipfile.BadZipFile):
        print("Audit did not finish. Check archive completeness, free disk space and required CSVs. Do not send a partial private cache or report.")
        return 1
    print("Complete archive scanned:", report["archive_end_reached"])
    print("Count tables verified:", report["count_tables_verified"])
    print("Send ONLY IR3_COMPLETE_AUDIT.json. Keep IR3_PRIVATE_CACHE.zip on this Windows computer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
