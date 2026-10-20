#!/usr/bin/env python3
"""Prepare a private, score-free PPMI candidate bundle from the prior IR3 cache.

Python 3.10+ standard library only. Run on the Windows computer holding
IR3_PRIVATE_CACHE.zip. Never send the private output through chat or commit it.
"""

from __future__ import annotations

import argparse
import array
import csv
import hashlib
import io
import json
import math
import os
import struct
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from audit_ir3_all_once import (AGE, GENES, PORTAL, STATUS, age_index, valid_age,
                                classify_portal, from_bytes, make_portal_index,
                                norm_id, release_summary)
from recover_ir3_links_from_cache import release_index, sample_candidates

REFERENCE_SHA256 = "ee90d4bb7e44ec507018b06ad8bf35a832291975f2f3b778288cbfc7940e15eb"
AUDIT_SCRIPT_SHA256 = "754ce5b9c6b9b25aeb256f7288d329da7d511aeecdbbe1be55b1313106627aeb"
LINK_SCRIPT_SHA256 = "a6b190e87387744c82c43c9c70cda4722f78a88f49c5f6037f164e8c637e6263"
EXPECTED_COUNT_FILES = 4871
EXPECTED_CANDIDATES = 472
GENE_COUNT = 16364
REPORT_NAME = "PPMI_WINDOWS_PREP_AGGREGATE.json"
BUNDLE_NAME = "PPMI_WINDOWS_PRIVATE_BUNDLE.zip"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def progress(done: int, total: int, start: float, phase: str) -> None:
    fraction = done / total if total else 1.0
    width = 24
    bar = "#" * round(width * fraction) + "-" * (width - round(width * fraction))
    elapsed = time.monotonic() - start
    eta = elapsed * (total - done) / done if done else 0
    print(f"\r{phase} [{bar}] {done}/{total} {fraction:5.1%} | "
          f"elapsed {elapsed/60:.1f}m | ETA {eta/60:.1f}m", end="", flush=True)
    if done == total:
        print(flush=True)


def read_csv_from_zip(cache: zipfile.ZipFile, name: str) -> list[dict[str, str]]:
    with cache.open("inputs/" + name) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")))


def load_reference(blob: bytes) -> tuple[list[str], list[int]]:
    if sha256_bytes(blob) != REFERENCE_SHA256:
        raise ValueError("The cached gene reference hash differs from the frozen reference")
    rows = list(csv.DictReader(io.StringIO(blob.decode("utf-8-sig")), delimiter="\t"))
    genes = [r["gene_id"] for r in rows]
    lengths = [int(r["exon_union_length_bp"]) for r in rows]
    if len(genes) != GENE_COUNT or len(set(genes)) != GENE_COUNT or min(lengths) <= 0:
        raise ValueError("Gene reference has the wrong size, duplicates, or invalid lengths")
    return genes, lengths


def rank_count_vector(blob: bytes, lengths: list[int]) -> tuple[bytes, dict[str, float | int]]:
    """Use fixed exon-union lengths, then GEO-equivalent within-sample midranks."""
    n = len(lengths)
    if len(blob) != 8 * n:
        raise ValueError("Cached count vector has the wrong byte length")
    counts = array.array("d")
    counts.frombytes(blob)
    if sys.byteorder != "little":
        counts.byteswap()
    if any(not math.isfinite(x) or x < 0 for x in counts):
        raise ValueError("Cached count vector contains missing or invalid values")
    adjusted = [count / length for count, length in zip(counts, lengths)]
    order = sorted(range(n), key=adjusted.__getitem__)
    ranks = array.array("f", [0.0]) * n
    lo = 0
    while lo < n:
        hi = lo + 1
        value = adjusted[order[lo]]
        while hi < n and adjusted[order[hi]] == value:
            hi += 1
        # 0-based equal-value positions [lo, hi): midrank=(lo+1+hi)/2.
        fraction = (lo + hi) / (2 * n)
        for at in range(lo, hi):
            ranks[order[at]] = fraction
        lo = hi
    if abs(math.fsum(ranks) / n - 0.5) > 2e-7:
        raise ValueError("Within-sample fractional-rank mean is not 0.5")
    if sys.byteorder != "little":
        ranks.byteswap()
    return ranks.tobytes(), {"zero_count_genes": sum(x == 0 for x in counts),
                              "mean_rank": math.fsum(ranks) / n}


def exact_links(index: list[dict], release: list[dict], meta_blob: bytes) -> tuple[dict[int, dict], dict[str, int]]:
    exact_index, _ = release_index(meta_blob)
    visit_index: dict[tuple[str, str], set[int]] = defaultdict(set)
    for i, row in enumerate(release):
        key = (norm_id(row.get("PATNO")), (row.get("CLINICAL_EVENT") or "").strip().upper())
        if key[0] and key[1]:
            visit_index[key].add(i)
    by_release: dict[int, dict] = {}
    methods = Counter()
    for rec in index:
        a, amethod = sample_candidates(rec["count_header"], exact_index, visit_index)
        b, bmethod = sample_candidates(rec["member_name"], exact_index, visit_index)
        if amethod == "identifier_vs_visit_conflict" or bmethod == "identifier_vs_visit_conflict":
            raise ValueError("Count-file identifier conflicts with participant/visit metadata")
        if a and b and a != b:
            raise ValueError("Count header and filename point to different release rows")
        hits = a or b
        if len(hits) != 1:
            raise ValueError("A cached count file has no unique exact metadata link")
        release_i = next(iter(hits))
        if release_i in by_release:
            raise ValueError("Multiple count files point to one release row")
        by_release[release_i] = rec
        methods["header_and_filename" if a and b else "header_only" if a else "filename_only"] += 1
    if len(by_release) != EXPECTED_COUNT_FILES or len(index) != EXPECTED_COUNT_FILES:
        raise ValueError("Cached count-file count differs from the audited 4,871")
    return by_release, dict(methods)


def candidate_metadata(release_i: int, label: str, pid: str, release: list[dict],
                       portal: list[dict], portal_index: dict, ages: dict) -> dict:
    r = release[release_i]
    names = {(r.get(c) or "").strip() for c in ("HudAlphaID", "Sample")}
    hits = set().union(*(portal_index.get(x, set()) for x in names if x))
    if len(hits) != 1:
        raise ValueError("Candidate portal metadata does not have one exact release link")
    p = portal[next(iter(hits))]
    age_values = {float(x) for x in ages.get((pid, "BL"), []) if valid_age(x)}
    if len(age_values) != 1:
        raise ValueError("Candidate has no single valid baseline age")
    if classify_portal(p) != label:
        raise ValueError("Candidate historical label changed")
    rplate = (r.get("Plate") or "").strip()
    pplate = (p.get("Plate") or "").strip()
    if not rplate:
        raise ValueError("Release RNA sequencing plate is missing")
    return {"release_row_index": release_i, "PATNO": pid, "EVENT_ID": "BL",
            "historical_label": "PD" if label.startswith("idiopathic") else "Control",
            "age_at_BL_visit": next(iter(age_values)),
            "sex": (r.get("GENDER") or "").strip(),
            "plate": rplate, "plate_source": "metaDataIR3_original",
            "portal_plate": pplate,
            "release_metadata": r, "portal_metadata": p}


def run(cache_path: Path, output_dir: Path) -> tuple[Path, Path]:
    if not cache_path.is_file():
        raise ValueError("IR3_PRIVATE_CACHE.zip was not found")
    for filename, expected in (("audit_ir3_all_once.py", AUDIT_SCRIPT_SHA256),
                               ("recover_ir3_links_from_cache.py", LINK_SCRIPT_SHA256)):
        if sha256_file(Path(__file__).with_name(filename)) != expected:
            raise ValueError("Bundled audit dependency hash differs from the audited version")
    output_dir.mkdir(parents=True, exist_ok=True)
    bundle = output_dir / BUNDLE_NAME
    report = output_dir / REPORT_NAME
    if bundle.exists() or report.exists():
        raise ValueError("Output already exists; choose a new empty output directory")
    work = output_dir / "PRIVATE_INPROGRESS"
    if work.exists():
        raise ValueError("PRIVATE_INPROGRESS already exists; inspect the previous failed run")
    work.mkdir()
    started = time.monotonic()
    with zipfile.ZipFile(cache_path) as cache:
        bad = cache.testzip()
        if bad:
            raise ValueError("Private cache ZIP integrity check failed")
        gene_blob = cache.read("inputs/" + GENES)
        genes, lengths = load_reference(gene_blob)
        meta_blob = cache.read("inputs/metaDataIR3_original")
        release = from_bytes(meta_blob)
        portal = read_csv_from_zip(cache, PORTAL)
        age_rows = read_csv_from_zip(cache, AGE)
        status_rows = read_csv_from_zip(cache, STATUS)
        _, candidates = release_summary(release, portal, age_rows, status_rows)
        if len(release) != EXPECTED_COUNT_FILES or len(candidates) != EXPECTED_CANDIDATES:
            raise ValueError("Release or candidate count differs from the previous aggregate audit")
        if len({pid for _, pid in candidates.values()}) != EXPECTED_CANDIDATES:
            raise ValueError("Candidate participant IDs are not one-to-one")
        index = json.loads(cache.read("PRIVATE_index.json"))
        linked, link_methods = exact_links(index, release, meta_blob)
        if any(rec["cache_entry"] not in cache.NameToInfo for rec in linked.values()):
            raise ValueError("At least one cached count vector is absent")
        portal_index = make_portal_index(portal)
        ages = age_index(age_rows)
        selected = sorted(candidates)
        plate_groups = defaultdict(Counter)
        missing_plates = Counter()
        portal_plate_comparison = Counter()
        zero_counts = []
        rank_path = work / "rank_features.f32"
        metadata_path = work / "candidate_metadata.jsonl"
        last_update = 0.0
        with rank_path.open("wb") as ranks_out, metadata_path.open("w", encoding="utf-8", newline="\n") as metadata_out:
            for done, release_i in enumerate(selected, 1):
                label, pid = candidates[release_i]
                item = candidate_metadata(release_i, label, pid, release, portal, portal_index, ages)
                blob = cache.read(linked[release_i]["cache_entry"])
                rank_blob, stats = rank_count_vector(blob, lengths)
                ranks_out.write(rank_blob)
                metadata_out.write(json.dumps(item, ensure_ascii=False, sort_keys=True) + "\n")
                zero_counts.append(stats["zero_count_genes"])
                if item["plate"]:
                    plate_groups[item["plate"]][item["historical_label"]] += 1
                else:
                    missing_plates[item["historical_label"]] += 1
                if not item["portal_plate"]:
                    portal_plate_comparison["portal_plate_missing"] += 1
                elif item["portal_plate"] == item["plate"]:
                    portal_plate_comparison["same_literal_value"] += 1
                else:
                    portal_plate_comparison["different_literal_value"] += 1
                now = time.monotonic()
                if done == 1 or done == len(selected) or now - last_update >= 10:
                    progress(done, len(selected), started, "Prepare")
                    last_update = now
        if rank_path.stat().st_size != EXPECTED_CANDIDATES * GENE_COUNT * 4:
            raise ValueError("Rank matrix byte count is wrong")
        ref_path = work / GENES
        ref_path.write_bytes(gene_blob)
        manifest = {
            "schema": "PPMI_Windows_candidate_bundle_v1", "status": "provisional_prelock_preparation_only",
            "candidate_count": len(selected), "gene_count": len(genes),
            "gene_order": "Exact row order of common_genes_gencode_v29.tsv, matching frozen GEO inputs",
            "feature_transformation": "featureCounts / fixed exon-union length, then within-sample fractional midrank (midrank-0.5)/16364",
            "feature_dtype": "little-endian float32", "feature_shape": [len(selected), len(genes)],
            "feature_row_order": "candidate_metadata.jsonl line order; sorted release metadata row index",
            "scores_computed": False, "models_fitted": False,
            "selection": "Provisional baseline Long RNA/QC/age/sex/historical-label candidates from prior aggregate audit; final eligibility pending",
            "source_cache_sha256": sha256_file(cache_path),
            "source_modules_sha256": {"audit_ir3_all_once.py": AUDIT_SCRIPT_SHA256,
                                      "recover_ir3_links_from_cache.py": LINK_SCRIPT_SHA256},
            "preparation_script_sha256": sha256_file(Path(__file__)),
            "files_sha256": {p.name: sha256_file(p) for p in (rank_path, metadata_path, ref_path)},
        }
        manifest_path = work / "MANIFEST.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        tmp_bundle = output_dir / (BUNDLE_NAME + ".inprogress")
        with zipfile.ZipFile(tmp_bundle, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=3, allowZip64=True) as out:
            for p in (rank_path, metadata_path, ref_path, manifest_path):
                out.write(p, p.name)
        with zipfile.ZipFile(tmp_bundle) as check:
            if check.testzip() is not None:
                raise ValueError("Output bundle ZIP verification failed")
        tmp_bundle.replace(bundle)
        aggregate = {
            "status": "complete", "purpose": "metadata and feature preparation before second public lock; no scores",
            "candidate_count": len(selected), "PD": sum(x.startswith("idiopathic") for x, _ in candidates.values()),
            "Control": sum(x.startswith("healthy") for x, _ in candidates.values()),
            "plate_by_group": [{"plate": plate, "PD": group["PD"], "Control": group["Control"]}
                               for plate, group in sorted(plate_groups.items())],
            "plate_source": "Release RNA sequencing metadata (metaDataIR3_original); portal Plate is retained separately and is not assumed to be the same field",
            "release_vs_portal_plate_literal_comparison": dict(portal_plate_comparison),
            "missing_plate_by_group": dict(missing_plates),
            "link_methods_all_count_files": link_methods,
            "zero_count_genes_per_candidate": {"min": min(zero_counts), "max": max(zero_counts)},
            "gene_reference_sha256": REFERENCE_SHA256,
            "source_cache_sha256": manifest["source_cache_sha256"],
            "bundle_sha256": sha256_file(bundle), "bundle_bytes": bundle.stat().st_size,
            "elapsed_minutes": round((time.monotonic() - started) / 60, 2),
            "privacy": "Aggregate only. The private bundle contains participant-level data; never upload it to chat, GitHub, or Zenodo.",
        }
        report.write_text(json.dumps(aggregate, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for p in (rank_path, metadata_path, ref_path, manifest_path):
        p.unlink()
    work.rmdir()
    return bundle, report


def self_test() -> None:
    b = struct.pack("<4d", 0, 1, 1, 2)
    ranks, info = rank_count_vector(b, [1, 1, 1, 1])
    got = struct.unpack("<4f", ranks)
    assert got == (0.125, 0.5, 0.5, 0.875), got
    assert info["zero_count_genes"] == 1
    b = struct.pack("<4d", 1, 2, 3, 4)
    ranks, _ = rank_count_vector(b, [1, 2, 3, 4])
    assert struct.unpack("<4f", ranks) == (0.5,) * 4
    print("SELF_TEST_PASSED")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, help="Existing IR3_PRIVATE_CACHE.zip on Windows")
    parser.add_argument("--output-dir", type=Path, help="New empty private output directory")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.cache is None or args.output_dir is None:
        parser.error("--cache and --output-dir are required")
    try:
        bundle, report = run(args.cache, args.output_dir)
    except (OSError, ValueError, KeyError, UnicodeError, json.JSONDecodeError, zipfile.BadZipFile, csv.Error) as exc:
        print(f"\nPREPARATION_FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        print("Keep the original cache. Do not transfer partial output.", file=sys.stderr)
        return 1
    print(f"PRIVATE_BUNDLE_READY: {bundle}")
    print(f"SHARE_ONLY_AGGREGATE: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
