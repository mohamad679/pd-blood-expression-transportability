#!/usr/bin/env python3
"""One-pass, score-free preparation of remaining PPMI inputs from the IR3 cache.

Python 3.10+ standard library only. The complete metadata and optional all-visit
rank vectors are private outputs. No model, cell-panel, or similarity score is
computed. stdout and the aggregate JSON contain counts and hashes only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import tempfile
import time
import zipfile
from collections import Counter
from pathlib import Path

from audit_ir3_all_once import AGE, GENES, PORTAL, STATUS, from_bytes, release_summary
from audit_ppmi_preqc_qc import (
    EXPECTED_BUNDLE_SHA256,
    audit as preqc_audit,
    candidate_keys as bundle_candidate_keys,
    patient,
    self_test as preqc_self_test,
)
from prepare_ppmi_windows_bundle import (
    exact_links,
    load_reference,
    rank_count_vector,
    self_test as rank_self_test,
)

RELEASE_MEMBER = "inputs/metaDataIR3_original"
INDEX_MEMBER = "PRIVATE_index.json"
SOURCE_MEMBERS = (
    RELEASE_MEMBER,
    "inputs/" + PORTAL,
    "inputs/" + AGE,
    "inputs/" + STATUS,
    INDEX_MEMBER,
    "inputs/" + GENES,
)
CAPSULE_NAME = "PPMI_COMPLETE_METADATA_PRIVATE.zip"
REPORT_NAME = "PPMI_COMPLETE_METADATA_AGGREGATE.json"
RANK_BUNDLE_NAME = "PPMI_ALL_VISIT_RANKS_PRIVATE.zip"
CAPSULE_SCHEMA = "PPMI_complete_metadata_capsule_v1"
EXPECTED_RELEASE_ROWS = 4871
EXPECTED_CANDIDATES = 472
EXPECTED_GROUPS = {"PD": 311, "Control": 161}


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def sha256_file(path: Path, show_progress: bool = False) -> str:
    digest = hashlib.sha256()
    size = path.stat().st_size
    done = 0
    started = time.monotonic()
    last_display = 0.0
    last_decile = -1
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
            done += len(chunk)
            if not show_progress:
                continue
            now = time.monotonic()
            decile = min(10, int(10 * done / size)) if size else 10
            should_show = done == size or (sys.stdout.isatty() and now - last_display >= 3) or (
                not sys.stdout.isatty() and decile > last_decile
            )
            if not should_show:
                continue
            fraction = done / size if size else 1.0
            width = 24
            filled = round(width * fraction)
            elapsed = now - started
            eta = elapsed * (size - done) / done if done else 0
            message = (f"Hashing cache [{'#' * filled}{'-' * (width - filled)}] "
                       f"{fraction:5.1%} | {done / 2**30:.2f}/{size / 2**30:.2f} GiB | "
                       f"elapsed {elapsed / 60:.1f}m | ETA {eta / 60:.1f}m")
            print(("\r" if sys.stdout.isatty() else "") + message,
                  end="" if sys.stdout.isatty() and done < size else "\n", flush=True)
            last_display = now
            last_decile = decile
    return digest.hexdigest()


def parse_csv(blob: bytes) -> list[dict[str, str]]:
    text = blob.decode("utf-8-sig")
    header = text.splitlines()[0] if text else ""
    delimiter = "\t" if header.count("\t") > header.count(",") else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("A source table has missing or duplicate column names")
    return list(reader)


def read_source_zip(path: Path, capsule: bool) -> tuple[dict[str, bytes], dict]:
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        if not set(SOURCE_MEMBERS) <= names:
            raise ValueError("A required cached metadata member is missing")
        blobs = {name: archive.read(name) for name in SOURCE_MEMBERS}
        if not capsule:
            return blobs, {}
        if "MANIFEST.json" not in names:
            raise ValueError("Metadata capsule manifest is missing")
        manifest = json.loads(archive.read("MANIFEST.json"))
        if manifest.get("schema") != CAPSULE_SCHEMA:
            raise ValueError("Metadata capsule schema differs")
        for name, blob in blobs.items():
            if sha256_bytes(blob) != manifest["members_sha256"].get(name):
                raise ValueError("A metadata capsule member hash differs")
        return blobs, manifest


def candidate_selection(
    release: list[dict[str, str]],
    portal: list[dict[str, str]],
    ages: list[dict[str, str]],
    statuses: list[dict[str, str]],
    index_blob: bytes,
    release_blob: bytes,
) -> tuple[set[tuple[str, str]], dict, dict[int, dict]]:
    _, eligible = release_summary(release, portal, ages, statuses)
    index = json.loads(index_blob)
    linked, methods = exact_links(index, release, release_blob)
    selected = sorted(set(eligible) & set(linked))
    keys: set[tuple[str, str]] = set()
    groups: Counter[str] = Counter()
    donors: set[str] = set()
    for release_i in selected:
        group, pid = eligible[release_i]
        if group not in {"idiopathic_PD_three_labels_agree", "healthy_control_three_labels_agree"}:
            raise ValueError("Selected candidate has an unexpected historical-group category")
        key = (patient(pid), (release[release_i].get("Sample") or "").strip())
        if not all(key) or key in keys:
            raise ValueError("Selected candidate identifiers are empty or repeated")
        keys.add(key)
        groups["PD" if group.startswith("idiopathic") else "Control"] += 1
        donors.add(key[0])
    if len(release) != EXPECTED_RELEASE_ROWS or len(selected) != EXPECTED_CANDIDATES:
        raise ValueError("Release or selected candidate count differs from the audited run")
    if len(donors) != EXPECTED_CANDIDATES:
        raise ValueError("Selected candidate participant linkage is not one-to-one")
    if dict(groups) != EXPECTED_GROUPS:
        raise ValueError("Selected candidate historical-group counts differ from the audited run")
    return keys, {
        "release_rows_with_exact_count_link": len(linked),
        "candidate_rows_with_exact_count_link": len(selected),
        "candidate_unique_participants": len(donors),
        "candidate_historical_groups": dict(sorted(groups.items())),
        "count_link_methods": methods,
    }, linked


def write_capsule(path: Path, blobs: dict[str, bytes], source_sha: str) -> None:
    manifest = {
        "schema": CAPSULE_SCHEMA,
        "purpose": "Complete source metadata and count-file linkage index; no expression vectors or scores",
        "source_cache_sha256": source_sha,
        "members_sha256": {name: sha256_bytes(blobs[name]) for name in SOURCE_MEMBERS},
        "members_bytes": {name: len(blobs[name]) for name in SOURCE_MEMBERS},
        "contains_participant_identifiers": True,
        "share_in_chat_or_public_repository": False,
    }
    pending = path.with_name(path.name + ".inprogress")
    with zipfile.ZipFile(pending, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name in SOURCE_MEMBERS:
            archive.writestr(name, blobs[name])
        archive.writestr("MANIFEST.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with zipfile.ZipFile(pending) as archive:
        if archive.testzip() is not None:
            raise ValueError("Metadata capsule ZIP integrity check failed")
    pending.replace(path)


def sample_progress(done: int, total: int, started: float, last_display: float,
                    last_decile: int, phase: str = "Preparing ranks") -> tuple[float, int]:
    now = time.monotonic()
    decile = min(10, int(10 * done / total))
    if done != total and not ((sys.stdout.isatty() and now - last_display >= 5) or
                              (not sys.stdout.isatty() and decile > last_decile)):
        return last_display, last_decile
    fraction = done / total
    filled = round(24 * fraction)
    elapsed = now - started
    eta = elapsed * (total - done) / done if done else 0
    amount = (f"{done / 2**20:.0f}/{total / 2**20:.0f} MiB" if phase == "Packaging ranks"
              else f"{done}/{total}")
    message = (f"{phase} [{'#' * filled}{'-' * (24 - filled)}] "
               f"{amount} {fraction:5.1%} | elapsed {elapsed / 60:.1f}m | ETA {eta / 60:.1f}m")
    print(("\r" if sys.stdout.isatty() else "") + message,
          end="" if sys.stdout.isatty() and done < total else "\n", flush=True)
    return now, decile


def write_all_visit_ranks(
    cache_path: Path,
    output_dir: Path,
    blobs: dict[str, bytes],
    release: list[dict[str, str]],
    linked: dict[int, dict],
    source_sha: str,
) -> tuple[Path, dict]:
    """Prepare all within-sample ranks, without centering or computing similarity."""
    genes, lengths = load_reference(blobs["inputs/" + GENES])
    if len(linked) != len(release):
        raise ValueError("Some release rows lack an exact cached count-file link")
    path = output_dir / RANK_BUNDLE_NAME
    pending = path.with_name(path.name + ".inprogress")
    rank_file = output_dir / "rank_features.f32.inprogress"
    metadata_file = output_dir / "all_visit_metadata.jsonl.inprogress"
    started = time.monotonic()
    last_display = 0.0
    last_decile = -1
    release_qc = Counter()
    rank_digest = hashlib.sha256()
    metadata_digest = hashlib.sha256()
    try:
        with zipfile.ZipFile(cache_path) as cache, rank_file.open("wb") as rank_out, \
                metadata_file.open("w", encoding="utf-8", newline="\n") as metadata_out:
            for i, row in enumerate(release):
                rank_blob, _ = rank_count_vector(cache.read(linked[i]["cache_entry"]), lengths)
                rank_out.write(rank_blob)
                rank_digest.update(rank_blob)
                record = {
                    "release_row_index": i,
                    "PATNO": patient(row.get("PATNO")),
                    "Sample": (row.get("Sample") or "").strip(),
                    "CLINICAL_EVENT": (row.get("CLINICAL_EVENT") or "").strip(),
                    "DIAGNOSIS": (row.get("DIAGNOSIS") or "").strip(),
                    "GENDER": (row.get("GENDER") or "").strip(),
                    "QCflagIR3": (row.get("QCflagIR3") or "").strip(),
                    "Plate": (row.get("Plate") or "").strip(),
                }
                metadata_line = json.dumps(record, sort_keys=True) + "\n"
                metadata_out.write(metadata_line)
                metadata_digest.update(metadata_line.encode("utf-8"))
                release_qc[record["QCflagIR3"].lower()] += 1
                last_display, last_decile = sample_progress(i + 1, len(release), started,
                                                             last_display, last_decile)
        expected_bytes = len(release) * len(genes) * 4
        if rank_file.stat().st_size != expected_bytes:
            raise ValueError("All-visit rank matrix has unexpected byte count")
        rank_manifest = {
            "schema": "PPMI_all_visit_rank_preparation_v1",
            "purpose": "Preparation for the registered post-lock longitudinal identity diagnostic",
            "source_cache_sha256": source_sha,
            "sample_rows": len(release),
            "gene_count": len(genes),
            "feature_shape": [len(release), len(genes)],
            "feature_dtype": "little-endian float32",
            "feature_order": "Release metadata row order, matching all_visit_metadata.jsonl",
            "gene_order": "Exact order in common_genes_gencode_v29.tsv",
            "transformation": "featureCounts divided by fixed exon-union length, then within-sample fractional midrank",
            "scores_computed": False,
            "models_fitted": False,
            "files_sha256": {
                "rank_features.f32": rank_digest.hexdigest(),
                "all_visit_metadata.jsonl": metadata_digest.hexdigest(),
                GENES: sha256_bytes(blobs["inputs/" + GENES]),
            },
        }
        print("Packaging private all-visit rank bundle", flush=True)
        with zipfile.ZipFile(pending, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=1, allowZip64=True) as archive:
            with rank_file.open("rb") as source_file, archive.open("rank_features.f32", "w", force_zip64=True) as sink:
                copied = 0
                last_pack_display = 0.0
                last_pack_decile = -1
                pack_started = time.monotonic()
                for block in iter(lambda: source_file.read(8 * 1024 * 1024), b""):
                    sink.write(block)
                    copied += len(block)
                    last_pack_display, last_pack_decile = sample_progress(
                        copied, expected_bytes, pack_started, last_pack_display, last_pack_decile,
                        phase="Packaging ranks"
                    )
            archive.write(metadata_file, "all_visit_metadata.jsonl")
            archive.writestr(GENES, blobs["inputs/" + GENES])
            archive.writestr("MANIFEST.json", json.dumps(rank_manifest, indent=2, sort_keys=True) + "\n")
        print("Checking private rank ZIP integrity", flush=True)
        with zipfile.ZipFile(pending) as archive:
            if archive.testzip() is not None:
                raise ValueError("All-visit rank bundle ZIP integrity check failed")
        pending.replace(path)
    finally:
        if rank_file.exists():
            rank_file.unlink()
        if metadata_file.exists():
            metadata_file.unlink()
    return path, {"sample_rows": len(release), "release_qc_counts": dict(sorted(release_qc.items())),
                  "rank_bundle_sha256": sha256_file(path), "rank_bundle_bytes": path.stat().st_size}


def run(args: argparse.Namespace) -> None:
    source = args.cache or args.capsule
    if not source.is_file():
        raise FileNotFoundError("Source cache or capsule was not found")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Output directory is not empty; choose a fresh directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    print("Stage 1/5: reading only metadata and linkage members", flush=True)
    blobs, manifest = read_source_zip(source, capsule=args.capsule is not None)
    release = from_bytes(blobs[RELEASE_MEMBER])
    portal = parse_csv(blobs["inputs/" + PORTAL])
    ages = parse_csv(blobs["inputs/" + AGE])
    statuses = parse_csv(blobs["inputs/" + STATUS])

    print("Stage 2/5: checking exact candidate linkage and historical groups", flush=True)
    keys, selection, linked = candidate_selection(
        release, portal, ages, statuses, blobs[INDEX_MEMBER], blobs[RELEASE_MEMBER]
    )
    frozen_bundle_comparison = "not_supplied"
    if args.candidate_bundle is not None:
        if sha256_file(args.candidate_bundle) != EXPECTED_BUNDLE_SHA256:
            raise ValueError("Candidate bundle hash differs from the frozen provisional bundle")
        if bundle_candidate_keys(args.candidate_bundle) != keys:
            raise ValueError("Frozen candidate bundle differs from the source metadata selection")
        frozen_bundle_comparison = "exact_match"

    print("Stage 3/5: counting complete pre-QC and post-QC groups", flush=True)
    qc = preqc_audit(release, portal, ages, keys)
    if not qc["reconciliation"]["candidate_subset_check_passed"]:
        raise ValueError("Candidate subset reconciliation failed")

    print("Stage 4/5: recording source hash and preparing private inputs", flush=True)
    source_sha = manifest["source_cache_sha256"] if args.capsule else sha256_file(source, show_progress=True)
    capsule_path = args.output_dir / CAPSULE_NAME
    if args.cache:
        write_capsule(capsule_path, blobs, source_sha)
        rank_path, rank_summary = write_all_visit_ranks(source, args.output_dir, blobs,
                                                        release, linked, source_sha)
    else:
        capsule_path = source
        rank_path, rank_summary = None, None

    print("Stage 5/5: writing aggregate report", flush=True)
    report = {
        "status": "complete_metadata_only_no_scores",
        "source_kind": "existing_ir3_private_cache" if args.cache else "metadata_capsule_verification",
        "source_cache_sha256": source_sha,
        "source_member_sha256": {name: sha256_bytes(blobs[name]) for name in SOURCE_MEMBERS},
        "source_table_rows": {
            "release": len(release), "portal": len(portal), "age": len(ages), "participant_status": len(statuses)
        },
        "selection": selection,
        "preqc": qc,
        "frozen_candidate_bundle_comparison": frozen_bundle_comparison,
        "metadata_capsule_sha256": sha256_file(capsule_path),
        "metadata_capsule_bytes": capsule_path.stat().st_size,
        "all_visit_rank_preparation": rank_summary,
        "elapsed_minutes": round((time.monotonic() - start) / 60, 2),
        "privacy": "Aggregate counts and hashes only. Never share the private metadata capsule in chat or public repositories.",
        "interpretation": "Metadata audit and all-visit feature preparation; provisional population until full eligibility and crosswalk review; no model fit, panel/model score, or longitudinal similarity computed.",
    }
    out = args.output_dir / REPORT_NAME
    pending = out.with_name(out.name + ".inprogress")
    pending.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    pending.replace(out)
    print("PPMI_COMPLETE_METADATA_AUDIT_DONE", flush=True)
    print("SHARE_AGGREGATE_JSON: " + str(out), flush=True)
    if args.cache:
        print("TRANSFER_PRIVATELY_TO_MAC_ONLY: " + str(capsule_path), flush=True)
        print("TRANSFER_PRIVATELY_TO_MAC_ONLY: " + str(rank_path), flush=True)


def self_test() -> None:
    preqc_self_test()
    rank_self_test()
    sample = b"PATNO,EVENT_ID\n1,BL\n"
    assert parse_csv(sample) == [{"PATNO": "1", "EVENT_ID": "BL"}]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / CAPSULE_NAME
        blobs = {name: sample for name in SOURCE_MEMBERS}
        write_capsule(path, blobs, "a" * 64)
        restored, manifest = read_source_zip(path, capsule=True)
        assert restored == blobs and manifest["schema"] == CAPSULE_SCHEMA
    print("ONE_SHOT_METADATA_SELF_TEST_PASSED", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--cache", type=Path, help="Existing IR3_PRIVATE_CACHE.zip on Windows")
    source.add_argument("--capsule", type=Path, help="Private metadata capsule for Mac verification")
    parser.add_argument("--candidate-bundle", type=Path, help="Optional frozen PPMI_WINDOWS_PRIVATE_BUNDLE.zip")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (args.cache is None and args.capsule is None) or args.output_dir is None:
        parser.error("Specify --cache or --capsule and a new --output-dir")
    try:
        run(args)
        return 0
    except Exception as error:
        # Never print exception text or a traceback: upstream libraries might
        # include a participant/sample value in either one.
        print(f"PPMI_COMPLETE_METADATA_AUDIT_FAILED: {type(error).__name__}", file=sys.stderr)
        print("No participant rows or identifiers were printed. Keep partial output private.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
