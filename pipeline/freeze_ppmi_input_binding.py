#!/usr/bin/env python3
"""Freeze aggregate PPMI input and population hashes without computing scores.

This program reads participant-level rows locally, but emits only counts, hashes,
dimensions, and integrity flags. Never add row-level data to stdout or reports.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np

from audit_ppmi_candidate_bundle import audit as audit_candidates


EXPECTED_BUNDLE_SHA256 = "1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c"
EXPECTED_RAW_MATRIX_SHA256 = "00eda6e7a3df0f5f3486cf8f3ed42c4ee2ad958bd93cfddcacf689e612c7f611"
EXPECTED_GENE_REFERENCE_SHA256 = "ee90d4bb7e44ec507018b06ad8bf35a832291975f2f3b778288cbfc7940e15eb"
N_ROWS = 472
N_GENES = 16364


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def geo_array_hash(array: np.ndarray) -> str:
    """The exact shape/dtype/bytes convention in geo_pipeline.array_hash."""
    value = np.ascontiguousarray(array)
    h = hashlib.sha256()
    h.update(str(value.shape).encode())
    h.update(value.dtype.str.encode())
    h.update(value.tobytes())
    return h.hexdigest()


def ordered_population_hash(rows: list[dict]) -> str:
    records = []
    for position, row in enumerate(rows):
        records.append([
            position,
            row["release_row_index"],
            str(row["PATNO"]),
            str(row["release_metadata"]["Sample"]),
            str(row["historical_label"]),
        ])
    canonical = json.dumps(records, ensure_ascii=False, separators=(",", ":"))
    return sha256_bytes(canonical.encode("utf-8"))


def freeze(bundle: Path, preprocessing: Path) -> dict:
    bundle_sha = sha256_file(bundle)
    if bundle_sha != EXPECTED_BUNDLE_SHA256:
        raise ValueError("Bundle SHA-256 differs from the audited Mac copy")
    candidate_audit = audit_candidates(bundle)
    if (candidate_audit["candidate_count"] != N_ROWS
            or candidate_audit["unique_participants"] != N_ROWS
            or candidate_audit["unique_release_samples"] != N_ROWS
            or candidate_audit["unique_release_rows"] != N_ROWS
            or candidate_audit["group_counts"] != {"PD": 311, "Control": 161}
            or not candidate_audit["internal_file_hashes_match"]):
        raise ValueError("Candidate metadata integrity checks failed")

    with zipfile.ZipFile(bundle) as archive:
        manifest = json.loads(archive.read("MANIFEST.json"))
        raw = archive.read("rank_features.f32")
        metadata_bytes = archive.read("candidate_metadata.jsonl")
        gene_bytes = archive.read("common_genes_gencode_v29.tsv")

    raw_sha = sha256_bytes(raw)
    if raw_sha != EXPECTED_RAW_MATRIX_SHA256 or raw_sha != manifest["files_sha256"]["rank_features.f32"]:
        raise ValueError("Raw input-matrix SHA-256 mismatch")
    if len(raw) != N_ROWS * N_GENES * 4 or manifest["feature_shape"] != [N_ROWS, N_GENES]:
        raise ValueError("Raw input-matrix dimensions mismatch")
    if sha256_bytes(gene_bytes) != EXPECTED_GENE_REFERENCE_SHA256:
        raise ValueError("Bundled gene-reference SHA-256 mismatch")
    rows = [json.loads(line) for line in metadata_bytes.decode("utf-8").splitlines()]
    if len(rows) != N_ROWS:
        raise ValueError("Candidate metadata and matrix rows differ")
    population_sha = ordered_population_hash(rows)

    gene_reader = csv.DictReader(io.StringIO(gene_bytes.decode("utf-8")), delimiter="\t")
    gene_ids = np.array([row["gene_id"] for row in gene_reader])
    if len(gene_ids) != N_GENES or len(set(gene_ids.tolist())) != N_GENES:
        raise ValueError("Gene reference length or uniqueness failed")
    with np.load(preprocessing, allow_pickle=False) as prep:
        if not np.array_equal(gene_ids, prep["genes"]):
            raise ValueError("PPMI gene order disagrees with frozen GEO preprocessing")

    raw_array = np.frombuffer(raw, dtype="<f4").reshape(N_ROWS, N_GENES)
    if not np.isfinite(raw_array).all():
        raise ValueError("Nonfinite raw PPMI matrix value")
    lattice = raw_array.astype(np.float64) * (2 * N_GENES)
    lattice_error = float(np.max(np.abs(lattice - np.rint(lattice))))
    if lattice_error >= 0.01:
        raise ValueError("PPMI ranks are not on the expected GEO half-integer lattice")
    scored = np.rint(lattice) / (2 * N_GENES)
    if (not np.isfinite(scored).all() or scored.min() <= 0 or scored.max() >= 1
            or np.max(np.abs(scored.mean(axis=1) - .5)) > 1e-12):
        raise ValueError("Recovered PPMI rank matrix violates fixed invariants")

    result = {
        "status": "population_and_matrix_binding_passed_no_scores",
        "population": "final_computational_baseline_population",
        "candidate_count": N_ROWS,
        "group_counts": candidate_audit["group_counts"],
        "unique_participants": N_ROWS,
        "unique_release_samples": N_ROWS,
        "bundle_sha256": bundle_sha,
        "raw_matrix_member": "rank_features.f32",
        "raw_matrix_sha256": raw_sha,
        "raw_matrix_bytes": len(raw),
        "raw_matrix_shape": [N_ROWS, N_GENES],
        "raw_matrix_dtype": "<f4",
        "scored_matrix_transform": "float32 -> float64 -> round(rank*32728)/32728",
        "scored_matrix_shape": [N_ROWS, N_GENES],
        "scored_matrix_dtype": "<f8",
        "scored_matrix_geo_array_hash_sha256": geo_array_hash(scored),
        "ordered_population_sha256": population_sha,
        "population_hash_schema": "canonical compact UTF-8 JSON array of [matrix_row,release_row_index,exact_PATNO,exact_release_Sample,historical_label]",
        "candidate_metadata_sha256": sha256_bytes(metadata_bytes),
        "gene_reference_sha256": sha256_bytes(gene_bytes),
        "geo_preprocessing_sha256": sha256_file(preprocessing),
        "gene_order_matches_geo": True,
        "finite_and_rank_invariants_pass": True,
        "metadata_integrity_pass": True,
        "scores_computed": False,
        "models_fitted": False,
        "limitation": "Computational sample linkage is not independent biological donor-identity proof; group-specific pre-QC rates remain descriptive/unavailable without complete release metadata.",
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--geo-preprocessing", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Aggregate output already exists")
    result = freeze(args.bundle, args.geo_preprocessing)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
