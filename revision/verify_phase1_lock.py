"""Read-only phase-1 checksum audit; emits hashes, never participant records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path.home() / "PPMI_private"
MANIFEST = ROOT / "docs/PPMI_SECOND_LOCK_MANIFEST_PUBLISHED.json"
ZENODO_COPY = ROOT / "revision/outputs/zenodo_manifest_remote_P1.json"


def digest_stream(handle) -> str:
    h = hashlib.sha256()
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        h.update(block)
    return h.hexdigest()


def digest_file(path: Path) -> str:
    with path.open("rb") as handle:
        return digest_stream(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", choices=("start", "end", "expanded_start", "expanded_end"), required=True)
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    expected = manifest["sha256"]
    pooled = ROOT / "private_data/geo_runs_20261018/pooled_final"
    bundle = PRIVATE / "PPMI_WINDOWS_PRIVATE_BUNDLE.zip"
    rows = []

    def add(name: str, path: Path, reference: str | None = None, member: str | None = None) -> None:
        actual = ""
        status = "MISSING"
        if path.is_file():
            if member:
                with ZipFile(path) as archive, archive.open(member) as handle:
                    actual = digest_stream(handle)
            else:
                actual = digest_file(path)
            status = "PASS" if reference is None or actual == reference else "FAIL"
        rows.append((name, str(path), member or "", reference or "", actual, status))

    add("published_manifest_local", MANIFEST, digest_file(ZENODO_COPY) if ZENODO_COPY.is_file() else None)
    if ZENODO_COPY.is_file():
        add("published_manifest_Zenodo_download", ZENODO_COPY, digest_file(MANIFEST))
    for model in ("M1", "M2", "M3", "M4", "M5"):
        add(f"{model}_model", pooled / f"{model}_model.npz", expected[f"{model}_model"])
    add("geo_preprocessing", pooled / "preprocessing.npz", expected["geo_preprocessing"])
    add("sex_artifact", ROOT / "private_data/geo_final_sex_model_20261020/sex_model.npz", expected["sex_artifact"])
    add("private_bundle", bundle, expected["bundle"])
    add("raw_matrix", bundle, manifest["raw_matrix_sha256"], "rank_features.f32")
    add("candidate_metadata", bundle, manifest["candidate_metadata_sha256"], "candidate_metadata.jsonl")
    add("gene_reference", bundle, manifest["gene_reference_sha256"], "common_genes_gencode_v29.tsv")
    code_map = {
        "candidate_audit_code": "pipeline/audit_ppmi_candidate_bundle.py",
        "frozen_geo_scorer_code": "pipeline/frozen_geo_scorer.py",
        "geo_decisions": "private_data/geo_runs_20261018/pooled_final/decisions.json",
        "paired_bootstrap_code": "pipeline/paired_bootstrap.py",
        "population_audit_code": "pipeline/freeze_ppmi_input_binding.py",
        "ppmi_executor_code": "pipeline/score_ppmi_after_lock.py",
        "marker_manifest": "protocol/v3_marker_manifest.json",
    }
    for key, relative in code_map.items():
        add(key, ROOT / relative, expected[key])
    for relative, reference in manifest["supporting_sha256"].items():
        add(f"supporting:{relative}", ROOT / relative, reference)
    score_summary = PRIVATE / "score_20261021/AGGREGATE_RESULTS.json"
    if score_summary.is_file():
        reported = json.loads(score_summary.read_text())
        for key, locked_key in (("ordered_population_sha256", "ordered_population_sha256"),
                                ("raw_matrix_sha256", "raw_matrix_sha256"),
                                ("scored_matrix_geo_array_hash_sha256", "scored_matrix_geo_array_hash_sha256")):
            value = reported.get(key, "")
            rows.append((f"prior_score_reported:{key}", str(score_summary), "",
                         manifest[locked_key], value,
                         "REPORTED_MATCH" if value == manifest[locked_key] else "FAIL"))
    add("PPMI_predictions_current", PRIVATE / "score_20261021/PRIVATE_ROW_ORDERED_PREDICTIONS.npz")
    add("PPMI_sex_predictions_current", PRIVATE / "score_20261021/PRIVATE_ROW_ORDERED_SEX_CONTROL.npz")
    out = ROOT / "revision/outputs" / f"lock_verification_P1_{args.snapshot}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("artifact", "source_path", "zip_member", "expected_sha256", "actual_sha256", "status"))
        writer.writerows(rows)
    print(f"{args.snapshot}: {sum(row[-1] == 'PASS' for row in rows)} PASS, "
          f"{sum(row[-1] == 'REPORTED_MATCH' for row in rows)} reported matches, "
          f"{sum(row[-1] == 'FAIL' for row in rows)} FAIL, "
          f"{sum(row[-1] == 'MISSING' for row in rows)} MISSING")


if __name__ == "__main__":
    main()
