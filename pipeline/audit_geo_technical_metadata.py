"""Metadata-only technical batch/chip audit; no expression values, fits, or scores."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import geo_pipeline as geo

INPUT_AUDIT = ROOT / "data/manifests/input_audit.json"
RELEASE = ROOT / "results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json"
OUT = ROOT / "results/geo_technical_metadata_audit_20261019"
SOURCES = ("GSE57475", "GSE99039")


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    audit = json.loads(INPUT_AUDIT.read_text(encoding="utf-8"))
    release = json.loads(RELEASE.read_text(encoding="utf-8"))
    metadata_by_source: dict[str, list[dict]] = {}
    all_rows: list[dict] = []
    # Read only prepared sample metadata. Do not load prepared expression arrays.
    for source in SOURCES:
        meta_path = ROOT / "prepared" / f"{source}_sample_metadata.json"
        expected = audit["input_sha256"][meta_path.name]
        observed = _sha(meta_path)
        if observed != expected:
            raise ValueError(f"Metadata hash mismatch for {source}")
        rows = json.loads(meta_path.read_text(encoding="utf-8"))
        if len({r["gsm"] for r in rows}) != len(rows):
            raise ValueError(f"Duplicate GSM in {source} metadata")
        for row in rows:
            row = dict(row)
            row["source"] = source
            metadata_by_source.setdefault(source, []).append(row)
            all_rows.append(row)

    # Reconstruct the published frozen splits from labels/donor IDs only. The
    # dummy X array is never inspected by split_manifest or donor_folds.
    meta_order = [r for source in SOURCES for r in metadata_by_source[source]]
    n = len(meta_order)
    samples = np.asarray([r["gsm"] for r in meta_order])
    labels = np.asarray([int(r["case"]) for r in meta_order], dtype=int)
    source_arr = np.asarray([r["source"] for r in meta_order])
    donors = np.asarray([r["source"] + ":" + str(r.get("subject_id") or r["gsm"])
                         for r in meta_order])
    # Metadata-only dummy arrays; X/age/sex are not used in split construction.
    data = geo.Data(np.zeros((n, 1)), np.zeros(n), np.zeros(n), labels,
                    source_arr, donors, samples, np.asarray(["metadata-only"]), [])
    splits = geo.split_manifest(data, audit["input_sha256"])
    split_sha = hashlib.sha256(json.dumps(
        splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()
    if split_sha != release["split_manifest_sha256"]:
        raise ValueError("Metadata-only reconstructed splits do not match frozen release")
    frozen = {entry["scope"]: entry for entry in splits["splits"]}

    fold_of_sample: dict[str, int] = {}
    for source in SOURCES:
        for fold in range(5):
            scope = f"{source}/outer/{fold}"
            for gsm in frozen[scope]["test"]:
                if gsm in fold_of_sample:
                    raise ValueError(f"Sample assigned to multiple outer test folds: {gsm}")
                fold_of_sample[gsm] = fold

    coverage_rows, class_batch_rows, fold_batch_rows = [], [], []
    for source in SOURCES:
        rows = metadata_by_source[source]
        batch_kind_counts = Counter(r.get("batch_kind", "") for r in rows)
        n_batch = sum(bool(str(r.get("batch", "")).strip()) for r in rows)
        coverage_rows.append({
            "source": source, "n_samples": len(rows),
            "technical_metadata_present": n_batch,
            "technical_metadata_missing": len(rows) - n_batch,
            "coverage_percent": round(100 * n_batch / len(rows), 4) if rows else None,
            "distinct_batch_or_chip_proxy_levels": len({str(r.get("batch", "")).strip()
                                                         for r in rows if str(r.get("batch", "")).strip()}),
            "batch_kind": "; ".join(f"{key or 'missing'}={value}"
                                     for key, value in sorted(batch_kind_counts.items())),
        })
        batches = sorted({str(r.get("batch", "")).strip() or "MISSING"
                          for r in rows})
        for batch in batches:
            group = [r for r in rows if (str(r.get("batch", "")).strip() or "MISSING") == batch]
            n_case = sum(int(r["case"]) == 1 for r in group)
            n_control = sum(int(r["case"]) == 0 for r in group)
            case_pct = 100 * n_case / len(group) if group else None
            class_batch_rows.append({
                "source": source, "batch_or_chip_proxy": batch, "n": len(group),
                "PD": n_case, "control": n_control,
                "PD_percent": round(case_pct, 2) if case_pct is not None else None,
                "class_purity_percent": round(max(n_case, n_control) * 100 / len(group), 2)
                if group else None,
                "single_class_batch": n_case == 0 or n_control == 0,
                "at_least_90_percent_one_class_descriptive_flag": bool(
                    case_pct is not None and (case_pct <= 10 or case_pct >= 90)),
            })
            for fold in range(5):
                fg = [r for r in group if fold_of_sample.get(r["gsm"]) == fold]
                fn = len(fg)
                fold_batch_rows.append({
                    "source": source, "batch_or_chip_proxy": batch,
                    "outer_test_fold": fold, "n": fn,
                    "PD": sum(int(r["case"]) == 1 for r in fg),
                    "control": sum(int(r["case"]) == 0 for r in fg),
                    "percent_of_batch_in_fold": round(100 * fn / len(group), 2)
                    if group else None,
                })

    OUT.mkdir(parents=True, exist_ok=False)
    _csv(OUT / "technical_metadata_coverage.csv", coverage_rows,
         list(coverage_rows[0]))
    _csv(OUT / "class_by_batch.csv", class_batch_rows, list(class_batch_rows[0]))
    _csv(OUT / "batch_by_outer_test_fold.csv", fold_batch_rows, list(fold_batch_rows[0]))
    result = {
        "status": "metadata_only_audit_complete",
        "input_metadata_sha256": {source: _sha(ROOT / "prepared" / f"{source}_sample_metadata.json")
                                   for source in SOURCES},
        "frozen_split_manifest_sha256": split_sha,
        "frozen_split_matches_published_release": True,
        "PPMI_loaded": False, "expression_values_loaded": False,
        "model_fits": False, "scores_or_performance_metrics_computed": False,
        "technical_metadata_used_in_M1_to_M5": False,
        "outer_fold_assignment": "Each sample is assigned to the unique frozen outer test fold in which it is held out.",
        "near_single_class_descriptive_threshold": "PD prevalence <=10% or >=90%; descriptive flag only, not a hypothesis test.",
        "coverage": coverage_rows,
        "batch_level_summary": class_batch_rows,
        "batch_by_outer_test_fold": fold_batch_rows,
        "outputs": ["technical_metadata_coverage.csv", "class_by_batch.csv",
                    "batch_by_outer_test_fold.csv"],
    }
    (OUT / "RESULTS.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                                        encoding="utf-8")
    lines = [
        "# GEO technical metadata audit (metadata only)", "",
        "No expression values were loaded; no model was fit and no score/performance metric was computed. Technical metadata was not used in M1–M5.", "",
        f"Frozen split manifest verified: `{split_sha}`.", "",
        "## Metadata coverage", "",
        "| Source | N | Metadata present | Missing | Coverage | Levels | Type |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in coverage_rows:
        lines.append(f"| {r['source']} | {r['n_samples']} | {r['technical_metadata_present']} | {r['technical_metadata_missing']} | {r['coverage_percent']:.1f}% | {r['distinct_batch_or_chip_proxy_levels']} | {r['batch_kind']} |")
    lines.extend(["", "## Class composition by batch/chip proxy", "",
                  "See `class_by_batch.csv` for every level. The 10%/90% flag is descriptive only.", "",
                  "| Source | Batch/chip | N | PD | Control | PD % | Purity % | Near single-class flag |",
                  "|---|---|---:|---:|---:|---:|---:|---|"])
    for r in class_batch_rows:
        lines.append(f"| {r['source']} | {r['batch_or_chip_proxy']} | {r['n']} | {r['PD']} | {r['control']} | {r['PD_percent']:.1f}% | {r['class_purity_percent']:.1f}% | {r['at_least_90_percent_one_class_descriptive_flag']} |")
    lines.extend(["", "## Batch distribution across frozen outer test folds", "",
                  "Each sample is counted once, in its assigned outer test fold. Full counts are in `batch_by_outer_test_fold.csv`.", ""])
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"status": result["status"], "coverage": coverage_rows,
                      "near_single_class_batches": [r for r in class_batch_rows
                          if r["at_least_90_percent_one_class_descriptive_flag"]],
                      "single_class_batches": [r for r in class_batch_rows if r["single_class_batch"]],
                      "frozen_split_manifest_sha256": split_sha,
                      "output": str(OUT)}, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
