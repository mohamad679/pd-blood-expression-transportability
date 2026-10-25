"""Metadata-only, frozen-fold batch prevalence diagnostic for GENEPARK."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from geo_pipeline import Data, split_manifest


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/geo_batch_only_diagnostic_20261019"
EXPECTED_SPLIT_SHA256 = "2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2"


def main() -> None:
    audit = json.loads((ROOT / "data/manifests/input_audit.json").read_text())
    records: list[dict] = []
    source_ids: dict[str, set[str]] = {}
    for source in ("GSE57475", "GSE99039"):
        source_records = json.loads((ROOT / f"prepared/{source}_sample_metadata.json").read_text())
        source_ids[source] = {r["gsm"] for r in source_records}
        records.extend(source_records)
    if any(not row.get("gsm") or row.get("case") not in (0, 1) or not row.get("batch") for row in records):
        raise ValueError("Missing sample ID, binary case label, or recorded batch")
    if len({r["gsm"] for r in records}) != len(records):
        raise ValueError("Duplicate sample IDs in metadata")

    # Construct only the metadata carrier required by the frozen split function.
    # X is a zero-column placeholder: no expression file is opened or loaded.
    data = Data(
        X=np.empty((len(records), 0), dtype=float),
        age=np.asarray([r["age"] for r in records], dtype=float),
        sex=np.asarray([r["sex"] for r in records], dtype=float),
        y=np.asarray([r["case"] for r in records], dtype=int),
        source=np.asarray(["GSE57475" if r["gsm"] in source_ids["GSE57475"] else "GSE99039" for r in records]),
        donor=np.asarray([
            ("GSE57475" if r["gsm"] in source_ids["GSE57475"] else "GSE99039") + ":" + str(r["subject_id"] or r["gsm"])
            for r in records
        ], dtype=str),
        sample=np.asarray([r["gsm"] for r in records], dtype=str),
        genes=np.asarray([], dtype=str),
        marker_indices=[],
    )
    manifest = split_manifest(data, audit["input_sha256"])
    payload = json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False).encode()
    split_hash = hashlib.sha256(payload).hexdigest()
    if split_hash != EXPECTED_SPLIT_SHA256:
        raise ValueError(f"Frozen split hash mismatch: {split_hash}")

    by_id = {r["gsm"]: r for r in records}
    # Primary fixed-fold M1-M5 comparison comes only from the existing GEO summary.
    original = json.loads((ROOT / "results/geo_development_20261018/aggregate_summary.json").read_text())
    folds = []
    unseen_total = 0
    for item in manifest["splits"]:
        if not item["scope"].startswith("GSE99039/outer/"):
            continue
        fold = int(item["scope"].rsplit("/", 1)[1])
        train = [by_id[sid] for sid in item["train"]]
        test = [by_id[sid] for sid in item["test"]]
        train_prevalence = float(np.mean([r["case"] for r in train]))
        batch_rates: dict[str, float] = {}
        for batch in sorted({r["batch"] for r in train}):
            members = [r for r in train if r["batch"] == batch]
            batch_rates[batch] = sum(r["case"] for r in members) / len(members)
        scores = []
        unseen = 0
        for row in test:
            rate = batch_rates.get(row["batch"])
            if rate is None:
                rate = train_prevalence
                unseen += 1
            scores.append(rate)
        y_test = [r["case"] for r in test]
        auc = float(roc_auc_score(y_test, scores))
        unseen_total += unseen
        orig_metrics = original["outer_folds"]["GSE99039"]["fold_results"][item["scope"]]
        folds.append({
            "scope": item["scope"], "n_train": len(train), "n_test": len(test),
            "train_pd_prevalence": train_prevalence, "batch_only_AUROC": auc,
            "unseen_test_batch_samples_fallback_to_train_prevalence": unseen,
            **{f"{model}_AUROC_original": orig_metrics[model]["AUROC"] for model in ("M1", "M2", "M3", "M4", "M5")},
        })
    if len(folds) != 5:
        raise ValueError(f"Expected exactly five fixed GENEPARK outer folds, got {len(folds)}")

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "fold_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(folds[0]))
        writer.writeheader()
        writer.writerows(folds)
    comparison_means = {
        model: float(np.mean([row[f"{model}_AUROC_original"] for row in folds]))
        for model in ("M1", "M2", "M3", "M4", "M5")
    }
    summary = {
        "status": "complete", "analysis": "metadata-only batch-specific training PD prevalence diagnostic",
        "source": "GSE99039", "fold_count": 5, "split_manifest_sha256": split_hash,
        "mean_fold_AUROC": float(np.mean([row["batch_only_AUROC"] for row in folds])),
        "same_fold_mean_AUROC_original_models": comparison_means,
        "fold_results": folds,
        "unseen_test_batch_samples_total": unseen_total,
        "expression_loaded": False, "PPMI_loaded": False, "test_labels_used_for_scores": False,
        "definition": "Within each frozen outer fold, assign each test sample its batch's PD prevalence in outer training; unseen batches fall back to overall outer-training PD prevalence; no smoothing.",
        "interpretation": "Predictive batch-label association only; not causal attribution or evidence that batch solely explains model performance.",
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    report = [
        "# GENEPARK batch-only diagnostic",
        "",
        "Computed from sample metadata only on the five frozen GSE99039 outer folds. No expression or PPMI files were loaded.",
        "",
        f"Mean fold AUROC, batch-only: **{summary['mean_fold_AUROC']:.3f}**.",
        "",
        "| Fold | Batch-only AUROC | M2 AUROC | M5 AUROC | Unseen-batch fallback samples |",
        "|---|---:|---:|---:|---:|",
    ]
    report.extend(
        f"| {row['scope'].rsplit('/', 1)[1]} | {row['batch_only_AUROC']:.3f} | {row['M2_AUROC_original']:.3f} | {row['M5_AUROC_original']:.3f} | {row['unseen_test_batch_samples_fallback_to_train_prevalence']} |"
        for row in folds
    )
    report.extend(["", "| Same-fold mean AUROC | Value |", "|---|---:|"])
    report.extend(f"| {model} | {value:.3f} |" for model, value in comparison_means.items())
    report.extend([
        "",
        "Interpretation: batch prevalence is predictively informative in this fixed-fold diagnostic and its mean is close to the original M5 mean. It does not establish that batch causally explains model performance or represents disease biology. One test sample had an unseen training batch and received the outer-training prevalence fallback.",
        "",
        f"Frozen split SHA-256 verified: `{split_hash}`.",
    ])
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
