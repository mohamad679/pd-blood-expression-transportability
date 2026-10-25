"""Within-batch AUROC diagnostic using saved GENEPARK outer-test predictions only."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from geo_pipeline import Data, split_manifest


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/geo_within_batch_auc_20261020"
EXPECTED_SPLIT_SHA256 = "2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2"
MINIMUM_PER_CLASS = 5


def show_progress(current: int, total: int, label: str, width: int = 24) -> None:
    """Show a compact progress bar, keeping redirected output readable."""
    filled = round(width * current / total)
    bar = "#" * filled + "-" * (width - filled)
    line = f"Progress [{bar}] {current}/{total} ({100 * current / total:5.1f}%) | {label}"
    if sys.stdout.isatty():
        sys.stdout.write("\r" + line + ("\n" if current == total else ""))
        sys.stdout.flush()
    else:
        print(line, flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _frozen_manifest(records: list[dict], input_hashes: dict) -> tuple[dict, str]:
    by_source: dict[str, list[dict]] = defaultdict(list)
    for row in records:
        by_source[row["source"]].append(row)
    samples = [row["gsm"] for row in records]
    data = Data(
        X=np.empty((len(records), 0), dtype=float),
        age=np.asarray([row["age"] for row in records], dtype=float),
        sex=np.asarray([row["sex"] for row in records], dtype=float),
        y=np.asarray([row["case"] for row in records], dtype=int),
        source=np.asarray([row["source"] for row in records], dtype=str),
        donor=np.asarray([
            row["source"] + ":" + str(row["subject_id"] or row["gsm"])
            for row in records
        ], dtype=str),
        sample=np.asarray(samples, dtype=str),
        genes=np.asarray([], dtype=str),
        marker_indices=[],
    )
    manifest = split_manifest(data, input_hashes)
    payload = json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False).encode()
    return manifest, hashlib.sha256(payload).hexdigest()


def main() -> None:
    audit_path = ROOT / "data/manifests/input_audit.json"
    audit = json.loads(audit_path.read_text())
    records: list[dict] = []
    input_hashes: dict[str, str] = {}
    for source in ("GSE57475", "GSE99039"):
        metadata_path = ROOT / f"prepared/{source}_sample_metadata.json"
        rows = json.loads(metadata_path.read_text())
        if any(not row.get("gsm") or row.get("case") not in (0, 1)
               or not row.get("batch") for row in rows):
            raise ValueError(f"Missing sample ID, binary diagnosis, or batch in {source}")
        if len({row["gsm"] for row in rows}) != len(rows):
            raise ValueError(f"Duplicate sample ID in {source} metadata")
        records.extend({**row, "source": source} for row in rows)
        input_hashes[metadata_path.name] = sha256(metadata_path)

    if len({row["gsm"] for row in records}) != len(records):
        raise ValueError("Duplicate sample ID across GEO metadata")
    manifest, split_hash = _frozen_manifest(records, audit["input_sha256"])
    if split_hash != EXPECTED_SPLIT_SHA256:
        raise ValueError(f"Frozen split hash mismatch: {split_hash}")

    rows_by_id = {row["gsm"]: row for row in records}
    collected: dict[str, dict] = {}
    seen: Counter[str] = Counter()
    output_rows: list[dict] = []
    folds = [s for s in manifest["splits"] if s["scope"].startswith("GSE99039/outer/")]
    if len(folds) != 5:
        raise ValueError(f"Expected 5 fixed GENEPARK outer folds, found {len(folds)}")

    for fold_index, split in enumerate(sorted(folds, key=lambda x: x["scope"]), start=1):
        fold = int(split["scope"].rsplit("/", 1)[1])
        prediction_path = ROOT / f"private_data/geo_runs_20261018/GSE99039_outer_{fold}/heldout_predictions.npz"
        if not prediction_path.is_file():
            raise FileNotFoundError(f"Missing saved held-out predictions for outer/{fold}")
        with np.load(prediction_path, allow_pickle=False) as prediction:
            required = {"samples", "y", "M4", "M5"}
            if not required.issubset(prediction.files):
                raise ValueError(f"Missing prediction arrays in outer/{fold}: {sorted(required - set(prediction.files))}")
            samples = prediction["samples"].astype(str).tolist()
            labels = prediction["y"].astype(int)
            expected_test = split["test"]
            if samples != expected_test:
                raise ValueError(f"Saved sample order does not match frozen outer/{fold} test split")
            if len(set(samples)) != len(samples):
                raise ValueError(f"Duplicate held-out sample in outer/{fold}")
            for i, sample_id in enumerate(samples):
                row = rows_by_id.get(sample_id)
                if row is None or row["source"] != "GSE99039":
                    raise ValueError(f"Unknown or wrong-source sample in saved outer/{fold} predictions")
                if int(row["case"]) != int(labels[i]):
                    raise ValueError(f"Saved label disagrees with metadata in outer/{fold}")
                seen[sample_id] += 1
                collected[sample_id] = {
                    "batch": str(row["batch"]),
                    "y": int(labels[i]),
                    "M4": float(prediction["M4"][i]),
                    "M5": float(prediction["M5"][i]),
                }
                if not np.isfinite([collected[sample_id]["M4"], collected[sample_id]["M5"]]).all():
                    raise ValueError(f"Nonfinite held-out prediction in outer/{fold}")
        input_hashes[prediction_path.name] = sha256(prediction_path)
        show_progress(fold_index, len(folds), f"loaded held-out predictions (outer/{fold})")

    expected_genepark = {row["gsm"] for row in records if row["source"] == "GSE99039"}
    if set(seen) != expected_genepark or any(count != 1 for count in seen.values()):
        raise ValueError("GENEPARK samples are not represented exactly once in outer held-out predictions")

    batch_names = sorted({row["batch"] for row in records if row["source"] == "GSE99039"})
    aucs: dict[str, list[float]] = {"M4": [], "M5": []}
    for batch in batch_names:
        members = [row for row in collected.values() if row["batch"] == batch]
        y = np.asarray([row["y"] for row in members], dtype=int)
        n_control = int(np.sum(y == 0))
        n_pd = int(np.sum(y == 1))
        evaluable = n_control >= MINIMUM_PER_CLASS and n_pd >= MINIMUM_PER_CLASS
        result = {
            "batch": batch,
            "n_control_oof": n_control,
            "n_pd_oof": n_pd,
            "evaluable": evaluable,
            "exclusion_reason": "" if evaluable else f"requires >= {MINIMUM_PER_CLASS} per class",
        }
        if evaluable:
            for model in ("M4", "M5"):
                score = np.asarray([row[model] for row in members], dtype=float)
                result[f"{model}_within_batch_AUROC"] = float(roc_auc_score(y, score))
                aucs[model].append(result[f"{model}_within_batch_AUROC"])
        else:
            result.update({"M4_within_batch_AUROC": "", "M5_within_batch_AUROC": ""})
        output_rows.append(result)

    if not aucs["M4"] or len(aucs["M4"]) != len(aucs["M5"]):
        raise ValueError("No eligible batches or model summaries are inconsistent")
    macro = {model: float(np.mean(values)) for model, values in aucs.items()}
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "batch_results.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    summary = {
        "status": "complete",
        "analysis": "within-batch AUROC on frozen GENEPARK outer out-of-fold predictions",
        "source": "GSE99039",
        "outer_fold_count": 5,
        "sample_count": len(collected),
        "batch_count": len(batch_names),
        "eligible_batch_count": len(aucs["M4"]),
        "minimum_per_class_per_batch": MINIMUM_PER_CLASS,
        "macro_mean_within_batch_AUROC": macro,
        "batch_results": output_rows,
        "split_manifest_sha256": split_hash,
        "input_sha256": input_hashes,
        "new_fits_performed": False,
        "PPMI_accessed": False,
        "interpretation": "Descriptive within-batch ranking diagnostic; not causal attribution and not proof of biological signal.",
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    report = [
        "# GENEPARK within-batch AUROC diagnostic",
        "",
        "Uses saved M4/M5 predictions from the five frozen outer test folds. No refitting or PPMI access.",
        "Each batch AUROC compares PD-control pairs within that batch over pooled out-of-fold samples. Batches with fewer than five samples in either class are listed but excluded from the unweighted macro mean.",
        "",
        f"Eligible batches: **{len(aucs['M4'])}/{len(batch_names)}**.",
        f"Macro mean within-batch AUROC: **M4 {macro['M4']:.3f}; M5 {macro['M5']:.3f}**.",
        "",
        "| Batch | Control n | PD n | Evaluable | M4 AUROC | M5 AUROC |",
        "|---|---:|---:|:---:|---:|---:|",
    ]
    for row in output_rows:
        m4 = f"{row['M4_within_batch_AUROC']:.3f}" if row["evaluable"] else "—"
        m5 = f"{row['M5_within_batch_AUROC']:.3f}" if row["evaluable"] else "—"
        report.append(f"| {row['batch']} | {row['n_control_oof']} | {row['n_pd_oof']} | {'Yes' if row['evaluable'] else 'No'} | {m4} | {m5} |")
    report.extend([
        "",
        "Interpretation: performance above 0.5 within a batch weakens a purely between-batch explanation, but does not establish a causal or biological disease signal. Batch-specific estimates are descriptive and should be read with their sample sizes.",
        "",
        f"Frozen split SHA-256 verified: `{split_hash}`.",
    ])
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n")
    print(f"Wrote aggregate results to {OUT}", flush=True)


if __name__ == "__main__":
    main()
