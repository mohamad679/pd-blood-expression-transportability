"""Descriptive paired donor bootstrap of saved GEO-to-GEO transfer scores.

The fixed decision is docs/GEO_TRANSFER_BOOTSTRAP_DECISION_20261020_EN.md.
No model is fitted and no PPMI data are accessed.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from paired_bootstrap import auc_multiplicities


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_REL = "pipeline/bootstrap_geo_transfer_descriptive.py"
DECISION_REL = "docs/GEO_TRANSFER_BOOTSTRAP_DECISION_20261020_EN.md"
OUTPUT = ROOT / "results/geo_transfer_descriptive_bootstrap_20261020.json"
REPETITIONS = 5000
MODELS = ("M1", "M2", "M3", "M4", "M5")
DIRECTIONS = (
    ("GSE57475", "GSE99039", "GSE57475_transfer", 20261020),
    ("GSE99039", "GSE57475", "GSE99039_transfer", 20261021),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=ROOT, text=True).strip()


def check_published_code() -> str:
    head = git_value("rev-parse", "HEAD")
    if head != git_value("rev-parse", "origin/main"):
        raise ValueError("Commit must be pushed to origin/main before calculation")
    for relative in (SCRIPT_REL, DECISION_REL):
        committed = subprocess.check_output(("git", "show", f"HEAD:{relative}"), cwd=ROOT)
        if committed != (ROOT / relative).read_bytes():
            raise ValueError(f"Uncommitted or mismatched decision/code: {relative}")
    return head


def bootstrap_direction(source: str, target: str, directory: str, seed: int,
                        original: dict) -> dict:
    predictions = ROOT / "private_data/geo_runs_20261018" / directory / "heldout_predictions.npz"
    metrics_file = predictions.with_name("RESULTS.json")
    metadata_file = ROOT / "prepared" / f"{target}_sample_metadata.json"
    expected_metadata_hash = original["input_sha256"][metadata_file.name]
    if sha256(metadata_file) != expected_metadata_hash:
        raise ValueError(f"Target metadata hash mismatch: {target}")

    with np.load(predictions, allow_pickle=False) as saved:
        if set(saved.files) != {"samples", "y", *MODELS}:
            raise ValueError(f"Unexpected saved prediction fields: {directory}")
        samples = [str(sample) for sample in saved["samples"]]
        y = np.asarray(saved["y"], dtype=int)
        scores = {model: np.asarray(saved[model], dtype=float) for model in MODELS}

    metadata_rows = json.loads(metadata_file.read_text())
    metadata = {str(row["gsm"]): row for row in metadata_rows}
    if (len(metadata) != len(metadata_rows) or len(samples) != len(metadata)
            or len(set(samples)) != len(samples) or set(samples) != set(metadata)
            or set(np.unique(y)) != {0, 1}):
        raise ValueError(f"Target sample or outcome coverage mismatch: {directory}")
    donors = [f"{target}:{metadata[sample]['subject_id'] or sample}" for sample in samples]
    if len(set(donors)) != len(donors):
        raise ValueError(f"Repeated operational donor in {directory}")
    if any(int(metadata[sample]["case"]) != int(label)
           for sample, label in zip(samples, y)):
        raise ValueError(f"Saved diagnosis differs from frozen metadata: {directory}")
    for model, values in scores.items():
        if (values.shape != y.shape or not np.isfinite(values).all()
                or np.any((values < 0) | (values > 1))):
            raise ValueError(f"Invalid saved probability in {directory}/{model}")

    prior_metrics = json.loads(metrics_file.read_text())["metrics"]
    observed = {model: float(auc_multiplicities(
        y, scores[model], np.ones((1, len(y)), dtype=int))[0]) for model in MODELS}
    for model in MODELS:
        if abs(observed[model] - float(prior_metrics[model]["AUROC"])) > 1e-12:
            raise ValueError(f"Saved AUROC does not reproduce: {directory}/{model}")

    rng = np.random.default_rng(seed)
    strata = [np.flatnonzero(y == value) for value in (0, 1)]
    draws = {model: np.empty(REPETITIONS, dtype=float) for model in MODELS}
    for start in range(0, REPETITIONS, 100):
        stop = min(REPETITIONS, start + 100)
        counts = np.zeros((stop - start, len(y)), dtype=int)
        for index in strata:
            counts[:, index] = rng.multinomial(
                len(index), np.full(len(index), 1 / len(index)), size=stop - start
            )
        for model in MODELS:
            draws[model][start:stop] = auc_multiplicities(y, scores[model], counts)

    differences = draws["M5"] - draws["M4"]
    intervals = {model: [float(value) for value in np.quantile(draws[model], [.025, .975])]
                 for model in MODELS}
    return {
        "direction": f"{source} -> {target}",
        "test_n": len(y),
        "test_PD_n": int(y.sum()),
        "test_control_n": int(len(y) - y.sum()),
        "unique_test_donors": len(donors),
        "seed": seed,
        "observed_AUROC": observed,
        "percentile_95_CI_AUROC": intervals,
        "observed_M5_minus_M4_AUROC": observed["M5"] - observed["M4"],
        "percentile_95_CI_M5_minus_M4_AUROC": [
            float(value) for value in np.quantile(differences, [.025, .975])
        ],
        "prediction_sha256": sha256(predictions),
        "original_result_sha256": sha256(metrics_file),
        "target_metadata_sha256": sha256(metadata_file),
    }


def main() -> None:
    if OUTPUT.exists():
        raise ValueError(f"Refusing to overwrite existing result: {OUTPUT}")
    commit = check_published_code()
    original_file = ROOT / "results/geo_development_20261018/aggregate_summary.json"
    original = json.loads(original_file.read_text())
    rows = []
    for number, (source, target, directory, seed) in enumerate(DIRECTIONS, start=1):
        row = bootstrap_direction(source, target, directory, seed, original)
        rows.append(row)
        print(f"Descriptive GEO transfer bootstrap {number}/2: {row['direction']} complete", flush=True)
    report = {
        "status": "complete",
        "created_UTC": datetime.now(timezone.utc).isoformat(),
        "purpose": "descriptive fixed-model GEO-to-GEO transfer uncertainty only",
        "method": "5000 paired donor-level, case/control-stratified resamples; percentile 95% CI",
        "repetitions": REPETITIONS,
        "conditional_on_fixed_fitted_models": True,
        "confirmatory_PPMI_analysis": False,
        "PPMI_accessed": False,
        "git_commit": commit,
        "code_sha256": sha256(ROOT / SCRIPT_REL),
        "decision_sha256": sha256(ROOT / DECISION_REL),
        "original_summary_sha256": sha256(original_file),
        "directions": rows,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(f"GEO_TRANSFER_DESCRIPTIVE_BOOTSTRAP_COMPLETE: {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
