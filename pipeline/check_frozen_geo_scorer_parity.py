#!/usr/bin/env python3
"""Compare inference-only frozen scorer against the original GEO pipeline.

Reads GEO inputs and saved GEO artifacts only. Does not open PPMI data or print
individual predictions. Emits aggregate numerical differences and hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.special import expit

import geo_pipeline as geo
from frozen_geo_scorer import FrozenGEOScorer, MODELS


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def reference_prep(final_dir: Path) -> geo.Prep:
    prep = geo.Prep.__new__(geo.Prep)
    with np.load(final_dir / "preprocessing.npz", allow_pickle=False) as z:
        for key in ("age_mean", "age_sd", "base_mean", "base_sd", "keep", "selected",
                    "gene_mean", "gene_sd", "genes", "fit_samples_hash"):
            setattr(prep, key, z[key])
    return prep


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-dir", required=True, type=Path)
    parser.add_argument("--marker-manifest", required=True, type=Path)
    parser.add_argument("--sex-artifact", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Parity report already exists")
    data, geo_input_hashes = geo.load_geo()
    scorer = FrozenGEOScorer(args.final_dir, args.marker_manifest, args.sex_artifact)
    actual = scorer.predict_models(data.X, data.age, data.sex, data.genes)

    prep = reference_prep(args.final_dir)
    fits = {}
    for name in MODELS:
        with np.load(args.final_dir / f"{name}_model.npz", allow_pickle=False) as z:
            fits[name] = SimpleNamespace(intercept_=z["intercept"], coef_=z["coef"])
    decisions = json.loads((args.final_dir / "decisions.json").read_text(encoding="utf-8"))
    expected = geo.predict_suite(prep, fits, decisions, data)
    model_differences = {
        name: float(np.max(np.abs(actual[name] - expected[name]))) for name in MODELS
    }

    with np.load(args.sex_artifact, allow_pickle=False) as z:
        sex_design = (data.X[:, z["selected"]] - z["mean"]) / z["sd"]
        sex_expected = expit(z["intercept"] + sex_design @ z["coef"])
    sex_actual = scorer.predict_sex_control(data.X, data.genes)
    sex_difference = float(np.max(np.abs(sex_actual - sex_expected)))

    subset = data.take(np.arange(0, len(data.y), 7))
    subset_actual = scorer.predict_models(subset.X, subset.age, subset.sex, subset.genes)
    subset_difference = max(
        float(np.max(np.abs(subset_actual[name] - actual[name][::7]))) for name in MODELS
    )
    tolerance = 1e-12
    if max(*model_differences.values(), sex_difference, subset_difference) > tolerance:
        raise ValueError("Frozen scorer differs from original GEO path")
    result = {
        "status": "geo_frozen_scorer_parity_passed",
        "geo_record_count": len(data.y),
        "models": list(MODELS),
        "maximum_absolute_probability_difference_by_model": model_differences,
        "sex_control_maximum_absolute_probability_difference": sex_difference,
        "subset_invariance_maximum_absolute_difference": subset_difference,
        "tolerance": tolerance,
        "geo_input_hashes": geo_input_hashes,
        "preprocessing_sha256": digest(args.final_dir / "preprocessing.npz"),
        "model_artifact_sha256": {name: digest(args.final_dir / f"{name}_model.npz") for name in MODELS},
        "decisions_sha256": digest(args.final_dir / "decisions.json"),
        "marker_manifest_sha256": digest(args.marker_manifest),
        "sex_artifact_sha256": digest(args.sex_artifact),
        "scorer_code_sha256": digest(Path(__file__).with_name("frozen_geo_scorer.py")),
        "parity_code_sha256": digest(Path(__file__)),
        "PPMI_opened": False,
        "PPMI_scored": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "geo_record_count": result["geo_record_count"],
        "maximum_absolute_difference": max(*model_differences.values(), sex_difference, subset_difference),
        "PPMI_scored": False,
    }))


if __name__ == "__main__":
    main()
