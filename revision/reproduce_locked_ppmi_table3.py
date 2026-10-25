#!/usr/bin/env python3
"""Independently rescore locked PPMI with the GEO reference path; aggregate only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from sklearn.metrics import log_loss, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from check_frozen_geo_scorer_parity import reference_prep
from score_ppmi_after_lock import load_locked_input
from postlock_common import private_output, require_execution_authorization, write_aggregate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization-record", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--final-dir", required=True, type=Path)
    parser.add_argument("--marker-manifest", required=True, type=Path)
    parser.add_argument("--stored-predictions", required=True, type=Path)
    parser.add_argument("--primary-aggregate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require_execution_authorization(args.authorization_record, analysis="TABLE3_AUDIT")
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    ranks, age, sex, y, donors, genes = load_locked_input(args.bundle, binding)
    manifest = json.loads(args.marker_manifest.read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    data = geo.Data(ranks, age, sex, y, np.repeat("PPMI", len(y)), donors,
                    donors.copy(), genes, markers)
    prep = reference_prep(args.final_dir)
    fits = {}
    for name in geo.MODELS:
        with np.load(args.final_dir / f"{name}_model.npz", allow_pickle=False) as model:
            fits[name] = SimpleNamespace(intercept_=model["intercept"], coef_=model["coef"])
    decisions = json.loads((args.final_dir / "decisions.json").read_text(encoding="utf-8"))
    new = geo.predict_suite(prep, fits, decisions, data)
    prior = json.loads(args.primary_aggregate.read_text(encoding="utf-8"))
    metrics = prior["disease"]["metrics"]
    differences = {}
    with np.load(args.stored_predictions, allow_pickle=False) as saved:
        if not np.array_equal(saved["y"], y):
            raise ValueError("Stored predictions differ in row order or diagnosis")
        for name in geo.MODELS:
            prediction = np.asarray(new[name], float)
            differences[name] = {
                "maximum_absolute_probability_difference": float(np.max(np.abs(prediction - saved[name]))),
                "absolute_AUROC_difference": float(abs(roc_auc_score(y, prediction) - metrics[name]["AUROC"])),
                "absolute_log_loss_difference": float(abs(log_loss(y, prediction, labels=[0, 1]) - metrics[name]["log_loss"])),
            }
    delta = float(roc_auc_score(y, new["M5"]) - roc_auc_score(y, new["M4"]))
    delta_difference = abs(delta - float(prior["disease"]["primary"]["delta_AUROC"]))
    maximum = max([delta_difference] + [v for row in differences.values() for v in row.values()])
    result = {"analysis": "TABLE3_AUDIT", "status": "PASS" if maximum <= 1e-10 else "FAIL",
              "independent_prediction_path": "geo_pipeline.Prep.transform and predict_suite",
              "input_matrix_sha256": binding["raw_matrix_sha256"],
              "model_differences": differences,
              "absolute_paired_delta_AUROC_difference": delta_difference,
              "tolerance": 1e-10,
              "note": "This verifies fixed-model scoring and aggregate point metrics; E7 separately reproduces bootstrap CI quantiles."}
    write_aggregate(output, result)
    print("TABLE3_INDEPENDENT_REPRODUCTION:", result["status"])
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
