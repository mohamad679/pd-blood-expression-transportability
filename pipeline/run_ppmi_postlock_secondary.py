#!/usr/bin/env python3
"""Post-lock, aggregate-only PPMI secondary checks on saved GEO predictions.

This program reads restricted rows locally. It writes one owner-only aggregate
JSON and never prints participant records, identifiers, predictions or CBC values.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, mean_absolute_error, roc_auc_score

from audit_ppmi_cbc_coverage import canonical_patno
from audit_ppmi_collection_linkage import parse_month_year
from frozen_geo_scorer import FrozenGEOScorer, MODELS, PANELS
from freeze_ppmi_input_binding import freeze
from paired_bootstrap import paired_bootstrap
from score_ppmi_after_lock import load_locked_input, progress, verify_lock


def calibration(y: np.ndarray, score: np.ndarray) -> dict:
    bins = np.minimum(np.floor(score * 10).astype(int), 9)
    return {
        "AUROC": float(roc_auc_score(y, score)),
        "AUPRC_average_precision": float(average_precision_score(y, score)),
        "log_loss": float(log_loss(y, score, labels=[0, 1])),
        "Brier": float(brier_score_loss(y, score)),
        "calibration_bins": [
            {"bin": i, "n": int(np.sum(bins == i)),
             "mean_predicted_PD": float(score[bins == i].mean()) if np.any(bins == i) else None,
             "observed_PD_fraction": float(y[bins == i].mean()) if np.any(bins == i) else None}
            for i in range(10)
        ],
    }


def matched_sensitivity(y: np.ndarray, age: np.ndarray, sex: np.ndarray,
                        donors: np.ndarray, predictions: dict[str, np.ndarray]) -> dict:
    cases, controls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    edges = [(abs(float(age[i] - age[j])), str(donors[i]), str(donors[j]), int(i), int(j))
             for i in cases for j in controls if sex[i] == sex[j] and abs(age[i] - age[j]) <= 3]
    edges.sort()
    used_case, used_control, pairs = set(), set(), []
    for gap, _case, _control, i, j in edges:
        if i not in used_case and j not in used_control:
            pairs.append((i, j, gap))
            used_case.add(i)
            used_control.add(j)
    result = {
        "rule": "same recorded sex, age difference <=3 years, sorted greedy no replacement",
        "eligible_edges": len(edges), "matched_pairs": len(pairs),
        "unmatched_PD": len(cases) - len(pairs), "unmatched_controls": len(controls) - len(pairs),
        "pairs_female": int(sum(sex[i] == 0 for i, _, _ in pairs)),
        "pairs_male": int(sum(sex[i] == 1 for i, _, _ in pairs)),
        "age_gap_mean": float(np.mean([p[2] for p in pairs])) if pairs else None,
        "age_gap_max": float(max(p[2] for p in pairs)) if pairs else None,
    }
    if len(pairs) < 5:
        result["status"] = "not_evaluable_fewer_than_five_pairs"
        return result
    ix = np.array([i for i, _, _ in pairs] + [j for _, j, _ in pairs], dtype=int)
    result["status"] = "evaluable"
    result["models"] = {name: {"AUROC": float(roc_auc_score(y[ix], pred[ix])),
                               "AUPRC_average_precision": float(average_precision_score(y[ix], pred[ix]))}
                        for name, pred in predictions.items()}
    rng = np.random.default_rng(20261014)
    draws = np.empty(5000)
    for k in range(5000):
        pick = rng.integers(0, len(pairs), len(pairs))
        sample = np.array([pairs[z][0] for z in pick] + [pairs[z][1] for z in pick], dtype=int)
        draws[k] = roc_auc_score(y[sample], predictions["M5"][sample]) - roc_auc_score(y[sample], predictions["M4"][sample])
    result["M5_minus_M4"] = {"delta_AUROC": result["models"]["M5"]["AUROC"] - result["models"]["M4"]["AUROC"],
                               "CI95": np.quantile(draws, [.025, .975]).tolist(),
                               "bootstrap_pairs": 5000, "seed": 20261014,
                               "conditional_on_fixed_models": True}
    return result


def age_control(ranks: np.ndarray, age: np.ndarray, y: np.ndarray,
                genes: np.ndarray, model_paths: list[Path], result_paths: list[Path]) -> dict:
    lookup = {str(g): k for k, g in enumerate(genes)}
    out = {}
    for source, model_path, result_path in zip(("GSE57475", "GSE99039"), model_paths, result_paths, strict=True):
        with np.load(model_path, allow_pickle=False) as z:
            selected_ids = [str(g) for g in z["gene_ids"]]
            if len(selected_ids) != len(set(selected_ids)) or not set(selected_ids).issubset(lookup):
                raise ValueError("Fixed age-model gene map does not match PPMI input")
            ix = np.array([lookup[g] for g in selected_ids])
            mean, sd = z["mean"], z["sd"]
            coef, intercept = z["coef"], float(z["intercept"])
            if len(ix) != len(mean) or len(ix) != len(sd) or len(ix) != len(coef) or np.any(sd <= 0):
                raise ValueError("Invalid fixed age-model dimensions")
            pred = intercept + ((ranks[y == 0][:, ix] - mean) / sd) @ coef
        record = json.loads(result_path.read_text(encoding="utf-8"))
        baseline = float(record["age"]["training_median_age"])
        truth = age[y == 0]
        rho = float(spearmanr(truth, pred).statistic)
        out[source] = {"n_controls": len(truth), "Spearman": rho if np.isfinite(rho) else None,
                       "MAE": float(mean_absolute_error(truth, pred)),
                       "MAE_GEO_training_median": float(mean_absolute_error(truth, np.full(len(truth), baseline))),
                       "GEO_development_direction_capacity": bool(record["age"]["development_direction_capacity"])}
    return out


ASSAYS = {"HMT8": 3, "HMT15": 3, "HMT10": 2, "HMT17": 2, "HMT13": 5,
          "HMT9": 0, "HMT16": 0}


def read_collection_months(path: Path, donors: np.ndarray) -> dict[str, tuple[int, int] | None]:
    wanted = set(map(str, donors))
    records: dict[str, list[tuple[int, int] | None]] = defaultdict(list)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"PATNO", "EVENT_ID", "BLDCLLDT"}.issubset(reader.fieldnames or []):
            raise ValueError("Collection export lacks required fields")
        for row in reader:
            donor = canonical_patno(row.get("PATNO"))
            if donor in wanted and (row.get("EVENT_ID") or "").strip().upper() == "BL":
                records[donor].append(parse_month_year(row.get("BLDCLLDT")))
    return {donor: vals[0] if len(vals) == 1 else None for donor, vals in records.items()}


def cbc_control(ranks: np.ndarray, y: np.ndarray, donors: np.ndarray,
                marker_indices: tuple[np.ndarray, ...], collection: Path, cbc: Path) -> dict:
    controls = {str(donors[k]): k for k in np.flatnonzero(y == 0)}
    months = read_collection_months(collection, donors)
    data: dict[tuple[str, str], list[tuple[int, str, str, float | None, str]]] = defaultdict(list)
    with cbc.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"PATNO", "EVENT_ID", "LTSTCODE", "LCOLLDT", "COLLTM", "LSIRES", "LSIUNIT"}.issubset(reader.fieldnames or []):
            raise ValueError("CBC export lacks required fields")
        for row in reader:
            donor = canonical_patno(row.get("PATNO"))
            code = (row.get("LTSTCODE") or "").strip().upper()
            if donor not in controls or code not in ASSAYS or months.get(donor) is None:
                continue
            if parse_month_year(row.get("LCOLLDT")) != months[donor]:
                continue
            event = (row.get("EVENT_ID") or "").strip().upper()
            tier = 0 if event == "BL" else 1 if event == "SC" else 2
            try:
                value = float((row.get("LSIRES") or "").strip().replace(",", ""))
                if not np.isfinite(value):
                    value = None
            except ValueError:
                value = None
            data[donor, code].append((tier, (row.get("LCOLLDT") or "").strip(),
                                      (row.get("COLLTM") or "").strip(), value,
                                      (row.get("LSIUNIT") or "").strip()))
    out = {"eligible_controls": len(controls), "controls_with_unique_BL_collection_month": sum(months.get(d) is not None for d in controls),
           "assays": {}}
    for code, panel in ASSAYS.items():
        usable, excluded_conflict, excluded_non_numeric, missing = [], 0, 0, 0
        for donor, k in controls.items():
            rows = data.get((donor, code), [])
            if not rows:
                missing += 1
                continue
            top = min(row[0] for row in rows)
            rows = [row for row in rows if row[0] == top]
            occasions = {(row[1], row[2]) for row in rows}
            signatures = {(row[3], row[4]) for row in rows}
            if len(occasions) != 1 or len(signatures) != 1:
                excluded_conflict += 1
                continue
            value, unit = next(iter(signatures))
            if value is None or not unit:
                excluded_non_numeric += 1
                continue
            usable.append((k, value, unit))
        units = sorted({entry[2] for entry in usable})
        summary = {"panel": PANELS[panel], "n_usable_before_unit_check": len(usable),
                   "missing_same_month_assay": missing, "excluded_top_tier_conflict": excluded_conflict,
                   "excluded_non_numeric_or_blank_unit": excluded_non_numeric,
                   "SI_unit_count": len(units)}
        if len(units) != 1:
            summary.update({"status": "not_evaluable_inconsistent_or_missing_SI_units", "n": 0})
            out["assays"][code] = summary
            continue
        indices = np.array([u[0] for u in usable], dtype=int)
        values = np.array([u[1] for u in usable], dtype=float)
        scores = ranks[indices][:, marker_indices[panel]].mean(axis=1)
        summary.update({"SI_unit": units[0], "n": len(values)})
        if len(values) < 3 or len(np.unique(values)) < 2 or len(np.unique(scores)) < 2:
            summary["status"] = "not_evaluable_too_few_or_constant"
        else:
            rho = float(spearmanr(scores, values).statistic)
            rng = np.random.default_rng(20261021)
            draws = []
            for _ in range(5000):
                pick = rng.integers(0, len(values), len(values))
                if len(np.unique(values[pick])) > 1 and len(np.unique(scores[pick])) > 1:
                    draws.append(float(spearmanr(scores[pick], values[pick]).statistic))
            ci = np.quantile(draws, [.025, .975]).tolist() if draws else None
            summary.update({"status": "evaluated" if len(values) >= 30 else "descriptive_below_n_30",
                            "Spearman": rho, "CI95": ci, "valid_bootstrap_draws": len(draws),
                            "supportive": bool(len(values) >= 30 and rho >= .4 and ci is not None and ci[0] > 0)})
        out["assays"][code] = summary
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("published-lock-manifest", "bundle", "final-dir", "marker-manifest", "sex-artifact",
                 "primary-aggregate", "saved-predictions", "age-model-57475", "age-model-99039",
                 "age-results-57475", "age-results-99039", "cbc-current", "collection", "output-dir"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    a = parser.parse_args()
    started = time.monotonic()
    private_root = (Path.home() / "PPMI_private").resolve()
    output = a.output_dir.resolve()
    if private_root not in output.parents or output.exists():
        raise ValueError("Use a new output directory inside ~/PPMI_private")
    progress(0, 8, started, "checking published lock")
    lock = json.loads(a.published_lock_manifest.read_text(encoding="utf-8"))
    verify_lock(lock, a.published_lock_manifest, a.bundle, a.final_dir, a.marker_manifest, a.sex_artifact)
    binding = freeze(a.bundle, a.final_dir / "preprocessing.npz")
    for key in ("ordered_population_sha256", "raw_matrix_sha256", "scored_matrix_geo_array_hash_sha256"):
        if binding[key] != lock[key]:
            raise ValueError(f"Locked input mismatch: {key}")
    ranks, age, sex, y, donors, genes = load_locked_input(a.bundle, binding)
    if binding["group_counts"] != lock["final_population_counts"] or len(y) != 472:
        raise ValueError("Locked population changed")
    progress(1, 8, started, "locked input and population verified")
    with np.load(a.saved_predictions, allow_pickle=False) as z:
        if set(z.files) != set(MODELS) | {"y", "sex", "sex_control"} or not np.array_equal(z["y"], y) or not np.array_equal(z["sex"], sex):
            raise ValueError("Saved prediction order differs from locked input")
        pred = {name: np.asarray(z[name], dtype=float) for name in MODELS}
    if any(p.shape != y.shape or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)) for p in pred.values()):
        raise ValueError("Invalid saved probabilities")
    primary = json.loads(a.primary_aggregate.read_text(encoding="utf-8"))
    if (primary["status"] != "completed_after_sex_gate" or not primary["sex_gate"]["pass"]
            or any(primary[k] != binding[k] for k in ("ordered_population_sha256", "raw_matrix_sha256", "scored_matrix_geo_array_hash_sha256"))):
        raise ValueError("Completed primary result and second lock disagree")
    for name in MODELS:
        if abs(float(roc_auc_score(y, pred[name])) - primary["disease"]["metrics"][name]["AUROC"]) > 1e-12:
            raise ValueError("Saved predictions differ from completed primary result")
    progress(2, 8, started, "saved predictions and primary result verified")
    result = {"status": "complete", "zenodo_second_lock_doi": lock["zenodo_doi"],
              "ordered_population_sha256": binding["ordered_population_sha256"],
              "population": {"PD": int(sum(y)), "Control": int(sum(y == 0)), "PD_fraction": float(y.mean())},
              "models": {name: calibration(y, pred[name]) for name in MODELS}}
    result["M3_minus_M1_secondary"] = paired_bootstrap(y, pred["M1"], pred["M3"], donors,
                                                           repetitions=5000, seed=20261014)[0]
    progress(3, 8, started, "AUPRC, calibration and M3-M1 complete")
    result["age_sex_matched"] = matched_sensitivity(y, age, sex, donors, pred)
    progress(4, 8, started, "age-sex matched sensitivity complete")
    result["age_control"] = age_control(ranks, age, y, genes,
                                         [a.age_model_57475, a.age_model_99039],
                                         [a.age_results_57475, a.age_results_99039])
    progress(5, 8, started, "pretrained GEO age controls complete")
    scorer = FrozenGEOScorer(a.final_dir, a.marker_manifest, a.sex_artifact)
    if not np.array_equal(genes, scorer.genes):
        raise ValueError("CBC marker gene order differs from frozen GEO")
    result["CBC_controls"] = cbc_control(ranks, y, donors, scorer.marker_indices, a.collection, a.cbc_current)
    progress(6, 8, started, "CBC controls complete")
    result["longitudinal_similarity"] = {"status": "unavailable",
                                         "reason": "locked baseline bundle has no independently verified nonbaseline autosomal ranks"}
    result["interpretation_scope"] = "All CIs descriptive or secondary and conditional on fixed GEO models; primary endpoint unchanged"
    progress(7, 8, started, "aggregate report assembled")
    os.umask(0o077)
    output.mkdir(mode=0o700)
    destination = output / "AGGREGATE_SECONDARY.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    destination.chmod(0o600)
    progress(8, 8, started, "owner-only aggregate report saved")
    print(json.dumps({"status": result["status"], "aggregate_output": str(destination)}))


if __name__ == "__main__":
    main()
