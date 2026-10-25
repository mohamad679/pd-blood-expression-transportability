#!/usr/bin/env python3
"""E3: post-lock exploratory nested M4/M5 CV; synthetic check only in phase 2.

Real PPMI execution fails closed without a separate phase-3 authorization record.
The program prints only aggregate progress, never donor IDs or row predictions.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
import warnings
import zipfile
from collections import Counter

import numpy as np
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from score_ppmi_after_lock import load_locked_input

from postlock_common import private_output, require_execution_authorization, rng_for, summary, write_aggregate

GRID = (.0001, .001, .01, .1, 1.)
RHOS = (0., .5, 1.)
PREDICTION_WARNINGS: Counter[str] = Counter()


def donor_folds(data: geo.Data, n_splits: int, seed: int):
    donors, inverse = np.unique(data.donor, return_inverse=True)
    labels = np.empty(len(donors), dtype=int)
    for donor_ix in range(len(donors)):
        values = np.unique(data.y[inverse == donor_ix])
        if len(values) != 1:
            raise ValueError("A donor has inconsistent diagnoses")
        labels[donor_ix] = values[0]
    counts = np.bincount(labels, minlength=2)
    if len(counts) != 2 or min(counts) < n_splits:
        raise ValueError("Too few donors per class for stratification")
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train_d, test_d in folds.split(donors, labels):
        train = np.flatnonzero(np.isin(inverse, train_d))
        test = np.flatnonzero(np.isin(inverse, test_d))
        if set(data.donor[train]) & set(data.donor[test]):
            raise AssertionError("Donor leakage")
        yield train, test


def plate_folds(data: geo.Data, plates: np.ndarray, n_splits: int, seed: int):
    if plates.shape != data.y.shape or len(set(plates)) < n_splits:
        raise ValueError("Insufficient recorded plates for grouped CV")
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                   random_state=seed)
    for train, test in splitter.split(np.zeros(len(data.y)), data.y, groups=plates):
        if set(plates[train]) & set(plates[test]):
            raise AssertionError("Plate leakage")
        if set(data.donor[train]) & set(data.donor[test]):
            raise AssertionError("Donor leakage in plate CV")
        if len(np.unique(data.y[train])) != 2 or len(np.unique(data.y[test])) != 2:
            raise ValueError("Plate-grouped fold lacks one diagnosis class")
        yield train, test


def plate_prevalence_score(train_y: np.ndarray, train_plates: np.ndarray,
                           test_plates: np.ndarray) -> np.ndarray:
    global_prevalence = float(np.mean(train_y))
    rates = {str(plate): float(np.mean(train_y[train_plates == plate]))
             for plate in np.unique(train_plates)}
    return np.array([rates.get(str(plate), global_prevalence)
                     for plate in test_plates], dtype=float)


def paired_oof_interval(y: np.ndarray, p4: np.ndarray, p5: np.ndarray,
                        draws: int, regime: int) -> tuple[float, list[float]]:
    if not (len(y) == len(p4) == len(p5)) or len(np.unique(y)) != 2:
        raise ValueError("Invalid OOF score arrays")
    point = float(roc_auc_score(y, p5) - roc_auc_score(y, p4))
    groups = [np.flatnonzero(y == label) for label in (0, 1)]
    rng = rng_for(3, 50000 + regime)
    values = np.empty(draws, float)
    for i in range(draws):
        ix = np.concatenate([rng.choice(group, len(group), replace=True)
                             for group in groups])
        values[i] = roc_auc_score(y[ix], p5[ix]) - roc_auc_score(y[ix], p4[ix])
    return point, np.quantile(values, [.025, .975]).astype(float).tolist()


def locked_plates(bundle: Path, donors: np.ndarray) -> np.ndarray:
    with zipfile.ZipFile(bundle) as archive:
        rows = [json.loads(line) for line in archive.read(
            "candidate_metadata.jsonl").decode("utf-8").splitlines()]
    if len(rows) != len(donors) or any(str(row["PATNO"]) != str(donor)
                                       for row, donor in zip(rows, donors)):
        raise ValueError("Plate metadata differs from locked donor order")
    plates = np.array([str(row["plate"]).strip() for row in rows])
    if any(not plate for plate in plates):
        raise ValueError("Recorded plate is missing")
    return plates


def score_model(prep, fit, data, name: str) -> np.ndarray:
    x, _ = prep.transform(data, name)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        p = np.asarray(fit.predict(x), dtype=float)
    for item in captured:
        PREDICTION_WARNINGS[f"{item.category.__name__}: {item.message}"] += 1
    if p.shape != data.y.shape or not np.isfinite(p).all():
        raise ValueError("Invalid validation prediction")
    return p


def validation_score(train: geo.Data, validation: geo.Data,
                     prediction: np.ndarray, objective: str) -> float:
    if objective == "AUROC":
        return float(roc_auc_score(validation.y, prediction))
    if objective != "weighted_deviance":
        raise ValueError("Unknown tuning objective")
    sources = np.unique(train.source)
    weights = np.zeros(len(validation.y), float)
    for source in sources:
        for label in (0, 1):
            train_n = int(((train.source == source) & (train.y == label)).sum())
            mask = (validation.source == source) & (validation.y == label)
            if train_n == 0:
                raise ValueError("Training source lacks a class")
            weights[mask] = len(train.y) / (2 * len(sources) * train_n)
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Validation source is absent from training")
    # `geo.choose` maximizes; negative weighted log-loss implements minimization.
    return -float(log_loss(validation.y, prediction, labels=[0, 1],
                           sample_weight=weights))


def checked_fit_with_warning_log(prep, data, model, lam, rho, context,
                                 fit_log: Path, warning_log: Path):
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = geo.checked_fit(prep, data, model, lam, rho, context, fit_log)
    if captured:
        counts = {}
        for item in captured:
            key = f"{item.category.__name__}: {item.message}"
            counts[key] = counts.get(key, 0) + 1
        with warning_log.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps({"stage": context["stage"], "model": model,
                                     "counts": counts}, sort_keys=True) + "\n")
    return result


def fit_pair(train: geo.Data, inner, output: Path, scope: str,
             max_genes: int = 1000, grid=GRID, rhos=RHOS,
             expansion_rounds: int = 2, prep_factory=None,
             objective: str = "AUROC"):
    """Reproduce v3.1 baseline-first tuning for only M4 and M5."""
    output.mkdir(parents=True, exist_ok=False)
    fit_log = output / "fits.jsonl"
    warning_log = output / "warnings.jsonl"
    if prep_factory is None:
        prep_factory = lambda data, max_genes: geo.Prep(data, max_genes=max_genes)
    fold_data = []
    baseline_scores = []
    for fold, (ti, vi) in enumerate(inner):
        a, b = train.take(ti), train.take(vi)
        if set(a.donor) & set(b.donor):
            raise AssertionError("Inner donor leakage")
        prep = prep_factory(a, max_genes)
        m4, _ = checked_fit_with_warning_log(prep, a, "M4", 0., 0.,
                                dict(scope=scope, stage="inner", fold=fold), fit_log,
                                warning_log)
        baseline_scores.append(validation_score(
            a, b, score_model(prep, m4, b, "M4"), objective))
        fold_data.append((prep, a, b))
    candidates = [dict(baseline=True, **{"lambda": 0.}, rho=0.,
                       score=float(np.mean(baseline_scores)), nonzero=False)]
    lambdas = list(grid)
    for round_i in range(expansion_rounds + 1):
        existing = {(c["lambda"], c["rho"]) for c in candidates if not c["baseline"]}
        for lam in sorted(lambdas, reverse=True):
            for rho in rhos:
                if (lam, rho) in existing:
                    continue
                scores, nonzero = [], False
                for fold, (prep, a, b) in enumerate(fold_data):
                    fit, records = checked_fit_with_warning_log(prep, a, "M5", lam, rho,
                        dict(scope=scope, stage="inner", fold=fold,
                             expansion_round=round_i), fit_log, warning_log)
                    scores.append(validation_score(
                        a, b, score_model(prep, fit, b, "M5"), objective))
                    nonzero |= records[-1]["diagnostics"]["nonzero_penalized_coefficients"] > 0
                candidates.append(dict(baseline=False, **{"lambda": float(lam)},
                                       rho=float(rho), score=float(np.mean(scores)),
                                       nonzero=bool(nonzero)))
        selected = geo.choose(candidates, tolerance=1e-6)
        edge = (not selected["baseline"] and selected["nonzero"]
                and selected["lambda"] in (min(lambdas), max(lambdas)))
        if not edge or round_i == expansion_rounds:
            break
        lambdas.append(min(lambdas) / 10 if selected["lambda"] == min(lambdas)
                       else max(lambdas) * 10)
    prep = prep_factory(train, max_genes)
    m4, _ = checked_fit_with_warning_log(prep, train, "M4", 0., 0.,
                            dict(scope=scope, stage="outer_refit"), fit_log,
                            warning_log)
    if selected["baseline"]:
        m5 = None
    else:
        m5, records = checked_fit_with_warning_log(prep, train, "M5", selected["lambda"],
                                      selected["rho"],
                                      dict(scope=scope, stage="outer_refit"), fit_log,
                                      warning_log)
        if records[-1]["diagnostics"]["nonzero_penalized_coefficients"] == 0:
            m5 = None
    decision = {"selected": selected, "persistent_boundary": bool(edge),
                "expansion_rounds": round_i, "candidate_count": len(candidates),
                "M5_reuses_M4": m5 is None}
    return prep, m4, m5, decision


def run_regime(data: geo.Data, plates: np.ndarray, output: Path, regime: str,
               repetitions: int, outer_folds: int, inner_folds: int,
               max_genes: int, grid, rhos, expansion_rounds: int,
               bootstrap_draws: int) -> dict:
    if regime not in ("donor", "plate_grouped"):
        raise ValueError("Unknown E3 CV regime")
    output = Path(output)
    if output.exists():
        raise ValueError("Use a fresh E3 output directory")
    output.mkdir(parents=True)
    started = time.monotonic()
    repetition_means = []
    plate_only_means = []
    plate_only_fold_values = []
    oof4 = np.zeros(len(data.y), float)
    oof5 = np.zeros(len(data.y), float)
    boundary_count = 0
    for rep in range(repetitions):
        stream_offset = 0 if regime == "donor" else 100000
        seed = int(rng_for(3, stream_offset + rep).integers(0, 2**32 - 1))
        paired = []
        plate_scores = []
        seen = np.zeros(len(data.y), dtype=np.int8)
        outer = (donor_folds(data, outer_folds, seed) if regime == "donor" else
                 plate_folds(data, plates, outer_folds, seed))
        for fold, (ti, vi) in enumerate(outer):
            train, valid = data.take(ti), data.take(vi)
            inner_seed = int(rng_for(3, stream_offset + rep * 100 + fold + 1)
                             .integers(0, 2**32 - 1))
            inner = list(donor_folds(train, inner_folds, inner_seed) if regime == "donor"
                         else plate_folds(train, plates[ti], inner_folds, inner_seed))
            scope = f"E3/{regime}/repetition/{rep}/outer/{fold}"
            prep, m4, m5, decision = fit_pair(train, inner,
                output / f"rep_{rep:02d}_fold_{fold}", scope, max_genes,
                grid, rhos, expansion_rounds)
            p4 = score_model(prep, m4, valid, "M4")
            p5 = p4 if m5 is None else score_model(prep, m5, valid, "M5")
            oof4[vi] += p4
            oof5[vi] += p5
            seen[vi] += 1
            paired.append(float(roc_auc_score(valid.y, p5) - roc_auc_score(valid.y, p4)))
            plate_p = plate_prevalence_score(train.y, plates[ti], plates[vi])
            plate_scores.append(float(roc_auc_score(valid.y, plate_p)))
            boundary_count += int(decision["persistent_boundary"])
        if not np.all(seen == 1):
            raise AssertionError("Each donor must have exactly one OOF score per repetition")
        repetition_means.append(float(np.mean(paired)))
        plate_only_means.append(float(np.mean(plate_scores)))
        plate_only_fold_values.append(plate_scores.copy())
        elapsed = time.monotonic() - started
        eta = elapsed / (rep + 1) * (repetitions - rep - 1)
        print(f"E3 {regime} progress {rep + 1}/{repetitions} repetitions; "
              f"elapsed {elapsed / 3600:.1f}h; ETA ~{eta / 3600:.1f}h", flush=True)
    p4 = oof4 / repetitions
    p5 = oof5 / repetitions
    point, ci = paired_oof_interval(data.y, p4, p5, bootstrap_draws,
                                     0 if regime == "donor" else 1)
    return {"status": "complete_exploratory", "CV_regime": regime,
            "seed_rule": "SeedSequence([20261014, 3, regime_offset + repetition])",
            "outer_repetitions": repetitions, "outer_folds": outer_folds,
            "inner_folds": inner_folds, "mean_paired_fold_delta_by_repetition": repetition_means,
            "repetition_summary": summary(np.array(repetition_means)),
            "range": [float(min(repetition_means)), float(max(repetition_means))],
            "plate_only_fold_AUROC_by_repetition": plate_only_fold_values,
            "plate_only_mean_fold_AUROC_by_repetition": plate_only_means,
            "plate_only_repetition_summary": summary(np.array(plate_only_means)),
            "mean_OOF_M4_AUROC": float(roc_auc_score(data.y, p4)),
            "mean_OOF_M5_AUROC": float(roc_auc_score(data.y, p5)),
            "mean_OOF_paired_delta_AUROC": point,
            "mean_OOF_paired_donor_bootstrap_CI95": ci,
            "bootstrap_draws": bootstrap_draws,
            "bootstrap_limitation": "Conditional on fitted CV models; omits training variability",
            "persistent_boundary_count": boundary_count,
            "prediction_warning_counts": dict(PREDICTION_WARNINGS)}


def interpretation_band(delta: float, ci: list[float]) -> str:
    if delta >= .05:
        return "within_PPMI_increment_not_biological_without_plate_control"
    if -.02 <= delta <= .02 and ci[1] < .05:
        return "increment_of_0.05_excluded_conditionally"
    if delta < -.02:
        return "adverse_point_estimate_uncertainty_reported_separately"
    return "inconclusive_increment"


def plate_comparison(ordinary_delta: float, grouped: dict) -> str:
    if grouped.get("status") != "complete_exploratory":
        return "plate_grouped_not_evaluable_no_attribution"
    grouped_delta = grouped["mean_OOF_paired_delta_AUROC"]
    if ordinary_delta >= .05 and grouped_delta < .02:
        return "ordinary_increment_plate_associated_not_biological_attribution"
    if ordinary_delta >= .05 and grouped_delta >= .05:
        return "increment_survives_plate_grouping_not_biological_attribution"
    if ordinary_delta >= .05:
        return "ordinary_increment_attenuated_plate_role_inconclusive"
    return "no_ordinary_increment_at_0.05_threshold_plate_results_descriptive"


def run_nested(data: geo.Data, plates: np.ndarray, output: Path,
               repetitions: int = 20, outer_folds: int = 5, inner_folds: int = 5,
               max_genes: int = 1000, grid=GRID, rhos=RHOS,
               expansion_rounds: int = 2, bootstrap_draws: int = 5000) -> dict:
    if len(np.unique(data.donor)) != len(data.y):
        raise ValueError("E3 requires one locked record per donor")
    output = Path(output)
    if output.exists():
        raise ValueError("Use a fresh E3 output directory")
    output.mkdir(parents=True)
    ordinary = run_regime(data, plates, output / "donor_cv", "donor",
                          repetitions, outer_folds, inner_folds, max_genes,
                          grid, rhos, expansion_rounds, bootstrap_draws)
    try:
        grouped = run_regime(data, plates, output / "plate_cv", "plate_grouped",
                             repetitions, outer_folds, inner_folds, max_genes,
                             grid, rhos, expansion_rounds, bootstrap_draws)
    except ValueError as error:
        grouped = {"status": "not_evaluable", "reason_type": type(error).__name__,
                   "reason": "Plate-grouped fold or fit not evaluable under fixed rules"}
    ordinary["interpretation_band"] = interpretation_band(
        ordinary["mean_OOF_paired_delta_AUROC"],
        ordinary["mean_OOF_paired_donor_bootstrap_CI95"])
    ordinary["plate_comparison"] = plate_comparison(
        ordinary["mean_OOF_paired_delta_AUROC"], grouped)
    return {"analysis": "E3", "status": "complete_exploratory",
            "donor_CV": ordinary, "plate_grouped_CV": grouped}


def synthetic_self_test() -> None:
    from postlock_common import synthetic_ranks
    rng = np.random.default_rng(17)
    n, g = 48, 40
    ranks = synthetic_ranks(n, g, seed=17)
    y = np.array([0, 1] * (n // 2), dtype=int)
    data = geo.Data(ranks, rng.uniform(50, 75, n), rng.integers(0, 2, n), y,
                    np.repeat("SYNTHETIC", n), np.array([f"synthetic-{i}" for i in range(n)]),
                    np.array([f"synthetic-{i}" for i in range(n)]),
                    np.array([f"G{i:04d}" for i in range(g)]),
                    [np.arange(i * 3, i * 3 + 3) for i in range(6)])
    folds = list(donor_folds(data, 2, 1234))
    assert len(folds) == 2
    assert all(not set(data.donor[a]) & set(data.donor[b]) for a, b in folds)
    arbitrary = np.where(data.y == 1, .7, .3)
    assert validation_score(data, data, arbitrary, "weighted_deviance") < 0
    plates = np.array([f"plate-{i//4}" for i in range(n)])
    grouped_folds = list(plate_folds(data, plates, 2, 1234))
    assert len(grouped_folds) == 2
    assert all(not set(plates[a]) & set(plates[b]) for a, b in grouped_folds)
    assert np.allclose(plate_prevalence_score(data.y[:4], plates[:4],
                                             np.array(["unseen"])), .5)
    point, ci = paired_oof_interval(data.y, arbitrary, arbitrary, 20, 0)
    assert point == 0. and ci == [0., 0.]
    assert plate_comparison(.05, {"status": "complete_exploratory",
                                  "mean_OOF_paired_delta_AUROC": .019}).startswith("ordinary_increment_plate_associated")
    assert plate_comparison(.05, {"status": "complete_exploratory",
                                  "mean_OOF_paired_delta_AUROC": .05}).startswith("increment_survives_plate_grouping")
    assert plate_comparison(.05, {"status": "not_evaluable"}) == "plate_grouped_not_evaluable_no_attribution"
    # The optimizer is exercised on a small artificial fold, not PPMI or GEO.
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        train = data.take(folds[0][0])
        inner = list(donor_folds(train, 2, 5678))
        prep, m4, m5, decision = fit_pair(train, inner, Path(temporary) / "fit",
            "E3/synthetic", max_genes=8, grid=(.1,), rhos=(0.,), expansion_rounds=0)
        valid = data.take(folds[0][1])
        assert np.isfinite(score_model(prep, m4, valid, "M4")).all()
        if m5 is not None:
            assert np.isfinite(score_model(prep, m5, valid, "M5")).all()
        assert decision["candidate_count"] == 2
        result = run_nested(data, plates, Path(temporary) / "nested", repetitions=1,
            outer_folds=2, inner_folds=2, max_genes=8, grid=(.1,),
            rhos=(0.,), expansion_rounds=0, bootstrap_draws=20)
        assert result["status"] == "complete_exploratory"
        assert len(result["donor_CV"]["mean_paired_fold_delta_by_repetition"]) == 1
        assert len(result["donor_CV"]["plate_only_fold_AUROC_by_repetition"][0]) == 2
    print("E3_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(x is None for x in (args.authorization_record, args.bundle, args.binding, args.output)):
        parser.error("Real run requires authorization record, bundle, binding and output")
    require_execution_authorization(args.authorization_record, analysis="E3")
    os.umask(0o077)
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    ranks, age, sex, y, donors, genes = load_locked_input(args.bundle, binding)
    plates = locked_plates(args.bundle, donors)
    manifest = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    data = geo.Data(ranks, age, sex, y, np.repeat("PPMI", len(y)), donors,
                    donors.copy(), genes, markers)
    result = run_nested(data, plates, output)
    write_aggregate(output / "E3_AGGREGATE.json", result)
    print("E3_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
