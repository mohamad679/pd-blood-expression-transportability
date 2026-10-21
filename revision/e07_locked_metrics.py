#!/usr/bin/env python3
"""E7: descriptive metrics and calibration from locked predictions only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time
import warnings

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             log_loss, roc_auc_score)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from paired_bootstrap import auc_multiplicities, paired_bootstrap
from score_ppmi_after_lock import load_locked_input
from postlock_common import private_output, require_execution_authorization, write_aggregate

MODELS = ("M1", "M2", "M3", "M4", "M5")
GRID = np.arange(.05, .951, .05)


def primary_counts(y: np.ndarray, repetitions: int = 5000) -> np.ndarray:
    """Exactly reconstruct paired_bootstrap.py's seed and 100-draw blocks."""
    y = np.asarray(y)
    rng = np.random.default_rng(20261014)
    classes = [np.flatnonzero(y == c) for c in (0, 1)]
    blocks = []
    for offset in range(0, repetitions, 100):
        b = min(100, repetitions - offset)
        counts = np.zeros((b, len(y)), dtype=int)
        for ix in classes:
            counts[:, ix] = rng.multinomial(len(ix), np.full(len(ix), 1 / len(ix)), size=b)
        blocks.append(counts)
    return np.vstack(blocks)


def basic_metrics(y: np.ndarray, score: np.ndarray) -> dict:
    return {"AP": float(average_precision_score(y, score)),
            "Brier": float(brier_score_loss(y, score)),
            "log_loss": float(log_loss(y, score, labels=[0, 1]))}


def calibration_fit(y: np.ndarray, score: np.ndarray,
                    weights: np.ndarray | None = None) -> tuple[float, float] | None:
    clipped = np.clip(score, 1e-12, 1 - 1e-12)
    logit = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            fit = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000,
                                     tol=1e-9).fit(logit, y, sample_weight=weights)
        if any("ConvergenceWarning" in item.category.__name__ for item in captured):
            return None
        intercept, slope = float(fit.intercept_[0]), float(fit.coef_[0, 0])
        return (intercept, slope) if np.isfinite([intercept, slope]).all() else None
    except (ValueError, FloatingPointError, OverflowError):
        return None


def loess_at(x: np.ndarray, y: np.ndarray, grid: np.ndarray,
             span: float = .5) -> tuple[np.ndarray, np.ndarray]:
    """Fixed-span local linear tricube smoother; effective n is Kish n."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    result, effective = np.full(len(grid), np.nan), np.zeros(len(grid))
    k = max(2, int(np.ceil(span * len(x))))
    for j, point in enumerate(grid):
        distance = np.abs(x - point)
        h = float(np.partition(distance, k - 1)[k - 1])
        if h == 0:
            weights = (distance == 0).astype(float)
        else:
            u = np.minimum(distance / h, 1.)
            weights = (1 - u**3)**3
        total = weights.sum()
        if total <= 0:
            continue
        effective[j] = total**2 / np.sum(weights**2)
        if effective[j] < 20:
            continue
        centered = x - point
        s0 = total
        s1 = float(np.sum(weights * centered))
        s2 = float(np.sum(weights * centered**2))
        t0 = float(np.sum(weights * y))
        t1 = float(np.sum(weights * centered * y))
        determinant = s0 * s2 - s1**2
        result[j] = ((t0 * s2 - t1 * s1) / determinant if determinant > 1e-12
                     else t0 / s0)
    return result, effective


def matched_pairs(age: np.ndarray, sex: np.ndarray, y: np.ndarray,
                  donors: np.ndarray) -> list[tuple[int, int]]:
    cases, controls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    edges = [(abs(float(age[i] - age[j])), str(donors[i]), str(donors[j]), int(i), int(j))
             for i in cases for j in controls if sex[i] == sex[j] and abs(age[i] - age[j]) <= 3]
    edges.sort()
    used_case, used_control, pairs = set(), set(), []
    for _gap, _case, _control, i, j in edges:
        if i not in used_case and j not in used_control:
            pairs.append((i, j))
            used_case.add(i)
            used_control.add(j)
    return pairs


def matched_auc_difference(y: np.ndarray, p1: np.ndarray, p4: np.ndarray,
                           pairs: list[tuple[int, int]], draws: int = 5000) -> dict:
    if len(pairs) < 5:
        return {"status": "not_evaluable_fewer_than_five_pairs", "n_pairs": len(pairs)}
    ix = np.array([i for i, _ in pairs] + [j for _, j in pairs], int)
    rng = np.random.default_rng(20261014)
    differences = np.empty(draws)
    for draw in range(draws):
        pick = rng.integers(0, len(pairs), len(pairs))
        sample = np.array([pairs[k][0] for k in pick] + [pairs[k][1] for k in pick], int)
        differences[draw] = (roc_auc_score(y[sample], p4[sample])
                             - roc_auc_score(y[sample], p1[sample]))
    return {"status": "complete_descriptive", "n_pairs": len(pairs),
            "M4_minus_M1_AUROC": float(roc_auc_score(y[ix], p4[ix])
                                       - roc_auc_score(y[ix], p1[ix])),
            "CI95_percentile": np.quantile(differences, [.025, .975]).tolist()}


def compute(y: np.ndarray, predictions: dict[str, np.ndarray],
            age: np.ndarray, sex: np.ndarray, donors: np.ndarray,
            primary_aggregate: dict | None = None, draws: int = 5000,
            quiet: bool = False) -> dict:
    started = time.monotonic()
    total_units = draws * 5
    completed_units = 0

    def tick(stage: str, increment: int = 1) -> None:
        nonlocal completed_units
        completed_units += increment
        if quiet or (completed_units % max(1, draws // 2) and completed_units != total_units):
            return
        elapsed = time.monotonic() - started
        eta = elapsed / completed_units * (total_units - completed_units)
        print(f"E7 progress {completed_units}/{total_units} units "
              f"({completed_units / total_units:.0%}) | elapsed {elapsed / 60:.1f}m "
              f"| ETA ~{eta / 60:.1f}m | {stage}", flush=True)
    y = np.asarray(y, int)
    if set(np.unique(y)) != {0, 1} or len(np.unique(donors)) != len(y):
        raise ValueError("Unique donor and two classes required")
    if set(predictions) != set(MODELS):
        raise ValueError("All five locked models are required")
    for name, p in predictions.items():
        if p.shape != y.shape or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
            raise ValueError(f"Invalid locked {name} predictions")
    counts = primary_counts(y, draws)
    p4, p5 = predictions["M4"], predictions["M5"]
    auc_delta_draws = auc_multiplicities(y, p5, counts) - auc_multiplicities(y, p4, counts)
    primary_ci = np.quantile(auc_delta_draws, [.025, .975])
    if primary_aggregate is not None:
        locked_ci = np.asarray(primary_aggregate["disease"]["primary"]["CI95"], float)
        if not np.allclose(primary_ci, locked_ci, atol=1e-12, rtol=0):
            raise ValueError("Primary bootstrap draws do not reproduce locked CI")
    prevalence = float(y.mean())
    entropy = -prevalence*np.log(prevalence) - (1-prevalence)*np.log(1-prevalence)
    result = {"analysis": "E7", "status": "complete_descriptive",
              "bootstrap_draws": draws, "seed": 20261014,
              "primary_CI_parity": primary_aggregate is not None,
              "null_references": {"AUROC": .5, "AP": prevalence,
                                  "Brier": prevalence*(1-prevalence),
                                  "log_loss": float(entropy)},
              "models": {}, "M5_minus_M4": {}, "calibration": {},
              "reliability": {}}
    metric_draws = {name: {m: np.empty(draws) for m in ("AP", "Brier", "log_loss")}
                    for name in MODELS}
    for name in MODELS:
        result["models"][name] = basic_metrics(y, predictions[name])
    for draw, weights in enumerate(counts):
        ix = np.repeat(np.arange(len(y)), weights)
        for name in MODELS:
            record = basic_metrics(y[ix], predictions[name][ix])
            for metric, value in record.items():
                metric_draws[name][metric][draw] = value
        tick("bootstrap metrics")
    for name in MODELS:
        for metric, values in metric_draws[name].items():
            result["models"][name][metric + "_CI95_percentile"] = np.quantile(values, [.025, .975]).tolist()
    for metric in ("AP", "Brier", "log_loss"):
        difference = metric_draws["M5"][metric] - metric_draws["M4"][metric]
        result["M5_minus_M4"][metric] = {
            "point": result["models"]["M5"][metric] - result["models"]["M4"][metric],
            "CI95_percentile": np.quantile(difference, [.025, .975]).tolist()}
    full_auc = (auc_multiplicities(y, p4, counts)
                - auc_multiplicities(y, predictions["M1"], counts))
    result["M4_minus_M1_full"] = {
        "point": float(roc_auc_score(y, p4) - roc_auc_score(y, predictions["M1"])),
        "CI95_percentile": np.quantile(full_auc, [.025, .975]).tolist()}
    pairs = matched_pairs(age, sex, y, donors)
    result["M4_minus_M1_matched"] = matched_auc_difference(
        y, predictions["M1"], p4, pairs, draws)
    for name in ("M4", "M5"):
        p = predictions[name]
        point = calibration_fit(y, p)
        fits = []
        for weights in counts:
            fits.append(calibration_fit(y, p, weights))
            tick(f"{name} calibration")
        valid = np.array([row for row in fits if row is not None], float).reshape(-1, 2)
        if point is None or len(valid) < .9 * draws:
            result["calibration"][name] = {"status": "not_evaluable_fit_failures",
                                            "successful_draws": len(valid)}
        else:
            result["calibration"][name] = {
                "status": "complete_descriptive", "intercept": point[0], "slope": point[1],
                "intercept_CI95_percentile": np.quantile(valid[:, 0], [.025, .975]).tolist(),
                "slope_CI95_percentile": np.quantile(valid[:, 1], [.025, .975]).tolist(),
                "successful_draws": len(valid)}
        curve, effective = loess_at(p, y, GRID)
        bands = np.full((draws, len(GRID)), np.nan)
        for i, weights in enumerate(counts):
            ix = np.repeat(np.arange(len(y)), weights)
            bands[i], _ = loess_at(p[ix], y[ix], GRID)
            tick(f"{name} reliability")
        rows = []
        for j, grid_value in enumerate(GRID):
            valid_draws = bands[np.isfinite(bands[:, j]), j]
            rows.append({"predicted_probability": float(grid_value),
                         "smoothed_observed": float(curve[j]) if np.isfinite(curve[j]) else None,
                         "effective_n": float(effective[j]),
                         "CI95_pointwise": (np.quantile(valid_draws, [.025, .975]).tolist()
                                            if len(valid_draws) >= .9 * draws else None)})
        result["reliability"][name] = {"span": .5, "grid": rows}
    return result


def synthetic_self_test() -> None:
    rng = np.random.default_rng(7)
    n = 60
    y = np.array([0, 1] * (n // 2))
    predictions = {name: np.clip(.5 + .13*y + rng.normal(0, .12, n), .01, .99)
                   for name in MODELS}
    age = rng.uniform(55, 75, n)
    sex = np.array([0, 1] * (n // 2))
    donors = np.array([f"synthetic-{i}" for i in range(n)])
    result = compute(y, predictions, age, sex, donors, draws=10, quiet=True)
    original, _ = paired_bootstrap(y, predictions["M4"], predictions["M5"],
                                    donors, repetitions=10, seed=20261014)
    reconstructed = (auc_multiplicities(y, predictions["M5"], primary_counts(y, 10))
                     - auc_multiplicities(y, predictions["M4"], primary_counts(y, 10)))
    assert np.allclose(np.quantile(reconstructed, [.025, .975]), original["CI95"])
    assert result["models"]["M4"]["AP"] > 0
    assert result["primary_CI_parity"] is False
    assert len(result["reliability"]["M4"]["grid"]) == len(GRID)
    print("E7_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--locked-predictions", type=Path)
    parser.add_argument("--primary-aggregate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(value is None for value in (args.authorization_record, args.bundle,
        args.binding, args.locked_predictions, args.primary_aggregate, args.output)):
        parser.error("Real E7 requires authorization and all locked inputs")
    require_execution_authorization(args.authorization_record, analysis="E7")
    output = private_output(args.output)
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    _ranks, age, sex, y, donors, _genes = load_locked_input(args.bundle, binding)
    with np.load(args.locked_predictions, allow_pickle=False) as z:
        if not np.array_equal(z["y"], y):
            raise ValueError("Prediction and bound-input row order differ")
        predictions = {name: np.asarray(z[name], float) for name in MODELS}
    primary = json.loads(args.primary_aggregate.read_text(encoding="utf-8"))
    write_aggregate(output, compute(y, predictions, age, sex, donors, primary))
    print("E7_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
