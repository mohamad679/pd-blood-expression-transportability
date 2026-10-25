"""Synthetic coverage validation for the registered paired donor-bootstrap CI.

This evaluates fixed-score sampling uncertainty only. It never reads participant
files or refits models. The provisional planned class sizes are used solely as
synthetic sample sizes; scenario truths are known from their generating laws.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri

from paired_bootstrap import auc_multiplicities, paired_bootstrap


N_CONTROL = 161
N_PD = 311
INNER_BOOTSTRAPS = 5000
BASE_SEED = 20261019


def continuous_auc(case_shift: float) -> float:
    """AUC for equally dispersed unit-normal control and case scores."""
    return float(ndtr(case_shift / math.sqrt(2.0)))


def ordinal_auc(case_shift: float, cutpoints=(-0.5, 0.5)) -> float:
    """Exact population AUC for normal scores rounded into three tied levels."""
    c0 = np.array([ndtr(cutpoints[0]), ndtr(cutpoints[1]) - ndtr(cutpoints[0]),
                   1.0 - ndtr(cutpoints[1])])
    q = np.array([ndtr(cutpoints[0] - case_shift),
                  ndtr(cutpoints[1] - case_shift) - ndtr(cutpoints[0] - case_shift),
                  1.0 - ndtr(cutpoints[1] - case_shift)])
    return float(q[2] * (c0[0] + c0[1]) + q[1] * c0[0]
                 + 0.5 * np.dot(q, c0))


SCENARIOS = {
    "continuous_null_distinct": {
        "shift_m4": math.sqrt(2.0) * ndtri(0.70),
        "shift_m5": math.sqrt(2.0) * ndtri(0.70),
        "correlation": 0.8,
        "ties": False,
        "identical": False,
        "purpose": "distinct continuous model scores with true delta zero",
    },
    "continuous_delta_005": {
        "shift_m4": math.sqrt(2.0) * ndtri(0.70),
        "shift_m5": math.sqrt(2.0) * ndtri(0.75),
        "correlation": 0.8,
        "ties": False,
        "identical": False,
        "purpose": "distinct continuous scores with a true AUROC gain of 0.05",
    },
    "tied_null_distinct": {
        "shift_m4": math.sqrt(2.0) * ndtri(0.70),
        "shift_m5": math.sqrt(2.0) * ndtri(0.70),
        "correlation": 0.8,
        "ties": True,
        "identical": False,
        "purpose": "distinct, three-level tied model scores with true delta zero",
    },
    "identical_model": {
        "shift_m4": math.sqrt(2.0) * ndtri(0.70),
        "shift_m5": math.sqrt(2.0) * ndtri(0.70),
        "correlation": 1.0,
        "ties": False,
        "identical": True,
        "purpose": "identical-model branch; degenerate zero-width interval",
    },
}


def generate_sample(rng: np.random.Generator, scenario: dict) -> tuple[np.ndarray, ...]:
    n = N_CONTROL + N_PD
    y = np.r_[np.zeros(N_CONTROL, dtype=np.int8), np.ones(N_PD, dtype=np.int8)]
    shared = rng.standard_normal(n)
    if scenario["identical"]:
        latent4 = shared
        latent5 = shared
    else:
        independent = rng.standard_normal(n)
        rho = scenario["correlation"]
        latent4 = shared
        latent5 = rho * shared + math.sqrt(1.0 - rho * rho) * independent
    s4 = latent4 + y * scenario["shift_m4"]
    s5 = latent5 + y * scenario["shift_m5"]
    if scenario["ties"]:
        s4 = np.digitize(s4, [-0.5, 0.5]).astype(float)
        s5 = np.digitize(s5, [-0.5, 0.5]).astype(float)
    donors = np.arange(n)
    return y, s4, s5, donors


def true_delta(scenario: dict) -> float:
    auc = ordinal_auc if scenario["ties"] else continuous_auc
    return auc(scenario["shift_m5"]) - auc(scenario["shift_m4"])


def auc_multiplicities_fast(y: np.ndarray, scores: np.ndarray,
                            counts: np.ndarray) -> np.ndarray:
    """Same tie-aware AUC as auc_multiplicities, grouping sorted score ties."""
    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    starts = np.r_[0, np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]) + 1]
    ordered_counts = counts[:, order]
    pos = np.add.reduceat(ordered_counts * y[order], starts, axis=1)
    neg = np.add.reduceat(ordered_counts * (1 - y[order]), starts, axis=1)
    denominator = pos.sum(axis=1) * neg.sum(axis=1)
    if np.any(denominator == 0):
        raise ValueError("Resample missing a class")
    return (pos * (np.cumsum(neg, axis=1) - 0.5 * neg)).sum(axis=1) / denominator


def fast_paired_bootstrap(y: np.ndarray, p4: np.ndarray, p5: np.ndarray,
                          repetitions: int = INNER_BOOTSTRAPS,
                          seed: int = 20261014) -> tuple[float, list[float], np.ndarray]:
    """Vectorized equivalent of the registered paired class-stratified CI."""
    rng = np.random.default_rng(seed)
    classes = [np.flatnonzero(y == c) for c in (0, 1)]
    draws = []
    for offset in range(0, repetitions, 100):
        b = min(100, repetitions - offset)
        counts = np.zeros((b, len(y)), dtype=np.int64)
        for ix in classes:
            counts[:, ix] = rng.multinomial(len(ix), np.full(len(ix), 1 / len(ix)), size=b)
        draws.extend(auc_multiplicities_fast(y, p5, counts)
                     - auc_multiplicities_fast(y, p4, counts))
    draws = np.asarray(draws)
    point = float(auc_multiplicities_fast(y, p5, np.ones((1, len(y))))[0]
                  - auc_multiplicities_fast(y, p4, np.ones((1, len(y))))[0])
    return point, np.quantile(draws, [0.025, 0.975]).tolist(), draws


def verify_fast_calculation() -> dict:
    """Compare 5,000 exact resamples to the public implementation on tied and untied data."""
    rng = np.random.default_rng(20261018)
    y = np.r_[np.zeros(20, dtype=np.int8), np.ones(30, dtype=np.int8)]
    x4 = rng.normal(size=len(y)) + y * 0.7
    x5 = 0.8 * x4 + 0.6 * rng.normal(size=len(y)) + y * 0.2
    tied = (np.digitize(x4, [-0.5, 0.5]).astype(float),
            np.digitize(x5, [-0.5, 0.5]).astype(float))
    results = {}
    for name, a, b in (("untied", x4, x5), ("tied", *tied)):
        reference, ref_draws = paired_bootstrap(
            y, a, b, np.arange(len(y)), repetitions=INNER_BOOTSTRAPS, seed=20261014
        )
        point, interval, draws = fast_paired_bootstrap(y, a, b)
        error = float(np.max(np.abs(ref_draws - draws)))
        if error > 5e-15 or not np.allclose(reference["CI95"], interval, atol=5e-15, rtol=0):
            raise AssertionError(f"Vectorized CI disagrees with public implementation ({name})")
        results[name] = {"max_resample_difference": error,
                         "CI95_matches": True,
                         "point_matches": math.isclose(reference["delta_AUROC"], point,
                                                       abs_tol=5e-15, rel_tol=0)}
    return results


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> list[float]:
    p = successes / trials
    denom = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / denom
    half = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials)) / denom
    return [max(0.0, center - half), min(1.0, center + half)]


def run(outer_repetitions: int) -> dict:
    output = {
        "status": "complete",
        "purpose": "synthetic fixed-score confidence-interval coverage; not model-training or PPMI analysis",
        "sample_sizes": {"control": N_CONTROL, "case": N_PD},
        "inner_bootstrap_repetitions": INNER_BOOTSTRAPS,
        "outer_sampling_repetitions_per_scenario": outer_repetitions,
        "outer_seed": BASE_SEED,
        "bootstrap_seed_per_dataset": 20261014,
        "interval": "paired, class-stratified donor percentile CI, 95%, fixed fitted scores",
        "vectorized_algorithm_verified_against_public_implementation": verify_fast_calculation(),
        "scenarios": {},
    }
    for index, (name, scenario) in enumerate(SCENARIOS.items()):
        rng = np.random.default_rng(BASE_SEED + index)
        target = true_delta(scenario)
        covered = 0
        zero_width = 0
        widths = []
        points = []
        for _ in range(outer_repetitions):
            y, s4, s5, donors = generate_sample(rng, scenario)
            point, interval, _ = fast_paired_bootstrap(
                y, s4, s5, repetitions=INNER_BOOTSTRAPS, seed=20261014
            )
            low, high = interval
            covered += int(low <= target <= high)
            zero_width += int(low == high)
            widths.append(high - low)
            points.append(point)
        output["scenarios"][name] = {
            "purpose": scenario["purpose"],
            "true_delta_AUROC": target,
            "empirical_coverage": covered / outer_repetitions,
            "coverage_monte_carlo_95pct_Wilson_CI": wilson_interval(covered, outer_repetitions),
            "covered_intervals": covered,
            "outer_repetitions": outer_repetitions,
            "mean_interval_width": float(np.mean(widths)),
            "zero_width_intervals": zero_width,
            "mean_estimated_delta": float(np.mean(points)),
        }
        print(f"Finished {name}: coverage={covered / outer_repetitions:.3f}", flush=True)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-repetitions", type=int, default=200)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.outer_repetitions < 1:
        parser.error("--outer-repetitions must be positive")
    report = run(args.outer_repetitions)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for name, result in report["scenarios"].items():
        print(f"{name}: coverage={result['empirical_coverage']:.3f} "
              f"MC95%CI={result['coverage_monte_carlo_95pct_Wilson_CI']} "
              f"true_delta={result['true_delta_AUROC']:.4f}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
