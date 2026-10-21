#!/usr/bin/env python3
"""E6: donor/bootstrap workflow refits on GEO, scored on fixed PPMI."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import geo_pipeline as geo
from e03_ppmi_nested_cv import fit_pair, score_model
from e04_transductive_universe import ppmi_data, save_fit
from postlock_common import (private_output, require_execution_authorization,
                             rng_for, sha256, synthetic_ranks, write_aggregate)


def stratified_donor_bootstrap(data: geo.Data, repetition: int) -> geo.Data:
    if len(np.unique(data.donor)) != len(data.y):
        raise ValueError("E6 requires one GEO record per donor")
    rng = rng_for(6, repetition)
    indices = []
    for source in sorted(set(map(str, data.source))):
        for label in (0, 1):
            stratum = np.flatnonzero((data.source == source) & (data.y == label))
            if len(stratum) < 5:
                raise ValueError("Source/class stratum has fewer than five donors")
            indices.extend(rng.choice(stratum, len(stratum), replace=True).tolist())
    result = data.take(np.asarray(indices, int))
    if len(result.y) != len(data.y):
        raise AssertionError("Bootstrap failed to preserve stratum sizes")
    return result


def attempt(data: geo.Data, target: geo.Data, output: Path, repetition: int,
            inner_folds: int = 5, max_genes: int = 1000,
            grid=(.0001, .001, .01, .1, 1.), rhos=(0., .5, 1.),
            expansion_rounds: int = 2) -> dict:
    train = stratified_donor_bootstrap(data, repetition)
    scope = f"E6/bootstrap/{repetition}"
    inner = geo.donor_folds(train, inner_folds, scope + "/inner")
    fit_dir = output / "fits"
    prep, m4, m5, decision = fit_pair(train, inner, fit_dir, scope,
        max_genes=max_genes, grid=grid, rhos=rhos,
        expansion_rounds=expansion_rounds)
    save_fit(fit_dir, prep, m4, m5, decision)
    p4 = score_model(prep, m4, target, "M4")
    p5 = p4 if m5 is None else score_model(prep, m5, target, "M5")
    return {"repetition": repetition, "status": "complete",
            "M4_AUROC": float(roc_auc_score(target.y, p4)),
            "M5_AUROC": float(roc_auc_score(target.y, p5)),
            "M5_minus_M4_AUROC": float(roc_auc_score(target.y, p5)
                                        - roc_auc_score(target.y, p4)),
            "M5_selection": decision}


def run(data: geo.Data, target: geo.Data, output: Path,
        max_attempts: int = 200, switch_after: int = 20,
        switch_seconds: float = 6*3600,
        inner_folds: int = 5, max_genes: int = 1000,
        grid=(.0001, .001, .01, .1, 1.), rhos=(0., .5, 1.),
        expansion_rounds: int = 2, synthetic: bool = False) -> dict:
    if output.exists():
        raise ValueError("Use a fresh E6 output directory; never overwrite attempts")
    output.mkdir(parents=True)
    start = time.monotonic()
    results = []
    cap = max_attempts
    for repetition in range(max_attempts):
        if repetition >= cap:
            break
        folder = output / f"attempt_{repetition:03d}"
        folder.mkdir()
        attempt_start = time.monotonic()
        try:
            result = attempt(data, target, folder, repetition, inner_folds,
                             max_genes, grid, rhos, expansion_rounds)
        except Exception as error:
            # Keep the attempt and fit logs; do not substitute another draw.
            result = {"repetition": repetition, "status": "failed",
                      "error_type": type(error).__name__}
        result["elapsed_seconds"] = time.monotonic() - attempt_start
        (folder / "attempt_aggregate.json").write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False), encoding="utf-8")
        results.append(result)
        if repetition + 1 == switch_after and max_attempts == 200:
            cap = 100 if time.monotonic() - start > switch_seconds else 200
        if (repetition + 1) % 5 == 0 or repetition + 1 == cap:
            elapsed = time.monotonic() - start
            eta = elapsed / (repetition + 1) * (cap - repetition - 1)
            print(f"E6 progress {repetition + 1}/{cap} attempts | "
                  f"elapsed {elapsed/3600:.1f}h | ETA ~{eta/3600:.1f}h", flush=True)
    successful = [r["M5_minus_M4_AUROC"] for r in results if r["status"] == "complete"]
    report = {"analysis": "E6", "exploratory": True,
              "attempts": len(results), "successful_attempts": len(successful),
              "failed_attempts": len(results) - len(successful),
              "runtime_threshold_seconds": switch_seconds,
              "first_20_attempt_runtime_seconds": sum(r["elapsed_seconds"] for r in results[:20]),
              "cap_after_runtime_rule": cap,
              "source_class_stratified_donor_bootstrap": True,
              "fixed_PPMI_population": True,
              "status": ("complete_exploratory" if len(successful) >= .9 * cap
                         else "not_evaluable_excess_fit_failures")}
    if report["status"] == "complete_exploratory":
        report["M5_minus_M4_AUROC_mean"] = float(np.mean(successful))
        report["M5_minus_M4_AUROC_CI95_percentile"] = np.quantile(
            successful, [.025, .975]).tolist()
    if synthetic:
        (output / "E6_AGGREGATE.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    else:
        write_aggregate(output / "E6_AGGREGATE.json", report)
    return report


def synthetic_self_test() -> None:
    rng = np.random.default_rng(61)
    n, g = 48, 40
    ranks = synthetic_ranks(n, g, 61)
    data = geo.Data(ranks, rng.uniform(50, 75, n), rng.integers(0, 2, n),
                    np.array([0, 1] * (n//2)), np.repeat("SYNTHETIC", n),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"s{i}" for i in range(n)]),
                    np.array([f"G{i:04d}" for i in range(g)]),
                    [np.arange(i*3, i*3+3) for i in range(6)])
    boot = stratified_donor_bootstrap(data, 0)
    assert len(boot.y) == n and int(boot.y.sum()) == int(data.y.sum())
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        report = run(data, data, Path(temporary) / "run", max_attempts=2,
            switch_after=1, switch_seconds=3600, inner_folds=2,
            max_genes=8, grid=(.1,), rhos=(0.,), expansion_rounds=0,
            synthetic=True)
        assert report["attempts"] == 2
        assert report["successful_attempts"] == 2
    print("E6_SYNTHETIC_TEST_PASSED")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-self-test", action="store_true")
    parser.add_argument("--authorization-record", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.synthetic_self_test:
        synthetic_self_test()
        return
    if any(v is None for v in (args.authorization_record, args.bundle,
                                args.binding, args.output_dir)):
        parser.error("Real E6 requires authorization, bundle, binding and output")
    require_execution_authorization(args.authorization_record, analysis="E6")
    os.umask(0o077)
    output = private_output(args.output_dir)
    data, _ = geo.load_geo()
    target = ppmi_data(args.bundle, args.binding)
    if not np.array_equal(data.genes, target.genes):
        raise ValueError("GEO/PPMI gene universe differs")
    result = run(data, target, output)
    print(f"E6_FINISHED: {result['status']}; aggregate output under ~/PPMI_private")


if __name__ == "__main__":
    main()
