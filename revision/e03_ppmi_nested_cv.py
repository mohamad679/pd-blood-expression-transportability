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
from collections import Counter

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

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
             expansion_rounds: int = 2):
    """Reproduce v3.1 baseline-first tuning for only M4 and M5."""
    output.mkdir(parents=True, exist_ok=False)
    fit_log = output / "fits.jsonl"
    warning_log = output / "warnings.jsonl"
    fold_data = []
    baseline_scores = []
    for fold, (ti, vi) in enumerate(inner):
        a, b = train.take(ti), train.take(vi)
        if set(a.donor) & set(b.donor):
            raise AssertionError("Inner donor leakage")
        prep = geo.Prep(a, max_genes=max_genes)
        m4, _ = checked_fit_with_warning_log(prep, a, "M4", 0., 0.,
                                dict(scope=scope, stage="inner", fold=fold), fit_log,
                                warning_log)
        auc = float(roc_auc_score(b.y, score_model(prep, m4, b, "M4")))
        baseline_scores.append(auc)
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
                    scores.append(float(roc_auc_score(b.y, score_model(prep, fit, b, "M5"))))
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
    prep = geo.Prep(train, max_genes=max_genes)
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


def run_nested(data: geo.Data, output: Path, repetitions: int = 20,
               outer_folds: int = 5, inner_folds: int = 5,
               max_genes: int = 1000, grid=GRID, rhos=RHOS,
               expansion_rounds: int = 2) -> dict:
    output = Path(output)
    if output.exists():
        raise ValueError("Use a fresh E3 output directory")
    output.mkdir(parents=True)
    started = time.monotonic()
    repetition_means = []
    boundary_count = 0
    for rep in range(repetitions):
        seed = int(rng_for(3, rep).integers(0, 2**32 - 1))
        paired = []
        for fold, (ti, vi) in enumerate(donor_folds(data, outer_folds, seed)):
            train, valid = data.take(ti), data.take(vi)
            inner_seed = int(rng_for(3, rep * 100 + fold + 1).integers(0, 2**32 - 1))
            inner = list(donor_folds(train, inner_folds, inner_seed))
            scope = f"E3/repetition/{rep}/outer/{fold}"
            prep, m4, m5, decision = fit_pair(train, inner,
                output / f"rep_{rep:02d}_fold_{fold}", scope, max_genes,
                grid, rhos, expansion_rounds)
            p4 = score_model(prep, m4, valid, "M4")
            p5 = p4 if m5 is None else score_model(prep, m5, valid, "M5")
            paired.append(float(roc_auc_score(valid.y, p5) - roc_auc_score(valid.y, p4)))
            boundary_count += int(decision["persistent_boundary"])
        repetition_means.append(float(np.mean(paired)))
        elapsed = time.monotonic() - started
        eta = elapsed / (rep + 1) * (repetitions - rep - 1)
        print(f"E3 progress {rep + 1}/{repetitions} repetitions; "
              f"elapsed {elapsed / 3600:.1f}h; ETA ~{eta / 3600:.1f}h", flush=True)
    mean = float(np.mean(repetition_means))
    interpretation = ("within_RNA_signal_not_transportable" if mean >= .05 else
                      "no_meaningful_increment_detected" if -.02 <= mean <= .02 else
                      "inconclusive_small_positive" if .02 < mean < .05 else
                      "adverse_within_platform_increment")
    return {"analysis": "E3", "status": "complete_exploratory",
            "seed_rule": "SeedSequence([20261014, 3, repetition])",
            "outer_repetitions": repetitions, "outer_folds": outer_folds,
            "inner_folds": inner_folds, "mean_paired_fold_delta_by_repetition": repetition_means,
            "mean_delta_AUROC": mean, "repetition_summary": summary(np.array(repetition_means)),
            "range": [float(min(repetition_means)), float(max(repetition_means))],
            "persistent_boundary_count": boundary_count,
            "prediction_warning_counts": dict(PREDICTION_WARNINGS),
            "interpretation_band": interpretation}


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
        result = run_nested(data, Path(temporary) / "nested", repetitions=1,
            outer_folds=2, inner_folds=2, max_genes=8, grid=(.1,),
            rhos=(0.,), expansion_rounds=0)
        assert result["status"] == "complete_exploratory"
        assert len(result["mean_paired_fold_delta_by_repetition"]) == 1
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
    manifest = json.loads((ROOT / "protocol/v3_marker_manifest.json").read_text(encoding="utf-8"))
    lookup = {str(g): i for i, g in enumerate(genes)}
    markers = [np.array([lookup[g] for g in manifest["cell_marker_panels"][name]["gene_ids"]], int)
               for name in geo.PANELS]
    data = geo.Data(ranks, age, sex, y, np.repeat("PPMI", len(y)), donors,
                    donors.copy(), genes, markers)
    result = run_nested(data, output)
    write_aggregate(output / "E3_AGGREGATE.json", result)
    print("E3_COMPLETE: aggregate output written under ~/PPMI_private")


if __name__ == "__main__":
    main()
