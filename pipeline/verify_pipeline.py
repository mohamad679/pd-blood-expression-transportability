"""Adversarial synthetic checks; no participant data are read by this module."""
import argparse
from pathlib import Path
import tempfile
import time
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

from geo_pipeline import (Data, Prep, array_hash, choose, donor_folds,
                          predict_suite, train_suite, write_json)
from paired_bootstrap import auc_multiplicities, paired_bootstrap


def synthetic(seed=1901, n=120, p=40):
    rng = np.random.default_rng(seed)
    source = np.repeat(['A', 'B'], n // 2)
    y = np.tile(np.repeat([0, 1], n // 4), 2)
    age = rng.normal(60, 8, n)
    sex = rng.binomial(1, .5, n)
    x = rng.normal(size=(n, p)) + (source == 'A')[:, None] * rng.normal(size=p)
    x[:, 0] += .25 * y
    X = (rankdata(x, axis=1) - .5) / p
    return Data(X, age, sex, y, source, np.array(['d'+str(i) for i in range(n)]),
                np.array(['s'+str(i) for i in range(n)]),
                np.array(['g%04d' % i for i in range(p)]),
                [np.arange(i*2, i*2+2) for i in range(6)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    started = time.monotonic(); checks = []
    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks.append(name)
    d = synthetic()
    # Repeated donors, including multiple rows with the same clinical label.
    repeated = d.take(np.repeat(np.arange(len(d.y)), 2))
    folds = donor_folds(repeated, 3, 'synthetic/repeated')
    check('repeated_donors_never_cross_outer_fold', all(
        not set(repeated.donor[a]) & set(repeated.donor[b]) for a, b in folds))
    check('all_rows_held_out_once', sorted(np.concatenate([b for a, b in folds])) == list(range(240)))
    for a, b in folds:
        sub = repeated.take(a)
        check('repeated_donors_never_cross_inner_fold_' + str(len(checks)), all(
            not set(sub.donor[c]) & set(sub.donor[e]) for c, e in donor_folds(sub, 3, 'synthetic/inner')))
    outer_train = np.arange(0, 100)
    # Ensure both sources/classes in training, with a separate unused holdout.
    train, test = d.take(outer_train), d.take(np.arange(100, 120))
    prep = Prep(train, 15)
    changed = synthetic()
    changed.X[100:] = np.random.default_rng(2).normal(1e6, 1e4, (20, 40))
    changed.age[100:] = 1e9; changed.y[100:] = 1 - changed.y[100:]
    check('heldout_extreme_X_age_and_relabel_do_not_change_training_prep',
          prep.fingerprint() == Prep(changed.take(outer_train), 15).fingerprint())
    check('gene_scaler_equals_training_only_mean', np.allclose(
        prep.gene_mean, train.X[:, prep.selected].mean(axis=0)))
    check('marker_scores_are_samplewise_fixed_panel_means', np.allclose(
        prep.baseline_raw(test)[:, 4], test.X[:, test.marker_indices[0]].mean(axis=1)))
    X4, b4 = prep.transform(train, 'M4'); X5, b5 = prep.transform(train, 'M5')
    check('M4_M5_identical_free_baseline', b4 == b5 and np.array_equal(X4, X5[:, :b5]))
    check('candidate_tie_prefers_baseline', choose([
        dict(score=.6, baseline=True, **{'lambda': 0}, rho=0),
        dict(score=.6000005, baseline=False, **{'lambda': .1}, rho=1)])['baseline'])
    # Whole training/tuning refit invariance to a heldout cohort perturbation.
    with tempfile.TemporaryDirectory() as td:
        inner = donor_folds(train, 3, 'synthetic/suite')
        p1, f1, c1 = train_suite(train, inner, Path(td)/'original', 'synthetic',
            max_genes=15, grid=(.1, 1.), rhos=(0., .5), expansion_rounds=0)
        p2, f2, c2 = train_suite(changed.take(outer_train), inner,
            Path(td)/'perturbed', 'synthetic', max_genes=15,
            grid=(.1, 1.), rhos=(0., .5), expansion_rounds=0)
        check('heldout_mutation_does_not_change_tuning_decisions', c1 == c2)
        check('heldout_mutation_does_not_change_fitted_coefficients', all(
            np.array_equal(f1[m].coef_, f2[m].coef_) and f1[m].intercept_ == f2[m].intercept_
            for m in f1))
        # Test the particular inner fold's fit-state against validation extremes.
        a, b = inner[0]
        corrupted = train.take(np.arange(len(train.y)))
        corrupted.X = corrupted.X.copy(); corrupted.age = corrupted.age.copy()
        corrupted.X[b] = 1e9; corrupted.age[b] = -1e9
        check('inner_validation_mutation_does_not_change_its_training_state',
              Prep(train.take(a), 15).fingerprint() == Prep(corrupted.take(a), 15).fingerprint())
        # Explicit baseline branch must use exactly saved baseline predictions.
        c1['M5']['selected']['reuses_baseline'] = 'M4'
        f1['M5'] = f1['M4']
        pred = predict_suite(p1, f1, c1, test)
        check('M5_baseline_branch_bit_identical_to_M4', np.array_equal(pred['M5'], pred['M4']))
    rng = np.random.default_rng(29)
    y = np.repeat([0, 1], 30)
    s4, s5 = rng.integers(0, 8, (2, 60))  # many ties, two separate rankings
    counts = rng.multinomial(60, np.ones(60)/60, size=80)
    check('bootstrap_multiplicities_equal_independent_sklearn_with_ties', np.allclose(
        auc_multiplicities(y, s4, counts),
        [roc_auc_score(y, s4, sample_weight=w) for w in counts], atol=1e-14))
    ids = np.arange(60)
    r, delta = paired_bootstrap(y, s4, s5, ids)
    reverse, rev = paired_bootstrap(y, s5, s4, ids)
    identical, zero = paired_bootstrap(y, s4, s4, ids)
    check('paired_bootstrap_swapping_models_exactly_reverses_difference', np.array_equal(delta, -rev))
    check('identical_scores_have_zero_resampling_difference', np.array_equal(zero, np.zeros(5000)))
    rejected = False
    try:
        paired_bootstrap(y, s4, s5, np.zeros(60))
    except ValueError:
        rejected = True
    check('bootstrap_rejects_repeated_visits_as_independent_donors', rejected)
    # Synthetic high-dimensional preprocessing without an expensive training fit.
    high = synthetic(p=1000)
    ph = Prep(high, 1000)
    check('1000_gene_design_finite_and_train_scaled', np.isfinite(ph.transform(high, 'M5')[0]).all())
    write_json(args.output, dict(status='passed', checks=len(checks), names=checks,
        participant_data_read=False, real_training_performed=False,
        seconds=time.monotonic()-started, bootstrap_repetitions=5000,
        bootstrap_coverage_validation='not_yet_run; algebra_and_resampling_checks_only',
        numerical_1000_gene_fit_validation='existing_v3.1_32_fit_stress_not_repeated_here'))
    print('PASSED', len(checks), 'synthetic pipeline checks', flush=True)


if __name__ == '__main__':
    main()
