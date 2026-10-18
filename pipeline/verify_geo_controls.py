"""Synthetic checks for auxiliary controls and clipping attribution; no study input."""
import argparse
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy.special import expit
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score

import geo_pipeline as geo
from geo_controls import FeaturePrep, bootstrap_control, control_folds
from audit_sigmoid import probability_diagnostic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20261018)
    n, p = 80, 45
    sex = np.tile([0, 1], n//2)
    age = rng.uniform(40, 80, n)
    X = rng.normal(size=(n, p))
    X[:, :5] += sex[:, None]*5
    X[:, 5:10] += age[:, None]/3
    data = geo.Data(X, age, sex, np.zeros(n, int), np.repeat('synthetic', n),
        np.array(['donor'+str(i) for i in range(n)]),
        np.array(['sample'+str(i) for i in range(n)]),
        np.array(['gene%03d' % i for i in range(p)]), [])
    a, b = data.take(np.arange(60)), data.take(np.arange(60, 80))
    prep = FeaturePrep(a, np.arange(5, p), 10)
    fingerprint = prep.fingerprint()
    test = b.take(np.arange(len(b.y)))
    test.X = test.X + 1000
    prep.transform(test)
    assert fingerprint == prep.fingerprint()
    assert np.max(abs(prep.transform(a).mean(axis=0))) < 1e-12
    assert np.allclose(prep.transform(a).std(axis=0), 1)
    # All disease labels are zero. Sex folds still contain both sex classes.
    folds = control_folds(data, 'sex')
    for ti, vi in folds:
        assert set(data.donor[ti]).isdisjoint(data.donor[vi])
        assert len(np.unique(data.sex[vi])) == 2
    assert all(set(data.donor[ti]).isdisjoint(data.donor[vi]) for ti,vi in control_folds(data, 'age'))
    sx = FeaturePrep(a, np.arange(5))
    model, records = geo.fit_checked(sx.transform(a), a.sex, 0, .01, 0.,
        sample_weight=np.ones(len(a.y)), context={'test': 'synthetic_sex'})
    assert records[-1]['passed']
    assert roc_auc_score(b.sex, model.predict(sx.transform(b))) > .99
    ag = FeaturePrep(a, np.arange(5, 10), 5)
    model_age = Ridge(alpha=.01, solver='svd').fit(ag.transform(a), a.age)
    assert np.mean(abs(model_age.predict(ag.transform(b))-b.age)) < 3
    linear = np.array([-40., -35., 35., 40.])
    y = np.array([0, 1, 0, 1])
    eps = 50*np.finfo(float).eps
    diagnostic = probability_diagnostic(y, linear, np.clip(expit(linear), eps, 1-eps), np.repeat('s',4))
    assert diagnostic['mixed_class_tail_pairs'] == 2
    assert diagnostic['AUROC']['linear'] != diagnostic['AUROC']['clipped']
    uncertainty = bootstrap_control('age', np.arange(10.), np.arange(10.), 4.5, 'synthetic/age', 100)
    assert uncertainty['percentile_95_CI']['MAE'] == [0., 0.]
    assert uncertainty['percentile_95_CI']['MAE_minus_median_baseline'][1] < 0
    result = dict(study_data_loaded=False, passed=True, checks=[
        'test_data_cannot_change_preprocessing', 'training_standardization',
        'donor_disjoint_sex_and_age_folds', 'sex_target_not_disease_label',
        'known_sex_signal_and_KKT', 'known_autosomal_age_signal',
        'mixed_class_tail_tie_detection', 'paired_age_bootstrap_difference'],
        synthetic_clipping_example=diagnostic)
    geo.write_json(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
