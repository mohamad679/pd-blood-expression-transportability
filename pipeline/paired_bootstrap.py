"""Fixed-score, case/control-stratified paired donor AUROC bootstrap.

Primary baseline population must contain exactly one prediction per donor.
Repeated visits must be collapsed by the locked eligibility rule before calling.
This function does not represent training or model-selection uncertainty.
"""
import numpy as np
from scipy.stats import rankdata


def auc_multiplicities(y, scores, counts):
    """AUROC for many resamples without materializing duplicated rows; exact ties."""
    y, scores, counts = np.asarray(y), np.asarray(scores), np.asarray(counts)
    _, group = np.unique(scores, return_inverse=True)
    membership = np.eye(group.max() + 1)[group]
    pos = (counts * y) @ membership
    neg = (counts * (1 - y)) @ membership
    denominator = pos.sum(axis=1) * neg.sum(axis=1)
    if np.any(denominator == 0):
        raise ValueError('Resample missing a class')
    return (pos * (np.cumsum(neg, axis=1) - .5 * neg)).sum(axis=1) / denominator


def paired_bootstrap(y, p4, p5, donors, repetitions=5000, seed=20261014):
    y, p4, p5, donors = [np.asarray(a) for a in (y, p4, p5, donors)]
    if (y.ndim != 1 or any(a.shape != y.shape for a in (p4, p5, donors))
            or set(np.unique(y)) != {0, 1} or len(np.unique(donors)) != len(y)
            or not np.isfinite(p4).all() or not np.isfinite(p5).all()):
        raise ValueError('One finite paired prediction per unique binary-outcome donor required')
    rng = np.random.default_rng(seed)
    classes = [np.flatnonzero(y == c) for c in (0, 1)]
    draws = []
    for offset in range(0, repetitions, 100):
        b = min(100, repetitions - offset)
        counts = np.zeros((b, len(y)))
        for ix in classes:
            counts[:, ix] = rng.multinomial(len(ix), np.full(len(ix), 1 / len(ix)), size=b)
        draws.extend(auc_multiplicities(y, p5, counts) - auc_multiplicities(y, p4, counts))
    draws = np.asarray(draws)
    point = float(auc_multiplicities(y, p5, np.ones((1, len(y))))[0]
                  - auc_multiplicities(y, p4, np.ones((1, len(y))))[0])
    return dict(delta_AUROC=point, CI95=np.quantile(draws, [.025, .975]).tolist(),
                repetitions=repetitions, seed=seed,
                conditional_on_fixed_fitted_models=True), draws
