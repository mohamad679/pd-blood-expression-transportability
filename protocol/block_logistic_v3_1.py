"""v3.1 weighted partial-penalty logistic model; no data loading or tuning."""

import numpy as np
from glum import GeneralizedLinearRegressor


def balanced_source_weights(y, cohort):
    """Equal total source mass and equal class mass within each training source.

    For K observed sources, w_i = n/(2*K*n[source_i,class_i]). Recompute using
    only the current training partition. With K=1, this balances the two classes.
    """
    y, cohort = np.asarray(y), np.asarray(cohort)
    if y.ndim != 1 or cohort.shape != y.shape or set(np.unique(y)) != {0, 1}:
        raise ValueError("Aligned binary labels and source vector required")
    sources = np.unique(cohort)
    weights = np.empty(len(y), dtype=float)
    for source in sources:
        for label in (0, 1):
            mask = (cohort == source) & (y == label)
            count = int(mask.sum())
            if count == 0:
                raise ValueError("Each training source must contain both classes")
            weights[mask] = len(y) / (2 * len(sources) * count)
    return weights


def kish_effective_n(weights):
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("Finite positive weights required")
    return float(w.sum() ** 2 / np.dot(w, w))


def fit_block_logistic(X, y, n_baseline, gene_alpha=0.01, gene_l1_ratio=0.5,
                       sample_weight=None):
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if X.ndim != 2 or not 0 <= n_baseline <= X.shape[1]:
        raise ValueError("Invalid baseline block size")
    if y.shape != (len(X),) or set(np.unique(y)) != {0.0, 1.0}:
        raise ValueError("Binary nonconstant outcome required")
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite model inputs")
    if gene_alpha < 0 or not 0 <= gene_l1_ratio <= 1:
        raise ValueError("Invalid gene penalty")
    if sample_weight is not None:
        sample_weight = np.asarray(sample_weight, dtype=float)
        if sample_weight.shape != y.shape or not np.isfinite(sample_weight).all() or np.any(sample_weight <= 0):
            raise ValueError("Finite positive sample weights aligned with y required")
    penalty = np.r_[np.zeros(n_baseline), np.ones(X.shape[1] - n_baseline)]
    model = GeneralizedLinearRegressor(
        family="binomial", link="logit", fit_intercept=True,
        alpha=gene_alpha, l1_ratio=gene_l1_ratio, P1=penalty, P2=penalty,
        solver="irls-cd", max_iter=1000, max_inner_iter=100000,
        gradient_tol=1e-7, step_size_tol=None, scale_predictors=False,
        selection="cyclic", random_state=20261014,
    )
    model.fit(X, y, sample_weight=sample_weight)
    return model
