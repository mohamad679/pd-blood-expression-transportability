"""Normalized weighted binomial objective and audited numerical fits, synthetic first.

P2 is a diagonal vector. X is the actual caller-standardized design. No data
loading, feature selection, target-cohort statistics, or tuning occurs here.
"""

import json
from pathlib import Path
import time
import warnings

import numpy as np
from scipy.special import expit
from glum import GeneralizedLinearRegressor


KKT_LIMIT = 1e-6
INITIAL = dict(max_iter=1000, max_inner_iter=100000, gradient_tol=1e-7)
RETRY = dict(max_iter=5000, max_inner_iter=1000000, gradient_tol=1e-8)
REFERENCE_GLM = dict(max_iter=5000, max_inner_iter=1000000, gradient_tol=1e-9)


def objective_and_kkt(X, y, weights, parameters, lam, rho, p1, p2,
                      loss_gradient_multiplier=1.0):
    """Return true objective and KKT in original input coefficient coordinates.

    loss_gradient_multiplier is exclusively for the scale-error negative check;
    it does not change the objective or any fit.
    """
    X, y = np.asarray(X, float), np.asarray(y, float)
    weights = np.ones(len(y)) if weights is None else np.asarray(weights, float)
    theta, p1, p2 = (np.asarray(v, float) for v in (parameters, p1, p2))
    if (X.ndim != 2 or y.shape != (len(X),) or weights.shape != y.shape
            or theta.shape != (X.shape[1]+1,) or p1.shape != (X.shape[1],)
            or p2.shape != p1.shape or np.any(weights <= 0)
            or np.any(p1 < 0) or np.any(p2 < 0) or lam < 0 or not 0 <= rho <= 1):
        raise ValueError('Invalid normalized binomial objective inputs')
    if not all(np.isfinite(v).all() for v in (X, y, weights, theta, p1, p2)):
        raise ValueError('Nonfinite objective inputs')
    a, beta = theta[0], theta[1:]
    normalized = weights / weights.sum()
    eta = a + X @ beta
    probability = expit(eta)
    objective = (np.dot(normalized, np.logaddexp(0, eta)-y*eta)
                 + lam*rho*np.dot(p1, np.abs(beta))
                 + 0.5*lam*(1-rho)*np.dot(p2, beta**2))
    errors = normalized * (probability-y) * loss_gradient_multiplier
    gradient = X.T @ errors
    smooth = gradient + lam*(1-rho)*p2*beta
    residual = smooth + lam*rho*p1*np.sign(beta)
    zero = beta == 0
    residual[zero] = np.sign(smooth[zero])*np.maximum(
        np.abs(smooth[zero])-lam*rho*p1[zero], 0)
    residual = np.r_[errors.sum(), residual]
    return {
        'objective': float(objective),
        'KKT_max': float(np.max(np.abs(residual))),
        'KKT_L1': float(np.sum(np.abs(residual))),
        'intercept_KKT': float(abs(residual[0])),
        'baseline_KKT_max': float(np.max(np.abs(residual[1:][(p1 == 0) & (p2 == 0)]), initial=0)),
        'nonzero_penalized_coefficients': int(np.sum((beta != 0) & ((p1 != 0) | (p2 != 0)))),
        'finite': bool(np.isfinite(objective) and np.isfinite(residual).all()),
    }


def append_record(log_path, record):
    if log_path is not None:
        with Path(log_path).open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(record, ensure_ascii=False, allow_nan=False)+'\n')


def solve_once(X, y, weights, lam, rho, p1, p2, settings, context,
               attempt, log_path=None):
    """One fully recorded fit; an inner warning alone is not an outer failure."""
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=float)
    weights = None if weights is None else np.asarray(weights, dtype=float)
    model = None
    error = None
    diagnostics = None
    start = time.monotonic()
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        try:
            # Documented numerical amendment: direct IRLS least squares for pure
            # L2, same objective/penalties. CD remains for L1 and baseline-only.
            solver = 'irls-ls' if rho == 0 and lam > 0 else 'irls-cd'
            model = GeneralizedLinearRegressor(
                family='binomial', link='logit', fit_intercept=True,
                alpha=lam, l1_ratio=rho, P1=np.asarray(p1), P2=np.asarray(p2),
                solver=solver, step_size_tol=None, scale_predictors=False,
                selection='cyclic', random_state=20261014, **settings)
            model.fit(X, y, sample_weight=weights)
            theta = np.r_[model.intercept_, model.coef_]
            diagnostics = objective_and_kkt(X, y, weights, theta, lam, rho, p1, p2)
        except (ArithmeticError, ValueError, RuntimeError) as exception:
            error = dict(type=type(exception).__name__, message=str(exception))
    passed = bool(diagnostics and diagnostics['finite'] and diagnostics['KKT_max'] <= KKT_LIMIT)
    record = {
        **context, 'lambda': lam, 'rho': rho, 'attempt': attempt,
        'n': len(y), 'p': X.shape[1], 'solver_settings': settings,
        'solver': solver, 'scale_predictors': False, 'sample_weight': 'uniform' if weights is None else 'explicit',
        'weight_sum': float(len(y) if weights is None else weights.sum()),
        'warnings': [dict(category=type(w.message).__name__, message=str(w.message)) for w in captured],
        'exception': error, 'diagnostics': diagnostics, 'KKT_limit': KKT_LIMIT,
        'passed': passed, 'n_iter': int(model.n_iter_) if model is not None and hasattr(model, 'n_iter_') else None,
        'seconds': round(time.monotonic()-start, 4),
    }
    append_record(log_path, record)
    return model, record


def fit_checked(X, y, n_baseline, lam, rho, sample_weight, context, log_path=None):
    """Fixed initial fit, one cold retry iff final KKT fails, never drop a candidate."""
    penalty = np.r_[np.zeros(n_baseline), np.ones(X.shape[1]-n_baseline)]
    records = []
    for attempt, settings in enumerate((INITIAL, RETRY), start=1):
        model, record = solve_once(X, y, sample_weight, lam, rho, penalty,
                                   penalty, settings, context, attempt, log_path)
        records.append(record)
        if record['passed']:
            return model, records
    raise RuntimeError('Unresolved numerical fit; affected tuning/analysis invalid: '+json.dumps(context))
