"""Preregistered local numerical scale/calibration checks; synthetic data only."""

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from numerical_kkt_v3_1 import KKT_LIMIT, REFERENCE_GLM, objective_and_kkt, solve_once


def synthetic_problem(seed, genes):
    rng = np.random.default_rng(np.random.SeedSequence([seed, genes]))
    n, k = 381, 10
    B = rng.normal(size=(n, k))
    G = rng.normal(size=(n, genes))
    G[:, 0] += 0.5*B[:, 0]
    X = np.column_stack((B, G))
    X = (X-X.mean(axis=0))/X.std(axis=0)
    y = rng.binomial(1, expit(-0.2+0.7*X[:, 0]-0.5*X[:, 1]+0.6*X[:, k])).astype(float)
    weights = np.exp(rng.normal(0, 0.7, n))
    p1 = np.r_[np.zeros(k), rng.uniform(0.6, 1.6, genes)]
    p2 = np.r_[np.zeros(k), rng.uniform(0.5, 1.5, genes)]
    return X, y, weights, p1, p2


def independent_solver(X, y, weights, lam, rho, p1, p2):
    """Separate SciPy implementation; L1 uses nonnegative split gene variables."""
    k = 10
    B = np.column_stack((np.ones(len(y)), X[:, :k]))
    G = X[:, k:]
    s = np.ones(len(y))/len(y) if weights is None else weights/weights.sum()
    size = B.shape[1]
    p = G.shape[1]
    start = time.monotonic()

    if rho == 0:
        D = np.column_stack((np.ones(len(y)), X))
        diagonal = np.r_[0, p2]

        def evaluate(theta):
            eta = D @ theta
            objective = np.dot(s, np.logaddexp(0, eta)-y*eta)+lam/2*np.dot(diagonal, theta**2)
            gradient = D.T @ (s*(expit(eta)-y))+lam*diagonal*theta
            return float(objective), gradient

        result = minimize(evaluate, np.zeros(X.shape[1]+1), method='L-BFGS-B', jac=True,
                          options=dict(ftol=1e-14, gtol=1e-10, maxiter=20000, maxls=100))
        theta = result.x
    else:
        def evaluate(split):
            beta = split[:size]
            positive, negative = split[size:size+p], split[size+p:]
            gamma = positive-negative
            eta = B @ beta + G @ gamma
            objective = (np.dot(s, np.logaddexp(0, eta)-y*eta)
                         + lam*rho*np.dot(p1[k:], positive+negative)
                         + lam*(1-rho)/2*np.dot(p2[k:], gamma**2))
            residual = s*(expit(eta)-y)
            gene_gradient = G.T @ residual+lam*(1-rho)*p2[k:]*gamma
            gradient = np.r_[B.T @ residual, gene_gradient+lam*rho*p1[k:],
                             -gene_gradient+lam*rho*p1[k:]]
            return float(objective), gradient

        result = minimize(evaluate, np.zeros(size+2*p), method='L-BFGS-B', jac=True,
                          bounds=[(None, None)]*size+[(0, None)]*(2*p),
                          options=dict(ftol=1e-14, gtol=1e-10, maxiter=20000, maxls=100))
        theta = np.r_[result.x[:size], result.x[size:size+p]-result.x[size+p:]]
    diagnostic = objective_and_kkt(X, y, weights, theta, lam, rho, p1, p2)
    return theta, {**diagnostic, 'scipy_success': bool(result.success),
                   'message': str(result.message), 'iterations': int(result.nit),
                   'seconds': round(time.monotonic()-start, 4)}


def main():
    root = Path(__file__).parent
    log = root/'kkt_scale_v3_1_solver_v2_fits.jsonl'
    out = root/'kkt_scale_v3_1_solver_v2.json'
    if log.exists() or out.exists():
        raise FileExistsError('Preserve prior run: choose a versioned output before rerunning')
    rows = []
    for phase, seed in (('calibration', 20261017), ('independent_validation', 20261018)):
        for genes in (40, 1000):
            X, y, w, p1, p2 = synthetic_problem(seed, genes)
            for weighted in (False, True):
                weights = w if weighted else None
                for rho in (0., 0.5, 1.):
                    for lam in (0.01, 0.2):
                        context = dict(phase=phase, stage='reference_scale', arm='nonuniform' if weighted else 'uniform',
                                       rep=0, fold=None, seed=seed, genes=genes)
                        model, fit = solve_once(X, y, weights, lam, rho, p1, p2,
                                                REFERENCE_GLM, context, 1, log)
                        reference, ref = independent_solver(X, y, weights, lam, rho, p1, p2)
                        if fit['diagnostics'] is None:
                            error, probability_error, wrong_scale = None, None, None
                        else:
                            error = abs(fit['diagnostics']['objective']-ref['objective'])
                            theta = np.r_[model.intercept_, model.coef_]
                            probability_error = float(np.max(np.abs(expit(theta[0]+X@theta[1:])-expit(reference[0]+X@reference[1:]))))
                            wrong_scale = objective_and_kkt(X, y, weights, theta, lam, rho, p1, p2,
                                                            loss_gradient_multiplier=2)['KKT_max']
                        passed = bool(fit['passed'] and ref['finite'] and ref['KKT_max'] <= KKT_LIMIT
                                      and error is not None and error <= 1e-8)
                        rows.append({**context, 'lambda': lam, 'rho': rho, 'glum': fit,
                                     'independent': ref, 'objective_difference': error,
                                     'max_probability_difference': probability_error,
                                     'wrong_double_loss_gradient_KKT': wrong_scale, 'passed': passed})
                        print(f'{phase} p={genes} weighted={weighted} lambda={lam} rho={rho}: '
                              f'pass={passed}, glum KKT={None if fit["diagnostics"] is None else fit["diagnostics"]["KKT_max"]}, '
                              f'SciPy KKT={ref["KKT_max"]}', flush=True)
                        out.with_suffix('.inprogress.json').write_text(json.dumps({'rows': rows, 'completed': len(rows)}, indent=2)+'\n')
    detectable = [r for r in rows if r['glum']['diagnostics'] and r['glum']['diagnostics']['nonzero_penalized_coefficients'] > 0]
    negative_pass = bool(detectable) and all(r['wrong_double_loss_gradient_KKT'] > KKT_LIMIT for r in detectable)
    result = dict(status='passed' if all(r['passed'] for r in rows) and negative_pass else 'review_required',
                  synthetic_only=True, participant_data_read=False, KKT_limit=KKT_LIMIT,
                  threshold_changed_after_outputs=False, checks=len(rows),
                  calibration_seed=20261017, validation_seed=20261018,
                  wrong_factor_two_detected=negative_pass,
                  glum_KKT_max=max(r['glum']['diagnostics']['KKT_max'] for r in rows if r['glum']['diagnostics']),
                  scipy_KKT_max=max(r['independent']['KKT_max'] for r in rows),
                  max_objective_difference=max(r['objective_difference'] for r in rows if r['objective_difference'] is not None),
                  failed_checks=sum(not r['passed'] for r in rows), rows=rows)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
