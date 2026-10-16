"""One prespecified full-grid, full-training-sized numerical stress fixture."""

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import rankdata

from block_logistic_v3_1 import balanced_source_weights, kish_effective_n
from numerical_kkt_v3_1 import fit_checked


def draw_problem():
    rng = np.random.default_rng(20261019)
    source = np.r_[np.zeros(142, dtype=int), np.ones(239, dtype=int)]
    y = np.r_[rng.permutation(np.r_[np.ones(93), np.zeros(49)]),
              rng.permutation(np.r_[np.ones(120), np.zeros(119)])]
    age = rng.normal(size=381)
    sex = rng.binomial(1, 0.5, 381)
    baseline = np.column_stack((age, age**2, sex, age*sex, rng.normal(size=(381, 6))))
    genes = rng.normal(size=(381, 1000))
    genes += (1-2*source)[:, None]*np.r_[np.full(500, 2.), np.full(500, -2.)][None, :]
    genes = (rankdata(genes, axis=1, method='average')-0.5)/1000
    X = np.column_stack((baseline, genes))
    X = (X-X.mean(axis=0))/X.std(axis=0)
    return np.ascontiguousarray(X), y, source


def main():
    root = Path(__file__).parent
    scale = json.loads((root/'kkt_scale_v3_1_solver_v2.json').read_text())
    if scale['status'] != 'passed':
        raise RuntimeError('Scale/calibration gate must pass before stress run')
    out, log = root/'stress_1000_kkt_v3_1_solver_v2.json', root/'stress_1000_kkt_v3_1_solver_v2_fits.jsonl'
    if out.exists() or log.exists():
        raise FileExistsError('Do not overwrite original stress outputs')
    X, y, source = draw_problem()
    rows = []
    start = time.monotonic()
    for weighted in (False, True):
        weights = balanced_source_weights(y, source) if weighted else None
        candidates = [None]+[(lam, rho) for lam in (1., 0.1, 0.01, 0.001, 0.0001) for rho in (0., 0.5, 1.)]
        for candidate in candidates:
            lam, rho = (0., 0.) if candidate is None else candidate
            data = np.ascontiguousarray(X[:, :10]) if candidate is None else X
            context = dict(phase='1000_gene_stress', stage='final_like',
                           arm='weighted' if weighted else 'unweighted', rep=0, fold=None,
                           baseline_only=candidate is None, seed=20261019)
            try:
                model, records = fit_checked(data, y, 10, lam, rho, weights, context, log)
                passed = True
                last = records[-1]
                row = dict(**context, **{'lambda': lam, 'rho': rho}, passed=passed,
                           attempts=len(records), final_KKT=last['diagnostics']['KKT_max'],
                           total_warnings=sum(len(r['warnings']) for r in records),
                           final_nonzero_genes=int(np.sum(model.coef_[10:] != 0)),
                           seconds=sum(r['seconds'] for r in records),
                           train_Kish=kish_effective_n(weights) if weighted else 381.)
            except RuntimeError as error:
                row = dict(**context, **{'lambda': lam, 'rho': rho}, passed=False,
                           exception=str(error))
            rows.append(row)
            print(json.dumps(row), flush=True)
            out.with_suffix('.inprogress.json').write_text(json.dumps(dict(completed=len(rows), rows=rows), indent=2)+'\n')
    attempts = [json.loads(line) for line in log.read_text().splitlines()]
    result = dict(status='passed' if all(r['passed'] for r in rows) else 'review_required',
                  synthetic_only=True, participant_data_read=False, seed=20261019,
                  n=381, genes=1000, baseline_features=10, total_final_like_fits=len(rows),
                  fits_with_initial_KKT_failure=sum(not r['passed'] for r in attempts if r['attempt'] == 1),
                  unresolved_fits=sum(not r['passed'] for r in rows),
                  warning_count=sum(len(r['warnings']) for r in attempts),
                  fits_with_warnings=sum(r.get('total_warnings', 0) > 0 for r in rows),
                  max_final_KKT=max(r['final_KKT'] for r in rows if r['passed']),
                  lambda_grid_changed=False, KKT_threshold_changed=False,
                  elapsed_seconds=round(time.monotonic()-start, 2),
                  limitations=['One numerical fixture, not validation of all possible training folds.',
                               'No new-platform transfer or statistical power assessed.'], rows=rows)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
