"""Same 200 paired source-null seeds, now with per-fit final KKT and cold retry.

Historical code/results remain unchanged. New output never overwrites the old
source_null_v3_1.json. This simulation still has two sources, not a new platform.
"""

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from block_logistic_v3_1 import balanced_source_weights, kish_effective_n
from numerical_kkt_v3_1 import fit_checked
from simulate_source_null_v3_1 import LAMBDAS, RHOS, draw_data, standardize, macro_auc, wilson


def select_and_evaluate(X, y, source, XT, yT, sourceT, folds, corrected, repeat, log_path):
    cached = []
    arm = 'weighted_macro_tuning' if corrected else 'unweighted_pooled_tuning'
    for train, valid in folds:
        xtr, xval = standardize(X[train], X[valid])
        weights = balanced_source_weights(y[train], source[train]) if corrected else None
        cached.append((xtr, xval, train, valid, weights))
    score_cache = {}

    def candidate_score(candidate):
        if candidate in score_cache:
            return score_cache[candidate]
        scores = []
        for fold, (xtr, xval, tr, va, weights) in enumerate(cached):
            lam, rho = (0., 0.) if candidate is None else candidate
            train_design = np.ascontiguousarray(xtr[:, :10]) if candidate is None else xtr
            valid_design = np.ascontiguousarray(xval[:, :10]) if candidate is None else xval
            context = dict(phase='source_null_rerun', stage='inner', arm=arm, rep=repeat,
                           fold=fold, baseline_only=candidate is None)
            model, _ = fit_checked(train_design, y[tr], 10, lam, rho, weights, context, log_path)
            pred = model.predict(valid_design)
            scores.append(macro_auc(y[va], pred, source[va]) if corrected else float(roc_auc_score(y[va], pred)))
        score_cache[candidate] = float(np.mean(scores))
        return score_cache[candidate]

    grid, expansions = list(LAMBDAS), 0
    while True:
        selected, best = None, candidate_score(None)
        for candidate in [(lam, rho) for lam in sorted(grid, reverse=True) for rho in sorted(RHOS, reverse=True)]:
            score = candidate_score(candidate)
            if score > best+1e-6:
                selected, best = candidate, score
        if selected is None:
            break
        lam, rho = selected
        if expansions >= 2 or lam not in (min(grid), max(grid)):
            break
        grid.append(min(grid)/10 if lam == min(grid) else max(grid)*10)
        expansions += 1
    xtr, xt = standardize(X, XT)
    weights = balanced_source_weights(y, source) if corrected else None
    lam, rho = (0., 0.) if selected is None else selected
    train_design = np.ascontiguousarray(xtr[:, :10]) if selected is None else xtr
    test_design = np.ascontiguousarray(xt[:, :10]) if selected is None else xt
    context = dict(phase='source_null_rerun', stage='final', arm=arm, rep=repeat,
                   fold=None, baseline_only=selected is None)
    model, _ = fit_checked(train_design, y, 10, lam, rho, weights, context, log_path)
    nonzero = int(np.sum(np.abs(model.coef_[10:]) > 1e-8))
    pred = model.predict(test_design)
    return dict(gene_selected=bool(selected is not None and nonzero > 0), selected=selected,
                inner_selection_score=best, independent_test_pooled_AUC=float(roc_auc_score(yT, pred)),
                independent_test_macro_AUC=macro_auc(yT, pred, sourceT), grid_expansions=expansions,
                final_grid=sorted(grid), selected_at_remaining_boundary=bool(selected is not None and nonzero > 0 and selected[0] in (min(grid), max(grid))),
                nonzero_genes=nonzero, train_kish=kish_effective_n(weights) if corrected else float(len(y)))


def main():
    root = Path(__file__).parent
    parser = argparse.ArgumentParser()
    parser.add_argument('--overlap-stress-no-release', action='store_true',
                        help='Compute independent null while stress finishes; never substitutes for stress/release gate.')
    args = parser.parse_args()
    stress = root/'stress_1000_kkt_v3_1_solver_v2.json'
    stress_passed = stress.exists() and json.loads(stress.read_text())['status'] == 'passed'
    if not stress_passed and not args.overlap_stress_no_release:
        raise RuntimeError('Numerical 1000-gene fixture must pass before full rerun')
    if not stress_passed:
        print('Computational overlap only: stress gate pending; no release or protocol lock authorized by this run.', flush=True)
    out, log = root/'source_null_kkt_v3_1.json', root/'source_null_kkt_v3_1_fits.jsonl'
    if out.exists() or log.exists():
        raise FileExistsError('Version the outputs rather than overwrite an earlier run')
    repetitions, rows = 200, []
    start = time.monotonic()
    try:
        for repeat in range(repetitions):
            rng = np.random.default_rng(np.random.SeedSequence([20261016, repeat]))
            X, y, s = draw_data(rng)
            XT, yT, sT = draw_data(rng)
            folds = list(StratifiedKFold(5, shuffle=True, random_state=20261016+repeat).split(X, s*2+y))
            old = select_and_evaluate(X, y, s, XT, yT, sT, folds, False, repeat, log)
            new = select_and_evaluate(X, y, s, XT, yT, sT, folds, True, repeat, log)
            rows.append(dict(repeat=repeat, unweighted_pooled_tuning=old, weighted_macro_tuning=new))
            if (repeat+1) % 10 == 0:
                print(f'Completed {repeat+1}/200 fully audited paired repetitions', flush=True)
                out.with_suffix('.inprogress.json').write_text(json.dumps(dict(completed=len(rows), rows=rows), indent=2)+'\n')
    except RuntimeError as exception:
        out.write_text(json.dumps(dict(status='invalid_tuning', completed=len(rows), error=str(exception), rows=rows), indent=2)+'\n')
        raise
    summary = {}
    for arm in ('unweighted_pooled_tuning', 'weighted_macro_tuning'):
        chosen = np.array([r[arm]['gene_selected'] for r in rows], dtype=float)
        metrics = {}
        for metric in ('independent_test_pooled_AUC', 'independent_test_macro_AUC'):
            draws = np.array([r[arm][metric] for r in rows])
            metrics[metric] = dict(mean=float(draws.mean()), MCSE=float(draws.std(ddof=1)/np.sqrt(repetitions)))
        summary[arm] = dict(gene_selection_count=int(chosen.sum()), gene_selection_rate=float(chosen.mean()),
                            gene_selection_Wilson_CI=wilson(int(chosen.sum()), repetitions),
                            final_boundary_count=sum(r[arm]['selected_at_remaining_boundary'] for r in rows),
                            expanded_repetitions=sum(r[arm]['grid_expansions'] > 0 for r in rows), **metrics)
    attempts = [json.loads(line) for line in log.read_text().splitlines()]
    paired = np.array([int(r['unweighted_pooled_tuning']['gene_selected'])-int(r['weighted_macro_tuning']['gene_selected']) for r in rows])
    diff = float(paired.mean())
    se = float(paired.std(ddof=1)/np.sqrt(repetitions))
    ci = [diff-norm.ppf(0.975)*se, diff+norm.ppf(0.975)*se]
    reasons = []
    if ci[0] <= 0:
        reasons.append('Selection-rate reduction not established; review before registration.')
    for arm, stats in summary.items():
        m = stats['independent_test_macro_AUC']
        if abs(m['mean']-0.5) > 3*m['MCSE']+0.01:
            reasons.append(arm+': independent null macro-AUC unexpectedly far from chance.')
    history = json.loads((root/'source_null_v3_1.json').read_text())
    changes = []
    for arm in summary:
        changes.append(dict(arm=arm, changed_gene_selection_repetitions=sum(a[arm]['gene_selected'] != b[arm]['gene_selected'] for a, b in zip(rows, history['rows'])),
                            max_test_macro_AUC_change=max(abs(a[arm]['independent_test_macro_AUC']-b[arm]['independent_test_macro_AUC']) for a, b in zip(rows, history['rows']))))
    result = dict(status='passed' if not reasons else 'review_required', synthetic_only=True,
                  computational_overlap_with_stress=bool(args.overlap_stress_no_release),
                  null_status_does_not_replace_stress_or_publication_gate=True,
                  participant_data_read=False, repetitions=200, seed=20261016,
                  n_train=381, n_test=381, genes=40, baseline_features=10,
                  summary=summary, paired_selection_rate_reduction=diff,
                  paired_selection_rate_reduction_MC_CI=ci, review_reasons=reasons,
                  fit_attempts=len(attempts), retry_attempts=sum(a['attempt'] == 2 for a in attempts),
                  initial_KKT_failures=sum(a['attempt'] == 1 and not a['passed'] for a in attempts),
                  warning_messages=sum(len(a['warnings']) for a in attempts),
                  warned_attempts=sum(bool(a['warnings']) for a in attempts),
                  warned_attempts_passing_final_KKT=sum(bool(a['warnings']) and a['passed'] for a in attempts),
                  final_accepted_KKT_max=max(a['diagnostics']['KKT_max'] for a in attempts if a['passed']),
                  historical_comparison=changes, rows=rows,
                  elapsed_seconds=round(time.monotonic()-start, 2),
                  limitations=history['limitations']+['No synthetic third-source transfer evaluated.'])
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
