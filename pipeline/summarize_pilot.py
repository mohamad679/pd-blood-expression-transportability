"""Read-only aggregation of the timing pilot; never refits or changes selection."""
import argparse
import json
from pathlib import Path
import numpy as np
from geo_pipeline import write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    p = args.directory
    if not (p / 'RESULTS.json').exists():
        raise RuntimeError('Pilot has not completed; no final budget report written')
    result = json.loads((p / 'RESULTS.json').read_text())
    started = json.loads((p / 'STARTED.json').read_text())
    rows = [json.loads(x) for x in (p / 'fits.jsonl').read_text().splitlines()]
    choices = json.loads((p / 'decisions.json').read_text())
    grouped = {}
    for model in ('M1', 'M2', 'M3', 'M4', 'M5'):
        fits = [r for r in rows if r['model'] == model]
        seconds = np.array([r['seconds'] for r in fits])
        grouped[model] = dict(attempts=len(fits), seconds=float(sum(seconds)),
            median_seconds=float(np.median(seconds)), max_seconds=float(max(seconds)),
            warned_attempts=sum(bool(r['warnings']) for r in fits),
            max_KKT=max(r['diagnostics']['KKT_max'] for r in fits if r['diagnostics']))
    summary = dict(status='complete', n_train=started['n_train'], n_test=started['n_test'],
        elapsed_wall_seconds=result['seconds'], process_CPU_seconds=result['process_CPU_seconds'],
        fit_attempts=len(rows), failed_attempts=sum(not r['passed'] for r in rows),
        retries=sum(r['attempt'] > 1 for r in rows),
        warning_messages=sum(len(r['warnings']) for r in rows),
        warned_attempts=sum(bool(r['warnings']) for r in rows),
        max_KKT=max(r['diagnostics']['KKT_max'] for r in rows if r['diagnostics']),
        by_model=grouped,
        selected={m: {k: v for k, v in d.items() if k != 'candidates'}
                  for m, d in choices.items()},
        same_fold_timing_scenarios=dict(
            core_13_analyses_wall_hours=13 * result['seconds'] / 3600,
            core_plus_four_full_sensitivity_runs_hours=65 * result['seconds'] / 3600,
            permutation_full_five_model_1000x10_outer_analyses_wall_days=
                10000 * result['seconds'] / 86400,
            permutation_M4_M5_only_1000x10_fit_time_days_excluding_preprocessing=
                10000 * (grouped['M4']['seconds'] + grouped['M5']['seconds']) / 86400),
        forecast_limitation='Scenario arithmetic only: one small PROBE fold does not establish larger-cohort, pooled or sensitivity runtimes. Backend internal threads are recorded; model fits are sequential. Permutation execution scope is still to be fixed before its execution.',
        second_pre_PPMI_lock_done=False, PPMI_scores_computed=False,
        single_fold_performance_is_study_result=False)
    write_json(p / 'TIMING_BUDGET.json', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
