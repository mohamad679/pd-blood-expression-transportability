"""Synthetic-first executor checks and resource benchmarks; no PPMI inputs.

No real-label permutation is performed here. The real observed-label identity
replay is a separate command, required before the public amendment is ready.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
from scipy.stats import rankdata
from threadpoolctl import threadpool_info, threadpool_limits

try:
    from . import permutation_executor as pe
except ImportError:
    import permutation_executor as pe


def synthetic(seed=7411, n=200, p=40, real_sizes=False):
    geo = pe.engine()
    rng = np.random.default_rng(seed)
    sizes = (142, 239) if real_sizes else (n//2, n//2)
    n = sum(sizes)
    source = np.repeat(['A', 'B'], sizes)
    cases = (92, 119) if real_sizes else tuple(size//2 for size in sizes)
    y = np.concatenate([np.r_[np.zeros(size-c), np.ones(c)] for size, c in zip(sizes, cases)]).astype(int)
    x = rng.normal(size=(n, p)) + (source == 'A')[:, None] * rng.normal(0, 2, p)
    X = (rankdata(x, axis=1) - .5) / p  # no within-source PD effect
    return geo.Data(X, rng.normal(60, 8, n), rng.binomial(1, .5, n), y, source,
                    np.array([f'd{i:04d}' for i in range(n)]),
                    np.array([f's{i:04d}' for i in range(n)]),
                    np.array([f'g{i:05d}' for i in range(p)]),
                    [np.arange(i*2, i*2+2) for i in range(6)])


def numeric_summary(results):
    return [{key: row[key] for key in ('rep_id', 'label_map_sha256', 'statistics', 'fold_results')}
            for row in results]


_BARRIER = None


def benchmark_init(threads, barrier):
    global _BARRIER
    _BARRIER = barrier
    for name in pe.THREAD_ENV:
        os.environ[name] = str(threads)
    barrier.wait(timeout=120)


def benchmark_task(job):
    geo = pe.engine()
    started = time.monotonic()
    with threadpool_limits(limits=2):
        data = synthetic(seed=8100+job, p=16364, real_sizes=True)
        source = data.take(np.flatnonzero(data.source == 'B'))
        tr, _ = geo.donor_folds(source, 5, 'B/outer')[0]
        train = source.take(tr)
        prep = geo.Prep(train, 1000)
        fit, records = geo.checked_fit(prep, train, 'M5', .01, 0.,
            dict(scope='synthetic/resource', stage='stress', fold=job), None)
        X, _ = prep.transform(source, 'M5')
        probability = fit.predict(X)
        pools = threadpool_info()
    return dict(job=job, pid=os.getpid(), seconds=time.monotonic()-started,
                peak_RSS_bytes=pe.peak_rss_bytes(), KKT=records[-1]['diagnostics']['KKT_max'],
                n_train=len(train.y), n_features=len(fit.coef_),
                thread_pools=pools, probability=probability.tolist())


def resource_benchmark(workers):
    for name in pe.THREAD_ENV:
        os.environ[name] = '2'
    context = mp.get_context('spawn')
    barrier = context.Barrier(workers)
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, mp_context=context,
            initializer=benchmark_init, initargs=(2, barrier)) as pool:
        rows = list(pool.map(benchmark_task, range(4)))
    peaks = {}
    for row in rows:
        peaks[row['pid']] = max(peaks.get(row['pid'], 0), row['peak_RSS_bytes'])
    return dict(workers=workers, threads_per_worker=2, wall_seconds=time.monotonic()-started,
                observed_worker_processes=len(peaks), summed_worker_peak_RSS_bytes=sum(peaks.values()),
                memory_measurement='sum of per-process high-water marks; conservative, not simultaneous RSS',
                rows=rows)


def physical_memory():
    try:
        return int(subprocess.check_output(['sysctl', '-n', 'hw.memsize'], text=True).strip())
    except (OSError, subprocess.CalledProcessError, ValueError):
        return None


def swap_counters():
    try:
        text = subprocess.check_output(['vm_stat'], text=True)
        return {line.split(':')[0]: int(line.split(':')[1].strip().rstrip('.'))
                for line in text.splitlines() if line.startswith(('Swapouts:', 'Pageouts:'))}
    except (OSError, subprocess.CalledProcessError, ValueError):
        return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--work', type=Path, required=True)
    ap.add_argument('--benchmark', action='store_true')
    args = ap.parse_args()
    settings = json.loads(pe.CONFIG.read_text())
    versions = pe.locked_environment()
    geo = pe.engine()
    checks = []
    started = time.monotonic()

    def check(name, condition):
        if not condition:
            report.update(status='failed', failed_check=name, checks=checks,
                          seconds=time.monotonic()-started)
            pe.atomic_json(args.output, report)
            raise AssertionError(name)
        checks.append(name)
        pe.atomic_json(args.output, report)
        print('PASS', name, flush=True)

    args.work.mkdir(parents=True, exist_ok=True)
    run = args.work / ('synthetic_' + str(time.time_ns()))
    run.mkdir()
    report = dict(status='in_progress', participant_data_read=False,
                  real_randomized_fits_performed=False, PPMI_accessed=False,
                  code_sha256=pe.code_hashes(), config_sha256=pe.sha_file(pe.CONFIG),
                  verifier_sha256=pe.sha_file(__file__),
                  package_versions=versions, checks=checks)
    pe.atomic_json(args.output, report)
    data = synthetic()
    manifest = geo.split_manifest(data, {})
    layout = pe.outer_layout(data, manifest)
    changed, label_hash, seeds = pe.shuffled_data(data, layout, 0)
    check('ten_disjoint_outer_tests_cover_all_samples', len(layout)==10)
    check('permutation_changes_at_least_one_label', not np.array_equal(changed.y, data.y))
    check('each_outer_block_preserves_class_totals', all(
        sorted(changed.y[te]) == sorted(data.y[te]) for _, _, te in layout))
    reverse = pe.shuffled_data(data, list(reversed(layout)), 0)
    check('seed_and_labels_independent_of_layout_scheduling',
          label_hash == reverse[1] and seeds == reverse[2] and np.array_equal(changed.y, reverse[0].y))
    permuted_order = np.random.default_rng(90).permutation(len(data.y))
    reordered = data.take(permuted_order)
    check('donor_sorting_makes_label_map_row_order_invariant',
          label_hash == pe.shuffled_data(reordered, pe.outer_layout(reordered, manifest), 0)[1])
    check('different_rep_ids_derive_different_maps', label_hash != pe.shuffled_data(data, layout, 1)[1])
    check('identity_map_keeps_original_labels', np.array_equal(
        pe.shuffled_data(data, layout, 0, identity=True)[0].y, data.y))
    _, tr, _ = layout[0]
    train = changed.take(tr)
    inner = geo.donor_folds(train, 5, layout[0][0]+'/inner')
    check('inner_folds_rebuilt_from_current_labels', all(
        not set(train.donor[a]) & set(train.donor[b]) for a, b in inner))
    old_inner = next(item['inner'] for item in manifest['splits'] if item['scope'] == layout[0][0])
    fresh_ids = [train.sample[b].tolist() for _, b in inner]
    check('shuffled_inner_validation_memberships_differ_from_observed_strata',
          fresh_ids != [item['validation'] for item in old_inner])
    unchanged_train = data.take(tr)
    unchanged_inner = geo.donor_folds(unchanged_train, 5, layout[0][0]+'/inner')
    check('identity_inner_memberships_equal_frozen_observed_manifest',
          [unchanged_train.sample[b].tolist() for _, b in unchanged_inner] ==
          [item['validation'] for item in old_inner])
    repeated = data.take(np.repeat(np.arange(len(data.y)), 2))
    repeated = replace(repeated, sample=np.array([f'repeated_{i}' for i in range(len(repeated.y))]))
    rm = geo.split_manifest(repeated, {})
    rl = pe.outer_layout(repeated, rm)
    rd = pe.shuffled_data(repeated, rl, 0)[0]
    check('repeated_rows_of_donor_receive_identical_labels', all(
        len(np.unique(rd.y[rd.donor == d])) == 1 for d in np.unique(rd.donor)))

    small = dict(settings, max_genes=15, lambda_grid=[.1, 1.], rhos=[0., .5], expansion_rounds=1)
    inner3 = geo.donor_folds(data, 3, 'synthetic/equivalence')
    with threadpool_limits(limits=2):
        p1, f1, d1 = geo.train_suite(data, inner3, run/'full_suite', 'synthetic/equivalence',
            max_genes=15, grid=(.1, 1.), rhos=(0., .5), expansion_rounds=1)
        p2, f2, d2 = pe.train_pair(data, inner3, run/'pair', 'synthetic/equivalence', small)
    check('M4_M5_candidates_and_decisions_equal_frozen_engine', all(d1[m] == d2[m] for m in ('M4','M5')))
    check('M4_M5_coefficients_equal_frozen_engine', all(
        np.array_equal(f1[m].coef_, f2[m].coef_) and f1[m].intercept_ == f2[m].intercept_ for m in ('M4','M5')))
    check('preprocessing_equal_frozen_engine', p1.fingerprint() == p2.fingerprint())
    pr1 = geo.predict_suite(p1, f1, d1, data)
    pr2 = pe.predict_pair(p2, f2, d2, data)
    check('M4_M5_predictions_equal_frozen_engine', all(np.array_equal(pr1[m], pr2[m]) for m in ('M4','M5')))

    fake = run/'checkpoint'
    attempt = fake/'attempt_test'
    attempt.mkdir(parents=True)
    pe.atomic_json(attempt/'RESULTS.json', {'answer': 1})
    pe.complete_attempt(fake, attempt, {'code': 'test'}, {'answer': 1})
    check('valid_completed_checkpoint_reused', pe.validate_completion(fake, {'code': 'test'}) is not None)
    rejected = False
    try:
        pe.validate_completion(fake, {'code': 'changed'})
    except ValueError:
        rejected = True
    check('changed_provenance_rejected', rejected)
    pe.atomic_json(attempt/'RESULTS.json', {'answer': 2})
    rejected = False
    try:
        pe.validate_completion(fake, {'code': 'test'})
    except ValueError:
        rejected = True
    check('tampered_artifact_rejected', rejected)
    check('partial_attempt_without_COMPLETE_not_reused', pe.validate_completion(run/'partial', {}) is None)
    with pe.controller_lock(run/'lock'):
        rejected = False
        try:
            with pe.controller_lock(run/'lock'):
                pass
        except ValueError:
            rejected = True
        check('concurrent_controller_rejected', rejected)

    smoke = dict(settings, max_genes=15, lambda_grid=[.1], rhos=[0.], expansion_rounds=0,
                 threads_per_worker=1)
    binding = dict(synthetic_only=True, config=pe.object_hash(smoke))
    serial = pe.execute(data, manifest, smoke, run/'serial', binding, [0,1], 1)
    parallel_settings = dict(smoke, threads_per_worker=2)
    parallel = pe.execute(data, manifest, parallel_settings, run/'parallel', binding, [1,0], 2)
    check('worker_thread_and_order_changes_preserve_label_maps_and_fold_AUROCs',
          numeric_summary(serial) == numeric_summary(parallel))
    resumed = pe.execute(data, manifest, smoke, run/'serial', binding, [0,1], 1)
    check('resume_reuses_two_complete_reps', all(r['resumed'] for r in resumed))
    check('resume_statistics_equal_uninterrupted_statistics', numeric_summary(serial) == numeric_summary(resumed))
    # A preserved incomplete attempt models interruption before atomic commit.
    partial = run/'interrupted/rep_000/attempt_interrupted'
    partial.mkdir(parents=True)
    (partial/'preserved.txt').write_text('interrupted; never resume partial fit state')
    restarted = pe.execute(data, manifest, smoke, run/'interrupted', binding, [0], 1)
    check('interruption_preserves_old_attempt_and_recomputes_same_rep',
          (partial/'preserved.txt').exists() and numeric_summary(restarted) == numeric_summary(serial[:1]))
    failed_gate = False
    try:
        pe.verify_receipt({}, dict(code_sha256={}, config_sha256='none'), args.output)
    except ValueError:
        failed_gate = True
    check('real_permutation_gate_rejects_unpublished_or_incomplete_preflight', failed_gate)
    records = [json.loads(line) for path in run.rglob('fits.jsonl') for line in path.read_text().splitlines()]
    check('all_synthetic_numerical_attempts_pass_KKT', bool(records) and all(r['passed'] for r in records))
    report.update(synthetic_fit_attempts=len(records),
                  max_KKT=max(r['diagnostics']['KKT_max'] for r in records),
                  synthetic_null_smoke_statistics=[r['statistics'] for r in serial],
                  smoke_is_not_a_statistical_null_calibration=True)

    if args.benchmark:
        before = swap_counters()
        one = resource_benchmark(1)
        four = resource_benchmark(settings['workers'])
        after = swap_counters()
        errors = [float(np.max(np.abs(np.array(a['probability'])-np.array(b['probability']))))
                  for a,b in zip(one['rows'], four['rows'])]
        memory = physical_memory()
        budget = four['summed_worker_peak_RSS_bytes'] + pe.peak_rss_bytes() + 1024**3
        memory_pass = memory is not None and budget <= .75 * memory
        swap_delta = {k: after[k] - before[k] for k in before.keys() & after.keys()}
        for b in (one, four):
            for row in b['rows']:
                row.pop('probability')
        report.update(resource_benchmark=dict(serial=one, parallel=four,
            max_abs_probability_errors=errors, physical_RAM_bytes=memory,
            conservative_memory_budget_bytes=budget, system_swap_counter_deltas=swap_delta,
            system_swap_counters_before=before, system_swap_counters_after=after,
            full_permutation_runtime_not_estimated_from_this_fit_microbenchmark=True))
        pe.atomic_json(args.output, report)
        check('1000_gene_stress_all_KKT_pass', all(r['KKT'] <= 1e-6 for b in (one, four) for r in b['rows']))
        check('1000_gene_worker_equivalence_max_probability_error_1e_minus_6', max(errors) <= 1e-6)
        check('benchmark_uses_requested_worker_count', four['observed_worker_processes'] == settings['workers'])
        check('benchmark_backend_threads_bounded_at_two', all(
            pool['num_threads'] <= 2 for row in four['rows'] for pool in row['thread_pools']))
        check('measured_memory_plus_guard_under_75_percent_physical_RAM', memory_pass)
        # macOS exposes pageouts and compressor swapouts as distinct counters.
        # Pageouts are retained diagnostically; they do not establish swap growth.
        check('no_system_compressor_swapout_increase_during_benchmark',
              'Swapouts' in swap_delta and swap_delta['Swapouts'] == 0)
    report.update(status='synthetic_passed_identity_pending', checks=checks,
                  check_count=len(checks), seconds=time.monotonic()-started,
                  observed_identity_replay='required_separately_before_release',
                  memory_benchmark_completed=args.benchmark)
    pe.atomic_json(args.output, report)
    print('SYNTHETIC CHECKS PASSED; observed identity replay still required.', flush=True)


if __name__ == '__main__':
    main()
