"""Resumable M4/M5 permutations. Real randomized fits require a public receipt.

The frozen GEO engine is imported lazily and remains unchanged. Checkpoints and
all predictions stay private. No PPMI loading or scoring exists in this module.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import replace
import hashlib
import importlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import resource
import socket
import sys
import time
import uuid
import warnings

import numpy as np
from scipy.special import expit
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CONFIG = ROOT / 'protocol/permutation_config_v1.json'
THREAD_ENV = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


def engine():
    return importlib.import_module('pipeline.geo_pipeline')


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    with tmp.open('w', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def code_hashes():
    return {name: sha_file(ROOT / name) for name in (
        'pipeline/permutation_executor.py', 'pipeline/geo_pipeline.py',
        'pipeline/verify_permutation_executor.py', 'pipeline/finalize_permutation_preflight.py',
        'protocol/numerical_kkt_v3_1.py', 'protocol/block_logistic_v3_1.py',
        'protocol/config_v3_1.json', 'protocol/v3_marker_manifest.json',
        'requirements_glum_v3.lock.txt')}


def locked_environment():
    actual = {}
    for line in (ROOT / 'requirements_glum_v3.lock.txt').read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        package, expected = line.split('==')
        actual[package] = importlib.metadata.version(package)
        if actual[package] != expected:
            raise ValueError(f'Locked environment mismatch: {package}')
    return actual


def outer_layout(data, manifest):
    """Validate fixed outer membership; never rebuild it with shuffled labels."""
    lookup = {str(sample): i for i, sample in enumerate(data.sample)}
    if len(lookup) != len(data.sample):
        raise ValueError('Duplicate sample IDs')
    result, seen, donor_block = [], set(), {}
    for item in manifest['splits']:
        if '/outer/' not in item['scope']:
            continue
        scope = item['scope']
        tr = np.array([lookup[s] for s in item['train']], int)
        te = np.array([lookup[s] for s in item['test']], int)
        source = scope.split('/')[0]
        members = set(np.flatnonzero(data.source == source))
        if set(tr) & set(te) or set(tr) | set(te) != members or set(te) & seen:
            raise ValueError('Invalid outer partition')
        if set(data.donor[tr]) & set(data.donor[te]):
            raise ValueError('Outer donor leakage')
        for donor in np.unique(data.donor[te]):
            if str(donor) in donor_block:
                raise ValueError('Donor appears in multiple outer blocks')
            donor_block[str(donor)] = scope
        seen.update(te)
        result.append((scope, tr, te))
    if seen != set(range(len(data.y))):
        raise ValueError('Outer tests do not cover every sample once')
    return result


def shuffled_data(data, layout, rep_id, identity=False):
    """One coherent donor map per rep, preserving every outer block's totals."""
    if rep_id < 0:
        raise ValueError('Negative replicate ID')
    y = np.array(data.y, copy=True)
    label_map = {}
    seeds = {}
    for source in sorted(np.unique(data.source)):
        seed = int.from_bytes(hashlib.sha256(
            f'20261014:permutation/{rep_id}/{source}'.encode()).digest()[:4], 'big')
        seeds[str(source)] = seed
        rng = np.random.default_rng(seed)
        for scope, _, te in sorted(layout, key=lambda entry: entry[0]):
            if scope.split('/')[0] != source:
                continue
            donors = np.sort(np.unique(data.donor[te]))
            original = []
            for donor in donors:
                values = np.unique(data.y[data.donor == donor])
                if len(values) != 1:
                    raise ValueError('Inconsistent donor label')
                original.append(int(values[0]))
            permuted = original if identity else rng.permutation(original).tolist()
            for donor, value in zip(donors, permuted):
                y[data.donor == donor] = value
                label_map[str(donor)] = int(value)
    return replace(data, y=y), object_hash(label_map), seeds


def train_pair(train, inner, out, scope, settings):
    """M4/M5-only equivalent of the locked train_suite; all fits remain cold."""
    geo = engine()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    log = out / 'fits.jsonl'
    if log.exists():
        raise ValueError('Attempt already contains fit records')
    fold_data = []
    for ti, vi in inner:
        a, b = train.take(ti), train.take(vi)
        if set(a.donor) & set(b.donor):
            raise ValueError('Inner donor leakage')
        fold_data.append((geo.Prep(a, settings['max_genes']), a, b))
    final_prep = geo.Prep(train, settings['max_genes'])
    final_prep.save(out / 'preprocessing.npz')
    scores = []
    for f, (prep, a, b) in enumerate(fold_data):
        fit, _ = geo.checked_fit(prep, a, 'M4', 0., 0.,
            dict(scope=scope, stage='inner', fold=f, baseline_only=True), log)
        X, _ = prep.transform(b, 'M4')
        scores.append(geo.macro_auc(b, fit.predict(X)))
    baseline_score = float(np.mean(scores))
    base, _ = geo.checked_fit(final_prep, train, 'M4', 0., 0.,
        dict(scope=scope, stage='outer_or_final_refit', fold=None, baseline_only=True), log)
    candidates = [dict(baseline=True, **{'lambda': 0.}, rho=0.,
                       score=baseline_score, nonzero=False)]
    lambdas = list(settings['lambda_grid'])
    for round_i in range(settings['expansion_rounds'] + 1):
        existing = {(c['lambda'], c['rho']) for c in candidates if not c['baseline']}
        for lam in sorted(lambdas, reverse=True):
            for rho in settings['rhos']:
                if (lam, rho) in existing:
                    continue
                scores, nonzero = [], False
                for f, (prep, a, b) in enumerate(fold_data):
                    fit, records = geo.checked_fit(prep, a, 'M5', lam, rho,
                        dict(scope=scope, stage='inner', fold=f, baseline_only=False,
                             expansion_round=round_i), log)
                    X, _ = prep.transform(b, 'M5')
                    scores.append(geo.macro_auc(b, fit.predict(X)))
                    nonzero |= records[-1]['diagnostics']['nonzero_penalized_coefficients'] > 0
                candidates.append(dict(baseline=False, **{'lambda': float(lam)},
                    rho=float(rho), score=float(np.mean(scores)), nonzero=bool(nonzero)))
                atomic_json(out / 'M5_tuning.inprogress.json', candidates)
        selected = geo.choose(candidates)
        edge = (not selected['baseline'] and selected['nonzero']
                and selected['lambda'] in (min(lambdas), max(lambdas)))
        if not edge or round_i == settings['expansion_rounds']:
            break
        lambdas.append(min(lambdas) / 10 if selected['lambda'] == min(lambdas)
                       else max(lambdas) * 10)
    chosen = dict(selected)
    if selected['baseline']:
        extended = base
        chosen['reuses_baseline'] = 'M4'
    else:
        extended, records = geo.checked_fit(final_prep, train, 'M5',
            selected['lambda'], selected['rho'],
            dict(scope=scope, stage='outer_or_final_refit', fold=None, baseline_only=False), log)
        if records[-1]['diagnostics']['nonzero_penalized_coefficients'] == 0:
            extended = base
            chosen['reuses_baseline'] = 'M4'
            chosen['final_gene_block_zero'] = True
    decisions = dict(M4=dict(hyperparameters=None, inner_macro_AUROC=baseline_score),
        M5=dict(selected=chosen, candidates=candidates, persistent_boundary=bool(edge),
                expansion_rounds=round_i))
    fits = dict(M4=base, M5=extended)
    for name, fit in fits.items():
        np.savez_compressed(out / (name + '_model.npz'), intercept=fit.intercept_, coef=fit.coef_)
    atomic_json(out / 'decisions.json', decisions)
    return final_prep, fits, decisions


def predict_pair(prep, fits, decisions, test):
    predictions = {}
    for name in ('M4', 'M5'):
        actual = decisions[name].get('selected', {}).get('reuses_baseline', name)
        if actual != name:
            predictions[name] = predictions[actual].copy()
        else:
            X, _ = prep.transform(test, name)
            predictions[name] = expit(fits[name].intercept_ + X @ fits[name].coef_)
    return predictions


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == 'darwin' else value * 1024)


def validate_completion(directory, binding):
    """Never trust a completed flag without checking binding and every artifact."""
    directory = Path(directory)
    complete = directory / 'COMPLETE.json'
    if not complete.exists():
        return None
    record = json.loads(complete.read_text())
    if record['binding_sha256'] != object_hash(binding):
        raise ValueError('Checkpoint provenance mismatch')
    if not record['artifacts'] or record['result_file'] not in record['artifacts']:
        raise ValueError('Incomplete checkpoint artifact manifest')
    for relative, expected in record['artifacts'].items():
        path = (directory / relative).resolve()
        if directory.resolve() not in path.parents or sha_file(path) != expected:
            raise ValueError('Checkpoint artifact checksum mismatch')
    if json.loads((directory / record['result_file']).read_text()) != record['result']:
        raise ValueError('Checkpoint result differs from its checked artifact')
    return record


def complete_attempt(directory, attempt, binding, result):
    artifacts = {str(path.relative_to(directory)): sha_file(path)
                 for path in sorted(attempt.rglob('*')) if path.is_file()}
    if not artifacts:
        raise ValueError('Cannot complete an empty attempt')
    result_file = str((attempt / 'RESULTS.json').relative_to(directory))
    record = dict(binding_sha256=object_hash(binding), artifacts=artifacts,
                  result_file=result_file, result=result)
    atomic_json(Path(directory) / 'COMPLETE.json', record)
    return record


def _run_repetition_unlocked(data, manifest, settings, root, binding, rep_id, identity=False):
    """A crash leaves an uncommitted attempt; retry keeps it and starts cold."""
    geo = engine()
    layout = outer_layout(data, manifest)
    permuted, label_hash, seeds = shuffled_data(data, layout, rep_id, identity)
    binding = dict(binding, rep_id=rep_id, identity=identity, label_map_sha256=label_hash)
    directory = Path(root) / ('identity' if identity else f'rep_{rep_id:03d}')
    cached = validate_completion(directory, binding)
    if cached is not None:
        return dict(cached['result'], resumed=True)
    if list(directory.glob('attempt_*/FAILURE.json')):
        raise ValueError('A dependent replicate was declared invalid; do not retry it as an interruption')
    attempt = directory / ('attempt_' + uuid.uuid4().hex)
    attempt.mkdir(parents=True)
    atomic_json(attempt / 'STARTED.json', dict(binding=binding, seeds=seeds))
    started = time.monotonic()
    stats, fold_results, prediction_warnings = {}, [], []
    try:
        for scope, tr, te in layout:
            train, test = permuted.take(tr), permuted.take(te)
            inner = geo.donor_folds(train, settings['inner_folds'], scope + '/inner')
            out = attempt / scope.replace('/', '_')
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                prep, fits, decisions = train_pair(train, inner, out, scope, settings)
                probabilities = predict_pair(prep, fits, decisions, test)
            prediction_warnings.extend(dict(scope=scope, category=type(w.message).__name__,
                                            message=str(w.message)) for w in caught)
            if any(not np.isfinite(p).all() for p in probabilities.values()):
                raise ValueError('Nonfinite prediction')
            np.savez_compressed(out / 'heldout_predictions.npz', samples=test.sample,
                                y=test.y, **probabilities)
            scores = {m: float(roc_auc_score(test.y, probabilities[m])) for m in ('M4', 'M5')}
            fold_results.append(dict(scope=scope, **scores))
            atomic_json(attempt / 'progress.json', dict(completed_folds=len(fold_results),
                                                       fold_results=fold_results))
        for source in sorted(np.unique(data.source)):
            results = [r for r in fold_results if r['scope'].split('/')[0] == source]
            means = {m: float(np.mean([r[m] for r in results])) for m in ('M4', 'M5')}
            stats[str(source)] = dict(**means, delta=means['M5'] - means['M4'])
        result = dict(rep_id=rep_id, identity=identity, label_map_sha256=label_hash,
                      statistics=stats, fold_results=fold_results,
                      seconds=time.monotonic()-started, peak_RSS_bytes=peak_rss_bytes(),
                      PPMI_accessed=False)
        atomic_json(attempt / 'prediction_warnings.json', prediction_warnings)
        atomic_json(attempt / 'RESULTS.json', result)
        return dict(complete_attempt(directory, attempt, binding, result)['result'], resumed=False)
    except Exception as exc:
        atomic_json(attempt / 'FAILURE.json', dict(type=type(exc).__name__, message=str(exc),
                    completed_folds=len(fold_results), dependent_replicate_invalid=True))
        raise


def run_repetition(data, manifest, settings, root, binding, rep_id, identity=False):
    directory = Path(root) / ('identity' if identity else f'rep_{rep_id:03d}')
    with controller_lock(directory):
        return _run_repetition_unlocked(data, manifest, settings, root, binding, rep_id, identity)


_WORKER = None


def worker_init(data, manifest, settings, root, binding, identity):
    global _WORKER
    _WORKER = (data, manifest, settings, root, binding, identity)


def worker_task(rep_id):
    data, manifest, settings, root, binding, identity = _WORKER
    with threadpool_limits(limits=settings['threads_per_worker']):
        return run_repetition(data, manifest, settings, root, binding, rep_id, identity)


def execute(data, manifest, settings, root, binding, rep_ids, workers, identity=False):
    for variable in THREAD_ENV:
        os.environ[variable] = str(settings['threads_per_worker'])
    if workers == 1:
        worker_init(data, manifest, settings, root, binding, identity)
        return [worker_task(rep_id) for rep_id in rep_ids]
    results = []
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'),
            initializer=worker_init, initargs=(data, manifest, settings, root, binding, identity)) as pool:
        futures = {pool.submit(worker_task, rep_id): rep_id for rep_id in rep_ids}
        try:
            for future in as_completed(futures):
                results.append(future.result())
                print('Completed replicate', futures[future], flush=True)
        except Exception:
            for future in futures:
                future.cancel()
            raise
    return sorted(results, key=lambda r: r['rep_id'])


@contextmanager
def controller_lock(root):
    """Single local controller; demonstrably dead PIDs leave an audit record."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    path = root / 'CONTROLLER.lock'
    if path.exists():
        old = json.loads(path.read_text())
        if old['host'] != socket.gethostname():
            raise ValueError('Controller lock belongs to another host')
        try:
            os.kill(old['pid'], 0)
        except ProcessLookupError:
            path.rename(root / ('STALE_LOCK_' + uuid.uuid4().hex + '.json'))
        else:
            raise ValueError('Another controller is active')
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    with os.fdopen(descriptor, 'w') as stream:
        json.dump(dict(pid=os.getpid(), host=socket.gethostname()), stream)
    try:
        yield
    finally:
        path.unlink(missing_ok=True)


def verify_receipt(receipt, binding, preflight):
    report = json.loads(Path(preflight).read_text())
    if (receipt.get('publicly_verified') is not True
            or not receipt.get('record_url', '').startswith('https://zenodo.org/records/')
            or receipt.get('binding_sha256') != object_hash(binding)
            or receipt.get('preflight_report_sha256') != sha_file(preflight)
            or report.get('status') != 'passed'
            or report.get('binding_sha256') != object_hash(binding)
            or report.get('code_sha256') != binding['code_sha256']
            or report.get('config_sha256') != binding['config_sha256']):
        raise ValueError('Verified public amendment matching code/config/inputs/preflight is required')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=('plan', 'identity', 'run'))
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--config', type=Path, default=CONFIG)
    ap.add_argument('--receipt', type=Path)
    ap.add_argument('--preflight', type=Path)
    ap.add_argument('--original-runs', type=Path, default=ROOT / 'private_data/geo_runs_20261018')
    args = ap.parse_args()
    settings = json.loads(args.config.read_text())
    environment = locked_environment()
    geo = engine()
    data, inputs = geo.load_geo()
    manifest = geo.split_manifest(data, inputs)
    layout = outer_layout(data, manifest)
    if len(layout) != 10 or settings['B'] != 200:
        raise ValueError('Production requires ten outer folds and B=200')
    binding = dict(code_sha256=code_hashes(), config_sha256=sha_file(args.config),
                   input_sha256=inputs, split_manifest_sha256=hashlib.sha256(json.dumps(
                       manifest, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest(),
                   package_versions=environment, python_version=platform.python_version())
    args.output.mkdir(parents=True, exist_ok=True)
    if args.command == 'plan':
        atomic_json(args.output / 'PLAN.json', dict(binding=binding, settings=settings,
                    n=len(data.y), genes=len(data.genes), outer_scopes=len(layout),
                    real_randomized_fits_performed=False, PPMI_accessed=False))
        return
    if args.command == 'run':
        if args.receipt is None or args.preflight is None:
            raise ValueError('Public amendment receipt and passing preflight report required')
        verify_receipt(json.loads(args.receipt.read_text()), binding, args.preflight)
    with controller_lock(args.output):
        runfile = args.output / 'RUN_BINDING.json'
        if runfile.exists() and json.loads(runfile.read_text()) != binding:
            raise ValueError('Existing run provenance differs; refuse reuse')
        atomic_json(runfile, binding)
        if args.command == 'identity':
            result = execute(data, manifest, settings, args.output, binding, [0], 1, identity=True)[0]
            comparisons = []
            for scope, _, _ in layout:
                original = args.original_runs / scope.replace('/', '_') / 'heldout_predictions.npz'
                if not original.exists():
                    raise ValueError(f'Original prediction file missing for {scope}')
                new_completion = json.loads((args.output / 'identity/COMPLETE.json').read_text())
                attempt = next(iter(new_completion['artifacts'])).split('/')[0]
                fresh = args.output / 'identity' / attempt / scope.replace('/', '_') / 'heldout_predictions.npz'
                old_decisions = json.loads((original.parent / 'decisions.json').read_text())
                new_decisions = json.loads((fresh.parent / 'decisions.json').read_text())
                if any(old_decisions[model] != new_decisions[model] for model in ('M4', 'M5')):
                    raise ValueError(f'Identity tuning candidates or selection differ: {scope}')
                with np.load(original, allow_pickle=False) as old, np.load(fresh, allow_pickle=False) as new:
                    if not all(np.array_equal(old[k], new[k]) for k in ('samples', 'y')):
                        raise ValueError('Identity prediction membership differs')
                    for model in ('M4', 'M5'):
                        error = float(np.max(np.abs(old[model] - new[model])))
                        auc_error = abs(roc_auc_score(old['y'], old[model]) - roc_auc_score(new['y'], new[model]))
                        if error > 1e-6 or auc_error > 1e-12:
                            raise ValueError(f'Identity probabilities differ: {scope} {model}: {error}')
                        comparisons.append(dict(scope=scope, model=model,
                            max_abs_probability_error=error, AUROC_error=auc_error,
                            all_tuning_candidates_and_selection_identical=True))
            atomic_json(args.output / 'IDENTITY_VERIFIED.json', dict(status='passed',
                        binding_sha256=object_hash(binding), comparisons=comparisons,
                        statistics=result['statistics'], real_randomized_fits_performed=False))
            return
        identity = json.loads((args.output / 'IDENTITY_VERIFIED.json').read_text())
        if identity['status'] != 'passed' or identity['binding_sha256'] != object_hash(binding):
            raise ValueError('Observed identity replay not verified for this binding')
        results = execute(data, manifest, settings, args.output, binding,
                          list(range(settings['B'])), settings['workers'])
        summary = {}
        for source, observed in identity['statistics'].items():
            summary[source] = {}
            for statistic, value in observed.items():
                null = np.array([r['statistics'][source][statistic] for r in results])
                p = (1 + int(np.sum(null >= value - settings['tail_tolerance']))) / (settings['B'] + 1)
                summary[source][statistic] = dict(observed=value, null_mean=float(null.mean()),
                    null_SD=float(null.std(ddof=1)), p_value=p,
                    Monte_Carlo_SE_approx=float(np.sqrt(p * (1-p) / settings['B'])))
        atomic_json(args.output / 'SUMMARY.json', dict(status='complete', B=settings['B'],
                    results=summary, PPMI_accessed=False, exploratory_only=True))


if __name__ == '__main__':
    main()
