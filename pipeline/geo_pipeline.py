"""Cold-fit, donor-separated GEO development pipeline for preregistration v3.1.

No PPMI loader or external-scoring command is provided. Prepared matrices are
samplewise transformed feasibility artifacts; every learned operation is fitted
again inside its training partition. Scientific settings come from frozen v3.1.
"""
import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys
import time
import platform

import numpy as np
from scipy.special import expit
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             log_loss, roc_auc_score)
from sklearn.model_selection import StratifiedKFold
from threadpoolctl import threadpool_info

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'protocol'))
from numerical_kkt_v3_1 import fit_checked
from block_logistic_v3_1 import balanced_source_weights, kish_effective_n

PANELS = ['T cells', 'B lineage', 'Monocytic lineage', 'Neutrophils',
          'Erythroid', 'Platelet']
MODELS = ('M1', 'M2', 'M3', 'M4', 'M5')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False),
                   encoding='utf-8')
    tmp.replace(path)


def array_hash(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.shape).encode()); h.update(a.dtype.str.encode())
        h.update(a.tobytes())
    return h.hexdigest()


def scope_seed(scope):
    return int.from_bytes(hashlib.sha256(
        ('20261014:' + scope).encode()).digest()[:4], 'big')


@dataclass
class Data:
    X: np.ndarray
    age: np.ndarray
    sex: np.ndarray
    y: np.ndarray
    source: np.ndarray
    donor: np.ndarray
    sample: np.ndarray
    genes: np.ndarray
    marker_indices: list

    def take(self, ix):
        return Data(*(getattr(self, k)[ix] for k in
                      ('X', 'age', 'sex', 'y', 'source', 'donor', 'sample')),
                    self.genes, self.marker_indices)


def load_geo():
    """Strict alignment, no models and no outcome-guided preprocessing."""
    manifest = json.loads((ROOT / 'protocol/v3_marker_manifest.json').read_text())
    parts, hashes = [], {}
    for accession in ('GSE57475', 'GSE99039'):
        p = ROOT / 'prepared' / (accession + '_development.npz')
        mp = ROOT / 'prepared' / (accession + '_sample_metadata.json')
        hashes[p.name], hashes[mp.name] = digest(p), digest(mp)
        z = np.load(p, allow_pickle=False)
        ids, genes = z['sample_ids'], z['gene_ids']
        if len(set(ids)) != len(ids) or len(set(genes)) != len(genes):
            raise ValueError('Duplicate sample or gene identifier')
        meta = json.loads(mp.read_text())
        by_id = {r['gsm']: r for r in meta}
        if len(by_id) != len(meta) or set(by_id) != set(ids):
            raise ValueError('Metadata/expression sample mismatch')
        rows = [by_id[s] for s in ids]
        age, sex, y = z['age'], z['sex'], z['y']
        for i, r in enumerate(rows):
            expected_sex = {1: 1, 0: 0, 'M': 1, 'Male': 1, 'F': 0, 'Female': 0}[r['sex']]
            if (float(r['age']) != age[i] or expected_sex != sex[i]
                    or int(r['case']) != y[i]):
                raise ValueError('Sample metadata not aligned')
        # Original ranks are exact half-integer midranks / G, saved as float32.
        # Recover that fixed lattice, rather than rerank rounded expression.
        raw = z['rank_features'].astype(float)
        lattice = raw * (2 * len(genes))
        if np.max(np.abs(lattice - np.rint(lattice))) >= .01:
            raise ValueError('Prepared ranks are not on the expected lattice')
        X = np.rint(lattice) / (2 * len(genes))
        if (not np.isfinite(X).all() or X.min() <= 0 or X.max() >= 1
                or not np.allclose(X.mean(axis=1), .5, atol=1e-12)):
            raise ValueError('Invalid samplewise rank features')
        donors = np.array([accession + ':' + str(r['subject_id'] or r['gsm'])
                           for r in rows])
        lookup = {g: i for i, g in enumerate(genes)}
        markers = [np.array([lookup[g] for g in
                    manifest['cell_marker_panels'][panel]['gene_ids']], int)
                   for panel in PANELS]
        parts.append(Data(X, age, sex, y, np.repeat(accession, len(y)), donors,
                          ids, genes, markers))
    if not np.array_equal(parts[0].genes, parts[1].genes):
        raise ValueError('Gene universes/order differ')
    pooled = Data(*(np.concatenate([getattr(d, k) for d in parts]) for k in
                    ('X', 'age', 'sex', 'y', 'source', 'donor', 'sample')),
                  parts[0].genes, parts[0].marker_indices)
    return pooled, hashes


def donor_folds(data, n_splits, scope):
    """Stratify unique donors by source and class, then expand to sample rows."""
    donors = np.unique(data.donor)
    strata = []
    for donor in donors:
        ix = np.flatnonzero(data.donor == donor)
        if len(np.unique(data.y[ix])) != 1 or len(np.unique(data.source[ix])) != 1:
            raise ValueError('Inconsistent donor labels/source')
        strata.append(str(data.source[ix[0]]) + ':' + str(int(data.y[ix[0]])))
    _, counts = np.unique(strata, return_counts=True)
    if min(counts) < n_splits:
        raise ValueError('Too few donors in a source/class stratum')
    result = []
    splitter = StratifiedKFold(n_splits, shuffle=True, random_state=scope_seed(scope))
    for train, val in splitter.split(donors, strata):
        ti = np.flatnonzero(np.isin(data.donor, donors[train]))
        vi = np.flatnonzero(np.isin(data.donor, donors[val]))
        assert not set(data.donor[ti]) & set(data.donor[vi])
        result.append((ti, vi))
    return result


def split_manifest(data, hashes):
    def entry(d, tr, te, scope):
        inner = donor_folds(d.take(tr), 5, scope + '/inner')
        return dict(scope=scope, train=d.sample[tr].tolist(),
                    test=d.sample[te].tolist(), inner=[
                        dict(train=d.sample[tr[a]].tolist(),
                             validation=d.sample[tr[b]].tolist()) for a, b in inner])
    items = []
    for source in np.unique(data.source):
        d = data.take(np.flatnonzero(data.source == source))
        for fold, (tr, te) in enumerate(donor_folds(d, 5, source + '/outer')):
            items.append(entry(d, tr, te, source + '/outer/' + str(fold)))
        tr = np.flatnonzero(data.source == source)
        te = np.flatnonzero(data.source != source)
        items.append(entry(data, tr, te, source + '/transfer'))
    items.append(entry(data, np.arange(len(data.y)), np.array([], int), 'pooled/final'))
    return dict(version='3.1.1', master_seed=20261014,
                seed_rule='first_4_bytes_big_endian_SHA256_of_20261014:scope',
                input_sha256=hashes, splits=items,
                identity_limit='PROBE uses GSM as operational donor; source prefixes do not prove cross-source independence')


class Prep:
    """Training-only unweighted scaling, baseline rank reduction, MAD selection."""
    def __init__(self, train, max_genes=1000):
        self.age_mean = float(train.age.mean())
        self.age_sd = float(train.age.std()) or 1.
        raw = self.baseline_raw(train)
        self.base_mean, self.base_sd = raw.mean(axis=0), raw.std(axis=0)
        safe = np.where(self.base_sd > 0, self.base_sd, 1.)
        z = (raw - self.base_mean) / safe
        self.keep = []
        design = np.ones((len(z), 1))
        for j in range(z.shape[1]):
            extended = np.column_stack((design, z[:, j]))
            if self.base_sd[j] > 0 and np.linalg.matrix_rank(extended) > design.shape[1]:
                self.keep.append(j); design = extended
        variability = []
        for source in np.unique(train.source):
            x = train.X[train.source == source]
            variability.append(np.median(np.abs(x - np.median(x, axis=0)), axis=0))
        mad = np.min(variability, axis=0)
        eligible = np.flatnonzero(mad > 0)
        order = np.lexsort((train.genes[eligible], -mad[eligible]))
        self.selected = eligible[order[:max_genes]]
        gene = train.X[:, self.selected]
        self.gene_mean, self.gene_sd = gene.mean(axis=0), gene.std(axis=0)
        if np.any(self.gene_sd <= 0):
            raise ValueError('Positive MAD selected constant gene')
        self.genes = train.genes
        self.fit_samples_hash = array_hash(train.sample)

    def baseline_raw(self, data):
        age = (data.age - self.age_mean) / self.age_sd
        markers = [data.X[:, ix].mean(axis=1) for ix in data.marker_indices]
        return np.column_stack((age, age**2, data.sex, age * data.sex, *markers))

    def transform(self, data, model):
        raw = self.baseline_raw(data)
        z = (raw - self.base_mean) / np.where(self.base_sd > 0, self.base_sd, 1.)
        keep = [i for i in self.keep if i < 4] if model in ('M1', 'M3') else self.keep
        baseline = z[:, keep] if model != 'M2' else np.empty((len(z), 0))
        gene = (data.X[:, self.selected] - self.gene_mean) / self.gene_sd
        if model in ('M1', 'M4'):
            return np.ascontiguousarray(baseline), baseline.shape[1]
        return np.ascontiguousarray(np.column_stack((baseline, gene))), baseline.shape[1]

    def fingerprint(self):
        return array_hash(np.array([self.age_mean, self.age_sd]), self.base_mean,
                          self.base_sd, np.array(self.keep), self.selected,
                          self.gene_mean, self.gene_sd)

    def save(self, path):
        np.savez_compressed(path, age_mean=self.age_mean, age_sd=self.age_sd,
            base_mean=self.base_mean, base_sd=self.base_sd, keep=self.keep,
            selected=self.selected, gene_mean=self.gene_mean, gene_sd=self.gene_sd,
            genes=self.genes, fit_samples_hash=self.fit_samples_hash)


def macro_auc(data, score):
    return float(np.mean([roc_auc_score(data.y[data.source == s], score[data.source == s])
                          for s in np.unique(data.source)]))


def metrics(y, p):
    return dict(n=len(y), cases=int(sum(y)), prevalence=float(np.mean(y)),
                AUROC=float(roc_auc_score(y, p)),
                AUPRC=float(average_precision_score(y, p)),
                Brier=float(brier_score_loss(y, p)),
                log_loss=float(log_loss(y, p, labels=[0, 1])))


def choose(candidates, tolerance=1e-6):
    best = max(c['score'] for c in candidates)
    tied = [c for c in candidates if c['score'] >= best - tolerance]
    return max(tied, key=lambda c: (c['baseline'], c['lambda'], c['rho']))


def checked_fit(prep, data, model, lam, rho, context, log_path):
    X, n_base = prep.transform(data, model)
    w = balanced_source_weights(data.y, data.source)
    fit, records = fit_checked(X, data.y, n_base, lam, rho, w,
        dict(**context, model=model, preprocessing_sha256=prep.fingerprint(),
             train_sample_sha256=array_hash(data.sample), Kish=kish_effective_n(w)), log_path)
    return fit, records


def train_suite(train, inner, out, scope, max_genes=1000,
                grid=(.0001, .001, .01, .1, 1.), rhos=(0., .5, 1.),
                expansion_rounds=2):
    """All candidates checked; no warm starts or silently discarded failures."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    log = out / 'fits.jsonl'
    # Explicitly reject append-to-old-run accidents.
    if log.exists():
        raise ValueError('Output already contains fits; use a fresh run directory')
    fold_data = []
    for f, (ti, vi) in enumerate(inner):
        a, b = train.take(ti), train.take(vi)
        if set(a.donor) & set(b.donor):
            raise ValueError('Inner donor leakage')
        fold_data.append((Prep(a, max_genes), a, b))
    final_prep = Prep(train, max_genes)
    final_prep.save(out / 'preprocessing.npz')
    final, decisions = {}, {}
    baseline_scores = {}
    for name in ('M1', 'M4'):
        scores = []
        for f, (prep, a, b) in enumerate(fold_data):
            fit, _ = checked_fit(prep, a, name, 0., 0.,
                dict(scope=scope, stage='inner', fold=f, baseline_only=True), log)
            X, _ = prep.transform(b, name)
            scores.append(macro_auc(b, fit.predict(X)))
        baseline_scores[name] = float(np.mean(scores))
        fit, _ = checked_fit(final_prep, train, name, 0., 0.,
            dict(scope=scope, stage='outer_or_final_refit', fold=None, baseline_only=True), log)
        final[name] = fit
        decisions[name] = dict(hyperparameters=None, inner_macro_AUROC=baseline_scores[name])
        print(scope, name, 'baseline complete', flush=True)
    for name in ('M2', 'M3', 'M5'):
        baseline_name = 'M1' if name == 'M3' else 'M4'
        candidates = []
        if name != 'M2':
            candidates.append(dict(baseline=True, **{'lambda': 0.}, rho=0.,
                                   score=baseline_scores[baseline_name], nonzero=False))
        lambdas = list(grid)
        for round_i in range(expansion_rounds + 1):
            existing = {(c['lambda'], c['rho']) for c in candidates if not c['baseline']}
            for lam in sorted(lambdas, reverse=True):
                for rho in rhos:
                    if (lam, rho) in existing:
                        continue
                    scores, nonzero = [], False
                    for f, (prep, a, b) in enumerate(fold_data):
                        fit, records = checked_fit(prep, a, name, lam, rho,
                            dict(scope=scope, stage='inner', fold=f, baseline_only=False,
                                 expansion_round=round_i), log)
                        X, _ = prep.transform(b, name)
                        scores.append(macro_auc(b, fit.predict(X)))
                        nonzero |= records[-1]['diagnostics']['nonzero_penalized_coefficients'] > 0
                    candidates.append(dict(baseline=False, **{'lambda': float(lam)},
                        rho=float(rho), score=float(np.mean(scores)), nonzero=bool(nonzero)))
                    write_json(out / (name + '_tuning.inprogress.json'), candidates)
            selected = choose(candidates)
            edge = (not selected['baseline'] and selected['nonzero']
                    and selected['lambda'] in (min(lambdas), max(lambdas)))
            if not edge or round_i == expansion_rounds:
                break
            lambdas.append(min(lambdas) / 10 if selected['lambda'] == min(lambdas)
                           else max(lambdas) * 10)
        chosen = dict(selected)
        if selected['baseline']:
            final[name] = final[baseline_name]
            chosen['reuses_baseline'] = baseline_name
        else:
            fit, records = checked_fit(final_prep, train, name, selected['lambda'], selected['rho'],
                dict(scope=scope, stage='outer_or_final_refit', fold=None, baseline_only=False), log)
            if name in ('M3', 'M5') and records[-1]['diagnostics']['nonzero_penalized_coefficients'] == 0:
                final[name] = final[baseline_name]
                chosen['reuses_baseline'] = baseline_name
                chosen['final_gene_block_zero'] = True
            else:
                final[name] = fit
        decisions[name] = dict(selected=chosen, candidates=candidates,
                               persistent_boundary=bool(edge), expansion_rounds=round_i)
        write_json(out / 'decisions.json', decisions)
        print(scope, name, 'complete', flush=True)
    for name, fit in final.items():
        np.savez_compressed(out / (name + '_model.npz'),
                            intercept=fit.intercept_, coef=fit.coef_)
    return final_prep, final, decisions


def predict_suite(prep, fits, decisions, test):
    p = {}
    for name in MODELS:
        actual = decisions.get(name, {}).get('selected', {}).get('reuses_baseline', name)
        if actual != name:
            p[name] = p[actual].copy()
        else:
            X, _ = prep.transform(test, name)
            p[name] = expit(fits[name].intercept_ + X @ fits[name].coef_)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('command', choices=['plan', 'pilot', 'analysis'])
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--registration-receipt', type=Path)
    ap.add_argument('--scope', help='Exact planned scope; analysis command only')
    args = ap.parse_args()
    data, hashes = load_geo()
    splits = split_manifest(data, hashes)
    args.output.mkdir(parents=True, exist_ok=True)
    if args.command == 'plan':
        write_json(args.output / 'splits.json', splits)
        write_json(args.output / 'input_audit.json', dict(
            input_sha256=hashes, n=len(data.y), genes=len(data.genes),
            sources={str(s): metrics(data.y[data.source == s],
                     np.full(sum(data.source == s), .5))['n'] for s in np.unique(data.source)},
            split_manifest_sha256=digest(args.output / 'splits.json'),
            real_training_performed=False))
        return
    if not args.registration_receipt:
        raise ValueError('Public pre-pilot implementation amendment receipt required')
    receipt = json.loads(args.registration_receipt.read_text())
    if (receipt.get('publicly_verified') is not True
            or receipt['pipeline_code_sha256'] != digest(__file__)
            or receipt['split_manifest_sha256'] != hashlib.sha256(
                json.dumps(splits, indent=2, ensure_ascii=False, allow_nan=False).encode()).hexdigest()):
        raise ValueError('Public lock receipt does not match code/splits')
    scope = 'GSE57475/outer/0' if args.command == 'pilot' else args.scope
    selected_split = [s for s in splits['splits'] if s['scope'] == scope]
    if len(selected_split) != 1:
        raise ValueError('An exact planned scope is required')
    s = selected_split[0]
    lookup = {sample: i for i, sample in enumerate(data.sample)}
    a = data.take(np.array([lookup[x] for x in s['train']], int))
    b = data.take(np.array([lookup[x] for x in s['test']], int))
    local = {sample: i for i, sample in enumerate(a.sample)}
    inner = [(np.array([local[x] for x in f['train']], int),
              np.array([local[x] for x in f['validation']], int)) for f in s['inner']]
    write_json(args.output / 'STARTED.json', dict(
        scope=scope, input_sha256=hashes, registration=receipt,
        utc_started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        n_train=len(a.y), n_test=len(b.y), PPMI_accessed=False,
        runtime=dict(python=sys.version, platform=platform.platform(),
                     thread_pools=threadpool_info(), concurrent_model_fits=1,
                     warm_start=False)))
    started = time.monotonic()
    cpu_started = time.process_time()
    try:
        prep, fits, decisions = train_suite(a, inner, args.output, scope)
        probabilities = predict_suite(prep, fits, decisions, b) if len(b.y) else {}
        if probabilities:
            np.savez_compressed(args.output / 'heldout_predictions.npz',
                                samples=b.sample, y=b.y, **probabilities)
        write_json(args.output / 'RESULTS.json', dict(status='complete', scope=scope,
            seconds=time.monotonic() - started,
            process_CPU_seconds=time.process_time() - cpu_started,
            metrics={m: metrics(b.y, p) for m, p in probabilities.items()},
            interpretation=('timing pilot; one outer fold, not a study performance estimate'
                            if args.command == 'pilot' else 'prespecified GEO development analysis'),
            PPMI_scores_computed=False))
    except Exception as e:
        write_json(args.output / 'FAILURE.json', dict(status='affected_analysis_invalid',
            error_type=type(e).__name__, message=str(e),
            seconds=time.monotonic() - started, PPMI_scores_computed=False))
        raise


if __name__ == '__main__':
    main()
