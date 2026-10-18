"""Reconcile saved fit membership/preprocessing hashes against frozen GEO splits.

No model fitting, label permutation or PPMI access. Only aggregate evidence is
written; donor/sample IDs, matrices and predictions remain private.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

try:
    from . import geo_pipeline as geo
except ImportError:
    import geo_pipeline as geo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runs', type=Path, default=geo.ROOT/'private_data/geo_runs_20261018')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    data, hashes = geo.load_geo()
    manifest = geo.split_manifest(data, hashes)
    lookup = {str(sample): i for i, sample in enumerate(data.sample)}
    summaries = []
    with threadpool_limits(limits=1):
        for item in manifest['splits']:
            scope = item['scope']
            out = args.runs / scope.replace('/', '_')
            started = json.loads((out/'STARTED.json').read_text())
            if started['input_sha256'] != hashes:
                raise ValueError('Original scope input hashes differ')
            tr = np.array([lookup[s] for s in item['train']], int)
            te = np.array([lookup[s] for s in item['test']], int)
            if set(data.donor[tr]) & set(data.donor[te]):
                raise ValueError('Outer operational-donor overlap')
            contexts = {}
            final = data.take(tr)
            prep = geo.Prep(final)
            contexts[('outer_or_final_refit', None)] = dict(
                n=len(tr), train_hash=geo.array_hash(final.sample), fingerprint=prep.fingerprint(),
                Kish=geo.kish_effective_n(geo.balanced_source_weights(final.y, final.source)))
            with np.load(out/'preprocessing.npz', allow_pickle=False) as saved:
                if str(saved['fit_samples_hash']) != geo.array_hash(final.sample):
                    raise ValueError('Saved final preprocessing membership differs')
                for name in ('base_mean', 'base_sd', 'keep', 'selected', 'gene_mean', 'gene_sd'):
                    if not np.array_equal(saved[name], np.asarray(getattr(prep, name))):
                        raise ValueError(f'Saved train-only preprocessing differs: {scope} {name}')
            for f, inner in enumerate(item['inner']):
                a = np.array([lookup[s] for s in inner['train']], int)
                b = np.array([lookup[s] for s in inner['validation']], int)
                if (set(a) & set(b) or set(a) | set(b) != set(tr)
                        or set(data.donor[a]) & set(data.donor[b])
                        or set(a) & set(te) or set(b) & set(te)):
                    raise ValueError('Invalid inner/outer membership')
                train = data.take(a)
                contexts[('inner', f)] = dict(n=len(a),
                    train_hash=geo.array_hash(train.sample), fingerprint=geo.Prep(train).fingerprint(),
                    Kish=geo.kish_effective_n(geo.balanced_source_weights(train.y, train.source)))
            records = [json.loads(line) for line in (out/'fits.jsonl').read_text().splitlines()]
            for record in records:
                expected = contexts[(record['stage'], record['fold'])]
                if (record['scope'] != scope or record['n'] != expected['n']
                        or record['train_sample_sha256'] != expected['train_hash']
                        or record['preprocessing_sha256'] != expected['fingerprint']
                        or not np.isclose(record['Kish'], expected['Kish'], atol=1e-10, rtol=0)
                        or not record['passed']):
                    raise ValueError('Logged fit state differs from its frozen training partition')
            if len(te):
                with np.load(out/'heldout_predictions.npz', allow_pickle=False) as saved:
                    if not np.array_equal(saved['samples'], data.sample[te]) or not np.array_equal(saved['y'], data.y[te]):
                        raise ValueError('Saved heldout membership/outcomes differ from frozen test')
            summaries.append(dict(scope=scope, n_train=len(tr), n_test=len(te),
                fit_attempts=len(records), verified_training_contexts=len(contexts),
                max_KKT=max(r['diagnostics']['KKT_max'] for r in records)))
            print('Audited', scope, flush=True)
    geo.write_json(args.output, dict(status='passed', scope_count=len(summaries),
        fit_attempts=sum(s['fit_attempts'] for s in summaries), scopes=summaries,
        input_sha256=hashes, new_model_fits=0, real_label_permutations=0, PPMI_accessed=False,
        limitations=['Operational donor separation only; PROBE uses GSM fallback.',
                     'Source prefixes do not establish biological cross-source donor independence.',
                     'This reconciles saved logs/state with locked inputs and split membership; it is not proof against every possible unlogged data-processing error.']))


if __name__ == '__main__':
    main()
