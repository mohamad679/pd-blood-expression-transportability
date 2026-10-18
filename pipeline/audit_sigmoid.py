"""Replay locked GEO inner fits to audit prediction clipping, without replacing results."""
import argparse
from collections import defaultdict
import inspect
import json
from pathlib import Path
import sys
import time
import warnings

import numpy as np
from scipy.special import expit
from sklearn.metrics import roc_auc_score
from glum._link import LogitLink

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geo_pipeline as geo


def probability_diagnostic(y, linear, clipped, source):
    y, linear, clipped = map(np.asarray, (y, linear, clipped))
    probability = expit(linear)
    eps50 = 50 * np.finfo(probability.dtype).eps
    expected = np.clip(probability, eps50, 1-eps50)
    if not np.array_equal(expected, clipped):
        raise ValueError('Installed glum prediction differs from documented clipping rule')
    aucs = defaultdict(list)
    mixed_tail_pairs = 0
    for group in np.unique(source):
        ix = source == group
        for name, score in [('clipped', clipped), ('expit', probability), ('linear', linear)]:
            aucs[name].append(float(roc_auc_score(y[ix], score[ix])))
        for tail in [probability < eps50, probability > 1-eps50]:
            mixed_tail_pairs += int(np.sum(ix & tail & (y == 0)) * np.sum(ix & tail & (y == 1)))
    return dict(n=len(y), clipped_low=int(np.sum(probability < eps50)),
                clipped_high=int(np.sum(probability > 1-eps50)),
                mixed_class_tail_pairs=mixed_tail_pairs,
                AUROC={name: float(np.mean(scores)) for name, scores in aucs.items()})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--pilot', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--aggregate', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    if geo.ROOT / 'private_data' not in out.parents:
        raise ValueError('Detailed replay outputs must remain in ignored private_data')
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    data, hashes = geo.load_geo()
    if hashes != json.loads((geo.ROOT / 'data/manifests/input_audit.json').read_text())['input_sha256']:
        raise ValueError('Frozen input hashes differ')
    manifest = json.loads((geo.ROOT / 'data/manifests/splits.json').read_text())
    context_events, comparisons, fit_count = [], [], 0
    max_score_error, max_KKT = 0., 0.
    for scope in manifest['splits']:
        name = scope['scope']
        original = args.pilot if name == 'GSE57475/outer/0' else args.runs / name.replace('/', '_')
        decisions = json.loads((original / 'decisions.json').read_text())
        historical = [json.loads(line) for line in (original / 'fits.jsonl').read_text().splitlines()]
        records = [r for r in historical if r['stage'] == 'inner' and r['passed']]
        if len(records) != len({(r['fold'], r['model'], r['lambda'], r['rho']) for r in records}):
            raise ValueError('Repeated successful fit keys require explicit review')
        dest = out / name.replace('/', '_')
        dest.mkdir()
        scores = defaultdict(lambda: defaultdict(list))
        for f, split in enumerate(scope['inner']):
            a = data.take(np.flatnonzero(np.isin(data.sample, split['train'])))
            b = data.take(np.flatnonzero(np.isin(data.sample, split['validation'])))
            if set(a.donor) & set(b.donor):
                raise ValueError('Replay donor leakage')
            prep = geo.Prep(a, 1000)
            for old in [r for r in records if r['fold'] == f]:
                if prep.fingerprint() != old['preprocessing_sha256']:
                    raise ValueError('Replay preprocessing differs from frozen fit')
                context = dict(scope=name, stage='inner', fold=f, model=old['model'])
                model, checked = geo.checked_fit(prep, a, old['model'], old['lambda'],
                    old['rho'], context, dest / 'fits.jsonl')
                max_KKT = max(max_KKT, checked[-1]['diagnostics']['KKT_max'])
                X, _ = prep.transform(b, old['model'])
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    pred = model.predict(X)
                linear = model.linear_predictor(X)
                diagnostic = probability_diagnostic(b.y, linear, pred, b.source)
                event = dict(**context, **{'lambda': old['lambda']}, rho=old['rho'],
                    warnings=[str(w.message) for w in caught], **diagnostic)
                with (dest / 'prediction_audit.jsonl').open('a') as stream:
                    stream.write(json.dumps(event) + '\n')
                if caught or diagnostic['clipped_low'] or diagnostic['clipped_high']:
                    context_events.append(event)
                key = (old['model'], old['lambda'], old['rho'])
                for kind, value in diagnostic['AUROC'].items():
                    scores[key][kind].append(value)
                fit_count += 1
                if fit_count % 100 == 0:
                    print('Replayed', fit_count, 'inner fits;', name, flush=True)
        for model_name, decision in decisions.items():
            candidates = decision.get('candidates')
            if candidates is None:
                expected = decision['inner_macro_AUROC']
                replayed = float(np.mean(scores[(model_name, 0., 0.)]['clipped']))
                max_score_error = max(max_score_error, abs(expected-replayed))
                continue
            alternatives = {kind: [] for kind in ['clipped', 'expit', 'linear']}
            for candidate in candidates:
                base = 'M1' if model_name == 'M3' else 'M4'
                key = ((base, 0., 0.) if candidate['baseline']
                       else (model_name, candidate['lambda'], candidate['rho']))
                for kind in alternatives:
                    if len(scores[key][kind]) != 5:
                        raise ValueError('Missing replay fold score')
                    item = dict(candidate, score=float(np.mean(scores[key][kind])))
                    alternatives[kind].append(item)
                    if kind == 'clipped':
                        max_score_error = max(max_score_error, abs(item['score']-candidate['score']))
            selections = {kind: geo.choose(rows) for kind, rows in alternatives.items()}
            def identity(c):
                return c['baseline'], c['lambda'], c['rho']
            original_choice = decision['selected']
            comparisons.append(dict(scope=name, model=model_name,
                original_selected=original_choice,
                replay_selected=selections,
                original_choice_reproduced=identity(original_choice) == identity(selections['clipped']),
                expit_choice_changed=identity(selections['expit']) != identity(selections['clipped']),
                linear_choice_changed=identity(selections['linear']) != identity(selections['clipped']),
                max_candidate_AUROC_change_expit=max(abs(a['score']-b['score']) for a,b in
                    zip(alternatives['clipped'], alternatives['expit'])),
                max_candidate_AUROC_change_linear=max(abs(a['score']-b['score']) for a,b in
                    zip(alternatives['clipped'], alternatives['linear']))))
        print('Replay completed:', name, 'total inner fits:', fit_count, flush=True)
        geo.write_json(out / 'progress.json', dict(completed_scope=name, inner_fits=fit_count,
            seconds=time.time()-start))
    summary = dict(created_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        purpose='Retrospective numerical diagnostic; original results and pipeline unchanged',
        PPMI_loaded=False, code_sha256=geo.digest(Path(__file__)),
        locked_pipeline_sha256=geo.digest(Path(geo.__file__)), input_sha256=hashes,
        split_manifest_sha256=geo.digest(geo.ROOT / 'data/manifests/splits.json'),
        glum_inverse_source=inspect.getsource(LogitLink.inverse),
        glum_link_file_sha256=geo.digest(inspect.getfile(LogitLink)),
        inner_fits_replayed=fit_count, scope_count=len(manifest['splits']),
        max_replay_KKT=max_KKT, max_original_candidate_score_error=max_score_error,
        all_original_choices_reproduced=all(c['original_choice_reproduced'] for c in comparisons),
        expit_choice_changes=sum(c['expit_choice_changed'] for c in comparisons),
        linear_choice_changes=sum(c['linear_choice_changed'] for c in comparisons),
        clipped_prediction_context_count=len(context_events), prediction_events=context_events,
        choice_comparisons=comparisons, seconds=time.time()-start)
    geo.write_json(args.aggregate, summary)
    if max_score_error > 1e-8 or not summary['all_original_choices_reproduced']:
        raise RuntimeError('Replay did not reproduce original tuning; review aggregate')
    print('Completed sigmoid audit:', args.aggregate, flush=True)


if __name__ == '__main__':
    main()
