"""Fit the preregistered final auxiliary sex model on pooled GEO controls only."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import time
import warnings

import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geo_pipeline as geo
from geo_controls import FeaturePrep


def progress(done, max_fits, start, detail, finished=False):
    fraction = 1.0 if finished else done / max_fits
    width = 24
    bar = '#' * round(width * fraction) + '-' * (width - round(width * fraction))
    elapsed = (time.monotonic() - start) / 60
    eta = elapsed * (max_fits - done) / done if done else 0
    label = 'complete' if finished else f'up to {max_fits} fits'
    print(f"Progress [{bar}] {fraction:5.1%} | {done} fits, {label}"
          f" | elapsed {elapsed:.1f}m | ETA upper ~{eta:.1f}m | {detail}", flush=True)


def sex_folds(controls):
    if len(np.unique(controls.donor)) != len(controls.donor):
        raise ValueError('Repeated control donors')
    pseudo = geo.Data(controls.X, controls.age, controls.sex, controls.sex,
                      controls.source, controls.donor, controls.sample,
                      controls.genes, controls.marker_indices)
    folds = geo.donor_folds(pseudo, 5, 'controls/pooled/sex/inner')
    for train, validation in folds:
        if set(controls.donor[train]) & set(controls.donor[validation]):
            raise ValueError('Inner donor leakage')
        if len(np.unique(controls.sex[validation])) != 2:
            raise ValueError('Validation contains only one sex')
    return folds


def ordered_markers(data):
    path = geo.ROOT / 'protocol/v3_marker_manifest.json'
    items = json.loads(path.read_text())['sex_markers']
    ids = [str(item['gene_id']) for item in items]
    if len(ids) != 5 or len(set(ids)) != 5:
        raise ValueError('Expected five unique registered sex markers')
    if len(set(map(str, data.genes))) != len(data.genes):
        raise ValueError('Nonunique common gene IDs')
    lookup = {str(g): i for i, g in enumerate(data.genes)}
    if any(g not in lookup for g in ids):
        raise ValueError('Missing registered sex marker')
    return np.array([lookup[g] for g in ids], dtype=int), np.array(ids)


def fit_final(data, eligible, grid, out, input_hashes):
    controls = data.take(np.flatnonzero(data.y == 0))
    if set(np.unique(controls.y)) != {0}:
        raise ValueError('Noncontrol record entered auxiliary fit')
    folds = sex_folds(controls)
    if not all(len(np.unique(controls.source[a])) == 2 for a, _ in folds):
        raise ValueError('Each training fold must contain both GEO sources')
    if not np.all(np.isfinite(controls.X[:, eligible])):
        raise ValueError('Nonfinite sex marker input')
    prep_blocks = [(FeaturePrep(controls.take(a), eligible), controls.take(a),
                    controls.take(b)) for a, b in folds]
    initial_grid = sorted({float(x) for x in grid})
    if initial_grid != [0.0001, 0.001, 0.01, 0.1, 1.0]:
        raise ValueError('Registered sex lambda grid changed')
    lambdas = list(initial_grid)
    candidates = []
    warnings_seen = []
    start = time.monotonic()
    done, max_fits = 0, (len(initial_grid) + 2) * len(folds) + 1
    for round_i in range(3):
        pending = [lam for lam in sorted(lambdas, reverse=True)
                   if not any(c['lambda_value'] == lam for c in candidates)]
        for lam in pending:
            scores = []
            for fold_i, (prep, train, validation) in enumerate(prep_blocks):
                model, _ = geo.fit_checked(prep.transform(train), train.sex, 0,
                    lam, 0., sample_weight=np.ones(len(train.y)),
                    context=dict(scope='controls/pooled/sex', stage='inner',
                                 fold=fold_i, candidate_lambda=lam,
                                 preprocessing_sha256=prep.fingerprint()),
                    log_path=out / 'sex_fits.jsonl')
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    prediction = np.asarray(model.predict(prep.transform(validation)))
                if not np.all(np.isfinite(prediction)):
                    raise ValueError('Nonfinite inner sex prediction')
                scores.append(float(roc_auc_score(validation.sex, prediction)))
                warnings_seen.extend(dict(fold=fold_i, lambda_value=lam,
                                          warning=str(item.message)) for item in caught)
                done += 1
            candidates.append(dict(lambda_value=lam, score=float(np.mean(scores)),
                                   fold_AUROC=scores))
            progress(done, max_fits, start, f'lambda={lam:g} complete across 5 inner folds')
        best = max(c['score'] for c in candidates)
        selected = max((c for c in candidates if c['score'] >= best - 1e-6),
                       key=lambda c: c['lambda_value'])
        edge = selected['lambda_value'] in (min(lambdas), max(lambdas))
        if not edge or round_i == 2:
            break
        if selected['lambda_value'] == min(lambdas):
            lambdas.append(min(lambdas) / 10)
        else:
            lambdas.append(max(lambdas) * 10)
        print(f'Boundary extension {round_i + 1}/2; next lambda={lambdas[-1]:g}', flush=True)
    prep = FeaturePrep(controls, eligible)
    final_model, final_records = geo.fit_checked(
        prep.transform(controls), controls.sex, 0, selected['lambda_value'], 0.,
        sample_weight=np.ones(len(controls.y)),
        context=dict(scope='controls/pooled/sex', stage='final', fold=None,
                     preprocessing_sha256=prep.fingerprint()),
        log_path=out / 'sex_fits.jsonl')
    done += 1
    progress(done, max_fits, start, 'pooled final fit complete', finished=True)
    artifact = out / 'sex_model.npz'
    np.savez_compressed(artifact, selected=prep.selected,
        gene_ids=controls.genes[prep.selected], mean=prep.mean, sd=prep.sd,
        intercept=np.asarray(final_model.intercept_), coef=np.asarray(final_model.coef_))
    with np.load(artifact) as saved:
        reloaded_linear = saved['intercept'] + prep.transform(controls) @ saved['coef']
    original_linear = final_model.intercept_ + prep.transform(controls) @ final_model.coef_
    if not np.allclose(reloaded_linear, original_linear, rtol=0, atol=1e-12):
        raise ValueError('Artifact reload parity failed')
    split_data = [dict(train=controls.sample[a].tolist(),
                       validation=controls.sample[b].tolist()) for a, b in folds]
    geo.write_json(out / 'private_inner_folds.json', split_data)
    fit_records = [json.loads(line) for line in (out / 'sex_fits.jsonl').read_text().splitlines()]
    if not all(record['passed'] for record in fit_records):
        raise ValueError('A sex fit failed KKT validation')
    summary = dict(status='complete', PPMI_loaded=False, n_controls=len(controls.y),
        source_counts={str(source): int(np.sum(controls.source == source))
                       for source in np.unique(controls.source)},
        sex_counts={str(sex): int(np.sum(controls.sex == sex)) for sex in [0, 1]},
        selected_lambda=selected['lambda_value'], candidates=candidates,
        expansion_rounds=round_i, persistent_boundary=bool(edge),
        prediction_warnings=warnings_seen, fit_count=len(fit_records),
        fit_warning_count=sum(len(record['warnings']) for record in fit_records),
        max_KKT=max(record['diagnostics']['KKT_max'] for record in fit_records),
        final_KKT=final_records[-1]['diagnostics']['KKT_max'],
        input_sha256=input_hashes,
        split_manifest_sha256=geo.digest(out / 'private_inner_folds.json'),
        artifact_sha256=geo.digest(artifact),
        code_sha256=geo.digest(Path(__file__)),
        core_sha256=geo.digest(Path(geo.__file__)),
        decision_sha256=geo.digest(geo.ROOT / 'docs/GEO_FINAL_SEX_MODEL_DECISION_20261020_EN.md'),
        marker_manifest_sha256=geo.digest(geo.ROOT / 'protocol/v3_marker_manifest.json'),
        config_sha256=geo.digest(geo.ROOT / 'protocol/config_v3_1.json'),
        elapsed_minutes=(time.monotonic() - start) / 60)
    geo.write_json(out / 'SUMMARY.json', summary)
    return summary


def self_test():
    rng = np.random.default_rng(20261020)
    n = 80
    source = np.repeat(np.array(['GSE57475', 'GSE99039']), n // 2)
    sex = np.tile(np.repeat([0, 1], n // 4), 2)
    X = rng.normal(size=(n, 5)) + 0.8 * sex[:, None]
    genes = np.array([item['gene_id'] for item in
                      json.loads((geo.ROOT / 'protocol/v3_marker_manifest.json').read_text())['sex_markers']])
    data = geo.Data(X, np.full(n, 60.), sex, np.zeros(n, dtype=int), source,
                    np.array([f'donor-{i:03d}' for i in range(n)]),
                    np.array([f'sample-{i:03d}' for i in range(n)]), genes, {})
    eligible, ordered = ordered_markers(data)
    if not np.array_equal(genes, ordered) or len(sex_folds(data)) != 5:
        raise ValueError('Synthetic marker/fold check failed')
    with tempfile.TemporaryDirectory() as temp:
        out = Path(temp)
        summary = fit_final(data, eligible, [0.0001, 0.001, 0.01, 0.1, 1],
                            out, {'synthetic': True})
        if summary['status'] != 'complete' or not (out / 'sex_model.npz').exists():
            raise ValueError('Synthetic artifact check failed')
    print('SELF_TEST_PASSED', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.output is None:
        parser.error('--output is required except for --self-test')
    out = args.output.resolve()
    if geo.ROOT / 'private_data' not in out.parents:
        raise ValueError('GEO control artifacts must be in ignored private_data')
    out.mkdir(parents=True, exist_ok=False)
    data, hashes = geo.load_geo()
    audit = json.loads((geo.ROOT / 'data/manifests/input_audit.json').read_text())
    if hashes != audit['input_sha256']:
        raise ValueError('Prepared GEO inputs differ from frozen audit')
    eligible, ids = ordered_markers(data)
    config = json.loads((geo.ROOT / 'protocol/config_v3_1.json').read_text())
    symbols = [item['symbol'] for item in
               json.loads((geo.ROOT / 'protocol/v3_marker_manifest.json').read_text())['sex_markers']]
    if symbols != config['sex_positive_control']['markers']:
        raise ValueError('Manifest and config sex markers differ')
    summary = fit_final(data, eligible, config['sex_positive_control']['lambda_grid'],
                        out, hashes)
    print('GEO_FINAL_SEX_MODEL_COMPLETE:', out, flush=True)
    print('selected_lambda:', summary['selected_lambda'],
          'artifact_sha256:', summary['artifact_sha256'], flush=True)


if __name__ == '__main__':
    main()
