"""Registered sex/age transfer controls; GEO healthy controls only, no PPMI."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import warnings

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import balanced_accuracy_score, mean_absolute_error, roc_auc_score
from sklearn.model_selection import KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geo_pipeline as geo


class FeaturePrep:
    """Fixed columns or training-only positive-MAD autosomal selection/scaling."""
    def __init__(self, data, eligible, max_genes=None):
        eligible = np.asarray(eligible, int)
        if max_genes is not None:
            x = data.X[:, eligible]
            mad = np.median(np.abs(x - np.median(x, axis=0)), axis=0)
            good = np.flatnonzero(mad > 0)
            order = np.lexsort((data.genes[eligible[good]], -mad[good]))
            eligible = eligible[good[order[:max_genes]]]
        self.selected = eligible
        self.mean = data.X[:, eligible].mean(axis=0)
        self.sd = data.X[:, eligible].std(axis=0)
        if not len(eligible) or np.any(self.sd <= 0):
            raise ValueError('Empty or constant auxiliary feature block')

    def transform(self, data):
        return np.ascontiguousarray((data.X[:, self.selected] - self.mean) / self.sd)

    def fingerprint(self):
        return geo.array_hash(self.selected, self.mean, self.sd)


def control_folds(data, target):
    if len(np.unique(data.donor)) != len(data.donor):
        raise ValueError('Auxiliary controls require one record per donor')
    scope = 'controls/' + str(data.source[0]) + '/' + target + '/inner'
    if target == 'sex':
        sex_data = geo.Data(data.X, data.age, data.sex, data.sex, data.source,
                            data.donor, data.sample, data.genes, data.marker_indices)
        return geo.donor_folds(sex_data, 5, scope)
    donors = np.unique(data.donor)
    splitter = KFold(5, shuffle=True, random_state=geo.scope_seed(scope))
    return [(np.flatnonzero(np.isin(data.donor, donors[a])),
             np.flatnonzero(np.isin(data.donor, donors[b])))
            for a, b in splitter.split(donors)]


def bootstrap_control(target, truth, score, baseline, scope, repetitions=5000):
    rng = np.random.default_rng(geo.scope_seed(scope))
    truth, score = np.asarray(truth), np.asarray(score)
    values = []
    undefined = 0
    groups = [np.flatnonzero(truth == k) for k in (0, 1)] if target == 'sex' else None
    for _ in range(repetitions):
        ix = (np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
              if groups is not None else rng.integers(0, len(truth), len(truth)))
        if target == 'sex':
            values.append([roc_auc_score(truth[ix], score[ix])])
        else:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                corr = float(spearmanr(truth[ix], score[ix]).statistic)
            if not np.isfinite(corr):
                undefined += 1
            mae = mean_absolute_error(truth[ix], score[ix])
            base_mae = mean_absolute_error(truth[ix], np.repeat(baseline, len(ix)))
            values.append([corr, mae, mae - base_mae])
    values = np.asarray(values)
    names = ['AUROC'] if target == 'sex' else ['Spearman', 'MAE', 'MAE_minus_median_baseline']
    return dict(repetitions=repetitions, undefined_Spearman_replicates=undefined,
                percentile_95_CI={name: np.nanpercentile(values[:, j], [2.5, 97.5]).tolist()
                                  for j, name in enumerate(names)})


def fit_sex(train, test, eligible, out, grid, expansions=2):
    folds = control_folds(train, 'sex')
    blocks = [(FeaturePrep(train.take(a), eligible), train.take(a), train.take(b))
              for a, b in folds]
    candidates = []
    lambdas = list(grid)
    scope = 'controls/' + str(train.source[0]) + '/sex'
    for round_i in range(expansions + 1):
        for lam in sorted(lambdas, reverse=True):
            if any(c['lambda'] == lam for c in candidates):
                continue
            scores = []
            for f, (prep, a, b) in enumerate(blocks):
                model, _ = geo.fit_checked(prep.transform(a), a.sex, 0, lam, 0.,
                    sample_weight=np.ones(len(a.y)),
                    context=dict(scope=scope, stage='inner', fold=f,
                                 preprocessing_sha256=prep.fingerprint()),
                    log_path=out / 'sex_fits.jsonl')
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    pred = model.predict(prep.transform(b))
                scores.append(float(roc_auc_score(b.sex, pred)))
                if caught:
                    with (out / 'sex_prediction_warnings.jsonl').open('a') as stream:
                        stream.write(json.dumps(dict(fold=f, **{'lambda': lam},
                            warnings=[str(w.message) for w in caught])) + '\n')
            candidates.append({'lambda': float(lam), 'score': float(np.mean(scores))})
        best_score = max(c['score'] for c in candidates)
        selected = max((c for c in candidates if c['score'] >= best_score - 1e-6),
                       key=lambda c: c['lambda'])
        edge = selected['lambda'] in (min(lambdas), max(lambdas))
        if not edge or round_i == expansions:
            break
        lambdas.append(min(lambdas) / 10 if selected['lambda'] == min(lambdas)
                       else max(lambdas) * 10)
    prep = FeaturePrep(train, eligible)
    model, records = geo.fit_checked(prep.transform(train), train.sex, 0,
        selected['lambda'], 0., sample_weight=np.ones(len(train.y)),
        context=dict(scope=scope, stage='final', fold=None,
                     preprocessing_sha256=prep.fingerprint()), log_path=out / 'sex_fits.jsonl')
    Xtest = prep.transform(test)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        pred = model.predict(Xtest)
    linear = model.intercept_ + Xtest @ model.coef_
    auc = float(roc_auc_score(test.sex, pred))
    result = dict(selected=selected, candidates=candidates, persistent_boundary=bool(edge),
                  expansion_rounds=round_i, AUROC=auc,
                  linear_score_AUROC_diagnostic=float(roc_auc_score(test.sex, linear)),
                  BA_at_0_5_descriptive=float(balanced_accuracy_score(test.sex, pred >= .5)),
                  prediction_warnings=[str(w.message) for w in caught],
                  final_KKT=records[-1]['diagnostics']['KKT_max'],
                  development_direction_pass=bool(auc >= .95))
    result['uncertainty'] = bootstrap_control('sex', test.sex, pred, None, scope + '/bootstrap')
    np.savez_compressed(out / 'sex_model.npz', selected=prep.selected, mean=prep.mean,
        sd=prep.sd, intercept=model.intercept_, coef=model.coef_, gene_ids=train.genes[prep.selected])
    return result, pred, folds


def fit_age(train, test, eligible, out, alphas):
    folds = control_folds(train, 'age')
    blocks = [(FeaturePrep(train.take(a), eligible, 1000), train.take(a), train.take(b))
              for a, b in folds]
    candidates = []
    for alpha in alphas:
        scores = []
        for prep, a, b in blocks:
            model = Ridge(alpha=alpha, solver='svd', fit_intercept=True)
            model.fit(prep.transform(a), a.age)
            scores.append(float(mean_absolute_error(b.age, model.predict(prep.transform(b)))))
        candidates.append(dict(alpha=float(alpha), MAE=float(np.mean(scores))))
    best = min(c['MAE'] for c in candidates)
    selected = max((c for c in candidates if c['MAE'] <= best + 1e-6), key=lambda c: c['alpha'])
    prep = FeaturePrep(train, eligible, 1000)
    model = Ridge(alpha=selected['alpha'], solver='svd', fit_intercept=True)
    model.fit(prep.transform(train), train.age)
    pred = model.predict(prep.transform(test))
    median = float(np.median(train.age))
    corr = float(spearmanr(test.age, pred).statistic)
    mae = float(mean_absolute_error(test.age, pred))
    base = float(mean_absolute_error(test.age, np.repeat(median, len(test.y))))
    if not np.isfinite(pred).all() or not np.isfinite(corr):
        raise ValueError('Nonfinite age predictions/correlation')
    result = dict(selected=selected, candidates=candidates, Spearman=corr, MAE=mae,
                  training_median_age=median, MAE_train_median_predictor=base,
                  MAE_minus_median_baseline=mae-base, selected_genes=len(prep.selected),
                  development_direction_capacity=bool(corr >= .4 and mae < base))
    scope = 'controls/' + str(train.source[0]) + '/age/bootstrap'
    result['uncertainty'] = bootstrap_control('age', test.age, pred, median, scope)
    np.savez_compressed(out / 'age_model.npz', selected=prep.selected, mean=prep.mean,
        sd=prep.sd, intercept=model.intercept_, coef=model.coef_, gene_ids=train.genes[prep.selected])
    return result, pred, folds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True,
                        help='Fresh ignored private_data output directory')
    parser.add_argument('--aggregate', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    if geo.ROOT / 'private_data' not in out.parents:
        raise ValueError('Individual outputs must be inside ignored private_data')
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    data, hashes = geo.load_geo()
    audit = json.loads((geo.ROOT / 'data/manifests/input_audit.json').read_text())
    if hashes != audit['input_sha256']:
        raise ValueError('Prepared inputs differ from frozen audit')
    config = json.loads((geo.ROOT / 'protocol/config_v3_1.json').read_text())
    markers = json.loads((geo.ROOT / 'protocol/v3_marker_manifest.json').read_text())['sex_markers']
    lookup = {g: i for i, g in enumerate(data.genes)}
    sex_ix = [lookup[m['gene_id']] for m in markers]
    with np.load(geo.ROOT / 'prepared/GSE57475_development.npz') as z:
        chrom = z['chromosomes']
    with np.load(geo.ROOT / 'prepared/GSE99039_development.npz') as z:
        if not np.array_equal(chrom, z['chromosomes']):
            raise ValueError('Chromosome mapping differs across inputs')
    age_ix = np.flatnonzero(np.isin(chrom, ['chr' + str(i) for i in range(1, 23)]))
    if len(age_ix) != 15720:
        raise ValueError('Unexpected autosomal universe')
    controls = data.take(np.flatnonzero(data.y == 0))
    if len(np.unique(controls.donor)) != len(controls.y):
        raise ValueError('Repeated control donors')
    results, private_splits = {}, {}
    for source in np.unique(controls.source):
        train = controls.take(np.flatnonzero(controls.source == source))
        test = controls.take(np.flatnonzero(controls.source != source))
        if set(train.donor) & set(test.donor):
            raise ValueError('Cross-source donor leakage')
        dest = out / str(source)
        dest.mkdir()
        print('Control source:', source, 'n_train:', len(train.y), 'n_test:', len(test.y), flush=True)
        sex, sex_pred, sex_folds = fit_sex(train, test, sex_ix, dest,
            config['sex_positive_control']['lambda_grid'])
        age, age_pred, age_folds = fit_age(train, test, age_ix, dest,
            config['age_positive_control']['alpha_grid'])
        results[str(source)] = dict(test_source=str(test.source[0]), n_train=len(train.y),
            n_test=len(test.y), train_female=int(np.sum(train.sex == 0)),
            train_male=int(np.sum(train.sex == 1)), test_female=int(np.sum(test.sex == 0)),
            test_male=int(np.sum(test.sex == 1)), train_age_range=[float(train.age.min()), float(train.age.max())],
            test_age_range=[float(test.age.min()), float(test.age.max())], sex=sex, age=age)
        private_splits[str(source)] = {name: [dict(train=train.sample[a].tolist(),
            validation=train.sample[b].tolist()) for a, b in folds]
            for name, folds in [('sex', sex_folds), ('age', age_folds)]}
        np.savez_compressed(dest / 'heldout_predictions.npz', sample=test.sample,
            sex=test.sex, age=test.age, sex_prediction=sex_pred, age_prediction=age_pred)
        geo.write_json(dest / 'RESULTS.json', results[str(source)])
        print(source, 'sex AUROC:', sex['AUROC'], 'age Spearman:', age['Spearman'], flush=True)
    geo.write_json(out / 'private_splits.json', private_splits)
    fits = [json.loads(line) for p in out.glob('*/sex_fits.jsonl')
            for line in p.read_text().splitlines()]
    summary = dict(created_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        input_sha256=hashes, core_pipeline_sha256=geo.digest(Path(geo.__file__)),
        auxiliary_code_sha256=geo.digest(Path(__file__)),
        execution_plan_sha256=geo.digest(geo.ROOT / 'docs/GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md'),
        config_sha256=geo.digest(geo.ROOT / 'protocol/config_v3_1.json'),
        split_manifest_sha256=geo.digest(out / 'private_splits.json'),
        PPMI_loaded=False, uniform_auxiliary_weights=True, autosomal_genes=len(age_ix),
        directions=results, sex_fit_count=len(fits), sex_all_fits_passed=all(r['passed'] for r in fits),
        sex_max_KKT=max(r['diagnostics']['KKT_max'] for r in fits),
        sex_fit_warning_count=sum(len(r['warnings']) for r in fits),
        development_sex_gate_pass=all(r['sex']['development_direction_pass'] for r in results.values()),
        development_age_capacity_demonstrated=all(r['age']['development_direction_capacity'] for r in results.values()),
        seconds=time.time()-started)
    geo.write_json(args.aggregate, summary)
    print('Aggregate:', args.aggregate, flush=True)


if __name__ == '__main__':
    main()
