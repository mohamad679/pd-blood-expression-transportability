"""Paired synthetic source-prior stress test. Never reads participant files."""

import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import argparse
import json
from pathlib import Path
import time
import warnings

import numpy as np
from scipy.special import expit
from scipy.stats import norm, rankdata
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from block_logistic_v3_1 import balanced_source_weights, fit_block_logistic, kish_effective_n


LAMBDAS=(0.0001,0.001,0.01,0.1,1.)
RHOS=(0.,0.5,1.)


def draw_data(rng):
    source=np.r_[np.zeros(142,dtype=int),np.ones(239,dtype=int)]
    y=np.empty(381,dtype=int)
    y[:142]=rng.permutation(np.r_[np.ones(93,dtype=int),np.zeros(49,dtype=int)])
    y[142:]=rng.permutation(np.r_[np.ones(120,dtype=int),np.zeros(119,dtype=int)])
    age=rng.normal(size=381); sex=rng.binomial(1,0.5,381)
    baseline=np.column_stack((age,age**2,sex,age*sex,rng.normal(size=(381,6))))
    platform=(1-2*source)[:,None]*np.r_[np.full(20,2.),np.full(20,-2.)][None,:]
    genes=rng.normal(size=(381,40))+platform
    genes=(rankdata(genes,axis=1,method='average')-0.5)/40
    return np.column_stack((baseline,genes)),y,source


def standardize(train,test):
    mean=train.mean(axis=0); sd=train.std(axis=0)
    if np.any(sd<=0): raise ValueError('Synthetic fold has constant features')
    return (train-mean)/sd,(test-mean)/sd


def macro_auc(y,pred,source):
    return float(np.mean([roc_auc_score(y[source==s],pred[source==s]) for s in np.unique(source)]))


def select_and_evaluate(X,y,source,XT,yT,sourceT,folds,corrected):
    cached=[]
    for train,valid in folds:
        xtr,xval=standardize(X[train],X[valid])
        weights=balanced_source_weights(y[train],source[train]) if corrected else None
        cached.append((xtr,xval,train,valid,weights))
    score_cache={}

    def candidate_score(candidate):
        if candidate in score_cache: return score_cache[candidate]
        scores=[]
        for xtr,xval,tr,va,weights in cached:
            if candidate is None:
                model=fit_block_logistic(xtr[:,:10],y[tr],10,0.,0.,weights)
                pred=model.predict(np.ascontiguousarray(xval[:,:10]))
            else:
                lam,rho=candidate
                model=fit_block_logistic(xtr,y[tr],10,lam,rho,weights)
                pred=model.predict(xval)
            scores.append(macro_auc(y[va],pred,source[va]) if corrected else float(roc_auc_score(y[va],pred)))
        score_cache[candidate]=float(np.mean(scores))
        return score_cache[candidate]

    grid=list(LAMBDAS); expansions=0
    while True:
        candidates=[None]+[(lam,rho) for lam in sorted(grid,reverse=True) for rho in sorted(RHOS,reverse=True)]
        selected=None; best=candidate_score(None)
        for candidate in candidates[1:]:
            score=candidate_score(candidate)
            # Candidate order implements baseline, then larger lambda/rho ties.
            if score>best+1e-6: selected,best=candidate,score
        if selected is None:
            break
        lam,rho=selected
        if expansions>=2 or lam not in (min(grid),max(grid)):
            break
        grid.append(min(grid)/10 if lam==min(grid) else max(grid)*10)
        expansions+=1
    xtr,xt=standardize(X,XT)
    weights=balanced_source_weights(y,source) if corrected else None
    if selected is None:
        model=fit_block_logistic(xtr[:,:10],y,10,0.,0.,weights)
        pred=model.predict(np.ascontiguousarray(xt[:,:10])); nonzero=0
    else:
        lam,rho=selected
        model=fit_block_logistic(xtr,y,10,lam,rho,weights)
        nonzero=int(np.sum(np.abs(model.coef_[10:])>1e-8))
        pred=model.predict(xt)
    boundary=bool(selected is not None and nonzero>0 and selected[0] in (min(grid),max(grid)))
    return {'gene_selected':bool(selected is not None and nonzero>0),'selected':selected,
            'inner_selection_score':best,'independent_test_pooled_AUC':float(roc_auc_score(yT,pred)),
            'independent_test_macro_AUC':macro_auc(yT,pred,sourceT),'grid_expansions':expansions,
            'final_grid':sorted(grid),'selected_at_remaining_boundary':boundary,
            'nonzero_genes':nonzero,'train_kish':kish_effective_n(weights) if corrected else float(len(y))}


def wilson(successes,n):
    z=norm.ppf(0.975); phat=successes/n; den=1+z*z/n
    center=(phat+z*z/(2*n))/den
    half=z*np.sqrt(phat*(1-phat)/n+z*z/(4*n*n))/den
    return [float(center-half),float(center+half)]


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--repetitions',type=int,default=200)
    args=parser.parse_args(); repetitions=args.repetitions
    output=Path(__file__).with_name('source_null_v3_1.json')
    rows=[]; caught=[]; start=time.monotonic()
    for repeat in range(repetitions):
        rng=np.random.default_rng(np.random.SeedSequence([20261016,repeat]))
        X,y,s=draw_data(rng); XT,yT,sT=draw_data(rng)
        folds=list(StratifiedKFold(5,shuffle=True,random_state=20261016+repeat).split(X,s*2+y))
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter('always')
            old=select_and_evaluate(X,y,s,XT,yT,sT,folds,False)
            new=select_and_evaluate(X,y,s,XT,yT,sT,folds,True)
        caught.extend({'repeat':repeat,'message':str(m.message)} for m in messages)
        rows.append({'repeat':repeat,'unweighted_pooled_tuning':old,'weighted_macro_tuning':new})
        if (repeat+1)%10==0:
            print(f'Completed {repeat+1}/{repetitions} paired synthetic repetitions',flush=True)
            output.with_suffix('.inprogress.json').write_text(json.dumps({'completed':len(rows),'repetitions':repetitions,'rows':rows,'warnings':caught},indent=2)+'\n')
    summary={}
    for arm in ('unweighted_pooled_tuning','weighted_macro_tuning'):
        chosen=np.array([r[arm]['gene_selected'] for r in rows],dtype=float)
        metrics={}
        for metric in ('independent_test_pooled_AUC','independent_test_macro_AUC'):
            draws=np.array([r[arm][metric] for r in rows])
            metrics[metric]={'mean':float(draws.mean()),'MCSE':float(draws.std(ddof=1)/np.sqrt(repetitions))}
        summary[arm]={'gene_selection_count':int(chosen.sum()),'gene_selection_rate':float(chosen.mean()),
                      'gene_selection_Wilson_CI':wilson(int(chosen.sum()),repetitions),
                      'final_boundary_count':sum(r[arm]['selected_at_remaining_boundary'] for r in rows),
                      'expanded_repetitions':sum(r[arm]['grid_expansions']>0 for r in rows),**metrics}
    paired=np.array([int(r['unweighted_pooled_tuning']['gene_selected'])-int(r['weighted_macro_tuning']['gene_selected']) for r in rows])
    se=float(paired.std(ddof=1)/np.sqrt(repetitions)); diff=float(paired.mean())
    ci=[diff-norm.ppf(0.975)*se,diff+norm.ppf(0.975)*se]
    reasons=[]
    if ci[0]<=0: reasons.append('Paired gene-selection-rate reduction not established; review diagnostic numbers before publication.')
    if caught: reasons.append('Model warnings require review; no failed fit silently dropped.')
    for arm,stats in summary.items():
        m=stats['independent_test_macro_AUC']
        if abs(m['mean']-0.5)>3*m['MCSE']+0.01:
            reasons.append(arm+': independent null macro-AUC unexpectedly far from chance.')
    result={'status':'review_required' if reasons else 'diagnostics_as_expected','synthetic_only':True,
            'participant_data_read':False,'repetitions':repetitions,'seed':20261016,'n_train':381,'n_test':381,
            'genes':40,'baseline_features':10,'folds':5,'lambda_grid_start':LAMBDAS,'rho_grid':RHOS,
            'max_boundary_expansion_rounds':2,'summary':summary,'paired_selection_rate_reduction':diff,
            'paired_selection_rate_reduction_MC_CI':ci,'review_reasons':reasons,
            'warnings':caught,'rows':rows,'elapsed_seconds':round(time.monotonic()-start,2),
            'limitations':['One source-prior stress scenario; not a forecast of clinical model performance.',
                          'Synthetic baseline covariates do not contain source information; real marker scores may do so.',
                          'Gene-selection rate is not a Type-I-error rate.','Weights and macro tuning change together; their separate effects are not identified.']}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','summary','paired_selection_rate_reduction','paired_selection_rate_reduction_MC_CI','review_reasons','elapsed_seconds']},indent=2))


if __name__=='__main__':
    main()
