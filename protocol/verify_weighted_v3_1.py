"""Synthetic checks of weighted partial penalties; no participant data read."""

import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import importlib.metadata
import json
from pathlib import Path
import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.metrics import roc_auc_score

from block_logistic_v3 import fit_block_logistic as original_fit
from block_logistic_v3_1 import balanced_source_weights, fit_block_logistic, kish_effective_n
from verify_block_penalty_v3 import main as run_original_checks


def independent_reference(X, y, w, n_baseline, lam, rho):
    """L-BFGS-B with positive/negative gene parts, exact L1 epigraph form."""
    B = np.column_stack((np.ones(len(X)), X[:, :n_baseline]))
    G = X[:, n_baseline:]
    k, p = B.shape[1], G.shape[1]
    s = w / w.sum()

    def evaluate(theta):
        beta = theta[:k]
        positive, negative = theta[k:k+p], theta[k+p:]
        gamma = positive - negative
        eta = B @ beta + G @ gamma
        residual = s * (expit(eta)-y)
        objective = (np.dot(s, np.logaddexp(0, eta)-y*eta)
                     + lam*rho*np.sum(positive+negative)
                     + 0.5*lam*(1-rho)*np.dot(gamma,gamma))
        gene_gradient = G.T @ residual + lam*(1-rho)*gamma
        gradient = np.r_[B.T @ residual, gene_gradient+lam*rho, -gene_gradient+lam*rho]
        return objective, gradient

    start = np.zeros(k+2*p)
    result = minimize(evaluate, start, method="L-BFGS-B", jac=True,
                      bounds=[(None,None)]*k+[(0,None)]*(2*p),
                      options={"ftol":1e-14,"gtol":1e-9,"maxiter":10000,"maxls":100})
    reference = np.r_[result.x[:k], result.x[k:k+p]-result.x[k+p:]]
    return reference, float(result.fun), bool(result.success)


def main():
    run_original_checks()
    rng = np.random.default_rng(20261015)
    n,k,p = 1200,4,24
    B = rng.normal(size=(n,k))
    G = rng.normal(size=(n,p)); G[:,0] += 0.5*B[:,0]
    X = np.column_stack((B,G)); X = (X-X.mean(axis=0))/X.std(axis=0)
    y = rng.binomial(1,expit(-0.2 + X[:,:k] @ np.array([0.9,-0.7,0.5,0.25])+0.6*X[:,k])).astype(float)
    records=[]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        # Reproduce all 21 model configurations underlying the original 23 checks.
        configurations = [(k,a,r,True) for r in (0.,0.5,1.) for a in (0.,0.01,1.,100.)]
        configurations += [(k,a,r,False) for r in (0.,0.5,1.) for a in (0.01,0.2,1.)]
        configurations += [(0,0.2,0.5,False)]  # Original deliberately wrong-mask fit.
        for nb,lam,rho,baseline_only in configurations:
            data = np.ascontiguousarray(X[:,:k] if baseline_only else X)
            reference = original_fit(data,y,nb,lam,rho)
            reference_parameters = np.r_[reference.intercept_,reference.coef_]
            for constant in (0.1,1.,10.):
                model=fit_block_logistic(data,y,nb,lam,rho,np.full(n,constant))
                error=float(np.max(np.abs(np.r_[model.intercept_,model.coef_]-reference_parameters)))
                assert error<2e-5,("uniform weight scale invariance",lam,rho,constant,error)
                records.append({"test":"uniform_weights_reproduce_original","baseline_only":baseline_only,
                                "lambda":lam,"rho":rho,"constant_weight":constant,"max_coefficient_error":error})
        # Nonuniform normalized weighted likelihood versus an independent solver.
        w=np.exp(rng.normal(0,0.7,n))
        for rho in (0.,0.5,1.):
            for lam in (0.01,0.2,1.):
                model=fit_block_logistic(X,y,k,lam,rho,w)
                ref,obj,success=independent_reference(X,y,w,k,lam,rho)
                fitted=np.r_[model.intercept_,model.coef_]
                error=float(np.max(np.abs(fitted-ref)))
                eta=model.intercept_+X@model.coef_
                own_obj=float(np.dot(w/w.sum(),np.logaddexp(0,eta)-y*eta)
                              +lam*rho*np.sum(np.abs(model.coef_[k:]))
                              +0.5*lam*(1-rho)*np.dot(model.coef_[k:],model.coef_[k:]))
                grad=np.column_stack((np.ones(n),X[:,:k])).T@(w/w.sum()*(model.predict(X)-y))
                assert success and error<3e-5,("independent weighted solver",success,error)
                assert abs(own_obj-obj)<1e-9,("weighted objective",own_obj,obj)
                assert np.max(np.abs(grad))<2e-6,("weighted baseline gradient",grad)
                scaled=fit_block_logistic(X,y,k,lam,rho,7*w)
                scale_error=float(np.max(np.abs(np.r_[scaled.intercept_,scaled.coef_]-fitted)))
                assert scale_error<2e-5,("nonuniform rescaling",scale_error)
                records.append({"test":"independent_nonuniform_weighted_solver","lambda":lam,"rho":rho,
                                "max_coefficient_error":error,"objective_error":abs(own_obj-obj),
                                "baseline_gradient_max":float(np.max(np.abs(grad))),"rescale_error":scale_error})
        # Analytic source-prior example, both one-source and pooled weights.
        cohorts=np.r_[np.zeros(142,dtype=int),np.ones(239,dtype=int)]
        labels=np.r_[np.ones(93),np.zeros(49),np.ones(120),np.zeros(119)].astype(int)
        weights=balanced_source_weights(labels,cohorts)
        assert np.isclose(weights.sum(),381)
        for source in (0,1):
            for label in (0,1):
                assert np.isclose(weights[(cohorts==source)&(labels==label)].sum(),381/4)
        naive_auc=float(roc_auc_score(labels,1-cohorts))
        weighted_auc=float(roc_auc_score(labels,1-cohorts,sample_weight=weights))
        assert abs(naive_auc-1951/3408)<1e-12 and abs(weighted_auc-0.5)<1e-12
        eff={"pooled":{"n":381,"kish":kish_effective_n(weights),"cell_weights":{}}}
        for source,name in ((0,"PROBE"),(1,"GENEPARK_TRO")):
            ix=cohorts==source; sw=balanced_source_weights(labels[ix],cohorts[ix])
            eff[name]={"n":int(ix.sum()),"kish":kish_effective_n(sw)}
            for label in (0,1):
                eff['pooled']['cell_weights'][name+('_PD' if label else '_HC')]=float(weights[ix&(labels==label)][0])
        messages=[str(w.message) for w in captured]
        assert not messages,messages
    payload={"status":"passed","data":"synthetic_only","original_checks_reexecuted":23,
             "uniform_comparisons":66,"nonuniform_independent_comparisons":9,
             "glum":importlib.metadata.version('glum'),"checks":records,
             "cohort_only_AUC_unweighted":naive_auc,"cohort_only_AUC_weighted":weighted_auc,
             "kish_and_weights":eff,"warnings":messages,"participant_data_read":False,
             "normalization":"likelihood divided by sum(weights); global weight scaling leaves lambda unchanged"}
    out=Path(__file__).with_name('weighted_v3_1_verification.json')
    out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:payload[k] for k in ['status','uniform_comparisons','nonuniform_independent_comparisons',
                      'cohort_only_AUC_unweighted','cohort_only_AUC_weighted','kish_and_weights']},indent=2))


if __name__=='__main__':
    main()
