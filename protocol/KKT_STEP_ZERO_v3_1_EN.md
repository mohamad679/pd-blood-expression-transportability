# Numerical step zero v3.1 — before the new checks

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the timestamped Persian v3.1 source (https://doi.org/10.5281/zenodo.23405371): `KKT_STEP_ZERO_v3_1_FA.md` (SHA-256 `4c59a9d6c95b223edfcffcaeaad91cad0e4a7139f42068d180e6628064db5cfb`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

Date: 15 October 2026. This is a local pre-execution plan, not public registration. The previous simulation results and all their warnings are retained without rewriting. This stage uses entirely synthetic data. Real metadata were previously viewed for feasibility, but study models and PPMI scores have not been run.

## Objective, scale and computational coordinates

`F(a,b) = sum(w_i*(log(1+exp(eta_i))-y_i*eta_i))/sum(w) + lambda*rho*sum(P1_j*abs(b_j)) + lambda*(1-rho)*sum(P2_j*b_j^2)/2`, with `eta=a+X*b`. X contains features standardized using training data only; the backend uses `scale_predictors=False`. P2 in this implementation is only the diagonal vector, not a general matrix. Baseline and intercept are unpenalized. Baseline P1/P2 are zero and gene weights in the study model are one; the scale check additionally examines positive nonuniform penalty weights. Normalized weighted NLL is used, not full deviance containing a factor of two. No test count or prevalence enters the loss.

`g0=sum(w*(sigmoid(eta)-y))/sum(w)` and `g=X.T@(w*(sigmoid(eta)-y))/sum(w)`.

KKT residual: g0 for the intercept; for each nonzero coefficient, `g_j+lambda*(1-rho)*P2_j*b_j+lambda*rho*P1_j*sign(b_j)`; for an exactly zero coefficient, the signed value `sign(g_j)*max(abs(g_j)-lambda*rho*P1_j,0)`. With P1=P2=0, the baseline reduces to g_j. Small coefficients are not set to zero merely to improve KKT. The outer criterion is `max(abs(residual)) <= 1e-6` and finite outputs; the residual L1 norm is also reported. This is distinct from glum's internal sum-norm stopping criterion and makes no claim of an AUROC-error bound.

## Solver settings and bounded retry

The initial solve remains glum 3.4.1, binomial/logit, irls-cd, max_iter=1000, max_inner_iter=100000, gradient_tol=1e-7, step_size_tol=None, scale_predictors=False, selection=cyclic, seed=20261014. All warnings are retained. If final KKT passes, an internal warning alone does not invalidate the fit.

If outputs are nonfinite, a numerical exception occurs or final KKT fails, only **one fixed retry from zero** is allowed: max_iter=5000, max_inner_iter=1000000 and gradient_tol=1e-8; all other model, data, weights, penalties and criteria remain fixed. Retry only addresses a numerical issue; lambda, rho, feature set or score is not changed to select it. Both attempts are recorded. Remaining failure invalidates the fit; a candidate is not discarded on the basis of an unreliable AUC.

Each attempt logs context, stage, arm, rep, fold, lambda/rho candidate, feature count, weights, attempt, solver settings, time, warning type/text, KKT_max, KKT_L1, status and final iteration count. A warning may concern one internal step; the warning count is not the number of nonconverged models.

## Consequences of failure

- Within-GEO CV or directional transfer: the tuning/fit procedure and the entire dependent analysis are invalid; a fold is not removed and the remaining-fold average is not that analysis's result. Other independent analyses may continue.
- Final M4 has no hyperparameter; both a valid final M4 fit and valid final M5 tuning/fit are required. Unresolved failure of either makes the confirmatory comparison unevaluable; no CI is constructed from an invalid model.
- Sex control: invalid tuning/fitting makes the control unevaluable and leaves the confirmatory gate unestablished. An invalid age control only limits the age-transfer claim; under the protocol, it is not a separate stopping gate.
- No threshold is relaxed after outputs to obtain a pass. Scale/code problems require a documented correction, version preservation and full rerun. Without a demonstrable technical explanation, method rescue on the same holdout is not allowed.

## Verification and calibration plan

1. Synthetic data with n=381, baseline=10 and p=40 or p=1000; calibration seed 20261017 and separate-check seed 20261018. Examine uniform/nonuniform weights, rho=0,0.5,1 and lambda=0.01 and 0.2. Calibration uses a stricter glum solve with gradient_tol=1e-9 and the retry iteration limits. Independent SciPy uses the objective above, with an L1 split-variable/bounds approach or a smooth solve for rho=0.
2. Scale agreement requires both solvers finite, both KKT_max <=1e-6 and objective difference <=1e-8; equal coefficients are not required for nonunique L1 solutions. Maximum probability difference is also reported. A deliberate factor-of-two loss-gradient error in a fit with nonzero genes must show disagreement; nonuniform P1/P2 weights are also examined. The independent solver is not assessed only from its success flag; independent KKT is computed too.
3. The 1e-6 criterion is proposed in advance. Residual distributions across all reference and separate-check fits are reported; this run confirms or rejects it rather than deriving a new threshold from the observed maximum residual. Rejection requires a numerically justified amendment before locking; registration does not proceed for now.
4. After passing scale and separate checks, run a real-dimension stress test with p=1000,n=381, seed=20261019, two weight settings, the entire lambda={1e-4,1e-3,1e-2,0.1,1} and rho={0,0.5,1} grid, a strong platform effect and within-sample ranking; baseline-only is also checked. This is a numerical check of final-training-like fits, not a power assessment, third-platform transfer validation or reconstruction of every inner fold. Chronic failure is reported with its candidate and KKT; the grid is not automatically restricted in this preliminary run.
5. If the 1000-gene output is unexpected or unresolved failure occurs, present numbers before continuing and keep registration on hold. After a documented resolution, rerun all previous 23/66/9 tests and the null with the same seeds, retaining both versions. Grid restriction is allowed only with a numerical rationale before locking and a full rerun of the program.

The previous null simulation has an independent test from the same two sources, not a third platform. This limitation is explicit; the 1000-gene numerical experiment also does not test new-source transfer. No individual PPMI score is allowed at this stage, and no Zenodo registration takes place.
