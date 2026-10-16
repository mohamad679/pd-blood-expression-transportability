# KKT numerical-stage results — 16 October 2026

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the timestamped Persian v3.1 source (https://doi.org/10.5281/zenodo.23405371): `KKT_NUMERICAL_RESULTS_20261016_FA.md` (SHA-256 `2cb61a9f549e82ea011027b68584690c584a99e36f4762dacf5f37359eebab5c`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

**Final numerical status: scale, regression checks, the 1000-gene numerical stress test and the complete null passed. This document is frozen before Zenodo publication; at freezing, public registration has not yet occurred and no individual PPMI score has been run.** Only synthetic data are reported here.

## Plan and documented change

The document written and hashed before execution was `KKT_STEP_ZERO_v3_1_FA.md` (SHA-256 `4c59a9d6c95b223edfcffcaeaad91cad0e4a7139f42068d180e6628064db5cfb`); `KKT_STEP_ZERO_PREEXEC_SHA256.json` records that original filename and hash. `KKT_STEP_ZERO_v3_1_EN.md` is a later translation. All P1/P2 weights, normalized loss, training-standardized inputs and the intercept enter KKT. The max-absolute=1e-6 criterion was not relaxed after outputs. One cold retry with more iterations and gradient_tol=1e-8 was defined in advance. Consequences were distinguished for CV, final M4/M5 and sex/age controls.

The subsequent amendment selected irls-ls for pure L2 because cyclic CD was extremely slow in the correlated 1000-gene fixture, not because of real AUC or a changed desired result. The original amendment is `KKT_NUMERICAL_AMENDMENT_20261016_FA.md` (SHA-256 `ffad1329e0493ddd2f83536a77f3726901eabe910909ad23020e71dcae0ebe5f`); `KKT_NUMERICAL_AMENDMENT_20261016_EN.md` is its later translation. For rho>0 and baseline-only, irls-cd remains in use. The lambda/rho grid and acceptance threshold were unchanged.

## Scale agreement with independent optimization

Both the first run and the post-solver-amendment run contained 48 comparisons, all passing: 2 independent seeds for calibration and a separate check, 40/1000 genes, uniform/nonuniform weights, lambda=0.01/0.2 and rho=0/0.5/1. Positive nonuniform gene P1/P2 weights were also used; the baseline was unpenalized.

- Maximum glum KKT_max: **1.8293e-10**.
- Maximum independent-SciPy KKT_max: **3.3167e-8**.
- Maximum objective difference between solvers: **6.4282e-14**.
- Deliberately doubling the loss gradient was detected in fits with nonzero genes.
- Maximum probability difference was about **7.35e-7**; uniqueness of all L1 coefficients was not required.
- One line-search warning occurred in the separate check at p=1000,rho=0.5,lambda=0.01; final KKT was **1.6410e-11**, with valid independent objective agreement. An internal warning alone therefore does not prove final nonconvergence.

Full file: `kkt_scale_v3_1_solver_v2.json`; fit-level log: `kkt_scale_v3_1_solver_v2_fits.jsonl`; residual distributions by seed/dimension: `kkt_residual_distribution_v3_1.json`. First code and report versions remain in `kkt_history_first_pass_20261016`.

## Rerun of 23, 66 and 9

The same seeds and previous assertions were run in the separate `kkt_regression_rerun_20261016` directory. The 23 v3 checks passed using the original helper as the regression reference; all 66 constant-weight comparisons and 9 nonuniform-weight comparisons against SciPy passed with the new KKT-checked, solver-dispatching fitter. In total, 84 new weighted/rescaled fits have numerical logs; all passed on the first attempt, without warnings, with maximum KKT_max **1.9641e-9**. The 84 is a fit count, not a count of independent research tests.

## 1000-gene stress test and limits of the claim

The fixture uses n=381, baseline=10, p=1000 and a strong platform effect. It checks the entire grid of 15 candidates plus baseline-only under two weight settings: 32 final-training-like fits. The first version was interrupted after four completed fits because the ridge candidate was slow; this is not called mathematical nonconvergence or KKT failure. No result from an unfinished fit was used for selection.

In the amended version, **all 32 fits passed**; initial KKT failures, unresolved failures and required retries were all zero. Across the grid, maximum final KKT_max was **6.1943e-9**. **15 internal warnings in three fits** were logged; all those fits passed final KKT. This fixture therefore supplied no evidence for removing a candidate or restricting the grid. Grid and threshold remain their prior values.

The full stress run took about **1454.45 seconds (24.24 minutes)**. The slowest fit was weighted rho=0.5 and lambda=0.0001: **532.88 seconds**, 8 internal warnings, KKT_max=**3.8582e-10**. This is computational cost in a strongly correlated fixture, not numerical failure or a prediction of all real-fold runtimes. The full log is preserved in `stress_1000_kkt_v3_1_solver_v2_fits.jsonl`.

This is only one numerical fixture; it validates neither all real folds, statistical power nor third-platform transfer. The earlier null also has an independent test from the same two sources; that limitation remains explicit.

## Complete null rerun and resolution of warning ambiguity

All **200 paired repetitions** were rerun with the same seed=20261016 and previous generator. In total, **35,100 unique fits/attempts** were recorded: initial KKT failures, unresolved failures and retry attempts were all **zero**. Maximum KKT_max was **4.2659e-8**, below the fixed 1e-6 limit.

There were 975 warning messages in 227 fits and 84 repetitions: 883 coordinate-descent and 92 line-search messages. All 227 warned fits passed final KKT. Of the warned fits, 52 were unweighted inner CV, 171 weighted inner CV and 4 selected final weighted fits. No selected final unweighted fit warned. Warnings were neither removed nor hidden; each message remains with arm, rep, fold, candidate, solver and residual in the log.

| Independent synthetic output | Unweighted/pooled tuning | Weighted/macro tuning |
|---|---:|---:|
| Gene-block selection | 200/200 = 100% | 130/200 = 65% |
| Wilson CI for selection rate | 98.12–100% | 58.16–71.27% |
| Mean pooled test AUC | 0.5520096 | 0.5051516 |
| MCSE of that mean | 0.0016767 | 0.0020330 |
| Mean test macro-AUC | 0.5037540 | 0.5019759 |
| MCSE of that mean | 0.0020746 | 0.0021922 |
| Bounded grid-expansion count | 115/200 | 47/200 |
| Selected gene candidate at a remaining boundary | 0/200 | 0/200 |

The paired selection-rate reduction was **35 percentage points**, with an approximate Monte Carlo CI of **28.37 to 41.63**. No repetition changed its gene-selection decision from the first run; the maximum change in test macro-AUC was also **zero**. The prespecified null diagnostic criteria passed; these results still do not predict the paper or establish a Type-I-error rate, and they do not imply a "35% baseline-only probability in the real study."

Full report: `source_null_kkt_v3_1.json`; individual-fit log: `source_null_kkt_v3_1_fits.jsonl`; count audit: `source_null_kkt_v3_1_audit_summary.json`. The new run took about **1075.26 seconds (17.92 minutes)**. Solely to save time, the independent null computation overlapped the last two stress fits; no acceptance condition changed. The initial code/report with 1,394 warnings and the first incomplete stress run are preserved separately. The new run resolves the numerical concern necessary to continue; no retrospective claim is made about unrecorded KKT in the first run.

## Decision and remaining work

**The numerical prerequisite for registering v3.1 is met.** The grid was not restricted, KKT threshold not changed, no candidate removed, and no real data informed these decisions. The L2-specific solver change is documented with rationale, independent verification and full history.

Reserved DOI: **10.5281/zenodo.23405371**. DOI reservation and this local file are not public registration; the record becomes citable registration only upon publication. The package contains only methods, code, public gene references and synthetic results.

The real project has not yet run. After public registration, remaining tasks are the full pipeline implementation and leakage/split/bootstrap checks, GEO training/evaluation, fixed transfer controls and finally conditional confirmatory PPMI evaluation. Passing this report's checks does not establish statistical power, correctness of every future fold, biological donor identity, third-platform transfer or PD performance.
