# Project status

**Snapshot date:** 2026-10-18 (Asia/Tehran)

| Workstream | Status | Evidence / limitation |
|---|---|---|
| Feasibility and cohort choice | Completed for the registered scope | Two GEO development cohorts and an approved PPMI access path were identified. Feasibility did not itself constitute external validation. |
| GEO input and mapping audit | Completed | 381 GEO records, 16,364 shared genes; all four prepared GEO input hashes match the frozen audit. |
| Protocol and implementation lock | Publicly versioned | Pipeline and split hashes matched Zenodo amendment 3.1.1 before GEO model execution. |
| GEO nested CV and directional transfer | Completed for all 13 registered scopes | Five outer folds per GEO source, both directional transfers, and pooled final tuning/fitting. Aggregate results and numerical audit: [GEO_CV_RESULTS_20261018_EN.md](GEO_CV_RESULTS_20261018_EN.md). |
| Numerical fit audit | Passed with a warning-attribution limitation | All 3,146 post-pilot fit records passed KKT ≤ 1e-6; zero retries/exceptions. Six sigmoid-clipping warnings appeared during prediction/tuning but were not attributed to individual candidates by the pipeline; review and any pre-registered resolution are required before PPMI scoring. |
| Registered GEO controls and diagnostics | Incomplete | Sex/age transfer controls, CBC validation, bounded sensitivities, leakage audit, bootstrap coverage, and permutation scope/execution remain. |
| PPMI linkage and final population | Incomplete | The provisional 311 PD / 161 control count is not final. The identity/crosswalk and eligibility set must be reconciled and locked. |
| PPMI external scoring | Not started | No PPMI disease scores have been calculated. A second public lock of code/environment, final models/preprocessing, and population reference must precede any PPMI score. |
| Manuscript | Not started | GEO results are now available; external validation and remaining diagnostics are not complete. Causal shortcut fractions and clinical utility are not inferred. |

## GEO milestone

All planned nested outer folds, directional transfers, and pooled final fitting completed on 17–18 October 2026. Across the five outer folds, mean AUROC was 0.591 for M5 versus 0.562 for M4 in GSE57475 and 0.686 versus 0.663 in GSE99039. In directional transfer, M5 was below M4 by 0.178 AUROC from GSE57475 to GSE99039 and by 0.009 in the reverse direction. These are development results for the registered procedure and two GEO sources; they are not PPMI performance or causal estimates.

The four input hashes match the frozen audit, and the split-manifest SHA-256 is `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. The public implementation manifest covers code SHA-256 `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`. No individual-level predictions or PPMI data are in Git.

## Next gates

1. Review the sigmoid-clipping warning pathway and decide whether a timestamped public amendment is needed before any affected re-run.
2. Complete the registered GEO controls, sensitivities, leakage checks, bootstrap-coverage validation, and permutation plan.
3. Finalize PPMI file-to-release linkage and eligibility using the guide and complete sample-name segments; document publication conditions.
4. Before PPMI scoring, publish the second lock of code/environment, splits, fitted M1–M5 models and preprocessing, and a non-identifying reference/hash for the private population/crosswalk.
