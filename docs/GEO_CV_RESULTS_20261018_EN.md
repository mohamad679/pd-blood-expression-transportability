# GEO development CV and transfer results

**Run date:** 17–18 October 2026 (Tehran time)  
**Protocol:** registered v3.1 / public amendment 3.1.1  
**Analysis status:** all 13 prespecified GEO scopes completed; external PPMI scoring not started.

## Scope and provenance

The four local GEO inputs matched the frozen input hashes. The analysis used 381 records across GSE57475 (142; 93 PD, 49 controls) and GSE99039 / GENEPARK-TRO (239; 120 PD, 119 controls), with 16,364 common genes. The pipeline code SHA-256 was `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`; the split-manifest SHA-256 was `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. Both matched the public implementation manifest in [Zenodo amendment 3.1.1](https://doi.org/10.5281/zenodo.23405802) before these fits.

The 13 scopes were five outer folds within each GEO source, two directional transfers, and one pooled final fit. All tuning and preprocessing followed the registered source/class weighting, within-validation-source macro-AUROC criterion, and training-only operations. The pooled final fit did not evaluate a test set. No PPMI data or outcome scores entered any run.

## Nested outer-fold performance

Values below are the unweighted mean and standard deviation of the five prespecified outer-fold metrics within each source. They are descriptive estimates from the registered development CV, not clinical performance estimates. Fold-level test sets were small and their estimates varied appreciably.

| Source | Model | AUROC, mean ± SD | AUPRC, mean ± SD |
|---|---|---:|---:|
| GSE57475 / PROBE | M1 | 0.475 ± 0.074 | 0.640 ± 0.048 |
|  | M2 | 0.556 ± 0.074 | 0.706 ± 0.085 |
|  | M3 | 0.537 ± 0.141 | 0.694 ± 0.096 |
|  | M4 | 0.561 ± 0.170 | 0.761 ± 0.086 |
|  | M5 | 0.591 ± 0.123 | 0.762 ± 0.075 |
| GSE99039 / GENEPARK-TRO | M1 | 0.584 ± 0.075 | 0.612 ± 0.096 |
|  | M2 | 0.704 ± 0.093 | 0.722 ± 0.072 |
|  | M3 | 0.670 ± 0.120 | 0.687 ± 0.083 |
|  | M4 | 0.663 ± 0.127 | 0.690 ± 0.098 |
|  | M5 | 0.686 ± 0.105 | 0.707 ± 0.100 |

The within-source mean M5−M4 AUROC differences were +0.030 in GSE57475 and +0.023 in GSE99039. These fold summaries are descriptive and do not replace the preregistered external comparison or its uncertainty analysis.

## Directional cohort transfers

| Training source → held-out source | Held-out n (cases) | M4 AUROC | M5 AUROC | M5−M4 AUROC | M4 AUPRC | M5 AUPRC |
|---|---:|---:|---:|---:|---:|---:|
| GSE57475 → GSE99039 | 239 (120) | 0.667 | 0.488 | −0.178 | 0.628 | 0.485 |
| GSE99039 → GSE57475 | 142 (93) | 0.635 | 0.627 | −0.009 | 0.754 | 0.761 |

The gene-augmented M5 did not improve on M4 in either directional transfer; the first direction showed a large deterioration. This is a development-transfer finding for these two microarray cohorts and this locked model procedure. It does not estimate the causal fraction of disease versus shortcut signal, establish that blood contains no PD biology, or predict PPMI performance.

## Final pooled development fit

The registered pooled inner-CV tuning selected a nonzero gene block for M5 (ridge-only, `lambda=0.001`, `rho=1.0`), so the final M5 is distinct from M4. M2 selected `lambda=0.1`, `rho=0.0`; M3 selected `lambda=0.01`, `rho=1.0`. These are training decisions only. The corresponding model and preprocessing artifacts are held locally and were not added to this public repository; they must be hashed in the second pre-PPMI lock before any PPMI score is generated.

## Numerical and execution audit

- All 12 post-pilot run directories completed, including the 10 remaining outer/transfer scopes and pooled final fit. Together with the earlier GSE57475 outer/0 pilot, this covers all 13 registered scopes.
- The 3,146 fit records from the 12 new scopes passed the locked `KKT_max ≤ 1e-6` rule. There were zero failed records, retries, or exceptions; the largest recorded residual was `4.1093e-8`. The pilot separately had 269/269 passing fit attempts and a maximum residual of `3.8924e-8`.
- No nonempty fit-level warning field was recorded. Six `glum` sigmoid-clipping warnings were observed on stderr during validation-prediction/tuning calls. They were not attached to candidate/fold records by this implementation, so their exact attribution is unavailable. They indicate numerical saturation/clipping, not a recorded solver/KKT failure; the warning and logging limitation are retained here and require review before any external scoring. No pipeline change was made after these locked GEO results.
- The new scopes summed to about 1.67 wall hours and 6.82 process-CPU hours; including the pilot, all scopes used about 1.79 summed wall hours and 7.34 process-CPU hours. Process CPU sums across the backend's threads and is not elapsed time.
- Run directories contained sample-level predictions and fit artifacts. They remain in the local temporary workspace and are intentionally excluded from Git. This report contains aggregate metrics only.

## Remaining work before any PPMI score

1. Review the sigmoid-clipping warning pathway and make a timestamped, public amendment before any affected re-run if a methodological change is warranted. Do not silently alter this completed analysis.
2. Complete the registered GEO sex/age transfer controls, CBC marker validation, bounded sensitivities, leakage checks, and bootstrap-coverage simulation; finalize the permutation scope before execution.
3. Reconcile and lock the final PPMI population and one-to-one file/metadata crosswalk. The preliminary 311 PD / 161 control count is not yet a final eligibility count.
4. Hash the locked code/environment, splits, fitted M1–M5 models and preprocessing, and a reference to the private population/crosswalk. Publish the second pre-PPMI Zenodo version before computing any PPMI score.

The project has moved beyond feasibility and completed its primary GEO CV/transfer milestone. It has **not** completed external validation, and the conditions for the confirmatory PPMI analysis are not yet satisfied.
