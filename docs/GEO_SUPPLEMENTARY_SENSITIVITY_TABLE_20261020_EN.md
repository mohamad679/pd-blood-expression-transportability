# Supplementary Table S2. GEO technical diagnostics and sensitivity analyses

All entries are descriptive held-out AUROC. A within-cohort value is the unweighted mean of the five original outer-fold AUROCs unless stated otherwise. Transfer arrows identify training → test cohort. Δ = M5−M4, calculated before rounding. PROBE is GSE57475; GENEPARK is GSE99039.

## A. Gene-block sensitivity results

M4 is unchanged in the X/Y, residualization, and panel-gene-exclusion arms. The panel-gene arm refits M5 only.

| Analysis | PROBE M5 | GENEPARK M5 | PROBE→GENEPARK M5 | GENEPARK→PROBE M5 |
|---|---:|---:|---:|---:|
| Original | 0.591 | 0.686 | 0.488 | 0.627 |
| X/Y genes excluded | 0.533 | 0.692 | 0.519 | 0.616 |
| Train-only age/sex residualization | 0.538 | 0.699 | 0.547 | 0.600 |
| 172 panel genes excluded | 0.563 | 0.684 | 0.539 | 0.623 |

Unchanged M4 AUROC in these arms: PROBE 0.561, GENEPARK 0.663, PROBE→GENEPARK 0.667, GENEPARK→PROBE 0.635.

## B. Age/sex matching

These are AUROCs in the *matched*, label-defined population. They are not directly comparable with full-population AUROCs. The GENEPARK folds are shown individually because outer/0 and all PROBE outer folds were non-evaluable under the prespecified minimum-size rule; **the mean of four GENEPARK folds is not a valid five-fold result**. The pooled/final matching fit has no held-out evaluation and is omitted.

| Scope | Matched test pairs | M4 | M5 | Δ |
|---|---:|---:|---:|---:|
| GENEPARK outer/1 | 14 | 0.735 | 0.526 | −0.209 |
| GENEPARK outer/2 | 14 | 0.510 | 0.633 | +0.122 |
| GENEPARK outer/3 | 15 | 0.556 | 0.538 | −0.018 |
| GENEPARK outer/4 | 12 | 0.667 | 0.736 | +0.069 |
| GENEPARK→PROBE transfer | 46 | 0.623 | 0.632 | +0.009 |

## C. Partial-ComBat sensitivity

The disease covariate in the middle arm was used in *training* only. The transductive arm used held-out expression to estimate correction and is not an inductive external-validation estimate. M4 and M5 were both reconstructed from corrected expression. Values with Δ = 0 correspond to no positive gene-block increment; the tuning rule can select the no-gene branch.

| Analysis | PROBE M4/M5 (Δ) | GENEPARK M4/M5 (Δ) | PROBE→GENEPARK M4/M5 (Δ) | GENEPARK→PROBE M4/M5 (Δ) |
|---|---:|---:|---:|---:|
| Original | 0.561/0.591 (+0.030) | 0.663/0.686 (+0.023) | 0.667/0.488 (−0.178) | 0.635/0.627 (−0.009) |
| Inductive, no disease covariate | 0.594/0.559 (−0.035) | 0.624/0.607 (−0.017) | 0.655/0.604 (−0.051) | 0.612/0.612 (0.000) |
| Inductive, disease covariate in training | 0.594/0.596 (+0.001) | 0.648/0.674 (+0.026) | 0.655/0.655 (0.000) | 0.604/0.601 (−0.003) |
| **Transductive**, no disease covariate | 0.596/0.594 (−0.002) | 0.614/0.611 (−0.003) | 0.600/0.555 (−0.045) | 0.595/0.595 (0.000) |

## D. GENEPARK technical-batch diagnostics

Recorded batch metadata covered 239/239 participants in 11 batches. Batch-only prediction used the PD fraction in the *outer-training* participants of each test sample's batch, giving mean outer-fold AUROC 0.691. Four batches had at least five PD and five controls; within-batch AUROC below uses saved, pooled out-of-fold predictions from five different fitted models, with an unweighted mean over the four evaluable batches.

| Batch | Controls | PD | M4 within-batch AUROC | M5 within-batch AUROC |
|---|---:|---:|---:|---:|
| 08/25/09_TRO | 9 | 21 | 0.608 | 0.640 |
| 08/26/09_TRO | 27 | 9 | 0.601 | 0.531 |
| 09/02/09_TRO | 27 | 16 | 0.799 | 0.850 |
| 09/03/09_TRO | 38 | 22 | 0.612 | 0.623 |
| **Unweighted mean, four batches** | — | — | **0.655** | **0.661** |

Seven other batches were below the minimum size in one or both classes. GENEPARK batch–diagnosis imbalance included 09/01/09_TRO (31 PD, 2 controls) and 08/26/09_TRO (9 PD, 27 controls). PROBE chip metadata covered 142/142 participants across 16 title-derived proxies.

Editorial provenance: completed GEO sensitivity `RESULTS.json` files under `private_data/geo_xy_sensitivity_full_20261019`, `private_data/geo_age_sex_resid_full_20261019`, `private_data/geo_panel_gene_exclusion_full_20261019`, and `private_data/geo_age_sex_matching_full_20261019`; ComBat `SUMMARY.json` under `private_data/geo_combat_sensitivity_20261020_retry1`; technical reports under `results/geo_technical_metadata_audit_20261019`, `results/geo_batch_only_diagnostic_20261019`, and `results/geo_within_batch_auc_20261020`. ComBat execution and numerical limitations remain in `docs/GEO_COMBAT_EXECUTION_AUDIT_20261020_EN.md`.
