# GEO partial-ComBat sensitivity: draft article tables (2026-10-20)

All entries are held-out AUROC. **Δ = M5 − M4**, calculated before rounding. The within-cohort values are arithmetic means over the **five original outer folds** (not an AUROC on pooled out-of-fold predictions). Transfer values are single-direction held-out evaluations. `GSE57475` is PROBE; `GSE99039` is GENEPARK. These are exploratory sensitivity estimates, not a replacement for the uncorrected primary results.

## Within-cohort nested CV

| Analysis | PROBE M4 | PROBE M5 | PROBE Δ | GENEPARK M4 | GENEPARK M5 | GENEPARK Δ |
|---|---:|---:|---:|---:|---:|---:|
| Original, no ComBat | 0.561 | 0.591 | +0.030 | 0.663 | 0.686 | +0.023 |
| Inductive, no disease covariate | 0.594 | 0.559 | −0.035 | 0.624 | 0.607 | −0.017 |
| Inductive, disease covariate in training only | 0.594 | 0.596 | +0.001 | 0.648 | 0.674 | +0.026 |
| **Transductive**, no disease covariate | 0.596 | 0.594 | −0.002 | 0.614 | 0.611 | −0.003 |

## Directional cross-cohort transfer

| Analysis | PROBE → GENEPARK M4 | M5 | Δ | GENEPARK → PROBE M4 | M5 | Δ |
|---|---:|---:|---:|---:|---:|---:|
| Original, no ComBat | 0.667 | 0.488 | −0.178 | 0.635 | 0.627 | −0.009 |
| Inductive, no disease covariate | 0.655 | 0.604 | −0.051 | 0.612 | 0.612 | +0.000 |
| Inductive, disease covariate in training only | 0.655 | 0.655 | +0.000 | 0.604 | 0.601 | −0.003 |
| **Transductive**, no disease covariate | 0.600 | 0.555 | −0.045 | 0.595 | 0.595 | +0.000 |

PROBE → GENEPARK tests 239 samples (120 PD); GENEPARK → PROBE tests 142 samples (93 PD). Transfer rows are directional: the arrow identifies the **training → test** cohort. The transductive rows use expression from the designated held-out scope to estimate correction, so they must not be presented as prospective, inductive external validation. `pooled/final` has no held-out metric and is not in either table.

## Interpretation for Results and Discussion

Among the corrected arms, only the GENEPARK within-cohort **disease-preserving** arm shows an appreciable positive M5−M4 difference (+0.026); the two arms without a disease covariate give Δ ≤ 0 in GENEPARK. This pattern should be read alongside the [pre-outcome interpretation cautions](GEO_COMBAT_AUDIT_DECISION_20261019_EN.md): with batch and diagnosis imbalanced, retaining disease as a correction covariate can distort apparent group differences. Nygaard, Rødland & Hovig ([*Biostatistics*, 2016](https://doi.org/10.1093/biostatistics/kxv027)) discuss that risk in downstream differential-expression inference. Applying that warning to classifier AUROC is a cautionary inference, **not proof** that the +0.026 is an artifact. It is also smaller than the between-fold AUROC variability in GENEPARK (original M5 fold SD ≈ 0.105).

Without a disease covariate, inductive partial ComBat lowers GENEPARK M5 mean AUROC from 0.686 to 0.607. This is compatible with some original performance reflecting batch-associated expression, but correction can also remove genuine PD-associated variation because batch and disease are confounded. The analysis cannot separate those explanations. M4 changes too, because its panel scores are recalculated from corrected expression. This M4 implementation detail was clarified **after outcomes** in the decision document and is marked there accordingly.

Neither inductive transfer direction shows a positive gene-block increment after partial ComBat. The no-disease PROBE → GENEPARK M5 rises in absolute AUROC relative to the uncorrected M5 (0.488 → 0.604), but remains below the corresponding corrected M4 (0.655). That comparison does not support added value from the M5 gene block. Equal M5 and M4 entries in some rows can arise when the tuning rule selects the no-gene branch; they should be described as “gene block not selected,” not as evidence that its biological effect is zero.

## Methods and numerical limitations to retain with the tables

- This is **partial ComBat**: batches with fewer than three training samples, or unseen held-out batches in inductive arms, remain on original ranks. The initial whole-partition exclusion rule was changed after a metadata-only feasibility count **before** real GEO ComBat outcomes. The [execution audit](GEO_COMBAT_EXECUTION_AUDIT_20261020_EN.md) gives the exact chronology, including that three first-run arm-scope AUROC outputs were later seen after the amendment and reused in the completed retry. The first-run runner hash was not saved, and the runner file was edited before retry; exact source-code identity across both starts cannot be established retrospectively.
- ComBat operates on **within-sample ranks**, not raw or normalized expression intensities. The correction and all downstream panel/gene feature construction are refit within the appropriate fold for inductive arms.
- The original full run did not retain `neuroCombat` stderr. A one-scope [numerical audit](../results/geo_combat_numeric_audit_20261020/SUMMARY.json) found 25/25 corrected matrices finite but captured 30 divide-by-zero warnings from the package's relative convergence-change calculation. Matrix finiteness does not prove solver convergence, and this check covers only `GSE99039/outer/0`. These sensitivity estimates therefore carry an unresolved numerical warning and should not be described as numerically validated across all 39 arm-scopes.
- The original results and the three ComBat arms use the same frozen GEO split-manifest SHA-256: `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. No PPMI data were accessed in this sensitivity.

Data sources: [original GEO aggregate summary](../results/geo_development_20261018/aggregate_summary.json), [completed ComBat summary](../private_data/geo_combat_sensitivity_20261020_retry1/SUMMARY.json), and [numerical diagnostic summary](../results/geo_combat_numeric_audit_20261020/SUMMARY.json).
