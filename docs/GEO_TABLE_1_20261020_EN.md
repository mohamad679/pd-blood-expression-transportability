# Table 1. GEO development performance of all five fixed model families

**A. Discrimination.** Within-cohort entries are the unweighted mean ± sample SD of the five prespecified outer-fold AUROCs. Transfer entries are AUROC with descriptive 95% percentile intervals from 5,000 paired, diagnosis-stratified donor bootstrap resamples of saved held-out predictions. The intervals condition on the fitted models; they do not include training or model-selection uncertainty.

| Model | PROBE within-cohort (n=142) | GENEPARK within-cohort (n=239) | PROBE → GENEPARK (n=239) | GENEPARK → PROBE (n=142) |
|---|---:|---:|---:|---:|
| M1: age/sex | 0.475 ± 0.074 | 0.584 ± 0.075 | 0.622 [0.549, 0.692] | 0.585 [0.483, 0.682] |
| M2: gene block | 0.556 ± 0.074 | 0.704 ± 0.093 | 0.460 [0.388, 0.534] | 0.663 [0.570, 0.754] |
| M3: age/sex + gene block | 0.537 ± 0.141 | 0.670 ± 0.120 | 0.524 [0.450, 0.597] | 0.667 [0.575, 0.757] |
| M4: age/sex + cell markers | 0.561 ± 0.170 | 0.663 ± 0.127 | 0.667 [0.597, 0.735] | 0.635 [0.537, 0.730] |
| M5: M4 + gene block | 0.591 ± 0.123 | 0.686 ± 0.105 | 0.488 [0.414, 0.563] | 0.627 [0.532, 0.718] |
| **M5−M4** | **+0.030; paired-fold SD 0.049** | **+0.023; paired-fold SD 0.074** | **−0.178 [−0.236, −0.121]** | **−0.009 [−0.089, +0.071]** |

PROBE contains 93 PD and 49 controls; GENEPARK contains 120 PD and 119 controls. The within-cohort M5−M4 entries are means of the five paired fold differences, with their *sample SD*, not differences of independently estimated SDs. These GEO-to-GEO bootstrap intervals are descriptive and distinct from the reserved PPMI confirmatory analysis.

**B. Transfer log loss.** Smaller is better. Values use the original held-out transfer predictions, without bootstrap intervals.

| Model | PROBE → GENEPARK | GENEPARK → PROBE |
|---|---:|---:|
| M1 | 0.683 | 0.668 |
| M2 | 1.286 | 0.966 |
| M3 | 1.785 | 0.956 |
| M4 | 0.664 | 0.822 |
| M5 | **6.414** | 1.297 |

Data provenance for editorial review: frozen GEO aggregate results `results/geo_development_20261018/aggregate_summary.json`; original transfer `RESULTS.json` files in `private_data/geo_runs_20261018`; descriptive transfer intervals `results/geo_transfer_descriptive_bootstrap_20261020.json` under the pre-computation decision `docs/GEO_TRANSFER_BOOTSTRAP_DECISION_20261020_EN.md`.
