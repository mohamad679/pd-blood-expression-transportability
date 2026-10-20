# Supplementary Table S3. GEO transfer-capacity controls

These auxiliary controls were defined separately from the M1–M5 disease comparisons. All tests used the opposite GEO cohort as held-out data. The bootstrap intervals describe fixed-model uncertainty on the observed held-out participants; they do not include training uncertainty.

| Training → test cohort | Sex-marker AUROC [95% CI] | Sex capacity gate | Age-predictor Spearman [95% CI] | Age capacity gate |
|---|---:|:---:|---:|:---:|
| PROBE → GENEPARK | 1.000 [1.000, 1.000] | Pass | −0.058 [−0.239, +0.121] | Not met |
| GENEPARK → PROBE | 0.993 [0.975, 1.000] | Pass | +0.259 [−0.051, +0.522] | Not met |

The sex control uses five registered Y-chromosome markers and was fit only in healthy controls. In both directions it selected lambda=100 at the upper grid edge after two permitted expansions; those two auxiliary choices are not included in Supplementary Table S1's disease-model boundary count. The age predictor is a separate autosomal-gene control and did not meet its registered transfer-capacity criterion.

Editorial provenance: `docs/GEO_TRANSFER_CONTROLS_20261018_EN.md`. No PPMI scores were generated.
