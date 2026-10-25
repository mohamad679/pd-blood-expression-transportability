# GENEPARK within-batch AUROC diagnostic

Uses saved M4/M5 predictions from the five frozen outer test folds. No refitting or PPMI access.
Each batch AUROC compares PD-control pairs within that batch over pooled out-of-fold samples. Batches with fewer than five samples in either class are listed but excluded from the unweighted macro mean.

Eligible batches: **4/11**.
Macro mean within-batch AUROC: **M4 0.655; M5 0.661**.

| Batch | Control n | PD n | Evaluable | M4 AUROC | M5 AUROC |
|---|---:|---:|:---:|---:|---:|
| 03/05/10_TRO | 4 | 8 | No | — | — |
| 03/10/10_TRO | 6 | 2 | No | — | — |
| 08/20/09_TRO | 2 | 5 | No | — | — |
| 08/25/09_TRO | 9 | 21 | Yes | 0.608 | 0.640 |
| 08/26/09_TRO | 27 | 9 | Yes | 0.601 | 0.531 |
| 08/27/09_TRO | 1 | 0 | No | — | — |
| 09/01/09_TRO | 2 | 31 | No | — | — |
| 09/02/09_TRO | 27 | 16 | Yes | 0.799 | 0.850 |
| 09/03/09_TRO | 38 | 22 | Yes | 0.612 | 0.623 |
| 09/04/09_TRO | 1 | 3 | No | — | — |
| 09/16/09_TRO | 2 | 3 | No | — | — |

Interpretation: performance above 0.5 within a batch weakens a purely between-batch explanation, but does not establish a causal or biological disease signal. Batch-specific estimates are descriptive and should be read with their sample sizes.

Frozen split SHA-256 verified: `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.
