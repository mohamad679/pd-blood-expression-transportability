# GENEPARK batch-only diagnostic

Computed from sample metadata only on the five frozen GSE99039 outer folds. No expression or PPMI files were loaded.

Mean fold AUROC, batch-only: **0.691**.

| Fold | Batch-only AUROC | M2 AUROC | M5 AUROC | Unseen-batch fallback samples |
|---|---:|---:|---:|---:|
| 0 | 0.749 | 0.552 | 0.576 | 0 |
| 1 | 0.619 | 0.698 | 0.722 | 0 |
| 2 | 0.734 | 0.747 | 0.714 | 0 |
| 3 | 0.709 | 0.724 | 0.589 | 0 |
| 4 | 0.645 | 0.801 | 0.830 | 1 |

| Same-fold mean AUROC | Value |
|---|---:|
| M1 | 0.584 |
| M2 | 0.704 |
| M3 | 0.670 |
| M4 | 0.663 |
| M5 | 0.686 |

Interpretation: batch prevalence is predictively informative in this fixed-fold diagnostic and its mean is close to the original M5 mean. It does not establish that batch causally explains model performance or represents disease biology. One test sample had an unseen training batch and received the outer-training prevalence fallback.

Frozen split SHA-256 verified: `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.
