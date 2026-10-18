# Completed GEO label-permutation control — 19 October 2026

The registered, exploratory B=200 control completed on the Hetzner CPX62
with eight workers. The pilot's registered replicate 0 was reused. The frozen
run summary reports `status: complete`, `B: 200`, `PPMI_accessed: false`, and
`exploratory_only: true`. The returned archive contains `COMPLETE.json` for
each fixed replicate ID 0–199 and has a readable gzip/tar structure. There
were no `FAILURE.json` files in the returned archive. These are structural
and summary checks; this report does not claim an independent rerun of every
fit or a complete second hash verification of all 17,498 archive entries.

| GEO cohort | Statistic | Observed mean | Null mean | Null SD | Exploratory p | Approx. Monte Carlo SE |
|---|---|---:|---:|---:|---:|---:|
| GSE57475 / PROBE | M4 AUROC | 0.561462 | 0.506547 | 0.068529 | 0.203980 | 0.028493 |
| GSE57475 / PROBE | M5 AUROC | 0.591111 | 0.504120 | 0.068100 | 0.104478 | 0.021629 |
| GSE57475 / PROBE | M5 − M4 AUROC | 0.029649 | -0.002427 | 0.042895 | 0.189055 | 0.027687 |
| GSE99039 / GENEPARK-TR | M4 AUROC | 0.662953 | 0.508113 | 0.052584 | 0.004975 | 0.004975 |
| GSE99039 / GENEPARK-TR | M5 AUROC | 0.686081 | 0.507257 | 0.049672 | 0.004975 | 0.004975 |
| GSE99039 / GENEPARK-TR | M5 − M4 AUROC | 0.023128 | -0.000856 | 0.044142 | 0.233831 | 0.029929 |

The cohort statistic is the unweighted mean of five within-source outer-fold
AUROCs. The published amendment fixed 200 coherent donor-label-map
replicates and the one-sided tail calculation
`(1 + count[T_b >= T_observed - 1e-12]) / 201`. The lowest attainable p-value
is 1/201 = 0.004975. These six comparisons are secondary and exploratory, not
the confirmatory external endpoint. Global label shuffling tests a global
association null. Its M5-minus-M4 null distribution does not test incremental
gene value conditional on the M4 baseline, nor does a non-small p-value prove
absence of information or leakage. No PPMI model score has been computed.

The server-side archive `/root/pdgeo-results-20261019.tar.gz` was copied to
the Mac, tested as a readable gzip/tar archive, and matched its server SHA-256:

```text
50d9b536f5f0f46a176605859134e23f143bd3994840d885074dc79ca5a7fead
```

The verified private local copy is
`private_data/server_return/pdgeo-results-20261019.tar.gz` (ignored by Git).
The user also retained a copy in `~/Downloads/`. The archive includes the
complete private GEO permutation tree, `SUMMARY.json`, pilot/preflight
evidence and execution logs. Individual predictions and fit artifacts remain
private; only the aggregate statistics above belong in this public-ready
documentation. Server deletion was confirmed in Hetzner Console, followed by
an empty Primary IP inventory. This stops ongoing charges for those two
resources; the final invoice has not been inspected. The unused firewall may
remain as a configuration object and is not a computation resource.

Next registered work is the bounded remaining diagnostics and sensitivities,
final PPMI identity/eligibility and model preparation, then the second public
pre-PPMI lock. This permutation result alone does not complete external
validation or the study.
