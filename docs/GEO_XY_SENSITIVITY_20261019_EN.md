# GEO X/Y-removal sensitivity — 19 October 2026

## Execution

The preregistered X/Y-removal sensitivity was run on all 13 frozen GEO scopes:
five outer folds per dataset, two directional transfers, and the pooled final
fit. The 15,720 autosomal genes were eligible for the added gene block, and the
same training-only source-stratified MAD rule selected up to 1,000 genes. The
published primary pipeline and split hashes matched exactly. No PPMI data were
loaded or scored.

The run completed 3,429 checked fits in 5,904.75 wall seconds (1 h 38 min 25 s)
and 25,063.10 process-CPU seconds. All fits passed the numerical gate; maximum
recorded KKT residual was 4.11e-8. One `ConvergenceWarning` was recorded for
GSE99039 transfer, inner fold 3, M3 (`lambda=0.01`, `rho=1.0`): line search
failed and the next iteration was reported as very close to the current one.
Several sigmoid-clipping warnings appeared in the terminal. They were not
attached to individual candidate/fold records by this wrapper, so exact
attribution is unavailable. They indicate numerical saturation/clipping, not
a recorded KKT failure; preserve this logging limitation in any manuscript
description and do not call the run warning-free.

Leakage and invariance checks passed in every scope: all fit training hashes
matched a frozen training partition, all selected genes were autosomal, and
the M1/M4 baseline transforms were bit-identical before and after the
sensitivity in training and held-out data. Held-out data were scored only
after inner-fold selection. The pooled final fit has no held-out metric by
design.

## Descriptive held-out AUROC comparison

Values below are unweighted means across the five frozen outer folds within
each GEO dataset. `Delta` is X/Y-removal minus the original published result.
These are descriptive cohort-specific summaries; no pooled estimate or new
uncertainty interval was computed here.

| Dataset | Model | Original mean AUROC | X/Y removed | Delta |
|---|---:|---:|---:|---:|
| GSE57475 | M2 | 0.5564 | 0.5349 | -0.0215 |
| GSE57475 | M3 | 0.5371 | 0.4943 | -0.0428 |
| GSE57475 | M5 | 0.5911 | 0.5332 | -0.0580 |
| GSE99039 | M2 | 0.7042 | 0.7062 | +0.0019 |
| GSE99039 | M3 | 0.6704 | 0.7150 | +0.0445 |
| GSE99039 | M5 | 0.6861 | 0.6920 | +0.0059 |

M1 and M4 were exactly unchanged, as expected, because their baseline block was
not altered. Directional-transfer AUROC deltas (X/Y removed minus original)
were +0.0203/+0.0116/+0.0308 for M2/M3/M5 when training on GSE57475, and
-0.0215/-0.0342/-0.0101 when training on GSE99039. The changes were smaller
than between-fold variability and did not form a systematic performance
change. The safest summary is: **X/Y removal did not produce a systematic
change in performance; the differences were smaller than between-fold
variability.** After X/Y removal, M5−M4 was -0.028 in PROBE and +0.029 in
GENEPARK; these small, opposite-signed contrasts do not establish reliable
incremental value from the added gene block.

## Interpretation and artifact

This completes the registered GEO X/Y-removal sensitivity. It does not resolve
the broader shortcut question: effects vary by dataset, model, and transfer
direction. Interpret alongside the original results and the other registered
sensitivities, not as a causal attribution. The run output remains local under
`private_data/geo_xy_sensitivity_full_20261019/`; sample-level predictions and
fit logs are not for public release. The aggregate result is
`private_data/geo_xy_sensitivity_full_20261019/RESULTS.json`.

Next: run the registered training-only age/sex residualization sensitivity on
the same frozen GEO scopes, retaining the same no-leakage checks and held-out
evaluation order. Benchmarking for X/Y removal is complete; do not rerun it
without a specific reproducibility defect.
