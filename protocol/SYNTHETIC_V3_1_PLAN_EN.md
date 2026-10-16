# Synthetic test plan v3.1 — before observing outputs

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the timestamped Persian v3.1 source (https://doi.org/10.5281/zenodo.23405371): `SYNTHETIC_V3_1_PLAN_FA.md` (SHA-256 `45d15dc4576b92301016250bdf5eb369eb2a3b0a56499e4ef1eb24e6544fab1a`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

Drafting date: 15 October 2026. This plan uses only synthetic data and aggregate source counts; no individual GEO/PPMI record, expression or score is read. This file is not public registration.

## Numerical correctness

1. Compare the v3 test-model settings with constant weights 0.1, 1 and 10 against the unweighted fit. Coefficient-difference tolerance is 2e−5; for nonuniform weights, multiplying all weights by 7 must also preserve this invariance.

2. With nonuniform weights, optimize in SciPy the binomial objective with loss divided by the sum of weights and the same gene-block penalty. L1 uses positive/negative variables and bounds. Tolerances are coefficient difference 3e−5, objective difference 1e−9 and unpenalized-baseline gradient 2e−6. Optimizer or convergence failure is recorded, not ignored.

3. For each training set with K sources, `w_i=n/(2*K*n[source_i,class_i])`; each source receives equal mass and each class within a source receives equal mass. The same formula balances the two classes for a single source. Check weight sum, source×class cell mass and Kish=`sum(w)^2/sum(w²)`. Source-identifier AUC using the real counts must be 1951/3408 without weights and 0.5 with these weights.

## Prespecified null simulation

- 200 paired repetitions, master seed 20261016. Each repetition trains on two sources with 142 and 239 people and PD/control counts 93/49 and 120/119. Generate an independent artificial test with the same source×class counts.
- 40 artificial gene features, with a strong platform shift across the two gene halves, followed by within-sample ranking; 10 artificial baseline features with the same age/sex template and six random covariates without a source effect. Class is independent of features within each source. Thus **within-source disease signal is zero**, but the marginal source-driven association is deliberately present. This generator is only a mechanistic stress test; it is not equivalent to all biological features or the real study's 1,000 genes.
- Two arms on the same data and the same 5 source×class inner folds: (a) unweighted fitting with pooled-AUC tuning; (b) source×class weights with mean within-source-AUC tuning. Baseline fits also use their arm's weights. M4 is the no-gene candidate and is preferred for ties within 1e−6.
- Initial lambda grid `{1e−4,1e−3,1e−2,1e−1,1}`, rho `{0,0.5,1}`. If the selected gene-containing candidate is at a boundary, add a decade in that direction only within the same training CV, for at most two rounds. Baseline-only or an entirely zero block is not a gene-boundary selection. Report repetitions, expansion counts and remaining boundaries; the grid is not expanded until a desired result is obtained.
- Outputs: nonzero-block selection rate with Wilson CI, paired rate difference and Monte Carlo CI, pooled and macro AUC in the independent test, boundary frequency and effective training sample size. Gene selection under the null is not a Type-I-error test; similar rates alone also do not establish that the correction is necessary or useless.

## Consequences of results

Numerical correctness under the tolerances above is go/no-go. For the null, a higher unweighted selection rate is not guaranteed in advance. If the difference CI does not show a reduction, or mean independent-test macro-AUC differs from 0.5 by more than `3*MCSE+0.01`, flag outputs for review and **do not proceed with public registration before examining the numbers**. A flag is not proof of a coding error or an invalid idea. Any convergence failure is reviewed before publication as well. Retain all outputs, not just the preferred arm.
