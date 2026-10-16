# Text prepared for public preregistration — v3.1

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the timestamped Persian v3.1 source (https://doi.org/10.5281/zenodo.23405371): `PREREG_v3_1_FA.md` (SHA-256 `2d964a9b6f85ea4da578d8ef244eb2b6a0338537712a8e2db8d3a4836ffd743f`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

**Status at freezing:** all preregistration-stage numerical checks have passed and the text is ready for publication; public registration has not yet occurred. Freeze date: 16 October 2026. **Author: mohammadJavad Asgari**; affiliation was not supplied and is not stated. **Reserved DOI: 10.5281/zenodo.23405371**; registration is completed only by publication of the public record.

**Proposed title:** Cross-platform transfer of a blood gene-expression model in Parkinson's disease: incremental value beyond age, sex and cell markers, with two development sources and independent PPMI testing.

## Question and hypothesis

Is PPMI AUROC higher for M5 (unpenalized demographic/cellular baseline + elastic-net-penalized genes) than M4 (the same unpenalized baseline, without the gene block)? Report the primary quantity `AUROC(M5)−AUROC(M4)` with a paired 95% CI. The research-relevance margin is 0.05; no clinical-benefit threshold is claimed. Positive, negative and inconclusive results are all accepted.

## Data and prior exposure

Development: GSE57475/PROBE, 142 records, and GSE99039/GENEPARK-TRO, 239 complete cases. Test: PPMI Project133 IR3-B38, 311 idiopathic PD and 161 baseline controls with linked counts according to the local execution report of 15 October. No real disease/sex/age model and none of their PPMI scores have been run when this package is prepared. Metadata and group labels were examined for feasibility, counts and linkage; this test has not been completely unknown from the outset. Checking the structure of 4,871 count files and gene coverage was not performance evaluation.

## Registered method

Exact population definitions, input/exclusion rules, five models, two auxiliary models, six marker panels, tuning grid, splits and CIs are in **PROTOCOL_v3_1_EN.md and config_v3_1.json attached to this registration**. Use processed GEO data, GENCODE v29 and within-sample ranks, with count-length correction for RNA-seq. All fitting/tuning is training-only; primary transfer is inductive, with no PPMI centering. M4 and M5 are both binomial/logit with unpenalized baselines; glum 3.4.1 uses P1=P2=0 for baseline and 1 for genes. Penalty/weight correctness was checked using independent synthetic tests. All M1–M5 weights in all fits are calculated from training data as n/(2*K*n[source,class]); tuning uses mean AUC within validation sources. Final Kish 334.0463 out of 381 measures weight dispersion. The lambda grid expands by one decade toward the selected boundary for at most two rounds, using development CV only. Synthetic details/results are in KKT_NUMERICAL_RESULTS_20261016_EN.md; the null rerun passed final KKT; complete details/history are in KKT_NUMERICAL_RESULTS_20261016_EN.md.

Development evaluation uses five×five nested CV within each source and two directional transfers. Final models are fitted in GEO and saved/hashed before external testing. Primary CI uses 5,000 paired, class-stratified donor-bootstrap repetitions with fixed models; DeLong is not an alternative decision. Sample size, AUROC, AUPRC and uncertainty are reported separately for each population. Sensitivities are restricted to the protocol's bounded list.

## Controls and consequences

PPMI linkage comes from identifiers and release documentation. Sex is only independent QC; concordance does not prove individual identity. The age model must first demonstrate capacity in GEO transfer. The sex gate in PPMI controls is AUROC≥0.95 and CI lower bound≥0.90; BA at 0.5 is descriptive only. Development must achieve AUROC≥0.95 in both GEO transfers. Linkage/transformation/gate failure causes suspension; continuation requires a demonstrable technical bug, a timestamped amendment and repetition of the same control. Without a demonstrable technical explanation, confirmatory testing is unevaluable. Controls and their decision consequences are registered before the first PPMI score. Primary CBC analysis is in healthy controls, with PD separately exploratory; weak/unavailable CBC does not replace models or stop the analysis, but the marker baseline is described as having insufficient/partial cell validation and claims are limited. Longitudinal similarity is reported only; there is no automatic exclusion or reassignment.

## Precision limitations and interpretation

In some prediction-correlation scenarios, n=472 is insufficient to exclude an increase of 0.05 with appropriate precision. An estimate near zero with CI half-width around 0.055 does not exclude that increase. In the M5≡M4 branch, the result is "genes not selected in development" and a tautologically zero CI is not evidence against potential gene value; no primary replacement model is introduced. For distinct models, exclusion of an increase ≥0.05 requires CI upper bound<0.05; if the CI includes both zero and 0.05, the result is inconclusive. Design simulation is not actual performance or guaranteed power. Conclusions concern this model, input and population; causal shortcut contributions or absence of biological PD signal are not inferred.

## Changes and publication

Scope is closed after v3.1. Only a documented material failure affecting confirmatory validity justifies a timestamped amendment and targeted check. An unfavorable result or a new critique alone does not authorize a new control, endpoint, margin or model choice. All deviations and post-result analyses are labeled.

This preregistration must be publicly registered in Zenodo with an identifier/DOI **before the first PPMI score, including sex, age and cell-marker scores**. Attach an exact Git commit/tag if used; a local tag/hash alone is not public registration. The attachment contains no individual PPMI data, cache, sample identifiers or individual counts. Comply with the PPMI publication policy, DPC administrative review and required acknowledgment before manuscript submission.

**Reserved identifier:** 10.5281/zenodo.23405371; not yet public at freezing. **Public registration time:** determined by the Zenodo publication timestamp. **Commit/tag referenced in registration:** not yet determined.
