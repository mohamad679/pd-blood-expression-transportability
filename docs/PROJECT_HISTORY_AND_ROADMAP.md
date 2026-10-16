# Project history and execution roadmap

**Project:** Cross-cohort and cross-platform assessment of shortcut risk in blood transcriptomic Parkinson's disease classifiers  
**Author:** Mohammad Javad Asgari  
**Status date:** 2026-10-16  
**Purpose of this document:** a complete, candid hand-off from the original idea through feasibility work, protocol design, implementation, the first pilot, and the remaining execution.

## 1. The original idea

The project began with a methodological question: when a blood-based classifier reports that it distinguishes Parkinson's disease (PD) from controls, how much apparent performance could come from demographic imbalance, cohort identity, assay platform, batch, or blood-cell composition? More importantly, what predictive signal remains when an entire cohort is held out?

The initial plan was to audit suitable cohorts before modeling, register a fixed protocol, compare metadata-only and gene-expression models, test within-cohort and cohort-held-out performance, and then reserve PPMI RNA-seq for external evaluation. The scientific target is not to create the highest-scoring classifier. It is to estimate transportability and report uncertainty without attributing an AUC difference causally to “disease” or “shortcut.”

The project was therefore **feasibility work first, followed by protocol and software development, and only then analytical execution**. It has now crossed into execution: the primary M1–M5 pipeline has been implemented and one real GEO fold has completed. It has not completed the main study analyses or external validation.

## 2. Feasibility work and cohort decisions

### GEO development sources

Two public whole-blood GEO studies were retained for development:

| Source | Planned role | Current protocol population | Platform / notes |
|---|---|---:|---|
| GSE57475 / PROBE | Development source and one direction of cohort-held-out transfer | 142 records: 93 PD, 49 controls | GPL6947 microarray; age and sex available for the planned complete-case set. The operational donor identifier is GSM; independent evidence of unique biological donors is limited. |
| GSE99039 / GENEPARK-TRO | Development source and the opposite transfer direction | 239 records: 120 idiopathic PD, 119 controls | GPL570 microarray; complete age/sex case-control subset. The analysis uses the prespecified TRO subset; 199 other study records were excluded, and missingness depends on batch. |

Together these sources yield 381 planned development records. They are both microarray data, but use different platforms and cohort workflows. This permits within-source and directional leave-one-cohort-out analyses; it does not provide a third independent development cohort. No cross-source donor overlap is claimed absent an adequate crosswalk.

### Candidate and external sources

Other GEO records and literature were reviewed during feasibility. They were not promoted to development cohorts because comparability, phenotype, age/sex completeness, or sample provenance was insufficient for the locked design. Candidate screening is documented in the original working materials, while this repository preserves the final decision and relevant English protocol rather than every exploratory file.

PPMI Project 133 IR3/B38 long-RNA sequencing was chosen as the intended external test because it is an independent RNA-seq source with clinical and cell-related metadata. The participant obtained PPMI access through IDA/LONI Tier 2 approval. The roughly 42 GB archive was downloaded on a separate Windows machine and reviewed there because transferring the full archive was impractical. Aggregate audit outputs were returned for review.

The PPMI audit supports the feasibility of linking archive filenames to release metadata, but identifier parsing alone is not proof of biological identity. The current source documents also flag a discrepancy between portal metadata row counts and the methods release sample count, possible alias/ambiguity resolution, clinical label reconciliation (especially healthy controls), and the need to confirm data-use/publication requirements. Therefore **the final external population and identity crosswalk are not yet locked**. The 311 PD / 161 control count appears in the preregistration-era linkage report but remains provisional until revalidated under the final locked crosswalk and eligibility rules. No individual PPMI records, identifiers, expression counts, or predictions are stored in this repository.

## 3. Feature mapping and preprocessing design

The shared gene universe was built from unique Entrez-to-GENCODE v29 mapping restricted to protein-coding genes on canonical chromosomes. Probe-to-gene ambiguity was not rescued through symbol-only heuristics. For each array, mapped probes for a gene are summarized by their median within sample. The common universe contains 16,364 genes, including 15,720 autosomal genes, 629 chrX genes, and 15 chrY genes.

For the arrays, expression is converted to within-sample fractional midranks. For PPMI featureCounts, the registered design divides counts by reference exon-union length before applying the same sample-wise ranking. Ranking provides a common scale, but it does not guarantee a faithful cross-platform transformation: zero ties and low-expression noise in RNA-seq may behave differently from arrays. Thus the preregistered outcome-free transfer gates (sex and age, plus cell-marker checks where available) are required before a disease-score interpretation. The PPMI count coverage of every mapped feature still needs to be verified against the finalized files.

All learned operations used by the disease pipeline—gene selection, scaling, residualization, tuning, coefficients—must be estimated only in the training partition. The primary PPMI disease transfer is inductive: no PPMI population statistics or outcomes enter the fitted PD pipeline.

## 4. Registered research design

The English protocol and preregistration translation are in `protocol/`. The parent record and amendments are listed in the repository README. The main question is whether adding genes improves external AUROC beyond a strong age/sex and cell-marker baseline when all development is done in GEO. The study does not identify causal fractions of disease and shortcut signals and does not establish clinical benefit.

### Model set

| Model | Features | Purpose |
|---|---|---|
| M1 | Age, age squared, recorded sex, age-by-sex interaction | Demographic-only baseline |
| M2 | Up to 1,000 selected genes | Gene-only comparator |
| M3 | M1 terms plus up to 1,000 genes | Increment over demographics |
| M4 | M1 terms plus six fixed cell-marker scores | Strong baseline accounting for cell-related expression |
| M5 | Exactly M4 terms plus up to 1,000 genes | Primary incremental-gene model |

All models use weighted logistic regression from the locked glum implementation. Baseline coefficients are unpenalized in M3/M5, and the gene block receives the registered penalty. Training weights balance source-by-class mass; the same rule applies to all models and fit scopes. Hyperparameters are chosen by average within-validation-source AUROC, not pooled AUROC. The zero-gene baseline branch is retained; if tuning selects it, the result is reported as “genes not selected in development,” not as a confidence interval proving zero gene value.

### Evaluation and uncertainty

The prespecified primary GEO analyses comprise nested donor-separated CV within each source, both directional cohort-held-out transfers, and final pooled GEO fits. The intended primary PPMI evaluation compares M5 against M4 with paired, class-stratified donor bootstrap uncertainty (5,000 repetitions), using fixed models and no PPMI tuning. AUROC is primary; AUPRC, calibration-related measures, sample sizes, prevalence, and uncertainty are reported as descriptive complements under the protocol. A prespecified delta of 0.05 is used for precision interpretation; a non-significant estimate is not treated as evidence of equivalence.

Controls and sensitivities include metadata baselines, sex/age transfer gates, cell-marker validation (CBC where sufficiently linked), chromosome-sex-gene removal, matched/subsample analyses, training-only residualization, and ComBat inside-fold versus out-of-fold comparisons. The protocol restricts scope; exploratory changes after seeing outcomes must be labeled and cannot replace the confirmatory analysis.

## 5. Numerical validation and software implementation

Before real-data execution, the weighting and block-penalty implementation was checked against independent numerical routines. KKT scaling and thresholds were calibrated on synthetic data, and the solver was stress-tested at 1,000 genes. The synthetic pipeline suite passed 19 checks, including repeated-donor split integrity, training-only preprocessing, held-out mutation invariance, the M4/M5 baseline branch, paired bootstrap behavior with ties, and finite 1,000-gene preprocessing. These checks establish software properties only; they do not establish biological signal or bootstrap confidence-interval coverage. Coverage validation remains on the roadmap.

The repository includes the executable M1–M5 GEO engine, paired donor bootstrap, frozen configuration, dependency lock, split manifest, aggregate audits, and English provenance materials. The pipeline consumes locally prepared GEO matrices that are intentionally not checked in. A run requires the matching public implementation receipt; no receipt should be fabricated or bypassed.

## 6. What has actually been executed on study data

One real-data pilot scope completed: `GSE57475/outer/0`. It used five inner folds, 113 training records, 29 held-out records, and all five models. It took 427.34 seconds wall time, made 269 fit attempts, produced zero KKT failures/retries/warnings, and had maximum KKT residual 3.8924e-8. Test AUROCs were M1 0.3658, M2 0.5105, M3 0.3658, M4 0.3947, and M5 0.4895. M3 reused M1 in this fold. The detailed aggregate report is in `results/geo_pilot/GEO_PILOT_REPORT_20261016_EN.md`.

This is an **operational pilot fold inside the planned CV**, not an independent test and not a final estimate. It was not used to change the locked method. The complete GEO nested CV, both directional transfers, final GEO refit, sensitivity suite, and PPMI analyses have not been run. The pilot's predicted values, fitted coefficients, preprocessing arrays, and per-fit logs are deliberately excluded from the repository; they remain local in the analysis workspace.

## 7. What remains: execution roadmap

### Phase A — complete GEO development and internal validation

1. Verify the local prepared inputs against the committed input and split hashes; document exact software and hardware versions.
2. Run all remaining outer folds for GSE57475 and GSE99039, preserving one fresh output directory per locked scope.
3. Run both directional cohort-held-out transfers and the prespecified final pooled GEO model fitting.
4. Execute the registered sex/age transfer gates and bounded sensitivities; retain per-fit convergence diagnostics.
5. Validate donor-bootstrap CI coverage under the registered design, and finalize the permutation scope/compute plan. Any required method change must be timestamped and publicly amended before affected results are viewed; do not silently reduce permutations or omit label-dependent refitting.
6. Produce aggregate tables, uncertainty intervals, failure accounting, and a reproducible GEO results report. Audit that no split, transformation, tuning, or code path leaks held-out data.

**Milestone:** after all GEO CV/transfer analyses and required diagnostics are complete and quality-controlled, create a Git release/tag and publish a Zenodo snapshot. No incremental Zenodo upload is planned before this milestone.

### Phase B — finalize external PPMI eligibility without scoring disease

1. Resolve the archive-to-release crosswalk using the release guide and all relevant sample-name segments, including PATNO and visit; reconcile aliases and ambiguous links.
2. Validate one-to-one file-to-metadata mapping, subject uniqueness at the selected baseline visit, label concordance, QC, age/sex completeness, gene coverage, and exclusion rates by diagnosis.
3. Verify the permitted publication and acknowledgement requirements under the current PPMI data-use terms.
4. Finalize the external population definition and all private PPMI crosswalk inputs locally. Do not upload individual data or identifiers to GitHub or Zenodo.

### Phase C — second public lock, then external validation

Before **any** PPMI score, publish the second locked package: code and environment hashes, exact split manifest, finalized M1–M5 models and preprocessing states, control thresholds, and a hash/reference for the finalized private population/crosswalk. The public package must not contain subject-level PPMI information. Update Zenodo only at this specified milestone, and preserve the matching Git commit/tag.

Then run the outcome-free sex, age, and permitted cell-marker transfer gates. If a gate fails without a demonstrable technical bug, report the confirmatory PPMI disease test as unevaluable; do not alter the transformation after seeing the failure and then present the result as confirmatory. If the gates pass, run the fixed inductive M1–M5 PPMI transfer once, calculate the prespecified metrics and paired uncertainty, and preserve the audit trail.

### Phase D — interpretation and manuscript

Only after quality control and analysis should the paper be written. The primary claim must be bounded to this input, population, and model family. Report uncertainty and inconclusive results transparently. Avoid saying that a measured AUC gap is a causal percentage of shortcut signal or that successful cross-platform prediction proves clinical utility. Discuss single-pilot scope, sample-size precision, donor provenance, treatment context, missingness, gene mapping, cell mixture, and platform shift as limitations.

## 8. Git and Zenodo practice from this point

GitHub is the living project record. Each meaningful code, protocol, or results milestone should be committed, with English documentation and aggregate outputs only. Do not commit raw expression matrices, individual PPMI data, personal credentials, or local result files containing predictions or fitted subject-level artifacts. Keep the full internal audit locally under access controls.

Zenodo is reserved for the two project milestones described above: (1) complete, quality-controlled GEO CV package; (2) final pre-PPMI-scoring lock. The existing Zenodo records document preregistration/provenance history; they are not evidence that the full analyses or paper are complete.

## 9. Plain-language answer: is the project now “being done”?

Yes, the project has moved beyond feasibility: the analysis protocol is registered, the primary software exists, synthetic and numerical checks passed, and one real GEO fold ran. But the scientific study is still at an early execution stage. Most planned GEO analyses and all PPMI external validation remain outstanding. The next real research work is completing the frozen GEO pipeline—not writing more feasibility documents and not claiming results from the single pilot fold.
