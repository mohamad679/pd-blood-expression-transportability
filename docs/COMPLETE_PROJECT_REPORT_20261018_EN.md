# Complete project report: PD blood-expression transportability

**Author:** Mohammad Javad Asgari  
**ORCID:** [0009-0002-2313-571X](https://orcid.org/0009-0002-2313-571X)  
**Status date:** 18 October 2026, Asia/Tehran  
**Repository:** `mohamad679/pd-blood-expression-transportability` (private GitHub repository)  
**Purpose:** a detailed account of the original idea, feasibility decisions, data access and audits, protocol evolution, implementation, numerical checks, completed study analyses, actual results, unresolved issues, and the remaining execution roadmap.

This is a retrospective project-status report. It is not a new preregistration, a protocol amendment, or a claim that external validation has been completed. Preparing it did not involve fitting new disease models or calculating any PPMI score. Numeric tables below summarize preserved results; the methods remain governed by the registered source documents and amendments.

## 1. Current position

The project began with feasibility work and has now entered substantive execution. The core GEO M1–M5 pipeline is implemented, registered, and executed across all **13 planned scopes**: five outer folds within each of two GEO sources, two directional cohort-held-out transfers, and final pooled development fitting.

The development dataset contains **381 records**: 142 in GSE57475/PROBE and 239 in GSE99039/GENEPARK-TRO. The shared feature universe contains **16,364 protein-coding genes**. All 3,146 post-pilot fit records passed the fixed numerical acceptance criterion; the earlier pilot separately passed 269 fits. Six prediction/tuning sigmoid-clipping warnings were observed and require a focused review.

Within-source mean outer-fold AUROC for M4/M5 was **0.5615/0.5911** in PROBE and **0.6630/0.6861** in GENEPARK-TRO. Transfer M5−M4 AUROC was **−0.1784** from PROBE to GENEPARK-TRO and **−0.0086** in the reverse direction. These are actual GEO results, not synthetic scenarios or PPMI results.

The core GEO milestone was publicly released as [Zenodo version 3.1.3-geo-cv, DOI 10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851). No PPMI predictive or cell-marker score has been computed. Development transfer controls, several diagnostics/sensitivities, final PPMI eligibility and the second public lock remain outstanding.

### 1.1 A reporting correction found while preparing this account

The previously released short GEO report correctly recorded final M5 as `lambda=0.001, rho=1.0` but incorrectly called it "ridge-only." **Rho is the L1 mixing fraction: rho=1 means lasso/L1, whereas rho=0 means ridge/L2.** The final M5 gene block is L1-penalized and has 210 nonzero coefficients. This is a terminology correction, with no change to the code, fitted models, hyperparameters or results. The original released bytes are preserved; see [REPORTING_CORRECTION_20261018_EN.md](REPORTING_CORRECTION_20261018_EN.md).

## 2. The original scientific idea

The initial question was whether high reported accuracy of blood gene-expression classifiers for Parkinson's disease reflects disease-related information, demographic imbalance, technical batch/platform differences, or other shortcuts. Proposed checks included metadata-only classifiers, within-cohort CV versus cohort-held-out evaluation, train-only confounder handling, a batch-correction leakage comparison, and a sealed external cohort.

Early hypotheses were that internal CV might be optimistic, metadata might discriminate cases from controls in some cohorts, transfer might deteriorate, and a strong baseline might absorb much of the apparent gene advantage. Those were hypotheses rather than prespecified desired results. They could not legitimately determine which results were accepted.

Two important refinements followed:

1. An AUROC difference cannot identify a causal percentage attributable to disease versus shortcut. The final question concerns predictive value and transportability of a fixed procedure.
2. A negative RNA-seq transfer result can reflect biological differences, measurement/processing differences, medication, sample selection or the transformation itself. Positive controls and limited claims are necessary for interpretation.

The final primary question is therefore: **when the models are developed exclusively in GEO, does adding a gene block improve PPMI AUROC beyond an otherwise identical age/sex and six-marker-panel baseline?** The project is a methodological evaluation, with positive, negative and inconclusive outcomes all legitimate.

## 3. Feasibility and novelty assessment

### 3.1 What the literature review changed

The initial discussion considered whether poor external replication itself would be novel. A central prior study was Dayan et al., [Inherent variability limits clinical utility of reproducible Parkinson's transcriptomics signatures](https://www.nature.com/articles/s41531-025-01238-y), published online on 19 December 2025 in *npj Parkinson's Disease*, volume 12 (2026). Its abstract reports evaluation of 13 published signatures, median PPMI AUC of 59.7%, and a prospective cohort of 30 participants. Thus simply demonstrating lower external AUC would have substantial precedent.

The proposed contribution became a fixed cross-platform procedure with a strong shared baseline, explicit source/class weighting, transfer controls, numerical validation, transparent failure rules and precision-based reporting. This report does not claim a fresh exhaustive novelty search. Reconstructing the original 13 published classifiers was not selected as the current study design, and no claim is made that their original coefficients were recovered.

A Q1/Q2 submission remains a possible target, not a result established by feasibility, sample count, a DOI or the present GEO findings. Publication prospects depend on the remaining validation, interpretation, novelty relative to the final literature review and peer review.

### 3.2 Initial GEO screening

The initial go/no-go criterion required at least three usable, comparable cohorts with individually linked age and sex. The local audit covered the following GEO series. "Complete age/sex" below refers to available metadata records, not necessarily eligible independent whole-blood donors.

| GEO series | All series records | Platform | Records with both age and sex | Decision / limitation |
|---|---:|---|---:|---|
| GSE18838 | 28 | GPL5175 | 0 | Did not meet covariate completeness criterion. |
| GSE202665 | 59 | GPL20844 | 59 | Activated CD4+ T cells with repeated time points; nine title-derived donors. Not a comparable whole-blood cohort. |
| GSE49126 | 50 | GPL4133 | 0 | Did not meet covariate completeness criterion. |
| GSE57475 | 142 | GPL6947 | 142 | Retained as PROBE development cohort. |
| GSE6613 | 105 | GPL96 | 0 | Initial sealed-cohort candidate was not retained for this covariate-adjusted design. |
| GSE72267 | 59 | GPL571 | 0 | Did not meet covariate completeness criterion. |
| GSE99039 | 558 | GPL570 | 284 across all diagnoses | Only 239 IPD/control complete cases met the final development definition; retained as GENEPARK-TRO. |

The original three-GEO-cohort criterion was not met. The design was explicitly changed to **two development sources plus an independent PPMI RNA-seq test**. It is not described as LOCO across three or more development cohorts. The 59 activated-cell records were not treated as 59 independent comparable whole-blood donors.

## 4. Selected GEO populations

| Property | GSE57475 / PROBE | GSE99039 / GENEPARK-TRO |
|---|---:|---:|
| Development records | 142 | 239 |
| PD / idiopathic PD | 93 | 120 |
| Controls | 49 | 119 |
| PD fraction | 65.49% | 50.21% |
| Female / male, total cohort | 54 / 88 | 147 / 92 |
| Platform | GPL6947 | GPL570 |
| Probes in processed source matrix | 48,766 | 54,675 |
| Source preprocessing inherited by the study | GenomeStudio average normalization | RMA/log2/quantile normalization |
| Technical metadata | 16 title-derived chip proxies | 11 provided batches in the selected set |
| Operational donor identifier | GSM, without independent public donor-ID confirmation | 239 explicit unique subject IDs |

The full GSE99039 series includes 205 IPD and 233 controls, for 438 target-diagnosis records before complete-case selection. Retaining 120 IPD and 119 controls excludes 199 target-diagnosis records. The remaining 120 series records have other phenotypes. All retained target complete cases belong to the TRO subset, so **GENEPARK-TRO** is the appropriate analysis label. Missingness is related to batch/site, and this complete-case population is not the full multicentre study.

PROBE GSM separation is an operational split safeguard, not independent proof of biological donor uniqueness. The absence of cross-source donor overlap has not been established by a full crosswalk. These limitations must remain visible in the manuscript.

## 5. PPMI access, download, archive audit and linkage recovery

### 5.1 Access route and why the large archive was used

PPMI was initially unavailable to the participant. The participant subsequently obtained IDA/LONI access and reviewed the separate Tier 2 request form and access guidelines after the coordinator's email directed them to those documents. This opened an access route separate from an AMP-PD institutional agreement. It does not replace the continuing duty to verify the relevant use/publication conditions before submitting the paper.

The participant supplied the RNA-seq release methods PDF, PPMI user guide, portal metadata CSV, age-at-visit CSV and participant-status CSV. Metadata export and download behavior were initially inconsistent. The approximately 42 GB IR3 archive was downloaded on the participant's remotely accessed Windows computer in Italy because the current local connection made the download impractical.

That Windows activity belongs to data acquisition and feasibility history. Current GEO analytical execution is on macOS, and future user-facing shell instructions are to use macOS/bash rather than PowerShell.

### 5.2 Audit sequence and what each stage established

| Stage | Actual evidence | What it did not establish |
|---|---|---|
| Portal metadata review | 4,850 portal rows; provisional baseline, phenotype and covariate filters | No count-file link and no external model result. Early age field was age at consent, not verified sampling age. |
| Early archive metadata audit | Release metadata found; 4,871 rows; technical/QC fields present | The stream stopped early; it did not verify the whole archive or gene coverage. |
| Structure preview | Three seven-column featureCounts-style tables with gene and annotation fields | Header equality alone did not resolve sample linkage. |
| Complete one-pass audit | Archive end reached; all 4,871 count tables examined; private cache created | The original matching rule still linked zero count files to release rows. Structural success was not linkage success. |
| Guide-informed recovery from cached data | 4,871 unique file-to-release-row links: 4,752 header/filename agreements and 119 exact header-only links | Computational matching did not prove biological identity or finalize eligibility. |

The report from the complete user-run archive audit states:

- Archive size: **42,893,354,202 bytes** (approximately 42.89 GB decimal / 39.95 GiB).
- Archive members seen: **14,615**; count tables: **4,871**; `quant.sf` files: **9,742**.
- Every count table had **58,780 gene rows**.
- Every table contained all **16,364** reference genes, with reported length agreement fraction **1.0**.
- Release QC: **4,756 pass** and **115 fail**.
- The private cache was approximately **187.7 MB**; its purpose was to avoid repeatedly scanning the large archive.
- The all-once scan took **2,846.5 seconds**, approximately 47.4 minutes, on the remote computer.

These are results reported by the participant-run aggregate audit. The complete raw archive was not independently rescanned by the current macOS GEO pipeline. The equality between the methods figure of 4,756 and the audit's QC-pass count is informative; it does not by itself reconcile every record across the 4,850-row portal and 4,871-row release.

### 5.3 The identifier-parsing problem

An early strategy considered only a prefix before a dot. The PPMI guide describes composite identifiers in which PATNO and visit information can occur in middle segments. Throwing away those segments prevented valid matching.

The revised recovery rule used exact dot-delimited prefixes against documented release identifiers and, when needed, the guide's PATNO/visit composite format. When both header and filename resolved, they had to agree. No expression-based reassignment was used. The revised aggregate report recorded one count table for every release row and no multiply-linked release rows.

### 5.4 Provisional external population

The count-linked candidate set is **311 idiopathic PD and 161 healthy controls**, with distinct PATNO counts equal to candidate row counts. It is based on baseline visit, Long RNA/QC, valid covariates and agreement of historical labels. The earlier portal-only healthy-control count of 162 became 161 after the fuller age/eligibility linkage.

These **472 candidates remain provisional**, rather than a finalized externally validated population. Final one-to-one eligibility, aliases/ambiguity, reference-version consistency, age at sampling, diagnosis concordance, QC exclusions by group, treatment at sampling and population-file hashes still need to be documented and locked. A current 2026 status must not silently replace the historical baseline diagnosis or create future-outcome exclusion.

The audit reports plate, box, RIN and supplied neutrophil-score metadata. A supplied neutrophil score is not a measured CBC differential. Individual PPMI data, counts and private crosswalks have not been published in GitHub or Zenodo. The PPMI private expression cache is not assumed to be currently available on the macOS analysis machine simply because an aggregate remote audit exists.

## 6. Gene mapping and prepared input construction

### 6.1 Fixed annotation rules

The gene universe was constructed from GEO platform annotation and GENCODE v29 primary-assembly annotation using a single-integer Entrez identifier that mapped uniquely to a protein-coding gene on a canonical chromosome. Ambiguous mappings, absent mappings, non-protein-coding genes and noncanonical genes were excluded. There was no post hoc symbol/alias rescue based on performance.

| Mapping quantity | GPL6947 / PROBE | GPL570 / GENEPARK-TRO |
|---|---:|---:|
| Source probes | 48,766 | 54,675 |
| Retained mapped probes | 27,168 | 37,730 |
| Unique retained genes before intersection | 17,905 | 17,039 |
| Genes represented by multiple retained probes | 6,130 | 10,018 |

The intersection is **16,364 genes**: **15,720 autosomal**, **629 chromosome X** and **15 chromosome Y**. Reference exon-union lengths are available for all common genes. Annotation/reference coverage was checked separately from actual IR3 count-table coverage; those are different evidence stages.

### 6.2 Transformation

For microarrays, retained probes for each gene are summarized by their median within the sample. Genes are then transformed to fractional midranks, `(midrank - 0.5) / G`, with the fixed common universe `G=16,364`.

For RNA-seq, the registered plan first divides featureCounts by fixed reference exon-union length, then uses the same within-sample ranking. No target-cohort mean or target label is required for that step. Ranking is a common numerical scale; it is not proof that microarray and RNA-seq measurements are biologically interchangeable, especially with zero ties and weakly expressed genes.

The prepared arrays store float32 ranks. The implementation amendment reconstructs their fixed integer/half-integer rank lattice with `round(rank_float32 * 2G) / (2G)`, checking the reconstruction tolerance and samplewise mean of 0.5. It does not rerank rounded values and manufacture new ties.

The study uses already processed source matrices. No new whole-dataset cross-sample normalization was added during preparation, but inherited GenomeStudio/RMA processing remains a limitation. The analysis does not claim to begin from raw CEL or FASTQ files.

### 6.3 Inputs and hashes

Four prepared GEO files are available locally in the ignored `prepared/` directory. Their SHA-256 values were checked against the frozen audit before model execution:

| Input file | SHA-256 |
|---|---|
| `GSE57475_development.npz` | `77e4900a42c6b1a34857548550b5deb3c401593644d5dc9d9cc280b01594cf14` |
| `GSE57475_sample_metadata.json` | `ae42df5e7066a8f66ebb57256886f1355e6f9419c2966c84681baf22c434a18e` |
| `GSE99039_development.npz` | `d1cfa7f4e8556d57099ce85b754e25e22b2c114335abb4eb16d9259eb26a5d49` |
| `GSE99039_sample_metadata.json` | `bebcf13f321214081b5195527a9d1a1f3100ce8c1e81138a54a1fab94d4e8414` |

Prepared expression matrices and sample-level metadata are intentionally excluded from the Git history and milestone upload. The four-file `plan` check validated 381 records, 16,364 genes and the frozen split hash; it did not fit models.

## 7. How the protocol evolved

The design underwent explicit review before the main GEO execution. The important changes were substantive decisions rather than additional classifiers chosen to maximize AUC:

| Initial issue | Final decision |
|---|---|
| Causal "percentage disease versus shortcut" wording | Restrict claims to predictive incremental value and transfer. |
| Three-GEO-cohort go/no-go criterion could not be met | Register two development cohorts plus one independent PPMI test. |
| Age/sex alone insufficient as a comparator | Add M4 with six fixed cell-marker scores and compare M5 against the same baseline. |
| Penalty on all features confounded baseline comparison | Keep every baseline coefficient unpenalized; explicitly penalize only the added gene block. |
| Platform identity could exploit differing cohort prevalences | Balance source×class mass in every fit and tune using within-source macro-AUROC. |
| Gene augmentation might fail to help even in tuning | Include a baseline-only candidate in M3/M5; preserve the identical-model branch. |
| Cross-platform negative result could be uninterpretable | Register sex/age controls and limited CBC/identity diagnostics with explicit consequences. |
| A fixed probability threshold could reject a well-ranked sex model | Use sex AUROC and CI for the gate; BA at 0.5 is descriptive. |
| Numerical warnings and solver behavior could alter selection | Define normalized KKT, one fixed cold retry and analysis-level failure rules. |
| Multiple CI options could permit post hoc switching | Use one confirmatory paired stratified bootstrap method. |
| Insufficient precision for small increments | Predefine delta=0.05 as a research-relevance margin and accept inconclusive results. |
| Iterative review could expand scope indefinitely | Freeze a bounded list of controls/sensitivities; new checks require a material documented failure. |

The initially proposed pathway-informed model is not part of the locked five-model study. ComBat is a bounded development diagnostic rather than the primary external transformation. The original Holm idea does not create multiple primary tests: the final protocol contains one confirmatory comparison.

## 8. Final model and fitting specification

### 8.1 Five disease models

| Model | Unpenalized terms | Penalized terms | Role |
|---|---|---|---|
| M1 | Intercept, standardized age basis, age squared, recorded sex, age×sex | None | Demographic baseline. |
| M2 | Intercept | Up to 1,000 selected genes | Gene-only comparator. |
| M3 | Exactly the M1 baseline | Up to 1,000 selected genes | Increment beyond demographic predictors. |
| M4 | M1 terms plus six fixed marker scores | None | Shared demographic/cellular baseline. |
| M5 | Exactly the M4 baseline | Up to 1,000 selected genes | Primary gene-increment model. |

All five use a binomial/logit family. M1 therefore contains a four-term demographic basis, rather than only two raw age/sex columns. The baseline has 10 columns in the final M4/M5 fit. An unpenalized intercept is additional to those columns. Constant/dependent baseline columns are removed using a deterministic training-only rank check; M4 and M5 retain the same baseline design within a partition.

### 8.2 Objective and hyperparameters

For gene coefficients gamma and baseline coefficients beta, the fitted objective is:

```text
sum_i w_i * logistic_NLL(y_i, intercept + B_i beta + G_i gamma) / sum_i w_i
    + lambda * rho * sum_j |gamma_j|
    + lambda * (1-rho) * sum_j gamma_j^2 / 2
```

`P1=P2=0` for baseline columns, and `P1=P2=1` for gene columns. No artificial multiplication of baseline columns is used. Adding genes can change fitted baseline coefficients because the predictors are correlated; "unpenalized" does not mean numerically fixed.

- Lambda grid: `0.0001, 0.001, 0.01, 0.1, 1`.
- Rho / L1 fraction: `0, 0.5, 1` (ridge, elastic net, lasso respectively).
- M3/M5 also include the no-gene candidate.
- AUROC tie tolerance: `1e-6`; ties prefer no genes, then larger lambda, then larger rho.
- A selected nonzero gene candidate at a grid boundary permits one-decade extension toward that boundary, for at most two rounds on the same inner splits and metric.
- A remaining boundary is reported; it does not authorize unlimited search or external-data tuning.

### 8.3 Source/class weights and effective sample size

For each current training partition, `w_i = n_train / (2 * K_train * n_train[source_i, class_i])`. This balances the total mass of each training source×class cell. Weights sum to the training sample count. The rule applies to M1–M5, baseline fits, inner fits, outer refits, single-source transfers and final pooled fitting. For one source it balances the two classes.

Tuning averages unweighted AUROC within each validation source and then averages across inner folds. Training weights are not reapplied to held-out evaluation metrics. Preprocessing still uses unweighted training means/SDs.

Without this correction, cohort prevalence alone yields an aggregate source-only AUROC of `1951/3408 = 0.5724765` for the observed development group counts. This is a count-based diagnostic, not an observed gene-model score or a causal shortcut fraction.

Kish weight-dispersion summaries are **128.3662/142** in PROBE, **238.9958/239** in GENEPARK-TRO and **334.0463/381** in the pooled fit. Kish is not an exact information-based sample size for logistic regression or AUROC.

### 8.4 Features, folds and leakage boundaries

Gene selection uses positive median absolute deviation computed within each current training partition. For multisource training, the minimum within-source MAD ranks each gene. Gene-ID ties are broken deterministically; at most 1,000 genes are selected. Selection and scaling are repeated inside every inner fold.

Age is standardized from training data before polynomial/interaction construction. Constructed baseline columns and selected genes are then centered/scaled using that training partition only. Cohort ID, clinical severity, olfaction and technical IDs are excluded from the five disease-model predictors.

There are five outer and five inner folds within each source, shared across models. Donors stay together; stratification uses class, or cohort×class for pooled inner CV. Master seed is **20261014**; scope-specific seeds are derived from the first four big-endian bytes of SHA-256 of the master seed plus scope. Realized splits, not only the seed, are frozen.

Primary external transfer is inductive: fitted statistics come from GEO only. No PPMI population centering, label or prevalence enters the PD model. Each run uses a fresh output directory and an exact code/split registration receipt.

## 9. The six fixed marker panels

| Panel | Usable shared genes | Reference role |
|---|---:|---|
| T cells | 14 | Fixed MCP-counter-derived panel. |
| B lineage | 8 | Fixed MCP-counter-derived panel. |
| Monocytic lineage | 7 | Fixed MCP-counter-derived panel. |
| Neutrophils | 15 | Fixed MCP-counter-derived panel. |
| Erythroid | 11 | Consensus from 12 independent xCell reference signatures. |
| Platelet | 118 | Mapped consensus from three xCell reference signatures. |

The score is the mean within-sample rank across each fixed panel. These are marker scores, not measured cell fractions and not full executions of the MCP-counter or xCell algorithms. They can reflect activation and specimen processing as well as cell abundance. Erythroid signal may depend on globin processing and act as a platform proxy. The marker baseline cannot be described as complete control of blood-cell composition.

## 10. Synthetic and numerical validation already performed

These tests used synthetic data. Their success established software/numerical properties, not real PD signal, actual statistical power or PPMI transfer.

| Verification | Recorded result |
|---|---|
| Original block-penalty checks | 23 passed; a deliberately wrong penalty mask was detected. |
| Uniform-weight regression comparisons | 66 passed, including common weights 0.1, 1 and 10. |
| Nonuniform weights versus independent SciPy objective | Nine comparisons passed. |
| KKT scale/independent objective checks | 48 passed across 40/1,000 genes and separate calibration/check seeds; deliberate factor-of-two loss-gradient error detected. |
| Maximum glum / independent-SciPy KKT in the scale report | Approximately 1.8293e-10 / 3.3167e-8. |
| Maximum independent-objective difference | Approximately 6.4282e-14. |
| 1,000-gene numerical stress fixture | 32 final-like fits passed, with no unresolved KKT failure or retry. |
| Stress warnings | 15 messages in three fits; all passed final KKT. |
| Stress maximum final KKT | 6.1943e-9. |
| Stress runtime / slowest fit | About 1,454.45 s / 532.88 s. |
| Re-executed numerical null | 200 paired repetitions; 35,100 fit attempts; no unresolved KKT failure or retry. |
| Null warnings | 975 messages in 227 fits across 84 repetitions; warned fits passed final KKT. |

The first null run had 1,394 warnings and was retained. An initial slow stress run was interrupted after four completed fits; an unfinished computation was not treated as a failed KKT result. Synthetic evidence motivated a documented solver dispatch: `irls-ls` for positive-lambda pure L2, and `irls-cd` for L1/mixed/baseline fits. The lambda grid and acceptance threshold were not restricted to obtain a pass.

### 10.1 Cohort-shortcut null results

The fixed synthetic null contained two sources with different prevalences and a strong platform effect, but no within-source PD information. It compared unweighted pooled tuning with source/class weighting plus macro tuning:

| Synthetic result | Unweighted / pooled criterion | Weighted / macro criterion |
|---|---:|---:|
| Gene block selected | 200/200 (100%) | 130/200 (65%) |
| Wilson CI for selection proportion | 98.12–100% | 58.16–71.27% |
| Mean pooled test AUC | 0.5520096 | 0.5051516 |
| MCSE of pooled test AUC mean | 0.0016767 | 0.0020330 |
| Mean within-source macro test AUC | 0.5037540 | 0.5019759 |
| Grid expansions triggered | 115/200 | 47/200 |
| Selected gene candidate remaining on a boundary | 0/200 | 0/200 |

The observed gene-selection-rate reduction was 35 percentage points, with an approximate Monte Carlo interval of 28.37–41.63 points. This does not mean a 35% probability of baseline-only selection in the real study, nor does it establish the complete procedure's Type-I error. The independent synthetic test came from the same two source distributions; a third-platform synthetic transfer experiment was not performed.

### 10.2 Numerical acceptance and failure rule

Each fit must have finite output and maximum absolute KKT residual **≤1e-6** in the input coefficient coordinates, with the loss divided by total weight and the same P1/P2 terms used by the fitted objective. Exact-zero gene coefficients use a subgradient condition; baseline and intercept gradients are unpenalized. The L1 norm is also recorded.

Initial settings are `max_iter=1000`, `max_inner_iter=100000`, `gradient_tol=1e-7`, `scale_predictors=False`. A failed KKT check permits one cold retry at `5000`, `1000000`, `1e-8`, with the same data, weights and hyperparameters. A persistent invalid candidate invalidates its dependent tuning/analysis; candidates or folds are not silently discarded. Final M4/M5 numerical failure makes the confirmatory comparison unevaluable. Warning messages alone do not necessarily imply numerical failure, and KKT acceptance is not an AUROC-error bound.

## 11. Pipeline implementation and synthetic leakage checks

The repository contains the M1–M5 loading/preprocessing/training/scoring engine, realized splits, paired bootstrap utilities, preparation utility, numerical backend, protocol/configuration, dependency lock and aggregate reports. `geo_pipeline.py` currently exposes `plan`, `pilot` and `analysis` commands. It has no PPMI loader/scoring command and does not implement every registered auxiliary analysis merely because those analyses appear in the protocol.

The pipeline synthetic suite passed **19 checks**, including:

1. Repeated-donor separation across outer folds and inner partitions; all rows held out once.
2. Held-out expression/age/label perturbations leaving training preprocessing unchanged.
3. Gene standardization equal to training-only mean/SD, and fixed samplewise marker calculation.
4. Identical unpenalized M4/M5 baseline design and baseline-preferring candidate tie behavior.
5. Held-out perturbations leaving tuning decisions and fitted coefficients unchanged.
6. Inner-validation perturbations leaving their own training state unchanged.
7. Bit-identical M5/M4 predictions in the baseline-only branch.
8. Bootstrap multiplicities agreeing with an independent AUROC calculation with ties.
9. Swapping model labels exactly reversing paired bootstrap differences; identical scores producing zero differences.
10. Repeated visits rejected as independent bootstrap donors; finite training-scaled 1,000-gene design.

These tests are distinct from an audit of every real-data fit's membership. The pilot had a real membership audit; a complete post-bulk membership/leakage audit remains on the roadmap. Bootstrap algebra passed, but confidence-interval coverage simulation is unfinished.

## 12. Registration, translations and transition into real execution

| Event | Public/local evidence and timing |
|---|---|
| Initial feasibility and metadata review | Began 13 October 2026; not study-model training. |
| Preparation and remote archive audits | Mapping/rank features and structural/linkage checks; no PPMI model scores. |
| Original v3.1 registration | [10.5281/zenodo.23405371](https://doi.org/10.5281/zenodo.23405371). |
| Implementation amendment 3.1.1 | [10.5281/zenodo.23405802](https://doi.org/10.5281/zenodo.23405802); public verification recorded at 2026-10-15T22:07:25.879272Z. |
| First real study-model fit | PROBE outer/0 began 2026-10-15T22:07:32Z, corresponding to 16 October in Tehran. This was the transition from feasibility/software validation to real analytical execution. |
| First English translation | [10.5281/zenodo.23406064](https://doi.org/10.5281/zenodo.23406064); published after the pilot began. |
| English provenance correction | [10.5281/zenodo.23411084](https://doi.org/10.5281/zenodo.23411084), version 3.1.2-en. |
| GitHub working repository | Created and populated with English code/documentation on 16 October; meaningful subsequent progress committed. |
| Input/environment preflight | Four input hashes verified; user successfully ran `plan` on macOS after package-download timeouts. |
| Bulk GEO execution | Remaining 12 scopes executed 17–18 October 2026. |
| Core GEO milestone release | [10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851), version 3.1.3-geo-cv, published 18 October. |

The earlier English translation did not precede the real-data pilot. Its status statements describe the original snapshot. Some first-English references incorrectly named later EN files as documents that had been saved/hashed before synthetic outputs; the provenance correction identifies the actual FA originals and their hashes. The original Persian preregistration and published amendments govern if translations conflict. The current report is English and introduces no silent change of reference method.

Frozen configuration booleans such as `training_performed=false` and `public_preregistration_done=false` describe the freeze-time state; they are not current status flags to edit retrospectively. Current execution/publication status is maintained in `docs/STATUS.md` and the execution/publication receipts.

## 13. Real GEO execution and audit

### 13.1 Environment and operational issues

Recorded runtime: Python **3.11.15**, macOS **15.7.9**, x86_64, `glum==3.4.1`. Core locked versions include NumPy 2.4.6, SciPy 1.17.1, scikit-learn 1.9.1, pandas 3.0.6, tabmat 4.2.1 and threadpoolctl 3.7.0. The complete 20-package version list is preserved in `requirements_glum_v3.lock.txt`.

The default PyPI file host timed out. The user successfully used the Tsinghua PyPI mirror for the no-fit plan check. The agent's environment initially had cache/network restrictions; after installing the same pinned versions, successful offline resolution was used for subsequent scopes. No dependency version, solver rule or model grid was changed to solve the network problem.

Models were fitted serially without warm starts or simultaneous model fits. The recorded OpenMP libraries used up to eight internal threads. Summed process CPU time therefore exceeds elapsed wall time; neither number is called single-core wall time. GPU acceleration was not used.

### 13.2 Completed scopes, fit counts and wall time

| Scope | Training n | Test n (PD) | Fit attempts | Wall seconds |
|---|---:|---:|---:|---:|
| `GSE57475/outer/0` | 113 | 29 (19) | 269 | 427.34 |
| `GSE57475/outer/1` | 113 | 29 (19) | 238 | 201.47 |
| `GSE57475/outer/2` | 114 | 28 (19) | 269 | 412.20 |
| `GSE57475/outer/3` | 114 | 28 (18) | 269 | 371.57 |
| `GSE57475/outer/4` | 114 | 28 (18) | 285 | 423.19 |
| `GSE57475/transfer` | 142 | 239 (120) | 330 | 605.88 |
| `GSE99039/outer/0` | 191 | 48 (24) | 270 | 413.41 |
| `GSE99039/outer/1` | 191 | 48 (24) | 240 | 357.59 |
| `GSE99039/outer/2` | 191 | 48 (24) | 255 | 776.65 |
| `GSE99039/outer/3` | 191 | 48 (24) | 255 | 417.66 |
| `GSE99039/outer/4` | 192 | 47 (24) | 255 | 454.22 |
| `GSE99039/transfer` | 239 | 142 (93) | 240 | 609.53 |
| `pooled/final` | 381 | 0 (no test) | 240 | 956.12 |

The earlier pilot is counted once, as outer/0 in PROBE. It was not treated as a separate independent test. Twelve new directories plus that pilot cover 13 scopes. The total is **3,415 fit attempts**: 269 pilot plus 3,146 post-pilot. Sum of scope wall times is approximately **6,426.81 seconds (1.79 hours)**; summed process CPU is approximately **26,427.55 seconds (7.34 hours)**. Setup, waiting and documentation time are not included in that sum.

### 13.3 Numerical outcome and warning limitation

All 3,146 new records passed KKT on the first attempt. There were zero logged failed records, retries or exceptions. Maximum KKT was **4.10923814021108e-8**, below 1e-6. The pilot separately passed with maximum **3.8924002795232004e-8**.

The fit-level warning fields were empty. Nevertheless, six `Sigmoid function too close to 0 or 1. Clipping.` messages were visible on stderr during validation-prediction/tuning calls. The pipeline's warning capture surrounds fitting/KKT, while inner validation calls `fit.predict()` outside that capture. Exact candidate attribution was therefore not persisted. Final held-out prediction uses the saved coefficients and an explicit logistic transform. This explains the logging gap; it does not establish that clipping had no effect on candidate rankings or metrics. A focused review remains necessary.

Separate Arrow messages about inability to read processor-cache `sysctl` information arose under the restricted runtime. They were environment notices, not logged model/KKT failures.

### 13.4 Grid boundaries and baseline reuse

The allowed bounded grid expansion was exercised. Two models remained at a boundary after the maximum two rounds:

- PROBE outer/2, M2: `lambda≈1e-6, rho=0.5`.
- PROBE→GENEPARK-TRO transfer, M5: `lambda≈1e-6, rho=1.0`.

These boundary results are retained and do not authorize further tuning. In particular, the small transfer penalty and poor probability performance may motivate discussion, but do not prove the cause of the transfer result or permit a replacement model chosen from the target outcomes.

M3 reused M1 in PROBE outer/0, outer/1 and outer/2. M5 reused M4 in PROBE outer/1 and outer/3. Such baseline branches are permitted and should not be represented as independently fitted identical models. The final pooled M5 selected a nonzero gene block and is distinct from final M4.

## 14. Actual development results

### 14.1 Five-fold nested-CV summaries

Each value below is the equal mean of the five outer-fold metrics within the named source; AUROC/AUPRC are accompanied by fold SD. **SD is not a confidence interval**, and outer-fold metrics are not independent replications. Sources remain separate; cross-source pooled AUROC is not an endpoint.

| Source | Model | AUROC mean ± fold SD | AUPRC mean ± fold SD | Mean Brier | Mean log loss |
|---|---|---:|---:|---:|---:|
| GSE57475 | M1 | 0.4749 ± 0.0737 | 0.6396 ± 0.0483 | 0.2650 | 0.7263 |
| GSE57475 | M2 | 0.5564 ± 0.0738 | 0.7064 ± 0.0845 | 0.2788 | 1.1634 |
| GSE57475 | M3 | 0.5371 ± 0.1410 | 0.6939 ± 0.0964 | 0.2709 | 0.9645 |
| GSE57475 | M4 | 0.5615 ± 0.1701 | 0.7610 ± 0.0860 | 0.2752 | 0.7755 |
| GSE57475 | M5 | 0.5911 ± 0.1234 | 0.7622 ± 0.0748 | 0.2931 | 1.0392 |
| GSE99039 | M1 | 0.5835 ± 0.0751 | 0.6123 ± 0.0960 | 0.2489 | 0.6923 |
| GSE99039 | M2 | 0.7042 ± 0.0931 | 0.7224 ± 0.0715 | 0.2423 | 0.7483 |
| GSE99039 | M3 | 0.6704 ± 0.1204 | 0.6868 ± 0.0826 | 0.2759 | 0.9952 |
| GSE99039 | M4 | 0.6630 ± 0.1266 | 0.6897 ± 0.0979 | 0.2404 | 0.6827 |
| GSE99039 | M5 | 0.6861 ± 0.1052 | 0.7070 ± 0.1001 | 0.2509 | 0.7527 |

Mean within-source M5−M4 AUROC is **+0.0296** in PROBE and **+0.0231** in GENEPARK-TRO. No confidence interval or significance claim for these descriptive differences has been established by the current aggregate summary. AUPRC should be interpreted against case fractions of about 65.5% and 50.2%, not compared across populations as if prevalence were identical.

### 14.2 Outer-fold AUROC detail

| Source / outer fold | Test n (PD) | M1 | M2 | M3 | M4 | M5 | M5−M4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE57475/outer/0 | 29 (19) | 0.3658 | 0.5105 | 0.3658 | 0.3947 | 0.4895 | +0.0947 |
| GSE57475/outer/1 | 29 (19) | 0.5368 | 0.5000 | 0.5368 | 0.7000 | 0.7000 | +0.0000 |
| GSE57475/outer/2 | 28 (19) | 0.4327 | 0.6550 | 0.4327 | 0.3626 | 0.4327 | +0.0702 |
| GSE57475/outer/3 | 28 (18) | 0.5083 | 0.5000 | 0.7000 | 0.6333 | 0.6333 | +0.0000 |
| GSE57475/outer/4 | 28 (18) | 0.5306 | 0.6167 | 0.6500 | 0.7167 | 0.7000 | -0.0167 |
| GSE99039/outer/0 | 48 (24) | 0.4792 | 0.5521 | 0.4826 | 0.4514 | 0.5764 | +0.1250 |
| GSE99039/outer/1 | 48 (24) | 0.5347 | 0.6979 | 0.7049 | 0.7118 | 0.7222 | +0.0104 |
| GSE99039/outer/2 | 48 (24) | 0.6649 | 0.7465 | 0.7500 | 0.7448 | 0.7135 | -0.0312 |
| GSE99039/outer/3 | 48 (24) | 0.6137 | 0.7240 | 0.6285 | 0.6441 | 0.5885 | -0.0556 |
| GSE99039/outer/4 | 47 (24) | 0.6250 | 0.8007 | 0.7862 | 0.7627 | 0.8297 | +0.0670 |

The small PROBE test folds contain 28–29 records and the GENEPARK-TRO folds 47–48. Substantial fold-to-fold variation is visible. The first pilot's result was preserved and was not used to change the method.

### 14.3 Directional transfer: all five models

| Training → test | Test n (PD) | Model | AUROC | AUPRC | Brier | Log loss |
|---|---:|---|---:|---:|---:|---:|
| GSE57475 → GSE99039 | 239 (120) | M1 | 0.6218 | 0.6288 | 0.2447 | 0.6834 |
| GSE57475 → GSE99039 | 239 (120) | M2 | 0.4604 | 0.4817 | 0.3926 | 1.2858 |
| GSE57475 → GSE99039 | 239 (120) | M3 | 0.5239 | 0.5451 | 0.3732 | 1.7854 |
| GSE57475 → GSE99039 | 239 (120) | M4 | 0.6665 | 0.6277 | 0.2341 | 0.6644 |
| GSE57475 → GSE99039 | 239 (120) | M5 | 0.4882 | 0.4849 | 0.5063 | 6.4140 |
| GSE99039 → GSE57475 | 142 (93) | M1 | 0.5854 | 0.7005 | 0.2364 | 0.6676 |
| GSE99039 → GSE57475 | 142 (93) | M2 | 0.6634 | 0.7905 | 0.3123 | 0.9655 |
| GSE99039 → GSE57475 | 142 (93) | M3 | 0.6669 | 0.7895 | 0.3089 | 0.9558 |
| GSE99039 → GSE57475 | 142 (93) | M4 | 0.6351 | 0.7537 | 0.3000 | 0.8216 |
| GSE99039 → GSE57475 | 142 (93) | M5 | 0.6265 | 0.7612 | 0.3686 | 1.2969 |

For PROBE→GENEPARK-TRO, M4 AUROC is **0.6665** and M5 is **0.4882**, a difference of **−0.1784**. For GENEPARK-TRO→PROBE, M4 is **0.6351** and M5 is **0.6265**, a difference of **−0.0086**. Thus M5 did not improve over M4 in either direction at the observed point estimate.

The first-direction M5 Brier score is **0.5063** and log loss is **6.4140**, compared with M4 values of **0.2341** and **0.6644**. This is a serious probability-performance problem in that transfer. AUROC alone would not describe it. No target-data recalibration has been performed to repair the result.

These are two distinct directions of transfer, not a single pooled LOCO estimate. They do not isolate platform, medication, cell biology, source preprocessing, selection or demographic distribution as a unique cause.

### 14.4 Final pooled fitting

| Model | Selected lambda | Rho | Gene penalty | Nonzero gene coefficients | Inner selection AUROC |
|---|---:|---:|---|---:|---:|
| M1 | — | — | None | 0 | 0.5883 |
| M2 | 0.1 | 0 | L2 / ridge | 1000 | 0.6758 |
| M3 | 0.01 | 1 | L1 / lasso | 132 | 0.6695 |
| M4 | — | — | None | 0 | 0.6368 |
| M5 | 0.001 | 1 | L1 / lasso | 210 | 0.6739 |

M1/M4 have no tuned penalty. M2 is a ridge gene model, while M3/M5 selected lasso gene blocks. The final baseline remains unpenalized. The tuning scores used to make these selections are selection criteria on development data, not unbiased estimates of final external performance.

The 210 nonzero M5 coefficients are a property of this fitted procedure. No gene-level significance, stable biomarker discovery or pathway enrichment claim follows from this count. All final fitted preprocessing and coefficient arrays must be included in the required pre-PPMI commitment; fitting them is not the same as completing that public lock.

## 15. What the results currently support

The completed analyses show modest average internal improvements for the gene-augmented model and no observed improvement over the marker baseline in the two directional transfers. The deterioration is asymmetric and especially large for PROBE-trained M5. This supports studying transportability and calibration rather than claiming that strong within-source discrimination necessarily transfers.

They do **not** yet support:

- A result for PPMI or a clinical diagnostic claim.
- A causal estimate of how much signal comes from disease versus shortcuts.
- Proof of no PD-related biological information in blood.
- Proof that within-sample ranks transfer all relevant genes faithfully between assay types.
- Equivalence, absence of a 0.05 benefit, or statistical significance based only on point estimates or fold SD.
- Selection of a replacement primary model, transformation or threshold after observing the transfer results.
- Guaranteed Q1/Q2 acceptance.

The current potential paper is a methodological transportability study with a strong baseline and an auditable procedure. Its final strength depends on the remaining controls, uncertainty analysis, external eligibility, interpretation and limitations.

## 16. Planned evaluation metrics and inference

The sole confirmatory comparison remains **PPMI AUROC(M5)−AUROC(M4)**, conditional on the registered gates. Primary uncertainty is a **95% percentile paired donor-bootstrap CI using 5,000 class-stratified resamples**, with common resampling indices for both fixed developed models. This covers test-sample variability with models held fixed; it does not cover all model-training uncertainty.

AUROC, AUPRC, sample size, case fraction, Brier score and log loss are reported by population. Model probabilities from these case/control and class-weighted analyses are not claimed to be calibrated clinical population risk. Balanced accuracy at 0.5 is optional/descriptive for PD models, not a confirmatory decision rule. DeLong is not an alternate primary method to choose after seeing results; a single primary comparison does not require a primary Holm family correction.

The research-relevance margin is **delta*=0.05 AUROC**, not an established clinical-benefit threshold. For distinct M4/M5 models, an upper CI bound below 0.05 excludes an increase at least that large under this fixed procedure. A CI spanning both zero and 0.05 is inconclusive. A positive lower bound would support a positive increment, while uncertainty and applicability still matter.

If development had selected M5≡M4, the exact zero difference/CI would have been tautological. The registered interpretation would be "genes not selected in development," not proof that possible gene value is zero. Final pooled M5 is distinct, so that final branch was not selected here.

Earlier precision scenarios showed that n≈472 may not exclude a 0.05 increment in some model-score-correlation settings. A roughly zero-centred interval with half-width 0.055 would not exclude that margin. Such planning scenarios are not observed PPMI performance, actual bootstrap-coverage validation, guaranteed power or a reason to narrow the interval post hoc.

## 17. Remaining controls, diagnostics and their consequences

### 17.1 Sex-transfer control

Train ridge logistic regression only in healthy GEO controls using the five fixed Y markers `KDM5D, DDX3Y, RPS4Y1, UTY, EIF1AY`. XIST is absent from the current common universe. Use the same checked backend, lambda grid, rho=0, inner-AUROC tuning and bounded boundary rule. Both GEO directions must reach **AUROC≥0.95** before the external gate has demonstrated development capacity.

The PPMI healthy-control gate requires **AUROC≥0.95 and 95% CI lower bound≥0.90**, with 5,000 sex-stratified donor-bootstrap repetitions. Balanced accuracy at 0.5 is descriptive. Passing sex concordance cannot prove correct identity among donors of the same sex, and predicted sex must not overwrite recorded labels or automatically exclude participants.

### 17.2 Age-transfer control

Train an autosomal-gene ridge regressor in healthy GEO controls only, selecting up to 1,000 genes by the registered rules. Alpha grid: `0.01, 0.1, 1, 10, 100`; tuning metric MAE; ties prefer larger alpha. Report Spearman correlation and MAE against a constant predictor using the training-median age.

Development capacity requires **Spearman≥0.4 in both GEO directions and better MAE than the training-median predictor**. Limited age capacity does not automatically halt the study, but prevents treating an external age failure as evidence that PPMI itself is defective. This control tests one autosomal relationship, not all potential PD-related autosomal signal.

### 17.3 CBC validation and longitudinal identity diagnostic

The guide identifies a Blood Chemistry/Hematology table, but its required file/linkage has not been completed in this project. Prefer same-day blood measurements; otherwise the registered window is nearest within **±30 days**. Report units, timing gaps, missingness and sample counts.

The primary panel/count association is within healthy controls; PD associations are separate exploratory summaries. For a corresponding measurement, supportive descriptive evidence requires **n≥30, Spearman≥0.4 and CI lower bound>0**. Total lymphocytes do not specifically validate T/B scores; RBC count does not substitute for reticulocytes. Weak or absent CBC evidence limits the claim and leaves the fixed models unchanged. It is not a stop gate.

Longitudinal within-PATNO versus other-PATNO similarity, matched on recorded sex and wave where feasible, is a diagnostic only. Its separate technical centering does not enter the inductive PD model. It cannot automatically exclude/reassign samples or replace genetic identity matching. It has not been executed.

### 17.4 Bounded sensitivities and negative controls

The registered list includes X/Y-gene removal, training-only age/sex residualization, age/sex matching with a **three-year caliper** and no replacement, exclusion of the six panels' genes from the added gene block, development technical-metadata diagnostics and the development ComBat audit. Matched PPMI evaluation uses fixed models and a matched-pair bootstrap unit.

ComBat is not the primary pipeline. Target-batch transductive correction, correction for known training batches and a genuinely new batch are different deployment scenarios. Test disease labels must not enter batch estimation. A diagnostic using all X is to be disclosed as transductive, not presented as inductive performance.

The registered secondary negative control still specifies **1,000 development permutations**, within cohort×fixed-outer-fold blocks, repeating every label-dependent stage. It has not been silently reduced or executed. Its scope/cost must be made explicit before execution; a change requires a timestamped amendment. Permutations do not create more participants or remedy low precision. The earlier pilot extrapolations suggested a potentially multi-day burden and were scenarios, not validated scheduling promises.

## 18. Remaining execution roadmap and deliverables

| Order | Work package | Concrete deliverable | Decision consequence |
|---|---|---|---|
| 1 | Prediction-warning review | Trace clipping to the relevant scoring calls; assess whether selected-model/tuning behavior is affected; preserve current outputs and write a technical finding. | No silent pipeline change. A material bug requires an amendment before affected re-execution. |
| 2 | GEO auxiliary controls | Implement and synthetically verify sex/age routines; execute both healthy-control transfer directions and save fit/metric evidence. | Failure of the sex development prerequisite prevents conditional confirmatory PPMI execution; age has a reporting-limitation role. |
| 3 | Real execution membership/leakage audit | Compare every real fit's membership/preprocessing fingerprint with the frozen split manifest and verify prediction membership. | A material failure invalidates the affected analysis; do not drop problematic folds. |
| 4 | Bootstrap coverage validation | Fixed synthetic scenarios and a coverage report for the registered interval procedure, including ties and distinct/identical model branches. | Material failure requires a public amendment; no CI-method switching based on PPMI outcomes. |
| 5 | Registered sensitivities/permutation plan | Bounded implementations, computational plan, and any necessary pre-execution amendment; then complete the authorized secondary runs. | Do not silently remove stages or introduce outcome-driven sensitivities. |
| 6 | PPMI final eligibility | Validated one-to-one file/release/visit/donor linkage; resolved age/diagnosis/QC rules; exclusion flow; treatment and publication-condition audit. | Unresolved unique-link ambiguity stops confirmatory interpretation. |
| 7 | Second public lock | Hash manifest of code/environment, splits, all fitted M1–M5/preprocessing states, controls, and a non-identifying hash/reference for the private final population/crosswalk. | Must be public and verified before PPMI scores in this workflow; the GEO milestone is not this lock. |
| 8 | PPMI outcome-free controls | Execute fixed sex/age/CBC/longitudinal diagnostics under their registered roles and report all outcomes. | Apply the explicit sex-gate and technical-bug exit rule; no pipeline rescue on the same holdout. |
| 9 | Conditional PPMI disease scoring | One fixed inductive application of saved M1–M5, aggregate metrics and paired uncertainty, plus the registered bounded evaluations. | If conditions are unmet, report confirmatory outcome as unevaluable rather than invent a replacement endpoint. |
| 10 | Manuscript and reproducibility package | Methods/results/limitations, population flow, model/feature documentation, figures, aggregate data, code and appropriate acknowledgements. | Claims remain bounded to the evaluated procedure and populations. |

The immediate next analytical work is warning review and the development sex/age transfer controls. Another general feasibility cycle is not required to recognize that the core GEO execution has already occurred. Conversely, the existence of a working engine and published GEO results does not make the PPMI analysis ready today.

## 19. Stop rules and unresolved risks

The following remain explicit risks rather than established causes of the observed results:

- PPMI count/participant identity and final sample eligibility; a unique computational link is necessary but not biological identity proof.
- Limited number of independent development sources and the selected TRO complete-case population.
- Operational GSM donor IDs in PROBE and incomplete cross-source overlap evidence.
- Different inherited source preprocessing, assay dynamic range, RNA-seq ties/zeros and count-length transformation.
- Demographic imbalance, medication at sampling and recruitment/processing differences.
- Imperfect cell-marker proxies and unavailable or incomplete CBC validation.
- Persistent hyperparameter boundaries and probability overconfidence in some transfers.
- Finite external sample size and potentially inconclusive incremental-AUROC uncertainty.
- Warning attribution outside fit logs and uncompleted full real-data leakage/coverage audits.
- The substantial cost of full retuned permutations and sensitivities.

If the sex gate or linkage/transformation condition fails, continuation is permitted only for an independently demonstrable technical bug, with a timestamped amendment and repetition of the same control while preserving all runs. Without such evidence, the conditional confirmatory result is **unevaluable under the registered protocol**. A different transformation or model chosen after a nontechnical failure requires a new study/registration and independent validation.

## 20. Repository, public records and preservation

GitHub is the living work record. Important progress was committed, including input preflight (`de5d7f2`), environment blocker (`fc02fc0`), successful plan (`30042d4`), completed core GEO results (`c80d3e8`) and the Zenodo publication receipt (`d9e8fe9`). The full Git history preserves the exact changes. This report adds documentation and an explicit descriptive correction; it does not edit the frozen analysis configuration or models.

Zenodo is used for registration and major milestones rather than every incremental status update. All versions remain under [concept DOI 10.5281/zenodo.23405370](https://doi.org/10.5281/zenodo.23405370). The latest published milestone contains four files: English results report, aggregate JSON, SHA-256 release manifest and source ZIP at commit `c80d3e84b2ea07e5aef8ca7b63b9c01268d6c728`. All uploaded MD5 checksums matched local files. It is a post-results snapshot, not a new preregistration or final external lock.

| Commitment | SHA-256 |
|---|---|
| Locked `pipeline/geo_pipeline.py` | `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c` |
| Frozen split manifest | `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2` |
| Published source snapshot ZIP | `bcc039668b1c4124b28c7b4ed66d11d2b4d1404d3cb70ac9cfbc907486851fba` |

A byte-for-byte verified local backup of **167 post-pilot output files** is preserved under ignored `private_data/geo_runs_20261018/`. It includes scope-level start/results/decision records, fit logs, preprocessing states, model arrays and held-out predictions. The original temporary execution directories also remain. Local artifacts are not public uploads. The earlier pilot's aggregate and hash evidence are in `results/geo_pilot/`; full pilot artifacts belong to the earlier analysis workspace and their continuing availability should be checked before a complete reproducibility handoff.

### 20.1 Artifact guide

| Path relative to repository | Contents / role |
|---|---|
| `docs/COMPLETE_PROJECT_REPORT_20261018_EN.md` | This complete current account. |
| `docs/STATUS.md` | Concise current status and next gates. |
| `docs/PROJECT_HISTORY_AND_ROADMAP.md` | Earlier living history and roadmap. |
| `docs/REPORTING_CORRECTION_20261018_EN.md` | Lasso/ridge wording correction, preserving the published snapshot. |
| `docs/GEO_CV_RESULTS_20261018_EN.md` | Original short GEO milestone report as released, with the separately disclosed terminology error. |
| `results/geo_development_20261018/aggregate_summary.json` | Exact fold, transfer, selected-model and numerical summaries. |
| `results/geo_development_20261018/ZENODO_PUBLICATION_RECEIPT.json` | Verified DOI/version, source commit and upload status. |
| `results/geo_development_20261018/GEO_CV_RELEASE_MANIFEST_20261018.json` | Exact public file hashes and source commit. |
| `results/project_history_20261018/report_evidence_manifest.json` | Evidence files/hashes used for this account and aggregate feasibility provenance. |
| `pipeline/geo_pipeline.py` | Primary GEO loading, preprocessing, split execution, fitting and prediction. |
| `pipeline/paired_bootstrap.py` | Fixed-score paired donor resampling utilities. |
| `pipeline/verify_pipeline.py` and synthetic report | Implementation invariants and 19-check evidence. |
| `protocol/PROTOCOL_v3_1_EN.md`, `PREREG_v3_1_FULL_EN.md` | English reading copies with original-source/timing caveats. |
| `protocol/config_v3_1.json` | Frozen design, grids and gates; freeze-time status booleans are historical. |
| `protocol/numerical_kkt_v3_1.py` | Solver dispatch, objective/KKT checks and fixed retry. |
| `protocol/*verification.json`, stress/null JSON and EN reports | Numerical/synthetic evidence. |
| `protocol/TRANSLATION_*`, `PROVENANCE_CORRECTION*` | Translation/source provenance and corrections. |
| `data/manifests/` | Gene mapping, fixed gene list/lengths, input hashes and realized splits. |
| `prepared/` | Ignored local prepared GEO inputs. |
| `private_data/geo_runs_20261018/` | Ignored, persistent local analytical outputs. |

## 21. Final status statement

**Feasibility was performed first. The project has subsequently been defined, registered, coded, numerically tested, and executed on the two GEO development sources. The core GEO nested CV, two transfers and final development models now exist.** The project is neither only an idea nor a completed external-validation paper. The next work is the bounded diagnostic/control sequence, final external eligibility, the second lock and conditional PPMI evaluation described above.

