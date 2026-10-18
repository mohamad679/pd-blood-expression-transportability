# Project Continuation Report: Work Since the Full Retrospective Report

**Author:** Mohammad Javad Asgari  
**ORCID:** [0009-0002-2313-571X](https://orcid.org/0009-0002-2313-571X)  
**Report date:** 18 October 2026 (Asia/Tehran)  
**Repository:** `mohamad679/pd-blood-expression-transportability`  
**Scope:** Work completed after the baseline in [`COMPLETE_PROJECT_REPORT_20261018_EN.md`](COMPLETE_PROJECT_REPORT_20261018_EN.md), through the archived whole-blood collection linkage audit.

This is a continuation of the previous retrospective report, not a replacement preregistration, amendment, or claim of completed external validation. It records the work, results, boundaries, and remaining gates since that report. The disease pipeline and GEO results were not changed. No PPMI disease score, marker score, or marker association was calculated during this interval.

## 1. Starting point

At the previous report's cutoff, the main GEO nested cross-validation and directional transfer analyses were complete and had been released as Zenodo version 3.1.3-geo-cv. The six sigmoid warnings still required review; GEO sex/age transfer controls had not yet run; CBC availability and temporal linkage were incomplete; PPMI identity and eligibility remained provisional; the permutation proposal was not executed; and no second public lock existed for PPMI scoring.

The work below addressed several of those diagnostics and data-linkage questions. It did not turn the provisional PPMI cohort into a locked test population.

## 2. Fixed GEO control plan and synthetic checks

Before running the real auxiliary controls, the implementation plan was committed in Git (`36db979`). It used healthy controls only, both directional GEO transfers, five donor folds, training-only preprocessing and tuning, and fixed outcome-specific criteria. No PPMI input was used.

Eight synthetic checks passed before the real controls: held-out mutation invariance; training-only standardization; donor separation; use of the intended sex target rather than PD labels; recoverable synthetic sex and age signals; fixed KKT acceptance; detection of a sigmoid-induced mixed-class tail tie; and paired bootstrap behavior. The diagnostic plans, split membership and detailed sample-level outputs remain local/private; only aggregate evidence and hashes are in the repository.

## 3. Sigmoid-warning investigation

The warnings arose because the pinned `glum` logit prediction path clips probabilities near zero and one. Since clipping can tie cases and controls, the warnings were investigated rather than assumed harmless.

The first broad replay driver stopped before fitting because it passed a duplicate context argument. The driver was corrected. A broader diagnostic was then stopped after 111 completed fits to bound the investigation to the six contexts that had actually reported warnings. These attempts did not alter the original scientific runs. The scope and driver history are documented in Git.

The completed targeted replay covered all candidates and five inner folds in the six reported contexts, plus their required baseline candidates: **490 fits**. Findings:

- All six warning events were reproduced. Each involved one upper-tail probability; none merged a differently labelled case/control pair at the clipped tail.
- Clipped probabilities, unclipped `expit` probabilities, and linear scores produced the same AUROC for every audited validation call.
- Every recorded candidate mean was reproduced exactly; maximum absolute difference was **0.0**. No hyperparameter choice changed, and none of the warning-producing candidates had been selected originally.
- All 490 replay fits passed the fixed KKT criterion on the first attempt. Maximum KKT residual was **2.2578e-8**, below the `1e-6` threshold. There were no replay retries.
- Wall-clock replay time was **2,584.24 seconds (43.07 minutes)**.

The audit therefore closed without changing the scoring method or rerunning the core GEO science. It is bounded to the six observed warning contexts; it does not claim to have replayed every candidate in every GEO scope. The locked disease pipeline hash and original GEO artifacts remain unchanged. See [`SIGMOID_AUDIT_20261018_EN.md`](SIGMOID_AUDIT_20261018_EN.md) and [`SIGMOID_REPLAY_SCOPE_20261018_EN.md`](SIGMOID_REPLAY_SCOPE_20261018_EN.md).

## 4. GEO sex and age transfer controls

The controls used only healthy controls: 49 PROBE participants and 119 GENEPARK-TRO participants. Both directions were evaluated; no PD case entered fitting or evaluation. The four frozen GEO inputs matched the previously recorded hashes.

| Target | PROBE → GENEPARK-TRO | GENEPARK-TRO → PROBE | Predefined development criterion | Outcome |
|---|---:|---:|---|---|
| Sex AUROC | 1.000000 (95% bootstrap CI 1.000000–1.000000) | 0.993311 (0.974875–1.000000) | AUROC ≥ 0.95 in both directions | Passed |
| Age Spearman | -0.058118 (-0.238535–0.121006) | 0.258599 (-0.051083–0.522189) | Spearman ≥ 0.4 and MAE better than median baseline in both directions | Failed to demonstrate capacity |
| Age MAE, years | 7.038 vs 7.025 median baseline | 17.482 vs 7.694 median baseline | Same rule as above | Failed |

The sex model used the five preregistered Y markers. The selected penalty was lambda=100, the upper edge after the two allowed extensions, in both directions. Balanced accuracy at threshold 0.5 was 0.5 in both directions; it was descriptive only and did not invalidate the AUROC-based gate. Across 72 checked sex fits, all passed KKT, the largest residual was `2.6747e-8`, and no fit or prediction warnings or retries occurred. The success validates this limited sex/Y-marker signal across the two GEO sources; it does not certify autosomal, cell-composition, or PD-signal transport.

For age, the fixed alpha grid selected 100, its upper endpoint, in both directions. The observed failure therefore describes this fixed age predictor; the result cannot separate platform shift from the estimator's fixed feature selection, tuning range, and training sample. No grid extension or rescue transformation was performed. The age result does not establish failure of all autosomal signals and is not a PPMI-specific result.

The final GEO-trained sex model for PPMI has **not** been fit or evaluated. The PPMI sex gate remains for the post-lock external stage. Full methods and outputs are in [`GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md`](GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md) and [`GEO_TRANSFER_CONTROLS_20261018_EN.md`](GEO_TRANSFER_CONTROLS_20261018_EN.md).

## 5. Permutation scope and compute planning

The originally stated 1,000 permutations with full retuning across all five models would be prohibitively expensive using measured runtimes. A bounded proposal was recorded: **200 permutations**, models **M4 and M5**, and the ten within-source outer folds only. The proposal retains within-cohort/class constraints and reruns all label-dependent operations; it excludes directional transfers and final pooled fitting. The reported p-value form is `(1 + count at least as extreme) / (B + 1)`.

The measured M4/M5 runtime was 18.38 minutes per original pass, corresponding to an estimated **61.27-hour serial baseline for B=200**, before overhead or null-run runtime changes. This is a planning estimate, not an observed permutation runtime. The amendment is still a draft, is **not published or effective**, and **no registered real-GEO permutation has been run**. The executor, seed/membership/resume behavior, and worker/memory limits still need preflight before publication and execution. See [`STATUS.md`](STATUS.md) and [`../protocol/PERMUTATION_AMENDMENT_DRAFT_20261018_EN.md`](../protocol/PERMUTATION_AMENDMENT_DRAFT_20261018_EN.md).

## 6. CBC availability and version concordance

After IDA access, the user downloaded the current and archived `Blood Chemistry & Hematology` exports. The aggregate-only audit was matched exactly on canonical participant ID and visit against the 2021 Long RNA-seq metadata:

- The RNA metadata had **4,738** Long RNA-seq rows; **106** did not have a parseable visit key, leaving **4,632** valid rows/unique participant-visit pairs.
- Of those RNA sample-visits, **1,984** had a record in the current CBC export and **2,060** in the archived export; **76** were archived-only and none current-only.
- Among numeric records for the nine prespecified core CBC markers, **56,104** shared exact measurement keys across the two CBC versions had identical values and units. There were **zero** differing-value/unit keys and **zero** duplicate-ambiguous shared keys.

This establishes substantial visit-level measurement availability and concordance for the shared marker rows. It does not establish that the complete CBC files are identical, validate the cell-marker panels, or demonstrate disease transfer. The raw CBC files and row-level joins remain outside Git; only aggregate outputs are stored under ignored `private_data/`.

## 7. Archived whole-blood collection linkage and date precision

The user then supplied the archived `Whole_Blood_Sample_Collection` export. A new script reuses the established RNA sample parser and performs exact normalized `PATNO + visit/EVENT_ID` matching. It writes aggregate counts only; it does not read outcomes or expression values.

- The archive contained **13,793** rows and **13,793** valid unique participant-visit keys, with no duplicated keys. It exactly overlapped **4,521** of the **4,632** valid RNA sample-visits.
- The field `BLDCLLDT` was month/year precision for **13,235** archive rows and blank for **558**. Among the 4,521 RNA-overlap visits, **4,316** had a month/year value.
- `BLDCOLL` values among those exact overlaps were: code `0`, **206**; code `1`, **4,313**; code `2`, **2**. Their meanings have not been verified against the applicable PPMI codebook, so the report deliberately does not label them “collected” or “not collected.”
- CBC `LCOLLDT` was also month/year precision. Among exact RNA/collection overlaps with numeric core-marker CBC data, **1,824** current-table visits and **1,899** archived-table visits had the same reported month as `BLDCLLDT`.

Same reported month is useful as a candidate temporal subset: if both fields represent the collection dates implied by their names, calendar dates in the same month cannot be more than 30 days apart. But month/year precision cannot establish exact day differences, identify a same-day match, or reliably select the nearest CBC when records fall in adjacent months. Therefore the complete registered “same day or nearest within 30 days” criterion has **not** been verified, and no temporal eligibility gate is marked passed. The archive improves the linkage evidence but does not remove this date-resolution limitation.

The audit did not compute CBC-derived cell-marker scores, correlations, any PPMI outcome metric, or any PPMI disease prediction. See [`PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md`](PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md).

## 8. Repository and publication record

The continuing work was recorded in Git, not as new Zenodo releases. Relevant commits after the prior full report include:

| Commit | Work recorded |
|---|---|
| `36db979` | Frozen GEO sex/age control implementation and sigmoid audit plan before execution. |
| `eb7b9c0` | Corrected a retrospective replay-driver argument error before any fit in that attempt. |
| `153ae18`, `c17ba61` | Bounded and documented the replay scope and preserved pre-control plan provenance. |
| `fb4b5b6` | Recorded the GEO sex/age control results. |
| `532d90a` | Recorded the completed targeted sigmoid audit and unchanged core results. |
| `ea65937` | Recorded the still-unpublished permutation proposal and compute estimate. |
| `b29dfe7`, `9512320` | Added the aggregate CBC coverage audit and current/archive marker-row concordance. |
| `a5c7be3`, `ca91856` | Added the aggregate whole-blood collection linkage audit, its code, and the final overlap count. |

The latest public Zenodo milestone remains version 3.1.3-geo-cv, DOI [10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851). No Zenodo amendment has been published for the permutation proposal, and the second pre-PPMI lock has not been made. No PPMI scoring has occurred.

## 9. Current project state and work remaining

### Completed in this interval

1. Reproduced and closed the six-warning sigmoid diagnostic without changing AUROCs, selections, the locked model, or the core GEO results.
2. Ran the prespecified GEO healthy-control sex and age transfer controls; the sex gate passed and the fixed age predictor did not meet its capacity criterion.
3. Received the current and archived CBC tables, measured exact visit-level marker availability, and confirmed shared numeric marker-row concordance.
4. Received the archived whole-blood collection export, linked it to RNA visits exactly, quantified its date precision, and recorded the remaining temporal limitation.
5. Added English aggregate reports and reproducible local audit code to Git; kept all restricted rows and person-level outputs outside Git.

### Still required before any PPMI disease score

1. Verify the `BLDCOLL` value mapping from the correct PPMI dictionary/codebook and document which date fields correspond to which collection events.
2. Finish the bounded GEO leakage/membership checks, bootstrap-coverage simulation, and prespecified sensitivity analyses. No additional controls should be added merely to rescue an observed result.
3. Complete the tested, documented permutation executor and amendment; publish the amendment before any permutation run, then execute only that approved scope.
4. Finalize count-file-to-person/visit linkage, eligibility, diagnosis/covariate reconciliation, QC exclusions by group, treatment-at-sampling reporting, and the non-identifying population reference/hash.
5. Fit and preserve the exact final GEO-trained sex control and all final M1–M5 preprocessing/model artifacts needed for transfer.
6. Publish the second public lock of the code, environment, realized splits, final models/preprocessing, and non-identifying population reference **before** the first PPMI disease score. The registered PPMI controls and scoring can then proceed without using PPMI outcomes to alter the pipeline.
7. Evaluate external results with the prespecified uncertainty and failure rules. Only then finalize the manuscript and publication package.

The project has progressed beyond feasibility and has a completed GEO development study, but the external-validation phase is still pending. The current evidence supports a careful methodological paper in principle; it does not yet establish that the PPMI gene increment will be positive or precise, nor does it guarantee a Q1/Q2 outcome.

## 10. Reproducible files

- Previous full history and roadmap: [`COMPLETE_PROJECT_REPORT_20261018_EN.md`](COMPLETE_PROJECT_REPORT_20261018_EN.md)
- Current status: [`STATUS.md`](STATUS.md)
- Sigmoid audit: [`SIGMOID_AUDIT_20261018_EN.md`](SIGMOID_AUDIT_20261018_EN.md)
- GEO control plan and results: [`GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md`](GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md), [`GEO_TRANSFER_CONTROLS_20261018_EN.md`](GEO_TRANSFER_CONTROLS_20261018_EN.md)
- CBC collection linkage: [`PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md`](PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md)
- Aggregate-only linkage audit implementation: [`pipeline/audit_ppmi_collection_linkage.py`](../pipeline/audit_ppmi_collection_linkage.py)
- Restricted-input aggregate output, intentionally ignored by Git: `private_data/ppmi_collection_crosswalk_20261018.json`
