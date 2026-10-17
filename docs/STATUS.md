# Project status

**Snapshot date:** 2026-10-16 (Asia/Tehran)

| Workstream | Status | Evidence / limitation |
|---|---|---|
| Research question and feasibility | Completed for the current scope | Two GEO development sources and an approved PPMI data-access path were identified. This does not imply that external validation has been executed. |
| Gene mapping and cross-platform input design | Completed as a preparation audit | 16,364 shared protein-coding genes were retained; 15,720 are autosomal. Mapping coverage does not prove that each gene is measurable in every PPMI count table. |
| PPMI IR3 archive review | Partial | Remote archive and aggregate audits were reviewed. File-to-release-row matching was reported, but individual biological identity, final eligibility, and the final locked external population remain unresolved. No PPMI scores were computed. |
| Statistical protocol | Publicly registered and versioned | See the linked Zenodo lineage and English provenance notes in `protocol/`. The English documents are translations, not replacements for the original source record. |
| Solver and synthetic checks | Completed at this stage | Weighted/block-penalty and KKT checks, stress runs, and 19 pipeline synthetic checks are documented in the protocol and pilot report. Bootstrap CI coverage validation remains outstanding. |
| GEO pipeline | Implemented; execution incomplete | One prespecified GSE57475 outer fold completed. No complete nested CV or directional transfer result exists yet. |
| PPMI external scoring | Not started | A second public lock of code/environment, splits, final models/preprocessing, and the finalized private population must precede all PPMI scores. |
| Manuscript | Not started | It is premature to claim external validation, quantify shortcut contribution, or infer clinical utility. |

## GEO input-hash preflight (2026-10-16)

The four prepared GEO inputs were found in the prior local working directory and their SHA-256 values exactly match `data/manifests/input_audit.json`. The checked-in split-manifest SHA-256 also matches the value in that input audit. The input files are exposed to this checkout through ignored local symlinks under `prepared/`; no expression or sample-level data were added to Git.

The no-fit pipeline `plan` command was run successfully by the user on macOS using an alternate package index after the default PyPI file host timed out. It validated the four local GEO inputs and produced `n=381`, `genes=16,364`, `GSE57475=142`, and `GSE99039=239`. All four input SHA-256 values matched the frozen audit, and the generated split-manifest SHA-256 matched `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. The output reports `real_training_performed=false`; no model fit was run. The repository still does not contain the required local receipt proving that the exact pipeline-code and split hashes are covered by the public registration. Do not start additional real-data fits until that receipt is verified.

## First real-data pilot (not a study result)

The single completed fold was `GSE57475/outer/0`, with 113 training and 29 held-out records. Five-model test AUROCs were M1 0.3658, M2 0.5105, M3 0.3658, M4 0.3947, and M5 0.4895. The fold took 427.34 seconds wall time and completed 269 fit attempts with no warning or KKT failure; the maximum KKT residual was 3.8924e-8, below the fixed 1e-6 threshold. M3 reused M1 in this fold. These values are retained as an execution pilot, not as final performance estimates, a separate test, or confirmatory evidence.

## Immediate next step

Verify the public registration receipt against the exact pipeline-code and split hashes. Only after this check passes, run the remaining prespecified GEO analyses; estimate compute time from the first full planned scope. Preserve every fit log and aggregate result locally; commit only non-sensitive reproducibility evidence and aggregate outputs.
