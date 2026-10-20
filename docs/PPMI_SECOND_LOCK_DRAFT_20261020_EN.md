# Draft second pre-PPMI lock — not effective

**Status:** planning document only. No PPMI score is authorized by this draft. A completed, timestamped Zenodo version with a DOI and the matching exact Git commit/tag is required before the first PPMI cell-panel or model/control score.

## Current execution site

- The researcher's personal MacBook is the sole site for all **remaining** PPMI acquisition, local storage, eligibility and metadata checks, control/disease scoring, and inference.
- The existing provisional 472-candidate local bundle may be used for checks supported by its frozen contents. Obtain missing metadata directly into local, nonsynced Mac storage. Do not run further PPMI processing on another computer or server.
- Keep participant-level PPMI inputs and outputs under `~/PPMI_private` or an equivalently reviewed local path, outside Git, Zenodo, chat, and cloud-synced Documents. External scoring and inference require frozen GEO-trained artifacts and the published second lock.
- The current site decision is [documented separately](PPMI_MAC_ONLY_EXECUTION_DECISION_20261020_EN.md). Earlier acquisition and preparation records remain historical provenance and must be reported accurately; they are not instructions for future execution.

## Required fields before publication

| Item | Current status |
|---|---|
| Final PPMI eligibility and one-to-one crosswalk | Pending; provisional counts are not final |
| Plate × diagnosis metadata count | Completed for provisional 472 candidates; final population still pending |
| Existing input preparation code, release, output schema and SHA-256 | Recorded; verify binding to final local input |
| Local private bundle SHA-256 | Checked; `1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c` |
| Mac-only storage and applicable PPMI data-use terms | Local secure path established; terms and earlier sync history still under review |
| Mac scoring code, environment, exact command and portability checks | Pending |
| Frozen final GEO M1–M5 and preprocessing artifact hashes | Pending public manifest |
| Final GEO-trained auxiliary sex model hash | Fit and independently audited on GEO controls; SHA-256 `aa92a1516328477b1c4ba1a61e7bcec43263cc38574f57ec99c5f26465348200`. Include this exact artifact in the public lock manifest; no PPMI score has been run. |
| Control gates, bootstrap rule and private output location | Sex gate and its 5,000-resample path specified below; complete the other registered controls, scorer code/hash, and final private output manifest before publication |
| Exact Git commit/tag and public Zenodo version/DOI | Not yet published |

## Metadata inputs still needed on the Mac

The local 472-candidate bundle contains only **selected QC-pass candidates**, so it cannot supply a pre-QC denominator or group-specific QC-fail rates. Those counts require the complete corresponding IR3/B38 RNA release metadata (`metaDataIR3`, including `QCflagIR3`), the matching full portal sample metadata (including `QC` and Long RNA status), and the diagnosis/status source used to classify the same pre-QC baseline pool. Define the denominator before applying QC, show each exclusion stage separately by historical PD/control group, and reconcile the surviving rows one-to-one with the locked candidate crosswalk. Aggregate archive totals are not substitutes.

The archived whole-blood collection export is already on the Mac and provides `PATNO`, `EVENT_ID`, and `BLDCLLDT`, but its recorded date has only month/year precision. The current candidate bundle contains no medication-at-draw history, and no medication export was found locally. Obtain the applicable medication record and field codebook directly on the Mac. Join records by participant and sampling visit/date under a rule frozen before examining group contrasts; report missingness and unresolved timing rather than inferring treatment status from diagnosis or a later visit. A month-only date cannot establish exact same-day medication use. These are metadata audits only; they must not calculate expression-derived scores.

## Frozen PPMI sex-control scoring path

After the completed second public lock, run the fixed GEO-control-trained auxiliary model (`sex_model.npz`, SHA-256 `aa92a1516328477b1c4ba1a61e7bcec43263cc38574f57ec99c5f26465348200`) on the **healthy controls in the locked final PPMI population**. The currently linked candidate set has 161 controls; this is provisional, so the final lock must record the actual control count before scoring. Require unique participant linkage, recorded sex, and all five marker values for each included control. An unresolved linkage or transformation discrepancy suspends the gate rather than changing the cohort after seeing scores.

Use the same five registered Y genes, in the saved model order: `KDM5D` (`ENSG00000012817`), `DDX3Y` (`ENSG00000067048`), `RPS4Y1` (`ENSG00000129824`), `UTY` (`ENSG00000183878`), and `EIF1AY` (`ENSG00000198692`). Match by exact Ensembl ID, not by a new symbol lookup. Read their values from the already prepared, frozen PPMI **within-sample rank-percentile matrix**. Apply only the five GEO-control means and population SDs stored with the final model, followed by its stored intercept and coefficients through the locked prediction implementation. Do not recalculate means, SDs, ranks across samples, coefficients, calibration, threshold, or any other fitting statistic from PPMI; the within-sample ranking rule itself is fixed in protocol v3.1. Check marker order, finite values, model hash, matrix/gene-manifest hash, and synthetic/GEO-only scorer parity before reading PPMI predictions.

Calculate AUROC of the fixed model's probabilities against **recorded sex** among final PPMI controls only. For its 95% percentile CI, perform exactly 5,000 fixed-model donor bootstrap replicates, sampling with replacement **within each recorded-sex stratum** and preserving each stratum's original count. Use `geo_pipeline.scope_seed("ppmi/controls/sex/bootstrap")`, which yields seed `3222988182` under master seed `20261014`, and take the 2.5th and 97.5th percentiles of replicate AUROCs. Record the exact bootstrap code hash in the final lock. No refitting occurs in bootstrap resamples. The operational gate passes only if the observed AUROC is at least 0.95 **and** the CI lower bound is at least 0.90. Balanced accuracy at threshold 0.5 is descriptive only. If the gate fails or is unevaluable, suspend confirmatory PD scoring and interpretation under protocol v3.1; do not retune, recalibrate, overwrite sex labels, exclude discordant controls, or substitute another marker model based on PPMI results.

Keep all participant-level inputs, probabilities, and bootstrap indices in local nonsynced private storage on the researcher's Mac. Publish only permitted aggregate gate results after the lock and data-use review. This paragraph specifies the future execution path; **no PPMI sex score has yet been calculated**.

## Storage and data-use evidence required in the final lock

The intended execution and private storage location is the researcher's personal Mac, with FileVault enabled and participant-level files under local, nonsynced `~/PPMI_private` (owner-only directory and files). The applicable signed Tier 2 terms, the already completed historical device transfer, and any earlier iCloud synchronization still require review and an accurate compliance statement before a public second lock. The current local storage setup alone is not evidence that all past handling complied with those terms.

No field marked pending may be presented as satisfied. The final lock will reference private participant-level files only by non-identifying hashes and controlled local paths; it will not publish those files.
