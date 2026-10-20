# PPMI post-lock secondary analyses: execution record

**Status:** Operational details written after the primary PPMI result and before calculating the secondary results below. Protocol v3.1 and its published CBC amendment govern. These analyses cannot replace or modify the sole confirmatory endpoint, its population, fitted models, or inference rule. No PPMI result from this stage may enter a public repository, preprint, or manuscript before DPC review.

## Inputs and integrity

Use the 472-row baseline bundle and row-ordered M1–M5 predictions from the completed first PPMI score. Verify the published second-lock manifest, all three locked input digests, recorded sex and historical group order, one row per donor, and the primary M4/M5 AUROCs against the completed aggregate. Never regenerate M1–M5 predictions or fit a disease model. All outputs stay under owner-only `~/PPMI_private`; print aggregates only.

## Secondary fixed-score metrics

Report the case fraction and average precision (AUPRC) for all five models. Describe calibration with the already reported log loss, Brier score, and ten fixed-width probability bins `[0,0.1),...,[0.9,1]` for each model. A bin is a descriptive count, mean probability, and observed PD fraction; no recalibration, threshold tuning, or inferential claim is attached to it. Report M3−M1 AUROC with the same case/control-stratified paired donor bootstrap algorithm and 5,000 draws used for the primary comparison; label it secondary and conditional on the fixed GEO fits. The seed for this secondary draw is 20261014, shared for a directly comparable resampling scheme.

## Matched PPMI sensitivity

Match PD and control donors one-to-one within recorded-sex strata. Enumerate every eligible edge with absolute baseline-age difference ≤3 years, sort by `(age gap, PD donor ID, control donor ID)`, and greedily accept an edge only when neither donor has already been used. Because each locked donor has exactly one baseline sample, donor IDs fully resolve ties. The matching uses historical labels solely to define the matched target population. Apply the saved M1–M5 predictions to the selected rows; do not refit, retune, rescale, or recalibrate. Report matched pairs, class/sex counts, unmatched counts, age-gap summary, AUROC and AUPRC for each model. Require at least five matched pairs for performance estimation, as in the published GEO matching decision. For M5−M4 uncertainty, resample whole matched pairs with replacement 5,000 times, seed 20261014; use a percentile CI. This sensitivity estimates performance in the matched population and is not directly comparable to the full-population primary estimate.

## Age control

Apply each of the two **already fitted** GEO-control age models (`GSE57475` and `GSE99039`) separately to the 161 locked PPMI controls, matching exact Ensembl IDs to the locked rank matrix and using only stored GEO means, SDs, coefficients and intercepts. Report Spearman correlation and MAE versus each model's stored training-median-age predictor. No PPMI fitting, age-target selection, or new pooled age model is allowed. The pre-lock GEO transfer report did not demonstrate age-control capacity; a poor PPMI age result is descriptive and is not an exclusion or a failure gate.

## CBC and longitudinal controls

For CBC, retain the published v3.1.4 rule: controls first; exact participant; the RNA baseline-collection month; CBC rows at any visit in the same calendar month; visit priority BL, then SC, then other visits; no lower-tier fallback following a conflict. Collapse only identical repeated measurement records. Multiple collection occasions or conflicting numeric values/units at the highest tier are excluded for that participant and assay. Use the current export, not the archived export, as the analysis source. Read numeric SI results with a nonempty identical SI unit within an assay; rows without an interpretable common unit are excluded and counted. Compare the neutrophil panel with HMT8/HMT15, monocytic panel with HMT10/HMT17, and platelet panel with HMT13. HMT9/HMT16 (total lymphocytes) may describe T/B panel association but cannot validate T or B cells separately. HMT3/RBC does not validate the erythroid-reticulocyte panel. The assay names are identified in the PPMI Data User Guide. For each assay, report eligible n, exclusions, units, Spearman and 5,000 donor-bootstrap percentile CI (seed 20261021), and the prespecified supportive rule `n ≥ 30`, `rho ≥ 0.4`, lower CI >0. If a usable assay has fewer than 30 pairs or constant values, report non-evaluable; never select a new assay based on correlation.

Longitudinal similarity requires nonbaseline sample-level autosomal ranks. The locked 472-row bundle contains baseline candidates only. If no independently verified nonbaseline rank matrix is locally available, report the diagnostic as unavailable, with the missing input and no inferred identity claim. This is a nonstopping control; do not redownload a large archive or alter the baseline population to force it.

## Publication boundary

The PPMI publication policy calls for DPC administrative review of external manuscripts via **Submit a Manuscript or Abstract** before journal submission. The project adopts the stricter rule that any public PPMI result, including a preprint or Zenodo version with PPMI numbers, waits for that review. Keep the review packet and all PPMI numerical results outside tracked public project files until then.
