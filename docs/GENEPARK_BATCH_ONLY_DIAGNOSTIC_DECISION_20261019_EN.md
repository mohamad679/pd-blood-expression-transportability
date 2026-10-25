# GENEPARK batch-only diagnostic decision (2026-10-19)

## Frozen definition

This is a metadata-only, cross-validated diagnostic on GSE99039 (GENEPARK/TRO). It uses the five already frozen `GSE99039/outer/0..4` test folds from the v3.1.1 split manifest; it does not create or alter any split.

For each outer fold independently, estimate the PD proportion for every batch from that fold's **outer-training samples only**:

`score(batch b) = number of outer-training PD samples in b / number of outer-training samples in b`.

Give each outer-test sample the score for its recorded batch. If a test batch is absent from that fold's training set, use the overall PD proportion in that outer-training set and record the number of affected samples. There is no smoothing, fitting/tuning, expression, model selection, sample exclusion, or use of test labels to form scores. Test labels are used only to calculate the fold AUROC. Tied scores are handled by the standard AUROC definition.

Report all five fold AUROCs and their unweighted mean, consistent with the existing within-cohort outer-fold aggregation. Show the existing M1–M5 results on those same five folds alongside this diagnostic. Do not use this diagnostic as a replacement endpoint or to modify the models, splits, or interpretation rules.

## Interpretation boundary

This estimates how much discrimination is available from the association between the recorded batch and case status in held-out samples whose batch labels are represented in training. It is a useful shortcut diagnostic, but it is not causal evidence that batch explains biological performance, and it does not establish that batch is the sole cause of any model's AUROC. In particular, same-batch train/test observations permit learning a batch-specific disease prevalence. Report it as a **batch-only metadata diagnostic**, not as a deployable disease classifier or evidence of disease biology.

The primary dataset summary is descriptive. No additional significance test or threshold for declaring “explained” is introduced after observing the result.

## Data separation and provenance

Read only the prepared GSE99039 sample metadata (sample ID, binary case label, recorded batch) and the frozen split-manifest inputs needed to reconstruct the published folds. Expression matrices, PPMI files, scores, and model predictions are not read by the batch-only computation. Verify the reconstructed split manifest SHA-256 against the published value before computing AUROCs. Preserve the generated fold table, summary, and script as analysis provenance.
