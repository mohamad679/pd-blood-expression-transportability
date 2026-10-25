# GEO within-batch AUROC diagnostic decision (2026-10-20)

## Question and scope

Assess whether the already-fitted GENEPARK (GSE99039) M4 and M5 models separate PD and control samples within the same recorded batch. This is a descriptive GEO development diagnostic. It does not use PPMI, fit or tune a model, or make a causal claim about batch.

## Frozen inputs and unit of analysis

- Use only the five original frozen GSE99039 outer-test prediction files in `private_data/geo_runs_20261018/GSE99039_outer_{0..4}/heldout_predictions.npz`.
- Join the saved held-out sample IDs to the existing GSE99039 sample metadata by GEO sample ID to obtain batch and diagnosis.
- Require each prediction file's ordered sample IDs and labels to match its frozen outer-test split and metadata exactly. Every GSE99039 sample must occur in exactly one outer test fold.
- Use only saved M4 and M5 held-out probabilities. No refit, retuning, recomputation of scores, or use of in-sample predictions is allowed.

## Estimand and aggregation

For each recorded batch, calculate the AUROC among all out-of-fold samples assigned to that batch, comparing only PD-control pairs within the batch. A batch is evaluable only when the pooled out-of-fold set contains at least 5 PD and 5 control samples. List counts and the reason for non-evaluability for every batch. The primary summary is the unweighted macro-average of the batch-specific AUROCs over evaluable batches, reported separately for M4 and M5; also show each eligible batch's AUROC. No pooled cross-batch AUROC, significance test, or selection based on results will be reported.

The eligibility threshold applies to the full outer out-of-fold cohort within each batch, not separately to each fold. Each sample contributes one held-out prediction from its own outer fold. This provides enough observations per batch for the prespecified descriptive diagnostic while preserving out-of-fold predictions.

## Interpretation

Within-batch AUROC above 0.5 indicates that saved model predictions rank cases above controls within that batch in these held-out samples. It weakens a purely between-batch explanation, but does not prove a biological disease signal or rule out other confounding. Small numbers of eligible batches and noisy batch-specific AUROCs limit generalization. Report the per-batch sample sizes alongside the estimates.

## Data handling

Inputs are GEO development data only. The script writes only aggregate batch counts and AUROCs; it does not export sample IDs or participant-level records. PPMI files and scores are never accessed.
