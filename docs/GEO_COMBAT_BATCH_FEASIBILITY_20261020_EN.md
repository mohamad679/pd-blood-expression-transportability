# GEO ComBat batch feasibility count (metadata only; 2026-10-20)

This metadata-only feasibility check reconstructs the frozen GEO split manifest from sample metadata. It does not read gene-expression values, fit a correction, compute scores, or access PPMI.

Frozen split-manifest SHA-256 verified: `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.

## Impact of the initial whole-partition exclusion rule

The initial proposed rule marked an entire fit/validation or fit/test partition non-evaluable if any training batch had fewer than three samples.

| Cohort | Outer train/test partitions invalid | Inner train/validation fits invalid | Inner invalid by outer fold |
|---|---:|---:|---|
| GSE57475 (PROBE) | 0/5 | 8/25 | outer/0: 1; outer/1: 3; outer/2: 3; outer/3: 1; outer/4: 0 |
| GSE99039 (GENEPARK) | 5/5 | 24/25 | outer/0: 5; outer/1: 5; outer/2: 5; outer/3: 5; outer/4: 4 |

This rule would make GENEPARK ComBat evaluation practically unavailable, including every outer test fold. The counts support changing the handling rule before examining any ComBat outcome.

## Registered handling after the count

Use **partial ComBat** in every partition:

1. Within a training partition, estimate ComBat parameters only from batches with at least three training samples.
2. Keep training samples from smaller batches in the model fit, but leave their expression at the original within-sample rank values, without ComBat correction.
3. Apply frozen estimates to held-out samples only when their batch was included in the corresponding training fit. Leave held-out samples from small or unseen batches at their original within-sample ranks. Never estimate batch parameters from held-out expression or labels in an inductive arm.
4. In the all-expression arm, use expression and batch metadata from the full designated development scope, with no labels, and leave full-scope batches smaller than three uncorrected. Label this arm transductive.
5. Report corrected/uncorrected sample counts separately for training and held-out data in every partition. This sensitivity is partial adjustment, not complete batch removal.

The disease-preserving arm must additionally pass a training-only design-rank check after small batches are omitted from its ComBat fit. If batch and disease are not identifiable in that training design, that arm/partition is non-evaluable. No batch merging, population change, or outcome-informed repair is permitted.

## Full-cohort batch inventory

| Cohort | Batch count | Batches with n ≥ 3 | Samples in batches with n ≥ 3 | Small batches |
|---|---:|---:|---:|---|
| GSE57475 (PROBE) | 16 | 16 | 142/142 | None |
| GSE99039 (GENEPARK) | 11 | 10 | 238/239 | `08/27/09_TRO` (n=1) |

The threshold was also checked across the five outer training folds: both cohorts have multiple batches meeting it in every fold. The 3-sample cutoff avoids estimating within-batch variance from singleton batches and avoids relying on a two-sample variance estimate for this diagnostic. It does not automatically make the disease-preserving model identifiable; the explicit design-rank check remains required.

## Scope and hashes

- Metadata sources: `prepared/GSE57475_sample_metadata.json` and `prepared/GSE99039_sample_metadata.json`.
- The frozen manifest hash above was verified before counting.
- No expression matrix, held-out prediction, ComBat fit, performance metric, or PPMI file was accessed.
- The outcome-independent rule change is documented in [GEO_COMBAT_AUDIT_DECISION_20261019_EN.md](GEO_COMBAT_AUDIT_DECISION_20261019_EN.md).
