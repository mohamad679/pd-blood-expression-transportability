# GEO panel-gene exclusion sensitivity — implementation decision

**Status:** fixed before the panel-exclusion pilot and full sensitivity run.

## Fixed implementation

1. Use the six exact gene-ID lists in `protocol/v3_marker_manifest.json`.
   Their union contains **172 unique genes** (173 panel memberships); `TECPR2`
   (`ENSG00000196663`) occurs in both the Neutrophils and Platelet panels. The
   manifest SHA-256 is
   `729a9c450deadd6c65d0949c1e8f09ee76b455989a91580a60d41a1a758a8bc1`.
2. In every inner-training and outer/final-training partition, remove those
   IDs from the M5-added gene universe **before** source-stratified positive
   MAD ranking. Select up to 1,000 eligible remaining genes using the locked
   MAD ordering and gene-ID tie break. Validation and held-out samples do not
   influence selection or scaling.
3. Refit/tune **M5 only**. Reuse the frozen published M4 inner-CV score as
   M5's no-gene candidate and reuse the published M1–M4 model artifacts;
   do not refit M1–M4. The six panel scores continue to use their full frozen
   gene sets. Thus M4 is unchanged and the intervention is confined to the
   added gene block.
4. If the M5 tuning rule selects its no-gene baseline candidate, reuse M4
   exactly and report **“gene block not selected”**. If a penalized candidate
   is selected but the final M5 fit has zero penalized gene coefficients,
   likewise reuse M4 and report **“gene block not selected (final gene block
   is zero)”**. Do not describe either branch as an observed zero effect.
5. Use the original frozen 13 scopes, split/code/input hashes, candidate grid,
   class/source weighting, KKT gates, retry rule, and held-out scoring order.
   A single-scope pilot checks implementation and timing only; it is not an
   effect estimate and does not alter this decision. No PPMI data or scores are
   used.

## Descriptive baseline audit, before exclusion

This audit is computed from the already saved original `preprocessing.npz` and
`M5_model.npz` artifacts. It required no model fit. “Selected panel genes” is
the intersection of the original final 1,000-gene M5 selection and the 172-ID
union. “Nonzero panel coefficients” counts exact nonzero final M5 coefficients
among those selected genes, using the numerical protocol's `beta != 0`
definition. These are descriptive pre-sensitivity quantities, not evidence of
performance impact.

| Frozen scope | Selected panel genes | With nonzero M5 coefficient |
|---|---:|---:|
| GSE57475 outer 0 | 13 | 0 |
| GSE57475 outer 1 | 13 | 0 |
| GSE57475 outer 2 | 14 | 14 |
| GSE57475 outer 3 | 12 | 0 |
| GSE57475 outer 4 | 11 | 11 |
| GSE57475 transfer | 13 | 1 |
| GSE99039 outer 0 | 16 | 16 |
| GSE99039 outer 1 | 14 | 14 |
| GSE99039 outer 2 | 16 | 1 |
| GSE99039 outer 3 | 16 | 0 |
| GSE99039 outer 4 | 15 | 15 |
| GSE99039 transfer | 15 | 15 |
| Pooled final | 13 | 1 |
| **Total across scope-specific fits** | **181** | **88** |

The 181 selected memberships are summed across 13 separate fits; they are not
181 distinct genes. Counts range from 11–16 per scope (1.1–1.6% of each
1,000-gene selection). This indicates that relatively few marker-union genes
were selected in the original fits, while some had nonzero coefficients. It
does not determine how their removal will affect refitting or held-out
performance.
