# Preregistered decision: GEO age/sex matching sensitivity

**Status:** This decision was fixed before conducting the metadata count audit and before any matching-sensitivity model fit. This document is the authoritative English version. The [Persian translation](GEO_AGE_SEX_MATCHING_DECISION_20261019_FA.md) is retained for convenience.

## Population and matching algorithm

Perform one-to-one matching of PD and control samples separately within each GEO source and recorded-sex stratum. Require an absolute age difference no greater than three years and use each donor at most once. Enumerate all eligible case-control edges, sort by `(absolute age gap, case donor ID, control donor ID)`, then accept an edge in order only if neither donor has already been matched. If those keys tie completely, use case sample ID and then control sample ID as deterministic tie-breakers.

This is a deterministic greedy algorithm. It is not guaranteed to maximize the number of pairs or minimize total age distance. This limitation is intentional; the method will not be changed after inspecting results.

## Fold-safe order

For every outer scope, first match outer-train and define the matched outer-training subset `S`. Restrict the frozen inner folds to members of `S`. Within each inner fold, rematch the restricted inner-training side and inner-validation side separately. A counterpart assigned to the opposite side cannot be used in that partition; such a sample remains unmatched and is excluded there. No sample outside `S` may enter inner tuning, and no matched pair may cross the train/validation boundary.

Match outer-test independently within source and recorded-sex strata. Test matching uses disease labels to define the target population; any resulting estimate must therefore be called **performance in the matched population** and must not be presented as directly comparable to the primary full-population estimate. For transfer scopes, match the source training cohort and target test cohort independently. For `pooled/final`, match the pooled training set to `S` and use it only for descriptive inner hyperparameter selection; it has no held-out test or performance estimate.

## Minimum size and evaluability

Every outer test and every inner validation must retain at least five participants in each class after matching (at least five matched pairs). If an outer test is below this threshold, its held-out estimate is unevaluable. If any inner validation in an outer scope is below threshold, the entire scope is unevaluable for tuning. Do not drop or replace folds, change the caliper, or fit models for that scope. Report counts and exclusions, plus age-gap summaries, for every partition, including inner-training and inner-validation partitions, before any model fit.

## Reporting rule after the metadata count audit

Do not report the mean of the four remaining eligible GENEPARK outer folds as “matched GENEPARK cross-validation.” Declare the complete within-cohort GENEPARK CV result unevaluable for this sensitivity; the four eligible folds may be shown individually for descriptive purposes only. If it passes the thresholds, the GENEPARK-to-PROBE transfer is a complete, reportable scope. The `pooled/final` branch records descriptive inner hyperparameter selection only and is not used in any PPMI evaluation. PPMI evaluation uses the original locked models without refitting.

## Metadata-only count audit

Before fitting, generate a count table for outer-train before/after matching, independently matched outer test, and each inner-train/validation after restriction to `S` and side-specific rematching. Include class counts, number of pairs, matched and unmatched sample/donor counts, mean/median/maximum absolute age gap, and threshold status. This audit uses only prepared GEO metadata and frozen split identifiers. It does not access PPMI, fit models, compute model scores, or calculate performance metrics.

The split source is the frozen v3.1.1 GEO split manifest. The eligible analysis set and reporting labels above are fixed by this decision.
