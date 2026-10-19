# GEO age/sex-matching pilot report

**Scope:** `GSE99039/outer/1`  
**Purpose:** implementation and runtime pilot only; no outer-test performance estimate.

## Outcome

The pilot completed successfully in 873.2 seconds wall time (14.55 minutes) and 2,854.4 process-CPU seconds. It matched 70 outer-training pairs and 14 outer-test pairs. The outer test was not scored. All five inner-validation partitions met the preregistered minimum of five pairs; their pair counts were 8, 9, 5, 11, and 10. Inner-training rematched pair counts were 50, 51, 49, 52, and 52.

The preprocessing fit-sample hashes exactly matched the expected four-/five-fold inner-training memberships plus the final matched outer-training set. No donor crossed an inner train/validation boundary. All 300 recorded fit attempts passed their numerical gates; the maximum recorded KKT value was `2.124e-8`, below the `1e-6` limit. No solver or prediction warnings were recorded. PPMI was not loaded, and outer-test performance metrics were not computed.

## Finding that needs review before the full run

M2 and M5 selected the lowest lambda on the final expanded grid (`1e-6`, `rho=0`), with nonzero gene blocks. Both therefore remained on the lower search boundary after the allowed two expansion rounds. M3 selected `lambda=0.01`, `rho=0`. The M2/M5 boundary selections are a tuning-instability signal, not a failed numerical fit; all KKT checks passed. The preregistered grid will not be changed based on this pilot.

## Runtime implication

The earlier one-hour estimate is too optimistic to rely on. Six tuning-eligible scopes remain under the registered count rule (four GENEPARK outer folds for separate descriptive reporting, the complete GENEPARK-to-PROBE transfer, and pooled/final descriptive tuning). At the pilot's runtime, six similarly sized scopes alone would take about 87 minutes; pooled/final has more matched training data and may take longer. This is a rough extrapolation, not a runtime guarantee.

Because the pilot showed persistent lower-bound selections for two models and ran substantially longer than expected, the condition “if nothing unusual” is not met. Do not start the full sensitivity run until this finding and revised time estimate have been reviewed. No fold, model, grid, caliper, or reporting rule was changed in response to the pilot.

Raw pilot fit records, model artifacts, warning log, and machine-readable result are under `private_data/geo_age_sex_matching_pilot_20261019_GSE99039_outer_1/`.
