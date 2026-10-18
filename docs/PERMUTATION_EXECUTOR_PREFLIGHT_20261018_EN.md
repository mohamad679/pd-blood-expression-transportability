# Permutation executor preflight

**Date:** 2026-10-18. The executable scientific scope is B=200, M4/M5, ten within-GEO outer folds. No PPMI inputs are supported. No real-GEO labels have been permuted or scored during this preflight.

## CBC interpretation fixed before per-marker inspection

The amendment explicitly retains a marker-specific minimum of 30 valid pairs. With only 31 metadata candidates, missing values, incompatible units or conflicting records can legitimately make a marker non-evaluable. Matching, population, time window and visit priority will not be reopened to increase n. A passing association at approximately 30 participants will be described as limited supportive evidence. Same-month matching omits some possible adjacent-month ±30-day pairs; monthly dates cannot determine their true eligibility. Marker completeness and associations have not been calculated.

## Synthetic tests and resource results

The final synthetic run passed **32 checks** and **650 fit attempts**, all meeting KKT <= 1e-6; the largest residual was **2.6839175860937647e-8**. Eight additional 1,000-gene numerical stress fits were used for the resource comparison. The full source universe in the final resource fixture is 16,364 genes, with 381 synthetic participants and the real source sizes; fits select 1,000 genes and use 191 training participants in the larger source. These are artificial null data and provide no clinical performance evidence.

- M4/M5 preprocessing, candidate scores, selected decisions, coefficients and predictions match the frozen five-model engine in the synthetic equivalence test.
- Donor-label maps are invariant to row order and scheduling, preserve fixed outer-block class totals and apply consistently to repeated donor rows. Inner validation membership is regenerated from current labels; the identity case reproduces the frozen original inner membership.
- One-worker/one-thread and two-worker/two-thread runs give identical label-map hashes and fold AUROCs in the synthetic smoke runs. Two null repetitions are an implementation smoke test, not statistical null calibration.
- Resume reuses checked complete repetitions, rejects changed provenance and tampered artifacts, and preserves incomplete attempts while rerunning the same repetition cold. Controller and per-repetition locks reject concurrent ownership.
- Real randomized execution rejects absent public-amendment/preflight evidence. Completion requires every dependent fit/fold; no shortened denominator or silent dropping is implemented.

The retained resource configuration is **two workers, two threads per worker**. On this 8 GiB Mac, the final benchmark's conservative budget (sum of worker peak RSS, parent peak RSS and a 1 GiB guard) was **2.033 GiB**. System compressor swapouts increased by **0** pages; general pageouts increased by **1,645** pages. Four stress jobs took **17.46 s** with one worker and **11.14 s** with two workers, including process startup. Stress-fit probabilities were identical across configurations. This microbenchmark does not establish a full-permutation ETA or guarantee future system load.

## Preserved preliminary attempts and resource-gate correction

The preliminary four-worker run passed its numerical/equivalence checks but failed an initially combined zero-pageout/zero-swapout assertion. Its raw counter values were not persisted before that assertion; no numerical swap increase is asserted for that run. The failed report and private log are retained.

A fully instrumented two-worker rerun passed 31 checks and measured zero additional swapouts but 1,154 additional general pageouts, causing the same combined assertion to fail. All its counter values and fit results are retained in `synthetic_preflight_v2.json`. Its conservative budget was 1,730,355,200 bytes with the smaller 1,200-gene raw fixture.

Before the final full-gene-universe rerun, the resource criterion was corrected to distinguish the counters: macOS defines general pageouts separately from compressor swapouts. See [Apple's XNU counter definitions](https://github.com/apple-oss-distributions/xnu/blob/main/osfmk/mach/vm_statistics.h). Pageouts remain a reported diagnostic, while the operational acceptance criterion is no growth in compressor swapouts plus the conservative memory budget below 75% of physical RAM. This is an execution-resource clarification before publication and before any real randomized fit; the statistical pipeline, KKT criterion and seeds were unchanged. A machine-wide counter does not attribute activity to a specific worker.

## Observed-label identity gate

**Running; not yet certified.** The separate identity command refits only M4/M5 on the ten original GEO outer folds, using the original labels and reconstructed original inner folds. It must reproduce all original tuning candidates and selected decisions, every held-out membership/outcome array, AUROC within 1e-12 and probabilities within absolute 1e-6. It permutes no real labels. A public aggregate certificate can be created only after this gate passes with the same executor/configuration/environment binding.

The original PROBE fold-0 pilot was recovered from the prior feasibility workspace and copied into the ignored complete-run backup. All 13 of its previously recorded artifact hashes were verified before use; no predictions or models are published.

## Saved-fit membership and preprocessing audit

A separate read-only audit passed for **all 13 original GEO scopes and 3,415 fit attempts**. It reconciled every logged training-membership hash, regenerated training-only preprocessing fingerprint, sample count and weight-derived Kish size with the frozen outer/inner splits. Saved final preprocessing arrays and held-out membership/outcome arrays also matched. It performed zero new model fits or real-label permutations. Evidence: `geo_membership_audit.json`; implementation: `pipeline/audit_geo_membership.py`.

This verifies separation under the operational donor IDs available to the pipeline. PROBE's GSM fallback and source prefixes do not establish biological identity or cross-source donor independence. The audit does not claim to rule out every possible unlogged processing error.

## Commands and publication gate

`bash pipeline/run_permutation_preflight.sh` performs synthetic/resource checks, plans the frozen real membership, executes the observed identity replay and creates an aggregate certificate if all gates pass. It cannot launch real permutations. The real `run` command additionally requires a receipt for the publicly verified Zenodo amendment, matching exact executor/config/input/environment/split hashes and the passed certificate hash.

The source engine `pipeline/geo_pipeline.py` remains byte-identical to its public lock, SHA-256 `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`. Checkpoints, fit logs, models and subject-level outputs remain under ignored private directories. Public files contain source code, configuration and aggregate evidence only.
