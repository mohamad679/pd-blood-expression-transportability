# Permutation executor preflight

**Date:** 2026-10-18. **Current status: resource gate not passed; publication and real permutations not started.** The executable scientific scope is B=200, M4/M5, ten within-GEO outer folds. No PPMI inputs are supported. No real-GEO labels have been permuted or scored during this preflight.

## CBC interpretation fixed before per-marker inspection

The amendment explicitly retains a marker-specific minimum of 30 valid pairs. With only 31 metadata candidates, missing values, incompatible units or conflicting records can legitimately make a marker non-evaluable. Matching, population, time window and visit priority will not be reopened to increase n. A passing association at approximately 30 participants will be described as limited supportive evidence. Same-month matching omits some possible adjacent-month ±30-day pairs; monthly dates cannot determine their true eligibility. Marker completeness and associations have not been calculated.

## Synthetic tests and resource results

Before the terminal-failure resume guard was hardened, a synthetic run passed **32 checks** and **650 fit attempts**, all meeting KKT <= 1e-6; the largest residual was **2.6839175860937647e-8**. Eight additional 1,000-gene numerical stress fits were used for the resource comparison. The full source universe in the resource fixture is 16,364 genes, with 381 synthetic participants and the real source sizes; fits select 1,000 genes and use 191 training participants in the larger source. These are artificial null data and provide no clinical performance evidence. The final release-candidate status is described below and is not certified as passed.

- M4/M5 preprocessing, candidate scores, selected decisions, coefficients and predictions match the frozen five-model engine in the synthetic equivalence test.
- Donor-label maps are invariant to row order and scheduling, preserve fixed outer-block class totals and apply consistently to repeated donor rows. Inner validation membership is regenerated from current labels; the identity case reproduces the frozen original inner membership.
- One-worker/one-thread and two-worker/two-thread runs give identical label-map hashes and fold AUROCs in the synthetic smoke runs. Two null repetitions are an implementation smoke test, not statistical null calibration.
- Resume reuses checked complete repetitions, rejects changed provenance and tampered artifacts, and preserves incomplete attempts while rerunning the same repetition cold. Controller and per-repetition locks reject concurrent ownership.
- Real randomized execution rejects absent public-amendment/preflight evidence. Completion requires every dependent fit/fold; no shortened denominator or silent dropping is implemented.

The candidate resource configuration is **two workers, two threads per worker**. On this 8 GiB Mac, the pre-hardening benchmark's conservative budget (sum of worker peak RSS, parent peak RSS and a 1 GiB guard) was **2.033 GiB**. System compressor swapouts increased by **0** pages; general pageouts increased by **1,645** pages. Four stress jobs took **17.46 s** with one worker and **11.14 s** with two workers, including process startup. Stress-fit probabilities were identical across configurations. This microbenchmark does not establish a full-permutation ETA or guarantee future system load, and it does not override subsequent failed resource measurements.

## Preserved preliminary attempts and resource-gate correction

The preliminary four-worker run passed its numerical/equivalence checks but failed an initially combined zero-pageout/zero-swapout assertion. Its raw counter values were not persisted before that assertion; no numerical swap increase is asserted for that run. The failed report and private log are retained.

A fully instrumented two-worker rerun passed 31 checks and measured zero additional swapouts but 1,154 additional general pageouts, causing the same combined assertion to fail. All its counter values and fit results are retained in `synthetic_preflight_v2.json`. Its conservative budget was 1,730,355,200 bytes with the smaller 1,200-gene raw fixture.

Before the final full-gene-universe rerun, the resource criterion was corrected to distinguish the counters: macOS defines general pageouts separately from compressor swapouts. See [Apple's XNU counter definitions](https://github.com/apple-oss-distributions/xnu/blob/main/osfmk/mach/vm_statistics.h). Pageouts remain a reported diagnostic, while the operational acceptance criterion is no growth in compressor swapouts plus the conservative memory budget below 75% of physical RAM. This is an execution-resource clarification before publication and before any real randomized fit; the statistical pipeline, KKT criterion and seeds were unchanged. A machine-wide counter does not attribute activity to a specific worker.

## Observed-label identity gate

**First replay passed; final binding still pending.** The separate identity command refit only M4/M5 on the ten original GEO outer folds, using the original labels and reconstructed original inner folds. All original candidate/selection dictionaries, held-out membership and outcome arrays, and 20 AUROCs matched. Maximum absolute probability difference was **1.1102230246251565e-16**, below the preflight tolerance of 1e-6 (AUROC tolerance 1e-12). It permuted no real labels. This replay belongs to the pre-hardening executor binding at Git commit `dbdc878`; it is retained and is not relabeled as evidence for the later whole-file hash. The final-binding replay remains required after the current resource gate passes. A public aggregate certificate cannot yet be created.

The original PROBE fold-0 pilot was recovered from the prior feasibility workspace and copied into the ignored complete-run backup. All 13 of its previously recorded artifact hashes were verified before use; no predictions or models are published.

## Final resume hardening and current execution blocker

Code review identified that a repetition with an explicit terminal `FAILURE.json` must not be treated like an ordinary interrupted attempt on restart. The executor now refuses to retry that invalid dependent analysis. An interruption without terminal failure remains resumable by a fresh cold attempt with the same seed. A new synthetic injected-failure test passes. The public receipt also binds the complete certificate input/config/environment/split identity, and verifier/finalizer source hashes are included in the bound code manifest. The fitting, tuning and preprocessing algorithms are unchanged.

The final release candidate has **33 checks**. Two runs passed the first **32**, including all 650 numerical fits, the terminal-failure guard, worker/seed/resume invariance and the 1,000-gene stress checks. Both failed the last resource criterion:

| Retained report | Competing observed replay | Memory budget, GiB | Additional Swapouts, pages | Additional Pageouts, pages |
|---|---|---:|---:|---:|
| `synthetic_preflight_release.json` | Active | 2.208 | 91,386 | 1,692 |
| `synthetic_preflight_idle.json` | Finished | 2.124 | 23,780 | 3,190 |

On this Mac's 4,096-byte pages, the second increase is approximately 92.9 MiB. These are system-wide counters and do not identify which process caused the activity. They nevertheless fail the predefined no-additional-compressor-swapout criterion. The memory budget itself is below 75% of physical RAM; the blocker is the observed resource condition, not an unresolved KKT failure. No acceptance threshold was relaxed to pass these runs. No amendment is published and no real permutation is launched. Next: rerun the prepared Bash preflight when the host has less competing memory pressure or in a suitable execution environment, retaining these attempts.

## Saved-fit membership and preprocessing audit

A separate read-only audit passed for **all 13 original GEO scopes and 3,415 fit attempts**. It reconciled every logged training-membership hash, regenerated training-only preprocessing fingerprint, sample count and weight-derived Kish size with the frozen outer/inner splits. Saved final preprocessing arrays and held-out membership/outcome arrays also matched. It performed zero new model fits or real-label permutations. Evidence: `geo_membership_audit.json`; implementation: `pipeline/audit_geo_membership.py`.

This verifies separation under the operational donor IDs available to the pipeline. PROBE's GSM fallback and source prefixes do not establish biological identity or cross-source donor independence. The audit does not claim to rule out every possible unlogged processing error.

## Commands and publication gate

`bash pipeline/run_permutation_preflight.sh` performs synthetic/resource checks, plans the frozen real membership, executes the observed identity replay and creates an aggregate certificate if all gates pass. It cannot launch real permutations. The real `run` command additionally requires a receipt for the publicly verified Zenodo amendment, matching exact executor/config/input/environment/split hashes and the passed certificate hash.

The source engine `pipeline/geo_pipeline.py` remains byte-identical to its public lock, SHA-256 `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`. Checkpoints, fit logs, models and subject-level outputs remain under ignored private directories. Public files contain source code, configuration and aggregate evidence only.
