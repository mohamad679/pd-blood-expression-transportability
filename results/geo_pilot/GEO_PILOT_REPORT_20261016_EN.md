# Pipeline implementation and first real GEO fold report

> **Provenance and status.** This English text was prepared after the first real PROBE fold began at 2026-10-15 22:07:32 UTC. Its source is the post-pilot Persian report, prepared after the pilot; it was not preregistered: `GEO_PILOT_REPORT_20261016_FA.md` (SHA-256 `1faa7aa5519ac2f9b37e7cb94c6c2970ab488e9e165488919aa83f0592772d8c`). It introduces no method change. If wording conflicts, the original Persian preregistration and its published amendment govern. Statements of registration, execution, and publication status in the translated body describe its original snapshot, not the reader's current date. See [TRANSLATION_SOURCE_MAP.json](../TRANSLATION_SOURCE_MAP.json).

Date: 16 October 2026. This is the first execution of real study models; earlier
work comprised preparation/feasibility and synthetic tests. No disease, sex, age
or cell-marker score has been calculated in PPMI.

## Public registration and author name

The original record's name was corrected to **Asgari, Mohammad Javad**, and ORCID
**0009-0002-2313-571X** was added. PPMI linkage wording was limited to computational
identifier linkage reported by the local audit; verified biological donor identity
and a finalized population were not claimed.

Amendment 3.1.1 was made public in the same Zenodo lineage before the first real fit:
https://doi.org/10.5281/zenodo.23405802 . Original version 3.1 and its files remain:
https://doi.org/10.5281/zenodo.23405371 .

Amendment public verification: `2026-10-15T22:07:25.879272+00:00`; pilot start:
`2026-10-15T22:07:32Z`. Both correspond to 16 October in Tehran.
The four site checksums matched the local files. Code, input and split hashes are
in the amendment; this implementation lock does not replace the second final-model/
PPMI-population lock.

## Implemented and checked

- Primary M1–M5 training/tuning/prediction engine with the locked objective,
  weighting, block penalties, grid, tie rule and KKT rules.
- Donor splits for ten outer folds, two directional transfers and pooled final
  fitting; MAD selection, standardization and basis construction use training
  data only inside each fit.
- Paired donor bootstrap for fixed scores, including ties and class
  stratification. CI coverage was not checked here and remains required before PPMI.
- The 19 synthetic checks passed, covering repeated donors, held-out data/label
  perturbations, M4/M5 baseline equality, the baseline branch and independent
  verification of bootstrap AUROC calculation.
- Audit of all **269** real fit logs found training membership consistent with
  the public splits; the outer test entered no inner fit. Outputs match the locked
  test membership.

PROBE GSM is the operational donor identifier; this differs from independent
proof of no duplicates or cross-source overlap, and that limitation remains.

## First-fold timing and convergence

Only **GSE57475/outer/0** was executed: 113 training records, 29 test records,
1000 genes selected in the final fit, five inner folds and all five models.
Predictions and all attempts were preserved.

| Quantity | Result |
|---|---:|
| Total elapsed time | 427.34 seconds, about 7 minutes and 7 seconds |
| Process CPU, summed across threads | 1884.47 seconds, about 31.4 minutes |
| Fit/attempt count | 269 |
| KKT failures / retries / warnings | Zero / zero / zero |
| Maximum KKT residual | 3.8924e-8, below the fixed 1e-6 limit |
| Slowest fit | 19.25 seconds |

Model fits were not concurrent; internal libraries defaulted to up to 8
threads. Wall time must therefore not be called single-core CPU time. No warm
start, new model parallelism or solver change from the locked version was applied.

| Model | Fits | Total fit seconds | Development selection in this fold |
|---|---:|---:|---|
| M1 | 6 | 0.09 | Age/sex baseline |
| M2 | 106 | 258.26 | λ=1e-5, ρ=0 |
| M3 | 75 | 82.15 | Baseline-only; exact reuse of M1 |
| M4 | 6 | 0.12 | Age/sex and six-panel baseline |
| M5 | 76 | 82.77 | λ=0.01, ρ=0.5 |

M2 expanded the grid twice under the registered rule; final selection was
interior. No selected model remained at a boundary. M3 and M1 had bit-identical
predictions. This branch only shows that the gene block was not selected in this
particular tuning procedure.

## Observed results, with explicit limitations

Test AUROCs in 29 people: M1=0.3658, M2=0.5105, M3=0.3658, M4=0.3947 and
M5=0.4895. These are retained for transparency; they are neither final study
performance estimates nor evidence of successful/failed PPMI transfer. This fold's
M5−M4 difference is not a confirmatory result. The fold belongs to the original
primary CV and will not later be counted again as independent testing. Neither
method nor scope was changed on the basis of its performance.

## Budget: the core analysis is practical; permutations dominate

The 13 primary analyses require at most 3120 initial fits, or 4290 with maximum
grid expansion. If all analyses took exactly as long as this fold, the core cost
would be about **1.54 wall hours**, or **7.72 hours** with four complete sensitivity
reruns. This is simple multiplication only; GENEPARK and pooled fitting have
larger samples and their correlation/convergence behavior has not been measured.
This is not a firm delivery time.

By contrast, 1000 permutations with full retuning may require hundreds of
thousands to several million fits. Repeating this fold's time for five models
and ten outer analyses gives a simple scenario of **49.5 wall days**; restricting
to M4/M5 gives about **9.6 days** for the two models' fits alone, excluding
preprocessing and other overhead. These are not validated forecasts or a chosen
permutation scope; they illustrate the magnitude. Silently reducing permutation
count, omitting y-dependent refitting or removing slow candidates is not allowed.

Bulk execution was not started in this turn so the secondary control's scope and
cost could be reported clearly. Any required acceleration must be documented
with synthetic testing and a public amendment before affected execution. The core
analysis did not stop because of a numerical problem; the pilot passed completely.

## Remaining work

Complete the other GEO analyses, sex/age controls, sensitivities, registered
diagnostics, bootstrap-coverage simulation, and precise permutation compute plan
and its execution; then PPMI controls and final population determination under
the protocol. The primary engine is implemented, but these tasks must not be
represented as complete.

Before the first PPMI disease score, a new public version must hash code/environment,
splits, final M1–M5 models with all fitted preprocessing states, and the final
private population/crosswalk. Individual data are not published. This second lock
has not yet occurred.

Local evidence: `pilot_PROBE_fold0_20261016/TIMING_BUDGET.json`, `fits.jsonl`,
`REAL_FOLD_MEMBERSHIP_AUDIT.json`, `RESULTS.json` and `heldout_predictions.npz`.
