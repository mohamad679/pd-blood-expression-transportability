# Retrospective audit of the six GEO sigmoid warnings

**Completed:** 2026-10-18, 05:27:54 UTC. **PPMI read or scored:** no.

All six previously observed warnings were reproduced and attributed. Explicit
sigmoid clipping changed **zero validation AUROCs and zero hyperparameter
choices** in the audited contexts. Every original candidate mean was reproduced
exactly. This bounded diagnostic found no reason to amend the scoring method or
replace the completed scientific GEO results. The locked engine and original
artifacts remain unchanged.

Exact evidence is in the [aggregate audit](../results/geo_controls_20261018/sigmoid_audit.json).
The [scope clarification](SIGMOID_REPLAY_SCOPE_20261018_EN.md) discloses the
preliminary driver error and interrupted diagnostic attempt. The completed
[sex/age controls](GEO_TRANSFER_CONTROLS_20261018_EN.md) are a separate workstream.

## Mechanism and why the warning required investigation

Inspection of the installed, pinned glum 3.4.1 implementation showed that
`GeneralizedLinearRegressor.predict` calls `LogitLink.inverse`. This computes
scipy `expit` and clips values outside `[50*eps, 1-50*eps]`. For float64,
`50*eps` is approximately 1.1102e-14. It emits the warning:

> Sigmoid function too close to 0 or 1. Clipping.

The locked GEO pipeline captures warnings during fitting and its KKT check, but
inner-validation `fit.predict()` calls occur outside that capture. Therefore,
empty original fit-warning fields did not establish that prediction warnings
were absent. Passed KKT checks also could not establish that clipping preserved
AUROC: clipping can merge differently labelled samples into a tie. The synthetic
check constructed exactly such an example, in which AUROC changed from 0.75 to
0.50. Warnings were investigated rather than dismissed based on convergence.

## Bounded replay

All recorded candidates and five inner folds were replayed for the six reported
model/scope contexts, together with the required M1/M4 no-gene baselines:

| Scope | Audited gene models | Inner fits, including required baselines |
|---|---|---:|
| GSE99039/outer/0 | M5 | 80 |
| GSE99039/outer/2 | M3 | 95 |
| GSE99039/transfer | M5 | 80 |
| pooled/final | M2, M3, M5 | 235 |
| Total | Six model/scope contexts | 490 |

The replay reconstructed partitions from the frozen split manifest and verified
each training preprocessing fingerprint against the original fit record. It used
the unchanged cold-fit wrapper, weights, feature selection, penalties, solvers,
KKT threshold and retry rules. No alternate grid was searched. For every fit it
captured prediction warnings and compared AUROC from glum's clipped probability,
unclipped expit probability and the linear predictor. The original recorded
candidate grids were used for all three diagnostic selection comparisons.

This is a retrospective numerical diagnostic, not a new preregistration, a
replacement scientific analysis or a rerun of every fit in all 13 GEO scopes.
Its findings apply to the six observed warning contexts and their complete
candidate sets. No external PPMI property is established by it.

## Exact warning attribution

Fold indices below are **zero-based**. All six events were in inner validation,
each affecting exactly **one upper-tail probability**; none affected a lower-tail
probability. No event merged a PD/control pair at the clipped tail. Each row's
AUROC was identical for clipped probability, unclipped expit and linear score.

| Scope | Model | Fold index | Lambda | Rho | Validation n | Unchanged AUROC |
|---|---|---:|---:|---:|---:|---:|
| GSE99039/outer/0 | M5 | 4 | 0.0001 | 1 | 38 | 0.614958449 |
| GSE99039/outer/2 | M3 | 3 | 0.00001 | 1 | 38 | 0.770083102 |
| GSE99039/transfer | M5 | 3 | 0.0001 | 1 | 48 | 0.534722222 |
| pooled/final | M2 | 2 | 0.0001 | 1 | 76 | 0.695586623 |
| pooled/final | M3 | 2 | 0.0001 | 1 | 76 | 0.730902778 |
| pooled/final | M5 | 4 | 0.0001 | 1 | 76 | 0.683680556 |

Rho is the L1 ratio; rho=1 denotes an L1/lasso gene penalty. None of these six
warning-producing candidates was the selected candidate in its model/scope.

## Reproduction and selection results

- All 490 replay fits passed on the first attempt; no captured fit warnings,
  exceptions or retries. Largest KKT residual: **2.2577973759693726e-8**, below
  the fixed 1e-6 limit.
- Maximum absolute difference between replayed and original candidate mean
  AUROC: **0.0**, including the required baseline means.
- All six original selections were reproduced. Unclipped-expit selection
  changes: **0**. Linear-score selection changes: **0**.
- Explicit clipping changed no per-fold or candidate-mean AUROC in any of the
  490 audited validation calls.
- Actual replay wall time: **2,584.24 seconds (43.07 minutes)**. The earlier
  28.5-minute estimate was the sum of the recorded original fit times, not a
  guaranteed runtime for this diagnostic.

| Model/scope | Unchanged selected lambda | Unchanged rho | Original and replayed inner mean AUROC |
|---|---:|---:|---:|
| GSE99039/outer/0, M5 | 0.1 | 0 | 0.733462604 |
| GSE99039/outer/2, M3 | 0.0001 | 0.5 | 0.699335180 |
| GSE99039/transfer, M5 | 0.1 | 0 | 0.715775966 |
| pooled/final, M2 | 0.1 | 0 | 0.675819429 |
| pooled/final, M3 | 0.01 | 1 | 0.669489973 |
| pooled/final, M5 | 0.001 | 1 | 0.673884678 |

### Separate finite-precision effect

Three validation calls **without clipping** had slightly different AUROC for a
linear predictor versus its expit probability. This is a floating-point tie
effect, distinct from glum's explicit tail clipping. It occurred in two
nonselected M3 candidates in GSE99039/outer/2 and one nonselected pooled M3
candidate. Maximum candidate-mean differences were 0.000263158 and 0.000173611,
respectively. No selection changed. Full aggregate contexts are preserved in
the [rounding diagnostic](../results/geo_controls_20261018/rounding_diagnostic.json).
The scoring representation was not switched after inspecting these results.

## Preservation, provenance and consequence

The core pipeline SHA-256 remains
`6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`, and the split
manifest SHA-256 remains
`2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.
The completed replay driver SHA-256 is
`d3339598c3b506d34cfcfda5ce0a35c50204dfedfe72229a8dd380b876af0732`.
The installed glum link-file hash is recorded in the aggregate.

All 167 post-pilot original execution files match the pre-existing backup
byte-for-byte. Thirteen pilot artifacts match the original pilot hash manifest.
The [preservation record](../results/geo_controls_20261018/artifact_preservation.json)
contains the verification result. Original scientific model arrays, predictions
and detailed logs were not replaced or committed. Detailed diagnostic logs live
only in ignored `private_data/sigmoid_audit_20261018_targeted/`.

For completeness, the first replay-driver attempt stopped before any fit due to
a duplicate context keyword. It was corrected in commit `eb7b9c0`. A subsequent
all-scope diagnostic was interrupted after 111 completed fits to bound this
investigation to the six reported warning contexts; its artifacts were retained.
The bounded scope was documented before the successful 490-fit replay. Neither
preliminary attempt changed or invalidated an original scientific GEO run.

**Decision:** close the six-warning review without changing the registered
pipeline or repeating the core scientific analyses. Preserve the observed
warnings and their attribution. No method amendment or new Zenodo release is
needed for this diagnostic finding. Remaining leakage checks, bootstrap coverage,
bounded sensitivities, permutation decisions, PPMI linkage/eligibility and the
second public lock still precede external scoring. No PPMI score has occurred.

## Reproduction on macOS

With the frozen prepared inputs and original local GEO run artifacts available,
run from the repository root using fresh output paths:

```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
UV_CACHE_DIR=/private/tmp/pdgeo_uv_cache \
uv run --offline \
  --default-index https://pypi.tuna.tsinghua.edu.cn/simple \
  --with-requirements requirements_glum_v3.lock.txt \
  python pipeline/audit_sigmoid.py --reported-contexts-only \
  --runs private_data/geo_runs_20261018 \
  --pilot ../pd_shortcut_feasibility_2026-10-13/pipeline/pilot_PROBE_fold0_20261016 \
  --output private_data/sigmoid_audit_NEW_RUN \
  --aggregate /private/tmp/sigmoid_audit_NEW_RUN.json
```

Offline execution requires cached pinned dependencies. This command audits
existing GEO results; it neither reads nor scores PPMI.
