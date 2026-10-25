# Synthetic paired-bootstrap coverage check

**Date:** 2026-10-19 (Asia/Tehran)

## Scope

This is a synthetic fixed-score check of the registered paired, class-stratified
donor percentile interval for `AUROC(M5) - AUROC(M4)`. It reads no participant
files, fits no models, and uses no PPMI scores. The provisional planning class
sizes (161 controls, 311 cases) are used only as synthetic sample sizes. Each
outer sample uses the registered 5,000 bootstrap resamples and fixed bootstrap
seed. The simulation seed is 20261019.

The coverage experiment uses 200 independent outer samples in each of four
scenarios: continuous distinct scores under a true zero difference; continuous
distinct scores with a true difference of 0.05; distinct three-level scores
with ties and a true zero difference; and the identical-model branch. The first
three use correlated fixed score-generating distributions (correlation 0.8).
Known population AUC differences are computed from the generating distributions.

## Results

| Scenario | True ΔAUROC | Covered / 200 | Coverage | 95% Monte Carlo interval | Mean CI width |
|---|---:|---:|---:|---:|---:|
| Continuous, distinct, null | 0 | 191 | 95.5% | 91.7–97.6% | 0.0657 |
| Continuous, distinct, Δ=0.05 | 0.05 | 191 | 95.5% | 91.7–97.6% | 0.0649 |
| Tied, distinct, null | 0 | 192 | 96.0% | 92.3–98.0% | 0.0783 |
| Identical model | 0 | 200 | 100.0% | 98.1–100.0% | 0 |

The 100% identical-model result is a degenerate branch: all 200 intervals have
zero width and equal `[0, 0]`. It is a consistency check, not evidence about
coverage for distinct models. In the three nondegenerate scenarios, the
Monte Carlo intervals include nominal 95% coverage. With 200 outer samples,
these intervals remain broad; this is a screening result and does not establish
high-precision coverage calibration.

## Implementation check and disposition

The vectorized calculation was compared with the published bootstrap
implementation on 5,000 resamples for both untied and tied synthetic scores.
All resample differences were exactly zero, the point estimates matched, and
the percentile intervals matched. The vectorized code changes only the way
tie groups are accumulated; the resampling law, seed, interval definition and
number of bootstrap repetitions remain fixed.

No material coverage failure was observed in these four synthetic scenarios,
so this check alone does not trigger a protocol amendment or a change of CI
method. The finite outer-simulation precision should be kept visible in the
interpretation. This is not model-training uncertainty, clinical validation,
or a guarantee for all score distributions.

Reproducibility artifacts:

- `pipeline/simulate_bootstrap_coverage.py`
- `results/bootstrap_coverage_20261019.json`
