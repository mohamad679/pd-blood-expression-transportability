# Reporting correction: final M5 penalty name

**Date:** 18 October 2026. **Type:** descriptive correction; no methodological change or re-estimation.

The final-pooled-fit paragraph of `GEO_CV_RESULTS_20261018_EN.md`, preserved in the Zenodo version [3.1.3-geo-cv](https://doi.org/10.5281/zenodo.23445851), incorrectly calls the selected M5 penalty "ridge-only" while correctly reporting `lambda=0.001, rho=1.0`.

In the registered objective and implementation, `rho` is the L1 mixing fraction (`glum`'s `l1_ratio`):

- `rho=0`: pure L2 / ridge gene penalty;
- `rho=0.5`: mixed L1/L2 / elastic net;
- `rho=1`: pure L1 / lasso gene penalty.

The correct description of final M5 is **L1/lasso on the gene block, lambda=0.001, rho=1.0, with an unpenalized baseline**. The final fit has 210 nonzero gene coefficients among 1,000 selected genes. Final M3 is likewise L1/lasso (`lambda=0.01, rho=1.0`); final M2 is L2/ridge (`lambda=0.1, rho=0.0`).

The source is the unchanged `solve_once` implementation in `protocol/numerical_kkt_v3_1.py`, which passes `l1_ratio=rho`, together with the preserved pooled-final decisions and fit logs. The code, configuration, saved models, hyperparameter values and performance numbers are unchanged. This note preserves the released report as a historical snapshot rather than silently rewriting its published bytes. It should accompany subsequent manuscript/reporting work and be included in the next scheduled public milestone.
