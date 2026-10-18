# Parkinson's disease blood transcriptome: shortcut audit and transportability

This repository contains the English-language research record and executable GEO development pipeline for a study of how blood transcriptomic PD classifiers behave when they are evaluated across cohorts and platforms. The study asks how much of apparent within-cohort performance survives strict cohort-held-out testing and external transfer, and whether age, sex, platform, or cell-composition proxies explain part of the signal.

**Current stage:** the protocol and primary M1–M5 engine are implemented; synthetic checks passed; all 13 prespecified GEO nested-CV, directional-transfer, and pooled-fit scopes have completed. Aggregate results are in [`docs/GEO_CV_RESULTS_20261018_EN.md`](docs/GEO_CV_RESULTS_20261018_EN.md) and [`results/geo_development_20261018/aggregate_summary.json`](results/geo_development_20261018/aggregate_summary.json). Registered GEO controls and diagnostics remain incomplete. No PPMI model score has been computed. This repository is a working research project, not a completed paper or a claim of validated PD biomarker performance.

The full project history and remaining roadmap are in [`docs/PROJECT_HISTORY_AND_ROADMAP.md`](docs/PROJECT_HISTORY_AND_ROADMAP.md). Start with [`docs/STATUS.md`](docs/STATUS.md) for the current execution state.

## Repository layout

- `pipeline/`: GEO loading, preprocessing, donor-level split execution, M1–M5 training, model scoring, paired bootstrap utilities, synthetic checks, and preparation utility.
- `protocol/`: fixed numerical solver, block-penalty implementation, English protocol, preregistration translation, amendments, lock configuration, dependency lock, and provenance notes.
- `data/manifests/`: prespecified GEO split manifest, aggregate input and mapping audits, and shared gene manifest. No expression matrices or individual-level PPMI data are stored here.
- `results/geo_pilot/`: aggregate metrics, timing, public split-membership audit, hashes, and the English pilot report. Individual held-out predictions, fitted coefficients, preprocessing arrays, and per-fit logs are excluded.
- `results/geo_development_20261018/`: aggregate-only machine-readable GEO CV and transfer summary. Per-sample predictions, fitted coefficients, preprocessing arrays, and per-fit logs are excluded.
- `docs/`: project history, status, data boundaries, and next-step roadmap.

## Reproducibility and data boundaries

The primary pipeline consumes two prepared GEO development matrices (`GSE57475_development.npz` and `GSE99039_development.npz`) and corresponding sample metadata under `prepared/`. These local inputs are intentionally excluded from Git; see [`data/README.md`](data/README.md) for their required schema and privacy boundary. The pipeline has no PPMI loader or PPMI scoring command.

Use Python 3.11 with the pinned dependencies in [`requirements_glum_v3.lock.txt`](requirements_glum_v3.lock.txt). The registered implementation, configuration, and dependency versions are preserved in `protocol/`; do not edit those files in-place to silently change a locked analysis. Any method change that could affect confirmatory results requires a timestamped amendment before the affected analysis is run.

The checked-in split manifest and public preregistration must agree before execution. Existing real-data runs also require the locally held public-registration receipt with matching code and split hashes. Do not fabricate or bypass that receipt. Use a fresh output directory for every run.

## Current protocol and public record

- Parent preregistration: [Zenodo DOI 10.5281/zenodo.23405371](https://doi.org/10.5281/zenodo.23405371)
- Implementation amendment 3.1.1: [Zenodo DOI 10.5281/zenodo.23405802](https://doi.org/10.5281/zenodo.23405802)
- English provenance-corrected reference package: [Zenodo DOI 10.5281/zenodo.23411084](https://doi.org/10.5281/zenodo.23411084)
- Completed core GEO CV/transfer milestone, version 3.1.3-geo-cv: [Zenodo DOI 10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851). This is a post-results snapshot; remaining diagnostics and the final pre-PPMI lock are outstanding.

The English protocol files are translations and reading aids. Their provenance notes identify their source and timing; if the translation conflicts with the source preregistration or a published amendment, the published source documents govern. The GitHub repository is the living code and work log. Per the project decision, Zenodo will be updated only at major milestones: after the full GEO CV package is complete, and again before any PPMI scoring under the required final lock.

## Main models

All models use the same locked weighted logistic-regression family. M1 is age/sex; M2 is genes only; M3 adds genes to M1; M4 is age/sex plus six cell-marker panels; M5 adds genes to M4. Baseline terms are unpenalized in the nested gene models. The primary comparison is the incremental AUROC of M5 over M4 in the external PPMI population, conditional on prespecified transfer and QC gates. The protocol limits interpretation to predictive performance of this procedure; it does not estimate causal percentages attributable to disease versus shortcuts.

## License and use

This repository is private. No open-source license has been granted. PPMI data-use terms prohibit unauthorized sharing, public database posting, and re-identification attempts; all controlled PPMI material and individual-level records must remain outside GitHub and other public repositories.
