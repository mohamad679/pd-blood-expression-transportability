# Complete English translation with provenance correction — v3.1.2-en

**Author:** Mohammad Javad Asgari. **ORCID:** https://orcid.org/0009-0002-2313-571X.

This edition supplies the full English protocol, numerical plans, solver-amendment
rationale, numerical results, historical warning-hold report, implementation
amendment and GEO pilot report. It is a faithful translation/presentation edition,
**not a methodological amendment, a new preregistration date for the methods,
or a claim that the translation preceded all real-data analysis**.

Original preregistration: https://doi.org/10.5281/zenodo.23405371 .
Pre-execution implementation amendment: https://doi.org/10.5281/zenodo.23405802 .
Unchanged all-versions DOI: https://doi.org/10.5281/zenodo.23405370 .
First English-edition DOI: https://doi.org/10.5281/zenodo.23406064 . Publication timing
is determined by its public Zenodo record, not DOI reservation.

## Timing, source priority and exposure disclosure

The original method and implementation amendment were public before the first
real study-model fit. **One real PROBE/GSE57475 outer fold had already been
completed before this English translation was prepared and published.** Its
five-model predictions, hyperparameter choices and timing results had been viewed.
The pilot report is included in English. The pilot began at 2026-10-15 22:07:32 UTC; this edition and the earlier English edition were published afterward. No additional
real-data training or scoring was performed to produce this translation. No PPMI
disease, sex, age or cell-marker score has been calculated.

The translated bodies preserve the original freeze-time wording. Historical statements such as "public registration has not yet occurred" or
"the real project has not yet run" in translated source documents describe
their original freeze-time status. They are not assertions of the current status.
Their preservation is deliberate; their original timestamps are not backdated.

The specific corrections to English-file pre-execution claims and their original hashes are listed in [PROVENANCE_CORRECTION.md](PROVENANCE_CORRECTION.md). The [difference table](TRANSLATION_DIFFERENCES.md) lists every presentation difference, including these corrections.

## Read in this order

1. [Complete method](protocol/PROTOCOL_v3_1_EN.md), with the unchanged
   [frozen configuration](protocol/config_v3_1.json).
2. [Complete preregistration text](protocol/PREREG_v3_1_FULL_EN.md), alongside
   the existing [English overview](protocol/PREREG_v3_1_EN.md).
3. [Numerical step zero](protocol/KKT_STEP_ZERO_v3_1_EN.md),
   [solver amendment](protocol/KKT_NUMERICAL_AMENDMENT_20261016_EN.md),
   [synthetic plan](protocol/SYNTHETIC_V3_1_PLAN_EN.md) and
   [complete numerical results](protocol/KKT_NUMERICAL_RESULTS_20261016_EN.md).
4. [Implementation amendment](amendments/AMENDMENT_v3_1_1_EN.md) and its
   [translated summary](amendments/AMENDMENT_v3_1_1_SUMMARY_EN.md).
5. [First real GEO fold report](pipeline/GEO_PILOT_REPORT_20261016_EN.md).
6. [Translation differences and precedence](TRANSLATION_DIFFERENCES.md),
   [machine audit](TRANSLATION_AUDIT.json) and [critical parameters](CRITICAL_PARAMETERS.json).

These English texts are translations for reading the registered rules. **In any discrepancy, the original Persian v3.1 record and the published v3.1.1 amendment take precedence.** The existing implementation amendment continues to clarify/extend the parent method exactly as before. No translation grants
permission to change a model, population, threshold, gate, endpoint, grid or
inference rule. If a substantive translation discrepancy is discovered, disclose
and correct it publicly before the affected step; do not silently choose a
different interpretation after results. The original timestamped sources, not these later English filenames, establish what was registered and hashed when.

## What has and has not changed

All original Python/R code, all original JSON/JSONL data, the dependency lock,
gene definitions and the frozen config are byte-identical to the original source
snapshot. Thus every config key/value, seed, threshold and grid is unchanged.
Original config operational flags and its historical author spelling are preserved
as frozen evidence; current citation metadata are in `CREATORS.json` and Zenodo.

Every previously Persian document in the original ZIP has a complete English
counterpart, including historical copies. Markdown links point to English files.
The frozen config, original hash manifests and historical JSON logs retain some
original source filenames containing `_FA`: those are English-text provenance
strings, not embedded Persian documents. `TRANSLATION_SOURCE_MAP.json` resolves
each original source document to its English counterpart and records both hashes.
The original source release hash manifest is under `provenance/`; it is not a hash
manifest for the translated files. Use `ENGLISH_RELEASE_SHA256.json` for this edition.

`TRANSLATION_AUDIT.json` compares numeric literal inventories after normalizing
Persian digits, decimal/group separators and scientific notation. It also checks
exact config/code preservation and recursively scans the whole package for Arabic
script and filenames marked `_FA`. Digits embedded in identifiers are not treated
as quantities; technical identifiers are checked separately. Matching numbers do
not prove semantic equivalence; the prose and decision rules were also reviewed
against the source. The build/audit script is supplied for inspection.

## Data boundary and remaining work

No individual PPMI metadata, counts, private crosswalk, cache, model predictions,
credentials or controlled-access documents are included. GEO split identifiers
come from public sources. Pilot outputs supplied here are aggregate reports and
hash commitments, not individual expression data.

The full auxiliary-control, sensitivity, bootstrap-coverage, permutation and
external-evaluation tasks are still outstanding as disclosed in the pilot report.
This translation is not the required second public lock of final models,
preprocessing and the final PPMI population. That lock and all existing transfer
gates remain mandatory before PPMI disease scoring. Future project documents,
code comments, logs and public releases will be written in English.
