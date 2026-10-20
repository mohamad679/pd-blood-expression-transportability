# Decision before calculation: descriptive GEO transfer uncertainty

**Status:** Fixed before running `pipeline/bootstrap_geo_transfer_descriptive.py`. This is a post-development, descriptive analysis of the two already completed GEO-to-GEO transfers. It is not the confirmatory PPMI bootstrap, does not change any model or select a favorable result, and computes no new fitted prediction.

## Inputs and checks

Use only the saved original (uncorrected) held-out prediction files for `GSE57475/transfer` (PROBE → GENEPARK) and `GSE99039/transfer` (GENEPARK → PROBE). For each file, read the stored binary outcome and M1–M5 probabilities; map its sample IDs to the prepared target-cohort metadata solely to verify one prediction per unique operational donor and agreement of labels. Require exact target-cohort sample coverage, finite probabilities in [0,1], original metadata hashes from `data/manifests/input_audit.json`, and observed AUROC for every model matching the saved original transfer `RESULTS.json` within 1e-12. Stop rather than repair any mismatch. Do not read or score PPMI.

## Frozen resampling rule

For each direction independently, generate **5,000** donor-level bootstrap resamples. Resample with replacement *separately* within the PD and control strata, preserving their original sizes. Apply each resample's identical donor multiplicities to all five saved model scores. Compute tied-score AUROC for each M1–M5 and the paired difference `AUROC(M5) − AUROC(M4)` in every resample. Report each observed value and the two-sided **95% percentile interval** (empirical 2.5th and 97.5th quantiles). Use NumPy `default_rng` with fixed seed **20261020** for PROBE → GENEPARK and **20261021** for GENEPARK → PROBE. No p-value, threshold optimization, model refit, model selection, or multiplicity claim is added.

The intervals condition on the two fitted transfer models and these observed held-out cohorts; they do not include model-development, cohort-selection, platform, or between-study uncertainty. The sampling design is paired across models and stratified by outcome, matching the project's existing fixed-score bootstrap logic. They are **descriptive GEO transfer intervals**, not confirmatory PPMI evidence.

## Output and provenance

Write one aggregate-only JSON file under `results/` with direction, sample/class counts, observed AUROCs, intervals, observed M5−M4 and its interval, input file hashes, method, seeds, repetitions, decision-document hash, code hash and Git commit SHA. Do not export participant IDs, per-person scores, or bootstrap draws. Before this calculation, commit and push this decision and the exact executable script. If Git provenance or any data check fails, stop without calculating intervals.

## Pre-calculation implementation correction

The first committed invocation stopped with `KeyError: 'input_sha256'` while looking for metadata hashes in `results/geo_development_20261018/aggregate_summary.json`; that summary has no such field. It stopped before creating resampling draws or an output file. The hashes were already frozen in `data/manifests/input_audit.json`; only this lookup path was corrected, and this correction is committed and pushed before calculation. The resampling rule, seeds, metrics, and reporting plan above are unchanged.
