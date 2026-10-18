# Proposed amendment: bounded GEO label-permutation negative control

**DRAFT — not published or effective.** Prepared after the observed GEO results
and before this permutation stage. Parent protocol: v3.1 and implementation
amendments; latest results snapshot: [10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851).
**No registered real-GEO label-permutation replicate has been executed to date.**
Earlier synthetic source-null simulations are separate experiments.

This draft also records a proposed CBC temporal matching rule for the planned
cell-score validation. It is not effective until this combined amendment is
tested and publicly versioned. The frozen set of 161 baseline healthy controls
currently has zero CBC records at the exact RNA-seq participant+visit in either
the current or archived CBC export. Consequently the baseline-control CBC
association is non-evaluable under the registered minimum `n >= 30`; no CBC
correlation has been calculated. The rule below does not authorize broadening
the population after this availability result.

**Change and rationale.** Replace 1,000 repetitions across M1–M5 with **200 fixed
repetitions of M4/M5**, covering only the ten within-source outer folds: five
per GEO. Exclude directional transfers and pooled final fitting. Original
ten-fold wall time was 70.92 minutes; recorded M4/M5 fit time was 18.38 minutes
(M5: 26.10% of all model-fit time). Linear scaling gives 61.27 serial hours for
200 repetitions, excluding overhead and potentially different null tuning costs.
The reduction is computational, selected before permutation outputs. This remains
a secondary exploratory negative control; it changes no PPMI endpoint or model.

**Randomization and refitting.** For each repetition, make one coherent donor-level
label map, shuffling within cohort × fixed outer-test-fold membership. Use that
same map across all five outer analyses; permute both training and evaluation
labels. Preserve outer membership and block class totals. Rebuild stratified
inner donor folds from the permuted training labels using the locked splitter
and scope-seed rule; do not retain inner strata constructed from original labels.
Apply this same split-generation procedure to the observed-label reference and
verify identity reproduction before execution.

Recompute class/source weights, tuning, allowed lambda expansions and final
fits. Retain training-only scaling, the positive-MAD rule, tie rules, exact M4
reuse for a zero/no-gene M5, KKT <= 1e-6 and the one fixed cold retry. MAD is
label-free, but selected genes/scaling must correspond to each current training
partition; caches require matching training-membership and implementation hashes.
Capture prediction warnings. Never drop a failed candidate, fold or replicate,
replace its seed, or calculate a completed-test p-value from a shortened denominator.

**Statistics and interpretation.** For each cohort separately, use the unweighted
mean of its five outer-fold AUROCs for M4, M5 and paired M5-minus-M4. Evaluate
the same statistic with observed and permuted labels. For each statistic report
`p = (1 + count[T_b >= T_observed - 1e-12]) / 201`; the fixed tolerance includes
numerical ties conservatively. Report null summaries and Monte Carlo uncertainty.
Minimum p is 1/201, approximately 0.00498; near a tail probability of 0.05,
Monte Carlo SE is approximately 0.015. These six descriptive p-values are not
confirmatory decisions. Global label shuffling tests a global association null;
its delta distribution is **not** a conditional test of added gene information
given the baseline, and chance-like results do not certify absence of leakage.

**Reproducibility and execution gate.** Use IDs 0–199, both cohorts per ID.
For each cohort, derive the RNG seed from the first four big-endian bytes of
SHA256(`20261014:permutation/{rep_id}/{cohort}`); sort donors before shuffling.
Seeds are independent of worker count, scheduling and restart. Preserve the
existing inner scope seeds. Atomic per-repetition completion records must verify
code/config/input/label-map hashes and all required successful fits before reuse.
Fix B in advance; no significance-driven early stopping. Report incomplete or
invalid dependent analyses explicitly.

Before publication, validate membership isolation, identity behavior, synthetic
null behavior, resume integrity and numerical equivalence across worker/thread
settings. Benchmark peak memory, swapping and throughput before accepting four
workers × two threads; existing times used eight threads and cannot be divided
by four. Freeze the tested executor/config hashes in the public amendment before
any real permutation, including a runtime pilot. Use persistent logs, sleep
prevention and resumable outputs; neither protects against shutdown or all sleep
conditions. No scoring-method amendment is needed for the closed sigmoid audit.

**CBC participant/visit and temporal matching rule (proposed).** A CBC row is
eligible only when normalized participant ID and visit exactly match the RNA-seq
sample and the CBC `LCOLLDT` month/year equals the whole-blood collection
`BLDCLLDT` month/year for that same visit. Exclude adjacent months and do not
choose among records using an inferred nearest day. Two dates known only to fall
within the same calendar month are at most 30 days apart, but these exports do
not reveal same-day status or exact day distance. Validate specific cell-score
proxies against numeric CBC markers only if the prespecified minimum sample
size is met; otherwise report the association as non-evaluable. Do not calculate
correlations on the currently available baseline-control set (`n=0` exact-visit
CBC rows). `BLDCOLL` values remain uninterpreted until verified against the
applicable official PPMI code list. CBC associations are a marker-proxy check,
not a general test of autosomal or PD-specific transfer.

References: [permutation null](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.permutation_test_score.html);
[randomized p-values](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).
