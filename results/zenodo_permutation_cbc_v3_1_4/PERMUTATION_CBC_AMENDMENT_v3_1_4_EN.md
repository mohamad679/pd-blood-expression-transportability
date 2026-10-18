# Combined amendment: GEO permutation control and CBC date matching

**Release version 3.1.4-permutation-cbc; effective on public publication.** Prepared after the observed GEO results
and before this permutation stage. Parent protocol: v3.1 and implementation
amendments; latest results snapshot: [10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851).
**No registered real-GEO label-permutation replicate has been executed to date.**
Earlier synthetic source-null simulations are separate experiments.

This amendment also records the CBC temporal matching rule for the planned
cell-score validation. It is not effective until this combined amendment is
tested and publicly versioned. The parent CBC rule allows same-day or nearest
measurement within ±30 days; it does not require the same visit code. The frozen
161 baseline controls have zero exact-visit CBC rows, but a metadata-only
all-visit audit found 31 controls with CBC rows in the same calendar month as
their BL blood collection, all at SC, in both exports. Numeric per-marker n is
not established. No CBC correlation has been calculated. Eligibility and marker
panels are unchanged; the previous draft's same-visit requirement is withdrawn.

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
settings. Execution is planned on an x86_64 Linux cloud VM with the exact locked
package versions and Python 3.11.9. Record machine, CPU, OS and environment in
outputs; do not transfer any PPMI participant data. Test a ceiling of eight
workers × two threads using synthetic full-dimensional fixtures before accepting
that ceiling. Linux reads RAM from `/proc/meminfo` and swapout pages from
`/proc/vmstat` (`pswpout`); retain diagnostic `pgpgout` separately. Keep the
75%-RAM conservative budget and zero-swapout-growth gate; missing measurements
cannot pass. macOS counters retain their original meaning. Existing times used
eight threads and cannot be divided by eight. Freeze the tested executor/config hashes in the public amendment before
any real permutation, including a runtime pilot. Use persistent logs, sleep
prevention and resumable outputs; neither protects against shutdown or all sleep
conditions. After publication, use fixed registered replicate 0 at one worker
as a runtime pilot, then reuse it in the full B=200 run. No p-value or seed/B
selection follows from that pilot. Runtime concurrency may be 1, 2, 4 or 8,
within the tested ceiling, based only on timing/resources, without changing the
bound config, labels, fits or statistics. A serial pilot does not establish
parallel throughput. No scoring-method amendment is needed for the closed sigmoid audit.

**CBC participant and temporal matching rule.** Keep the frozen
baseline healthy-control population. Use an available collection month from the
exact participant+BL whole-blood collection record. `BLDCOLL` is not an eligibility
filter or execution gate; observed codes remain uninterpreted. Match the same
participant's CBC rows from any visit, requiring CBC `LCOLLDT` month/year to
equal baseline `BLDCLLDT` month/year. Exclude adjacent months. Two collection
dates within one calendar month are at most 30 days apart, but same-day,
closest-day and earliest-day ordering are unidentifiable, and some true
within-30-day pairs across adjacent months will be omitted.

For each participant and fixed laboratory marker, prefer the RNA visit (BL),
then SC, then other visits. Within the highest available visit-priority tier,
collapse exact duplicate representations of an identified measurement only if
values and units agree. Exclude unresolved multiple collection occasions or
conflicting repeats for that participant-marker; report exclusions and do not
fall back to a lower tier after a preferred-tier conflict. Do not break ties by
CSV order, inferred day, or result magnitude. Use the current export as primary;
the archive provides a concordance check, not additional observations.

The metadata availability upper bound is 31 controls. After publication and
before association analysis, report marker-specific usable numeric n, units,
ambiguity exclusions and date precision. Apply the original per-marker
`n >= 30`, Spearman >= 0.4 and CI lower bound > 0 supportive-evidence criteria;
otherwise report insufficient/non-evaluable evidence. No marker panel or M4/M5
model changes follow from these results. CBC validates specific cell-composition
proxies, not general autosomal or PD-specific transfer, and is not a PD stopping
gate.

This subset is at the evaluability boundary: missing results, incompatible units
or unresolved repeats may reduce individual marker counts below 30. A marker
with fewer than 30 valid paired observations is **non-evaluable**, and neither
the matching rule nor the population, visit priority or time window will be
changed to increase n after marker completeness is inspected. Even a marker
meeting the descriptive evidence rule at approximately 30 participants will be
reported as **limited supportive evidence**, not as definitive biological or
cross-platform validation. The same-month rule may yield fewer eligible pairs
than a true day-resolved ±30-day rule; month-level dates cannot establish how
many adjacent-month pairs would actually satisfy that window. The reason for
the 91 unmatched known-month controls is not inferred from the count alone.

References: [permutation null](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.permutation_test_score.html);
[randomized p-values](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).

## Passing Linux preflight certificate

Evidence freeze: 2026-10-18T13:39:25.938528+00:00. Reserved version DOI: 10.5281/zenodo.23449296.
Author: Mohammad Javad Asgari, ORCID 0009-0002-2313-571X.

All 33 synthetic checks passed across 650 fit attempts; maximum KKT residual
was 2.68391758783e-08. All 20 original-GEO model/fold identity comparisons
passed, with zero AUROC error and a maximum probability error of
7.64943663967e-14.
The eight-worker, two-thread-per-worker fixture observed all eight processes,
used a conservative memory budget of 5,996,040,192 bytes against
32,848,162,816 bytes RAM, and had zero pswpout growth. The fixture's serial and
parallel times were 10.3020 and 2.6809 seconds; these are microbenchmark times,
not forecasts of permutation runtime. The machine was an x86_64 Linux cloud VM
with 16 logical CPUs, AMD EPYC-Genoa CPU model, Ubuntu 26.04.1 and Python 3.11.9.

The exact binding SHA-256 is `4be7bd1bc5e92da5d094fc6f59dbd45df033b254b3dc31043b4af8ed1ccc0310`.
The passing certificate SHA-256 is `2b21d17727adcbe71022a411e4db870b46d616d34d8839da562b0b5b7f80f0ac`.
Code/config hashes and package versions are supplied in the release manifest
and archive. No real randomized GEO replicate, CBC association or PPMI disease
score has been computed at this evidence freeze. The first registered randomized
execution remains replicate 0 at one worker, after public release.

A second attempted controller was rejected while the original preflight was
running; its lock was not manually removed. The original run subsequently
completed and supplied this passing certificate. Transport/authentication/DNS
setup issues did not change scientific code, data, seeds or settings.
