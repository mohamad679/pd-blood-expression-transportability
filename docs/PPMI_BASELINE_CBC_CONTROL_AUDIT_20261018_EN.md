# Baseline-control CBC availability audit

**Date:** 2026-10-18  
**Purpose:** Check whether the frozen baseline healthy-control set has a usable CBC measurement at its RNA-seq visit before calculating any cell-score/CBC association.  
**Outcome data and expression:** Not read. No scores, correlations, or PPMI disease predictions were calculated.

## Frozen candidate set

The audit reproduces the previously used eligibility definition: 2021 RNA metadata marked Healthy Control, baseline `BL` / `M0`, Long RNA-seq `Yes`, `HighQC`; exact Healthy Control cohort in the 2026 participant-status export; one valid baseline age; recorded Male/Female sex; and both portal diagnosis fields equal to Control. Matching keys were exact normalized participant plus visit (`PATNO` + `EVENT_ID`), without fuzzy or adjacent-visit matching.

The resulting set contains **161** eligible baseline control sample-visits.

## Aggregate availability results

| Check | Count |
|---|---:|
| Eligible baseline controls | 161 |
| Whole-blood collection record at exact visit | 126 |
| Collection record with reported month/year | 122 |
| Current CBC export: any row at exact visit | 0 |
| Archived CBC export: any row at exact visit | 0 |
| Current or archived CBC: numeric core marker at exact visit and same reported month | 0 |
| Same-month absolute neutrophil and absolute lymphocyte pair | 0 |

Across the 126 exact-visit collection records, the observed `BLDCOLL` code counts are 4 records with code `0` and 122 with code `1`; 122 of the 126 records have a reported month/year. The aggregate audit does not cross-tabulate code against date availability. These code values remain **uninterpreted** until matched to an authoritative PPMI code list. No meaning is inferred from their numeric values.

## Temporal rule for any future association

The proposed rule is: require exact participant and visit agreement, then retain only CBC measurements whose reported collection month/year equals the whole-blood collection month/year. Exclude adjacent months. Given month-level dates, two dates in the same calendar month are at most 30 days apart, but the available fields cannot establish the exact day difference, same-day status, or which of multiple records is nearest. This rule is to be included in the combined, not-yet-published amendment before any CBC association is computed.

Under this rule the baseline-control CBC validation is **not evaluable**: there are zero eligible exact-visit CBC rows, below the protocol's `n >= 30` minimum. Do not expand the visit window or replace the control set after seeing this result. Any alternative CBC validation population would need a separately justified, pre-outcome rule and an amendment before calculating associations.

## Interpretation

CBC can validate whether selected cell-composition proxy scores agree with measured blood counts. Its absence here does not establish that all autosomal expression signal fails to transfer, and a positive CBC association would not by itself prove transfer of PD-specific autosomal signal.

## Reproducibility and privacy

The audit was reproduced from the local PPMI exports using `pipeline/audit_ppmi_baseline_control_cbc.py`. The aggregate JSON is stored under ignored `private_data/`; the committed report contains only counts and no participant identifiers, dates, row-level records, or protected source data.

The authoritative PPMI [Guidance Resources page](https://www.ppmi-info.org/access-data-specimens/guidance-resources) states that its Data Dictionary and Code List provide variable/schema definitions and values. A usable official code-list file was not obtained during this audit; `BLDCOLL` meaning therefore remains unresolved.
