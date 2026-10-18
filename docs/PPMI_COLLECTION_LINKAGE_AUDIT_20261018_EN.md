# PPMI Whole-Blood Collection Crosswalk Audit

**Date:** 18 October 2026  
**Scope:** Exact visit-level linkage and date-precision audit only.

## Purpose and privacy

This audit checks whether the archived `Whole_Blood_Sample_Collection` export contains records at the same participant-visit keys as the Long RNA-seq sample metadata, and whether the collection and CBC exports contain enough date precision to assess temporal alignment. The code reads the restricted files locally and writes aggregate counts only. It does not inspect disease labels or gene expression, and it computes no scores, correlations, or disease associations. Raw participant IDs, row-level dates, and restricted files are not committed.

## Method

The RNA metadata uses the repository's established parser for `PATNO` plus the visit token embedded in `PATNO Visit`. The archived collection export is linked by exact normalized `PATNO` plus `EVENT_ID`; no fuzzy or date-based linkage is attempted. CBC rows are restricted to the nine prespecified numeric core hematology markers and matched to those exact visits. Date comparison is limited to the reported calendar month.

## Aggregate results

- Long RNA-seq metadata: 4,738 rows; 4,632 valid unique participant-visit pairs; 106 rows lacked a parseable visit key.
- Archived whole-blood collection export: 13,793 rows and 13,793 valid unique participant-visit pairs, with no duplicate participant-visit pairs. It overlapped 4,521 RNA-seq visits exactly.
- The collection export's `BLDCOLL` field contained codes 0, 1, and 2. Among the 4,521 exact overlaps, the counts were 206, 4,313, and 2, respectively. These codes are reported without semantic labels because the applicable data-dictionary code mapping has not yet been verified.
- `BLDCLLDT` had month/year precision for 13,235 rows and was blank for 558; 4,316 of the exact RNA-seq overlaps had a parseable month/year value.
- Both CBC exports' `LCOLLDT` values also had month/year precision. Among exact RNA-seq/collection overlaps, the current table had numeric core-marker records for 1,923 visits, of which 1,824 had the same reported month as `BLDCLLDT`. The archived table had records for 1,999 visits, of which 1,899 had the same reported month.

## Interpretation and limits

The new archive provides a substantial exact visit-level crosswalk and month-level temporal information. If `BLDCLLDT` and CBC `LCOLLDT` are both collection dates as their names imply, same-month pairs are within 30 calendar days. Month/year precision cannot establish exact day differences or choose the nearest record across adjacent months, so the full registered same-day/nearest-within-30-days rule has **not** been verified. No temporal eligibility rule is declared passed. The code meanings for `BLDCOLL` also remain unverified and are not interpreted here.

This is a data-linkage and availability audit only. It is not the positive-control CBC score validation, does not establish cell-composition transfer, and does not authorize disease scoring on PPMI. The second public lock remains required before any PPMI disease score.

## Reproduction

Run `pipeline/audit_ppmi_collection_linkage.py` with local paths for the Long RNA-seq metadata, archived whole-blood collection CSV, current CBC CSV, and archived CBC CSV. The aggregate JSON output should be written under ignored `private_data/`; do not add restricted inputs or row-level outputs to Git.

The local aggregate output is `private_data/ppmi_collection_crosswalk_20261018.json` and is intentionally excluded from version control.
