# CBC metadata availability across visits in baseline controls

**Date:** 2026-10-18. This corrects the interpretation of the earlier same-visit audit. No CBC measurement values, expression, cell scores, correlations or disease predictions were used in this audit. The existing metadata eligibility rule was reused to identify the same 161 baseline healthy controls.

## Protocol check and correction

The frozen protocol, `PROTOCOL_v3_1_EN.md`, CBC section, prefers the blood-draw day and otherwise the nearest measurement within ±30 days. It does **not** require the CBC and RNA sample to have the same visit code. Therefore the earlier exact-visit count of zero was correctly calculated, but its use to declare CBC unavailable was too restrictive. The earlier suggestion that any cross-visit join would be an improper post-hoc rescue is withdrawn. This decision is based on metadata availability before any CBC association or PPMI score, and the combined amendment remains unpublished.

## Metadata-only results

Only `PATNO`, `EVENT_ID` and `LCOLLDT` were retained from each CBC export. CBC table rows represent individual laboratory entries, not independent people or separate blood draws. No marker-specific completeness is inferred from these row counts.

| Availability | Current export | Archived export |
|---|---:|---:|
| Frozen baseline healthy controls | 161 | 161 |
| Controls with a CBC table row at any visit | 161 | 161 |
| Controls with a CBC table row at `SC` | 161 | 161 |
| CBC rows at `SC` | 5,406 | 5,406 |
| Controls with a reported BL blood-collection month | 122 | 122 |
| Controls with any CBC row in that BL collection month | **31** | **31** |
| Same-month matches from `SC` | 31 | 31 |
| Same-month matches from `BL` or other visits | 0 | 0 |
| Same-month CBC table rows | 980 | 980 |

The full aggregate `EVENT_ID` breakdown is in `results/ppmi_cbc_metadata_20261018/any_visit_coverage.json`. CBC rows at other visits are present, but none supplies additional same-month matches in this control set. The all-visit availability count reaches 31, whereas the earlier same-visit count was zero. This supports retaining the planned CBC validation subject to the unchanged per-marker minimum sample size and ambiguity checks. **31 is an availability upper bound, not a confirmed numeric sample size for every marker.**

## Proposed fixed operational rule

1. Keep the same baseline healthy-control population; require a usable collection month from its exact baseline blood-collection record. No `BLDCOLL` code filter is applied. Its values may be counted without interpretation, and obtaining their codebook is not an execution gate.
2. Search CBC rows for the same participant at **any visit**, retaining only the same reported calendar month as the baseline collection. Exclude adjacent months. Same-month matching guarantees at most 30 calendar days if the reported fields are collection months; it does not identify same-day or closest-day records and omits some genuinely eligible cross-month pairs.
3. For each participant and prespecified laboratory marker, prioritize the RNA visit (`BL`), then `SC`, then other visits. This priority is fixed before measurement values or associations are examined.
4. Month/year precision cannot identify the earliest day. Do not use CSV order or measurement magnitude to resolve ties. Collapse exact duplicate representations of the same identified measurement only when their values and units agree; exclude unresolved multiple collection occasions or conflicting repeats at the highest available visit-priority tier for that participant-marker, and report these exclusions. No fallback to a lower tier is made after a conflict at the preferred tier.
5. Use the current export as the primary source; the archived export is an availability/concordance check and is not concatenated as additional observations. Determine valid numeric sample size separately for every fixed marker after this rule is publicly versioned. Apply the original `n >= 30`, Spearman and confidence-interval criteria without changing panels or selecting markers by results. If a marker has fewer than 30 usable pairs, report it as non-evaluable.

The rule clarifies handling of month-level dates and cross-visit records. It will be published with the permutation amendment before any CBC association. It does not change M4/M5 or the main PPMI population, and CBC is not a stopping gate for PD scoring. CBC supports validation of specific cell-composition proxies, not certification of general autosomal or PD-specific transfer.

**Fixed interpretation before marker completeness:** the 31 candidates leave little margin above the per-marker minimum of 30. Missing values, unit problems or conflicting records may legitimately make a marker non-evaluable; the matching rule will not be reopened to increase n. Any passing association at approximately 30 participants will be described as limited supportive evidence. The same-month rule can omit genuine ±30-day pairs in adjacent months and therefore may retain fewer pairs than a day-resolved rule. We have not measured how many of the other 91 controls with known BL collection month fall within 30 days, or established why each did not match.

## Verification and next stage

`pipeline/audit_ppmi_cbc_any_visit.py` reproduces this audit. Synthetic tests verify a same-month screening match, adjacent-month exclusion, unique-person counting despite repeated assay rows, missing baseline-date exclusion, and operation without `BLDCOLL` or CBC result columns. Restricted exports remain outside Git. Next: test the permutation executor's membership, seed, resume and numerical behavior, benchmark worker/memory settings, then publish the combined amendment before running any real permutation or CBC association.
