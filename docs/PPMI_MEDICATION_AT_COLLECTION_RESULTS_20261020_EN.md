# PPMI recorded dopaminergic medication at collection: descriptive audit

The [rule and source-field mapping](PPMI_MEDICATION_AT_COLLECTION_DECISION_20261020_EN.md) were committed and pushed as `d60995c` before the candidate-level medication count. The Mac-local audit used the frozen 472-candidate bundle, the archived whole-blood collection export, and the `LEDD_Concomitant_Medication_Log_20Oct2026.csv` export (10,294 rows; SHA-256 `cdab712e45d2a57d4e81b27fa414b40f228a9d7216edf31d55891beb69b4e10b`). It emitted aggregate counts only. The private aggregate report is `private_data/ppmi_ledd_at_collection_20261020.json`.

| Historical group | Candidates | Collection month unknown | Known month, no LEDD log rows | Known month, LEDD rows but no recorded ongoing use | Known month, indeterminate LEDD timing | Known month, definitely ongoing LEDD use |
|---|---:|---:|---:|---:|---:|---:|
| PD | 311 | 54 | 16 | 236 | 5 | 0 |
| Control | 161 | 39 | 121 | 1 | 0 | 0 |
| Total | 472 | 93 | 137 | 237 | 5 | 0 |

The categories are mutually exclusive and reconcile to the frozen 472 candidates. Across LEDD records linked to candidates with a known collection month, 2,449 records were not ongoing by their dates, four had a start or stop in the collection month, and one had an internally inconsistent date interval. No record met the prespecified definitely ongoing rule. These results describe **recorded dopaminergic medication in the LEDD log only**, not all medication exposure or proof that participants were untreated. The 93 with no parseable collection month have unknown timing even if they have medication records. No candidate was excluded, relabeled, or scored; the confirmatory analysis and second-lock gates are unchanged.
