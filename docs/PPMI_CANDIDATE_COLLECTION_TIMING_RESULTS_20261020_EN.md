# PPMI candidate blood-collection timing audit results

The prespecified [metadata-only rule](PPMI_CANDIDATE_COLLECTION_TIMING_DECISION_20261020_EN.md) was run on the researcher's Mac after commit `9db7a1c` was pushed. Inputs were the frozen 472-candidate bundle (SHA-256 `1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c`) and the archived whole-blood collection export (SHA-256 `e609f9adb835e3654a5de2d5c02169549e095cec1d36a1cc7365b7d1027080bc`; 13,793 rows). Only aggregate results were emitted. The private aggregate report is `private_data/ppmi_candidate_collection_timing_20261020.json`.

| Historical group | Candidates | One baseline collection row with parseable MM/YYYY | One row, missing/unparseable MM/YYYY | No baseline collection row | Multiple rows |
|---|---:|---:|---:|---:|---:|
| PD | 311 | 257 | 26 | 28 | 0 |
| Control | 161 | 122 | 4 | 35 | 0 |
| Total | 472 | 379 | 30 | 63 | 0 |

Parseable month/year coverage is 379/472 (80.3%): 257/311 (82.6%) for PD and 122/161 (75.8%) for controls. The four mutually exclusive coverage categories reconcile exactly to 472. These are date-precision and linkage counts, not medication-at-draw classifications. The source field provides month/year only, so the exact blood-draw day remains unknown. The 63 candidates without a row in this collection export nonetheless have baseline RNA samples; absence from this export is not evidence that blood was not drawn and is not an exclusion reason. All 93 candidates without a parseable collection month remain in the population and have unknown medication timing under the [frozen medication rule](PPMI_MEDICATION_AT_COLLECTION_DECISION_20261020_EN.md). Medication exposure remains pending an appropriate medication export and its codebook; no assumption of untreated status was made. No PPMI expression, model, control or disease score was calculated.
