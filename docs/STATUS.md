# Project status

**Snapshot date:** 2026-10-19 (Asia/Tehran)

**Cloud permutation result:** The registered exploratory B=200 M4/M5 GEO
label-permutation run completed with eight workers. The returned summary has
`status: complete`, `B: 200`, and `PPMI_accessed: false`. All 200 fixed
replicate completion markers are in the returned private archive. The server
archive transferred to the Mac, passed tar/gzip inspection, and matched SHA-256
`50d9b536f5f0f46a176605859134e23f143bd3994840d885074dc79ca5a7fead`.
A verified ignored local copy is in `private_data/server_return/`. Hetzner
Console subsequently showed no servers and no Primary IPs. Read the
[aggregate result report](GEO_PERMUTATION_RESULTS_20261019_EN.md) and the
[updated new-chat handoff](NEW_CHAT_HANDOFF_20261018_EN.md). The second public
pre-PPMI lock and external scoring remain outstanding.

The complete retrospective account is [COMPLETE_PROJECT_REPORT_20261018_EN.md](COMPLETE_PROJECT_REPORT_20261018_EN.md); subsequent work is summarized in the [continuation report](PROJECT_CONTINUATION_REPORT_20261018_EN.md). Its [reporting correction](REPORTING_CORRECTION_20261018_EN.md) records a terminology error in the released short report: final M5 (`rho=1`) has an L1/lasso gene penalty, not ridge. Code, hyperparameters, fitted models and metrics are unchanged. The full report also records two persistent lambda boundaries from the completed runs; no extra tuning has been performed.

| Workstream | Status | Evidence / limitation |
|---|---|---|
| Feasibility and cohort choice | Completed for the registered scope | Two GEO development cohorts and an approved PPMI access path were identified. Feasibility did not itself constitute external validation. |
| GEO input and mapping audit | Completed | 381 GEO records, 16,364 shared genes; all four prepared GEO input hashes match the frozen audit. |
| Protocol and implementation lock | Publicly versioned | Pipeline and split hashes matched Zenodo amendment 3.1.1 before GEO model execution. |
| GEO nested CV and directional transfer | Completed for all 13 registered scopes | Five outer folds per GEO source, both directional transfers, and pooled final tuning/fitting. Aggregate results and numerical audit: [GEO_CV_RESULTS_20261018_EN.md](GEO_CV_RESULTS_20261018_EN.md). |
| Numerical fit and sigmoid audit | Completed for the six reported warning contexts | All 3,146 original post-pilot fit records passed KKT ≤ 1e-6. A separate 490-fit replay reproduced all six warnings, all original candidate means exactly, and all selections; clipping changed no AUROC or selection. The locked engine and original artifacts are unchanged. [Audit](SIGMOID_AUDIT_20261018_EN.md). |
| GEO sex/age transfer controls | Completed | Sex AUROC 1.000 and 0.993 passed the two development directions; age Spearman -0.058 and 0.259 did not establish development capacity. This limits age-control claims without changing the pipeline. [Report](GEO_TRANSFER_CONTROLS_20261018_EN.md). |
| Saved GEO fit membership and preprocessing | Passed for all 13 original scopes | All 3,415 fit records' training hashes, regenerated train-only preprocessing fingerprints and Kish weights matched the frozen splits; saved final preprocessing and held-out membership/outcomes matched. No new model fits. Operational donor checks do not prove biological cross-source independence. [Evidence](../results/permutation_preflight_20261018/geo_membership_audit.json). |
| Remaining registered diagnostics | Incomplete | CBC marker validation, bounded sensitivities, remaining leakage checks and bootstrap coverage remain; B=200 permutations have completed and are archived. A metadata-only cross-visit audit found 31 baseline controls with same-month CBC rows at SC; numeric n for each marker is pending. See [follow-up audit](PPMI_CBC_ANY_VISIT_AUDIT_20261018_EN.md). |
| PPMI linkage and final population | Incomplete | The provisional 311 PD / 161 control count is not final. The identity/crosswalk and eligibility set must be reconciled and locked. |
| PPMI external scoring | Not started | No PPMI disease scores have been calculated. A second public lock of code/environment, final models/preprocessing, and population reference must precede any PPMI score. |
| Manuscript | Not started | GEO results are now available; external validation and remaining diagnostics are not complete. Causal shortcut fractions and clinical utility are not inferred. |

## GEO milestone

The core GEO milestone was published on 18 October 2026 as version **3.1.3-geo-cv**, DOI [10.5281/zenodo.23445851](https://doi.org/10.5281/zenodo.23445851), in the existing concept-DOI lineage. The four uploaded file checksums matched the local files. The publication receipt and release manifest are in `results/geo_development_20261018/`. This is a post-results snapshot and does not satisfy the second pre-PPMI lock.

All planned nested outer folds, directional transfers, and pooled final fitting completed on 17–18 October 2026. Across the five outer folds, mean AUROC was 0.591 for M5 versus 0.562 for M4 in GSE57475 and 0.686 versus 0.663 in GSE99039. In directional transfer, M5 was below M4 by 0.178 AUROC from GSE57475 to GSE99039 and by 0.009 in the reverse direction. These are development results for the registered procedure and two GEO sources; they are not PPMI performance or causal estimates.

The four input hashes match the frozen audit, and the split-manifest SHA-256 is `2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`. The public implementation manifest covers code SHA-256 `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`. No individual-level predictions or PPMI data are in Git.

A byte-for-byte verified backup of the 167 post-pilot run files is retained locally in the ignored `private_data/geo_runs_20261018/` directory, in addition to the temporary execution outputs. These files were not included in the public upload.

## Next gates

### Cloud preflight and amendment publication completed — 18 October 2026

The registered pilot completed (23.77 minutes, peak RSS 0.52 GiB), and its
replicate 0 was reused in the eight-worker fixed B=200 run. The full run has
now completed. The returned archive has all replicate completion markers and
a complete summary; its source, protocol and public binding were frozen before
execution. The read-only monitor used during execution changed no scientific
file. Six aggregate statistics and interpretation are in the
[completed result report](GEO_PERMUTATION_RESULTS_20261019_EN.md).

The returned Linux certificate passed all 33 synthetic checks (650 fit attempts,
maximum KKT 2.6839e-8) and all 20 GEO model/fold identity comparisons (zero AUROC
error; maximum probability error 7.6494e-14). All eight tested worker processes
were observed. The conservative memory budget was 5,996,040,192 bytes against
32,848,162,816 bytes RAM, with zero pswpout growth. An overlapping controller
was rejected; the original controller subsequently completed without deleting
its lock or changing code. These are preflight results, not permutations.

The combined amendment was published as version **3.1.4-permutation-cbc**, DOI
[10.5281/zenodo.23449296](https://doi.org/10.5281/zenodo.23449296), at
2026-10-18T13:46:35.803828+00:00. All four public file contents were downloaded
without authentication and verified byte-for-byte against the frozen release.
The exact code/config/input/machine binding and public-verification receipt are
recorded in `results/zenodo_permutation_cbc_v3_1_4/` and
`results/cloud_preflight/ZENODO_PUBLICATION_RECEIPT.json`.

At publication verification, no randomized GEO replicate, CBC association or
PPMI disease score had been executed. That historical statement is preserved
in the release evidence; pilot and full-run execution followed publication.
The server setup, encrypted GEO-only transfer, exact environment, preflight,
publication and pilot are now complete. Full execution has since completed and been returned to the Mac. The server
and Primary IP were deleted. No PPMI files were included in the cloud bundle.

1. **Completed:** the six-warning review found no effect of explicit clipping on AUROC or parameter selection; no scoring change or core scientific rerun is warranted by this diagnostic.
2. Sex/age controls are complete; complete the bounded sensitivities, leakage checks, bootstrap-coverage validation, and permutation plan. Age failure does not authorize a rescue control or holdout retuning.
3. Finalize PPMI file-to-release linkage and eligibility using the guide and complete sample-name segments; document publication conditions.
4. Before PPMI scoring, publish the second lock of code/environment, splits, fitted M1–M5 models and preprocessing, and a non-identifying reference/hash for the private population/crosswalk.

## Historical planning notes: permutation and reporting

**Historical record, superseded for execution/publication status by the entries
above and the new-chat handoff.** References below to an unpublished draft, a
failed Mac resource gate, or no real permutation describe earlier preparation.
They do not authorize repeating publication, preflight or the completed cloud run.
CBC selection/interpretation is now governed by the published joint amendment.

- Proposed scope: B=200, M4/M5, ten within-GEO outer folds. The [amendment draft](../protocol/PERMUTATION_AMENDMENT_DRAFT_20261018_EN.md) is **not published or effective**; no registered real-GEO permutation has run. Executor, membership/seed/resume tests and worker/memory benchmarks precede publication and execution. Measured M4/M5 fit time is 18.38 minutes per original pass, implying a 61.27-hour serial baseline for B=200 before overhead/null-runtime changes.
- **Executor preflight — current resource gate not passed:** a pre-hardening run passed 32 checks and the observed-label ten-fold replay then matched all 20 M4/M5 comparisons (maximum probability difference 1.1102e-16). Review subsequently added a guard against treating terminal numerical failure as an interruption, strengthened receipt binding, and added the corresponding test. The final candidate passes 32 of 33 checks and all 650 fit attempts (maximum KKT 2.6839e-8), including worker/order invariance, seed derivation, inner re-stratification and checkpoint integrity. Its last resource check fails: after the competing observed replay finished, the two-worker/two-thread benchmark still recorded 23,780 additional system compressor-swapout pages (about 92.9 MiB), despite a conservative memory budget of 2.124 GiB on the 8 GiB Mac. No process-level cause is asserted and no threshold is relaxed. Final-binding observed replay/certification is pending; no amendment publication or real permutation has occurred. The one-command Bash preflight is prepared for a host with less memory pressure. All prior attempts are retained. See [report](PERMUTATION_EXECUTOR_PREFLIGHT_20261018_EN.md).
- CBC is a priority for the already planned marker-score validation. After user-completed IDA login, the current `Blood Chemistry & Hematology` and archived table were located. The browser download did not deliver them, but the user subsequently supplied both local exports. CBC supports specific marker/measurement associations, not certification of all autosomal or PD transfer. The final GEO-trained sex control remains to be fitted and included in the second lock.
- **Update:** the user supplied both CBC exports. A local aggregate-only exact `PATNO + visit` coverage audit against the 2021 Long RNA-seq metadata found 4,738 Long RNA-seq rows; 106 lacked a parseable visit key. Of 4,632 matched sample/visit rows, 1,984 had a current-table CBC record and 2,060 an archived-table record (76 archived-only; none current-only). Core CBC measurements were available for 2,001 archived sample/visit pairs. Matching is visit-level; same-day alignment cannot be checked from the available RNA metadata. Raw CBC stays outside Git; no outcomes, scores, or correlations were computed. Aggregate output is local under ignored `private_data/`.
- **CBC version concordance and temporal gate:** among the nine prespecified numeric CBC markers, 56,104 shared exact measurement keys (PATNO, EVENT_ID, test code, collection date, and collection time) had identical values and units in current and archived exports; zero differed and zero were duplicate-ambiguous. This is marker-row concordance only, not a claim that the full tables are identical. The 2021 RNA metadata itself has no collection-date field. The user subsequently supplied the archived `Whole_Blood_Sample_Collection` export. Exact visit-key matching linked 4,521 of 4,632 valid unique RNA visits; 4,316 linked records had month/year collection dates. Both that date and CBC `LCOLLDT` are month/year precision. For core numeric CBC markers, 1,824 current-table and 1,899 archived-table visit pairs share the collection month. This can support a candidate within-30-day subset if both fields denote collection dates, but does not establish exact day differences or nearest-record selection across adjacent months; the full temporal rule remains unverified. `BLDCOLL` codes are not semantically interpreted pending verification against the applicable codebook. See [PPMI collection linkage audit](PPMI_COLLECTION_LINKAGE_AUDIT_20261018_EN.md). No CBC association, marker score, or PPMI disease score has been run; raw CBC and collection files remain outside Git.
- **Baseline-control CBC correction:** the exact-visit audit yielded zero CBC rows, but same visit was an additional restriction absent from the registered same-day/nearest ±30-day rule. The same 161 controls were therefore audited across all visits using CBC participant/visit/date metadata only. All 161 have SC rows; 122 have a known BL collection month; **31** have CBC rows in that same month, all from SC, in both exports. This retains a candidate CBC validation subset; 31 is an availability upper bound and numeric n must still be established separately for each fixed marker. The earlier broad non-evaluability conclusion is withdrawn. The draft now permits same-participant, any-visit, same-month matching, with fixed BL then SC then other-visit priority and explicit exclusion of unresolved ties; month-level dates cannot identify earliest/nearest days. `BLDCOLL` is neither an eligibility filter nor an execution blocker. No CBC measurement value was used in this follow-up audit, and no correlation, expression score or PPMI disease score was computed. See [audit and rule](PPMI_CBC_ANY_VISIT_AUDIT_20261018_EN.md).
- **CBC interpretation lock:** per-marker n<30 is non-evaluable and does not reopen the matching rule to increase n. Passing results at approximately n=30 are limited supportive evidence. The same-month subset may be smaller than the unknown day-resolved ±30-day subset; the monthly fields cannot quantify that difference. These explicit qualifications are in the joint draft before marker completeness or correlations are inspected.
- Give a short report after each stage: what was done, key numbers, unexpected findings. Keep ordinary progress in Git plus a brief status update; full retrospective reports only when requested. Planned remaining Zenodo milestones are the permutation amendment, the pre-PPMI lock and the final release. An unforeseen material method change would still require a timely amendment rather than an undocumented deviation.
- Order: finish membership/leakage and permutation executor preflight; publish the combined permutation/CBC temporal amendment only after executor and worker/memory tests pass; then launch resumable permutations. CBC marker-specific completeness and association analysis follow the published selection rule. Official `BLDCOLL` codebook verification may proceed independently and is not blocking. Schedule bootstrap coverage and sensitivities within measured resource limits. No new control or enlarged age grid is introduced to rescue the observed age result.
