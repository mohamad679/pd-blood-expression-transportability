# PHASE 2 REPORT — local draft, registration pending

Date: 22 October 2026. This report and the amendment are private working documents. No E1–E10 analysis has run on actual GEO or PPMI data, and no PPMI outcome result is stated here.

## 1. Outcome

The methods, interpretation rules, execution order, exclusions, and synthetic-tested code for all ten post-lock diagnostics are prepared. Phase 2 has **not reached its registration objective**: the approved PPMI project description and written guidance on public methods-only registration are unavailable. The amendment is therefore unpublished, and Phase 3 is closed.

## 2. Work-item status

| ID | Status | Evidence |
|---|---|---|
| P2-01, I-33 | PARTIAL | Public policy was reviewed, but it does not establish whether the approved application covers E3 training or whether DPC review is required before a methods-only Zenodo amendment. See `PHASE_2_ACCESS_GATE.md`. Clarification text was drafted, not sent. |
| P2-02, I-01/I-03/I-05/I-06/I-10/I-13/I-14 | DONE as draft | E1–E10 specify inputs, methods, seed, output, interpretation, and status as post-primary diagnostics. |
| P2-03, I-15 | DONE as draft | Full-genome elastic net, supervised filtering, and 438-record GENEPARK M2 are explicitly excluded with reasons. |
| P2-04 | DONE as draft | E1/E3 interpretation thresholds, E4/E5 exploratory status, staged GEO-before-PPMI ordering, and reporting of all outcomes are fixed. |
| P2-05 | DONE for synthetic preparation | One script per E1–E10, shared gate helpers, locked environment, all ten self-tests passed. Software baseline commit below. |
| P2-06 | PARTIAL | English amendment drafted and independently reviewed; no Zenodo registration or DOI. |

## 3. Incomplete work and exact gate

The approved PPMI application wording is needed to determine whether E3's nested within-PPMI model training is in scope. Written PPMI/DPC guidance is needed before publishing a methods-only Zenodo amendment, unless the approved agreement itself clearly resolves that question. The author reported that neither text is currently available. Until this is resolved, the draft remains local, no registration is claimed, and the execution helper rejects real E1–E10 runs unless a separate Phase-3 authorization record supplies the relevant evidence hash, approved analysis, and registered DOI. E9 additionally depends on a reproducible genetic-status field; absent that field it is reported as not evaluable.

## 4. Deliberate operationalizations and deviations

The amendment is explicitly **post-primary**, not a prospective registration for the already observed primary result. E3 has an explicit adverse-increment branch below −0.02. E6's six-hour rule uses the first 20 attempted replicates, including failed fits, and caps attempted draws. E7 deterministically reconstructs the original full-population bootstrap stream, checks parity with the locked quantiles, and constructs a separate matched-pair stream; it does not claim that bootstrap index arrays were saved. E9 requires a defined sporadic de novo subgroup and at least ten donors per class. E10 reports nine predefined comparisons across six panels and treats lymphocyte percentage and RBC/hemoglobin as imperfect proxies. These choices were recorded before E1–E10 real-data execution. The review corrected E1's cross-population comparison to use the same frozen scope fit and limited E8 PPMI scoring to pooled GEO fits.

## 5. Files, SHA-256, and code baseline

The E1–E10 scripts and shared helper were committed by `9d70630845bcd6de08dfcdd59b8d5c9006256c13` (E1/E2/E9/helper originated in the earlier private commit). The report and amendment are subsequently committed as documentation. Hashes below are SHA-256 of the final working files at report creation:

| File | SHA-256 |
|---|---|
| `revision/PHASE_2_AMENDMENT_DRAFT.md` | `a245ce210b26cbc6d2fbd40693f729e6190f84a7ed4a6e752263c11ae3a0f02f` |
| `revision/postlock_common.py` | `b7509724a8173908d85a8795b931443a45b232ef7cb780741e696b8b61b592c1` |
| `revision/e01_logit_decomposition.py` | `a43250596f70bf150db22d948137498a2d6c23ab7452a518a93140a3d2cd9c60` |
| `revision/e02_rank_shift.py` | `493eedd524083b59179cdf017ebe41984ed4d81bb9f10cc32036088b703c55b8` |
| `revision/e03_ppmi_nested_cv.py` | `8f4abaa86a7e6badfe1c4a984e8657570d7a44b7d5a51940327406436a14bbc6` |
| `revision/e04_transductive_universe.py` | `1f40ac5f0288f50e2fbc5e3bea85765f1c6342c5fc8060c5763b9c479e30d538` |
| `revision/e05_weighted_deviance_panel_exclusion.py` | `d3e2c0ad18eaa442b31ccc454f85329621c517ea5adf2b5b53630eee4ce81765` |
| `revision/e06_workflow_bootstrap.py` | `19327514321c385be731de214366aab78f9ab6cce86c1e015e4e76aa037e6305` |
| `revision/e07_locked_metrics.py` | `f58eade35282ed9115702688cb85317429987f5ea2251562f3b0dcde5b0c816f` |
| `revision/e08_panel_ablation.py` | `57afd6228e41541aedb7584fd55fb38993f9e4ad8855f9a07313e43da399d19d` |
| `revision/e09_sporadic_subgroup.py` | `9edd48d587853c1364fb4b431e849b76d0d6509580e32bbf6ca5d0670bab3b14` |
| `revision/e10_cbc_all_groups.py` | `e31c8edb99c418fa9328ae9fb9250115cc240d1dcd2d190577c1bfac2423d5c8` |
| `revision/reports/PHASE_2_ACCESS_GATE.md` | `80ff7d5a21052452d2780476b1a6e80fb7ab00b4cc2ac87166ef4fe9169202bc` |

## 6. Double-check

| Check | Result |
|---|---|
| Phase-2 issue IDs map to an E analysis, exclusion, or access gate | PASS |
| E1–E10 script list is one-to-one with amendment | PASS |
| Specified primary-result number strings absent from amendment | PASS: zero matches in the disclosure scan |
| Ten scripts pass synthetic-only self-tests | PASS, locked Python 3.11 environment; `py_compile` also passed |
| No actual PPMI or GEO analysis data accessed by these self-tests | PASS by the synthetic entry paths; no real-data command was run |
| Interpretation rules are numeric and fixed where the roadmap requires | PASS; limitations are explicit |
| Independent fresh-context review | PASS; ambiguities and their corrections recorded in `PHASE_2_PROGRESS.md` |
| PPMI application and DPC publication permission confirmed | OPEN |
| Zenodo amendment registration | NOT DONE; separate author approval and policy clarification pending |

## 7. Unexpected findings

The prior primary bootstrap report stored quantiles, not draw indices, so E7 must regenerate its deterministic draw stream and verify parity before descriptive extensions. The E1 scope mapping needed direct checks against frozen split membership and fitted-sample hashes; those checks are now in code. No diagnostic outcome was inspected.

## 8. Decision needed from the author

When available, provide the approved project-description wording or written PPMI/DPC reply addressing whether E3 model training is in scope and whether a result-free amendment may be registered publicly before DPC manuscript review. A separate explicit instruction is required to publish in Zenodo. No decision about the scientific results is requested at this gate.

## 9. Lock status

PASS for changes made in this phase: no tracked frozen model, prediction, or matrix artifact was edited by the Phase-2 work; the Phase-1 lock comparison remains the last full reference verification. A fresh full private-artifact hash comparison is deferred until the authorized pre-execution gate and is not represented as having been repeated here.

## 10. Stop point

Phase 2 remains partial and Phase 3 has not started. This report and draft await review; after the access/publication questions are resolved, the phase can be finalized. Do not interpret this report as consent to register or to execute real E1–E10 analyses.
