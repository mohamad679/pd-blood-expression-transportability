# Phase 2 progress — draft only, not a registration

Date: 22 October 2026. The author approved Phase 1. Phase 2 is **in progress and not approved for real-data execution or Zenodo publication**. This report contains no PPMI outcome result.

## P2-01 access and publication gate

The approved PPMI application text and written DPC advice about a methods-only Zenodo record are unavailable. The public DUA and publication policy do not settle those project-specific questions; see `PHASE_2_ACCESS_GATE.md`. The author response permitted local drafting and synthetic-only checks while treating the amendment as unpublished. No E1–E10 analysis was run on real GEO or PPMI data. E3 PPMI model training in particular remains gated on approved project scope. The clarification message has been drafted, not sent.

## Drafted methods

`revision/PHASE_2_AMENDMENT_DRAFT.md` specifies the rationale, E1–E10, seed streams, inputs, outputs, interpretation rules, excluded variants, ordering, and PPMI embargo. Its status is explicitly post-primary and unregistered. The original confirmatory result, models, and input matrix are unchanged. The draft maps I-01 to E1–E4; I-03 to E6; I-05, I-06 and I-22 to E7; I-10 to E8/E10; I-13 and I-14 to E5; I-15 to the excluded-variant table; I-02 to E9; and I-33 to the access gate.

The draft operationalizes three points beyond the roadmap's exact wording and will flag them in any registered version: E3 adds an explicit adverse-increment branch for a mean paired AUROC difference below −0.02; E9 requires a reproducible genetic-status field and a defined sporadic de novo subgroup; E6 applies its six-hour runtime rule to the first 20 *attempts*, including failed fits, and limits the total attempted draws. None of these was chosen using E1–E10 results.

## Independent fresh-context review

A fresh-context second reader examined the draft without editing files or accessing PPMI data. The review identified a mixed-model comparison in E1, an undefined matched bootstrap in E7, and ambiguity about which E8 fits score PPMI. The draft was corrected: E1 cross-platform SMDs use the same scope fit on both populations; E7 regenerates the primary full-sample draws and separately reconstructs the prespecified pair draws; E8 scores PPMI only with the pooled GEO fits. The reviewer also confirmed that the specified embargoed primary-result number strings do not appear in the draft. Access/DPC gates and the software commit placeholder remain open.

## Code status and tests

Draft scripts now exist one-to-one for E1–E10, with shared private-output and execution-gate helpers. **All ten synthetic-only self-tests passed** under the locked Python 3.11 environment, and all eleven Python modules compiled. E3's test runs the actual baseline/gene fitting path with reduced artificial dimensions and donor-separated nested folds. E4, E5 and E8 test the separated GEO-training/PPMI-scoring path; E7 checks paired-bootstrap parity on synthetic scores; E10 checks same-month CBC matching. Tests do not read actual PPMI or GEO analysis data. The E1 scope manifest is supplied explicitly at execution, and each scope's sample list/order is checked against the frozen split manifest; each model's fitted-sample hash is also checked. A standalone scope-manifest *builder* was not needed for the verification contract.

Real-data execution in every script requires a separate Phase-3 authorization record containing an approved-application evidence hash, approved analysis name, and registered amendment DOI. No such record exists. The draft code baseline is commit `9d70630845bcd6de08dfcdd59b8d5c9006256c13`. Synthetic pass does not certify every real-data input schema or its runtime; those will require a gated dry run after registration. E9 remains conditional on a reproducible genetic-status field. E3 remains conditional on project-scope confirmation. The scripts are **not** a license to execute either one on PPMI.

## Phase gate

Phase 2 is **TECHNICALLY PREPARED BUT PARTIAL**. The application/DPC questions remain unresolved; the amendment is an unpublished local draft. It cannot be submitted to Zenodo until that guidance is checked and the author separately approves publication. Phase 3 remains closed, and no real E1–E10 analysis has run. The original locked result and model artifacts were not modified.
