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

Draft code exists for E1, E2, E3 and E9, with shared private-output and execution-gate helpers. Each has a synthetic-only self-test that passed under the locked Python 3.11 environment. E3's test runs the actual baseline/gene fitting path with reduced artificial dimensions, donor-separated nested folds, and aggregate-only progress; warnings are captured in diagnostics rather than flooding the terminal. Tests do not read PPMI or GEO data. Real-data execution in all four scripts requires a separate Phase-3 authorization record and a registered amendment DOI. No such record exists.

The E1 scope-manifest builder/hash verification and E4–E8/E10 scripts are **not yet complete**. The E1, E2, E3 and E9 scripts remain drafts until their real-data input contracts and all numerical choices receive one final code review. The one-to-one E1–E10 script check therefore **fails at present**. No code commit hash may replace the amendment placeholder yet.

## Phase gate

Phase 2 is **PARTIAL**. It cannot be marked complete or submitted to Zenodo until the application/DPC questions are resolved, all ten scripts and synthetic checks are complete, the software commit hash is inserted, disclosure scanning is rerun, and the author separately approves publication. Phase 3 remains closed.
