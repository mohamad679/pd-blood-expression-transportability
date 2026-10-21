# Phase 2 access and publication gate — P2-01

Date: 22 October 2026. Status: **UNRESOLVED; stop before E1–E10 scripts, synthetic executions, amendment registration, or PPMI exploratory analysis.** This memo contains no PPMI outcome result and makes no assertion that the investigator's approved application covers the proposed post-lock analyses.

## Documents checked

- [PPMI Data Use Agreement, version 5.0, April 2026](https://www.ppmi-info.org/sites/default/files/docs/ppmi-data-use-agreement.pdf), especially clauses 2–5 and 9–13 and Appendix A.
- [PPMI Publication Policy, September 2024](https://www.ppmi-info.org/sites/default/files/docs/ppmi-publication-policy.pdf), especially external publications, DPC review, and non-manuscript outputs.
- Local `PPMI_Data_Access_Guidelines_202604.pdf`, version 8.0, 1 April 2026: RNA-sequencing large datasets are listed as Tier 2; the guidelines concern data access tiers and do not approve a particular analysis plan.
- Local `PPMI_Tier_2_Request_Form_Rev_20260609.docx`: the available copy is a blank form and does not show the researcher's submitted project description, approval, or signatures.
- Existing private DPC manuscript request draft: a proposal for later administrative manuscript review, **not** evidence of approved access scope or prior DPC advice about a methods-only Zenodo amendment.

## What the public documents establish

The DUA permits scientific investigation and, in Appendix A, explicitly contemplates model training using PPMI data when the environment complies with the DUA. It also requires investigators to accurately describe the planned analyses in their data-use application and prohibits disclosure of participant-specific data to unauthorized parties or AI tools that cannot guarantee containment. A local, user-operated Mac and aggregate-only reporting reduce the disclosure path, but do not prove that a new within-PPMI training design is covered by the investigator's application. No participant-level PPMI data should be supplied to this assistant, Git, Zenodo, or the clarification message.

The publication policy requires DPC administrative review of an external **manuscript** before journal submission. The DUA asks investigators to endeavor to submit non-manuscript scholarly outputs to the DPC and a public platform. Neither public document unambiguously states whether a *methods-only amendment registered in Zenodo after the primary result, with no PPMI result numbers*, requires DPC review **before** the Zenodo post. The public texts therefore do not justify treating prior publication as automatically approved. Regardless, the roadmap independently requires explicit investigator approval before any Zenodo publication.

## Unresolved facts and decisions

1. The approved data-use application or project description is not among the locally located files. Confirm whether it covers **E3: new nested cross-validation and model training within PPMI** (and, if relevant, E4/E5/E6/E8 refits whose scoring uses PPMI). A general DUA permission to train models is not proof that this proposed analysis is within the approved application. If the description is narrow, ask PPMI whether an application amendment is needed before E3.
2. Obtain written PPMI/DPC guidance on whether a public, result-free methodological amendment may be posted to Zenodo before DPC review, and whether the corresponding code-only record needs prior review. Until clarified, prepare files locally only; do not publish them.
3. If the application scope is broader and the investigator can provide its non-sensitive analysis wording, record the evidence and resolve this gate before advancing P2-02 through P2-06. A verbal assumption without the approved wording is insufficient for a compliance claim.

## Draft clarification text for investigator review — not sent

> Our approved project evaluates transportability of blood-expression classifiers using PPMI Project 133 RNA-sequencing data on an investigator-controlled Mac. We propose post-lock, explicitly exploratory analyses including a nested cross-validation comparison of a cell-marker baseline with a gene-augmented model within PPMI. Only aggregate outputs would be reported; participant-level data, predictions, and fitted PPMI model weights would not be placed in public repositories or external AI services. Is this analysis within the scope of the existing approved data-use application, or should we update it before proceeding? Separately, may we post a methods-only amendment and code to Zenodo **without PPMI results or participant-level material** before these new analyses, or should the Data and Publications Committee review that protocol-only record first?

This is a draft question, not a representation that DPC has answered it. No communication was sent.

## Gate status

P2-01: **PARTIAL / requires investigator application wording and PPMI/DPC clarification on Zenodo timing.** P2-02–P2-06: **NOT STARTED**. E1–E10: **NOT RUN**. Second-lock files and stored predictions: **unchanged**. No amendment or code was published.
