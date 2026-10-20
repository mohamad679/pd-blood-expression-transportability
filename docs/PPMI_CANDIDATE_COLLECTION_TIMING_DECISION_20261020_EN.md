# PPMI candidate blood-collection timing metadata audit

**Status:** decision fixed before running this candidate-level metadata audit. This is descriptive, pre-second-lock, and score-free.

Read only `candidate_metadata.jsonl` from the frozen 472-candidate private bundle (SHA-256 `1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c`) and the archived whole-blood collection CSV on the researcher's Mac. The 472 candidates retain their existing historical PD/control labels. Join by exact normalized PATNO and `EVENT_ID=BL`; normalization removes only a leading `PP-` prefix. No fuzzy matching, date-based relabeling, or change in the candidate population is allowed.

For each group, count candidates with no collection row, exactly one matching row with a parseable `MM/YYYY` `BLDCLLDT`, exactly one row with a missing/unparseable date, and multiple matching rows. Among multiple rows, distinguish all parseable months equal from conflicting or missing months. Do not choose a row from a duplicate set or interpret undocumented `BLDCOLL` status codes. Report only aggregate counts, file hashes and total row counts, never identifiers, dates, individual rows or expression. Month/year precision cannot establish exact draw day or medication-at-draw exposure.

Medication exposure is a separate audit. No medication file is currently present among the known Mac inputs, and no diagnosis-based assumption of untreated status is allowed. This audit may document collection-date coverage now; it cannot declare the medication-at-sampling prerequisite complete. Commit and push this decision and script before running the real candidate audit.
