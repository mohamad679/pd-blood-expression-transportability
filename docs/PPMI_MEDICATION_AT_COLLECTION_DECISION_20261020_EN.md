# PPMI medication at blood collection: descriptive metadata rule

**Status:** fixed before opening a medication export or examining medication-by-diagnosis counts. This is a metadata-only descriptive audit, not an eligibility criterion, model input, control gate, or confirmatory outcome. No candidate is excluded or relabeled because of medication or missing medication information.

Apply the rule to the frozen 472 baseline candidates by exact normalized PATNO. Use the candidate-level baseline collection `BLDCLLDT` month/year already audited in [the collection report](PPMI_CANDIDATE_COLLECTION_TIMING_RESULTS_20261020_EN.md). There are 379 candidates with one parseable collection month, and 93 without one (30 with a collection row lacking a parseable month and 63 without a matching row). The latter 93 have **unknown medication timing** regardless of medication-log contents. A missing collection-export row does not imply that blood was not drawn: these candidates have a baseline RNA sample and remain in the population.

For each medication record with a usable start month and a usable stop month or blank stop date:

- **Recorded use at collection:** start month strictly before the collection month, and stop date blank or stop month strictly after the collection month.
- **Indeterminate timing:** start or stop is in the collection month. The exact day is unavailable, so do not infer use or nonuse. Missing or unparseable boundary dates that prevent classification are also indeterminate.
- **No use indicated by that record:** a known start month strictly after the collection month or a known stop month strictly before it.

For a participant with multiple records, any definitely active record establishes recorded exposure; otherwise any indeterminate record yields indeterminate timing; otherwise report **no recorded ongoing use**. Keep a separate count for participants with no medication-log rows, who must not be called untreated. The name and scope of each medication export and the field mapping from its codebook must be recorded before computing counts. In particular, an LEDD/dopaminergic-medication log only supports conclusions about drugs covered by that log, not all medication classes. Report the groups' counts and coverage, including unclassifiable records, without participant IDs or row-level dates.

If the export or its dates are incomplete, report the limitation and proceed with the second-lock preparation. Medication is descriptive and cannot delay the lock for days or change the confirmatory analysis. No expression values, fitted model, or score is accessed in this audit.

## Source mapping fixed before the medication count

The researcher supplied `LEDD_Concomitant_Medication_Log_20Oct2026.csv` on the Mac (SHA-256 `cdab712e45d2a57d4e81b27fa414b40f228a9d7216edf31d55891beb69b4e10b`). Header-only/schema inspection found `PATNO`, `LEDTRT`, `STARTDT`, and `STOPDT`, with start/stop fields in month/year form or a blank stop. Map these to participant, recorded dopaminergic treatment, start month, and stop month, respectively. The [PPMI Data User Guide](https://www.ppmi-info.org/sites/default/files/docs/PPMI%20Data%20User%20Guide.pdf) identifies this log as a history of **dopaminergic** medication and uses these field names; it is not a complete history of other medications. Do not infer untreated status from absence of a row. A blank treatment name or internally inconsistent start/stop order is unclassifiable, not proof of exposure or nonexposure. Count those records separately. This field mapping was documented before any candidate-level medication classification was computed.
