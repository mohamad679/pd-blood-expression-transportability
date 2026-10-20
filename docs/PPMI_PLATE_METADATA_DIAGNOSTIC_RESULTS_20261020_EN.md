# PPMI plate metadata diagnostic: aggregate results

The method was registered in `PPMI_PLATE_METADATA_DIAGNOSTIC_DECISION_20261020_EN.md` and pushed in commit `aac6302` before calculation. The input was the provisional 472-candidate private bundle (SHA-256 `1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c`); the audit script SHA-256 was `dfe7536c9c5dc1b60ebc35edb7e56392302c496fdb5929879c4d08e2444371df`. The aggregate output is `private_data/ppmi_plate_metadata_diagnostic_20261020.json`, which is excluded from Git.

| Cohort and technical field | PD | Controls | In-sample prevalence-only AUROC |
| --- | ---: | ---: | ---: |
| PPMI, release RNA-seq plate | 311 | 161 | 0.577670 |
| GENEPARK, recorded batch | 120 | 119 | 0.751576 |

For PPMI, the 50,071 PD–control pairs yielded a twice-favorable count of 57,849 and an AUROC of 57,849 / 100,142 = 0.5776697090. The GENEPARK comparison was calculated by the same exact tied-pair rule from the already saved `results/geo_technical_metadata_audit_20261019/class_by_batch.csv` counts, without revisiting expression or fitting a model. Both numbers use their own full-cohort labels to construct prevalence scores and are optimistic descriptive associations. GENEPARK's separately reported outer-fold batch-only AUROC was 0.691; it is not the in-sample number in this table.

All 472 candidate rows matched the deterministic `N` (release) to `pN` (portal) normalization. The 29 observed raw plate pairs were one-to-one in both directions. The release `Plate` remains the reference technical variable. Twenty-eight plates contain both classes; their PD proportions range from 7/13 (53.8%) to 12/16 (75.0%). Plate `201` ↔ `p201` has one provisional PD candidate and no control, was retained, and therefore has 100% PD; the 53.8–75.0% range does not include it.

The locally supplied two-page Project 133 IR3/B38 long-RNA methods PDF (`PPMI_RNAseq_IR3_Release_20210402.pdf`) and 92-page PPMI Data User Guide (`PPMI_Data_User_Guide_20260806.pdf`) contain no explanation found for a `201` or `p201` plate designation. Its purpose is unknown; a rerun or reprocessing label is only a hypothesis and was not used to change eligibility.

The recorded plate carries substantially less diagnosis discrimination in these PPMI candidates than the recorded batch does in GENEPARK. This narrows the specific label–plate shortcut concern; it does not establish that PPMI is free from batch effects or that the final PPMI population has the same numbers. No expression-derived panel, sex, age, or disease model score was computed. Under the current execution decision, all remaining PPMI checks and, after the second public lock, scoring occur locally on the researcher's personal Mac. Pre-QC group denominators, medication and blood-draw timing, the final GEO-trained sex artifact, final eligibility, and that public lock remain pending.
