# PPMI input inventory through external evaluation

**Current plan, 2026-10-20:** The researcher's personal Mac is the site for all remaining PPMI processing. The proposed new Windows one-shot export was [superseded before execution](PPMI_ONE_TIME_CACHE_EXPORT_DECISION_20261020_EN.md). This inventory consolidates remaining needs; it does not authorize scoring before the published second lock.

| Need | Current source | Action |
|---|---|---|
| Group-specific pre-QC denominators | Complete `metaDataIR3` is not in the current Mac candidate bundle or supplied aggregate. Overall release QC is 4,756 pass / 115 fail among 4,871 rows; baseline is 1,530 pass / 47 fail among 1,577 rows. | Look for a separately downloadable small file through approved LONI access on the Mac. If absent, report group-specific counts as unavailable; descriptive only and not a lock gate. Do not re-download the large archive. |
| Complete portal, baseline age and participant-status exports | Portal and age copies are already on the Mac; previous audits used participant status. | Use existing local versions and hashes; do not request them again for descriptive QC if release metadata remains unavailable. |
| Frozen 472 baseline candidates, gene manifest and within-sample ranks | Existing hashed private candidate bundle on the Mac. | Verify final eligibility, input hashes, and scorer compatibility locally before lock. |
| Nonbaseline release sample ranks for longitudinal similarity | Not supplied by the 472-candidate baseline bundle. | Record as unavailable for now. A later separate diagnostic may be deferred; it cannot block the primary PPMI score or trigger Windows processing or a large redownload. |
| Baseline blood-collection month | Existing whole-blood collection export on the Mac | Descriptive timing audit complete; 93 candidate months unknown |
| Recorded dopaminergic medication timing | Existing LEDD export on the Mac | Descriptive audit complete; no candidate exclusion. Register the conditional GEO-to-PPMI treatment caution in the second lock before scoring. |
| CBC marker-control measurements | Existing current and archived CBC exports on the Mac | Post-lock fixed-score association audit; no new source request presently identified |
| Final M1–M5 and auxiliary control models, frozen transformations, inference code | GEO project artifacts/code on the Mac | Finish code/model manifest and synthetic/GEO parity checks before second lock |
| DPC publication review and acknowledgement | PPMI administrative policy | Separate publication task, not an additional RNA data table |

Any unexpected source-version or integrity failure should be reported as one consolidated validation finding, rather than converted into serial source requests or an ad hoc change to the external cohort.
