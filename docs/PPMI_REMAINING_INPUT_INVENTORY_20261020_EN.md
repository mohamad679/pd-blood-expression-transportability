# PPMI input inventory through external evaluation

This inventory is a single request plan, not an authorization to score before the second public lock. It separates source data from code/model work so the researcher is not asked for one missing table at a time.

| Need | Current source | Action |
|---|---|---|
| Complete release rows, including `QCflagIR3`, for group-specific pre-QC denominators | Already inside the original `IR3_PRIVATE_CACHE.zip` on the researcher's Windows computer | One-time cache script writes aggregate QC audit and private metadata capsule |
| Complete portal, baseline age and participant-status exports, plus count-file linkage index | Already inside that same cache; portal and age copies also on the Mac | Same one-time script captures all versions and hashes together |
| Frozen 472 baseline candidates, gene manifest and within-sample ranks | Existing hashed private candidate bundle on the Mac | Reconcile exact selection against the cache export; no new download |
| Nonbaseline release sample ranks for longitudinal similarity | Count vectors already inside the same cache | Same one-time script prepares all 4,871 within-sample rank vectors as a private all-visit bundle; no similarity score before lock |
| Baseline blood-collection month | Existing whole-blood collection export on the Mac | Descriptive timing audit complete; 93 candidate months unknown |
| Recorded dopaminergic medication timing | Existing LEDD export on the Mac | Descriptive audit complete; no candidate exclusion |
| CBC marker-control measurements | Existing current and archived CBC exports on the Mac | Post-lock fixed-score association audit; no new source request presently identified |
| Final M1–M5 and auxiliary control models, frozen transformations, inference code | GEO project artifacts/code on the Mac | Finish code/model manifest and synthetic/GEO parity checks before second lock |
| DPC publication review and acknowledgement | PPMI administrative policy | Separate publication task, not an additional RNA data table |

The Windows run produces only the files named in [the one-time decision](PPMI_ONE_TIME_CACHE_EXPORT_DECISION_20261020_EN.md). Share its aggregate JSON; transfer the private bundles directly to local Mac storage and provide only their paths. Once the capsule and rank bundle pass hash/reconciliation checks on the Mac, the currently known PPMI source-input list is closed. Any unexpected source-version or integrity failure must be reported as one consolidated validation finding, rather than converted into serial requests or an ad hoc change to the external cohort.
