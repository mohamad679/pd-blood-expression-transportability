# Draft second pre-PPMI lock — not effective

**Status:** planning document only. No PPMI score is authorized by this draft. A completed, timestamped Zenodo version with a DOI and the matching exact Git commit/tag is required before the first PPMI cell-panel or model/control score.

## Frozen execution topology

- Source archive: researcher's Windows computer in Italy; the approximately 42 GB archive remains there.
- Pre-score preparation: run the frozen extraction script on Windows; produce a compact private PPMI bundle containing only required metadata and count/length-corrected, within-sample-ranked features. No panel, sex, age, or disease model score is computed on Windows.
- Transfer: copy the compact bundle to the researcher's personal Mac; verify its SHA-256 against the Windows manifest. The researcher supplies its local path for the Mac-side analysis. Do not transfer the full archive unless a separate documented need arises.
- External scoring and inference: run locally on that Mac using frozen GEO-trained artifacts, after this draft has been completed and publicly locked. No PPMI participant-level data are included in the public lock.

The detailed binding decision is in `docs/PPMI_SECOND_LOCK_EXECUTION_DECISION_20261020_EN.md`. The original Windows-scoring proposal was rejected before scoring and is superseded.

## Required fields before publication

| Item | Current status |
|---|---|
| Final PPMI eligibility and one-to-one crosswalk | Pending; provisional counts are not final |
| Plate × diagnosis metadata count | Pending |
| Windows extraction code, input release, output schema and SHA-256 | Pending |
| Compact private bundle SHA-256 and Mac transfer verification | Pending |
| Mac scoring code, environment, exact command and portability checks | Pending |
| Frozen final GEO M1–M5 and preprocessing artifact hashes | Pending public manifest |
| Final GEO-trained auxiliary sex model hash | Pending fit/lock |
| Control gates, bootstrap rule and private output location | Must be verified against protocol v3.1 |
| Exact Git commit/tag and public Zenodo version/DOI | Not yet published |

No field marked pending may be presented as satisfied. The final lock will reference private participant-level files only by non-identifying hashes and controlled local paths; it will not publish those files.
