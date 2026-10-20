# Draft second pre-PPMI lock — not effective

**Status:** planning document only. No PPMI score is authorized by this draft. A completed, timestamped Zenodo version with a DOI and the matching exact Git commit/tag is required before the first PPMI cell-panel or model/control score.

## Current execution site

- The researcher's personal MacBook is the sole site for all **remaining** PPMI acquisition, local storage, eligibility and metadata checks, control/disease scoring, and inference.
- The existing provisional 472-candidate local bundle may be used for checks supported by its frozen contents. Obtain missing metadata directly into local, nonsynced Mac storage. Do not run further PPMI processing on another computer or server.
- Keep participant-level PPMI inputs and outputs under `~/PPMI_private` or an equivalently reviewed local path, outside Git, Zenodo, chat, and cloud-synced Documents. External scoring and inference require frozen GEO-trained artifacts and the published second lock.
- The current site decision is [documented separately](PPMI_MAC_ONLY_EXECUTION_DECISION_20261020_EN.md). Earlier acquisition and preparation records remain historical provenance and must be reported accurately; they are not instructions for future execution.

## Required fields before publication

| Item | Current status |
|---|---|
| Final PPMI eligibility and one-to-one crosswalk | Pending; provisional counts are not final |
| Plate × diagnosis metadata count | Completed for provisional 472 candidates; final population still pending |
| Existing input preparation code, release, output schema and SHA-256 | Recorded; verify binding to final local input |
| Local private bundle SHA-256 | Checked; `1b191f042885eedf2eea00c69ce79a269cc647ef4d5148b37b678b97c7d0333c` |
| Mac-only storage and applicable PPMI data-use terms | Local secure path established; terms and earlier sync history still under review |
| Mac scoring code, environment, exact command and portability checks | Pending |
| Frozen final GEO M1–M5 and preprocessing artifact hashes | Pending public manifest |
| Final GEO-trained auxiliary sex model hash | Fit and independently audited on GEO controls; SHA-256 `aa92a1516328477b1c4ba1a61e7bcec43263cc38574f57ec99c5f26465348200`. Include this exact artifact in the public lock manifest; no PPMI score has been run. |
| Control gates, bootstrap rule and private output location | Must be verified against protocol v3.1 |
| Exact Git commit/tag and public Zenodo version/DOI | Not yet published |

No field marked pending may be presented as satisfied. The final lock will reference private participant-level files only by non-identifying hashes and controlled local paths; it will not publish those files.
