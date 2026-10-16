# Data inputs and handling

This repository contains manifests and aggregate QC summaries only. It intentionally excludes raw and processed expression matrices, sample-level phenotypes, individual PPMI metadata, private IDA crosswalks, and external held-out predictions.

## GEO development inputs

The execution pipeline expects the following files under a local, untracked `prepared/` directory at the repository root:

- `GSE57475_development.npz`
- `GSE57475_sample_metadata.json`
- `GSE99039_development.npz`
- `GSE99039_sample_metadata.json`

The required fields and validation checks are enforced by `src/pipeline/geo_pipeline.py`. The immutable input hashes and planned sample counts are recorded in `data/manifests/input_audit.json`. The actual expression inputs are public GEO data, but are omitted from this repository to keep the code repository small and to avoid treating derived matrices as source-controlled artifacts. Use the established project preparation procedure and re-check hashes before running.

## PPMI external data

PPMI controlled-use files and subject-level metadata must stay in the approved local environment under the PPMI data-use agreement. Do not commit, upload, publicly post, or use the data to attempt re-identification. The final population and count-file crosswalk are not finalized in this repository. A public aggregate summary may be added only after checking it contains no subject-level values or identifying information and is permitted by the data-use terms.

## Pilot output policy

Only aggregate pilot metrics, runtime, public input/split hashes, and non-identifying audit summaries are tracked. The following remain local and are ignored: held-out predictions, fitted coefficients, preprocessing arrays, per-fit JSONL logs, and any raw or prepared sample-level data.
