# GEO pipeline

The pipeline reads only the two locally prepared GEO development matrices from `prepared/`. It does not load PPMI data and exposes no PPMI scoring command. The exact model and preprocessing rules are in the frozen English reference under `protocol/`.

## Local setup

Use Python 3.11 and install the pinned packages from the repository-root `requirements_glum_v3.lock.txt`. Keep `prepared/` untracked; the expected file names, checks, and immutable hashes are documented in `data/README.md` and `data/manifests/input_audit.json`.

Run the synthetic software checks without real study inputs:

```bash
python pipeline/verify_pipeline.py --output /tmp/pd_synthetic_checks.json
```

Generate/check the split manifest from local GEO inputs:

```bash
python pipeline/geo_pipeline.py plan --output /tmp/pd_geo_plan
```

Any real model fit requires a local registration receipt that confirms the exact code and split hashes are publicly registered. The receipt is not included here and must never be fabricated. The preregistered one-fold pilot command, if the receipt matches, is:

```bash
python pipeline/geo_pipeline.py pilot \
  --output /path/to/fresh/local/run_directory \
  --registration-receipt /path/to/local/public_registration_receipt.json
```

Run an exact planned analysis scope only after confirming its name in `data/manifests/splits.json`:

```bash
python pipeline/geo_pipeline.py analysis \
  --scope 'EXACT_SCOPE_FROM_MANIFEST' \
  --output /path/to/fresh/local/run_directory \
  --registration-receipt /path/to/local/public_registration_receipt.json
```

Every run requires a new output directory. The pipeline rejects appending to an earlier run, verifies the registration receipt against current code and split hashes, and writes a failure record rather than silently dropping failed candidates. Do not upload predictions, per-fit logs, preprocessing arrays, or coefficient files to GitHub.

## Current execution status

The original locked pilot `GSE57475/outer/0` completed on 2026-10-16. Its aggregate report is in `results/geo_pilot/`. No full nested CV or external PPMI score has been run.
