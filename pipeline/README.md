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

All 13 core GEO scopes completed on 2026-10-17/18, including the original
`GSE57475/outer/0` pilot from 2026-10-16. The pilot is included once in the
five-fold estimate. Core aggregates are in `results/geo_development_20261018/`.
No PPMI score has been run.

## Auxiliary controls and numerical diagnostic

`geo_controls.py` implements the registered sex/age transfer controls in healthy
GEO controls. Its implementation details were committed before execution in
`docs/GEO_CONTROL_EXECUTION_PLAN_20261018_EN.md`. Both sex transfer directions
passed; the age control did not establish development capacity. See
`docs/GEO_TRANSFER_CONTROLS_20261018_EN.md` and the aggregate under
`results/geo_controls_20261018/`. Use fresh, ignored private output directories.

```bash
python pipeline/verify_geo_controls.py --output /private/tmp/geo_aux_synthetic.json
python pipeline/geo_controls.py \
  --output private_data/geo_controls_NEW_RUN \
  --aggregate /private/tmp/geo_controls_NEW_RUN.json
```

`audit_sigmoid.py` is a retrospective replay diagnostic that leaves the locked
core pipeline and its completed results unchanged. The bounded run replays all
candidates of the six model/scope contexts in which warnings were reported;
it is not a new disease analysis or permission to replace tuning scores. See
`docs/SIGMOID_REPLAY_SCOPE_20261018_EN.md` for the scope and retained attempts.
The completed 490-fit diagnostic reproduced all six reported warnings and all
original candidate means and selections. Explicit clipping changed no AUROC.
Results and limitations are in `docs/SIGMOID_AUDIT_20261018_EN.md`.
