# GEO primary analysis engine

Public implementation amendment: https://doi.org/10.5281/zenodo.23405802
Parent preregistration: https://doi.org/10.5281/zenodo.23405371

`geo_pipeline.py` loads only the two prepared GEO sources. It has no PPMI loader
or external-scoring command. Every fit uses `numerical_kkt_v3_1.fit_checked`.
All learned preprocessing is training-only; all model fits start cold.

The public supplement ZIP fixes the code/dependencies, exact split manifest and
input hashes used for the pilot. Code inside that published ZIP must not be
overwritten to make an old hash appear current. Future fixes require a new
documented version before affected re-execution.

Run from the workspace root, with the audited environment:

```sh
uv run --offline --no-project --python 3.11 \
  --default-index https://pypi.tuna.tsinghua.edu.cn/simple \
  --with glum==3.4.1 python \
  pd_shortcut_feasibility_2026-10-13/pipeline/verify_pipeline.py \
  --output /tmp/geo_synthetic_checks.json
```

`plan` validates local inputs and generates donor-level splits without fitting.
`pilot` executes only the first locked PROBE outer fold, with all five models.
`analysis --scope SCOPE` executes one exact planned outer, transfer or final
training analysis from the split manifest. Use a fresh output directory for each
execution. Every real fit requires the locally recorded, publicly verified
implementation-amendment receipt; code/split hashes must agree.

```sh
uv run --offline --no-project --python 3.11 \
  --default-index https://pypi.tuna.tsinghua.edu.cn/simple \
  --with glum==3.4.1 python \
  pd_shortcut_feasibility_2026-10-13/pipeline/geo_pipeline.py pilot \
  --output pd_shortcut_feasibility_2026-10-13/pipeline/NEW_RUN_DIRECTORY \
  --registration-receipt \
  pd_shortcut_feasibility_2026-10-13/protocol/ZENODO_REGISTRATION_RECEIPT_v3_1_1.json
```

Outputs retain preprocessing, final model coefficients, candidate decisions,
per-attempt numerical checks/warnings/timing, and held-out predictions. An
exception writes `FAILURE.json`; no approximate candidate is silently skipped.
The timing pilot is part of the planned analysis, not independent replication
or a basis to select a new method from its disease performance.

This delivery covers the primary M1-M5 development engine and paired bootstrap
algebra. Registered sex/age auxiliary controls, sensitivities, permutation runs,
ComBat diagnostics and bootstrap coverage validation remain separate tasks.
Before PPMI disease scoring, the final model/preprocessing/population second
public lock must be completed; the implementation DOI is not that final lock.
