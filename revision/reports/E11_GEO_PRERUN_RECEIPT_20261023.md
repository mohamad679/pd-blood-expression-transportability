# E11 GEO-only audit: pre-run local hash receipt

Recorded before the first E11 real-data execution on 23 October 2026. This is a **local hash record**, not an independent timestamp, public registration, or a modification of the sealed Phase-2 amendment. The post hoc analysis was specified after the E5 results had been seen and after E6 completed. No PPMI outcome result is used to choose the E11 conditions.

- E11 addendum SHA-256: `20ee61ba3f403af020b23897795863978eb5b9c153b5b4f9bc298f82a07f881f`
- E11 executable SHA-256: `d5c94a84d2ee006ccbb13303054bc32fd25f19ef95dc42902f6c0e06f7d4fc01`
- Imported fitting-code SHA-256 values: `revision/e03_ppmi_nested_cv.py` `0f6b8eb24aa81158652d769fc39495ed3163241904b6d2428969833d68816157`; `revision/e04_transductive_universe.py` `1f40ac5f0288f50e2fbc5e3bea85765f1c6342c5fc8060c5763b9c479e30d538`; `revision/e05_weighted_deviance_panel_exclusion.py` `bfb2df0f4cb94af0df2e634e53ee278c4bbd1f14c5ffab5e4d1572c7cc784c32`; `pipeline/geo_pipeline.py` `6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`.
- Requirement lock SHA-256: `dd6d9fc3d9b48da6740f63b8218fa1086f25b325fc2d832c3d8f3da218b9da00`.
- Sealed Phase-2 amendment SHA-256, unchanged: `6947539f57182a92a75bc7ecc71a267a2184e5a2013875342a6156db89d30d39`
- Execution environment: `uv run --offline --no-project --python 3.11 --with-requirements requirements_glum_v3.lock.txt`
- Planned GEO-only output: `results/geo_tuning_seed_stability_20261023/`
- The program imports E5 fitting helpers but does not load or score any PPMI data. The source GEO input hashes will be written to its aggregate output.
- A no-data aggregate self-test passed before this receipt; it checked the 20-seed no-gene counting path.

If the executable or addendum changes before or during execution, stop and create a new receipt. Report any deviation with its timing. The end result remains post hoc and descriptive.
