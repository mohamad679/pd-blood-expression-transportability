# Windows PPMI preparation (before the second public lock)

This tool reads the **existing** `IR3_PRIVATE_CACHE.zip` made by the prior complete archive audit. It does not reread the 42 GB archive. It selects the 472 **provisional** baseline Long RNA candidates using the same exact linkage and metadata rules as the earlier audit, applies fixed exon-union-length division and within-sample fractional midranks, and writes a compact **private** bundle for the researcher's Mac. It also writes a separate aggregate-only plate-by-group report. It does not fit a model or compute cell-marker, age, sex, or PD scores. Final eligibility and the second public lock remain pending.

Required: Python 3.10 or newer on Windows, and the complete `IR3_PRIVATE_CACHE.zip` from the previous audit. No Python packages or internet download are required. The `PPMI_WINDOWS_PREP_TOOLS.zip` file includes all three `.py` files with frozen audit dependencies.

1. Copy `PPMI_WINDOWS_PREP_TOOLS.zip` from the Mac to the Windows computer's Downloads directory.
2. In Windows PowerShell, run:

```powershell
$tools = 'D:\IR3audit\PPMI_Prepare'
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\PPMI_WINDOWS_PREP_TOOLS.zip" -DestinationPath $tools -Force
py -3 "$tools\prepare_ppmi_windows_bundle.py" --self-test
Test-Path 'D:\IR3audit\AllOnce\IR3_PRIVATE_CACHE.zip'
```

If `Test-Path` prints `False`, locate the cache and substitute its true path in the command below. If `True` and `SELF_TEST_PASSED` appear, run:

```powershell
py -3 "$tools\prepare_ppmi_windows_bundle.py" --cache 'D:\IR3audit\AllOnce\IR3_PRIVATE_CACHE.zip' --output-dir 'D:\IR3audit\PPMI_Prepared_20261020'
```

At completion, the progress line reaches `472/472` and `PRIVATE_BUNDLE_READY` is printed. Send **only** `D:\IR3audit\PPMI_Prepared_20261020\PPMI_WINDOWS_PREP_AGGREGATE.json` to the analyst (or paste its contents). Keep `PPMI_WINDOWS_PRIVATE_BUNDLE.zip` private and local. For the later transfer to the personal Mac, copy that bundle through a private channel and verify its SHA-256 against `bundle_sha256` in the aggregate report. Never upload it to public GitHub, Zenodo, this chat, or the former Hetzner server.

If a command fails, send the error text, but do not transfer `PRIVATE_INPROGRESS` or a partial bundle. The script refuses to overwrite existing completed outputs. A fresh output directory can be used after resolving the error.
