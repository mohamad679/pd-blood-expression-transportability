PPMI one-time export from the existing IR3_PRIVATE_CACHE.zip
===========================================================

Purpose
-------
The script uses the existing cache on the researcher's Windows computer.
It does not need the original ~42-GB download and downloads nothing.
It writes one aggregate JSON plus two PRIVATE ZIPs. The private ZIPs must
be moved directly to the researcher's private Mac storage; never attach
them to chat, email them to the analyst, put them in Git, or publish them.
No PPMI model, panel, sex, age, PD, or longitudinal similarity score is run.

Prerequisites
-------------
* Python 3.10 or newer (the earlier Windows run used Python 3.11).
* Existing D:\IR3audit\AllOnce\IR3_PRIVATE_CACHE.zip.
* Enough free disk for a temporary ~319-MB rank matrix, the compressed
  rank ZIP, and the small metadata ZIP. Allow at least 1 GB free.
* A new, empty output folder name for each attempt.

Copy PPMI_ONE_SHOT_TOOLS_20261020.zip into the Windows Downloads folder,
then paste this ENTIRE block into PowerShell in order:

    $ErrorActionPreference = 'Stop'
    $package = "$env:USERPROFILE\Downloads\PPMI_ONE_SHOT_TOOLS_20261020.zip"
    $tools = 'D:\IR3audit\PPMI_OneShot_Tools_20261020'
    $cache = 'D:\IR3audit\AllOnce\IR3_PRIVATE_CACHE.zip'
    $out = 'D:\IR3audit\PPMI_OneShot_Output_20261020'
    if (-not (Test-Path -LiteralPath $package)) { throw 'Tool ZIP missing from Downloads' }
    if (-not (Test-Path -LiteralPath $cache)) { throw 'IR3 cache missing at the expected path' }
    Expand-Archive -LiteralPath $package -DestinationPath $tools -Force
    py -3 "$tools\prepare_ppmi_remaining_inputs.py" --self-test
    if ($LASTEXITCODE -ne 0) { throw 'Self-test failed; do not run the real cache' }
    py -3 "$tools\prepare_ppmi_remaining_inputs.py" --cache $cache --output-dir $out
    if ($LASTEXITCODE -ne 0) { throw 'Preparation failed; keep partial output private' }

If the frozen PPMI_WINDOWS_PRIVATE_BUNDLE.zip is still on Windows, add
--candidate-bundle 'D:\path\to\PPMI_WINDOWS_PRIVATE_BUNDLE.zip' to the real-run
command for an immediate exact-set check. The same check can be done on
the Mac later using its already stored bundle, so this option is not required.

The progress display has a byte-based bar and ETA while hashing the cache,
and a sample-based bar and ETA while preparing all-visit ranks. It updates
one line in a live terminal rather than printing a line for every sample.

After success
-------------
Send only:

    D:\IR3audit\PPMI_OneShot_Output_20261020\PPMI_COMPLETE_METADATA_AGGREGATE.json

Move these two PRIVATE files directly to ~/PPMI_private on the Mac and
tell the analyst only their Mac paths:

    PPMI_COMPLETE_METADATA_PRIVATE.zip
    PPMI_ALL_VISIT_RANKS_PRIVATE.zip

The aggregate JSON includes SHA-256 and byte size for each private ZIP;
verify those hashes after transfer. Neither private ZIP belongs in chat.

The private metadata ZIP contains complete release/portal/age/status rows,
the count-file linkage index, and the fixed gene reference, but no count
vectors. The all-visit ZIP contains prepared within-sample ranks for all
4,871 release samples. Its later longitudinal diagnostic still requires
the public second lock before any similarity score is calculated.

If the run fails, do not delete or overwrite the original cache. Use a
fresh output folder name for a retry and send only the generic error text
from the terminal. The script deliberately does not print private rows.
