# GEO-only Hetzner execution handoff

Prepared 2026-10-18. No server has been purchased, accessed or benchmarked by
this package. The combined amendment remains unpublished. The archive contains
public GEO study inputs and previously computed GEO reference outputs only.
It contains **no PPMI participant data**. Do not copy the whole Mac workspace.

## 1. Create the server yourself in Hetzner Console

The supplied screenshots show CPX62: 16 AMD vCPUs, 32 GB RAM and 640 GB disk,
Nuremberg, Ubuntu 26.04, IPv4. CPX resources are shared, not dedicated physical
cores. Eight workers with two numerical threads is a candidate ceiling, not a
guarantee of eightfold acceleration. Dedicated-resource CCX is an alternative
if consistent sustained CPU throughput is important; compare the live quote.
The screenshot quote is EUR 0.254/hour for the server and EUR 159.20/month
including IPv4 and displayed VAT. Treat the live checkout as authoritative;
the quote is not a project runtime or total cost estimate.

Choose Linux x86/AMD, not ARM/CAX. Keep the selected Ubuntu image. Setup installs
Python 3.11.9 separately and never substitutes the system Python for the lock.
No GPU, extra volume, placement group or cloud-init content is required.
Select one server, name it `pdgeo-permutation`, and keep IPv4 for this guide.

On the Mac, generate a dedicated key only if this path is unused:

```bash
test -f "$HOME/.ssh/pdgeo_hetzner" || ssh-keygen -t ed25519 -f "$HOME/.ssh/pdgeo_hetzner" -C pdgeo-hetzner
cat "$HOME/.ssh/pdgeo_hetzner.pub"
```

Use **Add SSH key** in the server form and paste only the `.pub` content. Keep
the private key and passphrase local. Add a firewall permitting inbound TCP 22
from your current public IP; retain the default outbound access for downloads.
You can leave paid backups unchecked for this reproducible temporary worker,
provided you retrieve outputs before deleting it. Click **Create & Buy now**
yourself after reviewing the quote. Copy the resulting IPv4 address.

## 2. Transfer this single archive from the Mac

Run in the Mac Terminal, replacing `SERVER_IPV4` with the actual address:

```bash
server_ip='SERVER_IPV4'
bundle='/Users/mohsenshamsijazeb/Documents/ChatGPT/Potential paper/pd-blood-expression-transportability/private_data/server_transfer/GEO_SERVER_BUNDLE.tar.gz'
scp -i "$HOME/.ssh/pdgeo_hetzner" "$bundle" "$bundle.sha256" "root@$server_ip:/root/"
ssh -i "$HOME/.ssh/pdgeo_hetzner" "root@$server_ip"
```

Verify the server's host key through the provider console on first connection.
The transfer uses an explicit file allowlist: four frozen prepared GEO files,
ten folds' original GEO predictions/decisions, and the minimal execution source,
lock, configuration and guide. It excludes the repository's PPMI loaders/tables,
private caches, credentials, raw archives, models and unrelated workspace files.

## 3. Install and validate on the server

Commands below run **on the Linux server**, after SSH login:

```bash
cd /root
sha256sum --check GEO_SERVER_BUNDLE.tar.gz.sha256
tar -xzf GEO_SERVER_BUNDLE.tar.gz
cd /root/pdgeo
apt-get update
apt-get install -y curl ca-certificates tmux libgomp1 libstdc++6
curl -LsSf https://astral.sh/uv/install.sh -o /root/install-uv.sh
sh /root/install-uv.sh
export PATH="$HOME/.local/bin:$PATH"
bash pipeline/server_setup.sh
```

All installed scientific package versions must exactly match the existing lock.
If a locked version has no compatible Linux distribution, **stop and return the
error**. Do not upgrade/downgrade packages, loosen pins or regenerate the lock.
Linux wheel availability is not established by the Mac tests.

## 4. Preflight only: synthetic benchmark and observed identity replay

```bash
tmux new -s pdgeo
bash pipeline/server_stage.sh preflight 2>&1 | tee preflight.log
```

Detach without stopping work: press Ctrl-B, release, then D. Reconnect with
`tmux attach -t pdgeo`. Log inspection: `tail -n 30 /root/pdgeo/preflight.log`.

This step performs synthetic checks, a full 16,364-input/1,000-selected-gene
stress fixture at one and eight workers, and a real-GEO **observed-label identity
replay**. It is not a randomized GEO permutation and does not touch PPMI.
All eight spawned workers must be observed. KKT, selected candidate decisions,
probabilities and frozen membership must meet the existing tolerances.

Linux RAM comes from `/proc/meminfo`; the swap gate uses the change in
`pswpout` pages in `/proc/vmstat`. `pgpgout` is diagnostic writeout in KiB,
not evidence of swap growth. macOS retains its separate `Swapouts`/`Pageouts`.
Missing counters fail closed. Keep the existing conservative memory budget
(worker high-water RSS sum + parent + 1 GiB <= 75% of RAM) and zero swap-growth
gate. This host-level check cannot attribute activity to a particular process.
No swap is created or disabled to manufacture a pass. Machine/OS/CPU/package
information is recorded without IPs or credentials. Run directly on the VM,
not inside a container whose resource limits differ from host `/proc` totals.

On success, return only `results/cloud_preflight/PREFLIGHT_RETURN.tar.gz`.
From a **new Mac Terminal**:

```bash
server_ip='SERVER_IPV4'
scp -i "$HOME/.ssh/pdgeo_hetzner" "root@$server_ip:/root/pdgeo/results/cloud_preflight/PREFLIGHT_RETURN.tar.gz" "$HOME/Downloads/"
```

If it fails, return `preflight.log` and the latest
`results/cloud_preflight/synthetic_*.json`; retain failed attempts. Do not run
`pilot` or `run`. Cross-OS numerical identity is not presumed to pass.

## 5. Publication gate, then one runtime pilot

After inspecting the server evidence, publish the combined amendment and the
exact passing certificate/config/code binding on the existing Zenodo record.
Only then install the verified matching publication receipt at
`results/cloud_preflight/ZENODO_PUBLICATION_RECEIPT.json`.
The receipt needs `publicly_verified: true`, the public `/records/` URL,
`binding_sha256` and `preflight_report_sha256`. These must be verified from the
actual public release, not filled in as a workaround for the gate.

After that gate, inside tmux:

```bash
bash pipeline/server_stage.sh pilot 1 2>&1 | tee runtime_pilot.log
```

This runs registered replicate **0**, at one worker/two threads, preserving its
completion/checksums for reuse in B=200. It does not compute an inferential
p-value. Return `private_data/geo_permutations_cloud_v1/RUNTIME_PILOT.json` for
runtime review (GEO aggregate results only). Do not discard replicate 0, pick
another seed or shorten B according to its values. One serial replicate measures
serial cost, not eight-worker throughput. Concurrency can be reduced to 1, 2 or
4 within the tested ceiling of 8, solely for resource/runtime reasons; all seeds,
fits, configuration binding and statistics remain unchanged.

## 6. Full run and recovery

After reviewing pilot runtime/resource use, start the fixed-B run inside tmux:

```bash
bash pipeline/server_stage.sh run 8 2>&1 | tee permutation.log
```

The same command resumes verified completed replicates after an interruption.
Keep the output directory and original attempts. Incomplete fits restart cold
with the same seed; terminal numerical failures remain invalid and cannot be
retried as interruptions. Do not delete lock/failure files to force progress.
Summary p-values are emitted only after all 200 valid completions. The serial
61.27-hour estimate from prior observed runs is only a planning baseline; null
tuning, shared CPU contention and parallel overhead can change it substantially.

Retrieve `SUMMARY.json`, logs and the complete checkpoint directory to durable
storage before deleting the server. Do not publish participant-level outputs.
Powering off a Hetzner VM does not stop billing: delete it only after retrieving
and verifying the results, and delete an unneeded billable Primary IPv4 too.
Deletion is a separate user action; these scripts never delete infrastructure.

## Read-only progress display

`pipeline/monitor_permutations.py` is a standalone standard-library observer,
outside the frozen scientific code binding. Run it with `python3` and
`--interval 10` (or 30). It never starts fits, writes checkpoints, changes labels,
or modifies the executor. It excludes the observed identity replay from the
replicate count and includes registered replicate 0 only once. It displays
completed outer folds, committed replicate markers, live replicate locks and
terminal failure records. These are display counts; the executor remains the
authority for checkpoint checksum validation and final completeness.

Fold percentages count equally sized work units, not equal runtime. Pilot mode
shows its ten folds separately. Approximate full-run ETA is displayed only after
at least two newly completed replicates in the current controller, excluding
resumed checkpoints. A partial fold is not given a guessed runtime percentage.
Ctrl-C ends the observer only. Successful completion requires B=200 committed
replicates and the final complete summary; 200 markers alone are not sufficient.

## Sources

### HTTPS transfer fallback (18 October 2026)

After the user created the VM and logged in, direct SCP of the 37 MiB archive
stalled at a few KiB/s. An encrypted transport asset is available at the GitHub
prerelease `geo-server-transfer-20261018`. This is a transport artifact, not a
scientific release, preregistration, or amendment. It contains the exact original
archive encrypted with OpenSSL AES-256-CBC, PBKDF2 and 200,000 iterations.

Download `GEO_SERVER_BUNDLE.tar.gz.enc` and the original `.tar.gz.sha256` from
that release directly on the VM. Transfer only the 65-byte local
`private_data/server_transfer/GEO_SERVER_BUNDLE.transfer-pass` separately over
authenticated SSH. Never publish that password file or commit it to Git.
Decrypt using `openssl enc -d -aes-256-cbc -pbkdf2 -iter 200000` with
`-pass file:/root/GEO_SERVER_BUNDLE.transfer-pass`, then verify the original
archive checksum before extraction. A partially uploaded unencrypted archive
must not be used. All publication and preflight gates above still apply.

- [Hetzner resource types](https://docs.hetzner.com/cloud/servers/faq/)
- [Hetzner SSH connection](https://docs.hetzner.com/cloud/servers/getting-started/connecting-to-the-server/)
- [Hetzner billing](https://docs.hetzner.com/cloud/billing/faq/)
- [uv Python management](https://docs.astral.sh/uv/guides/install-python/)
- [Linux proc documentation](https://docs.kernel.org/filesystems/proc.html)
