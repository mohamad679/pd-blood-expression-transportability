# Linux cloud handoff — 18 October 2026

## Changes

- Added Linux `/proc/meminfo` RAM and `/proc/vmstat` swap counters. `pswpout`
  controls the swap gate; `pgpgout` is diagnostic. Missing counters cannot pass.
  The existing memory and zero-swap-growth acceptance rules were not relaxed.
- Added a separate eight-worker/two-thread cloud config; scientific settings,
  seeds, B=200, frozen GEO engine and the original Mac config remain unchanged.
  Stress benchmarking now supplies at least one full-dimensional job per worker
  and compares the same job IDs against the serial arm.
- Record machine characteristics in the execution binding and require synthetic,
  identity and certificate stages to describe the same machine.
- Added a receipt-gated runtime pilot: registered replicate 0, retained in the
  full B=200 run, without an interim p-value. Runtime concurrency may be reduced
  within the tested ceiling without changing seeds, methods or denominator.
- Added Bash setup/preflight/pilot/run commands, tmux instructions and an
  explicit-allowlist transfer archive builder. Only the four prepared GEO inputs,
  ten GEO reference folds and required execution files are packaged. PPMI
  participant data, loaders, archives, caches and credentials are excluded.
- Updated the unpublished combined amendment draft to describe cloud execution.

## Local validation

Six portability tests passed using synthetic Linux `/proc` fixtures and the
macOS counter path, including missing-counter/memory cases. Bash syntax checks
passed. The final local synthetic executor regression passed all **27** checks,
with **650** numerical fit attempts and maximum KKT residual
**2.6839175860937647e-8**. It covered seeds, repeated-donor labels, rebuilt inner
folds, exact M4/M5 agreement with the frozen engine, worker/order equivalence,
verified resume, corruption refusal and terminal-failure handling.

The local regression deliberately omitted the resource benchmark: it does not
certify eight workers, Linux package installation, Linux numerical identity or
server speed. Earlier intermediate local reports are retained under the results
directory and are not used as final release certificates. The frozen GEO engine
still has SHA256
`6f7f6cf42608002b185a63044387c9b9f5e624fc0b9ff5ae156030e3655e1a2c`.

An isolated extraction of the transfer archive verified all **39** payload
hashes and ran its GEO `plan` command without any model fit: **381** records,
**16,364** genes, **10** outer scopes and the unchanged split manifest hash
`2f72dfc7105aa19bdce3fc8fd40f8b6bfccf52cd2ba500347eb7698f912077c2`.
The archive is approximately 38.5 MB and includes a full payload checksum list
and source-commit provenance. The transfer archive stays ignored locally.

## Remaining gate and next action

Follow [the server runbook](HETZNER_GEO_RUNBOOK_EN.md), transfer the single archive,
install the exact lock and run server **preflight only**. Return its aggregate
`PREFLIGHT_RETURN.tar.gz`, or the failure log/report. No server was purchased or
accessed, no eight-worker/server benchmark was performed, no amendment was
published, and no real-GEO randomized fit or PPMI score was executed here.

Only after the server checks and observed-label identity pass may the matching
combined amendment/certificate be published. Then run the retained replicate-0
runtime pilot, review its cost and launch the fixed-B resumable run. A CPX VM has
shared vCPUs; the screenshots' specifications are not a throughput guarantee.
