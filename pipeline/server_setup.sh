#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")/.."
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  printf '%s\n' 'This setup targets Linux x86_64 only.' >&2
  exit 1
fi
sha256sum --check TRANSFER_SHA256.txt
if ! command -v uv >/dev/null 2>&1; then
  printf '%s\n' 'Install uv using the documented command, then rerun setup.' >&2
  exit 1
fi
uv python install 3.11.9
uv venv --python 3.11.9 .venv
# Fail on unavailable exact versions; never regenerate or loosen the lock.
uv pip sync --python .venv/bin/python requirements_glum_v3.lock.txt
.venv/bin/python -c 'from pipeline.permutation_executor import locked_environment; print(locked_environment())'
printf '%s\n' 'Exact environment installed. Next: bash pipeline/server_stage.sh preflight'
