#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")/.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export NUMEXPR_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 PYTHONUNBUFFERED=1
report="results/cloud_preflight"
run="private_data/geo_permutations_cloud_v1"
config="protocol/permutation_config_cloud_v1.json"
py=".venv/bin/python"
mkdir -p "$report"
case "${1:-}" in
  preflight)
    if [[ -e "$report/ZENODO_PUBLICATION_RECEIPT.json" ]]; then
      printf '%s\n' 'Receipt exists: preserve registered evidence; do not rerun preflight.' >&2
      exit 1
    fi
    # Timestamped attempts retain earlier failures; no real randomized fit here.
    stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    synthetic="$report/synthetic_${stamp}.json"
    "$py" -m pipeline.verify_permutation_executor --config "$config" \
      --output "$synthetic" --work private_data/cloud_synthetic --benchmark
    "$py" -m pipeline.permutation_executor plan --config "$config" --output "$run"
    "$py" -m pipeline.permutation_executor identity --config "$config" --output "$run"
    "$py" -m pipeline.finalize_permutation_preflight --config "$config" \
      --synthetic "$synthetic" --run-root "$run" --output "$report/preflight_passed.json"
    tar -czf "$report/PREFLIGHT_RETURN.tar.gz" "$synthetic" \
      "$report/preflight_passed.json" "$run/PLAN.json" "$run/IDENTITY_VERIFIED.json"
    printf '%s\n' 'STOP: return results/cloud_preflight/PREFLIGHT_RETURN.tar.gz. Amendment publication is still required.'
    ;;
  pilot|run)
    "$py" -m pipeline.permutation_executor "$1" --config "$config" --output "$run" \
      --receipt "$report/ZENODO_PUBLICATION_RECEIPT.json" \
      --preflight "$report/preflight_passed.json" --workers "${2:-8}"
    ;;
  *) printf '%s\n' 'Usage: bash pipeline/server_stage.sh preflight|pilot|run [workers: 1,2,4,8]' >&2; exit 2 ;;
esac
