#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$project_dir"
report_dir="results/permutation_preflight_20261018"
private_dir="private_data/permutation_preflight_20261018"
run_dir="private_data/geo_permutations_20261018"
mkdir -p "$report_dir" "$private_dir"

# Uses the previously successful index without changing the locked versions.
runner=(uv run --default-index https://pypi.tuna.tsinghua.edu.cn/simple --with-requirements requirements_glum_v3.lock.txt)
"${runner[@]}" python -m pipeline.verify_permutation_executor \
  --output "$report_dir/synthetic_preflight_final.json" --work "$private_dir" --benchmark
"${runner[@]}" python -m pipeline.permutation_executor plan --output "$run_dir"
"${runner[@]}" python -m pipeline.permutation_executor identity --output "$run_dir"
"${runner[@]}" python -m pipeline.finalize_permutation_preflight \
  --synthetic "$report_dir/synthetic_preflight_final.json" --run-root "$run_dir" \
  --output "$report_dir/preflight_passed.json"

printf '%s\n' 'Preflight finished. Real permutations still require the public amendment receipt.'
