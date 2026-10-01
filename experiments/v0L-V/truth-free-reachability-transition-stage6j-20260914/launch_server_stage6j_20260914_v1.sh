#!/usr/bin/env bash
set -euo pipefail

repo_root="${STAGE6J_REPO_ROOT:-/root/multiFSYNC-stage6b-20260903-fa10cbc}"
package_root="${STAGE6J_PACKAGE_SOURCE_ROOT:-/root/multiFSYNC-jaoua-paired-20260909}"
r_library="${STAGE6J_R_LIB:-/root/multiFSYNC-jaoua-paired-lib-20260909}"
stage6g_output="${STAGE6J_STAGE6G_OUTPUT_ROOT:-/root/multifsync-results/stage6g_20260911_v1}"
output_root="${STAGE6J_OUTPUT_ROOT:-/root/multifsync-results/stage6j_20260914_v1}"
workers="${STAGE6J_WORKERS:-8}"
experiment_root="$repo_root/experiments/v0L-V/truth-free-reachability-transition-stage6j-20260914"
runner="$experiment_root/run_stage6j_20260914_v1.R"

test -f "$runner"
test -d "$package_root/R"
test -d "$r_library"
test -f "$stage6g_output/STAGE6G_TRUTH_FREE_COMPLETE.txt"
test ! -e "$output_root"

export STAGE6J_PACKAGE_SOURCE_ROOT="$package_root"
export STAGE6J_R_LIB="$r_library"
export STAGE6J_STAGE6G_OUTPUT_ROOT="$stage6g_output"
export STAGE6J_OUTPUT_ROOT="$output_root"
export STAGE6J_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=preflight
Rscript "$runner" --action=prepare --output-root="$output_root"
Rscript "$runner" --action=smoke --output-root="$output_root"
Rscript "$runner" --action=run --output-root="$output_root"

printf 'status=STAGE6J_LAUNCH_COMPLETE\ncompleted_at=%s\ntruth_unsealed=FALSE\n' \
  "$(date --iso-8601=seconds)" > "$output_root/STAGE6J_LAUNCH_COMPLETE.txt"
