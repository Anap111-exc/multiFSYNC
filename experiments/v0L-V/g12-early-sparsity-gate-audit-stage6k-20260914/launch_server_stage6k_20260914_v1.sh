#!/usr/bin/env bash
set -euo pipefail

repo_root="${STAGE6K_REPO_ROOT:-/root/multiFSYNC-stage6b-20260903-fa10cbc}"
package_root="${STAGE6K_PACKAGE_SOURCE_ROOT:-/root/multiFSYNC-jaoua-paired-20260909}"
r_library="${STAGE6K_R_LIB:-/root/multiFSYNC-jaoua-paired-lib-20260909}"
stage6g_output="${STAGE6K_STAGE6G_OUTPUT_ROOT:-/root/multifsync-results/stage6g_20260911_v1}"
stage6j_output="${STAGE6K_STAGE6J_OUTPUT_ROOT:-/root/multifsync-results/stage6j_20260914_v1}"
output_root="${STAGE6K_OUTPUT_ROOT:-/root/multifsync-results/stage6k_20260914_v1}"
workers="${STAGE6K_WORKERS:-4}"
experiment_root="$repo_root/experiments/v0L-V/g12-early-sparsity-gate-audit-stage6k-20260914"
runner="$experiment_root/run_stage6k_20260914_v1.R"
evaluator="$experiment_root/evaluate_stage6k_20260914_v1.R"

test -f "$runner"
test -f "$evaluator"
test -d "$package_root/R"
test -d "$r_library"
test -f "$stage6g_output/STAGE6G_TRUTH_FREE_COMPLETE.txt"
test -f "$stage6j_output/STAGE6J_TRUTH_FREE_COMPLETE.txt"
test ! -e "$output_root"

export STAGE6K_PACKAGE_SOURCE_ROOT="$package_root"
export STAGE6K_R_LIB="$r_library"
export STAGE6K_STAGE6G_OUTPUT_ROOT="$stage6g_output"
export STAGE6K_STAGE6J_OUTPUT_ROOT="$stage6j_output"
export STAGE6K_OUTPUT_ROOT="$output_root"
export STAGE6K_WORKERS="$workers"
export STAGE6K_EVAL_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=preflight
Rscript "$runner" --action=prepare --output-root="$output_root"
Rscript "$runner" --action=smoke --output-root="$output_root"
Rscript "$runner" --action=run --output-root="$output_root"
Rscript "$evaluator" --output-root="$output_root" --workers="$workers"

printf 'status=STAGE6K_LAUNCH_COMPLETE\ncompleted_at=%s\ntruth_evaluated_after_freeze=TRUE\nnew_data=FALSE\ncontinuation=FALSE\n' \
  "$(date --iso-8601=seconds)" > "$output_root/STAGE6K_LAUNCH_COMPLETE.txt"
