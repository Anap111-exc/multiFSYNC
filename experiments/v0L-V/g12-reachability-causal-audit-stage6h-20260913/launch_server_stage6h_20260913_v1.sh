#!/usr/bin/env bash
set -euo pipefail

repo_root="${STAGE6H_REPO_ROOT:-/root/multiFSYNC-stage6b-20260903-fa10cbc}"
package_root="${STAGE6H_PACKAGE_SOURCE_ROOT:-/root/multiFSYNC-jaoua-paired-20260909}"
r_library="${STAGE6H_R_LIB:-/root/multiFSYNC-jaoua-paired-lib-20260909}"
stage6g_output="${STAGE6H_STAGE6G_OUTPUT_ROOT:-/root/multifsync-results/stage6g_20260911_v1}"
stage6g_evaluation="${STAGE6H_STAGE6G_EVALUATION_ROOT:-$stage6g_output/truth_evaluation_20260912_v1}"
output_root="${STAGE6H_OUTPUT_ROOT:-/root/multifsync-results/stage6h_20260913_v1}"
workers="${STAGE6H_WORKERS:-7}"
experiment_root="$repo_root/experiments/v0L-V/g12-reachability-causal-audit-stage6h-20260913"
runner="$experiment_root/run_stage6h_20260913_v1.R"
evaluator="$experiment_root/evaluate_stage6h_20260913_v1.R"

test -f "$runner"
test -f "$evaluator"
test -d "$package_root/R"
test -d "$r_library"
test -f "$stage6g_output/STAGE6G_TRUTH_FREE_COMPLETE.txt"
test -f "$stage6g_evaluation/EVALUATION_COMPLETE.txt"
test ! -e "$output_root"
mkdir -p "$(dirname "$output_root")"

export STAGE6H_PACKAGE_SOURCE_ROOT="$package_root"
export STAGE6H_R_LIB="$r_library"
export STAGE6H_STAGE6G_OUTPUT_ROOT="$stage6g_output"
export STAGE6H_STAGE6G_EVALUATION_ROOT="$stage6g_evaluation"
export STAGE6H_OUTPUT_ROOT="$output_root"
export STAGE6H_WORKERS="$workers"
export STAGE6H_EVAL_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=prepare --output-root="$output_root"
Rscript "$runner" --action=run --output-root="$output_root"
Rscript "$evaluator" --output-root="$output_root" --workers="$workers"

printf 'status=STAGE6H_WORKFLOW_COMPLETE\ncompleted_at=%s\n' \
  "$(date --iso-8601=seconds)" > "$output_root/STAGE6H_WORKFLOW_COMPLETE.txt"

