#!/usr/bin/env bash
set -euo pipefail

stage6a_repo="${STAGE6A_REPO:-/root/multiFSYNC-stage6a-20260831}"
stage6a_lib="${STAGE6A_R_LIB:-/root/multiFSYNC-stage6a-lib-20260831}"
stage6a_output="${STAGE6A_OUTPUT_ROOT:-/root/g12-dimension-overspec-stage6a-20260831-v1}"
stage6a_dir="$stage6a_repo/experiments/v0L-V/g12-dimension-overspec-stage6a-20260831"
runner="$stage6a_dir/stage6a_runner_20260831_v1.R"

test -f "$stage6a_output/SERVER_PREPARATION_COMPLETE.txt"
test ! -f "$stage6a_output/STAGE6A_LAUNCHED.txt"
test ! -f "$stage6a_output/TRUTH_UNSEAL_AUTHORIZATION.txt"

export STAGE6A_REPO="$stage6a_repo"
export STAGE6A_R_LIB="$stage6a_lib"
export STAGE6A_OUTPUT_ROOT="$stage6a_output"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

nohup bash -lc "set -euo pipefail; export STAGE6A_R_LIB='$stage6a_lib'; Rscript '$runner' --action=generate --output-root='$stage6a_output'; Rscript '$runner' --action=capacity --output-root='$stage6a_output'; Rscript '$runner' --action=run --output-root='$stage6a_output'; Rscript '$runner' --action=summarize --output-root='$stage6a_output'" \
  > "$stage6a_output/STAGE6A_LAUNCH.stdout.log" \
  2> "$stage6a_output/STAGE6A_LAUNCH.stderr.log" < /dev/null &
stage6a_pid=$!

{
  echo "status=STAGE6A_LAUNCHED"
  echo "launched_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pid=$stage6a_pid"
  echo "registered_data=1"
  echo "registered_fit_configs=3"
  echo "registered_fits=36"
  echo "registered_capacity_probes=3"
  echo "workers=3"
  echo "n_cpus_per_fit=1"
  echo "truth_unseal_authorized=FALSE"
  echo "evaluation_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$stage6a_output/STAGE6A_LAUNCHED.txt"

echo "STAGE6A_LAUNCHED pid=$stage6a_pid output_root=$stage6a_output"
