#!/usr/bin/env bash
set -euo pipefail

stage6d_repo="${STAGE6D_REPO:-/root/multiFSYNC-stage6d-20260906}"
stage6d_lib="${STAGE6D_R_LIB:-/root/multiFSYNC-stage6d-lib-20260906}"
stage6d_output="${STAGE6D_OUTPUT_ROOT:-/root/g12-dimension-selection-stage6d-20260906-v1}"
stage6d_dir="$stage6d_repo/experiments/v0L-V/g12-dimension-selection-stage6d-20260906"
runner="$stage6d_dir/run_stage6d_20260906_v1.R"

test -f "$stage6d_output/SERVER_PREPARATION_COMPLETE.txt"
test ! -f "$stage6d_output/STAGE6D_LAUNCHED.txt"
test ! -f "$stage6d_output/TRUTH_UNSEAL_AUTHORIZATION.txt"

export STAGE6D_REPO="$stage6d_repo"
export STAGE6D_R_LIB="$stage6d_lib"
export STAGE6D_OUTPUT_ROOT="$stage6d_output"
export STAGE6D_WORKERS="7"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

nohup bash -lc "set -euo pipefail; export STAGE6D_R_LIB='$stage6d_lib'; export STAGE6D_WORKERS=7; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1; Rscript '$runner' --action=prepare --output-root='$stage6d_output'; Rscript '$runner' --action=run --output-root='$stage6d_output'" \
  > "$stage6d_output/STAGE6D_LAUNCH.stdout.log" \
  2> "$stage6d_output/STAGE6D_LAUNCH.stderr.log" < /dev/null &
stage6d_pid=$!

{
  echo "status=STAGE6D_LAUNCHED"
  echo "launched_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pid=$stage6d_pid"
  echo "registered_data=3"
  echo "registered_configs=5"
  echo "registered_starts_per_config=12"
  echo "registered_fits=180"
  echo "outer_workers=7"
  echo "n_cpus_per_fit=1"
  echo "truth_unseal_authorized=FALSE"
  echo "continuation_started=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$stage6d_output/STAGE6D_LAUNCHED.txt"

echo "STAGE6D_LAUNCHED pid=$stage6d_pid output_root=$stage6d_output"
