#!/usr/bin/env bash
set -euo pipefail

stage6b_repo="${STAGE6B_REPO:-/root/multiFSYNC-stage6b-20260902}"
stage6b_lib="${STAGE6B_R_LIB:-/root/multiFSYNC-stage6a-lib-20260831}"
stage6b_output="${STAGE6B_OUTPUT_ROOT:-/root/g12-truth-free-dimension-selection-stage6b-20260902-v1}"
stage6b_dir="$stage6b_repo/experiments/v0L-V/g12-truth-free-dimension-selection-stage6b-20260901"
runner="$stage6b_dir/run_stage6b_predictive_pilot_20260902_v1.R"

test -f "$stage6b_output/SERVER_PREPARATION_COMPLETE.txt"
test -f "$stage6b_output/STAGE6B_PREPARED.txt"
test ! -f "$stage6b_output/STAGE6B_LAUNCHED.txt"
test ! -f "$stage6b_output/TRUTH_UNSEAL_AUTHORIZATION.txt"
test ! -f "$stage6b_output/STAGE6B_PILOT_COMPLETE.txt"

export STAGE6B_REPO="$stage6b_repo"
export STAGE6B_R_LIB="$stage6b_lib"
export STAGE6B_OUTPUT_ROOT="$stage6b_output"
export STAGE6B_WORKERS="${STAGE6B_WORKERS:-3}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

nohup bash -lc "set -euo pipefail; export STAGE6B_R_LIB='$stage6b_lib'; export STAGE6B_OUTPUT_ROOT='$stage6b_output'; export STAGE6B_WORKERS='$STAGE6B_WORKERS'; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1; Rscript '$runner' --action=run --output-root='$stage6b_output'" \
  > "$stage6b_output/STAGE6B_LAUNCH.stdout.log" \
  2> "$stage6b_output/STAGE6B_LAUNCH.stderr.log" < /dev/null &
stage6b_pid=$!

{
  echo "status=STAGE6B_LAUNCHED"
  echo "launched_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pid=$stage6b_pid"
  echo "source_data_reused=1"
  echo "new_data_generated=FALSE"
  echo "factor_phase_fits=12"
  echo "conditional_fpca_phase_fits=8"
  echo "maximum_fits=20"
  echo "workers=$STAGE6B_WORKERS"
  echo "n_cpus_per_fit=1"
  echo "truth_access_authorized=FALSE"
  echo "continuation_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$stage6b_output/STAGE6B_LAUNCHED.txt"

echo "STAGE6B_LAUNCHED pid=$stage6b_pid output_root=$stage6b_output"
