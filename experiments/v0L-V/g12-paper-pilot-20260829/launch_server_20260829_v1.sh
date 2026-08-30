#!/usr/bin/env bash
set -euo pipefail

pilot_repo="${PAPER_PILOT_REPO:-/root/multiFSYNC-g12-paper-pilot-20260829}"
pilot_lib="${PAPER_PILOT_R_LIB:-/root/multiFSYNC-g12-paper-pilot-lib-20260829}"
pilot_output="${PAPER_PILOT_OUTPUT_ROOT:-/root/g12-paper-pilot-20260829-v1}"
pilot_dir="$pilot_repo/experiments/v0L-V/g12-paper-pilot-20260829"
runner="$pilot_dir/paper_pilot_runner_20260829_v1.R"

test -f "$pilot_output/SERVER_PREPARATION_COMPLETE.txt"
test ! -f "$pilot_output/PILOT_LAUNCHED.txt"
test ! -f "$pilot_output/TRUTH_UNSEAL_AUTHORIZATION.txt"

export PAPER_PILOT_REPO="$pilot_repo"
export PAPER_PILOT_R_LIB="$pilot_lib"
export PAPER_PILOT_OUTPUT_ROOT="$pilot_output"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

nohup bash -lc "set -euo pipefail; export PAPER_PILOT_R_LIB='$pilot_lib'; Rscript '$runner' --action=generate --output-root='$pilot_output'; Rscript '$runner' --action=run --output-root='$pilot_output'; Rscript '$runner' --action=summarize --output-root='$pilot_output'" \
  > "$pilot_output/PILOT_LAUNCH.stdout.log" \
  2> "$pilot_output/PILOT_LAUNCH.stderr.log" < /dev/null &
pilot_pid=$!

{
  echo "status=PILOT_LAUNCHED"
  echo "launched_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "pid=$pilot_pid"
  echo "registered_data=2"
  echo "registered_fits=26"
  echo "workers=3"
  echo "n_cpus_per_fit=1"
  echo "truth_unseal_authorized=FALSE"
  echo "evaluation_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$pilot_output/PILOT_LAUNCHED.txt"

echo "PILOT_LAUNCHED pid=$pilot_pid output_root=$pilot_output"
