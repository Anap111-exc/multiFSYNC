#!/usr/bin/env bash
set -euo pipefail

stage6d_repo="${STAGE6D_REPO:-/root/multiFSYNC-stage6d-20260906}"
stage6d_lib="${STAGE6D_R_LIB:-/root/multiFSYNC-stage6d-lib-20260906}"
stage6d_output="${STAGE6D_OUTPUT_ROOT:-/root/g12-dimension-selection-stage6d-20260906-v1}"
workers="${STAGE6D_EVAL_WORKERS:-7}"
script="$stage6d_repo/experiments/v0L-V/g12-dimension-selection-stage6d-20260906/evaluate_stage6d_20260908_v1.R"

test -f "$stage6d_output/STAGE6D_COMPLETE.txt"
test -f "$stage6d_output/TRUTH_UNSEAL_AUTHORIZATION.txt"
test ! -e "$stage6d_output/evaluation"
test "$workers" -ge 1
test "$workers" -le 7

export STAGE6D_R_LIB="$stage6d_lib"
export STAGE6D_EVAL_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

nohup Rscript "$script" --output-root="$stage6d_output" --workers="$workers" \
  > "$stage6d_output/STAGE6D_EVALUATION.stdout.log" \
  2> "$stage6d_output/STAGE6D_EVALUATION.stderr.log" < /dev/null &
evaluation_pid=$!

printf 'status=EVALUATION_LAUNCHED\nlaunched_utc=%s\npid=%s\nworkers=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$evaluation_pid" "$workers" \
  > "$stage6d_output/STAGE6D_EVALUATION_LAUNCHED.txt"
echo "STAGE6D_EVALUATION_LAUNCHED pid=$evaluation_pid workers=$workers"
