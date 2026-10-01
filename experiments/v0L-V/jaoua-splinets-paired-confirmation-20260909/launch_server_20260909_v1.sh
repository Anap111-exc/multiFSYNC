#!/usr/bin/env bash
set -euo pipefail

repo="/root/multiFSYNC-jaoua-paired-20260909"
library="/root/multiFSYNC-jaoua-paired-lib-20260909"
output_root="/root/jspc-output-20260909-v1"
experiment="$repo/experiments/v0L-V/jaoua-splinets-paired-confirmation-20260909"
runner="$experiment/run_paired_confirmation_20260909_v1.R"
watcher="$experiment/server_autostop_20260909_v1.sh"

if [ -e "$output_root" ]; then
  printf 'Refusing pre-existing output root: %s\n' "$output_root" >&2
  exit 2
fi
test -f "$runner"
test -f "$watcher"
test -d "$library"
command -v Rscript >/dev/null
command -v aliyun >/dev/null

mkdir -m 700 "$output_root"
export R_LIBS="$library"
export STAGE6A_R_LIB="$library"
export JSPC_OUTER_WORKERS=7
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

nohup Rscript "$runner" --action=full --output-root="$output_root" \
  > "$output_root/SUPERVISOR.log" 2>&1 &
runner_pid=$!
printf '%s\n' "$runner_pid" > "$output_root/SUPERVISOR.pid"

nohup bash "$watcher" > "$output_root/AUTO_STOP_WATCHER.nohup" 2>&1 &
watcher_pid=$!
printf '%s\n' "$watcher_pid" > "$output_root/AUTO_STOP_WATCHER.pid"

sleep 2
kill -0 "$runner_pid"
kill -0 "$watcher_pid"
printf 'LAUNCH_PASS runner_pid=%s watcher_pid=%s workers=7 per_fit_n_cpus=1\n' \
  "$runner_pid" "$watcher_pid"

