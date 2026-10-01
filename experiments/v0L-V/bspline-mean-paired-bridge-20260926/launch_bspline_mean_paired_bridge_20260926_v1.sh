#!/usr/bin/env bash
set -euo pipefail

repo_root=/root/multiFSYNC-mean-bridge-20260926-v1
runner="$repo_root/experiments/v0L-V/bspline-mean-paired-bridge-20260926/run_bspline_mean_paired_bridge_20260926_v1.R"
output_root=/root/multifsync-results/bspline_mean_paired_bridge_20260926_v1
workers=8

export R_LIBS_USER=/root/multiFSYNC-mean-bridge-lib-20260926:/root/rb22lib:/root/multiFSYNC-jaoua-paired-lib-20260909
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

finish_and_schedule_stop() {
  status=$?
  trap - EXIT
  printf 'supervisor_exit_status=%s\nfinished_utc=%s\ntruth_unseal_authorized=FALSE\nscientific_evaluation_started=FALSE\n' \
    "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    >"$output_root/SUPERVISOR_EXIT_STATUS.txt" || true
  sync
  nohup bash -c 'sleep 60; /usr/local/bin/aliyun ecs StopInstance --profile multifsync-self-stop --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks --StoppedMode StopCharging --ForceStop false --DryRun false' \
    >"$output_root/AUTO_STOP.log" 2>&1 </dev/null &
  exit "$status"
}
trap finish_and_schedule_stop EXIT

test -f "$runner"
test -d /root/multiFSYNC-mean-bridge-lib-20260926

if [[ ! -f "$output_root/DATA_GENERATION_COMPLETE.txt" ]]; then
  Rscript "$runner" --action=generate --output-root="$output_root"
fi

if [[ ! -f "$output_root/preflight/BENCHMARK_COMPLETE.txt" ]]; then
  /usr/bin/time -v Rscript "$runner" --action=benchmark \
    --output-root="$output_root"
fi

/usr/bin/time -v Rscript "$runner" --action=batch \
  --output-root="$output_root" --workers="$workers"
