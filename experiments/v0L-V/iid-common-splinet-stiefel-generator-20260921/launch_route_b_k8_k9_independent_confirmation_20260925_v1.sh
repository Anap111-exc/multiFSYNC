#!/usr/bin/env bash
set -euo pipefail

repo_root=/root/multiFSYNC-routeb-bridge-20260922-v1
runner="$repo_root/experiments/v0L-V/iid-common-splinet-stiefel-generator-20260921/run_route_b_k8_k9_independent_confirmation_20260925_v1.R"
output_root=/root/multifsync-results/route_b_k8_k9_independent_confirmation_20260925_v1
workers=8

export R_LIBS_USER=/root/rb22lib:/root/multiFSYNC-jaoua-paired-lib-20260909
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

finish_and_schedule_stop() {
  status=$?
  trap - EXIT
  printf 'supervisor_exit_status=%s\nfinished_utc=%s\n' \
    "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    >"$output_root/SUPERVISOR_EXIT_STATUS.txt" || true
  sync
  nohup bash -c 'sleep 60; /usr/local/bin/aliyun ecs StopInstance --profile multifsync-self-stop --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks --StoppedMode StopCharging --ForceStop false --DryRun false' \
    >"$output_root/AUTO_STOP.log" 2>&1 </dev/null &
  exit "$status"
}
trap finish_and_schedule_stop EXIT

if [[ ! -f "$output_root/DATA_GENERATION_COMPLETE.txt" ]]; then
  Rscript "$runner" --action=generate --output-root="$output_root"
fi

if [[ ! -f "$output_root/preflight/BENCHMARK_COMPLETE.txt" ]]; then
  /usr/bin/time -v Rscript "$runner" --action=benchmark \
    --output-root="$output_root"
fi

/usr/bin/time -v Rscript "$runner" --action=batch \
  --output-root="$output_root" --workers="$workers"
