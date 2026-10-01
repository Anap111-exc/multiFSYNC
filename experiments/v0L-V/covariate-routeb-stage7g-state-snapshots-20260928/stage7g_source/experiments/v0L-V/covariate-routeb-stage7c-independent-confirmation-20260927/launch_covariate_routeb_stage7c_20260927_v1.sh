#!/usr/bin/env bash
set -euo pipefail

repo_root=/root/multiFSYNC-stage7c-20260927-v1
runner="$repo_root/experiments/v0L-V/covariate-routeb-stage7c-independent-confirmation-20260927/run_covariate_routeb_stage7c_20260927_v1.R"
output_root=/root/multifsync-results/covariate_routeb_stage7c_20260927_v1

export R_LIBS_USER=/root/multiFSYNC-stage7c-lib-20260927-v1:/root/rb22lib:/root/multiFSYNC-jaoua-paired-lib-20260909:/root/multiFSYNC-mean-bridge-lib-20260926
export COVARIATE_ROUTEB_STAGE7C_OUTER_WORKERS=8
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
test -f "$output_root/DATA_GENERATION_COMPLETE.txt"
test -f "$output_root/preflight/BENCHMARK_COMPLETE.txt"
test -d /root/multiFSYNC-stage7c-lib-20260927-v1

/usr/bin/time -v Rscript "$runner" --action=batch \
  --output-root="$output_root"
