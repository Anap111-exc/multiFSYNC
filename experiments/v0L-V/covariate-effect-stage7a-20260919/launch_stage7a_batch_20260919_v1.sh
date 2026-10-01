#!/usr/bin/env bash
set -euo pipefail

source_root="${1:?Pass the independent Stage 7A source directory}"
output_root="${2:?Pass the existing generated and benchmarked result root}"
model_library="${3:?Pass the independent Stage 7A installed R library}"
dependency_library="${4:?Pass the existing Splinets dependency library}"
workers="${STAGE7A_OUTER_WORKERS:-8}"
runner="$source_root/experiments/v0L-V/covariate-effect-stage7a-20260919/run_covariate_stage7a_20260919_v1.R"

test -f "$runner"
test -d "$model_library"
test -d "$dependency_library"
test -f "$output_root/preflight/BENCHMARK_COMPLETE.txt"
test ! -e "$output_root/BATCH_PROCESS_STATUS.csv"
export R_LIBS="$model_library:$dependency_library${R_LIBS:+:$R_LIBS}"
export STAGE7A_OUTER_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
set +e
/usr/bin/time -v Rscript "$runner" --action=batch \
  --output-root="$output_root"
run_status=$?
set -e

printf 'exit_status=%s\ncompleted_at=%s\ntruth_unseal_authorized=FALSE\nscientific_evaluation_started=FALSE\n' \
  "$run_status" "$(date --iso-8601=seconds)" \
  > "$output_root/SUPERVISOR_FINISHED.txt"
sync

# Reuse the existing instance-role self-stop mechanism after the supervisor
# has finished. Do not stop the instance if another R workload is present.
if [ "${STAGE7A_AUTO_STOP:-TRUE}" = TRUE ] && \
   command -v aliyun >/dev/null 2>&1 && \
   ! pgrep -f '[/]usr/lib64/R/bin/exec/R|[R]script' >/dev/null; then
  aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun false
fi
exit "$run_status"
