#!/usr/bin/env bash
set -euo pipefail

source_root="${1:?Pass the independent Route-B source directory}"
output_root="${2:?Pass a new result directory}"
model_library="${3:?Pass the independent installed R library}"
dependency_library="${4:?Pass the Splinets dependency library}"
workers="${ROUTEB_OUTER_WORKERS:-8}"
runner="$source_root/experiments/v0L-V/iid-common-splinet-stiefel-generator-20260921/run_route_b_paired_fit_bridge_20260922_v1.R"

test -f "$runner"
test -d "$model_library"
test -d "$dependency_library"
test ! -e "$output_root"
export R_LIBS="$model_library:$dependency_library${R_LIBS:+:$R_LIBS}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=generate --output-root="$output_root"
/usr/bin/time -v -o "$output_root/K9_BENCHMARK_RESOURCE_USAGE.txt" \
  Rscript "$runner" --action=benchmark --output-root="$output_root"

available_kib="$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)"
if [ -z "$available_kib" ] || [ "$available_kib" -lt 8388608 ]; then
  printf 'status=PRECHECK_BLOCKED\nMemAvailable_kib=%s\n' "${available_kib:-NA}" \
    > "$output_root/PRECHECK_BLOCKED.txt"
  exit 20
fi

set +e
/usr/bin/time -v -o "$output_root/BATCH_RESOURCE_USAGE.txt" \
  Rscript "$runner" --action=batch --output-root="$output_root" \
    --workers="$workers"
run_status=$?
set -e

printf 'exit_status=%s\ncompleted_at=%s\ntruth_unseal_authorized=FALSE\nscientific_evaluation_started=FALSE\n' \
  "$run_status" "$(date --iso-8601=seconds)" \
  > "$output_root/SUPERVISOR_FINISHED.txt"
sync

# Reuse the instance-role self-stop mechanism only after the supervisor ends.
if [ "${ROUTEB_AUTO_STOP:-TRUE}" = TRUE ] && \
   command -v aliyun >/dev/null 2>&1 && \
   ! pgrep -f '[/]usr/lib64/R/bin/exec/R|[R]script' >/dev/null; then
  aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun false
fi
exit "$run_status"
