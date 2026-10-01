#!/usr/bin/env bash
set -euo pipefail

repo_root="${1:?Pass the independent Stage 6L source directory}"
r_library="${2:?Pass the independent installed R library}"
output_root="${3:?Pass a new output directory outside the source tree}"
dependency_library="${4:-}"
workers="${STAGE6L_OUTER_WORKERS:-8}"
runner="$repo_root/experiments/v0L-V/g12-loading-dgp-robustness-stage6l-20260915/run_stage6l_loading_bridge_20260915_v1.R"

test -f "$runner"
test -d "$r_library"
test ! -e "$output_root"
if [ -n "$dependency_library" ]; then test -d "$dependency_library"; fi
export STAGE6A_R_LIB="$r_library" STAGE6L_OUTER_WORKERS="$workers"
export R_LIBS="$r_library${dependency_library:+:$dependency_library}${R_LIBS:+:$R_LIBS}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
set +e
/usr/bin/time -v Rscript "$runner" --action=full --output-root="$output_root"
run_status=$?
set -e
if [ -d "$output_root" ]; then
  printf 'exit_status=%s\ncompleted_at=%s\ntruth_unsealed=FALSE\nscientific_evaluation_started=FALSE\n' \
    "$run_status" "$(date --iso-8601=seconds)" \
    > "$output_root/SUPERVISOR_FINISHED.txt"
fi
sync

# Reuse the existing instance-role self-stop permission, not a new watcher.
# Never stop while another R workload is present.
if [ "${STAGE6L_AUTO_STOP:-TRUE}" = TRUE ] && \
   command -v aliyun >/dev/null 2>&1 && \
   ! pgrep -f '[/]usr/lib64/R/bin/exec/R|[R]script' >/dev/null; then
  aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun false
fi
exit "$run_status"
