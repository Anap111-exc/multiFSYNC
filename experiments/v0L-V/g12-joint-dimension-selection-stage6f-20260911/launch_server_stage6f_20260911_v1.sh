#!/usr/bin/env bash
set -euo pipefail

: "${STAGE6F_OUTPUT_ROOT:?STAGE6F_OUTPUT_ROOT is required}"
: "${STAGE6F_STAGE6B_ROOT:?STAGE6F_STAGE6B_ROOT is required}"
: "${STAGE6F_OBSERVATION_BUNDLE:?STAGE6F_OBSERVATION_BUNDLE is required}"
: "${STAGE6F_TRUTH_BUNDLE:?STAGE6F_TRUTH_BUNDLE is required}"
: "${STAGE6F_R_LIB:?STAGE6F_R_LIB is required}"

workers="${STAGE6F_WORKERS:-7}"
case "$workers" in
  1|2|3|4|5|6|7) ;;
  *) echo "STAGE6F_WORKERS must be 1 through 7" >&2; exit 2 ;;
esac

script_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner="$script_root/run_stage6f_joint_dimension_selection_20260911_v1.R"
evaluator="$script_root/evaluate_stage6f_joint_dimension_selection_20260911_v1.R"
common_args=(
  "--output-root=$STAGE6F_OUTPUT_ROOT"
  "--stage6b-root=$STAGE6F_STAGE6B_ROOT"
  "--observation-bundle=$STAGE6F_OBSERVATION_BUNDLE"
)

mkdir -p "$STAGE6F_OUTPUT_ROOT/logs"
Rscript "$runner" --action=check "${common_args[@]}" \
  >"$STAGE6F_OUTPUT_ROOT/logs/check.log" 2>&1
Rscript "$runner" --action=prepare "${common_args[@]}" \
  >"$STAGE6F_OUTPUT_ROOT/logs/prepare.log" 2>&1

fit_ids="$STAGE6F_OUTPUT_ROOT/NEW_FIT_IDS.txt"
STAGE6F_FIT_IDS="$fit_ids" Rscript -e 'x <- read.csv(file.path(Sys.getenv("STAGE6F_OUTPUT_ROOT"), "FIT_MANIFEST.csv"), stringsAsFactors=FALSE); stopifnot(sum(x$execute_new_fit)==20L); writeLines(x$fit_id[x$execute_new_fit], Sys.getenv("STAGE6F_FIT_IDS"))'

export script_root runner
export STAGE6F_OUTPUT_ROOT STAGE6F_STAGE6B_ROOT STAGE6F_OBSERVATION_BUNDLE
run_one_fit() {
  fit_id="$1"
  safe_id="${fit_id//[^A-Za-z0-9_.-]/_}"
  Rscript "$runner" --action=run-fit \
    "--output-root=$STAGE6F_OUTPUT_ROOT" \
    "--stage6b-root=$STAGE6F_STAGE6B_ROOT" \
    "--observation-bundle=$STAGE6F_OBSERVATION_BUNDLE" \
    "--fit-id=$fit_id" \
    >"$STAGE6F_OUTPUT_ROOT/logs/${safe_id}.log" 2>&1
}
export -f run_one_fit
xargs -r -n 1 -P "$workers" bash -c 'run_one_fit "$1"' _ <"$fit_ids"

Rscript "$runner" --action=freeze-starts "${common_args[@]}" \
  >"$STAGE6F_OUTPUT_ROOT/logs/freeze-starts.log" 2>&1
Rscript "$runner" --action=score-winners "${common_args[@]}" \
  >"$STAGE6F_OUTPUT_ROOT/logs/score-winners.log" 2>&1
Rscript "$runner" --action=freeze-dimensions "${common_args[@]}" \
  >"$STAGE6F_OUTPUT_ROOT/logs/freeze-dimensions.log" 2>&1
Rscript "$evaluator" \
  "--output-root=$STAGE6F_OUTPUT_ROOT" \
  "--truth-path=$STAGE6F_TRUTH_BUNDLE" \
  "--authorization-token=USER_AUTHORIZED_STAGE6F_20260911" \
  >"$STAGE6F_OUTPUT_ROOT/logs/evaluate.log" 2>&1

printf '%s\n' \
  'status=STAGE6F_COMPLETE' \
  'registered_endpoints=32' \
  'reused_endpoints=12' \
  'new_fits=20' \
  'continuations=0' \
  'formal_paper_mc_result=FALSE' \
  >"$STAGE6F_OUTPUT_ROOT/STAGE6F_COMPLETE.txt"
echo "STAGE6F_COMPLETE output_root=$STAGE6F_OUTPUT_ROOT"

if [ "${STAGE6F_AUTO_STOP:-false}" = "true" ]; then
  : "${STAGE6F_INSTANCE_ID:?STAGE6F_INSTANCE_ID is required for auto-stop}"
  : "${STAGE6F_REGION_ID:?STAGE6F_REGION_ID is required for auto-stop}"
  : "${STAGE6F_ALIYUN_PROFILE:?STAGE6F_ALIYUN_PROFILE is required for auto-stop}"
  command -v aliyun >/dev/null 2>&1
  printf '%s\n' \
    'status=AUTO_STOP_REQUESTED' \
    "instance_id=$STAGE6F_INSTANCE_ID" \
    'stopped_mode=StopCharging' \
    >"$STAGE6F_OUTPUT_ROOT/AUTO_STOP_REQUESTED.txt"
  sync
  aliyun ecs StopInstance \
    --profile "$STAGE6F_ALIYUN_PROFILE" \
    --RegionId "$STAGE6F_REGION_ID" \
    --InstanceId "$STAGE6F_INSTANCE_ID" \
    --StoppedMode StopCharging \
    --ForceStop false \
    --DryRun false \
    >"$STAGE6F_OUTPUT_ROOT/AUTO_STOP_API_RESPONSE.txt" 2>&1
  printf '%s\n' 'status=STOP_API_ACCEPTED' \
    >"$STAGE6F_OUTPUT_ROOT/AUTO_STOP_API_ACCEPTED.txt"
  sync
fi
