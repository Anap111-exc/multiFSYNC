#!/usr/bin/env bash
# One-shot cost-saving stop after the authorized Stage 7I read-only analysis.
set -euo pipefail

parent=/root/multifsync-results
freeze="$parent/stage7i_blind_freeze_20260929_v1/BLIND_ENDPOINTS_FROZEN.txt"
science="$parent/stage7i_scientific_evaluation_20260929_v5/EVALUATION_COMPLETE.txt"
trace="$parent/stage7i_readonly_transition_20260930_v1/TRACE_REDUCTION_COMPLETE.txt"
snapshot="$parent/stage7i_readonly_snapshot_transition_20260930_v2/SNAPSHOT_RECONSTRUCTION_COMPLETE.txt"
residual="$parent/stage7i_readonly_fixed_direction_20260930_v2/FIXED_DIRECTION_RESIDUAL_COMPLETE.txt"
record="$parent/stage7i_readonly_autostop_20260930_v1"

test -f "$freeze" && grep -qx 'STAGE7I_BLIND_ENDPOINTS_FROZEN' <(head -n 1 "$freeze")
test -f "$science" && grep -qx 'status=STAGE7I_ENDPOINT_SCIENCE_COMPLETE' "$science"
for marker in "$trace" "$snapshot" "$residual"; do
  test -f "$marker"
  grep -qx 'status=COMPLETE' "$marker"
  grep -qx 'truth_opened=FALSE' "$marker"
  grep -qx 'new_fit_started=FALSE' "$marker"
done
test ! -e "$record"
if pgrep -x Rscript >/dev/null || pgrep -x R >/dev/null; then
  echo 'A live R or Rscript process prevents automatic stop.' >&2
  exit 1
fi
mkdir -m 0700 "$record"
sha256sum "$freeze" "$science" "$trace" "$snapshot" "$residual" \
  > "$record/GATE_HASHES.sha256"
date --iso-8601=seconds > "$record/GATES_PASS_AT.txt"
sync

if aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun true \
    > "$record/STOP_DRY_RUN_RESPONSE.txt" 2>&1; then
  dry_status=0
else
  dry_status=$?
fi
printf 'dry_run_exit_status=%s\n' "$dry_status" > "$record/STOP_DRY_RUN_STATUS.txt"
grep -q 'DryRunOperation' "$record/STOP_DRY_RUN_RESPONSE.txt"
printf 'requested_at=%s\nmode=StopCharging\nforce_stop=false\n' \
  "$(date --iso-8601=seconds)" > "$record/STOP_REQUEST_INTENDED.txt"
sync

if aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun false \
    > "$record/STOP_API_RESPONSE.txt" 2>&1; then
  api_status=0
else
  api_status=$?
fi
printf 'api_exit_status=%s\nresponse_recorded_at=%s\n' \
  "$api_status" "$(date --iso-8601=seconds)" > "$record/STOP_API_STATUS.txt"
sync
test "$api_status" -eq 0
