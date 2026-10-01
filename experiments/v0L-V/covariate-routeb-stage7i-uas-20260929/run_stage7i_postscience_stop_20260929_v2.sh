#!/usr/bin/env bash
# One-shot recovery for the original watcher, which lost its transient fit unit.
# Run only after Stage 7I science and trajectory artifacts have been inspected.
set -euo pipefail

parent=/root/multifsync-results
freeze="$parent/stage7i_blind_freeze_20260929_v1/BLIND_ENDPOINTS_FROZEN.txt"
science="$parent/stage7i_scientific_evaluation_20260929_v5/EVALUATION_COMPLETE.txt"
trajectory="$parent/stage7i_trajectory_readout_20260929_v1/TRAJECTORY_COMPLETE.txt"
fit="$parent/stage7i_as_20260929_v1/AS_FITS_COMPLETE.txt"
u_parity="$parent/stage7i_u_20260929_v1/u_parity_v1/PARITY_PASS.txt"
record="$parent/stage7i_postscience_stop_20260929_v2"

test -f "$freeze" && grep -qx 'STAGE7I_BLIND_ENDPOINTS_FROZEN' <(head -n 1 "$freeze")
test -f "$science" && grep -qx 'status=STAGE7I_ENDPOINT_SCIENCE_COMPLETE' "$science"
test -f "$trajectory" && grep -qx 'STAGE7I_BLIND_TRAJECTORY_READOUT_COMPLETE' <(head -n 1 "$trajectory")
test -f "$fit" && grep -qx 'STAGE7I_AS_FITS_COMPLETE' <(head -n 1 "$fit")
test -f "$u_parity"
test ! -e "$record"
mkdir -m 0700 "$record"
sha256sum "$freeze" "$science" "$trajectory" "$fit" "$u_parity" > "$record/GATE_HASHES.sha256"
printf 'status=ALL_POSTSCIENCE_GATES_PASS\nchecked_at=%s\n' \
  "$(date --iso-8601=seconds)" > "$record/GATES_PASS.txt"
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
printf 'status=STOP_REQUEST_INTENDED\nrequested_at=%s\nmode=StopCharging\nforce_stop=false\n' \
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
