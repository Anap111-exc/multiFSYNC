#!/usr/bin/env bash
# Existing ECS RAM role profile; request StopCharging, never OS-only poweroff.
set -euo pipefail
test "$#" -eq 1
finished=$(realpath "$1")
record=/root/multifsync-results/stage7l_autostop_20261001_v1
test -f "$finished"
test ! -e "$record"
if pgrep -x Rscript >/dev/null || pgrep -x R >/dev/null; then
  echo 'Another R job is active; skipping automatic stop.' >&2
  exit 1
fi
mkdir -m 0700 "$record"
cp "$finished" "$record/SUPERVISOR_FINISHED.txt"
sha256sum "$finished" > "$record/FINISHED_HASH.sha256"
printf 'mode=StopCharging\nforce_stop=false\nrequested_at=%s\n' \
  "$(date --iso-8601=seconds)" > "$record/STOP_REQUEST_INTENDED.txt"
sync
if aliyun ecs StopInstance --profile multifsync-self-stop \
  --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
  --StoppedMode StopCharging --ForceStop false --DryRun true \
  > "$record/STOP_DRY_RUN_RESPONSE.txt" 2>&1; then dry_status=0; else dry_status=$?; fi
printf 'dry_run_exit_status=%s\n' "$dry_status" > "$record/STOP_DRY_RUN_STATUS.txt"
grep -q 'DryRunOperation' "$record/STOP_DRY_RUN_RESPONSE.txt"
sync
if aliyun ecs StopInstance --profile multifsync-self-stop \
  --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
  --StoppedMode StopCharging --ForceStop false --DryRun false \
  > "$record/STOP_API_RESPONSE.txt" 2>&1; then api_status=0; else api_status=$?; fi
printf 'api_exit_status=%s\nresponse_recorded_at=%s\n' "$api_status" \
  "$(date --iso-8601=seconds)" > "$record/STOP_API_STATUS.txt"
sync
test "$api_status" -eq 0
