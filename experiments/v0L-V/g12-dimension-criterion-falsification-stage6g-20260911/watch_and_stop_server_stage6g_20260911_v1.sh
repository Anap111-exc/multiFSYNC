#!/usr/bin/env bash
set -u

output_root="${STAGE6G_OUTPUT_ROOT:-/root/multifsync-results/stage6g_20260911_v1}"
pid_file="${STAGE6G_PID_FILE:-/root/stage6g_20260911_v1.pid}"
runner_path="${STAGE6G_RUNNER_PATH:-/root/multiFSYNC-stage6b-20260903-fa10cbc/experiments/v0L-V/g12-dimension-criterion-falsification-stage6g-20260911/run_stage6g_20260911_v1.R}"
instance_id="i-bp1egmsairexm9w1rnks"
region_id="cn-hangzhou"
profile_name="multifsync-self-stop"
log_file="/root/stage6g-autostop-20260911-v1.log"
lock_dir="/run/stage6g-autostop-20260911-v1.lock"

log_message() {
  printf '%s %s\n' "$(date --iso-8601=seconds)" "$*" >> "$log_file"
}

if ! mkdir "$lock_dir" 2>/dev/null; then
  log_message "watcher_not_started reason=lock_exists"
  exit 1
fi
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT

for _ in $(seq 1 60); do
  [ -s "$pid_file" ] && break
  sleep 1
done
if [ ! -s "$pid_file" ]; then
  log_message "watcher_aborted reason=pid_file_missing"
  exit 2
fi

supervisor_pid=$(tr -cd '0-9' < "$pid_file")
if [ -z "$supervisor_pid" ]; then
  log_message "watcher_aborted reason=pid_invalid"
  exit 3
fi
log_message "watcher_armed pid=$supervisor_pid output_root=$output_root"

while kill -0 "$supervisor_pid" 2>/dev/null; do
  sleep 60
done
while pgrep -f "$runner_path" >/dev/null 2>&1; do
  sleep 30
done

complete_count=$(find "$output_root" -type f -name FIT_COMPLETE.txt 2>/dev/null | wc -l)
error_count=$(find "$output_root" -type f -name FIT_ERROR.txt 2>/dev/null | wc -l)
terminal_count=$(find "$output_root" -type f -name terminal_record.rds 2>/dev/null | wc -l)
if [ -f "$output_root/STAGE6G_TRUTH_FREE_COMPLETE.txt" ] && \
   [ -f "$output_root/truth_free_selection/DIMENSION_SELECTION_FROZEN.txt" ] && \
   [ "$terminal_count" -eq 120 ]; then
  outcome="complete_truth_remains_sealed"
else
  outcome="incomplete_retained_for_diagnosis"
fi

mkdir -p "$output_root"
printf 'status=AUTO_STOP_REQUESTED\noutcome=%s\ncomplete=%s\nerrors=%s\nterminal_records=%s\nrequested_at=%s\ninstance_id=%s\nstopped_mode=StopCharging\n' \
  "$outcome" "$complete_count" "$error_count" "$terminal_count" \
  "$(date --iso-8601=seconds)" "$instance_id" \
  > "$output_root/AUTO_STOP_REQUESTED.txt"
log_message "stop_requested outcome=$outcome complete=$complete_count errors=$error_count terminals=$terminal_count"
sync

if aliyun ecs StopInstance \
  --profile "$profile_name" \
  --RegionId "$region_id" \
  --InstanceId "$instance_id" \
  --StoppedMode StopCharging \
  --ForceStop false \
  --DryRun false > "$output_root/AUTO_STOP_API_RESPONSE.txt" 2>&1; then
  log_message "stop_api_accepted stopped_mode=StopCharging"
  printf 'status=STOP_API_ACCEPTED\naccepted_at=%s\noutcome=%s\n' \
    "$(date --iso-8601=seconds)" "$outcome" \
    > "$output_root/AUTO_STOP_API_ACCEPTED.txt"
  sync
  exit 0
else
  api_status=$?
fi

log_message "auto_stop_failed exit_status=$api_status outcome=$outcome"
printf 'status=AUTO_STOP_FAILED\nexit_status=%s\nfailed_at=%s\noutcome=%s\n' \
  "$api_status" "$(date --iso-8601=seconds)" "$outcome" \
  > "$output_root/AUTO_STOP_FAILED.txt"
exit "$api_status"
