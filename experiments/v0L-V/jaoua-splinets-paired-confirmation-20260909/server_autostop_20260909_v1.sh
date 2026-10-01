#!/usr/bin/env bash
set -u

output_root="/root/jspc-output-20260909-v1"
runner_path="/root/multiFSYNC-jaoua-paired-20260909/experiments/v0L-V/jaoua-splinets-paired-confirmation-20260909/run_paired_confirmation_20260909_v1.R"
instance_id="i-bp1egmsairexm9w1rnks"
region_id="cn-hangzhou"
profile_name="multifsync-self-stop"
log_file="$output_root/AUTO_STOP_WATCHER.log"
lock_dir="/run/jspc-autostop-20260909.lock"

log_message() {
  printf '%s %s\n' "$(date --iso-8601=seconds)" "$*" >> "$log_file"
}

if ! mkdir "$lock_dir" 2>/dev/null; then
  log_message "watcher_not_started reason=lock_exists"
  exit 1
fi
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT
log_message "watcher_armed instance=$instance_id output_root=$output_root"

while true; do
  complete_count=$(find "$output_root/fits" -type f -name FIT_COMPLETE.txt \
    2>/dev/null | wc -l)
  error_count=$(find "$output_root/fits" -type f -name FIT_ERROR.txt \
    2>/dev/null | wc -l)
  if [ "$error_count" -gt 0 ]; then
    log_message "auto_stop_aborted reason=fit_error complete=$complete_count errors=$error_count"
    printf 'status=AUTO_STOP_ABORTED_FIT_ERROR\ncomplete=%s\nerrors=%s\n' \
      "$complete_count" "$error_count" > "$output_root/AUTO_STOP_ABORTED.txt"
    exit 2
  fi

  selection_frozen=false
  experiment_complete=false
  runner_active=false
  [ -f "$output_root/truth_free_selection/TRUTH_FREE_SELECTION_FROZEN.txt" ] &&
    selection_frozen=true
  [ -f "$output_root/FITS_AND_SELECTION_COMPLETE.txt" ] &&
    experiment_complete=true
  pgrep -f "$runner_path" >/dev/null 2>&1 && runner_active=true

  if [ "$complete_count" -eq 72 ] &&
     [ "$selection_frozen" = true ] &&
     [ "$experiment_complete" = true ] &&
     [ "$runner_active" = false ]; then
    log_message "completion_verified complete=72 errors=0 selection_frozen=true runner_active=false"
    printf 'status=AUTO_STOP_REQUESTED\nrequested_at=%s\ninstance_id=%s\nstopped_mode=StopCharging\n' \
      "$(date --iso-8601=seconds)" "$instance_id" \
      > "$output_root/AUTO_STOP_REQUESTED.txt"
    sync
    if aliyun ecs StopInstance \
      --profile "$profile_name" \
      --RegionId "$region_id" \
      --InstanceId "$instance_id" \
      --StoppedMode StopCharging \
      --ForceStop false \
      --DryRun false > "$output_root/AUTO_STOP_API_RESPONSE.txt" 2>&1; then
      log_message "stop_api_accepted stopped_mode=StopCharging"
      printf 'status=STOP_API_ACCEPTED\naccepted_at=%s\n' \
        "$(date --iso-8601=seconds)" \
        > "$output_root/AUTO_STOP_API_ACCEPTED.txt"
      sync
      exit 0
    else
      api_status=$?
      log_message "auto_stop_failed reason=stop_api exit_status=$api_status"
      printf 'status=AUTO_STOP_FAILED\nexit_status=%s\nfailed_at=%s\n' \
        "$api_status" "$(date --iso-8601=seconds)" \
        > "$output_root/AUTO_STOP_FAILED.txt"
      exit "$api_status"
    fi
  fi
  sleep 60
done

