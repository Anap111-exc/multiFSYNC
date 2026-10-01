#!/usr/bin/env bash
# Finalize and request StopCharging only after all truth-blind Stage 7I gates pass.
# Run once as a separate systemd service while the fit service is active.
set -euo pipefail

fit_unit=stage7i-uas-20260929-v1.service
script_dir=/root/stage7i-exec-20260929-v1
results_parent=/root/multifsync-results
freeze_script="$script_dir/freeze_stage7i_blind_endpoints_20260929_v1.R"
expected_freeze_sha=2748eca82ec0f4c10c2e95edd077dcbe2cc01069148ef0483afca38fbd7622b5
record_root="$results_parent/stage7i_autostop_20260929_v1"
freeze_marker="$results_parent/stage7i_blind_freeze_20260929_v1/BLIND_ENDPOINTS_FROZEN.txt"

log() { printf '%s %s\n' "$(date --iso-8601=seconds)" "$*"; }
property() {
  local value
  value="$(systemctl show "$fit_unit" -p "$1")"
  printf '%s\n' "${value#*=}"
}

test -f "$freeze_script"
actual_freeze_sha="$(sha256sum "$freeze_script")"
actual_freeze_sha="${actual_freeze_sha%% *}"
test "$actual_freeze_sha" = "$expected_freeze_sha" || {
  log 'Freeze script hash mismatch; leaving ECS running.' >&2
  exit 3
}
test ! -e "$record_root" || {
  log 'Auto-stop record already exists; refusing a second watcher.' >&2
  exit 3
}
log "Waiting for $fit_unit to complete."

while :; do
  load_state="$(property LoadState)"
  active_state="$(property ActiveState)"
  result="$(property Result)"
  main_status="$(property ExecMainStatus)"
  test "$load_state" = loaded || {
    log "Fit unit unavailable ($load_state); leaving ECS running." >&2
    exit 4
  }
  case "$active_state" in
    active|activating)
      sleep 60
      ;;
    inactive)
      test "$result" = success && test "$main_status" = 0 || {
        log "Fit unit did not succeed (Result=$result, status=$main_status); leaving ECS running." >&2
        exit 4
      }
      break
      ;;
    *)
      log "Fit unit ended in $active_state (Result=$result); leaving ECS running." >&2
      exit 4
      ;;
  esac
done

test -f "$results_parent/stage7i_u_20260929_v1/u_parity_v1/PARITY_PASS.txt"
test -f "$results_parent/stage7i_as_20260929_v1/AS_FITS_COMPLETE.txt"
test ! -e "$results_parent/stage7i_blind_freeze_20260929_v1"
mkdir "$record_root"
watcher_sha="$(sha256sum "$0")"
watcher_sha="${watcher_sha%% *}"
printf 'fit_unit=%s\nfit_result=%s\nfit_exit_status=%s\nfreeze_script_sha256=%s\nwatcher_script_sha256=%s\n' \
  "$fit_unit" "$result" "$main_status" "$expected_freeze_sha" "$watcher_sha" \
  > "$record_root/START_GATES.txt"

log 'Rechecking U/A/S and freezing blind endpoints.'
Rscript "$freeze_script" "$results_parent" \
  > "$record_root/BLIND_FREEZE_STDOUT.txt" \
  2> "$record_root/BLIND_FREEZE_STDERR.txt"
test -f "$freeze_marker"
grep -qx 'STAGE7I_BLIND_ENDPOINTS_FROZEN' "$freeze_marker"
marker_sha="$(sha256sum "$freeze_marker")"
marker_sha="${marker_sha%% *}"
printf 'status=BLIND_FREEZE_COMPLETE\ncompleted_at=%s\nfreeze_marker_sha256=%s\n' \
  "$(date --iso-8601=seconds)" "$marker_sha" \
  > "$record_root/BLIND_FREEZE_COMPLETE.txt"
sync

log 'Checking StopInstance authorization with DryRun.'
if aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun true \
    > "$record_root/STOP_DRY_RUN_RESPONSE.txt" 2>&1; then
  dry_status=0
else
  dry_status=$?
fi
printf 'dry_run_exit_status=%s\n' "$dry_status" \
  > "$record_root/STOP_DRY_RUN_STATUS.txt"
grep -q 'DryRunOperation' "$record_root/STOP_DRY_RUN_RESPONSE.txt" || {
  log 'DryRun did not confirm the requested StopInstance; leaving ECS running.' >&2
  exit 5
}

printf 'status=STOP_REQUEST_INTENDED\nrequested_at=%s\ninstance_id=i-bp1egmsairexm9w1rnks\nregion=cn-hangzhou\nstopped_mode=StopCharging\nforce_stop=false\n' \
  "$(date --iso-8601=seconds)" > "$record_root/STOP_REQUEST_INTENDED.txt"
sync
log 'Requesting cost-saving StopCharging mode.'
if aliyun ecs StopInstance --profile multifsync-self-stop \
    --RegionId cn-hangzhou --InstanceId i-bp1egmsairexm9w1rnks \
    --StoppedMode StopCharging --ForceStop false --DryRun false \
    > "$record_root/STOP_API_RESPONSE.txt" 2>&1; then
  api_status=0
else
  api_status=$?
fi
printf 'api_exit_status=%s\nresponse_recorded_at=%s\n' \
  "$api_status" "$(date --iso-8601=seconds)" \
  > "$record_root/STOP_API_STATUS.txt"
sync
test "$api_status" -eq 0 || {
  log "StopInstance API failed (status=$api_status); check the saved response." >&2
  exit "$api_status"
}
log 'StopInstance API accepted; final instance state requires independent verification.'
