#!/usr/bin/env bash
set -euo pipefail
exec_dir=/root/stage7l-exec-20261001-v1
record=/root/multifsync-results/stage7l_supervisor_20261001_v1
result=/root/multifsync-results/stage7l_prospective_20261001_v1
test ! -e "$record"
test ! -e "$result"
mkdir -m 0700 "$record"
exec > >(tee -a "$record/RUN.log") 2>&1
finish() {
  status=$?
  trap - EXIT
  printf 'status=%s\nfinished_at=%s\nresult=%s\n' "$status" "$(date --iso-8601=seconds)" "$result" \
    > "$record/SUPERVISOR_FINISHED.txt"
  sync
  systemd-run --unit=stage7l-cost-stop-20261001-v1 --setenv=HOME=/root \
    /bin/bash "$exec_dir/autostop_stage7l.sh" "$record/SUPERVISOR_FINISHED.txt" || true
  exit "$status"
}
trap finish EXIT
cd "$exec_dir"
sha256sum --check SOURCE_FILES.sha256
test -f PREFLIGHT_PASS.txt
test -f ADAPTER_MICRO_PASS.txt
python3 schedule_stage7l.py run
test -f "$result/EXPERIMENT_COMPLETE.txt"
printf 'STAGE7L_COMPLETE\n' > "$record/COMPLETE.txt"
