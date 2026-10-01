#!/usr/bin/env bash
# Sequential truth-blind U gate -> A replication -> blind freeze.
set -euo pipefail

exec_dir=/root/stage7j-exec-20260930-v1
parent=/root/multifsync-results
u_root="$parent/stage7j_seed02_u_20260930_v1"
a_root="$parent/stage7j_seed02_a_20260930_v1"
freeze_root="$parent/stage7j_seed02_blind_freeze_20260930_v1"
record="$parent/stage7j_seed02_supervisor_20260930_v1"

test ! -e "$record"
mkdir -m 0700 "$record"
exec > >(tee -a "$record/RUN.log") 2>&1
last_step=initial_preflight

finish() {
  status=$?
  trap - EXIT
  printf 'status=%s\nlast_step=%s\nfinished_at=%s\ntruth_opened=FALSE\n' \
    "$status" "$last_step" "$(date --iso-8601=seconds)" \
    > "$record/SUPERVISOR_FINISHED.txt"
  sync
  systemd-run --unit=stage7j-seed02-cost-stop-20260930-v1 \
    --setenv=HOME=/root /bin/bash \
    "$exec_dir/run_stage7j_seed02_autostop_20260930_v1.sh" \
    "$record/SUPERVISOR_FINISHED.txt" || true
  exit "$status"
}
trap finish EXIT

cd "$exec_dir"
sha256sum --check SOURCE_FILES.sha256
test ! -e "$u_root" && test ! -e "$a_root" && test ! -e "$freeze_root"

u_runner="$exec_dir/run_stage7j_u_resume_20260930_v1.R"
a_runner="$exec_dir/run_stage7j_a_resume_20260930_v1.R"
freeze_runner="$exec_dir/freeze_stage7j_blind_endpoints_20260930_v1.R"
stage7c="$parent/covariate_routeb_stage7c_20260927_v1"
stage7g="$parent/stage7g_20260928_v1"
stage7g_runner=/root/stage7i-exec-20260929-v1/run_stage7g_state_snapshots_20260928_v1.R
stage7g_archive=/root/multiFSYNC-stage7g-20260928-v1.tar.gz
stage7h="$parent/stage7h_20260929_v1"
stage7i_library=/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1
stage7i_archive=/root/stage7i-exec-20260929-v1/multiFSYNC-stage7i-20260929-v1.tar.gz

u_args=("$stage7c" "$stage7g" "$stage7g_runner" "$stage7g_archive"
        "$stage7h" "$stage7i_library" "$stage7i_archive" "$u_root")

last_step=u_prepare
Rscript "$u_runner" prepare "${u_args[@]}"
last_step=u_449_sweeps
Rscript "$u_runner" u "${u_args[@]}"
last_step=u_full_parity
Rscript "$u_runner" verify "${u_args[@]}"

last_step=a_fixed_initialization
Rscript "$a_runner" prepare "$u_root" "$a_root"
last_step=a_449_sweeps
Rscript "$a_runner" fit "$u_root" "$a_root"

last_step=blind_endpoint_freeze
Rscript "$freeze_runner" "$u_root" "$a_root" "$freeze_root"
test -f "$freeze_root/BLIND_ENDPOINTS_FROZEN.txt"
printf 'STAGE7J_SEED02_BLIND_EXPERIMENT_COMPLETE\ntruth_opened=FALSE\n' \
  > "$record/COMPLETE.txt"
last_step=complete
