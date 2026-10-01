#!/usr/bin/env bash
# One S continuation, blind freeze, fixed scientific evaluation, cost stop.
set -euo pipefail

exec_dir=/root/stage7k-exec-20260930-v1
parent=/root/multifsync-results
s_root="$parent/stage7k_seed02_s_20260930_v1"
freeze_root="$parent/stage7k_seed02_blind_freeze_20260930_v1"
record="$parent/stage7k_seed02_supervisor_20260930_v1"
evaluation_root="$parent/stage7k_seed02_scientific_evaluation_20260930_v1"
u_root="$parent/stage7j_seed02_u_20260930_v1"
a_root="$parent/stage7j_seed02_a_20260930_v1"

test ! -e "$record" && test ! -e "$s_root" &&
  test ! -e "$freeze_root" && test ! -e "$evaluation_root"
mkdir -m 0700 "$record"
exec > >(tee -a "$record/RUN.log") 2>&1
last_step=initial_preflight
truth_opened=FALSE

finish() {
  status=$?
  trap - EXIT
  printf 'status=%s\nlast_step=%s\nfinished_at=%s\ntruth_opened=%s\n' \
    "$status" "$last_step" "$(date --iso-8601=seconds)" "$truth_opened" \
    > "$record/SUPERVISOR_FINISHED.txt"
  sync
  systemd-run --unit=stage7k-seed02-cost-stop-20260930-v1 \
    --setenv=HOME=/root /bin/bash \
    "$exec_dir/run_stage7k_seed02_autostop_20260930_v1.sh" \
    "$record/SUPERVISOR_FINISHED.txt" || true
  exit "$status"
}
trap finish EXIT

cd "$exec_dir"
sha256sum --check SOURCE_FILES.sha256
test -x /usr/bin/timeout
test -f "$a_root/STATE_S_SWEEP50.rds"
test -f "$u_root/u_parity_v1/PARITY_PASS.txt"

runner="$exec_dir/run_stage7k_s_resume_20260930_v1.R"
freeze="$exec_dir/freeze_stage7k_blind_endpoints_20260930_v1.R"
evaluator="$exec_dir/evaluate_stage7k_seed02_endpoints_20260930_v1.R"
support_repo=/root/multiFSYNC-covariate-routeb-pilot-20260927-v1
package_lib=/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1
stage7c_evaluator=/root/multiFSYNC-stage7c-20260927-v1/experiments/v0L-V/covariate-routeb-stage7c-independent-confirmation-20260927/evaluate_covariate_routeb_stage7c_20260928_v1.R
eval_args=("$parent" "$support_repo" "$package_lib" "$stage7c_evaluator")

last_step=s_preflight
Rscript "$runner" preflight "$parent" "$s_root"
last_step=s_449_sweeps
timeout --signal=TERM --kill-after=60s 10800s \
  Rscript "$runner" fit "$parent" "$s_root"
last_step=blind_freeze
Rscript "$freeze" "$u_root" "$a_root" "$s_root" "$freeze_root"
test -f "$freeze_root/BLIND_ENDPOINTS_FROZEN.txt"
printf 'status=0\nlast_step=blind_freeze\ntruth_opened=FALSE\n' \
  > "$record/BLIND_PHASE_FINISHED.txt"
printf 'STAGE7K_SEED02_BLIND_COMPLETE\n' > "$record/BLIND_COMPLETE.txt"
sync

last_step=science_preflight
Rscript "$evaluator" "${eval_args[@]}" preflight
last_step=scientific_evaluation
truth_opened=TRUE
Rscript "$evaluator" "${eval_args[@]}"
test -f "$evaluation_root/EVALUATION_COMPLETE.txt"
printf 'STAGE7K_SEED02_EXPERIMENT_COMPLETE\n' > "$record/COMPLETE.txt"
last_step=complete
