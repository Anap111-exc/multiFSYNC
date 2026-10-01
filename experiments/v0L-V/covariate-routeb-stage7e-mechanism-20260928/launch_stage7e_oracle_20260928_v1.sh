#!/usr/bin/env bash
set -euo pipefail

cd /root/multifsync-results/stage7e_20260928_v1
export R_LIBS=/root/multiFSYNC-stage7c-lib-20260927-v1
stage7c=../covariate_routeb_stage7c_20260927_v1
output=oracle_v1

test -f "$output/PREPARE_COMPLETE.txt"
test -f "$output/MICRO_SMOKE_PASS.txt"
test ! -e "$output/FITS_STARTED.txt"
test ! -e "$output/SUPERVISOR_PID.txt"
Rscript run_stage7e_oracle_20260928_v1.R preflight "$stage7c" "$output"
nohup Rscript run_stage7e_oracle_20260928_v1.R \
  supervisor "$stage7c" "$output" \
  > "$output/SUPERVISOR.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$output/SUPERVISOR_PID.txt"
printf 'STAGE7E_SUPERVISOR_LAUNCHED pid=%s\n' "$(cat "$output/SUPERVISOR_PID.txt")"
