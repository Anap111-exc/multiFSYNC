#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4C_ROOT:-/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1}"
runner_commit="${G12SS4C_RUNNER_COMMIT:?G12SS4C_RUNNER_COMMIT is required}"
cd "$experiment_root"

test ! -e STAGE4C_COMPLETE.txt
test ! -e FITS_COMPLETE.txt
test ! -e RUNNER.pid

Rscript s4c_generate_and_seal_20260826_v1.R
Rscript s4c_nonformal_smoke_20260826_v1.R
Rscript s4c_supervisor_20260826_v1.R --check-only=true

nohup bash run_stage4c_20260826_v1.sh > STAGE4C_RUN.log 2>&1 &
runner_pid=$!
printf '%s\n' "$runner_pid" > RUNNER.pid
printf 'status=LAUNCHED\nrunner_pid=%s\nlaunched_utc=%s\nrunner_commit=%s\npackage_code_commit=72d9a53f0ef5e9d3e5f49d9cf837958207b3c980\nfixed_T1=400\nfits=12\ntruth_unseal_authorized=FALSE\nautomatic_800=FALSE\n' \
  "$runner_pid" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$runner_commit" > STAGE4C_LAUNCH_RECORD.txt
printf 'STAGE4C_LAUNCHED pid=%s fits=12 fixed_t1=400 truth=0\n' "$runner_pid"
