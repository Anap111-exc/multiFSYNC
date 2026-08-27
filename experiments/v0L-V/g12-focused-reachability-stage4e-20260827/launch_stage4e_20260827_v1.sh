#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4E_ROOT:-/root/v0lv-g12-focused-reachability-stage4e-20260827-v1}"
runner_commit="${G12SS4E_RUNNER_COMMIT:?G12SS4E_RUNNER_COMMIT is required}"
cd "$experiment_root"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

test ! -e STAGE4E_COMPLETE.txt
test ! -e FITS_COMPLETE.txt
test ! -e RUNNER.pid

if [[ ! -f DATA_REUSE_COMPLETE.txt ]]; then
  bash prepare_stage4e_data_20260827_v1.sh
fi
Rscript s4e_nonformal_smoke_20260827_v1.R
Rscript s4e_supervisor_20260827_v1.R --check-only=true

nohup bash run_stage4e_20260827_v1.sh > STAGE4E_RUN.log 2>&1 &
runner_pid=$!
printf '%s\n' "$runner_pid" > RUNNER.pid
printf 'status=LAUNCHED\nrunner_pid=%s\nlaunched_utc=%s\nrunner_commit=%s\npackage_code_commit=72d9a53f0ef5e9d3e5f49d9cf837958207b3c980\nper_fit_n_cpus=1\nblas_threads=1\nouter_workers=4\nfixed_T1=400\nnew_fits=18\ntotal_endpoints_after_reuse=24\ntruth_free_selection_before_evaluation=TRUE\nformal_v0lv_result=FALSE\n' \
  "$runner_pid" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$runner_commit" \
  > STAGE4E_LAUNCH_RECORD.txt
printf 'STAGE4E_LAUNCHED pid=%s new_fits=18 workers=4 fixed_t1=400\n' "$runner_pid"
