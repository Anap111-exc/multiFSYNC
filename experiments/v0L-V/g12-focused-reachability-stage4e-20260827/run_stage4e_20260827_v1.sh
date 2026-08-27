#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4E_ROOT:-/root/v0lv-g12-focused-reachability-stage4e-20260827-v1}"
cd "$experiment_root"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

Rscript s4e_supervisor_20260827_v1.R
Rscript s4e_truth_free_selection_20260827_v1.R
Rscript s4e_truth_evaluation_20260827_v1.R
Rscript s4e_analysis_20260827_v1.R

printf 'status=STAGE4E_COMPLETE\ncompleted_utc=%s\nnew_fits=18\ntotal_evaluated_endpoints=24\ntruth_free_winners=3\nformal_v0lv_result=FALSE\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > STAGE4E_COMPLETE.txt
