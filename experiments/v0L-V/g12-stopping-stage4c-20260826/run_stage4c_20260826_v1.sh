#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4C_ROOT:-/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1}"
cd "$experiment_root"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

Rscript s4c_supervisor_20260826_v1.R
Rscript s4c_analysis_20260826_v1.R
