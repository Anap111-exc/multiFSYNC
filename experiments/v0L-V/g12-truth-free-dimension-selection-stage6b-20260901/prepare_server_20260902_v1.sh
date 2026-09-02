#!/usr/bin/env bash
set -euo pipefail

stage6b_repo="${STAGE6B_REPO:-/root/multiFSYNC-stage6b-20260902}"
stage6b_lib="${STAGE6B_R_LIB:-/root/multiFSYNC-stage6a-lib-20260831}"
stage6b_output="${STAGE6B_OUTPUT_ROOT:-/root/g12-truth-free-dimension-selection-stage6b-20260902-v1}"
stage6b_observation="${STAGE6B_OBSERVATION_BUNDLE:-/root/g12-dimension-overspec-stage6a-20260831-v1/data/g12s6a_01/observation_bundle.rds}"
stage6b_commit="${STAGE6B_SOURCE_COMMIT:-UNKNOWN}"
stage6b_dir="$stage6b_repo/experiments/v0L-V/g12-truth-free-dimension-selection-stage6b-20260901"
runner="$stage6b_dir/run_stage6b_predictive_pilot_20260902_v1.R"

test -d "$stage6b_repo/R"
test -d "$stage6b_lib"
test -f "$stage6b_observation"
test -f "$runner"
test -f "$stage6b_dir/FIT_CONFIGS_STAGE6B_PILOT_V1.csv"
test -f "$stage6b_dir/START_SEEDS_STAGE6B_PILOT_V1.csv"

available_kb=$(df -Pk /root | awk 'NR==2 {print $4}')
if [ "$available_kb" -lt 5242880 ]; then
  echo "ERROR: at least 5 GiB free disk is required; available_kb=$available_kb" >&2
  exit 2
fi

mkdir -p "$stage6b_output"
export STAGE6B_REPO="$stage6b_repo"
export STAGE6B_R_LIB="$stage6b_lib"
export STAGE6B_OUTPUT_ROOT="$stage6b_output"
export STAGE6B_OBSERVATION_BUNDLE="$stage6b_observation"
export STAGE6B_WORKERS="${STAGE6B_WORKERS:-3}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=prepare \
  --output-root="$stage6b_output" \
  --observation-bundle="$stage6b_observation"

Rscript -e '.libPaths(c(Sys.getenv("STAGE6B_R_LIB"), .libPaths())); cat("R=", R.version.string, "\n", sep=""); cat("multiFSYNC=", as.character(packageVersion("multiFSYNC")), "\n", sep=""); for (p in c("digest","abind","ellipse","gtools","magic","matrixcalc","matrixStats","pracma","MASS","lattice","splines")) cat(p,"=",as.character(packageVersion(p)),"\n",sep="")' \
  > "$stage6b_output/SERVER_R_PACKAGE_VERSIONS.txt"

{
  echo "status=SERVER_PREPARATION_COMPLETE"
  echo "prepared_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "source_commit=$stage6b_commit"
  echo "repo=$stage6b_repo"
  echo "library=$stage6b_lib"
  echo "output_root=$stage6b_output"
  echo "source_observation=$stage6b_observation"
  echo "workers=$STAGE6B_WORKERS"
  echo "n_cpus_per_fit=1"
  echo "new_data_generated=FALSE"
  echo "truth_access_authorized=FALSE"
  echo "continuation_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$stage6b_output/SERVER_PREPARATION_COMPLETE.txt"

echo "SERVER_PREPARATION_COMPLETE output_root=$stage6b_output"
