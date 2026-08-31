#!/usr/bin/env bash
set -euo pipefail

stage6a_repo="${STAGE6A_REPO:-/root/multiFSYNC-stage6a-20260831}"
stage6a_lib="${STAGE6A_R_LIB:-/root/multiFSYNC-stage6a-lib-20260831}"
stage6a_output="${STAGE6A_OUTPUT_ROOT:-/root/g12-dimension-overspec-stage6a-20260831-v1}"
stage6a_dir="$stage6a_repo/experiments/v0L-V/g12-dimension-overspec-stage6a-20260831"

test -d "$stage6a_repo/R"
test -f "$stage6a_dir/stage6a_runner_20260831_v1.R"

available_kb=$(df -Pk /root | awk 'NR==2 {print $4}')
if [ "$available_kb" -lt 15728640 ]; then
  echo "ERROR: at least 15 GiB free disk is required; available_kb=$available_kb" >&2
  exit 2
fi

mkdir -p "$stage6a_lib" "$stage6a_output"
R CMD INSTALL --library="$stage6a_lib" "$stage6a_repo"

export STAGE6A_REPO="$stage6a_repo"
export STAGE6A_R_LIB="$stage6a_lib"
export STAGE6A_OUTPUT_ROOT="$stage6a_output"
Rscript "$stage6a_dir/stage6a_runner_20260831_v1.R" --action=check

Rscript -e '.libPaths(c(Sys.getenv("STAGE6A_R_LIB"), .libPaths())); cat("R=", R.version.string, "\n", sep=""); cat("multiFSYNC=", as.character(packageVersion("multiFSYNC")), "\n", sep=""); for (p in c("digest","abind","ellipse","gtools","magic","matrixcalc","matrixStats","pracma","MASS","lattice","splines")) cat(p,"=",as.character(packageVersion(p)),"\n",sep="")' \
  > "$stage6a_output/SERVER_R_PACKAGE_VERSIONS.txt"

{
  echo "status=SERVER_PREPARATION_COMPLETE"
  echo "prepared_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo=$stage6a_repo"
  echo "library=$stage6a_lib"
  echo "output_root=$stage6a_output"
  echo "workers=3"
  echo "n_cpus_per_fit=1"
  echo "truth_unseal_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$stage6a_output/SERVER_PREPARATION_COMPLETE.txt"

echo "SERVER_PREPARATION_COMPLETE output_root=$stage6a_output"
