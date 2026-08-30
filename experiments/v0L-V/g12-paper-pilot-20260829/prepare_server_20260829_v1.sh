#!/usr/bin/env bash
set -euo pipefail

pilot_repo="${PAPER_PILOT_REPO:-/root/multiFSYNC-g12-paper-pilot-20260829}"
pilot_lib="${PAPER_PILOT_R_LIB:-/root/multiFSYNC-g12-paper-pilot-lib-20260829}"
pilot_output="${PAPER_PILOT_OUTPUT_ROOT:-/root/g12-paper-pilot-20260829-v1}"
pilot_dir="$pilot_repo/experiments/v0L-V/g12-paper-pilot-20260829"

test -d "$pilot_repo/R"
test -f "$pilot_dir/paper_pilot_runner_20260829_v1.R"

available_kb=$(df -Pk /root | awk 'NR==2 {print $4}')
if [ "$available_kb" -lt 15728640 ]; then
  echo "ERROR: at least 15 GiB free disk is required; available_kb=$available_kb" >&2
  exit 2
fi

mkdir -p "$pilot_lib" "$pilot_output"
R CMD INSTALL --library="$pilot_lib" "$pilot_repo"

export PAPER_PILOT_REPO="$pilot_repo"
export PAPER_PILOT_R_LIB="$pilot_lib"
export PAPER_PILOT_OUTPUT_ROOT="$pilot_output"
Rscript "$pilot_dir/paper_pilot_runner_20260829_v1.R" --action=check

Rscript -e '.libPaths(c(Sys.getenv("PAPER_PILOT_R_LIB"), .libPaths())); cat("R=", R.version.string, "\n", sep=""); cat("multiFSYNC=", as.character(packageVersion("multiFSYNC")), "\n", sep=""); for (p in c("digest","abind","ellipse","gtools","magic","matrixcalc","matrixStats","pracma","MASS","lattice","splines")) cat(p,"=",as.character(packageVersion(p)),"\n",sep="")' \
  > "$pilot_output/SERVER_R_PACKAGE_VERSIONS.txt"

{
  echo "status=SERVER_PREPARATION_COMPLETE"
  echo "prepared_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo=$pilot_repo"
  echo "library=$pilot_lib"
  echo "output_root=$pilot_output"
  echo "workers=3"
  echo "n_cpus_per_fit=1"
  echo "truth_unseal_authorized=FALSE"
  echo "formal_paper_mc_result=FALSE"
} > "$pilot_output/SERVER_PREPARATION_COMPLETE.txt"

echo "SERVER_PREPARATION_COMPLETE output_root=$pilot_output"
