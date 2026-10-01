#!/usr/bin/env bash
set -euo pipefail

repo_root="${STAGE6G_REPO_ROOT:-/root/multiFSYNC-stage6b-20260903-fa10cbc}"
package_root="${STAGE6G_PACKAGE_SOURCE_ROOT:-/root/multiFSYNC-jaoua-paired-20260909}"
r_library="${STAGE6G_R_LIB:-/root/multiFSYNC-jaoua-paired-lib-20260909}"
output_root="${STAGE6G_OUTPUT_ROOT:-/root/multifsync-results/stage6g_20260911_v1}"
workers="${STAGE6G_WORKERS:-7}"
experiment_root="$repo_root/experiments/v0L-V/g12-dimension-criterion-falsification-stage6g-20260911"
runner="$experiment_root/run_stage6g_20260911_v1.R"

test -f "$runner"
test -d "$package_root/R"
test -d "$r_library"
test ! -e "$output_root"
mkdir -p "$(dirname "$output_root")"

export STAGE6G_PACKAGE_SOURCE_ROOT="$package_root"
export STAGE6G_R_LIB="$r_library"
export STAGE6G_WORKERS="$workers"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1

Rscript "$runner" --action=check
Rscript "$runner" --action=prepare --output-root="$output_root"
exec Rscript "$runner" --action=run --output-root="$output_root"
