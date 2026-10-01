#!/usr/bin/env bash
# Stage 7I execution chain. Run only after explicit authorization and ECS input
# preflight. Every R action refuses existing outputs; no truth is opened here.
set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 STAGE7C_ROOT STAGE7G_ROOT STAGE7H_ROOT STAGE7G_SOURCE_ARCHIVE RESULTS_PARENT" >&2
  exit 2
fi

stage7c_root="$1"
stage7g_root="$2"
stage7h_root="$3"
stage7g_archive="$4"
results_parent="$5"
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
stage7i_archive="$script_dir/multiFSYNC-stage7i-20260929-v1.tar.gz"
stage7g_runner="$script_dir/run_stage7g_state_snapshots_20260928_v1.R"
stage7i_lib="$script_dir/stage7i_lib_ecs_20260929_v1"
u_root="$results_parent/stage7i_u_20260929_v1"
as_root="$results_parent/stage7i_as_20260929_v1"

for required in "$stage7c_root" "$stage7g_root" "$stage7h_root" \
                "$results_parent"; do
  test -d "$required" || { echo "Missing directory: $required" >&2; exit 3; }
done
for required in "$stage7g_archive" "$stage7i_archive" "$stage7g_runner"; do
  test -f "$required" || { echo "Missing file: $required" >&2; exit 3; }
done
test ! -e "$stage7i_lib" || { echo "Refusing existing package library" >&2; exit 3; }
test ! -e "$u_root" || { echo "Refusing existing U output" >&2; exit 3; }
test ! -e "$as_root" || { echo "Refusing existing A/S output" >&2; exit 3; }

mkdir "$stage7i_lib"
R CMD INSTALL --library="$stage7i_lib" "$stage7i_archive"

u_runner="$script_dir/run_stage7i_u_resume_20260929_v1.R"
Rscript "$u_runner" prepare "$stage7c_root" "$stage7g_root" \
  "$stage7g_runner" "$stage7g_archive" "$stage7h_root" \
  "$stage7i_lib" "$stage7i_archive" "$u_root"
Rscript "$u_runner" smoke "$stage7c_root" "$stage7g_root" \
  "$stage7g_runner" "$stage7g_archive" "$stage7h_root" \
  "$stage7i_lib" "$stage7i_archive" "$u_root"
Rscript "$u_runner" u "$stage7c_root" "$stage7g_root" \
  "$stage7g_runner" "$stage7g_archive" "$stage7h_root" \
  "$stage7i_lib" "$stage7i_archive" "$u_root"
Rscript "$u_runner" verify "$stage7c_root" "$stage7g_root" \
  "$stage7g_runner" "$stage7g_archive" "$stage7h_root" \
  "$stage7i_lib" "$stage7i_archive" "$u_root"

as_runner="$script_dir/run_stage7i_as_resume_20260929_v1.R"
Rscript "$as_runner" prepare "$u_root" "$as_root"
Rscript "$as_runner" fit "$u_root" "$as_root"

test -f "$u_root/u_parity_v1/PARITY_PASS.txt"
test -f "$as_root/AS_FITS_COMPLETE.txt"
printf 'STAGE7I_BLIND_FITS_COMPLETE\nU_ROOT=%s\nAS_ROOT=%s\n' \
  "$u_root" "$as_root"
