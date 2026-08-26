#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4C_ROOT:-/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1}"
cd "$experiment_root"

complete=$(find fits -name FIT_COMPLETE.txt -type f 2>/dev/null | wc -l)
errors=$(find fits -name FIT_ERROR.txt -type f 2>/dev/null | wc -l)
printf 'complete=%s error=%s total=12\n' "$complete" "$errors"
if [[ -f CURRENT_STATUS.csv ]]; then tail -n 2 CURRENT_STATUS.csv; fi
if [[ -f STAGE4C_COMPLETE.txt ]]; then cat STAGE4C_COMPLETE.txt; fi
pgrep -af '[R]script.*s4c_|[b]ash run_stage4c_20260826_v1.sh' || true
