#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4E_ROOT:-/root/v0lv-g12-focused-reachability-stage4e-20260827-v1}"
cd "$experiment_root"

complete=$(find fits -name FIT_COMPLETE.txt -type f 2>/dev/null | wc -l)
errors=$(find fits -name FIT_ERROR.txt -type f 2>/dev/null | wc -l)
printf 'complete=%s error=%s total_new=18\n' "$complete" "$errors"
if [[ -f CURRENT_STATUS.csv ]]; then tail -n 2 CURRENT_STATUS.csv; fi
if [[ -f TRUTH_FREE_SELECTION_COMPLETE.txt ]]; then
  cat TRUTH_FREE_SELECTION_COMPLETE.txt
fi
if [[ -f STAGE4E_COMPLETE.txt ]]; then cat STAGE4E_COMPLETE.txt; fi
pgrep -af '[R]script.*s4e_|[b]ash run_stage4e_20260827_v1.sh' || true
