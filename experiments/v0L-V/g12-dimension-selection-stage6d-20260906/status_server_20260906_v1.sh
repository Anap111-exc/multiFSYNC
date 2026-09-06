#!/usr/bin/env bash
set -euo pipefail

stage6d_output="${STAGE6D_OUTPUT_ROOT:-/root/g12-dimension-selection-stage6d-20260906-v1}"

echo "=== stage6d markers ==="
for marker in SERVER_PREPARATION_COMPLETE.txt STAGE6D_LAUNCHED.txt \
  STAGE6D_PREPARED.txt STAGE6D_COMPLETE.txt; do
  if [[ -f "$stage6d_output/$marker" ]]; then
    echo "$marker=present"
  else
    echo "$marker=absent"
  fi
done

complete=0
errors=0
fit_dirs=0
if [[ -d "$stage6d_output/fits" ]]; then
  complete=$(find "$stage6d_output/fits" -name FIT_COMPLETE.txt -type f | wc -l)
  errors=$(find "$stage6d_output/fits" -name FIT_ERROR.txt -type f | wc -l)
  fit_dirs=$(find "$stage6d_output/fits" -mindepth 1 -maxdepth 1 -type d | wc -l)
fi
echo "fit_directories=$fit_dirs"
echo "complete_fits=$complete"
echo "error_fits=$errors"
echo "in_progress_or_incomplete=$((fit_dirs-complete-errors))"

echo "=== supervisor ==="
pgrep -af '[r]un_stage6d_20260906_v1.R|[S]TAGE6D_LAUNCH' || true
echo "=== resources ==="
uptime
free -h
df -h / "$stage6d_output" 2>/dev/null | awk '!seen[$1]++'
du -sh "$stage6d_output" 2>/dev/null || true

if [[ -f "$stage6d_output/CURRENT_STATUS.csv" ]]; then
  echo "=== current status ==="
  cat "$stage6d_output/CURRENT_STATUS.csv"
fi
if [[ -s "$stage6d_output/STAGE6D_LAUNCH.stderr.log" ]]; then
  echo "=== stderr tail ==="
  tail -n 20 "$stage6d_output/STAGE6D_LAUNCH.stderr.log"
fi
