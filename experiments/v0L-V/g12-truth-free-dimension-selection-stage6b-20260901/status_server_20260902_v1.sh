#!/usr/bin/env bash
set -euo pipefail

stage6b_output="${STAGE6B_OUTPUT_ROOT:-/root/g12-truth-free-dimension-selection-stage6b-20260902-v1}"

echo "=== timestamp ==="
date -u
echo "=== process ==="
pgrep -af '[R](script| --| CMD)' || true
echo "=== current status ==="
if [ -f "$stage6b_output/CURRENT_STATUS.csv" ]; then
  cat "$stage6b_output/CURRENT_STATUS.csv"
else
  echo "CURRENT_STATUS_MISSING"
fi
echo "=== counts ==="
if [ -d "$stage6b_output/fits" ]; then
  complete=$(find "$stage6b_output/fits" -name FIT_COMPLETE.txt -type f | wc -l)
  errors=$(find "$stage6b_output/fits" -name FIT_ERROR.txt -type f | wc -l)
else
  complete=0
  errors=0
fi
echo "complete_fits=$complete"
echo "error_fits=$errors"
echo "=== selections ==="
if [ -d "$stage6b_output/truth_free_selection" ]; then
  find "$stage6b_output/truth_free_selection" -maxdepth 2 -type f \
    \( -name '*FROZEN.txt' -o -name '*SELECTION.csv' \) -print
else
  echo "NO_SELECTION_YET"
fi
echo "=== launch stdout tail ==="
tail -n 20 "$stage6b_output/STAGE6B_LAUNCH.stdout.log" 2>/dev/null || true
echo "=== launch stderr tail ==="
tail -n 20 "$stage6b_output/STAGE6B_LAUNCH.stderr.log" 2>/dev/null || true
echo "=== resources ==="
free -h
df -h /root
