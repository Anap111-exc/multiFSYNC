#!/usr/bin/env bash
set -euo pipefail

stage6a_output="${STAGE6A_OUTPUT_ROOT:-/root/g12-dimension-overspec-stage6a-20260831-v1}"

echo "=== timestamp ==="
date -u +%Y-%m-%dT%H:%M:%SZ
echo "=== current status ==="
if [ -f "$stage6a_output/CURRENT_STATUS.csv" ]; then
  cat "$stage6a_output/CURRENT_STATUS.csv"
else
  echo "CURRENT_STATUS_NOT_CREATED"
fi
echo "=== capacity probe ==="
if [ -f "$stage6a_output/CAPACITY_PROBE_COMPLETE.txt" ]; then
  cat "$stage6a_output/CAPACITY_PROBE_COMPLETE.txt"
elif [ -f "$stage6a_output/CAPACITY_PROBE_ERROR.txt" ]; then
  cat "$stage6a_output/CAPACITY_PROBE_ERROR.txt"
else
  echo "CAPACITY_PROBE_PENDING"
fi
echo "=== terminal counts ==="
complete=$(find "$stage6a_output/fits" -name FIT_COMPLETE.txt -type f 2>/dev/null | wc -l)
errors=$(find "$stage6a_output/fits" -name FIT_ERROR.txt -type f 2>/dev/null | wc -l)
echo "complete=$complete"
echo "errors=$errors"
echo "=== supervisor process ==="
if [ -f "$stage6a_output/STAGE6A_LAUNCHED.txt" ]; then
  stage6a_pid=$(awk -F= '$1=="pid" {print $2}' "$stage6a_output/STAGE6A_LAUNCHED.txt")
  if kill -0 "$stage6a_pid" 2>/dev/null; then
    ps -o pid,ppid,etime,%cpu,%mem,rss,vsz,cmd -p "$stage6a_pid"
  else
    echo "supervisor_not_running pid=$stage6a_pid"
  fi
else
  echo "STAGE6A_NOT_LAUNCHED"
fi
echo "=== active R workers ==="
pgrep -af '[R]script.*stage6a_runner_20260831_v1.R.*action=fit' || true
echo "=== disk ==="
df -h /root
du -sh "$stage6a_output" 2>/dev/null || true
echo "=== truth boundary ==="
if [ -f "$stage6a_output/TRUTH_UNSEAL_AUTHORIZATION.txt" ]; then
  echo "TRUTH_UNSEAL_AUTHORIZATION_PRESENT"
else
  echo "truth_unseal_authorized=FALSE"
fi
if [ -f "$stage6a_output/WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt" ]; then
  cat "$stage6a_output/WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt"
fi
