#!/usr/bin/env bash
set -euo pipefail

pilot_output="${PAPER_PILOT_OUTPUT_ROOT:-/root/g12-paper-pilot-20260829-v1}"

echo "=== timestamp ==="
date -u +%Y-%m-%dT%H:%M:%SZ
echo "=== current status ==="
if [ -f "$pilot_output/CURRENT_STATUS.csv" ]; then
  cat "$pilot_output/CURRENT_STATUS.csv"
else
  echo "CURRENT_STATUS_NOT_CREATED"
fi
echo "=== terminal counts ==="
complete=$(find "$pilot_output/fits" -name FIT_COMPLETE.txt -type f 2>/dev/null | wc -l)
errors=$(find "$pilot_output/fits" -name FIT_ERROR.txt -type f 2>/dev/null | wc -l)
echo "complete=$complete"
echo "errors=$errors"
echo "=== supervisor process ==="
if [ -f "$pilot_output/PILOT_LAUNCHED.txt" ]; then
  pilot_pid=$(awk -F= '$1=="pid" {print $2}' "$pilot_output/PILOT_LAUNCHED.txt")
  if kill -0 "$pilot_pid" 2>/dev/null; then
    ps -o pid,ppid,etime,%cpu,%mem,rss,vsz,cmd -p "$pilot_pid"
  else
    echo "supervisor_not_running pid=$pilot_pid"
  fi
else
  echo "PILOT_NOT_LAUNCHED"
fi
echo "=== active R workers ==="
pgrep -af '[R]script.*paper_pilot_runner_20260829_v1.R.*action=fit' || true
echo "=== disk ==="
df -h /root
du -sh "$pilot_output" 2>/dev/null || true
echo "=== truth boundary ==="
if [ -f "$pilot_output/TRUTH_UNSEAL_AUTHORIZATION.txt" ]; then
  echo "TRUTH_UNSEAL_AUTHORIZATION_PRESENT"
else
  echo "truth_unseal_authorized=FALSE"
fi
if [ -f "$pilot_output/WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt" ]; then
  cat "$pilot_output/WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt"
fi
