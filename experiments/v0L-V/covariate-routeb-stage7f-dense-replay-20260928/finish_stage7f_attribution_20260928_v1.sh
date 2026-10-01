#!/usr/bin/env bash
set -euo pipefail

readonly OUT=/root/multifsync-results/stage7f_20260928_v1
readonly ANALYZER=/root/multifsync-stage7f-20260928-v1/analyze_stage7f_trace_20260928_v1.R

# File gate only; no truth is opened and no additional fit is launched.
for attempt in $(seq 1 480); do
  if [[ -f "$OUT/FITS_AND_ENDPOINT_PARITY_COMPLETE.txt" ]]; then
    Rscript "$ANALYZER" "$OUT"
    exit 0
  fi
  if [[ -f "$OUT/FITS_INCOMPLETE.txt" ||
        -f "$OUT/endpoint_parity_v1/PARITY_FAILED.txt" ]]; then
    printf 'STAGE7F_ATTRIBUTION_STOPPED_AT_FIT_OR_PARITY_GATE\n' > \
      "$OUT/ATTRIBUTION_STOPPED.txt"
    exit 1
  fi
  sleep 30
done
printf 'STAGE7F_ATTRIBUTION_WAIT_TIMEOUT\n' > "$OUT/ATTRIBUTION_WAIT_TIMEOUT.txt"
exit 1
