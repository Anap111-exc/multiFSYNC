#!/usr/bin/env bash
set -euo pipefail

readonly STAGE7C=/root/multifsync-results/covariate_routeb_stage7c_20260927_v1
readonly OUT=/root/multifsync-results/stage7g_20260928_v1
readonly ANALYZER=/root/multifsync-stage7g-20260928-v1/analyze_stage7g_snapshots_20260928_v1.R
readonly ARCHIVE=/root/multifsync-results/STAGE7G_RESULTS_20260928_v1.tar.gz

for _ in $(seq 1 180); do
  if test -f "$OUT/FITS_AND_ENDPOINT_PARITY_COMPLETE.txt"; then
    break
  fi
  if test -f "$OUT/FITS_INCOMPLETE.txt" || \
     test -f "$OUT/endpoint_parity_v1/PARITY_FAILED.txt"; then
    printf 'STAGE7G_FINISH_ABORTED_FIT_OR_PARITY_FAILURE\n' > \
      "$OUT/FINISH_ABORTED.txt"
    exit 1
  fi
  sleep 60
done
test -f "$OUT/FITS_AND_ENDPOINT_PARITY_COMPLETE.txt"
test -f "$OUT/endpoint_parity_v1/PARITY_PASS.txt"
test ! -e "$OUT/analysis_v1"
test ! -e "$ARCHIVE"

Rscript "$ANALYZER" "$STAGE7C" "$OUT"
test -f "$OUT/analysis_v1/ANALYSIS_COMPLETE.txt"
tar -czf "$ARCHIVE" -C /root/multifsync-results stage7g_20260928_v1
sha256sum "$ARCHIVE" > "$OUT/STAGE7G_RESULTS_ARCHIVE_SHA256.txt"
printf 'STAGE7G_FULL_DIAGNOSTIC_COMPLETE\n' > "$OUT/FULL_DIAGNOSTIC_COMPLETE.txt"
