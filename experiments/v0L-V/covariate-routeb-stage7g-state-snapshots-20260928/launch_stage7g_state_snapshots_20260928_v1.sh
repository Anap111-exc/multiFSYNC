#!/usr/bin/env bash
set -euo pipefail

export R_LIBS=/root/multiFSYNC-stage7g-lib-20260928-v1
readonly RUNNER=/root/multifsync-stage7g-20260928-v1/run_stage7g_state_snapshots_20260928_v1.R
readonly STAGE7C=/root/multifsync-results/covariate_routeb_stage7c_20260927_v1
readonly OUT=/root/multifsync-results/stage7g_20260928_v1

Rscript "$RUNNER" supervisor "$STAGE7C" "$OUT"
Rscript "$RUNNER" verify "$STAGE7C" "$OUT"
test -f "$OUT/endpoint_parity_v1/PARITY_PASS.txt"
printf 'STAGE7G_FITS_AND_ENDPOINT_PARITY_COMPLETE\n' > \
  "$OUT/FITS_AND_ENDPOINT_PARITY_COMPLETE.txt"
