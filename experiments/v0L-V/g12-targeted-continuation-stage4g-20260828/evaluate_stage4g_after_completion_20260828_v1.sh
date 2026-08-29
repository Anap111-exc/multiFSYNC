#!/usr/bin/env bash
set -u

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$root"

while [[ ! -f CONTINUATIONS_COMPLETE.txt && ! -f CONTINUATION_ERROR.txt ]]; do
  if ! pgrep -f '[R]script s4g_targeted_continuation_20260828_v1[.]R' \
      >/dev/null 2>&1; then
    printf '%s\n' 'status=ERROR_RUNNER_EXITED_WITHOUT_TERMINAL_MARKER' \
      > STAGE4G_EVALUATION_WRAPPER_ERROR.txt
    exit 2
  fi
  sleep 30
done

if [[ -f CONTINUATION_ERROR.txt ]]; then
  printf '%s\n' 'status=SKIPPED_CONTINUATION_ERROR_RETAINED' \
    > STAGE4G_EVALUATION_WRAPPER_ERROR.txt
  exit 2
fi

/usr/bin/time -v Rscript s4g_evaluate_continuations_20260828_v1.R \
  > STAGE4G_EVALUATION.log 2>&1
status=$?
if [[ $status -ne 0 ]]; then
  printf 'status=ERROR\nexit_status=%d\n' "$status" \
    > STAGE4G_EVALUATION_WRAPPER_ERROR.txt
  exit "$status"
fi
printf 'status=PASS\ntruth_read_after_continuations_complete=TRUE\n' \
  > STAGE4G_EVALUATION_WRAPPER_COMPLETE.txt
