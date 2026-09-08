#!/usr/bin/env bash
set -euo pipefail

stage6d_output="${STAGE6D_OUTPUT_ROOT:-/root/g12-dimension-selection-stage6d-20260906-v1}"
complete=0
errors=0
if [[ -d "$stage6d_output/evaluation" ]]; then
  complete=$(find "$stage6d_output/evaluation" -name EVALUATION_COMPLETE.txt \
    -type f | wc -l)
  errors=$(find "$stage6d_output/evaluation" -name EVALUATION_ERROR.txt \
    -type f | wc -l)
fi
echo "evaluation_fit_complete=$complete"
echo "evaluation_errors=$errors"
if [[ -f "$stage6d_output/evaluation/EVALUATION_COMPLETE.txt" ]]; then
  echo "aggregate_evaluation_complete=TRUE"
else
  echo "aggregate_evaluation_complete=FALSE"
fi
pgrep -af '[e]valuate_stage6d_20260908_v1.R' || true
if [[ -s "$stage6d_output/STAGE6D_EVALUATION.stderr.log" ]]; then
  echo "=== stderr tail ==="
  tail -n 20 "$stage6d_output/STAGE6D_EVALUATION.stderr.log"
fi
