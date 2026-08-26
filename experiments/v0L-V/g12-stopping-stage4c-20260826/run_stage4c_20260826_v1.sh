#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4C_ROOT:-/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1}"
cd "$experiment_root"

Rscript s4c_supervisor_20260826_v1.R
Rscript s4c_analysis_20260826_v1.R
