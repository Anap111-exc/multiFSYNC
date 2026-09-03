# Stage 6B truth-free dimension-selection development

Date: 2026-09-01

Status: Stage 6B-A complete; Stage 6B-B holdout interface implemented and
software-tested. The first full-size predictive-pilot execution stopped after
all 12 factor-screen fits because eight zero-warning fits failed during result
persistence. A minimal persistence bugfix was prepared on 2026-09-03 for a
clean rerun in a new output directory. Truth access remains unauthorized.

## Contents

- `PROTOCOL_G12_TRUTH_FREE_DIMENSION_SELECTION_STAGE6B_V1_DRAFT_20260901.md`:
  development protocol and evidence boundary.
- `PROTOCOL_G12_STAGE6B_PREDICTIVE_PILOT_V1_20260902.md`: registered sequential
  factor-count/FPCA-cap prediction pilot.
- `STAGE6B_PREDICTIVE_PILOT_V2_BUGFIX_RERUN_20260903.md`: failure diagnosis and
  unchanged-design clean-rerun boundary.
- `FIT_CONFIGS_STAGE6B_PILOT_V1.csv` and
  `START_SEEDS_STAGE6B_PILOT_V1.csv`: fixed candidate grid and four paired
  start seeds.
- `stage6b_predictive_pilot_tools_20260902_v1.R`: conditional manifest and
  stratified one-standard-error selection contracts.
- `run_stage6b_predictive_pilot_20260902_v1.R`: prepare, fit, sequential
  truth-free selection and summary entry point.
- `test_stage6b_predictive_pilot_20260902_v1.R`: targeted registration,
  selection and truth-path guards.
- `prepare_server_20260902_v1.sh`, `launch_server_20260902_v1.sh` and
  `status_server_20260902_v1.sh`: minimal target-Linux execution scripts.
- `stage6b_truth_free_tools_20260901_v1.R`: truth-free Stage 6A audit plus the
  experiment-local holdout splitter, posterior-mean predictor and scorer.
- `run_stage6ba_readonly_audit_20260901_v1.R`: read-only audit entry point.
- `summarize_stage6ba_readonly_audit_20260901_v1.R`: descriptive aggregation
  of audit CSVs only.
- `test_stage6b_truth_free_tools_20260901_v1.R`: targeted contract tests.
- `run_stage6bb_holdout_micro_smoke_20260901_v1.R`: nonformal software smoke.
- `stage6ba_readonly_audit_20260901_v1/`: versioned machine-readable read-only
  audit evidence.
- `validation_windows_20260901_v1/`: small successful targeted-test and smoke
  evidence; no fit object is retained.
- `validation_linux_20260901_v1/`: target-server micro-smoke evidence; no fit
  object is retained.
- `REPORT_STAGE6B_A_READONLY_AUDIT_AND_6B_B_INTERFACE_20260901.md`: results and
  interpretation.

## Interface boundary

`s6b_make_middle_time_holdout()` removes one interior whole-time vector from
every subject and rejects bundles containing truth-named fields. Version 1 is
restricted to the registered no-covariate route (`d=0`).

`s6b_predict_subject_time()` reconstructs a fitted subject at an arbitrary time
by linearly interpolating the fitted mean and eigenfunctions and combining them
with the fitted shared and study-specific scores/loadings. The primary mode uses
all fitted factors; `ppi_0.5` is a sensitivity mode.

`s6b_score_holdout()` averages loss first within each subject-time vector and
then over subjects. RMSE, NRMSE and MAE are usable point-prediction diagnostics.
The reported NLPD uses residual noise only; it is not a calibrated full
posterior predictive density and has
`primary_dimension_score_ready=FALSE`.

## Execution boundary

The completed audit read existing Stage 6A fit endpoints only. The micro smoke
used `S=2`, `n_s=(4,4)`, `p=20`, one CPU and six short sweeps and passed on both
Windows and the target Linux server. Neither operation generated a new
full-size data set, started a full-size Stage 6B fit, accessed truth, continued
a fit, or changed the multiFSYNC package source.

The 2026-09-02 pilot authorization covers only the registered observation-only
split, at most 20 G12 fits, and truth-free selection. It does not authorize new
data generation, truth evaluation, continuation, a full-data refit, or any
paper-level claim.
