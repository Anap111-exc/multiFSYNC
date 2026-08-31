# Stage 6A: G12 fitted-dimension overspecification

This directory contains the registered development-only Stage 6A screening
experiment. It uses one new baseline-strong data seed and three separate
fitted-dimension arms. The factor-count and FPCA-rank perturbations are not
combined. The single-seed result cannot establish cross-data generality.

Core files:

- `PROTOCOL_G12_DIMENSION_OVERSPEC_STAGE6A_V1_20260831.md`: scientific and
  execution contract;
- `DATA_MANIFEST.csv`, `FIT_CONFIGS.csv`, `FIT_MANIFEST.csv`: 1 data set,
  3 configurations, 36 fits and 3 within-configuration winners;
- `METRIC_MANIFEST.csv` and `MANIFEST_QC.csv`: registered outputs and manifest
  checks;
- `stage6a_common_20260831_v1.R` and `stage6a_runner_20260831_v1.R`: minimal
  generator, fit, capacity, truth-free selection and separately authorized
  evaluation runner;
- `prepare_server_20260831_v1.sh`, `launch_server_20260831_v1.sh`, and
  `status_server_20260831_v1.sh`: target-Linux commands.

The launch sequence generates and seals data, runs the three registered capacity
probes, and proceeds to all 36 fits only if the capacity gate passes. It freezes
three truth-free winners and then stops. It does not unseal truth, evaluate,
continue a fit, or create a formal result.

The server output root is outside Git by default:
`/root/g12-dimension-overspec-stage6a-20260831-v1`.
