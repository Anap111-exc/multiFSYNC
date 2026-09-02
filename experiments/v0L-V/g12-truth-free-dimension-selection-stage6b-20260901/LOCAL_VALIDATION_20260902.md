# Stage 6B predictive-pilot local validation

Date: 2026-09-02

Environment: Windows, R 4.5.1. No full-size Stage 6B fit was started.

## Results

- all R files in the Stage 6B directory parsed successfully;
- predictive-pilot registration/selection tests passed `7/7`;
- the existing truth-free holdout-interface tests passed `8/8` without a
  full-size observation bundle;
- the registered runner environment check passed with 36 conditional manifest
  rows, 12 Phase-1 rows, at most 20 executed fits and one CPU per fit;
- the existing nonformal `S=2`, `n_s=(4,4)`, `p=20`, six-sweep
  split-fit-predict-score smoke passed;
- smoke all-factor held-out RMSE/NRMSE were `0.4955000/0.8167709`;
- smoke truth access was false and no full-size fit was started.

The micro smoke was run from an ASCII-only temporary source copy because the
current Windows `pkgload/desc` stack mis-decoded the Chinese project path when
resolving `DESCRIPTION`. The source copy was byte-for-byte copied from the
working tree; this is a Windows loading-path workaround, not a separate model
source. The target Linux check and full observation-only split contract remain
required before launch.

The multiFSYNC package source under `R/`, `DESCRIPTION` and `NAMESPACE` was not
modified in this Stage 6B change, so the previously validated installed package
is reused on the server.
