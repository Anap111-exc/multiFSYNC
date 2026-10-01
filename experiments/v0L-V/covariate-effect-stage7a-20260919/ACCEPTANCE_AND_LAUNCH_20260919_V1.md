# Stage 7A preflight and launch record (2026-09-19)

- Source-package-only change: hoist the study/subject mean cross-product
  cache outside the `j=1:p` loop in `R/update_beta.R`; no CAVI equation or
  ELBO change.
- Pre-change tiny `d>0` fit baseline saved locally in ignored
  `acceptance/beta_baseline_20260919_v1.rds`. Post-change comparison passed
  for beta posterior means, beta posterior covariances, mean block, ordinary
  ELBO, ordinary ELBO components, and sweep counts.
- New direct-matrix beta-update test: 12/12 expectations passed.
- Full Windows `testthat::test_local`: passed; 15 existing
  maximum-iteration-style warnings were reported, no failures.
- Windows two-arm micro smoke: 2/2 fits completed, 2/2 truth-free selections.
- Independent Alibaba Cloud Linux library installation: passed. Linux two-arm
  micro smoke: 2/2 fits completed, 2/2 truth-free selections.
- Uploaded source archive SHA256:
  `e85b6113d258464394cb898fa0237ea9d64511032d272711dd9284ff146ab00b`;
  the ECS copy matched byte-for-byte. A local immutable copy and the tiny
  pre-change baseline are retained under
  `对比实验/Stage7A协变量效应_20260919/`, outside the nested Git repository.
- ECS full-size preflight (`p=500,d=1`): 3 ordinary sweeps took 12.989 s,
  8 took 32.673 s; first 3 ordinary ELBO points were numerically identical.
  Estimated marginal sweep time: 3.9368 s. The `/usr/bin/time -v` preflight
  process peaked at 162532 KiB RSS. Serial-equivalent projection gives about
  9.82 hours for 144 fits at 8-way task concurrency before contention and
  setup overhead; it is not a guaranteed completion time.
- ECS data generation and pairing QC: 6 logical datasets, 12 sealed arm
  bundles, 144 registered fits. Truth unsealing and scientific evaluation
  were not authorized or started.
- Independent source: `/root/multiFSYNC-stage7a-20260919`.
  Independent installed library: `/root/multiFSYNC-stage7a-lib-20260919`.
  Full result root: `/root/multifsync-results/stage7a_covariate_20260919_v1`.
- Batch launcher PID 2600 observed running on 2026-09-19 at 20:31 CST, with
  eight independent fit tasks active, no terminal errors then. The launcher
  uses the pre-existing instance-role self-stop after the batch finishes and
  writes `SUPERVISOR_FINISHED.txt` first.

The batch is **running**, not completed. The truth-free selection has not yet
been frozen for the 144 main fits. No scientific results are claimed here.
The source changes and experiment scripts remain local/uncommitted; no GitHub
push or historical frozen-code replacement was performed.
