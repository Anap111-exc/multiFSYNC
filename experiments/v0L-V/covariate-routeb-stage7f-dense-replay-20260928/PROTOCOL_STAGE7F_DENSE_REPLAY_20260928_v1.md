# Stage 7F: frozen C1 winner dense replay

Authorization: the user explicitly authorized three dense diagnostic replays on 2026-09-28. This experiment is separate from the Stage 7C formal candidate pool. No new start, model update, prior, CAVI formula, ELBO formula, or truth-guided fit choice is authorized here.

## Inputs and binding

- Stage 7C result root: `/root/multifsync-results/covariate_routeb_stage7c_20260927_v1`.
- Stage 7C source archive: `/root/multiFSYNC-stage7c-20260927-v1.tar.gz`, SHA256 `93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924`.
- Isolated installed package: `/root/multiFSYNC-stage7c-lib-20260927-v1/multiFSYNC`, version `0.3.0.9000`.
- Frozen truth-free winner CSV SHA256: `1e1e7791b3e975287ab195c87a801ab0fa24f97d26558a85fccc3b11af70a877`.
- Replays: `covrb7c_01__C1__start06` (seed `99716106`), `covrb7c_02__C1__start11` (seed `99716211`), `covrb7c_03__C1__start06` (seed `99716306`). The runner also verifies individual observation and reference fit hashes before creating an output directory.
- New result root: `/root/multifsync-results/stage7f_20260928_v1`.

## Single intervention

Copy the Stage 7C main `bayesSYNC_multi_pre_score` call: `random`, one pre-score sweep, Gram unit-energy function initialization, `L_f=1`, `L_s=(1,1)`, `M_f=2`, `M_s=(2,2)`, `K=9`, `n_g=51`, `anneal=(1,1.9,100)`, `maxit=499`, `n_cpus=1`, `bool_scale=FALSE`, `bool_var_spec_prob=FALSE`, `d_0=500`, `lambda_orth=0`, practical G12 stopping with checkpoints 380 and 400, and the original seeds. Only change `trace_sweeps` from `(1,380,400)` to `sort(unique(c(1:20, seq(25,90,5), 91:110, seq(115,200,5), 250,300,350,380,400,450,499)))`. The dense schedule has 79 total-sweep indices. The update formulas and objective are unchanged.

## Ordered gates

1. Run local `selftest` for the registered IDs and trace schedule.
2. On the ECS, run `prepare`: check source, winner, observation, and reference fit hashes, verify the three C1 winners, then create the isolated result root and manifest.
3. Run a two-sweep micro smoke on the third observation to verify the diagnostic trace records both sweeps. The smoke is a separate short fit, with no annealing and no formal-candidate status.
4. Run exactly the three full replays, one CPU per fit, concurrently. Retain each `fit.rds`, driver trace, ordinary ELBO path, warnings, logs, and completion/error marker. No automatic extra start or rerun.
5. Compare each replay to its Stage 7C reference before mechanism interpretation. Require matching configuration and 99+400 sweep counts, max absolute difference at most `1e-6` in the full ordinary ELBO path and raw endpoint loading PPI, loadings, slab moments, scores, function coefficients, selected covariance blocks, and noise precision. At the three originally recorded iterations, require matching trace row keys and numeric metrics within `1e-6`. The full record is in `endpoint_parity_v1`.
6. If any check fails, write `PARITY_FAILED.txt`, stop interpretation, and audit the first discrepancy. If all pass, write `PARITY_PASS.txt`, combine traces, and analyze the timing of B2 loss against B1 and the third-data control. The old truth labels may be used only after the parity gate for scientific interpretation, never for fit choice or stopping.

The independent post-processing helper waits for the fit and parity completion marker, then runs the truth-blind trace table script once. It exits on a fit or parity failure and never starts another fit.

## Outputs and decision

The attribution table will show, by total sweep, annealing phase/temperature, block/study/component, loading PPI and posterior mean norm, score norm, time-function coefficient norm, data driver, and driver/precision ratio. A separate script runs only after the parity marker, compares PPI crossings at `1e-2`, `1e-4`, and `1e-6`, records adjacent observed changes, and marks *stable joint closure* when mean loading PPI stays below `1e-4` and posterior mean loading norm stays below `1e-3` through the terminal record. These thresholds locate a recorded interval; they do not establish causality or replace the full trajectory. The trace has coefficient norms rather than actual function L2 norms and has no per-sweep loading direction. If these limits prevent a clear mechanism decision, the next task is a separate diagnostic snapshot design, not an algorithm change.

Estimated full replay time is 75–100 minutes wall time for three parallel fits, based on the Stage 7C reference runtimes of 4309–4526 seconds each, plus setup, smoke, and analysis. Stop on failed binding, smoke, fit completion, or endpoint parity. Rollback boundary is this new Stage 7F directory only; existing Stage 7C/7E results and tracked model code remain untouched.
