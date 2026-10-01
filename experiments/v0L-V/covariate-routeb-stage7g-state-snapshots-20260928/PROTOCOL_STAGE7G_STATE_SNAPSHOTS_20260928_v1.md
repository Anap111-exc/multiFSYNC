# Stage 7G: full state snapshot replay of three frozen C1 winners

## Question and boundary

Stage 7F established exact endpoint parity and found that seed 01/02 lose the
study 2 specific block during the early annealing phase, but its aggregate
trace cannot identify the first directional deviation or disentangle score,
loading, and time-function scale changes. Stage 7G replays the same three
truth-blind Stage 7C C1 winners only to record raw variational states.

No new candidate is selected. The replay is mechanism evidence, not part of
the formal Stage 7C candidate pool. No truth file is opened by the runner.

## Source and fit contract

- Base: Stage 7C source archive SHA256
  `93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924`.
- Stage 7G source archive SHA256
  `b3e9050f1de7154baff28b7dd986254df0b634bf85d6f50f464d00dbd4c8f17f`.
- The archives contain the same 84 files; only `R/multi_core.R` differs. Its
  new private option sink copies state at iteration 0 and after the last
  coordinate update of each requested sweep. It does not read or change any
  model parameter, prior, CAVI update, ELBO, or stopping rule.
- Inputs, seeds, model dimensions, pre-score setting, G12 control, annealing,
  499-sweep budget, and the 79 Stage 7F driver-trace sweeps are identical to
  Stage 7F. Package version remains `0.3.0.9000` in an isolated library.
- Three registered fits: seed 01 C1 start 06; seed 02 C1 start 11; seed 03 C1
  start 06. Observation, Stage 7C fit, Stage 7F fit, and winner hashes are
  checked before fitting.

## Recorded state

The sink saves **all 500 end-of-sweep states, including iteration 0**, rather
than a sparse early subset. Each state includes raw means and covariances for
mean, beta, shared and specific time functions, shared and specific scores,
loading slab means/covariances, effective loadings, support probabilities and
second moments, noise/regularization precision moments, omega parameters,
temperature, and annealing status. The full fit also retains its basis/grid,
ELBO, diagnostics, original-scale outputs, and the existing 79-sweep
within-sweep driver trace. Input observations remain bound by SHA256. A CSV
indexes every snapshot. If a fit errors, its partial snapshots and warnings
are retained; there is no automatic rerun.

## Gates and analysis

1. Parse and selftest the runner, compare both source archives, install the
   isolated package, and run a two-sweep smoke. The smoke must obtain states
   0, 1, 2 and match final raw score/support arrays exactly.
2. Run only the three registered fits. Stop on any fit error, warning needing
   interpretation, missing state, nonfinite ELBO, or run-contract failure.
3. **Before opening truth**, compare all 400 ordinary ELBO values and 18 raw
   fields with frozen Stage 7C, the 9 old Stage 7C trace checkpoints, all 79
   Stage 7F driver-trace sweeps, and all available final snapshot fields with
   the fitted object. Require exact equality for snapshot and dense trace,
   and at most `1e-6` for the Stage 7C endpoint comparison. Stop mechanism
   interpretation if this gate fails.
4. Only after a parity pass, compute within-fit loading-direction changes,
   score second moments, fitted-basis function L2 energy, support, and block
   signal contribution across sweeps. Truth-assisted angles can then be
   calculated as a separate evaluation, never as a fitting/selection input.

The key discrimination is whether study 2 loadings rotate away from the true
direction before their norm and score second moment collapse, or whether a
correct direction first loses magnitude while the score-loading feedback
shrinks. Estimate the first departure window, report cross-seed contrasts,
and keep unsupported causal explanations conditional.

## Resource and rollback boundary

Three fits can run in parallel at one thread each on the already-started ECS.
Each raw snapshot was estimated at about 0.6 MiB, so retaining 500 per fit
adds roughly 0.3 GiB of in-memory records per worker before serialization.
The isolated source package, runner, and new output directory can be ignored
without changing the project package or historical evidence. No file is
deleted, reset, or overwritten.
