# multiFSYNC

`multiFSYNC` fits a Bayesian multi-study factor model for high-dimensional
continuous functional data.  It separates shared and study-specific factor
processes and permits time-varying coefficient functions for scalar baseline
covariates.

## Statistical conventions

- The CAVI objects named `zeta`, `xi`, `nu_phi`, and `nu_psi` are legacy names
  for the computational parameters `eta`, `chi`, `theta`, and `kappa`.
- Canonical FPCA functions, scores, eigenvalues, loading scales, factor order,
  and signs are produced by `orthonormalise_multi()`.
- Complete posterior second moments are used in FPCA post-processing, while the
  reported point estimate is capped at the fitted truncation rank `M`.
- Factor-level `omega` (`bool_var_spec_prob = FALSE`) is the thesis main model;
  variable-factor-level `omega` is an implemented sensitivity analysis.
- `lambda_orth = 0` is the default. Positive values are optional numerical
  regularisation, not hard identification constraints.

## ELBO and stopping

`compute_elbo_multi()` evaluates the complete ordinary-model ELBO only after a
full `T = 1` sweep. Its result contains the total, named components, direct
current expected-RSS cache, finite-value checks, and decline diagnostics. Both
factor-level and variable-factor-level `omega` hierarchies are included.

Three stopping rules are available:

- `"parameters"` monitors all main variational coordinates and covariance
  summaries. It is always defined and remains the public default for legacy
  script compatibility, but slowly moving Half-Cauchy auxiliary coordinates
  can make it unnecessarily conservative.
- `"elbo"` monitors the adjacent complete ordinary-ELBO increment. It is a
  valid stopping objective only when `lambda_orth = 0` and no adaptive
  Cholesky jitter is used during the `T = 1` phase.
- `"practical"` is recommended explicitly for production fits. Starting at
  the 80th `T = 1` sweep, it requires five consecutive passes of an ELBO
  safeguard and fitted-value, expected-RSS, and PPI stability checks over both
  20- and 60-sweep windows. A fit that has not passed by sweep 150 is returned
  as `slow_case = TRUE`, not reported as converged.

The practical defaults use an absolute ELBO rate threshold of `0.03` and
scientific-output thresholds of `1e-3`, `1e-3`, and `1e-2` for fitted NRMSE,
relative expected-RSS change, and maximum PPI change. Total-ELBO relative
change is retained as a diagnostic rather than a gate because positive and
negative components can cancel near zero. If `lambda_orth > 0` or a `T = 1`
solve requires jitter, an ELBO-based request is downgraded to parameter
stopping and the ordinary ELBO is diagnostic only; a requested practical fit
still respects its 150-sweep slow-case cap.

## Reproducible use

Observation times must be finite and pre-normalised to `[0,1]`.  If a custom
`time_g` is supplied it must be unique and strictly increasing.  With
study-specific mean functions, scalar covariates must retain full-rank
within-study variation; study-constant covariates are not identifiable from the
means.

`L_s` can be either the historical scalar (recycled to all studies) or a
length-`S` vector. For example, `L_s = c(0, 1, 2)` is paired with an `M_s`
list whose element lengths are `0`, `1`, and `2`. Fitted objects always expose
the explicit vector as `L_s_by_study`; `L_s` remains scalar when all study
counts are equal for backward compatibility.

```r
devtools::load_all(".")
dat <- simulate_multi_study_osullivan(
  S = 2, n_s = c(10, 10), p = 6, d = 1,
  L_f = 1, L_s = 1, M_f = 2,
  M_s = list(2, 2), K = 7, n_obs = 25,
  use_explicit = TRUE, mean_structure = "shared",
  specific_time_heterogeneity = TRUE,
  identified_loadings = TRUE, seed = 1
)
fit <- bayesSYNC_multi(
  dat$Y, Z = dat$Z, time_obs = dat$time_obs,
  L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
  K = 7, anneal = NULL, maxit = 150,
  convergence_rule = "practical",
  initialization = "random",
  bool_scale = TRUE, verbose = FALSE, seed = 101
)
```

Here `maxit = 150` is sufficient because annealing is disabled. With an
annealing ladder, `maxit` counts both annealing and `T = 1` sweeps and must
allow the planned annealing sweeps plus the practical limit of 200.
The primary workflow uses the public default `initialization = "random"`.
Residual FPCA remains an optional initialization diagnostic rather than the
main analysis path. It removes study means and covariate effects and initializes
only iteration-zero means; it uses no simulated truth and does not alter the
prior, CAVI updates, annealing objective, or ordinary ELBO. The practical rule
keeps strict 20-sweep fitted/RSS/PPI tolerances of
0.001/0.001/0.01, while its separately calibrated 60-sweep cumulative
tolerances are 0.003/0.006/0.01. The ELBO safeguard accepts either a total
absolute change rate no larger than 0.03 or a per-scalar-response rate no
larger than 1e-4, avoiding a sample-size-dependent stopping penalty.

For the current conservative thesis analysis, pre-specify twelve distinct
random seeds, assess after eight, and expand to twelve only if unresolved:

```r
fixed_seeds_12 <- c(
  101L, 202L, 303L, 404L, 505L,
  606L, 707L, 808L, 909L, 1010L, 1111L, 1212L
)
multi_fit <- bayesSYNC_multi_multistart(
  dat$Y, Z = dat$Z, time_obs = dat$time_obs,
  L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
  K = 7, anneal = NULL, maxit = 200,
  bool_scale = TRUE, verbose = FALSE,
  start_seeds = fixed_seeds_12,
  start_initializations = "random",
  stage_sizes = c(8L, 12L)
)
fit <- multi_fit$fit
```

The wrapper first requires a finite valid ordinary ELBO, no material ELBO
decline, and no `T = 1` jitter. It then selects the largest ELBO among those
valid endpoints. A selected endpoint that reaches the 200-sweep budget without
practical convergence is returned with `selected_endpoint_unfinished = TRUE`
and can never be production-ready. `initialization_stable = TRUE` additionally
requires at
least two starts in the selected ELBO basin, agreement of their fitted values,
expected RSS, and PPI, practical convergence of the selected endpoint, no
higher unreliable competitor, and the minimum exploration count. Always inspect
`production_ready`, `initialization_status`, `unresolved_reasons`,
`stage_history`, `higher_ineligible_competitor`, and
`higher_unreliable_competitor`. A successful finite but invalid fit with an
anomalously higher ELBO is never selected, but it blocks `production_ready`
until investigated.

`stage_sizes = NULL` runs every supplied seed. The calibrated production policy
assesses after eight starts and expands to twelve when unresolved. A fit
with `production_ready = FALSE` must remain explicitly diagnostic. Prefixes
with 1, 3, or 5 starts are diagnostic only. Staged fitting has a default
minimum-exploration gate of
`min(8, length(start_seeds))`, so `stage_sizes = c(3L, 5L, 12L)` cannot stop at
3 or 5 under the default. The frozen formal policy uses 8 -> 12 after
calibration; `conditional_initialization_stable` records what the older
two-support rule would have decided before this gate. Setting
`min_exploration_starts = 3L` reproduces that conditional behaviour for
diagnosis or compatibility, not for the formal analysis. The fixed seed order
and attempted count must be reported. `retry_only = TRUE` is also a diagnostic
fast path and must not be used as the production stability policy. No finite
multi-start policy is a global-optimum guarantee.

Basin support is counted by initialization independence ID, not raw replica
count. With `perturb_sd = 0`, repeated deterministic residual-FPCA starts count
only once. The objective-basin tolerance combines relative tolerance with a
0.05 absolute floor to avoid artificial splitting when positive and negative
ELBO components nearly cancel. Earlier ten-start calibration found only 21/40
practically converged starts and no independently resolved selected optimum.
That evidence motivates the 8 -> 12 exploration gate and the explicit
unfinished-endpoint status; it does not justify selecting with simulated truth.
The variable-factor omega extension remains a separately reported sensitivity
analysis.

Post-processed means are indexed as `fit$list_mu_hat[[study]][[variable]]`.
Original-scale coefficient functions are in `fit$list_beta_hat`; canonical
loadings are in `fit$mu_q_a_hat` and `fit$mu_q_b_specific_hat`; original-scale
noise variance estimates are in `fit$sigsq_eps_hat`. Raw unsuffixed variational
objects remain on the working scale for backward compatibility, with explicit
`*_original` copies where relevant.

## Tests

Run the standard package tests with `devtools::test()`. Lightweight manual
regression scripts are in `inst/manual-tests/`; run them from the package root
so that they resolve paths without depending on the former project directory.
The complete-ELBO regression coverage is in `test_elbo_components.R` and
`test_elbo_monotonicity.R`; residual-FPCA initialization, practical stopping,
and multi-start policy coverage is in `test_initialization_quality.R`,
`test_practical_stopping.R`, and `test_initialization_policy.R`.
The frozen small-scale stability experiment is under
`../../对比实验/实用停止规则与初始化稳定性_20260724`.
The component-wise local-basin diagnosis is under
`../../对比实验/初始化盆地诊断_20260724`.
The paired initialization-quality baseline and calibrated experiment are under
`../../对比实验/初始化质量修复_20260724`.
