test_that("practical control defaults and validation are explicit", {
  control <- multiFSYNC:::.validate_practical_control()
  expect_equal(
    control[c(
      "min_t1", "window", "long_window", "consecutive", "max_t1",
      "elbo_abs_rate", "elbo_per_response_rate",
      "fitted_nrmse", "rss_rel", "ppi_max_abs",
      "long_fitted_nrmse", "long_rss_rel", "long_ppi_max_abs"
    )],
    list(
      min_t1 = 80L, window = 20L, long_window = 60L,
      consecutive = 5L, max_t1 = 200L, elbo_abs_rate = 0.03,
      elbo_per_response_rate = 1e-4,
      fitted_nrmse = 1e-3,
      rss_rel = 1e-3, ppi_max_abs = 1e-2,
      long_fitted_nrmse = 3e-3,
      long_rss_rel = 6e-3, long_ppi_max_abs = 1e-2
    )
  )
  legacy <- multiFSYNC:::.validate_practical_control(list(
    fitted_nrmse = 2, rss_rel = 3, ppi_max_abs = 4
  ))
  expect_identical(legacy$long_fitted_nrmse, 2)
  expect_identical(legacy$long_rss_rel, 3)
  expect_identical(legacy$long_ppi_max_abs, 4)
  expect_error(
    multiFSYNC:::.validate_practical_control(
      list(window = 20L, min_t1 = 20L)),
    "long_window"
  )
  expect_error(
    multiFSYNC:::.validate_practical_control(list(unknown = 1)),
    "Unknown"
  )
})

test_that("practical windows reject objective drift and scientific-output drift", {
  control <- multiFSYNC:::.validate_practical_control(list(
    min_t1 = 21L, window = 10L, long_window = 20L,
    consecutive = 1L, max_t1 = 30L
  ))
  snapshot <- list(
    fitted = c(1, 2), rss = c(3, 4), ppi = numeric())

  stable_elbo <- -100 + (0:20) * 0.001
  stable_diag <- multiFSYNC:::.elbo_history_diagnostics(stable_elbo)
  stable <- multiFSYNC:::.practical_window_diagnostic(
    21L, stable_elbo, stable_diag,
    snapshot, snapshot, snapshot, control)
  expect_true(stable$objective_pass)
  expect_true(stable$output_pass)
  expect_true(stable$pass)

  large_scale_drift <- -1e6 + 0:20
  drift_diag <- multiFSYNC:::.elbo_history_diagnostics(large_scale_drift)
  drift <- multiFSYNC:::.practical_window_diagnostic(
    21L, large_scale_drift, drift_diag,
    snapshot, snapshot, snapshot, control)
  expect_lt(drift$elbo_rel_rate, 5e-5)
  expect_gt(drift$elbo_abs_rate, control$elbo_abs_rate)
  expect_false(drift$objective_pass)
  expect_false(drift$pass)

  scaled_drift <- multiFSYNC:::.practical_window_diagnostic(
    21L, large_scale_drift, drift_diag,
    snapshot, snapshot, snapshot, control,
    objective_size = 20000)
  expect_equal(scaled_drift$elbo_abs_rate, drift$elbo_abs_rate)
  expect_equal(scaled_drift$elbo_per_response_rate,
               drift$elbo_abs_rate / 20000)
  expect_true(scaled_drift$objective_pass)

  expect_error(
    multiFSYNC:::.practical_window_diagnostic(
      21L, stable_elbo, stable_diag,
      snapshot, snapshot, snapshot, control,
      objective_size = 0),
    "objective_size"
  )

  moving_snapshot <- snapshot
  moving_snapshot$fitted <- moving_snapshot$fitted + 0.02
  moving <- multiFSYNC:::.practical_window_diagnostic(
    21L, stable_elbo, stable_diag,
    moving_snapshot, snapshot, snapshot, control)
  expect_true(moving$objective_pass)
  expect_false(moving$short_output_pass)
  expect_false(moving$output_pass)
  expect_false(moving$pass)

  long_moving_reference <- snapshot
  long_moving_reference$rss <- long_moving_reference$rss + 0.1
  cumulative_drift <- multiFSYNC:::.practical_window_diagnostic(
    21L, stable_elbo, stable_diag,
    snapshot, snapshot, long_moving_reference, control)
  expect_true(cumulative_drift$short_output_pass)
  expect_false(cumulative_drift$long_output_pass)
  expect_false(cumulative_drift$pass)

  decreasing_elbo <- stable_elbo
  decreasing_elbo[11L] <- decreasing_elbo[10L] - 1
  decreasing_diag <- multiFSYNC:::.elbo_history_diagnostics(decreasing_elbo)
  decreasing <- multiFSYNC:::.practical_window_diagnostic(
    21L, decreasing_elbo, decreasing_diag,
    snapshot, snapshot, snapshot, control)
  expect_true(decreasing$monotone)
  expect_false(decreasing$long_monotone)
  expect_false(decreasing$pass)
})

test_that("empty PPI and non-finite or mismatched outputs are handled safely", {
  control <- multiFSYNC:::.validate_practical_control()
  empty_ppi <- list(fitted = c(1, 2), rss = c(3, 4), ppi = numeric())
  comparison <- multiFSYNC:::.practical_compare(
    empty_ppi, empty_ppi, control)
  expect_equal(comparison$ppi_max_abs, 0)
  expect_equal(comparison$ppi_quantile_abs, 0)
  expect_equal(comparison$factor_ppi_max_abs, 0)
  expect_true(comparison$pass)

  bad <- empty_ppi
  bad$fitted <- c(1, Inf)
  expect_false(
    multiFSYNC:::.practical_compare(bad, empty_ppi, control)$pass)
  short <- empty_ppi
  short$rss <- 3
  expect_false(
    multiFSYNC:::.practical_compare(short, empty_ppi, control)$pass)
})

test_that("quantile-factor PPI gate ignores isolated loading flips but not factor drift", {
  reference <- list(
    fitted = c(1, 2), rss = c(3, 4),
    ppi = numeric(100), factor_ppi = c(1, 0)
  )
  isolated_flip <- reference
  isolated_flip$ppi[1L] <- 1

  max_control <- multiFSYNC:::.validate_practical_control(list(
    ppi_gate = "max"
  ))
  max_comparison <- multiFSYNC:::.practical_compare(
    isolated_flip, reference, max_control
  )
  expect_equal(max_comparison$ppi_max_abs, 1)
  expect_equal(max_comparison$ppi_quantile_abs, 0)
  expect_false(max_comparison$ppi_pass)

  quantile_control <- multiFSYNC:::.validate_practical_control(list(
    ppi_gate = "quantile_factor",
    ppi_quantile = 0.95
  ))
  quantile_comparison <- multiFSYNC:::.practical_compare(
    isolated_flip, reference, quantile_control
  )
  expect_true(quantile_comparison$ppi_pass)
  expect_true(quantile_comparison$pass)

  factor_drift <- isolated_flip
  factor_drift$factor_ppi <- c(1, 1)
  factor_comparison <- multiFSYNC:::.practical_compare(
    factor_drift, reference, quantile_control
  )
  expect_false(factor_comparison$ppi_pass)
  expect_false(factor_comparison$pass)

  expect_error(
    multiFSYNC:::.validate_practical_control(list(ppi_gate = "bad")),
    "ppi_gate"
  )
  expect_error(
    multiFSYNC:::.validate_practical_control(list(ppi_quantile = 0)),
    "ppi_quantile"
  )
})

test_that("practical stopping counts only T=1 and distinguishes slow cases", {
  dat <- make_practical_test_data()
  common <- list(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, n_g = 21L, bool_scale = FALSE,
    verbose = FALSE, seed = 8102L,
    convergence_rule = "practical"
  )

  fit <- do.call(bayesSYNC_multi, c(common, list(
    anneal = c(1, 1.25, 2), maxit = 30L,
    practical_control = small_practical_control(FALSE)
  )))
  expect_true(fit$converged)
  expect_true(fit$practical_converged)
  expect_identical(fit$convergence_status, "converged_practical")
  expect_equal(fit$t1_sweeps, length(fit$ELBO))
  expect_equal(fit$annealing_sweeps, fit$i_iter - fit$t1_sweeps)
  expect_equal(fit$annealing_sweeps, 1L)
  expect_gte(fit$t1_sweeps, 4L)
  expect_equal(length(fit$practical_checkpoints), 2L)

  expect_warning(
    slow_fit <- do.call(bayesSYNC_multi, c(common, list(
      anneal = NULL, maxit = 20L,
      practical_control = small_practical_control(TRUE)
    ))),
    "slow case"
  )
  expect_false(slow_fit$converged)
  expect_false(slow_fit$practical_converged)
  expect_true(slow_fit$slow_case)
  expect_identical(slow_fit$convergence_status, "slow_case")
  expect_equal(slow_fit$t1_sweeps, 5L)
})

test_that("practical stopping validates its full T=1 budget and checkpoint horizon", {
  dat <- make_practical_test_data(8151L)
  common <- list(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, anneal = NULL, n_g = 21L,
    convergence_rule = "practical",
    bool_scale = FALSE, verbose = FALSE, seed = 8152L
  )

  expect_error(
    do.call(bayesSYNC_multi, c(common, list(
      maxit = 9L,
      practical_control = small_practical_control(FALSE)
    ))),
    "planned annealing sweeps"
  )

  checkpoint_control <- small_practical_control(FALSE)
  checkpoint_control$checkpoints <- 11L
  expect_error(
    do.call(bayesSYNC_multi, c(common, list(
      maxit = 10L,
      practical_control = checkpoint_control
    ))),
    "separate forced reference fit"
  )
})

test_that("positive lambda downgrades practical stopping", {
  dat <- make_practical_test_data(8201L)
  expect_warning(
    fit <- bayesSYNC_multi(
      Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
      L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      K = 3L, anneal = NULL, n_g = 21L, maxit = 5L,
      tol_abs = 1e100, tol_rel = 1e100,
      convergence_rule = "practical",
      practical_control = small_practical_control(FALSE),
      lambda_orth = 0.01,
      bool_scale = FALSE, verbose = FALSE, seed = 8202L
    ),
    "downgrading"
  )
  expect_identical(fit$convergence_rule_requested, "practical")
  expect_identical(fit$convergence_rule, "parameters")
  expect_identical(fit$elbo_role, "diagnostic_only")
  expect_false(fit$elbo_objective_valid)
})

test_that("a downgraded requested practical fit keeps its T=1 budget", {
  dat <- make_practical_test_data(8221L)
  control <- small_practical_control(FALSE)
  control$max_t1 <- 4L
  control$checkpoints <- c(3L, 4L)
  warnings_seen <- character()
  fit <- withCallingHandlers(
    bayesSYNC_multi(
      Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
      L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      K = 3L, anneal = NULL, n_g = 21L, maxit = 20L,
      tol_abs = 0, tol_rel = 0,
      convergence_rule = "practical",
      practical_control = control,
      lambda_orth = 0.01,
      bool_scale = FALSE, verbose = FALSE, seed = 8222L
    ),
    warning = function(condition) {
      warnings_seen <<- c(warnings_seen, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("downgrading", warnings_seen, fixed = TRUE)))
  expect_true(any(grepl(
    "eligible result within 4 T = 1 sweeps", warnings_seen, fixed = TRUE
  )))
  expect_false(fit$converged)
  expect_true(fit$slow_case)
  expect_identical(fit$convergence_rule, "parameters")
  expect_identical(fit$convergence_status, "slow_case")
  expect_equal(fit$t1_sweeps, 4L)
})

test_that("T=1 jitter downgrades practical stopping without a false practical stop", {
  dat <- make_practical_test_data(8251L)
  original_inverse_spd <- multiFSYNC:::.inverse_spd
  injected_jitter <- FALSE

  testthat::local_mocked_bindings(
    .inverse_spd = function(
      precision, context = "unspecified",
      jitter_relative = c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4)
    ) {
      inverse <- original_inverse_spd(
        precision, context = context,
        jitter_relative = jitter_relative
      )
      if (!injected_jitter &&
          !identical(context, "linear_prior_covariance")) {
        multiFSYNC:::.record_spd_diagnostic(
          context = paste0(context, "::test-injected-practical-jitter"),
          jitter = 1e-12,
          relative_jitter = 1e-12,
          attempts = 2L
        )
        injected_jitter <<- TRUE
      }
      inverse
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    fit <- bayesSYNC_multi(
      Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
      L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      K = 3L, anneal = NULL, n_g = 21L, maxit = 10L,
      tol_abs = 1e100, tol_rel = 1e100,
      convergence_rule = "practical",
      practical_control = small_practical_control(FALSE),
      bool_scale = FALSE, verbose = FALSE, seed = 8252L
    ),
    regexp = "Adaptive Cholesky jitter.*downgrading"
  )

  expect_true(injected_jitter)
  expect_true(fit$converged)
  expect_identical(fit$convergence_status, "converged_parameters")
  expect_identical(fit$convergence_rule_requested, "practical")
  expect_identical(fit$convergence_rule, "parameters")
  expect_identical(fit$elbo_role, "diagnostic_only")
  expect_false(fit$elbo_objective_valid)
  expect_identical(
    fit$elbo_invalid_reasons,
    "adaptive_cholesky_jitter_at_T1"
  )
  expect_identical(fit$elbo_t1_jitter_count, 1L)
  expect_identical(fit$linear_solver_diagnostics$jitter_count, 1L)
  expect_true(any(grepl(
    "test-injected-practical-jitter",
    fit$linear_solver_diagnostics$events$context,
    fixed = TRUE
  )))

  expect_false(fit$practical_converged)
  expect_false(fit$objective_converged)
  expect_equal(nrow(fit$practical_diagnostics), 0L)
  expect_null(fit$practical_final)
})

test_that("sparse reference checkpoints build fitted snapshots only when needed", {
  dat <- make_practical_test_data(8271L)
  original_snapshot <- multiFSYNC:::.practical_snapshot
  snapshot_calls <- 0L
  testthat::local_mocked_bindings(
    .practical_snapshot = function(...) {
      snapshot_calls <<- snapshot_calls + 1L
      original_snapshot(...)
    },
    .package = "multiFSYNC"
  )

  fit <- suppressWarnings(bayesSYNC_multi(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, anneal = NULL, n_g = 21L, maxit = 4L,
    tol_abs = 0, tol_rel = 0,
    convergence_rule = "parameters",
    practical_control = list(checkpoints = c(2L, 4L)),
    bool_scale = FALSE, verbose = FALSE, seed = 8272L
  ))

  expect_equal(snapshot_calls, 2L)
  expect_identical(names(fit$practical_checkpoints), c("2", "4"))
})

test_that("RSS fitted cache reuses the posterior mean without changing RSS", {
  dat <- make_practical_test_data(8301L)
  fit <- suppressWarnings(bayesSYNC_multi(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, anneal = NULL, n_g = 21L, maxit = 4L,
    tol_abs = 0, tol_rel = 0,
    convergence_rule = "parameters",
    bool_scale = FALSE, verbose = FALSE, seed = 8302L
  ))
  args <- list(
    Y = dat$Y, C = dat$C, list_cp_C = fit$list_cp_C,
    mu_q_nu_mu = fit$mu_q_nu_mu,
    Sigma_q_nu_mu = fit$Sigma_q_nu_mu,
    mu_q_nu_beta = fit$mu_q_nu_beta,
    Sigma_q_nu_beta = fit$Sigma_q_nu_beta, Z = dat$Z,
    mu_q_zeta = fit$mu_q_zeta, Sigma_q_zeta = fit$Sigma_q_zeta,
    mu_q_nu_phi = fit$mu_q_nu_phi,
    Sigma_q_nu_phi = fit$Sigma_q_nu_phi,
    mu_q_xi = fit$mu_q_xi, Sigma_q_xi = fit$Sigma_q_xi,
    mu_q_nu_psi = fit$mu_q_nu_psi,
    Sigma_q_nu_psi = fit$Sigma_q_nu_psi,
    mu_q_a = fit$mu_q_a, term_a = fit$term_a,
    mu_q_b_specific = fit$mu_q_b_specific,
    term_b_specific = fit$term_b_specific,
    S = fit$S, n_s = fit$n_s, p = fit$p,
    L_f = fit$L_f, L_s = fit$L_s
  )
  plain <- do.call(
    multiFSYNC:::compute_rss_cache,
    c(args, list(return_fitted = FALSE)))
  cached <- do.call(
    multiFSYNC:::compute_rss_cache,
    c(args, list(return_fitted = TRUE)))
  expect_equal(cached$expected_rss, plain$expected_rss, tolerance = 0)
  expect_equal(
    cached$expected_rss_sum, plain$expected_rss_sum, tolerance = 0)
  expect_null(plain$fitted_values)

  manual <- numeric()
  for (i in seq_len(fit$n_s[1L])) {
    for (j in seq_len(fit$p)) {
      manual <- c(
        manual,
        as.vector(dat$C[[1L]][[i]] %*% fit$mu_q_nu_mu[[1L]][[j]])
      )
    }
  }
  expect_equal(cached$fitted_values, manual, tolerance = 1e-12)
})
