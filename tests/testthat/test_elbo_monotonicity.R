# Complete-ELBO integration tests on deterministic, dependency-free toy data.
#
# The dimensions are deliberately small: these tests exercise CAVI bookkeeping
# and stopping semantics, not statistical recovery.

.complete_elbo_component_names <- c(
  "data_likelihood",
  "osullivan_mu", "osullivan_beta",
  "osullivan_phi", "osullivan_psi",
  "score_zeta", "score_xi",
  "spike_slab_shared", "spike_slab_specific",
  "omega_shared", "omega_specific",
  "half_cauchy_eps", "half_cauchy_mu",
  "half_cauchy_beta", "half_cauchy_phi",
  "half_cauchy_psi"
)

.make_tiny_elbo_data <- function(S = 1L, n = 3L, p = 2L,
                                 n_obs = 6L, with_covariate = FALSE,
                                 n_covariates = if (with_covariate) 1L else 0L,
                                 seed = 2101L) {
  set.seed(seed)

  time_obs <- lapply(seq_len(S), function(s) {
    lapply(seq_len(n), function(i) {
      seq(0.05, 0.95, length.out = n_obs)
    })
  })

  Z <- if (n_covariates > 0L) {
    lapply(seq_len(S), function(s) {
      x <- seq(-1, 1, length.out = n)
      vapply(seq_len(n_covariates), function(r) x^r, numeric(n))
    })
  } else {
    NULL
  }

  Y <- lapply(seq_len(S), function(s) {
    lapply(seq_len(n), function(i) {
      lapply(seq_len(p), function(j) {
        tt <- time_obs[[s]][[i]]
        z_effect <- if (n_covariates > 0L) {
          sum(Z[[s]][i, ] *
                rep(c(0.12, -0.07), length.out = n_covariates)) *
            (tt - 0.5)
        } else {
          0
        }
        0.08 * s + 0.05 * j +
          (0.14 + 0.02 * i) * sin(2 * pi * tt) +
          0.04 * cos(pi * (i + j) * tt) +
          z_effect + stats::rnorm(n_obs, sd = 0.015)
      })
    })
  })

  list(Y = Y, Z = Z, time_obs = time_obs)
}

.fit_tiny_elbo <- function(S = 1L, L_f = 0L, L_s = 0L,
                           with_covariate = FALSE,
                           n_covariates = if (with_covariate) 1L else 0L,
                           n = 3L, p = 2L, n_obs = 6L,
                           M_f = NULL, M_s_each = NULL,
                           bool_var_spec_prob = FALSE,
                           anneal = NULL, maxit = 8L,
                           lambda_orth = 0,
                           convergence_rule = "elbo",
                           tol_abs = 0, tol_rel = 0,
                           seed = 2201L) {
  dat <- .make_tiny_elbo_data(
    S = S, n = n, p = p, n_obs = n_obs,
    with_covariate = with_covariate,
    n_covariates = n_covariates, seed = seed)
  set.seed(seed + 1L)

  bayesSYNC_multi(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = L_f, L_s = L_s,
    M_f = if (L_f > 0L) {
      if (is.null(M_f)) rep(1L, L_f) else M_f
    } else {
      integer()
    },
    M_s = lapply(seq_len(S), function(s) {
      if (L_s > 0L) {
        if (is.null(M_s_each)) rep(1L, L_s) else M_s_each
      } else {
        integer()
      }
    }),
    K = 2L, n_g = 21L,
    anneal = anneal, maxit = maxit, n_cpus = 1L,
    tol_abs = tol_abs, tol_rel = tol_rel,
    convergence_rule = convergence_rule,
    lambda_orth = lambda_orth,
    bool_scale = FALSE,
    bool_var_spec_prob = bool_var_spec_prob,
    verbose = FALSE
  )
}

.manual_sparse_elbo_components <- function(fit) {
  spike_shared <- 0
  spike_specific <- 0
  omega_shared <- 0
  omega_specific <- 0

  if (fit$L_f > 0L) {
    for (l in seq_len(fit$L_f)) {
      for (j in seq_len(fit$p)) {
        if (fit$bool_var_spec_prob) {
          elog_omega <- fit$mu_q_log_omega_a[j, l]
          elog_one_minus <- fit$mu_q_log_1_omega_a[j, l]
        } else {
          elog_omega <- fit$mu_q_log_omega_a[l]
          elog_one_minus <- fit$mu_q_log_1_omega_a[l]
        }
        spike_shared <- spike_shared +
          multiFSYNC:::.elbo_spike_slab_scalar(
            fit$mu_q_gamma_a[j, l],
            fit$mu_q_normal_a[j, l],
            fit$Sigma_q_normal_a[j, l],
            elog_omega, elog_one_minus
          )
      }
    }

    if (fit$bool_var_spec_prob) {
      for (l in seq_len(fit$L_f)) {
        for (j in seq_len(fit$p)) {
          omega_shared <- omega_shared +
            multiFSYNC:::.elbo_beta_prior_entropy(
              fit$c_1_omega_a[j, l], fit$d_1_omega_a[j, l],
              fit$c_0, fit$d_0
            )
        }
      }
    } else {
      for (l in seq_len(fit$L_f)) {
        omega_shared <- omega_shared +
          multiFSYNC:::.elbo_beta_prior_entropy(
            fit$c_1_omega_a[l], fit$d_1_omega_a[l],
            fit$c_0, fit$d_0
          )
      }
    }
  }

  if (fit$L_s > 0L) {
    for (s in seq_len(fit$S)) {
      for (l in seq_len(fit$L_s)) {
        for (j in seq_len(fit$p)) {
          if (fit$bool_var_spec_prob) {
            elog_omega <- fit$mu_q_log_omega_b[[s]][j, l]
            elog_one_minus <- fit$mu_q_log_1_omega_b[[s]][j, l]
          } else {
            elog_omega <- fit$mu_q_log_omega_b[[s]][l]
            elog_one_minus <- fit$mu_q_log_1_omega_b[[s]][l]
          }
          spike_specific <- spike_specific +
            multiFSYNC:::.elbo_spike_slab_scalar(
              fit$mu_q_gamma_b[[s]][j, l],
              fit$mu_q_normal_b[[s]][j, l],
              fit$Sigma_q_normal_b[[s]][j, l],
              elog_omega, elog_one_minus
            )
        }
      }

      if (fit$bool_var_spec_prob) {
        for (l in seq_len(fit$L_s)) {
          for (j in seq_len(fit$p)) {
            omega_specific <- omega_specific +
              multiFSYNC:::.elbo_beta_prior_entropy(
                fit$c_1_omega_b[[s]][j, l],
                fit$d_1_omega_b[[s]][j, l],
                fit$c_0, fit$d_0
              )
          }
        }
      } else {
        for (l in seq_len(fit$L_s)) {
          omega_specific <- omega_specific +
            multiFSYNC:::.elbo_beta_prior_entropy(
              fit$c_1_omega_b[[s]][l],
              fit$d_1_omega_b[[s]][l],
              fit$c_0, fit$d_0
            )
        }
      }
    }
  }

  c(
    spike_slab_shared = spike_shared,
    spike_slab_specific = spike_specific,
    omega_shared = omega_shared,
    omega_specific = omega_specific
  )
}

.rebuild_tiny_design <- function(dat, fit) {
  time_obs_flat <- unlist(
    lapply(seq_len(fit$S), function(s) dat$time_obs[[s]]),
    recursive = FALSE
  )
  grid_obj <- multiFSYNC:::get_grid_objects(
    time_obs_flat, fit$K,
    n_g = fit$n_g, time_g = fit$time_g,
    format_univ = TRUE
  )

  C <- vector("list", fit$S)
  first <- 1L
  for (s in seq_len(fit$S)) {
    last <- first + fit$n_s[s] - 1L
    C[[s]] <- grid_obj$C[first:last]
    first <- last + 1L
  }
  C
}

.direct_final_rss <- function(dat, fit) {
  C <- .rebuild_tiny_design(dat, fit)
  direct_rss <- lapply(seq_len(fit$S), function(s) {
    matrix(NA_real_, nrow = fit$n_s[s], ncol = fit$p)
  })

  for (s in seq_len(fit$S)) {
    for (i in seq_len(fit$n_s[s])) {
      for (j in seq_len(fit$p)) {
        direct_rss[[s]][i, j] <- multiFSYNC:::compute_rss_single(
          s = s, i = i, j = j,
          Y = dat$Y, C = C, list_cp_C = fit$list_cp_C,
          mu_q_nu_mu = fit$mu_q_nu_mu,
          Sigma_q_nu_mu = fit$Sigma_q_nu_mu,
          mu_q_nu_beta = fit$mu_q_nu_beta,
          Sigma_q_nu_beta = fit$Sigma_q_nu_beta,
          Z = dat$Z,
          mu_q_zeta = fit$mu_q_zeta,
          Sigma_q_zeta = fit$Sigma_q_zeta,
          mu_q_nu_phi = fit$mu_q_nu_phi,
          Sigma_q_nu_phi = fit$Sigma_q_nu_phi,
          mu_q_xi = fit$mu_q_xi,
          Sigma_q_xi = fit$Sigma_q_xi,
          mu_q_nu_psi = fit$mu_q_nu_psi,
          Sigma_q_nu_psi = fit$Sigma_q_nu_psi,
          mu_q_a = fit$mu_q_a, term_a = fit$term_a,
          mu_q_b_specific = fit$mu_q_b_specific,
          term_b_specific = fit$term_b_specific,
          L_f = fit$L_f, L_s = fit$L_s
        )
      }
    }
  }
  direct_rss
}

.expect_complete_elbo_trace <- function(fit, label,
                                        check_monotonicity = TRUE,
                                        minimum_length = 8L) {
  expect_true(is.numeric(fit$ELBO), info = label)
  expect_true(length(fit$ELBO) >= minimum_length, info = label)
  expect_equal(length(fit$ELBO_components), length(fit$ELBO),
               info = label)
  expect_equal(length(fit$ELBO_finite), length(fit$ELBO),
               info = label)
  expect_equal(length(fit$ELBO_temperature), length(fit$ELBO),
               info = label)
  expect_true(all(fit$ELBO_temperature == 1), info = label)
  expect_true(is.list(fit$linear_solver_diagnostics), info = label)
  expect_true(fit$linear_solver_diagnostics$total_calls > 0L, info = label)
  expect_identical(fit$elbo_t1_jitter_count, 0L, info = label)

  component_names <- names(fit$ELBO_components[[1L]])
  expect_identical(
    component_names, .complete_elbo_component_names,
    info = label
  )
  expect_true(
    length(component_names) > 0L &&
      all(nzchar(component_names)) &&
      !anyDuplicated(component_names),
    info = label
  )
  expect_true(all(vapply(
    fit$ELBO_components,
    function(x) identical(names(x), component_names),
    logical(1)
  )), info = label)
  expect_true(all(vapply(
    fit$ELBO_components,
    function(x) all(is.finite(x)),
    logical(1)
  )), info = label)
  expect_true(all(vapply(
    fit$ELBO_finite,
    function(x) isTRUE(x$all) &&
      isTRUE(x$total) &&
      all(x$components) &&
      length(x$nonfinite_components) == 0L,
    logical(1)
  )), info = label)

  component_totals <- vapply(
    fit$ELBO_components, sum, numeric(1))
  expect_equal(
    unname(fit$ELBO), unname(component_totals),
    tolerance = 1e-12, info = label
  )
  expect_equal(
    fit$elbo_result$total,
    sum(fit$elbo_result$components),
    tolerance = 1e-12, info = label
  )

  if (check_monotonicity) {
    elbo_diff <- diff(fit$ELBO)
    relative_tolerance <- fit$ELBO_diagnostics$tolerance *
      (1 + pmax(
        abs(head(fit$ELBO, -1L)),
        abs(tail(fit$ELBO, -1L))
      ))

    expect_true(
      all(is.finite(elbo_diff) &
            elbo_diff >= -relative_tolerance),
      info = label
    )
    expect_equal(
      fit$ELBO_diagnostics$differences, elbo_diff,
      tolerance = 1e-14, info = label
    )
    expect_equal(
      fit$ELBO_diagnostics$scaled_tolerance, relative_tolerance,
      tolerance = 1e-14, info = label
    )
    expect_identical(fit$ELBO_diagnostics$decrease_count, 0L,
                     info = label)
    expect_identical(fit$ELBO_decrease_count, 0L, info = label)
    expect_true(fit$ELBO_diagnostics$monotone_within_tolerance,
                info = label)
    expect_equal(
      fit$ELBO_diagnostics$minimum_difference, min(elbo_diff),
      tolerance = 1e-14, info = label
    )
    expect_equal(
      fit$ELBO_min_difference, min(elbo_diff),
      tolerance = 1e-14, info = label
    )
  }
}

test_that("complete T=1 ELBO is finite and monotone in all small model cases", {
  cases <- list(
    pure_mean = list(
      S = 1L, L_f = 0L, L_s = 0L, with_covariate = FALSE,
      bool_var_spec_prob = FALSE, seed = 2301L
    ),
    single_shared_factor = list(
      S = 1L, L_f = 1L, L_s = 0L, with_covariate = FALSE,
      bool_var_spec_prob = FALSE, seed = 2302L
    ),
    shared_and_specific_factors = list(
      S = 2L, L_f = 1L, L_s = 1L, with_covariate = FALSE,
      bool_var_spec_prob = FALSE, seed = 2303L
    ),
    one_covariate = list(
      S = 2L, L_f = 0L, L_s = 0L, with_covariate = TRUE,
      bool_var_spec_prob = FALSE, seed = 2304L
    ),
    variable_factor_omega = list(
      S = 2L, L_f = 1L, L_s = 1L, with_covariate = FALSE,
      bool_var_spec_prob = TRUE, seed = 2305L
    )
  )

  for (case_name in names(cases)) {
    expect_warning(
      fit <- do.call(.fit_tiny_elbo, cases[[case_name]]),
      regexp = "Max iterations reached before convergence"
    )
    .expect_complete_elbo_trace(fit, case_name)

    expected_hierarchy <- if (cases[[case_name]]$bool_var_spec_prob) {
      "variable_factor"
    } else {
      "factor"
    }
    expect_identical(fit$elbo_result$omega_hierarchy,
                     expected_hierarchy, info = case_name)
    expect_identical(fit$bool_var_spec_prob,
                     cases[[case_name]]$bool_var_spec_prob,
                     info = case_name)
    manual_sparse <- .manual_sparse_elbo_components(fit)
    expect_equal(
      fit$elbo_result$components[names(manual_sparse)],
      manual_sparse,
      tolerance = 1e-12,
      info = case_name
    )
    expect_identical(fit$elbo_role, "objective", info = case_name)
    expect_true(fit$elbo_objective_valid, info = case_name)
    expect_identical(
      fit$elbo_invalid_reasons, character(), info = case_name
    )
  }
})

test_that("high-coupling Gauss-Seidel sweep has a monotone T=1 ELBO", {
  expect_warning(
    fit <- .fit_tiny_elbo(
      S = 2L, n = 4L, p = 2L,
      L_f = 2L, M_f = c(2L, 2L),
      L_s = 2L, M_s_each = c(2L, 2L),
      with_covariate = TRUE, n_covariates = 2L,
      bool_var_spec_prob = TRUE,
      maxit = 8L, seed = 2351L
    ),
    regexp = "Max iterations reached before convergence"
  )

  .expect_complete_elbo_trace(fit, "high-coupling Gauss-Seidel")
  expect_identical(fit$d, 2L)
  expect_identical(fit$L_f, 2L)
  expect_identical(fit$M_f, c(2L, 2L))
  expect_identical(fit$L_s, 2L)
  expect_true(all(vapply(
    fit$M_s, identical, logical(1), c(2L, 2L)
  )))
  expect_identical(fit$elbo_result$omega_hierarchy, "variable_factor")

  dat <- .make_tiny_elbo_data(
    S = 2L, n = 4L, p = 2L, n_obs = 6L,
    with_covariate = TRUE, n_covariates = 2L,
    seed = 2351L
  )
  direct_rss <- .direct_final_rss(dat, fit)
  direct_rss_sum <- do.call(rbind, lapply(direct_rss, colSums))
  expect_equal(direct_rss, fit$expected_rss, tolerance = 1e-11)
  expect_equal(
    direct_rss_sum, fit$expected_rss_sum, tolerance = 1e-11
  )
  expect_equal(
    fit$elbo_result$rss$expected_rss,
    direct_rss,
    tolerance = 1e-11
  )
  expect_equal(
    fit$elbo_result$rss$expected_rss_sum,
    direct_rss_sum,
    tolerance = 1e-11
  )
  expect_identical(
    fit$elbo_result$rss$source, "variance_update_cache"
  )

  total_obs_sj <- matrix(0, nrow = fit$S, ncol = fit$p)
  for (s in seq_len(fit$S)) {
    for (j in seq_len(fit$p)) {
      total_obs_sj[s, j] <- sum(vapply(
        seq_len(fit$n_s[s]),
        function(i) length(dat$Y[[s]][[i]][[j]]),
        integer(1)
      ))
    }
  }
  manual_data_likelihood <- -0.5 * sum(
    total_obs_sj * log(2 * pi) +
      total_obs_sj * fit$mu_q_log_sigsq_eps +
      fit$mu_q_recip_sigsq_eps * direct_rss_sum
  )
  expect_equal(
    fit$elbo_result$components[["data_likelihood"]],
    manual_data_likelihood,
    tolerance = 1e-11
  )
})

test_that("annealing records the ordinary ELBO only after reaching T=1", {
  expect_warning(
    fit <- .fit_tiny_elbo(
      S = 1L, L_f = 1L, L_s = 0L,
      anneal = c(1, 1.25, 2), maxit = 9L,
      bool_var_spec_prob = FALSE, seed = 2401L
    ),
    regexp = "Max iterations reached before convergence"
  )

  .expect_complete_elbo_trace(fit, "annealing")
  expect_true(fit$annealing_completed)
  expect_equal(fit$final_temperature, 1, tolerance = 1e-14)
  expect_true(length(fit$ELBO) < fit$i_iter)
  expect_identical(
    fit$ELBO_diagnostics$history_length,
    length(fit$ELBO)
  )
})

test_that("valid ELBO stopping can terminate a no-jitter T=1 fit", {
  fit <- .fit_tiny_elbo(
    S = 1L, L_f = 1L, L_s = 0L,
    convergence_rule = "elbo",
    tol_abs = 1e20, tol_rel = 0,
    maxit = 8L,
    bool_var_spec_prob = FALSE,
    seed = 2451L
  )

  expect_true(fit$converged)
  expect_identical(fit$convergence_rule_requested, "elbo")
  expect_identical(fit$convergence_rule, "elbo")
  expect_identical(fit$elbo_role, "objective")
  expect_true(fit$elbo_objective_valid)
  expect_identical(fit$elbo_invalid_reasons, character())
  expect_identical(fit$elbo_t1_jitter_count, 0L)
  expect_identical(fit$i_iter, 2L)
  expect_length(fit$ELBO, 2L)
})

test_that("positive lambda_orth downgrades ELBO stopping to diagnostic mode", {
  expect_warning(
    fit <- .fit_tiny_elbo(
      S = 1L, L_f = 1L, L_s = 0L,
      lambda_orth = 0.05,
      convergence_rule = "elbo",
      tol_abs = 1e20,
      bool_var_spec_prob = FALSE,
      seed = 2501L
    ),
    regexp = "downgrading the requested objective-based convergence rule"
  )

  .expect_complete_elbo_trace(
    fit, "positive lambda_orth", check_monotonicity = FALSE,
    minimum_length = 2L)
  expect_identical(fit$convergence_rule_requested, "elbo")
  expect_identical(fit$convergence_rule, "parameters")
  expect_identical(fit$elbo_role, "diagnostic_only")
  expect_false(fit$elbo_objective_valid)
  expect_identical(
    fit$elbo_invalid_reasons,
    "positive_soft_orthogonality_penalty"
  )
  expect_equal(fit$lambda_orth, 0.05)
})

test_that("a recorded T=1 jitter event downgrades ELBO stopping", {
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
      # The linear-prior covariance is inverted during initialisation, before
      # a T=1 sweep begins. Inject only into the first coordinate solve so the
      # event is correctly classified as T=1 jitter by the main loop.
      if (!injected_jitter &&
          !identical(context, "linear_prior_covariance")) {
        multiFSYNC:::.record_spd_diagnostic(
          context = paste0(context, "::test-injected-jitter"),
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
    fit <- .fit_tiny_elbo(
      S = 1L, L_f = 1L, L_s = 0L,
      convergence_rule = "elbo",
      tol_abs = 1e20, tol_rel = 0,
      maxit = 8L,
      bool_var_spec_prob = FALSE,
      seed = 2551L
    ),
    regexp = "Adaptive Cholesky jitter.*downgrading"
  )

  expect_true(injected_jitter)
  expect_true(fit$converged)
  expect_identical(fit$convergence_rule_requested, "elbo")
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
    "test-injected-jitter",
    fit$linear_solver_diagnostics$events$context,
    fixed = TRUE
  )))
})
