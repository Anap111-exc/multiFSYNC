.make_scale_trace_test_data <- function(seed = 9311L) {
  multiFSYNC::simulate_multi_study_osullivan(
    S = 1L, n_s = 3L, p = 3L, d = 1L,
    L_f = 1L, L_s = 1L,
    M_f = 1L, M_s = list(1L),
    K = 3L, n_obs = 6L, common_grid = TRUE,
    sigma_eps = 0.08, seed = seed, use_explicit = FALSE
  )
}

.fit_scale_trace_test <- function(dat, anneal = NULL, maxit = 3L,
                                  trace_control = NULL,
                                  traced = FALSE, seed = 9312L) {
  arguments <- list(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L),
    K = 2L, n_g = 17L,
    anneal = anneal, maxit = maxit,
    tol_abs = 0, tol_rel = 0,
    convergence_rule = "parameters",
    bool_scale = FALSE, verbose = FALSE, seed = seed
  )
  suppressWarnings(if (traced) {
    do.call(
      multiFSYNC:::.bayesSYNC_multi_scale_trace,
      c(arguments, list(trace_control = trace_control))
    )
  } else {
    do.call(multiFSYNC::bayesSYNC_multi, arguments)
  })
}

testthat::test_that("private scale tracing is numerically and RNG inert", {
  testthat::expect_false(
    "trace_control" %in% names(formals(multiFSYNC::bayesSYNC_multi))
  )
  testthat::expect_false(
    "scale_trace_control" %in% names(formals(multiFSYNC::bayesSYNC_multi))
  )

  dat <- .make_scale_trace_test_data()
  plain <- .fit_scale_trace_test(dat, traced = FALSE)
  seed_after_plain <- .Random.seed
  traced <- .fit_scale_trace_test(
    dat, traced = TRUE,
    trace_control = list(
      include_initial = TRUE,
      include_annealing = FALSE,
      t1_sweeps = 1:3
    )
  )
  seed_after_traced <- .Random.seed

  testthat::expect_false("scale_trace" %in% names(plain))
  testthat::expect_s3_class(traced$scale_trace, "data.frame")
  testthat::expect_identical(seed_after_traced, seed_after_plain)
  traced_without_trace <- traced
  traced_without_trace$scale_trace <- NULL
  testthat::expect_equal(
    traced_without_trace, plain, tolerance = 1e-12
  )
  testthat::expect_identical(traced_without_trace, plain)

  core_objects <- c(
    "mu_q_nu_mu", "Sigma_q_nu_mu",
    "mu_q_nu_beta", "Sigma_q_nu_beta",
    "mu_q_nu_phi", "Sigma_q_nu_phi",
    "mu_q_nu_psi", "Sigma_q_nu_psi",
    "mu_q_zeta", "Sigma_q_zeta",
    "mu_q_xi", "Sigma_q_xi",
    "mu_q_normal_a", "Sigma_q_normal_a", "mu_q_gamma_a", "mu_q_a",
    "mu_q_normal_b", "Sigma_q_normal_b", "mu_q_gamma_b",
    "mu_q_b_specific",
    "mu_q_recip_sigsq_eps", "mu_q_recip_sigsq_mu",
    "mu_q_recip_sigsq_beta", "mu_q_recip_sigsq_phi",
    "mu_q_recip_sigsq_psi",
    "c_1_omega_a", "d_1_omega_a",
    "c_1_omega_b", "d_1_omega_b",
    "expected_rss", "expected_rss_sum",
    "ELBO", "ELBO_components"
  )
  for (object in core_objects) {
    testthat::expect_equal(
      traced[[object]], plain[[object]], tolerance = 1e-12,
      info = object
    )
  }
})

testthat::test_that("scale trace has the requested annealing and T1 states", {
  dat <- .make_scale_trace_test_data(seed = 9321L)
  anneal <- c(1, 1.5, 4)
  fit <- .fit_scale_trace_test(
    dat, anneal = anneal, maxit = 7L, traced = TRUE, seed = 9322L,
    trace_control = list(
      include_initial = TRUE,
      include_annealing = TRUE,
      t1_sweeps = c(1L, 2L)
    )
  )
  trace <- fit$scale_trace

  testthat::expect_identical(
    trace$stage,
    c("initial", rep("annealing", 3L), rep("t1", 2L))
  )
  testthat::expect_identical(trace$iteration, 0:5)
  testthat::expect_identical(
    trace$t1_sweep,
    c(0L, rep(NA_integer_, 3L), 1L, 2L)
  )

  ladder <- multiFSYNC:::get_annealing_ladder_(anneal, verbose = FALSE)
  testthat::expect_equal(
    trace$temperature[1:4],
    c(1 / ladder[1L], 1 / ladder[1:3]),
    tolerance = 1e-14
  )
  testthat::expect_equal(trace$temperature[5:6], c(1, 1))

  testthat::expect_true(is.finite(trace$elbo[1L]))
  testthat::expect_identical(
    trace$objective_scope[1L], "prevariance_conditional"
  )
  testthat::expect_identical(
    trace$elbo_reason[1L],
    paste(
      "conditional objective; variance-only terms unavailable before",
      "the first variance update"
    )
  )
  annealing_rows <- trace$stage == "annealing"
  testthat::expect_true(all(is.finite(trace$elbo[annealing_rows])))
  testthat::expect_true(all(
    trace$objective_scope[annealing_rows] == "full"
  ))
  testthat::expect_true(all(
    trace$elbo_reason[annealing_rows] ==
      "complete tempered objective"
  ))
  t1_rows <- trace$stage == "t1"
  testthat::expect_equal(
    trace$elbo[t1_rows], fit$ELBO[1:2], tolerance = 1e-12
  )
  testthat::expect_equal(
    trace$elbo_data_likelihood[t1_rows],
    vapply(
      fit$ELBO_components[1:2],
      function(x) unname(x[["data_likelihood"]]),
      numeric(1)
    ),
    tolerance = 1e-12
  )
  testthat::expect_true(all(is.finite(trace$expected_rss_total)))
  testthat::expect_lt(
    max(trace$fitted_cache_max_abs_diff, na.rm = TRUE), 1e-12
  )

  required_columns <- c(
    "objective_scope", "observed_rms", "fitted_rms",
    "mean_contribution_rms", "beta_contribution_rms",
    "shared_contribution_rms", "specific_contribution_rms",
    "raw_observation_count", "raw_observation_nonfinite_count",
    "raw_prediction_value_count",
    "raw_fitted_nonfinite_count",
    "raw_mean_contribution_nonfinite_count",
    "raw_beta_contribution_nonfinite_count",
    "raw_shared_contribution_nonfinite_count",
    "raw_specific_contribution_nonfinite_count",
    "raw_factor_contribution_nonfinite_count",
    "raw_all_prediction_nonfinite_count",
    "noise_precision_min", "noise_precision_median",
    "noise_precision_max",
    "mu_rms", "beta_rms", "phi_rms", "psi_rms",
    "zeta_rms", "xi_rms", "a_rms", "b_rms",
    "mu_q_second_moment", "beta_q_second_moment",
    "score_q_second_moment",
    "slab_loading_q_second_moment",
    "effective_loading_q_second_moment",
    "dense_phi_mean_rms", "dense_phi_qsecond_rms",
    "dense_psi_mean_rms", "dense_psi_qsecond_rms",
    "ppi_min", "ppi_mean", "ppi_max", "ppi_gt_half",
    "raw_ppi_count", "raw_ppi_nonfinite_count",
    "raw_ppi_out_of_range_count",
    "omega_min", "omega_mean", "omega_max",
    "shared_omega_mean", "specific_omega_mean",
    "smooth_precision_min", "smooth_precision_median",
    "smooth_precision_max",
    "half_cauchy_aux_min", "half_cauchy_aux_median",
    "half_cauchy_aux_max"
  )
  testthat::expect_true(all(required_columns %in% names(trace)))

  raw_zero_columns <- c(
    "raw_observation_nonfinite_count",
    "raw_fitted_nonfinite_count",
    "raw_mean_contribution_nonfinite_count",
    "raw_beta_contribution_nonfinite_count",
    "raw_shared_contribution_nonfinite_count",
    "raw_specific_contribution_nonfinite_count",
    "raw_factor_contribution_nonfinite_count",
    "raw_all_prediction_nonfinite_count",
    "raw_ppi_nonfinite_count",
    "raw_ppi_out_of_range_count"
  )
  testthat::expect_true(all(
    as.matrix(trace[raw_zero_columns]) == 0
  ))
  testthat::expect_true(all(trace$raw_observation_count > 0L))
  testthat::expect_equal(
    trace$raw_prediction_value_count,
    5L * trace$raw_observation_count,
    tolerance = 0
  )
  testthat::expect_true(all(trace$raw_ppi_count > 0L))
})

testthat::test_that("within-sweep block tracing is complete and inert", {
  dat <- .make_scale_trace_test_data(seed = 9331L)
  anneal <- c(1, 1.5, 4)
  plain <- .fit_scale_trace_test(
    dat, anneal = anneal, maxit = 7L, traced = FALSE, seed = 9332L
  )
  seed_after_plain <- .Random.seed
  traced <- .fit_scale_trace_test(
    dat, anneal = anneal, maxit = 7L, traced = TRUE, seed = 9332L,
    trace_control = list(
      include_initial = FALSE,
      include_annealing = FALSE,
      t1_sweeps = integer(),
      block_annealing_sweeps = 1:2
    )
  )
  seed_after_traced <- .Random.seed

  testthat::expect_identical(seed_after_traced, seed_after_plain)
  traced_without_trace <- traced
  traced_without_trace$scale_trace <- NULL
  testthat::expect_identical(traced_without_trace, plain)

  trace <- traced$scale_trace
  expected_blocks <- c(
    "before_updates", "after_mu", "after_beta", "after_phi",
    "after_psi", "after_zeta", "after_xi",
    "after_shared_loadings", "after_specific_loadings",
    "after_variances", "after_omega"
  )
  testthat::expect_identical(
    trace$stage, rep("annealing_block", 2L * length(expected_blocks))
  )
  testthat::expect_identical(
    trace$iteration,
    rep(1:2, each = length(expected_blocks))
  )
  testthat::expect_identical(
    trace$block, rep(expected_blocks, times = 2L)
  )
  first_sweep <- trace$iteration == 1L
  prevariance_blocks <- c(
    "before_updates", "after_mu", "after_beta", "after_phi",
    "after_psi", "after_zeta", "after_xi",
    "after_shared_loadings", "after_specific_loadings"
  )
  testthat::expect_true(all(
    trace$objective_scope[
      first_sweep & trace$block %in% prevariance_blocks
    ] == "prevariance_conditional"
  ))
  testthat::expect_true(all(
    trace$objective_scope[
      first_sweep & trace$block %in% c("after_variances", "after_omega")
    ] == "full"
  ))
  testthat::expect_true(all(
    trace$objective_scope[trace$iteration == 2L] == "full"
  ))
  testthat::expect_true(all(is.finite(trace$elbo)))
  testthat::expect_true(all(is.finite(trace$expected_rss_total)))
  testthat::expect_lt(
    max(trace$fitted_cache_max_abs_diff, na.rm = TRUE), 1e-12
  )
  testthat::expect_true(all(trace$raw_all_prediction_nonfinite_count == 0L))
  testthat::expect_true(all(trace$raw_ppi_nonfinite_count == 0L))
  testthat::expect_true(all(trace$raw_ppi_out_of_range_count == 0L))
})

testthat::test_that("block trace sweep validation rejects invalid values", {
  testthat::expect_error(
    multiFSYNC:::.validate_scale_trace_control(list(
      block_annealing_sweeps = c(0, 1)
    )),
    "block_annealing_sweeps must contain positive integers"
  )
  testthat::expect_error(
    multiFSYNC:::.validate_scale_trace_control(list(
      block_annealing_sweeps = 1.5
    )),
    "block_annealing_sweeps must contain positive integers"
  )
})

testthat::test_that("private trace option is restored after fit errors", {
  option_name <- multiFSYNC:::.multiFSYNC_scale_trace_option
  original_options <- options()
  originally_existed <- option_name %in% names(original_options)
  original_value <- original_options[[option_name]]
  on.exit({
    options(stats::setNames(
      list(if (originally_existed) original_value else NULL),
      option_name
    ))
  }, add = TRUE)

  sentinel <- structure(
    list(marker = "pre-existing private option"),
    class = "scale_trace_test_sentinel"
  )
  options(stats::setNames(list(sentinel), option_name))
  testthat::expect_error(
    multiFSYNC:::.bayesSYNC_multi_scale_trace(
      Y = NULL, L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      trace_control = list(
        include_initial = TRUE,
        include_annealing = FALSE,
        t1_sweeps = 1L
      ),
      verbose = FALSE
    ),
    "Y must be a non-empty list"
  )
  testthat::expect_identical(getOption(option_name), sentinel)

  options(stats::setNames(list(NULL), option_name))
  testthat::expect_false(option_name %in% names(options()))
  testthat::expect_error(
    multiFSYNC:::.bayesSYNC_multi_scale_trace(
      Y = NULL, L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      trace_control = list(
        include_initial = FALSE,
        include_annealing = TRUE,
        t1_sweeps = integer()
      ),
      verbose = FALSE
    ),
    "Y must be a non-empty list"
  )
  testthat::expect_false(option_name %in% names(options()))
})
