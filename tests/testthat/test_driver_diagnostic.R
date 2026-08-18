.make_driver_diagnostic_test_data <- function(seed = 9611L) {
  multiFSYNC::simulate_multi_study_osullivan(
    S = 1L, n_s = 4L, p = 3L, d = 0L,
    L_f = 1L, L_s = 1L,
    M_f = 1L, M_s = list(1L),
    K = 3L, n_obs = 7L, common_grid = TRUE,
    sigma_eps = 0.08, seed = seed, use_explicit = FALSE
  )
}

.fit_driver_diagnostic_test <- function(dat, diagnostic = FALSE,
                                        dense_gate_sweeps = 0L,
                                        maxit = 4L, seed = 9612L) {
  arguments <- list(
    Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L),
    K = 2L, n_g = 19L,
    anneal = NULL, maxit = maxit,
    tol_abs = 0, tol_rel = 0,
    convergence_rule = "parameters",
    bool_scale = FALSE, verbose = FALSE, seed = seed
  )
  suppressWarnings(if (diagnostic) {
    do.call(
      multiFSYNC:::.bayesSYNC_multi_driver_diagnostic,
      c(
        arguments,
        list(control = list(
          trace_sweeps = seq_len(maxit),
          dense_gate_sweeps = dense_gate_sweeps
        ))
      )
    )
  } else {
    do.call(multiFSYNC::bayesSYNC_multi, arguments)
  })
}

testthat::test_that("private driver tracing is numerically and RNG inert", {
  testthat::expect_false(
    "driver_diagnostic_control" %in%
      names(formals(multiFSYNC::bayesSYNC_multi))
  )
  dat <- .make_driver_diagnostic_test_data()
  plain <- .fit_driver_diagnostic_test(dat, diagnostic = FALSE)
  seed_after_plain <- .Random.seed
  traced <- .fit_driver_diagnostic_test(dat, diagnostic = TRUE)
  seed_after_traced <- .Random.seed

  testthat::expect_false("driver_trace" %in% names(plain))
  testthat::expect_s3_class(traced$driver_trace, "data.frame")
  testthat::expect_identical(seed_after_traced, seed_after_plain)
  traced$driver_trace <- NULL
  traced$driver_diagnostic_control <- NULL
  traced$random_scale_calibration_diagnostics <- NULL
  testthat::expect_identical(traced, plain)
})

testthat::test_that("driver trace covers every coupled factor block", {
  dat <- .make_driver_diagnostic_test_data(seed = 9621L)
  fit <- .fit_driver_diagnostic_test(dat, diagnostic = TRUE, maxit = 3L)
  trace <- fit$driver_trace
  expected <- c(
    "shared_time_function", "specific_time_function",
    "shared_score", "specific_score",
    "shared_loading", "specific_loading"
  )
  testthat::expect_setequal(unique(trace$block), expected)
  testthat::expect_setequal(unique(trace$iteration), 1:3)
  testthat::expect_true(all(diff(trace$sequence) == 1L))
  testthat::expect_true(all(is.finite(trace$linear_driver_norm)))
  testthat::expect_true(all(is.finite(trace$precision_trace)))
  testthat::expect_true(all(trace$precision_trace > 0))
  loading <- grepl("_loading$", trace$block)
  testthat::expect_true(all(is.finite(trace$unrestricted_ppi_mean[loading])))
  testthat::expect_true(all(trace$unrestricted_ppi_mean[loading] >= 0))
  testthat::expect_true(all(trace$unrestricted_ppi_mean[loading] <= 1))
})

testthat::test_that("dense gate is bounded and releases the original update", {
  dat <- .make_driver_diagnostic_test_data(seed = 9631L)
  fit <- .fit_driver_diagnostic_test(
    dat, diagnostic = TRUE, dense_gate_sweeps = 2L, maxit = 5L)
  loading <- fit$driver_trace[
    grepl("_loading$", fit$driver_trace$block), , drop = FALSE]

  testthat::expect_true(all(loading$gate_active[loading$iteration <= 2L]))
  testthat::expect_true(all(loading$ppi_mean[loading$iteration <= 2L] == 1))
  testthat::expect_false(any(loading$gate_active[loading$iteration > 2L]))
  testthat::expect_equal(
    loading$ppi_mean[loading$iteration > 2L],
    loading$unrestricted_ppi_mean[loading$iteration > 2L],
    tolerance = 0
  )
  testthat::expect_true(
    isTRUE(attr(fit$driver_trace, "dense_gate_released")))
  testthat::expect_gte(attr(fit$driver_trace, "post_release_sweeps"), 1L)

  differences <- diff(fit$ELBO)
  tolerances <- fit$ELBO_diagnostics$scaled_tolerance
  testthat::expect_true(all(differences >= -tolerances))
})

testthat::test_that("dense gate requires a post-release sweep", {
  dat <- .make_driver_diagnostic_test_data(seed = 9641L)
  testthat::expect_error(
    .fit_driver_diagnostic_test(
      dat, diagnostic = TRUE, dense_gate_sweeps = 3L, maxit = 3L),
    "maxit greater than"
  )
})

testthat::test_that("data-free calibration fixes moments and preserves directions", {
  helper_formals <- names(formals(
    multiFSYNC:::.calibrate_random_factor_state
  ))
  testthat::expect_false(any(
    c("Y", "Z", "true_params", "truth") %in% helper_formals
  ))

  time_g <- seq(0, 1, length.out = 21L)
  C_g <- cbind(1, time_g, sin(2 * pi * time_g), cos(2 * pi * time_g))
  phi_before <- matrix(c(2, -1, 0.5, 3), ncol = 1L)
  psi_before <- matrix(c(-1, 2, 1.5, 0.25), ncol = 1L)
  zeta_before <- matrix(c(-2, -1, 0.5, 3), ncol = 1L)
  xi_before <- matrix(c(1, -3, 2, 0.5), ncol = 1L)
  slab_a_before <- matrix(c(-2, 1, 4), ncol = 1L)
  slab_b_before <- matrix(c(3, -1, 0.5), ncol = 1L)

  calibrated <- multiFSYNC:::.calibrate_random_factor_state(
    C_g = C_g, time_g = time_g,
    mu_q_nu_phi = list(phi_before),
    mu_q_zeta = list(list(zeta_before)),
    mu_q_normal_a = slab_a_before,
    Sigma_q_normal_a = matrix(1, 3L, 1L),
    mu_q_gamma_a = matrix(0.5, 3L, 1L),
    mu_q_nu_psi = list(list(psi_before)),
    mu_q_xi = list(list(xi_before)),
    mu_q_normal_b = list(slab_b_before),
    Sigma_q_normal_b = list(matrix(1, 3L, 1L)),
    mu_q_gamma_b = list(matrix(0.5, 3L, 1L)),
    S = 1L, L_f = 1L, L_s = 1L,
    M_f = 1L, M_s = list(1L), mode = "all"
  )

  direction_cosine <- function(before, after) {
    sum(before * after) / sqrt(sum(before^2) * sum(after^2))
  }
  testthat::expect_equal(
    direction_cosine(phi_before, calibrated$mu_q_nu_phi[[1L]]), 1,
    tolerance = 1e-14
  )
  testthat::expect_equal(
    direction_cosine(zeta_before, calibrated$mu_q_zeta[[1L]][[1L]]), 1,
    tolerance = 1e-14
  )
  testthat::expect_equal(
    direction_cosine(slab_a_before, calibrated$mu_q_normal_a), 1,
    tolerance = 1e-14
  )
  weights <- multiFSYNC:::.trap_weights(time_g)
  testthat::expect_equal(
    sum(weights * as.vector(C_g %*% calibrated$mu_q_nu_phi[[1L]])^2),
    1, tolerance = 1e-13
  )
  testthat::expect_equal(
    sqrt(mean(calibrated$mu_q_zeta[[1L]][[1L]]^2)),
    1, tolerance = 1e-14
  )
  testthat::expect_equal(
    sqrt(mean(calibrated$mu_q_normal_a^2)),
    1, tolerance = 1e-14
  )
  testthat::expect_true(all(calibrated$diagnostics$success))
  testthat::expect_equal(
    calibrated$diagnostics$after,
    rep(1, nrow(calibrated$diagnostics)),
    tolerance = 1e-13
  )
})

testthat::test_that("development-only function initializers have registered scales", {
  time_g <- seq(0, 1, length.out = 31L)
  C_g <- cbind(1, time_g, sin(2 * pi * time_g), cos(2 * pi * time_g))
  phi_before <- cbind(
    c(2, -1, 0.5, 3),
    c(-0.25, 1.5, -2, 0.75)
  )
  psi_before <- cbind(
    c(-1, 2, 1.5, 0.25),
    c(0.5, -0.75, 2.5, -1)
  )
  call_calibration <- function(mode) {
    multiFSYNC:::.calibrate_random_factor_state(
      C_g = C_g, time_g = time_g,
      mu_q_nu_phi = list(phi_before),
      mu_q_zeta = list(list(matrix(1, 4L, 2L))),
      mu_q_normal_a = matrix(c(-2, 1, 4), ncol = 1L),
      Sigma_q_normal_a = matrix(1, 3L, 1L),
      mu_q_gamma_a = matrix(0.5, 3L, 1L),
      mu_q_nu_psi = list(list(psi_before)),
      mu_q_xi = list(list(matrix(1, 4L, 2L))),
      mu_q_normal_b = list(matrix(c(3, -1, 0.5), ncol = 1L)),
      Sigma_q_normal_b = list(matrix(1, 3L, 1L)),
      mu_q_gamma_b = list(matrix(0.5, 3L, 1L)),
      S = 1L, L_f = 1L, L_s = 1L,
      M_f = 2L, M_s = list(2L), mode = mode
    )
  }

  jaoua <- call_calibration("coefficient_iid")
  testthat::expect_equal(jaoua$mu_q_nu_phi[[1L]][, 1L],
                         phi_before[, 1L], tolerance = 0)
  testthat::expect_equal(jaoua$mu_q_nu_phi[[1L]][, 2L],
                         phi_before[, 2L] * sqrt(2), tolerance = 1e-14)
  testthat::expect_equal(jaoua$mu_q_nu_psi[[1L]][[1L]][, 2L],
                         psi_before[, 2L] * sqrt(2), tolerance = 1e-14)

  gram_decay <- call_calibration("function_m_decay")
  weights <- multiFSYNC:::.trap_weights(time_g)
  shared_energy <- vapply(1:2, function(component) {
    dense <- C_g %*% gram_decay$mu_q_nu_phi[[1L]][, component]
    sum(weights * dense^2)
  }, numeric(1L))
  specific_energy <- vapply(1:2, function(component) {
    dense <- C_g %*% gram_decay$mu_q_nu_psi[[1L]][[1L]][, component]
    sum(weights * dense^2)
  }, numeric(1L))
  testthat::expect_equal(shared_energy, c(1, 1 / 2), tolerance = 1e-13)
  testthat::expect_equal(specific_energy, c(1, 1 / 2), tolerance = 1e-13)
  testthat::expect_true(all(gram_decay$diagnostics$success))
})

testthat::test_that("pre-score path is bounded, traced, and ELBO monotone", {
  dat <- .make_driver_diagnostic_test_data(seed = 9651L)
  arguments <- list(
    Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L),
    K = 2L, n_g = 19L, anneal = NULL, maxit = 5L,
    tol_abs = 0, tol_rel = 0,
    convergence_rule = "parameters",
    bool_scale = FALSE, verbose = FALSE, seed = 9652L,
    control = list(
      trace_sweeps = 1:5,
      dense_gate_sweeps = 0L,
      random_scale_calibration = "all",
      pre_score_sweeps = 1L
    )
  )
  fit <- suppressWarnings(do.call(
    multiFSYNC:::.bayesSYNC_multi_driver_diagnostic, arguments
  ))
  pre <- fit$driver_trace$phase == "pre_score"
  testthat::expect_true(any(pre))
  testthat::expect_true(all(fit$driver_trace$iteration[pre] == 1L))
  testthat::expect_true(all(
    fit$driver_trace$block[pre] %in% c("shared_score", "specific_score")
  ))
  testthat::expect_false(any(
    fit$driver_trace$phase[fit$driver_trace$iteration > 1L] == "pre_score"
  ))
  testthat::expect_true(all(
    diff(fit$ELBO) >= -fit$ELBO_diagnostics$scaled_tolerance
  ))
  testthat::expect_true(all(
    fit$random_scale_calibration_diagnostics$success
  ))
})
