.make_pre_score_interface_data <- function(seed = 9761L) {
  multiFSYNC::simulate_multi_study_osullivan(
    S = 1L, n_s = 4L, p = 3L, d = 0L,
    L_f = 1L, L_s = 1L,
    M_f = 1L, M_s = list(1L),
    K = 3L, n_obs = 7L, common_grid = TRUE,
    sigma_eps = 0.08, seed = seed, use_explicit = FALSE
  )
}

.pre_score_interface_args <- function(data, seed = 9762L) {
  list(
    Y = data$Y, Z = NULL, time_obs = data$time_obs,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L),
    K = 2L, n_g = 19L, anneal = NULL, maxit = 4L,
    tol_abs = 0, tol_rel = 0, convergence_rule = "parameters",
    bool_scale = FALSE, verbose = FALSE, seed = seed
  )
}

testthat::test_that("public pre-score interface fixes private controls", {
  data <- .make_pre_score_interface_data()
  fit <- suppressWarnings(do.call(
    multiFSYNC::bayesSYNC_multi_pre_score,
    c(.pre_score_interface_args(data), list(
      pre_score_sweeps = 1L, trace_sweeps = 1:2
    ))
  ))

  testthat::expect_identical(
    fit$driver_diagnostic_control$pre_score_sweeps, 1L
  )
  testthat::expect_identical(
    fit$driver_diagnostic_control$dense_gate_sweeps, 0L
  )
  testthat::expect_identical(
    fit$driver_diagnostic_control$random_scale_calibration, "none"
  )
  testthat::expect_identical(
    fit$pre_score_interface$interface, "bayesSYNC_multi_pre_score"
  )
  pre <- fit$driver_trace$phase == "pre_score"
  testthat::expect_true(any(pre))
  testthat::expect_true(all(fit$driver_trace$iteration[pre] == 1L))
  testthat::expect_true(all(
    fit$driver_trace$block[pre] %in% c("shared_score", "specific_score")
  ))
  testthat::expect_false(any(
    fit$driver_trace$phase[fit$driver_trace$iteration > 1L] == "pre_score"
  ))
})

testthat::test_that("zero pre-score is numerically identical to the core path", {
  data <- .make_pre_score_interface_data(seed = 9771L)
  arguments <- .pre_score_interface_args(data, seed = 9772L)
  plain <- suppressWarnings(do.call(multiFSYNC::bayesSYNC_multi, arguments))
  wrapped <- suppressWarnings(do.call(
    multiFSYNC::bayesSYNC_multi_pre_score,
    c(arguments, list(pre_score_sweeps = 0L, trace_sweeps = 1L))
  ))
  wrapped$driver_trace <- NULL
  wrapped$driver_diagnostic_control <- NULL
  wrapped$random_scale_calibration_diagnostics <- NULL
  wrapped$pre_score_interface <- NULL
  testthat::expect_identical(wrapped, plain)
})

testthat::test_that("public pre-score interface validates its narrow contract", {
  data <- .make_pre_score_interface_data(seed = 9781L)
  arguments <- .pre_score_interface_args(data, seed = 9782L)
  testthat::expect_error(
    do.call(
      multiFSYNC::bayesSYNC_multi_pre_score,
      c(arguments, list(pre_score_sweeps = -1L))
    ),
    "pre_score_sweeps"
  )
  testthat::expect_error(
    do.call(
      multiFSYNC::bayesSYNC_multi_pre_score,
      c(arguments, list(trace_sweeps = 0L))
    ),
    "trace_sweeps"
  )
  testthat::expect_error(
    do.call(
      multiFSYNC::bayesSYNC_multi_pre_score,
      c(arguments, list(control = list(pre_score_sweeps = 2L)))
    ),
    "does not accept a diagnostic control"
  )
  public_formals <- names(formals(multiFSYNC::bayesSYNC_multi_pre_score))
  testthat::expect_false(any(
    c("dense_gate_sweeps", "random_scale_calibration", "truth", "true_params") %in%
      public_formals
  ))
})
