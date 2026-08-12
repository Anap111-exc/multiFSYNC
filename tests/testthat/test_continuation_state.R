continuation_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    data <- simulate_multi_study_osullivan(
      S = 2L, n_s = c(3L, 3L), p = 3L, d = 0L,
      L_f = 1L, L_s = 1L,
      M_f = 1L, M_s = list(1L, 1L),
      K = 3L, n_obs = 7L, common_grid = TRUE,
      sigma_eps = 0.10, seed = 99101L, use_explicit = FALSE
    )
    fit <- suppressWarnings(bayesSYNC_multi(
      Y = data$Y, Z = NULL, time_obs = data$time_obs,
      L_f = 2L, L_s = 1L, M_f = c(1L, 1L),
      M_s = list(1L, 1L),
      K = 3L, anneal = NULL, n_g = 21L,
      tol_abs = 0, tol_rel = 0, maxit = 3L,
      n_cpus = 1L, verbose = FALSE, seed = 99102L,
      bool_scale = FALSE, bool_var_spec_prob = FALSE, d_0 = 3L,
      convergence_rule = "parameters", lambda_orth = 0,
      initialization = "random"
    ))
    cache <<- list(data = data, fit = fit)
    cache
  }
})

testthat::test_that("reduced continuation state maps reported shared columns", {
  fixture <- continuation_fixture()
  fit <- fixture$fit
  state <- multiFSYNC:::.make_reduced_continuation_state(fit, 1L)

  testthat::expect_s3_class(state, "multiFSYNC_continuation_state")
  testthat::expect_identical(state$L_f, 1L)
  testthat::expect_identical(
    state$diagnostics$shared_working_indices,
    as.integer(fit$factor_order_shared[[1L]])
  )
  testthat::expect_equal(dim(state$parameters$mu_q_a), c(fit$p, 1L))
  testthat::expect_length(state$parameters$mu_q_nu_phi, 1L)
  testthat::expect_length(state$parameters$mu_q_zeta[[1L]], 1L)
  testthat::expect_equal(
    state$parameters$mu_q_a[, 1L],
    fit$mu_q_a[, fit$factor_order_shared[[1L]]]
  )
})

testthat::test_that("zero-shared continuation state has valid empty shapes", {
  fit <- continuation_fixture()$fit
  state <- multiFSYNC:::.make_reduced_continuation_state(fit, integer())

  testthat::expect_identical(state$L_f, 0L)
  testthat::expect_identical(state$M_f, integer())
  testthat::expect_equal(dim(state$parameters$mu_q_a), c(fit$p, 0L))
  testthat::expect_length(state$parameters$mu_q_nu_phi, 0L)
  testthat::expect_length(state$parameters$mu_q_zeta[[1L]], 0L)
  testthat::expect_identical(state$L_s, fit$L_s)
})

testthat::test_that("T1 reduced continuation performs complete finite CAVI sweeps", {
  fixture <- continuation_fixture()
  state <- multiFSYNC:::.make_reduced_continuation_state(fixture$fit, 1L)
  fit <- suppressWarnings(bayesSYNC_multi(
    Y = fixture$data$Y, Z = NULL, time_obs = fixture$data$time_obs,
    L_f = state$L_f, L_s = state$L_s,
    M_f = state$M_f, M_s = state$M_s,
    K = state$K, anneal = NULL, n_g = length(state$time_g),
    time_g = state$time_g,
    tol_abs = 0, tol_rel = 0, maxit = 3L,
    n_cpus = 1L, verbose = FALSE, seed = 99103L,
    bool_scale = FALSE,
    bool_var_spec_prob = state$bool_var_spec_prob, d_0 = 3L,
    convergence_rule = "parameters", lambda_orth = 0,
    initialization = "random", continuation_state = state
  ))

  testthat::expect_true(fit$continuation_diagnostics$used)
  testthat::expect_identical(
    fit$initialization_diagnostics$used_method, "continuation_state"
  )
  testthat::expect_identical(fit$L_f, 1L)
  testthat::expect_length(fit$ELBO, 3L)
  testthat::expect_true(all(vapply(
    fit$ELBO_finite, function(x) isTRUE(x$all), logical(1)
  )))
  testthat::expect_equal(fit$ELBO_decrease_count, 0L)
})

testthat::test_that("continuation validation rejects altered state or objective", {
  fixture <- continuation_fixture()
  state <- multiFSYNC:::.make_reduced_continuation_state(fixture$fit, 1L)
  altered <- state
  altered$time_g[[1L]] <- altered$time_g[[1L]] + 0.01

  call_fit <- function(candidate_state, anneal = NULL) {
    bayesSYNC_multi(
      Y = fixture$data$Y, Z = NULL, time_obs = fixture$data$time_obs,
      L_f = candidate_state$L_f, L_s = candidate_state$L_s,
      M_f = candidate_state$M_f, M_s = candidate_state$M_s,
      K = candidate_state$K, anneal = anneal,
      n_g = length(fixture$fit$time_g), time_g = fixture$fit$time_g,
      maxit = 2L, n_cpus = 1L, verbose = FALSE,
      bool_scale = FALSE, bool_var_spec_prob = FALSE, d_0 = 3L,
      convergence_rule = "parameters", lambda_orth = 0,
      continuation_state = candidate_state
    )
  }
  testthat::expect_error(
    suppressWarnings(call_fit(altered)),
    "scaling or dense time grid differs"
  )
  testthat::expect_error(
    suppressWarnings(call_fit(state, c(1, 1.2, 2))),
    "anneal = NULL"
  )
})
