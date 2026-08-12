test_that("all generators support study-specific factor counts", {
  L_s <- c(0L, 1L, 2L)
  M_s <- list(integer(0), 1L, c(1L, 2L))
  common <- list(
    S = 3L, n_s = c(2L, 2L, 2L), p = 5L, d = 0L,
    L_f = 1L, L_s = L_s, M_f = 1L, M_s = M_s,
    K = 4L, n_obs = 8L, common_grid = TRUE,
    bool_sparse_loadings = FALSE, identified_loadings = TRUE,
    seed = 9401L
  )

  generated <- list(
    do.call(simulate_multi_study_data, common),
    do.call(simulate_multi_study_osullivan, common),
    do.call(simulate_multi_study_structured, c(common, list(n_dense = 60L)))
  )

  for (dat in generated) {
    truth <- dat$true_params
    expect_identical(truth$L_s_by_study, L_s)
    expect_identical(truth$L_s, L_s)
    specific_functions <- if (!is.null(truth$nu_psi_true)) {
      truth$nu_psi_true
    } else {
      truth$psi_dense_list
    }
    expect_identical(lengths(specific_functions), L_s)
    expect_identical(lengths(truth$xi_true), L_s)
    expect_identical(vapply(truth$b_true, ncol, integer(1)), L_s)
    expect_identical(vapply(truth$gamma_b_true, ncol, integer(1)), L_s)
  }
})

test_that("scalar L_s retains the historical public representation", {
  dat <- simulate_multi_study_data(
    S = 2L, n_s = c(2L, 2L), p = 4L, d = 0L,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L, 1L),
    K = 4L, n_obs = 8L, bool_sparse_loadings = FALSE,
    seed = 9402L
  )

  expect_identical(dat$true_params$L_s, 1L)
  expect_identical(dat$true_params$L_s_by_study, c(1L, 1L))
  expect_identical(vapply(dat$true_params$b_true, ncol, integer(1)),
                   c(1L, 1L))
})

test_that("CAVI and the complete ELBO preserve jagged specific blocks", {
  L_s <- c(0L, 1L, 2L)
  M_s <- list(integer(0), 1L, c(1L, 1L))
  dat <- simulate_multi_study_data(
    S = 3L, n_s = c(2L, 2L, 2L), p = 5L, d = 0L,
    L_f = 1L, L_s = L_s, M_f = 1L, M_s = M_s,
    K = 4L, n_obs = 12L, bool_sparse_loadings = FALSE,
    identified_loadings = TRUE, seed = 9403L
  )

  settings <- list(
    list(variable_omega = FALSE, initialization = "random"),
    list(variable_omega = TRUE, initialization = "random"),
    list(variable_omega = FALSE, initialization = "residual_fpca")
  )
  for (setting_index in seq_along(settings)) {
    variable_omega <- settings[[setting_index]]$variable_omega
    fit <- suppressWarnings(bayesSYNC_multi(
      Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
      L_f = 1L, L_s = L_s, M_f = 1L, M_s = M_s,
      K = 4L, anneal = NULL, n_g = 41L, maxit = 1L,
      convergence_rule = "parameters", bool_scale = FALSE,
      bool_var_spec_prob = variable_omega,
      initialization = settings[[setting_index]]$initialization,
      verbose = FALSE, seed = 9403L + setting_index
    ))

    expect_identical(fit$L_s, L_s)
    expect_identical(fit$L_s_by_study, L_s)
    expect_identical(fit$L_s_max, 2L)
    expect_identical(lengths(fit$mu_q_nu_psi), L_s)
    expect_identical(lengths(fit$mu_q_xi), L_s)
    expect_identical(vapply(fit$mu_q_b_specific, ncol, integer(1)), L_s)
    expect_identical(vapply(fit$mu_q_b_specific_hat, ncol, integer(1)), L_s)
    expect_true(length(fit$ELBO) == 1L && is.finite(fit$ELBO))
    expect_true(isTRUE(fit$elbo_result$finite$all))

    state <- .make_reduced_continuation_state(fit, 1L)
    expect_identical(state$L_s, L_s)
    expect_identical(lengths(state$parameters$mu_q_nu_psi), L_s)
  }
})

test_that("M_s is validated against each study's own count", {
  dat <- simulate_multi_study_data(
    S = 2L, n_s = c(2L, 2L), p = 4L, d = 0L,
    L_f = 1L, L_s = c(0L, 1L), M_f = 1L,
    M_s = list(integer(0), 1L), K = 4L, n_obs = 8L,
    seed = 9405L
  )

  expect_error(
    bayesSYNC_multi(
      Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
      L_f = 1L, L_s = c(0L, 1L), M_f = 1L,
      M_s = list(1L, 1L), K = 4L, anneal = NULL,
      maxit = 1L, verbose = FALSE
    ),
    "must be empty because L_s\\[1\\] = 0"
  )
})
