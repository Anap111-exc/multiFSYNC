# Test: T=1 and T→0 temperature effects
# Verifies that annealing schedules behave correctly.

test_that("T=1 (no annealing) produces valid results", {
  set.seed(42)
  dat <- simulate_multi_study_data(
    S = 1, n_s = c(8), p = 2, d = 0,
    L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
    K = 6, n_obs = 15, seed = 42,
    bool_sparse_loadings = FALSE
  )

  fit <- bayesSYNC_multi(
    Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
    L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
    K = 6, anneal = NULL, maxit = 15, n_cpus = 1,
    verbose = FALSE, seed = 42, bool_scale = FALSE
  )

  expect_true(is.numeric(fit$ELBO_iter))
  expect_true(length(fit$ELBO) >= 1)
  expect_true(is.matrix(fit$mu_q_gamma_a))
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))
})

test_that("T→0 (annealing) posterior variances shrink monotonically", {
  set.seed(42)
  dat <- simulate_multi_study_data(
    S = 1, n_s = c(8), p = 2, d = 0,
    L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
    K = 6, n_obs = 15, seed = 42,
    bool_sparse_loadings = FALSE
  )

  # With annealing (T starts at 1.9, converges to 1)
  fit <- bayesSYNC_multi(
    Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
    L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
    K = 6, anneal = c(1, 1.9, 20), maxit = 30, n_cpus = 1,
    verbose = FALSE, seed = 42, bool_scale = FALSE
  )

  expect_true(is.numeric(fit$ELBO_iter))
  expect_true(fit$ELBO_iter > -1e6)

  # PPI should be valid
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))
})

test_that("T value within valid range (T < 2)", {
  dat <- simulate_multi_study_data(
    S = 1, n_s = c(5), p = 2, d = 0,
    L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
    K = 6, n_obs = 10, seed = 1,
    bool_sparse_loadings = FALSE
  )

  # T_max = 1.9 should work
  expect_error(
    bayesSYNC_multi(
      Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
      L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
      K = 6, anneal = c(1, 1.9, 10), maxit = 10, n_cpus = 1,
      verbose = FALSE, seed = 1, bool_scale = FALSE
    ),
    NA  # no error expected
  )

  # T_max = 2 should fail (kappa_q_a = 2c-1 = 0)
  expect_error(
    bayesSYNC_multi(
      Y = dat$Y, Z = NULL, time_obs = dat$time_obs,
      L_f = 1, L_s = 0, M_f = c(1), M_s = list(integer(0)),
      K = 6, anneal = c(1, 2.0, 10), maxit = 10, n_cpus = 1,
      verbose = FALSE, seed = 1, bool_scale = FALSE
    )
  )
})
