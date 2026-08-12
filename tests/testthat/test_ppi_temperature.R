# Test: PPI temperature effect verification under the fixed mixed reference
# measure used by the tempered objective.

test_that("tempered PPI logit retains the continuous-slab normaliser", {
  mu_val <- 2; sigma2_val <- 1
  log_omega <- 0; log_1_omega <- 0
  entropy_term <- 0.5 * (mu_val^2 / sigma2_val + log(sigma2_val))

  c1 <- 1
  logit_c1 <- multiFSYNC:::.tempered_spike_slab_logit(
    log_omega - log_1_omega, mu_val, sigma2_val, c1)
  expect_equal(logit_c1, 2.0)

  c05 <- 0.5
  logit_c05 <- multiFSYNC:::.tempered_spike_slab_logit(
    log_omega - log_1_omega, mu_val, sigma2_val, c05)
  expect_equal(logit_c05, entropy_term + 0.25 * log(2 * pi))

  logit_without_measure_term <-
    c05 * (log_omega - log_1_omega) + entropy_term
  expect_equal(logit_c05 - logit_without_measure_term,
               0.25 * log(2 * pi))
})

test_that("tempered PPI logit equals direct Gaussian integration", {
  c_val <- 0.55
  quadratic_precision <- 2.3
  linear_term <- 0.7
  expected_log_prior_odds <- -0.4

  slab_variance <- 1 / (c_val * quadratic_precision)
  slab_mean <- linear_term / quadratic_precision
  direct_integral <- integrate(
    function(a) {
      exp(c_val * (-0.5 * quadratic_precision * a^2 +
                     linear_term * a - 0.5 * log(2 * pi)))
    },
    lower = -Inf, upper = Inf, rel.tol = 1e-13
  )$value
  direct_logit <- c_val * expected_log_prior_odds + log(direct_integral)

  expect_equal(
    multiFSYNC:::.tempered_spike_slab_logit(
      expected_log_prior_odds, slab_mean, slab_variance, c_val),
    direct_logit,
    tolerance = 1e-12
  )
})

test_that("PPI temperature integration test", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(8), p=2, d=0,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6, n_obs=15, seed=42,
    bool_sparse_loadings=FALSE)

  fit_c1 <- suppressWarnings(bayesSYNC_multi(
    Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6,
    anneal=NULL, maxit=10, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE))

  expect_true(all(fit_c1$mu_q_gamma_a >= 0 & fit_c1$mu_q_gamma_a <= 1))
})

test_that("PPI stays in (0,1) across iterations", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(10), p=3, d=0,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=6, n_obs=10, seed=42,
    bool_sparse_loadings=FALSE)

  fit <- suppressWarnings(bayesSYNC_multi(
    Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=6,
    anneal=NULL, maxit=10, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE))

  expect_true(all(fit$mu_q_gamma_a >= 0))
  expect_true(all(fit$mu_q_gamma_a <= 1))
})
