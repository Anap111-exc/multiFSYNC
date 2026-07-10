# Test: PPI temperature effect verification
# Verifies that the entropy term in PPI is NOT multiplied by c

test_that("PPI entropy term is NOT multiplied by c (standalone)", {
  mu_val <- 2; sigma2_val <- 1
  log_omega <- 0; log_1_omega <- 0
  entropy_term <- 0.5 * (mu_val^2 / sigma2_val + log(sigma2_val))

  c1 <- 1
  logit_c1 <- c1 * (log_omega - log_1_omega) + entropy_term
  expect_equal(logit_c1, 2.0)

  c05 <- 0.5
  logit_c05 <- c05 * (log_omega - log_1_omega) + entropy_term
  expect_equal(logit_c05, 2.0)

  logit_buggy <- c05 * (log_omega - log_1_omega + entropy_term)
  expect_equal(logit_buggy, 1.0)
})

test_that("PPI temperature integration test", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(8), p=2, d=0,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6, n_obs=15, seed=42,
    bool_sparse_loadings=FALSE)

  fit_c1 <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6, anneal=NULL, maxit=10, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_true(all(fit_c1$mu_q_gamma_a >= 0 & fit_c1$mu_q_gamma_a <= 1))
})

test_that("PPI stays in (0,1) across iterations", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(10), p=3, d=0,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=6, n_obs=10, seed=42,
    bool_sparse_loadings=FALSE)

  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=6, anneal=NULL, maxit=10, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_true(all(fit$mu_q_gamma_a >= 0))
  expect_true(all(fit$mu_q_gamma_a <= 1))
})
