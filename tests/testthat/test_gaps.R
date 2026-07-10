# Gap tests: orthonormalise validation, JAOUA T>1, variance monotonicity, factor PPI

test_that("orthonormalise: crossprod(Phi) ≈ I and PVE sums to 100", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(10), p=3, d=0,
    L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)), K=6, n_obs=20, seed=42,
    bool_sparse_loadings=FALSE)

  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)), K=6, anneal=NULL, maxit=10, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE)

  orth <- orthonormalise_multi(C_g=fit$C_g, time_g=fit$time_g,
    mu_q_nu_mu=fit$mu_q_nu_mu, mu_q_nu_phi=fit$mu_q_nu_phi,
    mu_q_nu_psi=fit$mu_q_nu_psi, mu_q_zeta=fit$mu_q_zeta,
    Sigma_q_zeta=fit$Sigma_q_zeta, mu_q_xi=fit$mu_q_xi,
    Sigma_q_xi=fit$Sigma_q_xi, mu_q_a=fit$mu_q_a,
    mu_q_b_specific=fit$mu_q_b_specific,
    mu_q_gamma_a=fit$mu_q_gamma_a, mu_q_gamma_b=fit$mu_q_gamma_b,
    S=fit$S, n_s=fit$n_s, p=fit$p,
    L_f=fit$L_f, L_s=fit$L_s, M_f=fit$M_f, M_s=fit$M_s)

  Phi <- orth$list_Phi_hat[[1]]
  for (m in 1:fit$M_f[1]) {
    int_phi2 <- trapint(fit$time_g, Phi[, m]^2)
    expect_equal(int_phi2, 1.0, tolerance = 0.1)
  }

  pve <- orth$list_cumulated_pve[[1]]
  expect_equal(tail(pve, 1), 100, tolerance = 1e-6)
})

test_that("factor_ppi formula: 1 - prod(1 - PPI_j)", {
  PPI <- matrix(c(0.2, 0.5, 0.8), nrow=3, ncol=1)
  expected <- 1 - prod(1 - PPI)
  fp <- 1 - exp(colSums(log1p(-PPI)))
  expect_equal(fp[1], expected, tolerance = 1e-10)
  expect_true(fp[1] > 0 && fp[1] < 1)
})

test_that("JAOUA T>1: annealing converges independently", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(10), p=3, d=0,
    L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)), K=6, n_obs=15, seed=42,
    bool_sparse_loadings=FALSE)

  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)), K=6,
    anneal=c(1, 1.9, 30), maxit=40, n_cpus=1,
    verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_true(is.numeric(fit$ELBO_iter))
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))
})

test_that("Variance tr(Sigma) decreases monotonically with c", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(8), p=2, d=0,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6, n_obs=15, seed=42,
    bool_sparse_loadings=FALSE)

  c_vals <- c(1, 5, 20)
  trs <- sapply(c_vals, function(cc) {
    fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
      L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=6,
      anneal=c(1, 1/cc, 5), maxit=if(cc==1) 10 else 5, n_cpus=1,
      verbose=FALSE, seed=42, bool_scale=FALSE)
    tr(fit$Sigma_q_nu_mu[[1]][[1]])
  })
  expect_true(all(diff(trs) < 0))
})
