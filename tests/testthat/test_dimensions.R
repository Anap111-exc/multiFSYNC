# Test: dimension checks for all variational parameters

test_that("All parameter dimensions match expected sizes", {
  set.seed(42)
  dat <- simulate_multi_study_data(
    S = 2, n_s = c(5, 4), p = 3, d = 1,
    L_f = 2, L_s = 1, M_f = c(2, 2), M_s = list(c(1), c(1)),
    K = 6, n_obs = 15, seed = 42,
    bool_sparse_loadings = FALSE
  )

  fit <- bayesSYNC_multi(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 2, L_s = 1, M_f = c(2, 2), M_s = list(c(1), c(1)),
    K = 6, anneal = NULL, maxit = 5, n_cpus = 1,
    verbose = FALSE, seed = 42, bool_scale = FALSE
  )

  K_total <- fit$K_total
  S <- fit$S; n_s <- fit$n_s; p <- fit$p; d <- fit$d
  L_f <- fit$L_f; L_s <- fit$L_s; M_f <- fit$M_f; M_s <- fit$M_s

  # nu_mu: (K+2) x 1 per (s,j)
  expect_equal(length(fit$mu_q_nu_mu), S)
  expect_equal(length(fit$mu_q_nu_mu[[1]]), p)
  expect_equal(length(fit$mu_q_nu_mu[[1]][[1]]), K_total)

  # nu_beta: (K+2) x 1 per (j,r)
  expect_equal(length(fit$mu_q_nu_beta), p)
  expect_equal(length(fit$mu_q_nu_beta[[1]]), d)

  # nu_phi: (K+2) x M_f[l] per l
  expect_equal(length(fit$mu_q_nu_phi), L_f)
  expect_equal(dim(fit$mu_q_nu_phi[[1]]), c(K_total, M_f[1]))

  # nu_psi: (K+2) x M_s[[s]][l] per (s,l)
  expect_equal(length(fit$mu_q_nu_psi), S)
  expect_equal(length(fit$mu_q_nu_psi[[1]]), L_s)
  expect_equal(dim(fit$mu_q_nu_psi[[1]][[1]]), c(K_total, M_s[[1]][1]))

  # zeta scores: n_s[s] x M_f[l]
  expect_equal(length(fit$mu_q_zeta), S)
  expect_equal(length(fit$mu_q_zeta[[1]]), L_f)
  expect_equal(dim(fit$mu_q_zeta[[1]][[1]]), c(n_s[1], M_f[1]))

  # xi scores: n_s[s] x M_s[[s]][l]
  expect_equal(length(fit$mu_q_xi), S)
  expect_equal(length(fit$mu_q_xi[[1]]), L_s)
  expect_equal(dim(fit$mu_q_xi[[1]][[1]]), c(n_s[1], M_s[[1]][1]))

  # Shared loadings a: p x L_f
  expect_equal(dim(fit$mu_q_a), c(p, L_f))
  expect_equal(dim(fit$mu_q_gamma_a), c(p, L_f))

  # Specific loadings b: p x L_s per study
  expect_equal(length(fit$mu_q_b_specific), S)
  expect_equal(dim(fit$mu_q_b_specific[[1]]), c(p, L_s))

  # Variance parameters
  expect_equal(dim(fit$mu_q_recip_sigsq_eps), c(S, p))
  expect_equal(length(fit$mu_q_recip_sigsq_phi), L_f)
  expect_equal(length(fit$mu_q_recip_sigsq_phi[[1]]), M_f[1])

  # Precomputation
  expect_equal(length(fit$list_cp_C), S)
  expect_equal(length(fit$list_cp_C[[1]]), n_s[1])
  expect_equal(dim(fit$list_cp_C[[1]][[1]]), c(K_total, K_total))
  expect_equal(dim(fit$sum_list_cp_C[[1]]), c(K_total, K_total))
})
