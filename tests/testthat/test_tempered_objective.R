test_that("Gaussian, inverse-Gamma, and Beta powers match canonical shapes", {
  c_val <- 1 / 1.9

  gaussian_precision <- 3.7
  gaussian_linear <- -0.8
  gaussian_variance <- 1 / (c_val * gaussian_precision)
  gaussian_mean <- c_val * gaussian_variance * gaussian_linear
  expect_equal(gaussian_variance, 1 / (c_val * gaussian_precision))
  expect_equal(gaussian_mean, gaussian_linear / gaussian_precision)

  conditional_ig_shape <- 4.2
  conditional_ig_scale <- 1.6
  tempered_ig_shape <- c_val * (conditional_ig_shape + 1) - 1
  tempered_ig_scale <- c_val * conditional_ig_scale
  x <- c(0.4, 0.9, 2.1)
  powered_ig_kernel <-
    c_val * (-(conditional_ig_shape + 1) * log(x) -
               conditional_ig_scale / x)
  tempered_ig_kernel <-
    -(tempered_ig_shape + 1) * log(x) - tempered_ig_scale / x
  expect_equal(diff(powered_ig_kernel), diff(tempered_ig_kernel))

  conditional_beta_shape1 <- 2.4
  conditional_beta_shape2 <- 5.3
  tempered_beta_shape1 <- 1 + c_val * (conditional_beta_shape1 - 1)
  tempered_beta_shape2 <- 1 + c_val * (conditional_beta_shape2 - 1)
  omega <- c(0.2, 0.45, 0.8)
  powered_beta_kernel <- c_val * (
    (conditional_beta_shape1 - 1) * log(omega) +
      (conditional_beta_shape2 - 1) * log1p(-omega))
  tempered_beta_kernel <-
    (tempered_beta_shape1 - 1) * log(omega) +
      (tempered_beta_shape2 - 1) * log1p(-omega)
  expect_equal(diff(powered_beta_kernel), diff(tempered_beta_kernel))
})

test_that("Half-Cauchy auxiliary and spline variance shapes use c(alpha+1)-1", {
  c_val <- 0.6
  aux <- multiFSYNC:::update_a_eps(
    mu_q_recip_sigsq_eps = matrix(c(0.7, 1.1), nrow = 1),
    A = 10,
    c_val = c_val
  )
  expect_equal(aux$kappa_q_a, 2 * c_val - 1)
  expect_equal(aux$lambda_q_a_eps,
               c_val * (matrix(c(0.7, 1.1), nrow = 1) + 1 / 10^2))

  K <- 3
  expected_spline_shape <- c_val * (K + 3) / 2 - 1
  expect_equal(c_val * (K + 1) / 2 + c_val - 1,
               expected_spline_shape)
})

test_that("both omega hierarchies use the exact tempered Beta update", {
  gamma <- matrix(c(0.2, 0.7, 0.4, 0.9), nrow = 2)
  c_val <- 0.6
  c0 <- 1.3
  d0 <- 2.4

  factor_level <- multiFSYNC:::update_omega_shared(
    gamma, p = 2, c_0 = c0, d_0 = d0,
    c_val = c_val, bool_var_spec_prob = FALSE)
  expect_equal(factor_level$c_1_omega_a,
               1 + c_val * (c0 + colSums(gamma) - 1))
  expect_equal(factor_level$d_1_omega_a,
               1 + c_val * (d0 + 2 - colSums(gamma) - 1))

  variable_level <- multiFSYNC:::update_omega_shared(
    gamma, p = 2, c_0 = c0, d_0 = d0,
    c_val = c_val, bool_var_spec_prob = TRUE)
  expect_equal(variable_level$c_1_omega_a,
               1 + c_val * (c0 + gamma - 1))
  expect_equal(variable_level$d_1_omega_a,
               1 + c_val * (d0 + 1 - gamma - 1))
})

test_that("shared and specific loading updates use the exact tempered logit", {
  c_val <- 0.6
  elog_omega <- -0.2
  elog_one_minus <- -1.0
  common <- list(
    Y = list(list(matrix(0, 1, 2))),
    C = list(list(matrix(c(1, 0), 1, 2))),
    list_cp_C = list(list(diag(c(1, 0)))),
    list_cp_C_Y = list(list(matrix(0, 2, 2))),
    mu_q_nu_mu = list(list(c(0, 0), c(0, 0))),
    mu_q_nu_beta = NULL,
    Z = NULL,
    mu_q_recip_sigsq_eps = matrix(1, 1, 2),
    S = 1L, n_s = 1L, p = 2L, d = 0L,
    K_total = 2L, c_val = c_val, n_cpus = 1L
  )

  shared <- do.call(multiFSYNC:::update_a_loadings, c(common, list(
    mu_q_zeta = list(list(matrix(0, 1, 1))),
    mu_q_nu_phi = list(matrix(c(1, 0), 2, 1)),
    mu_q_xi = NULL,
    mu_q_nu_psi = NULL,
    mu_q_a = matrix(0, 2, 1),
    term_a = matrix(1, 2, 1),
    mu_q_gamma_a = matrix(0.5, 2, 1),
    mu_q_normal_a = matrix(0, 2, 1),
    Sigma_q_normal_a = matrix(1, 2, 1),
    mu_q_b_specific = NULL,
    mu_q_log_omega_a = elog_omega,
    mu_q_log_1_omega_a = elog_one_minus,
    tr_qi_shared = array(1, dim = c(1, 1, 1)),
    L_f = 1L, L_s = 0L
  )))
  expected_shared_logit <- multiFSYNC:::.tempered_spike_slab_logit(
    elog_omega - elog_one_minus,
    shared$mu_q_normal_a[1, 1], shared$Sigma_q_normal_a[1, 1], c_val)
  expect_equal(shared$mu_q_gamma_a[1, 1], plogis(expected_shared_logit))

  specific <- do.call(multiFSYNC:::update_b_loadings, c(common, list(
    mu_q_zeta = list(list()),
    mu_q_nu_phi = list(),
    mu_q_xi = list(list(matrix(0, 1, 1))),
    mu_q_nu_psi = list(list(matrix(c(1, 0), 2, 1))),
    mu_q_a = matrix(numeric(), 2, 0),
    mu_q_b_specific = list(matrix(0, 2, 1)),
    term_b_specific = list(matrix(1, 2, 1)),
    mu_q_gamma_b = list(matrix(0.5, 2, 1)),
    mu_q_normal_b_specific = list(matrix(0, 2, 1)),
    Sigma_q_normal_b = list(matrix(1, 2, 1)),
    mu_q_log_omega_b = list(elog_omega),
    mu_q_log_1_omega_b = list(elog_one_minus),
    tr_xi_specific = array(1, dim = c(1, 1, 1)),
    L_f = 0L, L_s = 1L
  )))
  expected_specific_logit <- multiFSYNC:::.tempered_spike_slab_logit(
    elog_omega - elog_one_minus,
    specific$mu_q_normal_b[[1]][1, 1],
    specific$Sigma_q_normal_b[[1]][1, 1], c_val)
  expect_equal(specific$mu_q_gamma_b[[1]][1, 1],
               plogis(expected_specific_logit))
})
