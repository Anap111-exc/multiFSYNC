testthat::test_that("O'Sullivan ELBO keeps linear and nonlinear Gaussian blocks", {
  K <- 2L
  D <- K + 2L
  mean_q <- c(0.4, -0.3, 0.2, 0.7)
  Sigma_q <- matrix(
    c(
      0.70, 0.08, 0.03, 0.00,
      0.08, 0.50, 0.02, 0.04,
      0.03, 0.02, 0.40, 0.06,
      0.00, 0.04, 0.06, 0.60
    ),
    nrow = D, byrow = TRUE
  )
  Sigma_0 <- matrix(c(1.8, 0.25, 0.25, 1.1), nrow = 2)
  inv_Sigma_0 <- solve(Sigma_0)
  logdet_Sigma_0 <- as.numeric(
    determinant(Sigma_0, logarithm = TRUE)$modulus
  )

  shape_sigsq <- 2.6
  rate_sigsq <- 1.4
  expected_recip_sigsq <- shape_sigsq / rate_sigsq
  expected_log_sigsq <- log(rate_sigsq) - digamma(shape_sigsq)

  alpha_index <- 1:2
  u_index <- 3:4
  alpha_second <- drop(
    crossprod(mean_q[alpha_index],
              inv_Sigma_0 %*% mean_q[alpha_index])
  ) + sum(diag(
    inv_Sigma_0 %*% Sigma_q[alpha_index, alpha_index, drop = FALSE]
  ))
  u_second <- sum(mean_q[u_index]^2) +
    sum(diag(Sigma_q[u_index, u_index, drop = FALSE]))

  expected_log_prior_alpha <- -0.5 * (
    2 * log(2 * pi) + logdet_Sigma_0 + alpha_second
  )
  expected_log_prior_u <- -0.5 * (
    K * log(2 * pi) +
      K * expected_log_sigsq +
      expected_recip_sigsq * u_second
  )
  entropy_q <- 0.5 * (
    D * (1 + log(2 * pi)) +
      as.numeric(determinant(Sigma_q, logarithm = TRUE)$modulus)
  )
  expected <- expected_log_prior_alpha +
    expected_log_prior_u + entropy_q

  actual <- multiFSYNC:::.elbo_osullivan_block(
    mean = mean_q,
    covariance = Sigma_q,
    expected_recip_sigsq = expected_recip_sigsq,
    expected_log_sigsq = expected_log_sigsq,
    inv_Sigma_0 = inv_Sigma_0,
    logdet_Sigma_0 = logdet_Sigma_0,
    K = K
  )

  testthat::expect_equal(actual, expected, tolerance = 1e-12)
})

testthat::test_that("standard-normal score ELBO includes all Gaussian constants", {
  mean_q <- c(0.3, -0.5, 0.2)
  Sigma_q <- matrix(
    c(
      0.8, 0.1, 0.0,
      0.1, 0.6, 0.05,
      0.0, 0.05, 0.9
    ),
    nrow = 3, byrow = TRUE
  )
  M <- length(mean_q)

  expected_log_prior <- -0.5 * (
    M * log(2 * pi) + sum(mean_q^2) + sum(diag(Sigma_q))
  )
  entropy_q <- 0.5 * (
    M * (1 + log(2 * pi)) +
      as.numeric(determinant(Sigma_q, logarithm = TRUE)$modulus)
  )

  actual <- multiFSYNC:::.elbo_standard_normal_score(mean_q, Sigma_q)
  testthat::expect_equal(
    actual, expected_log_prior + entropy_q, tolerance = 1e-12
  )
  testthat::expect_equal(
    multiFSYNC:::.elbo_standard_normal_score(rep(0, M), diag(M)),
    0,
    tolerance = 1e-14
  )
})

testthat::test_that("Half-Cauchy ELBO contains both IG priors and entropies", {
  shape_sigsq <- 2.3
  rate_sigsq <- 1.7
  shape_aux <- 0.9
  rate_aux <- 0.65
  A <- 4

  elog_sigsq <- log(rate_sigsq) - digamma(shape_sigsq)
  erecip_sigsq <- shape_sigsq / rate_sigsq
  elog_aux <- log(rate_aux) - digamma(shape_aux)
  erecip_aux <- shape_aux / rate_aux

  log_prior_sigsq_given_aux <- -0.5 * elog_aux - lgamma(0.5) -
    1.5 * elog_sigsq - erecip_aux * erecip_sigsq
  log_prior_aux <- 0.5 * log(1 / A^2) - lgamma(0.5) -
    1.5 * elog_aux - erecip_aux / A^2
  entropy_ig <- function(shape, rate) {
    log(rate) + lgamma(shape) -
      (shape + 1) * digamma(shape) + shape
  }
  expected <- log_prior_sigsq_given_aux + log_prior_aux +
    entropy_ig(shape_sigsq, rate_sigsq) +
    entropy_ig(shape_aux, rate_aux)

  actual <- multiFSYNC:::.elbo_half_cauchy_pair(
    shape_sigsq = shape_sigsq,
    rate_sigsq = rate_sigsq,
    shape_aux = shape_aux,
    rate_aux = rate_aux,
    A = A
  )

  testthat::expect_equal(actual, expected, tolerance = 1e-12)
  testthat::expect_false(isTRUE(all.equal(
    actual,
    log_prior_sigsq_given_aux + log_prior_aux +
      entropy_ig(shape_sigsq, rate_sigsq)
  )))
})

testthat::test_that("Beta prior-minus-variational term uses positive entropy", {
  shape1_q <- 2.4
  shape2_q <- 3.7
  shape1_prior <- 0.8
  shape2_prior <- 5.2
  elog_omega <- digamma(shape1_q) -
    digamma(shape1_q + shape2_q)
  elog_one_minus <- digamma(shape2_q) -
    digamma(shape1_q + shape2_q)

  expected_log_prior <- -lbeta(shape1_prior, shape2_prior) +
    (shape1_prior - 1) * elog_omega +
    (shape2_prior - 1) * elog_one_minus
  beta_entropy <- lbeta(shape1_q, shape2_q) -
    (shape1_q - 1) * elog_omega -
    (shape2_q - 1) * elog_one_minus
  expected <- expected_log_prior + beta_entropy

  actual <- multiFSYNC:::.elbo_beta_prior_entropy(
    shape1_q, shape2_q, shape1_prior, shape2_prior
  )

  testthat::expect_equal(actual, expected, tolerance = 1e-12)
  testthat::expect_gt(
    abs(actual - (expected_log_prior - beta_entropy)),
    1e-3
  )
  testthat::expect_equal(
    multiFSYNC:::.elbo_beta_prior_entropy(
      shape1_q, shape2_q, shape1_q, shape2_q
    ),
    0,
    tolerance = 1e-14
  )
})

testthat::test_that("spike-and-slab ELBO includes slab and indicator entropies", {
  ppi <- 0.3
  slab_mean <- -0.7
  slab_variance <- 0.4
  omega_shape1 <- 1.8
  omega_shape2 <- 3.1
  expected_log_omega <- digamma(omega_shape1) -
    digamma(omega_shape1 + omega_shape2)
  expected_log_one_minus <- digamma(omega_shape2) -
    digamma(omega_shape1 + omega_shape2)

  expected_log_slab_prior <- ppi * -0.5 * (
    log(2 * pi) + slab_mean^2 + slab_variance
  )
  slab_entropy <- ppi * 0.5 * (
    log(2 * pi * slab_variance) + 1
  )
  expected_log_indicator_prior <- ppi * expected_log_omega +
    (1 - ppi) * expected_log_one_minus
  indicator_entropy <- -ppi * log(ppi) -
    (1 - ppi) * log(1 - ppi)
  expected <- expected_log_slab_prior + slab_entropy +
    expected_log_indicator_prior + indicator_entropy

  actual <- multiFSYNC:::.elbo_spike_slab_scalar(
    ppi, slab_mean, slab_variance,
    expected_log_omega, expected_log_one_minus
  )

  testthat::expect_equal(actual, expected, tolerance = 1e-12)
  testthat::expect_equal(
    multiFSYNC:::.elbo_spike_slab_scalar(
      0, slab_mean, slab_variance,
      expected_log_omega, expected_log_one_minus
    ),
    expected_log_one_minus,
    tolerance = 1e-14
  )
  testthat::expect_equal(
    multiFSYNC:::.elbo_spike_slab_scalar(
      1, slab_mean, slab_variance,
      expected_log_omega, expected_log_one_minus
    ),
    0.5 * (log(slab_variance) + 1 -
             slab_mean^2 - slab_variance) + expected_log_omega,
    tolerance = 1e-14
  )
})

testthat::test_that("ELBO history diagnostics count only material decreases", {
  history <- c(10, 12, 12 - 1e-9, 11, 13)
  diagnostics <- multiFSYNC:::.elbo_history_diagnostics(
    history, tolerance = 1e-8
  )

  testthat::expect_equal(diagnostics$differences, diff(history))
  testthat::expect_equal(diagnostics$decrease_count, 1)
  testthat::expect_equal(
    diagnostics$minimum_difference, min(diff(history))
  )
  testthat::expect_equal(diagnostics$last_difference, 2)
  testthat::expect_false(diagnostics$monotone_within_tolerance)

  tolerated <- multiFSYNC:::.elbo_history_diagnostics(
    c(1, 1 - 1e-9), tolerance = 1e-8
  )
  testthat::expect_equal(tolerated$decrease_count, 0)
  testthat::expect_true(tolerated$monotone_within_tolerance)
})

testthat::test_that("adaptive Cholesky uses jitter only after failure", {
  multiFSYNC:::.reset_spd_diagnostics()

  precision <- matrix(c(2.0, 0.3, 0.3, 1.5), nrow = 2)
  inverse_no_jitter <- multiFSYNC:::.inverse_spd(
    precision, context = "component-test-spd"
  )
  diagnostic_no_jitter <- attr(
    inverse_no_jitter, "solver_diagnostic", exact = TRUE
  )

  testthat::expect_false(diagnostic_no_jitter$used_jitter)
  testthat::expect_equal(diagnostic_no_jitter$jitter, 0)
  testthat::expect_equal(diagnostic_no_jitter$attempts, 1L)
  testthat::expect_equal(
    as.numeric(inverse_no_jitter),
    as.numeric(solve(precision)),
    tolerance = 1e-12
  )

  singular_precision <- diag(c(1, 0))
  inverse_with_jitter <- multiFSYNC:::.inverse_spd(
    singular_precision, context = "component-test-singular"
  )
  diagnostic_with_jitter <- attr(
    inverse_with_jitter, "solver_diagnostic", exact = TRUE
  )
  jittered_precision <- singular_precision +
    diagnostic_with_jitter$jitter * diag(2)

  testthat::expect_true(diagnostic_with_jitter$used_jitter)
  testthat::expect_gt(diagnostic_with_jitter$jitter, 0)
  testthat::expect_equal(diagnostic_with_jitter$attempts, 2L)
  testthat::expect_equal(
    jittered_precision %*% inverse_with_jitter,
    diag(2),
    tolerance = 1e-10
  )

  aggregate_diagnostics <- multiFSYNC:::.get_spd_diagnostics()
  testthat::expect_equal(aggregate_diagnostics$total_calls, 2L)
  testthat::expect_equal(aggregate_diagnostics$jitter_count, 1L)
  testthat::expect_true(aggregate_diagnostics$used_jitter)
  testthat::expect_equal(
    aggregate_diagnostics$events$context,
    "component-test-singular"
  )

  multiFSYNC:::.reset_spd_diagnostics()
})

testthat::test_that("error-variance update returns the direct current RSS cache", {
  S <- 1L
  n_s <- 2L
  p <- 2L
  L_f <- 0L
  L_s <- 0L

  C <- list(list(
    cbind(1, c(0, 0.5, 1), c(-0.2, 0.1, 0.4)),
    cbind(1, c(0.25, 0.75), c(0.3, -0.1))
  ))
  list_cp_C <- lapply(C, function(study) {
    lapply(study, crossprod)
  })
  Y <- list(list(
    list(c(0.2, 0.1, -0.1), c(-0.4, -0.2, 0.3)),
    list(c(0.35, -0.05), c(-0.15, 0.25))
  ))

  mu_q_nu_mu <- list(list(
    c(0.1, -0.2, 0.05),
    c(-0.3, 0.25, 0.1)
  ))
  Sigma_q_nu_mu <- list(list(
    diag(c(0.04, 0.03, 0.02)),
    matrix(
      c(
        0.05, 0.01, 0.00,
        0.01, 0.04, 0.005,
        0.00, 0.005, 0.03
      ),
      nrow = 3, byrow = TRUE
    )
  ))

  empty_shared <- matrix(numeric(), nrow = p, ncol = 0L)
  empty_specific <- list(matrix(numeric(), nrow = p, ncol = 0L))
  mu_q_recip_a_eps <- matrix(c(0.8, 1.1), nrow = S)
  total_obs_sj <- matrix(c(5, 5), nrow = S)

  rss_cache <- multiFSYNC:::compute_rss_cache(
    Y = Y, C = C, list_cp_C = list_cp_C,
    mu_q_nu_mu = mu_q_nu_mu,
    Sigma_q_nu_mu = Sigma_q_nu_mu,
    mu_q_nu_beta = NULL, Sigma_q_nu_beta = NULL, Z = NULL,
    mu_q_zeta = NULL, Sigma_q_zeta = NULL,
    mu_q_nu_phi = list(), Sigma_q_nu_phi = list(),
    mu_q_xi = NULL, Sigma_q_xi = NULL,
    mu_q_nu_psi = list(), Sigma_q_nu_psi = list(),
    mu_q_a = empty_shared, term_a = empty_shared,
    mu_q_b_specific = empty_specific,
    term_b_specific = empty_specific,
    S = S, n_s = n_s, p = p, L_f = L_f, L_s = L_s
  )

  updated <- multiFSYNC:::update_sigsq_eps(
    Y = Y, C = C, list_cp_C = list_cp_C,
    mu_q_nu_mu = mu_q_nu_mu,
    Sigma_q_nu_mu = Sigma_q_nu_mu,
    mu_q_nu_beta = NULL, Sigma_q_nu_beta = NULL, Z = NULL,
    mu_q_zeta = NULL, Sigma_q_zeta = NULL,
    mu_q_nu_phi = list(), Sigma_q_nu_phi = list(),
    mu_q_xi = NULL, Sigma_q_xi = NULL,
    mu_q_nu_psi = list(), Sigma_q_nu_psi = list(),
    mu_q_a = empty_shared, term_a = empty_shared,
    mu_q_b_specific = empty_specific,
    term_b_specific = empty_specific,
    mu_q_recip_a_eps = mu_q_recip_a_eps,
    S = S, n_s = n_s, p = p, L_f = L_f, L_s = L_s,
    total_obs_sj = total_obs_sj,
    c_val = 1
  )

  direct_rss <- matrix(NA_real_, nrow = n_s, ncol = p)
  manual_rss <- matrix(NA_real_, nrow = n_s, ncol = p)
  for (i in seq_len(n_s)) {
    for (j in seq_len(p)) {
      direct_rss[i, j] <- multiFSYNC:::compute_rss_single(
        s = 1L, i = i, j = j,
        Y = Y, C = C, list_cp_C = list_cp_C,
        mu_q_nu_mu = mu_q_nu_mu,
        Sigma_q_nu_mu = Sigma_q_nu_mu,
        mu_q_nu_beta = NULL, Sigma_q_nu_beta = NULL, Z = NULL,
        mu_q_zeta = NULL, Sigma_q_zeta = NULL,
        mu_q_nu_phi = list(), Sigma_q_nu_phi = list(),
        mu_q_xi = NULL, Sigma_q_xi = NULL,
        mu_q_nu_psi = list(), Sigma_q_nu_psi = list(),
        mu_q_a = empty_shared, term_a = empty_shared,
        mu_q_b_specific = empty_specific,
        term_b_specific = empty_specific,
        L_f = L_f, L_s = L_s
      )
      residual <- Y[[1]][[i]][[j]] -
        as.vector(C[[1]][[i]] %*% mu_q_nu_mu[[1]][[j]])
      manual_rss[i, j] <- sum(residual^2) + sum(diag(
        list_cp_C[[1]][[i]] %*% Sigma_q_nu_mu[[1]][[j]]
      ))
    }
  }

  testthat::expect_equal(direct_rss, manual_rss, tolerance = 1e-12)
  testthat::expect_equal(
    rss_cache$expected_rss[[1]], direct_rss, tolerance = 1e-12
  )
  testthat::expect_equal(
    updated$expected_rss[[1]], direct_rss, tolerance = 1e-12
  )
  testthat::expect_equal(
    updated$expected_rss_sum,
    matrix(colSums(direct_rss), nrow = S),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    updated$expected_rss_sum,
    rss_cache$expected_rss_sum,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    updated$lambda_q_sigsq_eps,
    mu_q_recip_a_eps + 0.5 * rss_cache$expected_rss_sum,
    tolerance = 1e-12
  )
})
