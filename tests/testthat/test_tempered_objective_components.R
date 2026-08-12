testthat::test_that("tempered Gaussian blocks equal energy plus T entropy", {
  temperature <- 1.7
  K <- 2L
  D <- K + 2L
  mean_q <- c(0.3, -0.2, 0.4, -0.1)
  Sigma_q <- diag(c(0.7, 1.1, 0.5, 0.9))
  inv_Sigma0 <- matrix(c(1.4, 0.2, 0.2, 0.8), 2, 2)
  logdet_Sigma0 <- -multiFSYNC:::.elbo_logdet_positive(inv_Sigma0)
  expected_precision <- 1.3
  expected_log_variance <- -0.25

  alpha <- mean_q[1:2]
  u <- mean_q[3:4]
  alpha_second <- drop(crossprod(alpha, inv_Sigma0 %*% alpha)) +
    sum(diag(inv_Sigma0 %*% Sigma_q[1:2, 1:2, drop = FALSE]))
  u_second <- drop(crossprod(u)) + sum(diag(Sigma_q[3:4, 3:4]))
  expected_log_prior <- -0.5 * (
    D * log(2 * pi) + logdet_Sigma0 + alpha_second +
      K * expected_log_variance + expected_precision * u_second
  )
  entropy_q <- 0.5 * (
    D * (1 + log(2 * pi)) +
      multiFSYNC:::.elbo_logdet_positive(Sigma_q)
  )

  actual <- multiFSYNC:::.elbo_osullivan_block(
    mean_q, Sigma_q, expected_precision, expected_log_variance,
    inv_Sigma0, logdet_Sigma0, K, temperature = temperature
  )
  testthat::expect_equal(
    actual, expected_log_prior + temperature * entropy_q,
    tolerance = 1e-12
  )

  score_mean <- c(0.2, -0.5)
  score_covariance <- diag(c(0.6, 1.2))
  score_energy <- -0.5 * (
    2 * log(2 * pi) + sum(score_mean^2) +
      sum(diag(score_covariance))
  )
  score_entropy <- 0.5 * (
    2 * (1 + log(2 * pi)) +
      multiFSYNC:::.elbo_logdet_positive(score_covariance)
  )
  testthat::expect_equal(
    multiFSYNC:::.elbo_standard_normal_score(
      score_mean, score_covariance, temperature = temperature
    ),
    score_energy + temperature * score_entropy,
    tolerance = 1e-12
  )
})

testthat::test_that("tempered spike-and-slab block uses the mixed-measure entropy", {
  temperature <- 1.8
  ppi <- 0.37
  slab_mean <- -0.45
  slab_variance <- 0.72
  elog_omega <- -0.6
  elog_one_minus <- -0.9

  expected_log_prior <-
    ppi * (elog_omega - 0.5 * (
      log(2 * pi) + slab_mean^2 + slab_variance
    )) +
    (1 - ppi) * elog_one_minus
  entropy_q <-
    -multiFSYNC:::.elbo_xlogx(ppi) -
    multiFSYNC:::.elbo_xlogx(1 - ppi) +
    0.5 * ppi * (log(2 * pi) + 1 + log(slab_variance))

  actual <- multiFSYNC:::.elbo_spike_slab_scalar(
    ppi, slab_mean, slab_variance,
    elog_omega, elog_one_minus,
    temperature = temperature
  )
  testthat::expect_equal(
    actual, expected_log_prior + temperature * entropy_q,
    tolerance = 1e-12
  )

  included_energy <- elog_omega - 0.5 * (
    log(2 * pi) + slab_mean^2 + slab_variance
  )
  excluded_energy <- elog_one_minus
  slab_entropy <- 0.5 * (
    log(2 * pi) + 1 + log(slab_variance)
  )
  expected_ppi <- stats::plogis(
    (included_energy - excluded_energy) / temperature + slab_entropy
  )
  objective <- function(probability) {
    multiFSYNC:::.elbo_spike_slab_scalar(
      probability, slab_mean, slab_variance,
      elog_omega, elog_one_minus,
      temperature = temperature
    )
  }
  numerical_ppi <- stats::optimize(
    objective, interval = c(1e-8, 1 - 1e-8), maximum = TRUE,
    tol = 1e-12
  )$maximum
  testthat::expect_equal(numerical_ppi, expected_ppi, tolerance = 1e-7)
})

testthat::test_that("tempered Beta and Half-Cauchy blocks equal energy plus T entropy", {
  temperature <- 1.6
  aq <- 2.7
  bq <- 4.1
  a0 <- 1.3
  b0 <- 2.2
  elog_omega <- digamma(aq) - digamma(aq + bq)
  elog_one_minus <- digamma(bq) - digamma(aq + bq)
  beta_energy <-
    (a0 - 1) * elog_omega + (b0 - 1) * elog_one_minus -
    lbeta(a0, b0)
  beta_entropy <-
    lbeta(aq, bq) - (aq - 1) * elog_omega -
    (bq - 1) * elog_one_minus
  testthat::expect_equal(
    multiFSYNC:::.elbo_beta_prior_entropy(
      aq, bq, a0, b0, temperature = temperature
    ),
    beta_energy + temperature * beta_entropy,
    tolerance = 1e-12
  )

  shape_sigsq <- 3.4
  rate_sigsq <- 1.7
  shape_aux <- 0.8
  rate_aux <- 2.1
  A <- 5
  elog_sigsq <- log(rate_sigsq) - digamma(shape_sigsq)
  erecip_sigsq <- shape_sigsq / rate_sigsq
  elog_aux <- log(rate_aux) - digamma(shape_aux)
  erecip_aux <- shape_aux / rate_aux
  hc_energy <-
    -0.5 * elog_aux - lgamma(0.5) -
    1.5 * elog_sigsq - erecip_aux * erecip_sigsq +
    0.5 * log(1 / A^2) - lgamma(0.5) -
    1.5 * elog_aux - erecip_aux / A^2
  hc_entropy <-
    multiFSYNC:::.elbo_ig_entropy(shape_sigsq, rate_sigsq) +
    multiFSYNC:::.elbo_ig_entropy(shape_aux, rate_aux)
  testthat::expect_equal(
    multiFSYNC:::.elbo_half_cauchy_pair(
      shape_sigsq, rate_sigsq, shape_aux, rate_aux, A,
      temperature = temperature
    ),
    hc_energy + temperature * hc_entropy,
    tolerance = 1e-12
  )
})

testthat::test_that("all tempered helpers reduce exactly to ordinary ELBO at T1", {
  testthat::expect_identical(
    multiFSYNC:::.elbo_standard_normal_score(c(0.2), matrix(0.7, 1, 1)),
    multiFSYNC:::.elbo_standard_normal_score(
      c(0.2), matrix(0.7, 1, 1), temperature = 1
    )
  )
  testthat::expect_identical(
    multiFSYNC:::.elbo_spike_slab_scalar(0.4, 0.2, 0.8, -0.7, -0.9),
    multiFSYNC:::.elbo_spike_slab_scalar(
      0.4, 0.2, 0.8, -0.7, -0.9, temperature = 1
    )
  )
  testthat::expect_identical(
    multiFSYNC:::.elbo_beta_prior_entropy(2.1, 3.2, 1.4, 2.3),
    multiFSYNC:::.elbo_beta_prior_entropy(
      2.1, 3.2, 1.4, 2.3, temperature = 1
    )
  )
  testthat::expect_identical(
    multiFSYNC:::.elbo_half_cauchy_pair(2.3, 1.1, 0.9, 1.4, 3),
    multiFSYNC:::.elbo_half_cauchy_pair(
      2.3, 1.1, 0.9, 1.4, 3, temperature = 1
    )
  )
})
