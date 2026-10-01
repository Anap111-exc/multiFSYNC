test_that("beta update matches direct small-matrix coordinates", {
  set.seed(971901L)
  S <- 2L
  n_s <- c(2L, 3L)
  p <- 3L
  d <- 2L
  K <- 3L
  K_total <- K + 2L
  C <- lapply(n_s, function(n) {
    lapply(seq_len(n), function(i) matrix(rnorm(6L * K_total), 6L))
  })
  Y <- lapply(n_s, function(n) {
    lapply(seq_len(n), function(i) {
      lapply(seq_len(p), function(j) rnorm(6L))
    })
  })
  Z <- lapply(n_s, function(n) matrix(rnorm(n * d), n, d))
  list_cp_C <- lapply(C, function(study) lapply(study, crossprod))
  list_cp_C_Y <- lapply(seq_len(S), function(s) {
    lapply(seq_len(n_s[[s]]), function(i) {
      vapply(seq_len(p), function(j) {
        as.vector(crossprod(C[[s]][[i]], Y[[s]][[i]][[j]]))
      }, numeric(K_total))
    })
  })
  mu_q_nu_mu <- lapply(seq_len(S), function(s) {
    lapply(seq_len(p), function(j) rnorm(K_total))
  })
  mu_q_nu_beta <- lapply(seq_len(p), function(j) {
    lapply(seq_len(d), function(r) rnorm(K_total))
  })
  Sigma_q_nu_beta <- lapply(seq_len(p), function(j) {
    lapply(seq_len(d), function(r) diag(K_total))
  })
  eps_prec <- matrix(runif(S * p, 1, 3), S, p)
  beta_prec <- matrix(runif(p * d, 0.5, 1.5), p, d)
  inv_Sigma_beta <- diag(c(1.4, 2.1))
  temperature <- 1.2

  actual <- update_nu_beta(
    Y = Y, C = C, list_cp_C = list_cp_C,
    list_cp_C_Y = list_cp_C_Y, mu_q_nu_mu = mu_q_nu_mu,
    mu_q_nu_beta = mu_q_nu_beta,
    Sigma_q_nu_beta = Sigma_q_nu_beta, Z = Z,
    mu_q_zeta = NULL, mu_q_nu_phi = NULL,
    mu_q_xi = NULL, mu_q_nu_psi = NULL,
    mu_q_a = NULL, mu_q_b_specific = NULL,
    mu_q_recip_sigsq_eps = eps_prec,
    mu_q_recip_sigsq_beta = beta_prec,
    inv_Sigma_beta = inv_Sigma_beta,
    S = S, n_s = n_s, p = p, d = d,
    L_f = 0L, L_s = c(0L, 0L), K = K, K_total = K_total,
    c_val = temperature, n_cpus = 1L
  )

  expected_mu <- mu_q_nu_beta
  expected_sigma <- Sigma_q_nu_beta
  prior <- matrix(0, K_total, K_total)
  prior[1:2, 1:2] <- inv_Sigma_beta
  for (j in seq_len(p)) for (r in seq_len(d)) {
    precision <- prior
    precision[3:K_total, 3:K_total] <- beta_prec[j, r] * diag(K)
    linear <- numeric(K_total)
    for (s in seq_len(S)) for (i in seq_len(n_s[[s]])) {
      gram <- list_cp_C[[s]][[i]]
      z <- Z[[s]][i, r]
      precision <- precision + eps_prec[s, j] * z^2 * gram
      residual <- list_cp_C_Y[[s]][[i]][, j] -
        as.vector(gram %*% mu_q_nu_mu[[s]][[j]])
      for (other in setdiff(seq_len(d), r)) {
        residual <- residual - Z[[s]][i, other] *
          as.vector(gram %*% expected_mu[[j]][[other]])
      }
      linear <- linear + eps_prec[s, j] * z * residual
    }
    expected_sigma[[j]][[r]] <- solve(temperature * precision)
    expected_mu[[j]][[r]] <- as.vector(
      temperature * expected_sigma[[j]][[r]] %*% linear)
  }
  for (j in seq_len(p)) for (r in seq_len(d)) {
    expect_equal(as.vector(actual$Sigma_q_nu_beta[[j]][[r]]),
                 as.vector(expected_sigma[[j]][[r]]), tolerance = 1e-10)
    expect_equal(actual$mu_q_nu_beta[[j]][[r]],
                 expected_mu[[j]][[r]], tolerance = 1e-10)
  }
})
