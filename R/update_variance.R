# =============================================================================
# update_variance.R — Variance parameter update functions
#
# Covers parameter blocks:
#   4.9  sigma^2_eps     (measurement error variance)
#   4.10 sigma^2_mu, sigma^2_beta  (mean & regression coefficient variances)
#   4.11 sigma^2_phi, sigma^2_psi  (eigenfunction variances)
#   4.12 a_*             (Half-Cauchy auxiliary variables)
#
# Each function follows the 5-step template:
#   1. Function signature with inputs/outputs
#   2. Dimension assertions
#   3. Formula translation (LaTeX -> R)
#   4. Residual assembly (where applicable)
#   5. Return updated parameters
#
# Key rules:
#   - kappa_q formulas include c as per derivations §4.9-4.12
#   - Penalized part only: exclude linear coefficients [1:2] for sigma^2_mu etc.
#   - kappa_q_a = 2c - 1 > 0 => c > 1/2 => T < 2 (critical constraint)
#   - Indexing: study [[s]] outermost
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

# ---- 4.12: Half-Cauchy auxiliary variable a_{eps,sj} ----

#' Update auxiliary variable for measurement error variance
#'
#' @param mu_q_recip_sigsq_eps S x p matrix, E_q[1/sigma^2_eps]
#' @param A Half-Cauchy scale hyperparameter (default 1e5)
#' @param c_val Temperature constant c = 1/T
#' @return List with updated mu_q_recip_a_eps, kappa_q_a, lambda_q_a_eps
update_a_eps <- function(mu_q_recip_sigsq_eps, A = 1e5, c_val = 1) {
  kappa_q_a <- 2 * c_val - 1
  lambda_q_a_eps <- c_val * (mu_q_recip_sigsq_eps + 1 / A^2)
  mu_q_recip_a_eps <- kappa_q_a / lambda_q_a_eps
  create_named_list(mu_q_recip_a_eps, kappa_q_a, lambda_q_a_eps)
}

# ---- 4.12: Half-Cauchy auxiliary variable a_{mu,sj} ----

#' Update auxiliary variable for mean function variance
#'
#' @param mu_q_recip_sigsq_mu S x p matrix, E_q[1/sigma^2_mu]
#' @param A Half-Cauchy scale
#' @param c_val Temperature constant
#' @return List with updated parameters
update_a_mu <- function(mu_q_recip_sigsq_mu, A = 1e5, c_val = 1) {
  kappa_q_a <- 2 * c_val - 1
  lambda_q_a_mu <- c_val * (mu_q_recip_sigsq_mu + 1 / A^2)
  mu_q_recip_a_mu <- kappa_q_a / lambda_q_a_mu
  create_named_list(mu_q_recip_a_mu, kappa_q_a, lambda_q_a_mu)
}

# ---- 4.12: Half-Cauchy auxiliary variable a_{beta,jr} (NEW) ----

#' Update auxiliary variable for regression coefficient variance
#'
#' @param mu_q_recip_sigsq_beta p x d matrix, E_q[1/sigma^2_beta]
#' @param A Half-Cauchy scale
#' @param c_val Temperature constant
#' @return List with updated parameters
update_a_beta <- function(mu_q_recip_sigsq_beta, A = 1e5, c_val = 1) {
  kappa_q_a <- 2 * c_val - 1
  lambda_q_a_beta <- c_val * (mu_q_recip_sigsq_beta + 1 / A^2)
  mu_q_recip_a_beta <- kappa_q_a / lambda_q_a_beta
  create_named_list(mu_q_recip_a_beta, kappa_q_a, lambda_q_a_beta)
}

# ---- 4.12: Half-Cauchy auxiliary variable a_{phi,ml} ----

#' Update auxiliary variable for shared eigenfunction variance
#'
#' @param mu_q_recip_sigsq_phi List of vectors, E_q[1/sigma^2_phi], indexed [[l]][m]
#' @param A Half-Cauchy scale
#' @param c_val Temperature constant
#' @return List with updated parameters (mu_q_recip_a_phi returned as list of vectors)
update_a_phi <- function(mu_q_recip_sigsq_phi, A = 1e5, c_val = 1) {
  kappa_q_a <- 2 * c_val - 1
  lambda_q_a_phi <- lapply(mu_q_recip_sigsq_phi, function(v) c_val * (v + 1 / A^2))
  mu_q_recip_a_phi <- lapply(lambda_q_a_phi, function(x) kappa_q_a / x)
  create_named_list(mu_q_recip_a_phi, kappa_q_a, lambda_q_a_phi)
}

# ---- 4.12: Half-Cauchy auxiliary variable a_{psi,slm} (NEW) ----

#' Update auxiliary variable for specific eigenfunction variance
#'
#' @param mu_q_recip_sigsq_psi List of list of vectors, E_q[1/sigma^2_psi], indexed [[s]][[l]][m]
#' @param A Half-Cauchy scale
#' @param c_val Temperature constant
#' @return List with updated parameters (mu_q_recip_a_psi returned as list of list of vectors)
update_a_psi <- function(mu_q_recip_sigsq_psi, A = 1e5, c_val = 1) {
  kappa_q_a <- 2 * c_val - 1
  lambda_q_a_psi <- lapply(mu_q_recip_sigsq_psi, function(s_list) lapply(s_list, function(v) c_val * (v + 1 / A^2)))
  mu_q_recip_a_psi <- lapply(lambda_q_a_psi, function(s_list) lapply(s_list, function(x) kappa_q_a / x))
  create_named_list(mu_q_recip_a_psi, kappa_q_a, lambda_q_a_psi)
}

# ---- 4.9: Measurement error variance sigma^2_{eps,sj} ----

#' Compute expected residual sum of squares for variable j, study s, individual i
#'
#' @param s Study index
#' @param i Individual index within study s
#' @param j Variable index
#' @keywords internal
compute_rss_single <- function(s, i, j, Y, C, list_cp_C,
                                mu_q_nu_mu, Sigma_q_nu_mu,
                                mu_q_nu_beta, Sigma_q_nu_beta, Z,
                                mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi,
                                mu_q_xi, Sigma_q_xi, mu_q_nu_psi,
                                mu_q_a, term_a,
                                mu_q_b_specific, term_b_specific,
                                L_f, L_s) {

  C_si <- C[[s]][[i]]
  y_sij <- Y[[s]][[i]][[j]]

  # Expected contributions using current variational means
  y_hat <- as.vector(C_si %*% mu_q_nu_mu[[s]][[j]])

  # Beta contribution
  if (!is.null(mu_q_nu_beta) && !is.null(Z)) {
    d <- ncol(Z[[s]])
    for (r in 1:d) {
      y_hat <- y_hat + Z[[s]][i, r] * as.vector(C_si %*% mu_q_nu_beta[[j]][[r]])
    }
  }

  # Shared factor contribution
  for (l in seq_len(L_f)) {
    phi_zeta <- as.vector(mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ])
    y_hat <- y_hat + mu_q_a[j, l] * as.vector(C_si %*% phi_zeta)
  }

  # Specific factor contribution
  if (L_s > 0) {
    for (l in seq_len(L_s)) {
      psi_xi <- as.vector(mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ])
      y_hat <- y_hat + mu_q_b_specific[[s]][j, l] * as.vector(C_si %*% psi_xi)
    }
  }

  residual <- y_sij - y_hat
  rss <- as.numeric(crossprod(residual))

  # Variance contribution from mean function posterior
  rss <- rss + tr(list_cp_C[[s]][[i]] %*% Sigma_q_nu_mu[[s]][[j]])

  # Variance contributions from shared factors
  for (l in seq_len(L_f)) {
    phi_zeta <- as.vector(mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ])
    mean_norm2 <- as.numeric(crossprod(C_si %*% phi_zeta))
    f_norm2_expect <- mean_norm2 +
      tr(list_cp_C[[s]][[i]] %*%
         (mu_q_nu_phi[[l]] %*% Sigma_q_zeta[[s]][[l]][[i]] %*% t(mu_q_nu_phi[[l]])))
    rss <- rss + term_a[j, l] * f_norm2_expect - mu_q_a[j, l]^2 * mean_norm2
  }

  # Variance contributions from specific factors
  if (L_s > 0) {
    for (l in seq_len(L_s)) {
      psi_xi <- as.vector(mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ])
      mean_norm2_spec <- as.numeric(crossprod(C_si %*% psi_xi))
      g_norm2_expect <- mean_norm2_spec +
        tr(list_cp_C[[s]][[i]] %*%
           (mu_q_nu_psi[[s]][[l]] %*% Sigma_q_xi[[s]][[l]][[i]] %*% t(mu_q_nu_psi[[s]][[l]])))
      rss <- rss + term_b_specific[[s]][j, l] * g_norm2_expect -
             mu_q_b_specific[[s]][j, l]^2 * mean_norm2_spec
    }
  }

  # Variance contributions from beta coefficients
  if (!is.null(mu_q_nu_beta) && !is.null(Sigma_q_nu_beta) && !is.null(Z)) {
    for (r in 1:ncol(Z[[s]])) {
      rss <- rss + Z[[s]][i, r]^2 * tr(list_cp_C[[s]][[i]] %*% Sigma_q_nu_beta[[j]][[r]])
    }
  }

  rss
}

#' Update measurement error variance sigma^2_{eps,sj}
#'
#' @param Y Observation data list
#' @param C Spline design matrices
#' @param list_cp_C Precomputed crossprod(C) per study-individual
#' @param mu_q_nu_mu Variational mean of nu_mu
#' @param Sigma_q_nu_mu Variational covariance of nu_mu
#' @param mu_q_nu_beta Variational mean of nu_beta (or NULL)
#' @param Z Covariate matrices per study (or NULL)
#' @param mu_q_zeta Variational mean of zeta
#' @param Sigma_q_zeta Variational covariance of zeta
#' @param mu_q_nu_phi Variational mean of nu_phi
#' @param mu_q_xi Variational mean of xi
#' @param Sigma_q_xi Variational covariance of xi
#' @param mu_q_nu_psi Variational mean of nu_psi
#' @param mu_q_a Spike-and-slab mean of shared loadings a
#' @param term_a E[a^2] for shared loadings
#' @param mu_q_b_specific Spike-and-slab mean of specific loadings b
#' @param term_b_specific E[b^2] for specific loadings
#' @param mu_q_recip_a_eps Current E[1/a_eps]
#' @param S Number of studies
#' @param n_s Vector of individuals per study
#' @param p Number of variables
#' @param L_f Number of shared factors
#' @param L_s Number of specific factors
#' @param c_val Temperature constant c = 1/T
#' @param n_cpus Number of CPU cores
#' @return List with updated kappa, lambda, mu_q_recip_sigsq_eps
#'
#' @export
update_sigsq_eps <- function(Y, C, list_cp_C,
                              mu_q_nu_mu, Sigma_q_nu_mu,
                              mu_q_nu_beta, Z,
                              mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi,
                              mu_q_xi, Sigma_q_xi, mu_q_nu_psi,
                              mu_q_a, term_a,
                              mu_q_b_specific, term_b_specific,
                              mu_q_recip_a_eps,
                               S, n_s, p, L_f, L_s,
                               total_obs_sj = NULL,
                               c_val = 1, n_cpus = 1) {

  # Dimensions
  stopifnot(nrow(mu_q_recip_a_eps) == S, ncol(mu_q_recip_a_eps) == p)

  # kappa: c * (total_obs_sj[s,j] + 1)/2 + c - 1
  kappa_q_sigsq_eps <- matrix(NA, nrow = S, ncol = p)
  for (s in 1:S) {
    for (j in 1:p) {
      nobs_sj <- if (!is.null(total_obs_sj)) total_obs_sj[s, j] else length(Y[[1]][[1]][[1]]) * n_s[s]
      kappa_q_sigsq_eps[s, j] <- c_val * (nobs_sj + 1) / 2 + c_val - 1
    }
  }

  # lambda: c * (E[a^{-1}_eps] + 0.5 * Σ_i RSS_{sij})
  lambda_q_sigsq_eps <- matrix(NA, nrow = S, ncol = p)

  for (s in 1:S) {
    for (j in 1:p) {
      rss_sum <- sum(sapply(1:n_s[s], function(i) {
        compute_rss_single(s, i, j, Y, C, list_cp_C,
                           mu_q_nu_mu, Sigma_q_nu_mu,
                           mu_q_nu_beta, Sigma_q_nu_beta, Z,
                           mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi,
                           mu_q_xi, Sigma_q_xi, mu_q_nu_psi,
                           mu_q_a, term_a,
                           mu_q_b_specific, term_b_specific,
                           L_f, L_s)
      }))
      lambda_q_sigsq_eps[s, j] <- c_val * (mu_q_recip_a_eps[s, j] + 0.5 * rss_sum)
    }
  }

  mu_q_recip_sigsq_eps <- kappa_q_sigsq_eps / lambda_q_sigsq_eps

  # Expected log sigma^2 for ELBO: E[log sigma^2] = log(lambda) - digamma(kappa)
  mu_q_log_sigsq_eps <- log(lambda_q_sigsq_eps) - digamma(kappa_q_sigsq_eps)

  create_named_list(kappa_q_sigsq_eps, lambda_q_sigsq_eps,
                    mu_q_recip_sigsq_eps, mu_q_log_sigsq_eps)
}

# ---- 4.10: Mean function variance sigma^2_{mu,sj} ----

#' Update mean function variance sigma^2_{mu,sj}
#'
#' @param mu_q_nu_mu Variational mean, list[[s]][[j]] (K+2)x1
#' @param Sigma_q_nu_mu Variational covariance, list[[s]][[j]] (K+2)x(K+2)
#' @param mu_q_recip_a_mu Current E[1/a_mu], S x p matrix
#' @param S Number of studies
#' @param p Number of variables
#' @param K Spline basis dimension (excl. linear part)
#' @param c_val Temperature constant
#' @param n_cpus Number of CPU cores
#' @return List with updated parameters
#'
#' @export
update_sigsq_mu <- function(mu_q_nu_mu, Sigma_q_nu_mu,
                             mu_q_recip_a_mu,
                             S, p, K, c_val = 1, n_cpus = 1) {

  kappa_q_sigsq_mu <- c_val * (K + 1) / 2 + c_val - 1

  lambda_q_sigsq_mu <- matrix(NA, nrow = S, ncol = p)

  for (s in 1:S) {
    for (j in 1:p) {
      nu_pen <- mu_q_nu_mu[[s]][[j]][-c(1:2)]
      Sigma_pen <- Sigma_q_nu_mu[[s]][[j]][-c(1:2), -c(1:2)]
      lambda_q_sigsq_mu[s, j] <- c_val * (
        mu_q_recip_a_mu[s, j] + 0.5 * (as.numeric(crossprod(nu_pen)) + tr(Sigma_pen))
      )
    }
  }

  mu_q_recip_sigsq_mu <- kappa_q_sigsq_mu / lambda_q_sigsq_mu
  mu_q_log_sigsq_mu <- log(lambda_q_sigsq_mu) - digamma(kappa_q_sigsq_mu)

  create_named_list(kappa_q_sigsq_mu, lambda_q_sigsq_mu,
                    mu_q_recip_sigsq_mu, mu_q_log_sigsq_mu)
}

# ---- 4.10: Regression coefficient variance sigma^2_{beta,jr} (NEW) ----

#' Update regression coefficient variance sigma^2_{beta,jr}
#'
#' @param mu_q_nu_beta Variational mean, list[[j]][[r]] (K+2)x1
#' @param Sigma_q_nu_beta Variational covariance, list[[j]][[r]] (K+2)x(K+2)
#' @param mu_q_recip_a_beta Current E[1/a_beta], p x d matrix
#' @param p Number of variables
#' @param d Number of covariates
#' @param K Spline basis dimension
#' @param c_val Temperature constant
#' @param n_cpus Number of CPU cores
#' @return List with updated parameters
#'
#' @export
update_sigsq_beta <- function(mu_q_nu_beta, Sigma_q_nu_beta,
                               mu_q_recip_a_beta,
                               p, d, K, c_val = 1, n_cpus = 1) {

  if (is.null(mu_q_nu_beta) || d == 0) return(NULL)

  kappa_q_sigsq_beta <- c_val * (K + 1) / 2 + c_val - 1

  lambda_q_sigsq_beta <- matrix(NA, nrow = p, ncol = d)

  for (j in 1:p) {
    for (r in 1:d) {
      nu_pen <- mu_q_nu_beta[[j]][[r]][-c(1:2)]
      Sigma_pen <- Sigma_q_nu_beta[[j]][[r]][-c(1:2), -c(1:2)]
      lambda_q_sigsq_beta[j, r] <- c_val * (
        mu_q_recip_a_beta[j, r] + 0.5 * (as.numeric(crossprod(nu_pen)) + tr(Sigma_pen))
      )
    }
  }

  mu_q_recip_sigsq_beta <- kappa_q_sigsq_beta / lambda_q_sigsq_beta
  mu_q_log_sigsq_beta <- log(lambda_q_sigsq_beta) - digamma(kappa_q_sigsq_beta)

  create_named_list(kappa_q_sigsq_beta, lambda_q_sigsq_beta,
                    mu_q_recip_sigsq_beta, mu_q_log_sigsq_beta)
}

# ---- 4.11: Shared eigenfunction variance sigma^2_{phi,ml} ----

#' Update shared eigenfunction variance sigma^2_{phi,ml}
#'
#' @param mu_q_nu_phi Variational mean, list[[l]][,m] (K+2)x1
#' @param Sigma_q_nu_phi Variational covariance, list[[l]][[m]] (K+2)x(K+2)
#' @param mu_q_recip_a_phi Current E[1/a_phi], list of vectors indexed [[l]][m]
#' @param L_f Number of shared factors
#' @param M_f Vector of FPCA component counts per shared factor
#' @param K Spline basis dimension
#' @param c_val Temperature constant
#' @return List with updated parameters (mu_q_recip_sigsq_phi, mu_q_log_sigsq_phi as lists)
#'
#' @export
update_sigsq_phi <- function(mu_q_nu_phi, Sigma_q_nu_phi,
                              mu_q_recip_a_phi,
                              L_f, M_f, K, c_val = 1) {

  kappa_q_sigsq_phi <- c_val * (K + 1) / 2 + c_val - 1

  lambda_q_sigsq_phi <- lapply(seq_len(L_f), function(l) rep(NA, M_f[l]))

  for (l in seq_len(L_f)) {
    for (m in 1:M_f[l]) {
      nu_pen <- mu_q_nu_phi[[l]][-c(1:2), m]
      Sigma_pen <- Sigma_q_nu_phi[[l]][[m]][-c(1:2), -c(1:2)]
      lambda_q_sigsq_phi[[l]][m] <- c_val * (
        mu_q_recip_a_phi[[l]][m] + 0.5 * (as.numeric(crossprod(nu_pen)) + tr(Sigma_pen))
      )
    }
  }

  mu_q_recip_sigsq_phi <- lapply(lambda_q_sigsq_phi, function(x) kappa_q_sigsq_phi / x)
  mu_q_log_sigsq_phi <- lapply(lambda_q_sigsq_phi, function(x) log(x) - digamma(kappa_q_sigsq_phi))

  create_named_list(kappa_q_sigsq_phi, lambda_q_sigsq_phi,
                    mu_q_recip_sigsq_phi, mu_q_log_sigsq_phi)
}

# ---- 4.11: Specific eigenfunction variance sigma^2_{psi,slm} (NEW) ----

#' Update specific eigenfunction variance sigma^2_{psi,slm}
#'
#' @param mu_q_nu_psi Variational mean, list[[s]][[l]][,m] (K+2)x1
#' @param Sigma_q_nu_psi Variational covariance, list[[s]][[l]][[m]] (K+2)x(K+2)
#' @param mu_q_recip_a_psi Current E[1/a_psi], list of list of vectors indexed [[s]][[l]][m]
#' @param S Number of studies
#' @param L_s Number of specific factors
#' @param M_s List of FPCA component counts per study per specific factor
#' @param K Spline basis dimension
#' @param c_val Temperature constant
#' @return List with updated parameters (mu_q_recip_sigsq_psi, mu_q_log_sigsq_psi as list of lists)
#'
#' @export
update_sigsq_psi <- function(mu_q_nu_psi, Sigma_q_nu_psi,
                              mu_q_recip_a_psi,
                              S, L_s, M_s, K, c_val = 1) {

  if (is.null(mu_q_nu_psi) || L_s == 0) return(NULL)

  kappa_q_sigsq_psi <- c_val * (K + 1) / 2 + c_val - 1

  lambda_q_sigsq_psi <- lapply(1:S, function(s) lapply(1:L_s, function(l) rep(NA, M_s[[s]][l])))

  for (s in 1:S) {
    for (l in 1:L_s) {
      for (m in 1:M_s[[s]][l]) {
        nu_pen <- mu_q_nu_psi[[s]][[l]][-c(1:2), m]
        Sigma_pen <- Sigma_q_nu_psi[[s]][[l]][[m]][-c(1:2), -c(1:2)]
        lambda_q_sigsq_psi[[s]][[l]][m] <- c_val * (
          mu_q_recip_a_psi[[s]][[l]][m] + 0.5 * (as.numeric(crossprod(nu_pen)) + tr(Sigma_pen))
        )
      }
    }
  }

  mu_q_recip_sigsq_psi <- lapply(lambda_q_sigsq_psi, function(s_list) lapply(s_list, function(x) kappa_q_sigsq_psi / x))
  mu_q_log_sigsq_psi <- lapply(lambda_q_sigsq_psi, function(s_list) lapply(s_list, function(x) log(x) - digamma(kappa_q_sigsq_psi)))

  create_named_list(kappa_q_sigsq_psi, lambda_q_sigsq_psi,
                    mu_q_recip_sigsq_psi, mu_q_log_sigsq_psi)
}

# ---- Main variance update wrapper ----

#' Update all variance and auxiliary variable parameters (calls all sub-functions)
#'
#' @param ... All current variational parameters
#' @param A Half-Cauchy scale
#' @param c_val Temperature constant
#' @param n_cpus Number of CPU cores
#' @return List with all updated variance-related parameters
#'
#' @export
update_all_variances <- function(Y, C, list_cp_C,
                                  mu_q_nu_mu, Sigma_q_nu_mu,
                                  mu_q_nu_beta, Sigma_q_nu_beta, Z,
                                  mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi, Sigma_q_nu_phi,
                                  mu_q_xi, Sigma_q_xi, mu_q_nu_psi, Sigma_q_nu_psi,
                                  mu_q_a, term_a,
                                  mu_q_b_specific, term_b_specific,
                                  mu_q_recip_sigsq_eps, mu_q_recip_a_eps,
                                  mu_q_recip_sigsq_mu, mu_q_recip_a_mu,
                                  mu_q_recip_sigsq_beta, mu_q_recip_a_beta,
                                  mu_q_recip_sigsq_phi, mu_q_recip_a_phi,
                                  mu_q_recip_sigsq_psi, mu_q_recip_a_psi,
                                   S, n_s, p, d, L_f, L_s, M_f, M_s, K,
                                   total_obs_sj = NULL,
                                   A = 1e5, c_val = 1, n_cpus = 1) {

  # 1. Update sigma^2_eps (measurement error)
  res_eps <- update_sigsq_eps(Y, C, list_cp_C,
                               mu_q_nu_mu, Sigma_q_nu_mu,
                               mu_q_nu_beta, Z,
                               mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi,
                               mu_q_xi, Sigma_q_xi, mu_q_nu_psi,
                               mu_q_a, term_a,
                               mu_q_b_specific, term_b_specific,
                               mu_q_recip_a_eps,
                               S, n_s, p, L_f, L_s,
                               total_obs_sj, c_val, n_cpus)

  # 2. Update sigma^2_mu (mean function variance)
  res_mu <- update_sigsq_mu(mu_q_nu_mu, Sigma_q_nu_mu,
                             mu_q_recip_a_mu,
                             S, p, K, c_val, n_cpus)

  # 3. Update sigma^2_beta (regression coefficient variance)
  res_beta <- update_sigsq_beta(mu_q_nu_beta, Sigma_q_nu_beta,
                                 mu_q_recip_a_beta,
                                 p, d, K, c_val, n_cpus)

  # 4. Update sigma^2_phi (shared eigenfunction variance)
  res_phi <- update_sigsq_phi(mu_q_nu_phi, Sigma_q_nu_phi,
                               mu_q_recip_a_phi,
                               L_f, M_f, K, c_val)

  # 5. Update sigma^2_psi (specific eigenfunction variance)
  res_psi <- update_sigsq_psi(mu_q_nu_psi, Sigma_q_nu_psi,
                               mu_q_recip_a_psi,
                               S, L_s, M_s, K, c_val)

  # 6-10. Update auxiliary variables a_*
  res_a_eps  <- update_a_eps(res_eps$mu_q_recip_sigsq_eps, A, c_val)
  res_a_mu   <- update_a_mu(res_mu$mu_q_recip_sigsq_mu, A, c_val)

  res_a_beta <- NULL
  if (!is.null(res_beta) && d > 0) {
    res_a_beta <- update_a_beta(res_beta$mu_q_recip_sigsq_beta, A, c_val)
  }

  res_a_phi  <- update_a_phi(res_phi$mu_q_recip_sigsq_phi, A, c_val)

  res_a_psi <- NULL
  if (!is.null(res_psi) && L_s > 0) {
    res_a_psi <- update_a_psi(res_psi$mu_q_recip_sigsq_psi, A, c_val)
  }

  create_named_list(
    # sigma^2_eps
    kappa_q_sigsq_eps = res_eps$kappa_q_sigsq_eps,
    lambda_q_sigsq_eps = res_eps$lambda_q_sigsq_eps,
    mu_q_recip_sigsq_eps = res_eps$mu_q_recip_sigsq_eps,
    mu_q_log_sigsq_eps = res_eps$mu_q_log_sigsq_eps,
    mu_q_recip_a_eps = res_a_eps$mu_q_recip_a_eps,

    # sigma^2_mu
    kappa_q_sigsq_mu = res_mu$kappa_q_sigsq_mu,
    lambda_q_sigsq_mu = res_mu$lambda_q_sigsq_mu,
    mu_q_recip_sigsq_mu = res_mu$mu_q_recip_sigsq_mu,
    mu_q_log_sigsq_mu = res_mu$mu_q_log_sigsq_mu,
    mu_q_recip_a_mu = res_a_mu$mu_q_recip_a_mu,

    # sigma^2_beta
    kappa_q_sigsq_beta = if (!is.null(res_beta)) res_beta$kappa_q_sigsq_beta else NULL,
    lambda_q_sigsq_beta = if (!is.null(res_beta)) res_beta$lambda_q_sigsq_beta else NULL,
    mu_q_log_sigsq_beta = if (!is.null(res_beta)) res_beta$mu_q_log_sigsq_beta else NULL,
    mu_q_recip_sigsq_beta = if (!is.null(res_beta)) res_beta$mu_q_recip_sigsq_beta else NULL,
    mu_q_recip_a_beta = if (!is.null(res_a_beta)) res_a_beta$mu_q_recip_a_beta else NULL,

    # sigma^2_phi
    kappa_q_sigsq_phi = res_phi$kappa_q_sigsq_phi,
    lambda_q_sigsq_phi = res_phi$lambda_q_sigsq_phi,
    mu_q_recip_sigsq_phi = res_phi$mu_q_recip_sigsq_phi,
    mu_q_log_sigsq_phi = res_phi$mu_q_log_sigsq_phi,
    mu_q_recip_a_phi = res_a_phi$mu_q_recip_a_phi,

    # sigma^2_psi
    kappa_q_sigsq_psi = if (!is.null(res_psi)) res_psi$kappa_q_sigsq_psi else NULL,
    lambda_q_sigsq_psi = if (!is.null(res_psi)) res_psi$lambda_q_sigsq_psi else NULL,
    mu_q_log_sigsq_psi = if (!is.null(res_psi)) res_psi$mu_q_log_sigsq_psi else NULL,
    mu_q_recip_sigsq_psi = if (!is.null(res_psi)) res_psi$mu_q_recip_sigsq_psi else NULL,
    mu_q_recip_a_psi = if (!is.null(res_a_psi)) res_a_psi$mu_q_recip_a_psi else NULL,

    # shared kappa_q_a
    kappa_q_a = res_a_eps$kappa_q_a
  )
}
