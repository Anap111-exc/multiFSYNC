# =============================================================================
# update_mu.R — Mean function spline coefficients nu_{mu,sj} update
#
# Covers parameter block 4.1 (derivations.md §3.1).
#
# Step 1: Function signature
# Step 2: Dimension assertions
# Step 3: Formula translation (LaTeX -> R)
# Step 4: Residual assembly (r^mu excludes mu's own contribution)
# Step 5: Return updated parameters
#
# Key formulas:
#   (Sigma^q)^{-1} = c * (sigma^{-2}_eps * sum_i C_i^T C_i + sigma^{-2}_mu * I_K)
#                   BUT in code: blkdiag(inv_Sigma_beta, sigma^{-2}_mu * diag(K))
#   mu^q = c * Sigma^q * sigma^{-2}_eps * sum_i C_i^T (y_i - r^mu_i)
#
#   r^mu_{sij} = sum_r z_{sir} beta_{jr}
#              + sum_l a_{jl} f^{(l)}_{si}
#              + sum_l b_{sjl} g^{(l)}_{si}
#
# Indexing: study [[s]] outermost, [[j]] variable, [[i]] individual
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update mean function spline coefficients nu_{mu,sj}
#'
#' @param Y Observation data: list[[s]][[i]][[j]] = vector of length n_obs
#' @param C Spline design matrices: C[[s]][[i]] = n_obs x (K+2) matrix
#' @param list_cp_C Precomputed crossprod(C): list_cp_C[[s]][[i]]
#' @param list_cp_C_Y Precomputed crossprod(C, Y): list_cp_C_Y[[s]][[i]][,j] = (K+2) x 1
#' @param sum_list_cp_C Precomputed sum_i crossprod(C): sum_list_cp_C[[s]] = (K+2)x(K+2)
#' @param mu_q_recip_sigsq_eps S x p matrix, E_q[1/sigma^2_eps]
#' @param mu_q_recip_sigsq_mu S x p matrix, E_q[1/sigma^2_mu]
#' @param inv_Sigma_beta 2 x 2 matrix, prior precision for linear part
#' @param mu_q_nu_beta Variational mean for nu_beta (or NULL if d=0)
#' @param Z Covariate matrices per study (or NULL)
#' @param mu_q_zeta mu_q_zeta[[s]][[l]] = n_s[s] x L_f matrix of score means
#' @param mu_q_nu_phi mu_q_nu_phi[[l]] = (K+2) x L_f matrix of eigenfunction coeffs
#' @param mu_q_xi mu_q_xi[[s]][[l]] = n_s[s] x L_s matrix of specific score means
#' @param mu_q_nu_psi mu_q_nu_psi[[s]][[l]] = (K+2) x L_s matrix of specific eigenfunctions
#' @param mu_q_a p x L_f matrix, spike-and-slab mean of shared loadings
#' @param mu_q_b_specific List of S matrices, p x L_s, spike-and-slab mean of specific loadings
#' @param S Number of studies
#' @param n_s Vector of individuals per study
#' @param p Number of variables
#' @param d Number of covariates
#' @param L_f Number of shared factors
#' @param L_s Number of specific factors
#' @param K Spline basis dimension (penalized part)
#' @param K_total Total spline dimension (= K + 2)
#' @param c_val Temperature constant c = 1/T
#' @param n_cpus Number of CPU cores
#'
#' @return List with updated mu_q_nu_mu, Sigma_q_nu_mu, inv_Sigma_q_nu_mu
#'
#' @noRd
update_nu_mu <- function(Y, C, list_cp_C, list_cp_C_Y,
                          mu_q_nu_mu, Sigma_q_nu_mu,
                          sum_list_cp_C,
                          mu_q_recip_sigsq_eps, mu_q_recip_sigsq_mu,
                          inv_Sigma_beta,
                          mu_q_nu_beta, Z,
                          mu_q_zeta, mu_q_nu_phi,
                          mu_q_xi, mu_q_nu_psi,
                          mu_q_a, mu_q_b_specific,
                          S, n_s, p, d, L_f, L_s, K, K_total,
                          c_val = 1, n_cpus = 1) {

  # ---- Step 2: Dimension assertions ----
  stopifnot(length(Y) == S)
  stopifnot(dim(mu_q_recip_sigsq_eps) == c(S, p))
  stopifnot(dim(mu_q_recip_sigsq_mu) == c(S, p))
  stopifnot(dim(inv_Sigma_beta) == c(2, 2))

  # ---- Step 3 & 4: Per-study, per-variable update ----

  mu_q_nu_mu <- vector("list", S)
  Sigma_q_nu_mu <- vector("list", S)
  inv_Sigma_q_nu_mu <- vector("list", S)

  for (s in 1:S) {
    L_ss <- .L_s_at(L_s, s)

    mu_q_nu_mu[[s]] <- vector("list", p)
    Sigma_q_nu_mu[[s]] <- vector("list", p)
    inv_Sigma_q_nu_mu[[s]] <- vector("list", p)

    for (j in 1:p) {

      # ---- Prior precision: blkdiag(inv_Sigma_beta, sigma_mu * diag(K)) ----
      inv_prior <- blkdiag(inv_Sigma_beta, mu_q_recip_sigsq_mu[s, j] * diag(K))

      # ---- Posterior precision: c * (prior + sigma_eps * sum C^T C) ----
      prec <- c_val * (inv_prior + mu_q_recip_sigsq_eps[s, j] * sum_list_cp_C[[s]])

      # ---- Posterior covariance ----
      Sigma_q_nu_mu[[s]][[j]] <- .inverse_spd(
        prec, context = sprintf("nu_mu[s=%d,j=%d]", s, j))
      inv_Sigma_q_nu_mu[[s]][[j]] <- inv_prior

      # ---- Residual assembly: r^mu_{sij} (excludes mu's own contribution) ----
      sum_term_mu_sj <- rep(0, K_total)

      for (i in 1:n_s[s]) {

        # Start from C_i^T y_{sij}
        residual_contrib <- rep(0, K_total)

        # --- Beta contribution: sum_r z_{sir} * C_i nu_beta_{jr} ---
        if (!is.null(mu_q_nu_beta) && !is.null(Z) && d > 0) {
          for (r in 1:d) {
            beta_val <- as.vector(C[[s]][[i]] %*% mu_q_nu_beta[[j]][[r]])
            residual_contrib <- residual_contrib +
              Z[[s]][i, r] * as.vector(crossprod(C[[s]][[i]], beta_val))
          }
        }

        # --- Shared factor contribution: sum_l a_{jl} * f^{(l)}_{si} ---
        if (L_f > 0) {
          for (l in 1:L_f) {
            f_sil <- as.vector(
              C[[s]][[i]] %*% mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ]
            )
            residual_contrib <- residual_contrib +
              mu_q_a[j, l] * as.vector(crossprod(C[[s]][[i]], f_sil))
          }
        }

        # --- Specific factor contribution: sum_l b_{sjl} * g^{(l)}_{si} ---
        if (L_ss > 0L && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
          for (l in seq_len(L_ss)) {
            g_sil <- as.vector(
              C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ]
            )
            residual_contrib <- residual_contrib +
              mu_q_b_specific[[s]][j, l] * as.vector(crossprod(C[[s]][[i]], g_sil))
          }
        }

        # sum_i C_i^T y_{sij} - sum_i C_i^T (residual)
        sum_term_mu_sj <- sum_term_mu_sj +
          list_cp_C_Y[[s]][[i]][, j] - residual_contrib
      }

      # ---- Step 5: Posterior mean ----
      mu_q_nu_mu[[s]][[j]] <- as.vector(
        c_val * Sigma_q_nu_mu[[s]][[j]] %*%
          (mu_q_recip_sigsq_eps[s, j] * sum_term_mu_sj)
      )
    }
  }

  create_named_list(mu_q_nu_mu, Sigma_q_nu_mu, inv_Sigma_q_nu_mu)
}
