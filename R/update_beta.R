# =============================================================================
# update_beta.R — Regression coefficient spline coefficients nu_{beta,jr}
#
# Covers parameter block 4.2 (derivations.md §3.2). NEW module.
#
# nu_beta is SHARED across all studies. Precision and mean aggregate
# data from all S studies.
#
# (Sigma^q_beta)^{-1} = c * (sum_s sigma^{-2}_eps[s,j] *
#                          sum_i z^2_{sir} * C_i^T C_i +
#                          sigma^{-2}_beta[j,r] * I_K) + blkdiag prior
#
# mu^q_beta = c * Sigma * sum_s sigma^{-2}_eps[s,j] *
#             sum_i z_{sir} * C_i^T * residual
#
# residual r^beta: y - mu - sum_{r'≠r} z_{sir'} beta_{jr'}
#                  - sum_l a_{jl} f^{(l)}_{si} - sum_l b_{sjl} g^{(l)}_{si}
#
# Indexing: study [[s]] outermost
# =============================================================================

#' Update regression coefficient spline coefficients nu_{beta,jr}
#'
#' Cross-study aggregation: beta is shared, uses ALL studies' data.
#'
#' @export
update_nu_beta <- function(Y, C, list_cp_C, list_cp_C_Y,
                            mu_q_nu_mu,
                            mu_q_nu_beta, Sigma_q_nu_beta, Z,
                            mu_q_zeta, mu_q_nu_phi,
                            mu_q_xi, mu_q_nu_psi,
                            mu_q_a, mu_q_b_specific,
                            mu_q_recip_sigsq_eps, mu_q_recip_sigsq_beta,
                            inv_Sigma_beta,
                            S, n_s, p, d, L_f, L_s, K, K_total,
                            c_val = 1, n_cpus = 1) {

  if (d == 0 || is.null(Z) || is.null(mu_q_nu_beta)) return(NULL)

  mu_q_nu_beta_new <- vector("list", p)
  Sigma_q_nu_beta_new <- vector("list", p)

  for (j in 1:p) {

    mu_q_nu_beta_new[[j]] <- vector("list", d)
    Sigma_q_nu_beta_new[[j]] <- vector("list", d)

    # Precompute list_cp_C_nu_mu for all (s,i) — needed for mu subtraction
    list_cp_C_nu_mu <- vector("list", S)
    for (s in 1:S) {
      list_cp_C_nu_mu[[s]] <- lapply(1:n_s[s], function(i) {
        sapply(1:p, function(jj) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[jj]])
      })
    }

    for (r in 1:d) {

      # ---- Cross-study precision accumulation ----
      prec_data <- matrix(0, K_total, K_total)
      for (s in 1:S) {
        eps_prec <- mu_q_recip_sigsq_eps[s, j]
        for (i in 1:n_s[s]) {
          z_sir2 <- Z[[s]][i, r]^2
          prec_data <- prec_data + eps_prec * z_sir2 * list_cp_C[[s]][[i]]
        }
      }

      # Prior precision
      inv_prior <- blkdiag(inv_Sigma_beta,
                           mu_q_recip_sigsq_beta[j, r] * diag(K))

      # Posterior covariance
      Sigma_q_nu_beta_new[[j]][[r]] <- solve(c_val * (prec_data + inv_prior))

      # ---- Linear term: cross-study sum ----
      sum_mu <- rep(0, K_total)

      for (s in 1:S) {
        eps_prec <- mu_q_recip_sigsq_eps[s, j]
        for (i in 1:n_s[s]) {

          z_sir <- Z[[s]][i, r]
          if (abs(z_sir) < 1e-15) next
          C_si <- C[[s]][[i]]
          cp_C_si <- list_cp_C[[s]][[i]]

          # Start: C^T y - C^T C mu
          ct_res <- list_cp_C_Y[[s]][[i]][, j] -
                    list_cp_C_nu_mu[[s]][[i]][, j]

          # Subtract OTHER covariates r' != r
          for (r_other in setdiff(1:d, r)) {
            ct_beta_other <- cp_C_si %*% mu_q_nu_beta[[j]][[r_other]]
            ct_res <- ct_res - Z[[s]][i, r_other] * as.vector(ct_beta_other)
          }

          # Subtract shared factor contributions
          if (L_f > 0) {
            for (l in 1:L_f) {
              f_sil <- as.vector(
                C_si %*% mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ]
              )
              ct_res <- ct_res - mu_q_a[j, l] * crossprod(C_si, f_sil)
            }
          }

          # Subtract specific factor contributions
          if (L_s > 0 && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
            for (l in 1:L_s) {
              g_sil <- as.vector(
                C_si %*% mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ]
              )
              ct_res <- ct_res - mu_q_b_specific[[s]][j, l] *
                                crossprod(C_si, g_sil)
            }
          }

          sum_mu <- sum_mu + eps_prec * z_sir * as.vector(ct_res)
        }
      }

      # ---- Posterior mean ----
      mu_q_nu_beta_new[[j]][[r]] <- as.vector(
        c_val * Sigma_q_nu_beta_new[[j]][[r]] %*% sum_mu
      )
    }
  }

  create_named_list(mu_q_nu_beta = mu_q_nu_beta_new,
                    Sigma_q_nu_beta = Sigma_q_nu_beta_new)
}
