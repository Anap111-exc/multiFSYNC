# =============================================================================
# update_xi.R — Study-specific factor FPCA scores xi^{(sl)}_{si}
#
# Covers parameter block 4.4 (derivations.md §3.4).
#
# Structure parallel to update_zeta.R, but:
#   - Within study s only (no cross-study aggregation)
#   - Uses b_{sjl} loadings instead of a_{jl}
#   - Uses nu_psi / xi instead of nu_phi / zeta
#   - H_psi instead of H_phi
#   - Residual excludes current specific factor l,
#     includes ALL shared factors
#
# Indexing: study [[s]] outermost
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update study-specific factor FPCA scores xi
#'
#' @export
update_xi <- function(Y, C, list_cp_C, list_cp_C_Y,
                       mu_q_nu_mu,
                       mu_q_nu_beta, Z,
                       mu_q_zeta, mu_q_nu_phi,
                       mu_q_xi, Sigma_q_xi,
                       mu_q_nu_psi, Sigma_q_nu_psi,
                       mu_q_a,
                       mu_q_b_specific, term_b_specific,
                       mu_q_recip_sigsq_eps,
                       S, n_s, p, d, L_f, L_s, M_s,
                       c_val = 1, n_cpus = 1) {

  if (L_s == 0 || is.null(mu_q_nu_psi)) return(NULL)

  mu_q_xi_new <- vector("list", S)
  Sigma_q_xi_new <- vector("list", S)
  tr_xi_specific <- array(NA, dim = c(S, L_s, max(n_s)))

  for (s in 1:S) {

    mu_q_xi_new[[s]] <- vector("list", L_s)
    Sigma_q_xi_new[[s]] <- vector("list", L_s)

    # Precompute C_i^T C_i * mu_q_nu_mu_j
    list_cp_C_nu_mu <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(j) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[j]])
    })

    for (l in 1:L_s) {

      M_sl <- M_s[[s]][l]

      # Prior precision: I_{M_sl}
      inv_Sigma_xi <- diag(M_sl)

      mu_q_xi_new[[s]][[l]] <- matrix(0, nrow = n_s[s], ncol = M_sl)
      Sigma_q_xi_new[[s]][[l]] <- vector("list", n_s[s])

      # Scalar precision: sum_j sigma^{-2}_eps[s,j] * E[b^2_{sjl}]
      sum_sigma_b_sl <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l])

      for (i in 1:n_s[s]) {

        C_si <- C[[s]][[i]]
        cp_C_si <- list_cp_C[[s]][[i]]

        # ---- H_psi matrix (M_sl x M_sl) ----
        H_mean <- crossprod(mu_q_nu_psi[[s]][[l]], cp_C_si) %*%
                  mu_q_nu_psi[[s]][[l]]
        H_var <- diag(sapply(1:M_sl, function(m) {
          tr(cp_C_si %*% Sigma_q_nu_psi[[s]][[l]][[m]])
        }), nrow = M_sl, ncol = M_sl)
        H_psi <- H_mean + H_var

        # ---- Posterior covariance ----
        prec_xi <- c_val * (sum_sigma_b_sl * H_psi + inv_Sigma_xi)
        Sigma_q_xi_new[[s]][[l]][[i]] <- solve(prec_xi + 1e-8 * diag(M_sl))

        # ---- Linear term ----
        sum_mu <- rep(0, ncol(C_si))

        # Positive: sum_j b_{sjl} * sigma^{-2}_eps * (C^T y - C^T C nu_mu)
        for (jj in 1:p) {
          coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_b_specific[[s]][jj, l]
          if (abs(coef) < 1e-15) next
          ct_res <- list_cp_C_Y[[s]][[i]][, jj] - list_cp_C_nu_mu[[i]][, jj]
          sum_mu <- sum_mu + coef * as.vector(ct_res)
        }

        # Subtract beta
        if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
          for (jj in 1:p) {
            coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_b_specific[[s]][jj, l]
            if (abs(coef) < 1e-15) next
            for (r in 1:d) {
              ct_beta <- cp_C_si %*% mu_q_nu_beta[[jj]][[r]]
              sum_mu <- sum_mu - coef * Z[[s]][i, r] * as.vector(ct_beta)
            }
          }
        }

        # Subtract ALL shared factor contributions
        if (L_f > 0) {
          for (l_shared in 1:L_f) {
            f_other <- as.vector(
              C_si %*% mu_q_nu_phi[[l_shared]] %*% mu_q_zeta[[s]][[l_shared]][i, ]
            )
            ct_f <- crossprod(C_si, f_other)
            for (jj in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, jj] *
                      mu_q_b_specific[[s]][jj, l] * mu_q_a[jj, l_shared]
              if (abs(coef) < 1e-15) next
              sum_mu <- sum_mu - coef * as.vector(ct_f)
            }
          }
        }

        # Subtract OTHER specific factors (l_other != l)
        if (L_s > 1) {
          for (l_other in setdiff(1:L_s, l)) {
            g_other <- as.vector(
              C_si %*% mu_q_nu_psi[[s]][[l_other]] %*% mu_q_xi[[s]][[l_other]][i, ]
            )
            ct_g <- crossprod(C_si, g_other)
            for (jj in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, jj] *
                      mu_q_b_specific[[s]][jj, l] *
                      mu_q_b_specific[[s]][jj, l_other]
              if (abs(coef) < 1e-15) next
              sum_mu <- sum_mu - coef * as.vector(ct_g)
            }
          }
        }

        # ---- Posterior mean ----
        mu_q_xi_new[[s]][[l]][i, ] <- as.vector(
          c_val * Sigma_q_xi_new[[s]][[l]][[i]] %*%
          crossprod(mu_q_nu_psi[[s]][[l]], sum_mu)
        )

        # ---- Trace for loadings ----
        tr_xi_specific[s, l, i] <- tr(
          H_psi %*% (Sigma_q_xi_new[[s]][[l]][[i]] +
                      tcrossprod(mu_q_xi_new[[s]][[l]][i, ]))
        )
      }
    }
  }

  create_named_list(mu_q_xi = mu_q_xi_new,
                    Sigma_q_xi = Sigma_q_xi_new,
                    tr_xi_specific = tr_xi_specific)
}
