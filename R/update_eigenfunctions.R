# =============================================================================
# update_eigenfunctions.R — Eigenfunction spline coefficients update
#
# Covers parameter blocks:
#   4.5  nu_{phi,ml}  — shared eigenfunction coefficients
#   4.6  nu_{psi,slm} — study-specific eigenfunction coefficients
#
# Key formulas (derivations.md §3.5-3.6):
#   (Sigma^q)^{-1} = c * (sum_vec * list_sum + inv_prior)
#   mu^q            = c * Sigma^q * sum_term
#
#   sum_vec:  Σ_{s,j} sigma^{-2}_eps * E[a^2]   (shared)
#             Σ_j sigma^{-2}_eps * E[b^2]        (specific, per study)
#   list_sum: Σ_{s,i} E[zeta^2] * C_i^T C_i     (weighted sum of C^T C)
#
# Residual assembly (Appendix C.4): must exclude
#   - current factor l (cross-factor terms q≠l)
#   - current component m (same-factor other-component m'≠m)
#   - beta / specific factor contributions
#
# Indexing: study [[s]] outermost
#
# Notation (after rename):
#   L_f    : number of shared factors
#   L_s    : number of study-specific factors
#   M_f[l] : number of components (spline basis) for shared factor l
#   M_s[[s]][l] : number of components for study-s factor l
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update shared eigenfunction coefficients nu_{phi,ml}
#'
#' Cross-study aggregation: uses data from ALL studies to update
#' shared eigenfunctions.  Each shared factor l can have a different
#' number of components M_f[l].
#'
#' @export
update_nu_phi <- function(Y, C, list_cp_C, list_cp_C_Y,
                           mu_q_nu_mu, Sigma_q_nu_mu,
                           mu_q_nu_beta, Z,
                           mu_q_zeta, Sigma_q_zeta, mu_q_nu_phi,
                           mu_q_xi, mu_q_nu_psi,
                           mu_q_a, term_a,
                           mu_q_b_specific,
                           term_b_specific = NULL,
                           mu_q_recip_sigsq_eps, mu_q_recip_sigsq_phi,
                           inv_Sigma_beta,
                           S, n_s, p, d, L_f, L_s, M_f, K, K_total,
                           c_val = 1, n_cpus = 1) {

  # ---- Step 2: Dimension assertions ----
  stopifnot(length(mu_q_nu_phi) == L_f)
  # per-factor component counts may differ, no longer assert uniform ncol

  mu_q_nu_phi_new <- vector("list", L_f)
  Sigma_q_nu_phi_new <- vector("list", L_f)
  inv_Sigma_q_nu_phi_new <- vector("list", L_f)

  for (l in seq_len(L_f)) {

    M_l <- M_f[l]

    mu_q_nu_phi_new[[l]] <- matrix(0, nrow = K_total, ncol = M_l)
    Sigma_q_nu_phi_new[[l]] <- vector("list", M_l)
    inv_Sigma_q_nu_phi_new[[l]] <- vector("list", M_l)

    for (m in 1:M_l) {

      # --- Prior precision ---
      inv_prior <- blkdiag(inv_Sigma_beta,
                           mu_q_recip_sigsq_phi[[l]][m] * diag(K))

      # --- Per-study weighted sum: Σ_s w_s * Σ_i E[zeta^2] * C_i^T C_i ---
      # w_s = Σ_j sigma^{-2}_{eps,sj} * E[a²_{jl}] (study-specific scalar)
      prec_data <- matrix(0, K_total, K_total)
      for (s in 1:S) {
        w_s <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l])
        for (i in 1:n_s[s]) {
          Ez2 <- Sigma_q_zeta[[s]][[l]][[i]][m, m] +
                 mu_q_zeta[[s]][[l]][i, m]^2
          prec_data <- prec_data + w_s * Ez2 * list_cp_C[[s]][[i]]
        }
      }

      # --- Soft orthogonal penalty: penalize overlap with other shared factors ---
      if (L_f > 1) {
        lambda_orth <- 0.1  # penalty strength
        for (l_other in setdiff(seq_len(L_f), l)) {
          w_other <- sum(mu_q_recip_sigsq_eps[, 1] * term_a[, l_other])  # rough scalar per study
          for (s in 1:S) {
            w_s_other <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l_other])
            for (i in 1:n_s[s]) {
              f_other <- as.vector(C[[s]][[i]] %*% mu_q_nu_phi[[l_other]] %*% mu_q_zeta[[s]][[l_other]][i, ])
              f_ct <- crossprod(C[[s]][[i]], f_other)
              prec_data <- prec_data + lambda_orth * w_s_other * tcrossprod(f_ct)
            }
          }
        }
      }

      # --- 3) Shared vs specific (prevent shared from absorbing specific signal) ---
      if (L_s > 0 && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
        for (s in 1:S) {
          for (l_spec in seq_len(L_s)) {
            w_s_spec <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l_spec])
            for (i in 1:n_s[s]) {
              g_spec <- as.vector(C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l_spec]] %*% mu_q_xi[[s]][[l_spec]][i, ])
              g_ct <- crossprod(C[[s]][[i]], g_spec)
              prec_data <- prec_data + lambda_orth * w_s_spec * tcrossprod(g_ct)
            }
          }
        }
      }

      # --- Posterior precision ---
      prec <- c_val * (prec_data + inv_prior)
      Sigma_q_nu_phi_new[[l]][[m]] <- solve(prec + 1e-8 * diag(K_total))
      inv_Sigma_q_nu_phi_new[[l]][[m]] <- inv_prior

      # --- Linear term (Appendix C.4): cross-study aggregation ---
      sum_term <- rep(0, K_total)

      for (s in 1:S) {
        for (i in 1:n_s[s]) {

          zeta_silm <- mu_q_zeta[[s]][[l]][i, m]
          C_si <- C[[s]][[i]]
          cp_C_si <- list_cp_C[[s]][[i]]

          # ---- Positive data term (minus mu contribution) for all j ----
          for (j in 1:p) {
            coef <- mu_q_recip_sigsq_eps[s, j] * mu_q_a[j, l] * zeta_silm
            if (abs(coef) < 1e-15) next

            # C_i^T * y_{sij} - C_i^T C_i * nu_mu_{sj}
            ct_y <- list_cp_C_Y[[s]][[i]][, j]
            ct_mu <- cp_C_si %*% mu_q_nu_mu[[s]][[j]]
            sum_term <- sum_term + coef * as.vector(ct_y - ct_mu)
          }

          # ---- Subtract beta contribution ----
          if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
            for (j in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, j] * mu_q_a[j, l] * zeta_silm
              if (abs(coef) < 1e-15) next
              for (r in 1:d) {
                ct_beta <- cp_C_si %*% mu_q_nu_beta[[j]][[r]]
                sum_term <- sum_term - coef * Z[[s]][i, r] * as.vector(ct_beta)
              }
            }
          }

          # ---- Subtract same-factor other-components (m' != m) ----
          if (M_l > 1) {
            for (m_tilde in setdiff(1:M_l, m)) {
              zeta_silm_t <- mu_q_zeta[[s]][[l]][i, m_tilde]
              cov_m_mt <- Sigma_q_zeta[[s]][[l]][[i]][m, m_tilde]
              Ez_prod <- zeta_silm * zeta_silm_t + cov_m_mt

              if (abs(Ez_prod) < 1e-15) next

              ct_phi_mt <- cp_C_si %*% mu_q_nu_phi[[l]][, m_tilde]

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] * term_a[j, l] * Ez_prod
                sum_term <- sum_term - coef * as.vector(ct_phi_mt)
              }
            }
          }

          # ---- Subtract other shared factors (l_tilde != l) ----
          if (L_f > 1) {
            for (l_tilde in setdiff(1:L_f, l)) {
              f_other <- as.vector(
                C_si %*% mu_q_nu_phi[[l_tilde]] %*%
                mu_q_zeta[[s]][[l_tilde]][i, ]
              )
              ct_f_other <- crossprod(C_si, f_other)

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        mu_q_a[j, l] * mu_q_a[j, l_tilde] * zeta_silm
                if (abs(coef) < 1e-15) next
                sum_term <- sum_term - coef * as.vector(ct_f_other)
              }
            }
          }

          # ---- Subtract all specific factors ----
          if (L_s > 0 && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
            for (l_spec in 1:L_s) {
              g_other <- as.vector(
                C_si %*% mu_q_nu_psi[[s]][[l_spec]] %*%
                mu_q_xi[[s]][[l_spec]][i, ]
              )
              ct_g_other <- crossprod(C_si, g_other)

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        mu_q_a[j, l] * mu_q_b_specific[[s]][j, l_spec] *
                        zeta_silm
                if (abs(coef) < 1e-15) next
                sum_term <- sum_term - coef * as.vector(ct_g_other)
              }
            }
          }
        }
      }

      # --- Posterior mean ---
      mu_q_nu_phi_new[[l]][, m] <- as.vector(
        c_val * Sigma_q_nu_phi_new[[l]][[m]] %*% sum_term
      )

    } # end m loop
  } # end l loop

  create_named_list(mu_q_nu_phi = mu_q_nu_phi_new,
                    Sigma_q_nu_phi = Sigma_q_nu_phi_new,
                    inv_Sigma_q_nu_phi = inv_Sigma_q_nu_phi_new)
}


#' Update study-specific eigenfunction coefficients nu_{psi,slm}
#'
#' Within-study only: each study's eigenfunctions are updated
#' independently using only that study's data.  The number of components
#' can vary by study and factor: M_s[[s]][l].
#'
#' @export
update_nu_psi <- function(Y, C, list_cp_C, list_cp_C_Y,
                           mu_q_nu_mu, Sigma_q_nu_mu,
                           mu_q_nu_beta, Z,
                           mu_q_zeta, mu_q_nu_phi,
                           mu_q_xi, Sigma_q_xi, mu_q_nu_psi,
                           mu_q_a,
                           term_a = NULL,
                           mu_q_b_specific, term_b_specific,
                           mu_q_recip_sigsq_eps, mu_q_recip_sigsq_psi,
                           inv_Sigma_beta,
                           S, n_s, p, d, L_f, L_s, M_s, K, K_total,
                           c_val = 1, n_cpus = 1) {

  if (is.null(mu_q_nu_psi) || L_s == 0) return(NULL)

  mu_q_nu_psi_new <- vector("list", S)
  Sigma_q_nu_psi_new <- vector("list", S)
  inv_Sigma_q_nu_psi_new <- vector("list", S)

  for (s in 1:S) {

    mu_q_nu_psi_new[[s]] <- vector("list", L_s)
    Sigma_q_nu_psi_new[[s]] <- vector("list", L_s)
    inv_Sigma_q_nu_psi_new[[s]] <- vector("list", L_s)

    for (l in 1:L_s) {

      M_sl <- M_s[[s]][l]

      mu_q_nu_psi_new[[s]][[l]] <- matrix(0, nrow = K_total, ncol = M_sl)
      Sigma_q_nu_psi_new[[s]][[l]] <- vector("list", M_sl)
      inv_Sigma_q_nu_psi_new[[s]][[l]] <- vector("list", M_sl)

      # --- Scalar precision: Σ_j sigma^{-2}_eps * term_b (within study s) ---
      sum_vec_b_sl <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l])

      for (m in 1:M_sl) {

        # --- Prior precision ---
        inv_prior <- blkdiag(inv_Sigma_beta,
                             mu_q_recip_sigsq_psi[[s]][[l]][m] * diag(K))

        # --- Weighted C^T C: Σ_i E[xi^2_{sim}] * C_i^T C_i ---
        list_sum_psi_slm <- matrix(0, K_total, K_total)
        for (i in 1:n_s[s]) {
          Ex2 <- Sigma_q_xi[[s]][[l]][[i]][m, m] +
                 mu_q_xi[[s]][[l]][i, m]^2
          list_sum_psi_slm <- list_sum_psi_slm + Ex2 * list_cp_C[[s]][[i]]
        }

        # --- Soft orthogonal penalties ---
        lambda_orth <- 0.1
        
        # 1) Specific vs other specific factors in the same study
        if (L_s > 1) {
          for (l_other in setdiff(seq_len(L_s), l)) {
            w_s_other <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l_other])
            for (i in 1:n_s[s]) {
              g_other <- as.vector(C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l_other]] %*% mu_q_xi[[s]][[l_other]][i, ])
              g_ct <- crossprod(C[[s]][[i]], g_other)
              list_sum_psi_slm <- list_sum_psi_slm + lambda_orth * w_s_other * tcrossprod(g_ct)
            }
          }
        }
        
        # 2) Specific vs all shared factors (prevent shared from absorbing specific signal)
        if (L_f > 0) {
          for (l_shared in seq_len(L_f)) {
            w_sh <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l_shared])
            for (i in 1:n_s[s]) {
              f_sh <- as.vector(C[[s]][[i]] %*% mu_q_nu_phi[[l_shared]] %*% mu_q_zeta[[s]][[l_shared]][i, ])
              f_ct <- crossprod(C[[s]][[i]], f_sh)
              list_sum_psi_slm <- list_sum_psi_slm + lambda_orth * w_sh * tcrossprod(f_ct)
            }
          }
        }

        # --- Posterior precision ---
        prec <- c_val * (sum_vec_b_sl * list_sum_psi_slm + inv_prior)
        Sigma_q_nu_psi_new[[s]][[l]][[m]] <- solve(prec + 1e-8 * diag(K_total))
        inv_Sigma_q_nu_psi_new[[s]][[l]][[m]] <- inv_prior

        # --- Linear term: within-study s only ---
        sum_term <- rep(0, K_total)

        for (i in 1:n_s[s]) {

          xi_silm <- mu_q_xi[[s]][[l]][i, m]
          C_si <- C[[s]][[i]]
          cp_C_si <- list_cp_C[[s]][[i]]

          # ---- Positive data term (minus mu) ----
          for (j in 1:p) {
            coef <- mu_q_recip_sigsq_eps[s, j] *
                    mu_q_b_specific[[s]][j, l] * xi_silm
            if (abs(coef) < 1e-15) next

            ct_y <- list_cp_C_Y[[s]][[i]][, j]
            ct_mu <- cp_C_si %*% mu_q_nu_mu[[s]][[j]]
            sum_term <- sum_term + coef * as.vector(ct_y - ct_mu)
          }

          # ---- Subtract beta ----
          if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
            for (j in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, j] *
                      mu_q_b_specific[[s]][j, l] * xi_silm
              if (abs(coef) < 1e-15) next
              for (r in 1:d) {
                ct_beta <- cp_C_si %*% mu_q_nu_beta[[j]][[r]]
                sum_term <- sum_term - coef * Z[[s]][i, r] * as.vector(ct_beta)
              }
            }
          }

          # ---- Subtract same-factor other-components (m' != m) ----
          if (M_sl > 1) {
            for (m_tilde in setdiff(1:M_sl, m)) {
              xi_silm_t <- mu_q_xi[[s]][[l]][i, m_tilde]
              cov_m_mt <- Sigma_q_xi[[s]][[l]][[i]][m, m_tilde]
              Ex_prod <- xi_silm * xi_silm_t + cov_m_mt

              if (abs(Ex_prod) < 1e-15) next

              ct_psi_mt <- cp_C_si %*% mu_q_nu_psi[[s]][[l]][, m_tilde]

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        term_b_specific[[s]][j, l] * Ex_prod
                sum_term <- sum_term - coef * as.vector(ct_psi_mt)
              }
            }
          }

          # ---- Subtract other specific factors (l_tilde != l) ----
          if (L_s > 1) {
            for (l_tilde in setdiff(1:L_s, l)) {
              g_other <- as.vector(
                C_si %*% mu_q_nu_psi[[s]][[l_tilde]] %*%
                mu_q_xi[[s]][[l_tilde]][i, ]
              )
              ct_g_other <- crossprod(C_si, g_other)

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        mu_q_b_specific[[s]][j, l] *
                        mu_q_b_specific[[s]][j, l_tilde] *
                        xi_silm
                if (abs(coef) < 1e-15) next
                sum_term <- sum_term - coef * as.vector(ct_g_other)
              }
            }
          }

          # ---- Subtract all shared factors ----
          if (L_f > 0) {
            for (l_shared in 1:L_f) {
              f_other <- as.vector(
                C_si %*% mu_q_nu_phi[[l_shared]] %*%
                mu_q_zeta[[s]][[l_shared]][i, ]
              )
              ct_f_other <- crossprod(C_si, f_other)

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        mu_q_b_specific[[s]][j, l] *
                        mu_q_a[j, l_shared] *
                        xi_silm
                if (abs(coef) < 1e-15) next
                sum_term <- sum_term - coef * as.vector(ct_f_other)
              }
            }
          }
        }

        # --- Posterior mean ---
        mu_q_nu_psi_new[[s]][[l]][, m] <- as.vector(
          c_val * Sigma_q_nu_psi_new[[s]][[l]][[m]] %*% sum_term
        )
      }
    }
  }

  create_named_list(mu_q_nu_psi = mu_q_nu_psi_new,
                    Sigma_q_nu_psi = Sigma_q_nu_psi_new,
                    inv_Sigma_q_nu_psi = inv_Sigma_q_nu_psi_new)
}
