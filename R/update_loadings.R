# =============================================================================
# update_loadings.R — Loadings a_{jl} and b_{sjl} with Spike-and-Slab PPI
#
# Covers parameter blocks 4.7 (shared a) and 4.8 (specific b).
#
# CRITICAL: PPI formula uses correct temperature scaling:
#   logit(PPI) = c * (E[log ω] - E[log(1-ω)])
#               + 0.5 * (mu^2/sigma^2 + log(sigma^2))
#   The Gaussian entropy term is NOT multiplied by c!
#
# Shared loadings a: cross-study aggregation
# Specific loadings b: within study only
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update shared loadings a_{jl} and PPI gamma_{jl}
#'
#' Cross-study aggregation: uses all studies' data.
#'
#' @param tr_qi_shared L_f x S x max(n_s) array of E[||f^{(l)}_{si}||^2]
#'
#' @export
update_a_loadings <- function(Y, C, list_cp_C, list_cp_C_Y,
                               mu_q_nu_mu,
                               mu_q_nu_beta, Z,
                               mu_q_zeta, mu_q_nu_phi,
                               mu_q_xi, mu_q_nu_psi,
                               mu_q_a, term_a, mu_q_gamma_a,
                               mu_q_normal_a, Sigma_q_normal_a,
                               mu_q_b_specific,
                               mu_q_recip_sigsq_eps,
                               mu_q_log_omega_a, mu_q_log_1_omega_a,
                               tr_qi_shared,
                               S, n_s, p, d, L_f, L_s, K_total,
                               c_val = 1, n_cpus = 1) {

  # Precompute list_tcp_nu_phi_zeta: (K+2) x n_s[s] per study per factor
  list_tcp_nu_phi_zeta <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    list_tcp_nu_phi_zeta[[l]] <- vector("list", S)
    for (s in 1:S) {
      list_tcp_nu_phi_zeta[[l]][[s]] <-
        tcrossprod(mu_q_nu_phi[[l]], mu_q_zeta[[s]][[l]])
    }
  }

  # Precompute C_i^T C_i * mu_q_nu_mu per study
  list_cp_C_nu_mu <- vector("list", S)
  for (s in 1:S) {
    list_cp_C_nu_mu[[s]] <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(jj) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[jj]])
    })
  }

  # Precompute specific factor contributions for residual
  # g_factor_contrib[[s]][[l]][,i] = C_i * nu_psi_sl * xi_sli  (n_obs × 1)
  # We need C^T * g for residual, so precompute in (K+2) space
  g_ct_contrib <- vector("list", S)
  if (L_s > 0 && !is.null(mu_q_nu_psi)) {
    for (s in 1:S) {
      g_ct_contrib[[s]] <- lapply(1:n_s[s], function(i) {
        contrib <- rep(0, K_total)
        for (l in seq_len(L_s)) {
          g_ct <- crossprod(C[[s]][[i]], as.vector(
            C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ]))
          contrib <- contrib + mu_q_b_specific[[s]][, l] * as.vector(g_ct)
        }
        # contrib is K_total x p: column j = sum_l b_{sjl} * g^{(l)}_si contribution
        matrix(contrib, nrow = K_total, ncol = p)
      })
    }
  }

  # ---- Update per factor l (accumulate into matrices) ----
  mu_q_normal_a_new <- matrix(0, nrow = p, ncol = L_f)
  Sigma_q_normal_a_new <- matrix(1, nrow = p, ncol = L_f)
  mu_q_gamma_a_new <- matrix(0.5, nrow = p, ncol = L_f)
  mu_q_a_new <- matrix(0, nrow = p, ncol = L_f)
  term_a_new <- matrix(0, nrow = p, ncol = L_f)

  for (l in seq_len(L_f)) {

    # Aggregate E[||f||^2] across all individuals in all studies
    rs_tr_qi_l <- 0
    for (s in 1:S) {
      for (i in 1:n_s[s]) {
        rs_tr_qi_l <- rs_tr_qi_l + tr_qi_shared[l, s, i]
      }
    }

    # sum_mu_q per variable j
    sum_mu_a <- rep(0, p)

    for (s in 1:S) {
      eps_prec_j <- mu_q_recip_sigsq_eps[s, ]

      for (i in 1:n_s[s]) {

        tcp_phi_zeta_li <- list_tcp_nu_phi_zeta[[l]][[s]][, i]  # (K+2) × 1

        for (jj in 1:p) {

          # f̄^T * (C^T y - C^T C mu)
          term <- as.numeric(crossprod(tcp_phi_zeta_li,
            list_cp_C_Y[[s]][[i]][, jj] - list_cp_C_nu_mu[[s]][[i]][, jj]))

          # Subtract beta contribution
          if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
            for (r in 1:d) {
              ct_beta <- list_cp_C[[s]][[i]] %*% mu_q_nu_beta[[jj]][[r]]
              term <- term - Z[[s]][i, r] *
                as.numeric(crossprod(tcp_phi_zeta_li, as.vector(ct_beta)))
            }
          }

          # Subtract OTHER shared factors (l_other != l)
          if (L_f > 1) {
            for (l_other in setdiff(seq_len(L_f), l)) {
              ct_f_other <- list_cp_C[[s]][[i]] %*%
                list_tcp_nu_phi_zeta[[l_other]][[s]][, i]
              term <- term - mu_q_a[jj, l_other] *
                as.numeric(crossprod(tcp_phi_zeta_li, as.vector(ct_f_other)))
            }
          }

          # Subtract specific factors
          if (L_s > 0 && !is.null(g_ct_contrib[[s]])) {
            term <- term - as.numeric(
              crossprod(tcp_phi_zeta_li, g_ct_contrib[[s]][[i]][, jj]))
          }

          sum_mu_a[jj] <- sum_mu_a[jj] + eps_prec_j[jj] * term
        }
      }
    }

    # ---- Slab variance (per-variable, cross-study aggregation) ----
    for (jj in 1:p) {
      denom <- 1 + sum(sapply(1:S, function(ss) {
        mu_q_recip_sigsq_eps[ss, jj] *
          sum(tr_qi_shared[l, ss, 1:n_s[ss]])
      }))
      Sigma_q_normal_a_new[jj, l] <- 1 / (c_val * denom)
    }

    # ---- Slab mean ----
    mu_q_normal_a_new[, l] <- c_val * Sigma_q_normal_a_new[, l] * sum_mu_a

    # ---- PPI (CRITICAL: entropy term NOT multiplied by c) ----
    log_odds_prior <- c_val * (mu_q_log_omega_a[l] - mu_q_log_1_omega_a[l])
    entropy_term <- 0.5 * (mu_q_normal_a_new[, l]^2 / Sigma_q_normal_a_new[, l] +
                           log(Sigma_q_normal_a_new[, l]))

    mu_q_gamma_a_new[, l] <- 1 / (1 + exp(-(log_odds_prior + entropy_term)))

    # ---- Spike-and-slab moments ----
    mu_q_a_new[, l] <- mu_q_gamma_a_new[, l] * mu_q_normal_a_new[, l]
    term_a_new[, l] <- (Sigma_q_normal_a_new[, l] + mu_q_normal_a_new[, l]^2) *
                       mu_q_gamma_a_new[, l]
  }

  create_named_list(mu_q_normal_a = mu_q_normal_a_new,
                    Sigma_q_normal_a = Sigma_q_normal_a_new,
                    mu_q_gamma_a = mu_q_gamma_a_new,
                    mu_q_a = mu_q_a_new,
                    term_a = term_a_new)
}


#' Update study-specific loadings b_{sjl} and PPI
#'
#' Within-study only: each study's loadings are independent.
#'
#' @param tr_xi_specific S x L_s x max(n_s) array of E[||g^{(l)}_{si}||^2]
#'
#' @export
update_b_loadings <- function(Y, C, list_cp_C, list_cp_C_Y,
                               mu_q_nu_mu,
                               mu_q_nu_beta, Z,
                               mu_q_zeta, mu_q_nu_phi,
                               mu_q_xi, mu_q_nu_psi,
                               mu_q_a,
                               mu_q_b_specific, term_b_specific, mu_q_gamma_b,
                               mu_q_normal_b_specific, Sigma_q_normal_b,
                               mu_q_recip_sigsq_eps,
                               mu_q_log_omega_b, mu_q_log_1_omega_b,
                               tr_xi_specific,
                               S, n_s, p, d, L_f, L_s, K_total,
                               c_val = 1, n_cpus = 1) {

  if (L_s == 0 || is.null(mu_q_nu_psi)) return(NULL)

  mu_q_normal_b_new <- Sigma_q_normal_b_new <- vector("list", S)
  mu_q_gamma_b_new <- mu_q_b_spec_new <- term_b_spec_new <- vector("list", S)

  for (s in 1:S) {

    mu_q_normal_b_new[[s]] <- matrix(0, p, L_s)
    Sigma_q_normal_b_new[[s]] <- matrix(1, p, L_s)
    mu_q_gamma_b_new[[s]] <- matrix(0.5, p, L_s)

    # Precompute within study
    list_cp_C_nu_mu_s <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(jj) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[jj]])
    })

    # Precompute shared factor contributions for residual
    f_ct_contrib_s <- lapply(1:n_s[s], function(i) {
      contrib <- rep(0, K_total)
      for (l in seq_len(L_f)) {
        f_ct <- crossprod(C[[s]][[i]], as.vector(
          C[[s]][[i]] %*% mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ]))
        contrib <- contrib + mu_q_a[, l] * as.vector(f_ct)
      }
      matrix(contrib, nrow = K_total, ncol = p)
    })

    for (l in seq_len(L_s)) {

      # Aggregate E[||g||^2] within study s
      rs_tr_xi_sl <- 0
      for (i in 1:n_s[s]) {
        rs_tr_xi_sl <- rs_tr_xi_sl + tr_xi_specific[s, l, i]
      }

      tcp_psi_xi_sl <- tcrossprod(mu_q_nu_psi[[s]][[l]], mu_q_xi[[s]][[l]])

      sum_mu_b <- rep(0, p)
      eps_prec_j <- mu_q_recip_sigsq_eps[s, ]

      for (i in 1:n_s[s]) {

        tcp_psi_xi_li <- tcp_psi_xi_sl[, i]

        for (jj in 1:p) {

          term <- as.numeric(crossprod(tcp_psi_xi_li,
            list_cp_C_Y[[s]][[i]][, jj] - list_cp_C_nu_mu_s[[i]][, jj]))

          # Subtract beta
          if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
            for (r in 1:d) {
              ct_beta <- list_cp_C[[s]][[i]] %*% mu_q_nu_beta[[jj]][[r]]
              term <- term - Z[[s]][i, r] *
                as.numeric(crossprod(tcp_psi_xi_li, as.vector(ct_beta)))
            }
          }

          # Subtract ALL shared factors
          if (L_f > 0) {
            term <- term - as.numeric(
              crossprod(tcp_psi_xi_li, f_ct_contrib_s[[i]][, jj]))
          }

          # Subtract OTHER specific factors
          if (L_s > 1) {
            for (l_other in setdiff(seq_len(L_s), l)) {
              ct_g_other <- list_cp_C[[s]][[i]] %*%
                tcrossprod(mu_q_nu_psi[[s]][[l_other]],
                           mu_q_xi[[s]][[l_other]])[, i]
              term <- term - mu_q_b_specific[[s]][jj, l_other] *
                as.numeric(crossprod(tcp_psi_xi_li, as.vector(ct_g_other)))
            }
          }

          sum_mu_b[jj] <- sum_mu_b[jj] + eps_prec_j[jj] * term
        }
      }

      # ---- Slab variance and mean ----
      Sigma_q_normal_b_new[[s]][, l] <-
        1 / (c_val * (eps_prec_j * rs_tr_xi_sl + 1))
      mu_q_normal_b_new[[s]][, l] <-
        c_val * Sigma_q_normal_b_new[[s]][, l] * sum_mu_b

      # ---- PPI (entropy term NOT multiplied by c) ----
      log_odds_prior <- c_val *
        (mu_q_log_omega_b[[s]][l] - mu_q_log_1_omega_b[[s]][l])
      entropy_term <- 0.5 *
        (mu_q_normal_b_new[[s]][, l]^2 / Sigma_q_normal_b_new[[s]][, l] +
         log(Sigma_q_normal_b_new[[s]][, l]))

      mu_q_gamma_b_new[[s]][, l] <-
        1 / (1 + exp(-(log_odds_prior + entropy_term)))
    }

    mu_q_b_spec_new[[s]] <- mu_q_gamma_b_new[[s]] * mu_q_normal_b_new[[s]]
    term_b_spec_new[[s]] <-
      (Sigma_q_normal_b_new[[s]] + mu_q_normal_b_new[[s]]^2) *
      mu_q_gamma_b_new[[s]]
  }

  create_named_list(mu_q_normal_b = mu_q_normal_b_new,
                    Sigma_q_normal_b = Sigma_q_normal_b_new,
                    mu_q_gamma_b = mu_q_gamma_b_new,
                    mu_q_b_specific = mu_q_b_spec_new,
                    term_b_specific = term_b_spec_new)
}
