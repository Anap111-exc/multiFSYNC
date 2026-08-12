# =============================================================================
# update_zeta.R — Shared standardised computational scores eta^{(l)}_{si}
#
# Covers parameter block 4.3 (derivations.md §3.3).
# The function/object name `zeta` is retained only for backward compatibility.
# Statistically, every zeta object in this file represents eta with unit prior
# variance.  It becomes the theoretical FPCA score zeta only after the joint
# score/time-function post-processing in orthonormalise_multi().
#
# (Sigma^q_zeta)^{-1} = c * (sum_sigma * H_phi + I_{M_f[l]})
# mu^q_zeta = c * Sigma_zeta * crossprod(nu_phi, sum_mu)
#
# H_phi = nu_phi^T * C^T C * nu_phi + diag(tr(C^T C * Sigma_nu_phi_m))
# sum_sigma = sum_j sigma^{-2}_eps[s,j] * E[a^2_{jl}]
#
# Residual r^zeta: excludes current shared factor l,
# includes all specific factors and all other shared factors.
#
# Indexing: study [[s]] outermost
#
# Symbol convention:
#   L_f  number of shared factors
#   L_s  number of specific factors
#   M_f  vector, M_f[l] = number of FPCA basis functions for shared factor l
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update shared standardised computational scores eta
#'
#' The exported function name and argument names retain `zeta` for backward
#' compatibility; these objects are working scores eta, not canonical FPCA
#' scores zeta.
#'
#' @noRd
update_zeta <- function(Y, C, list_cp_C, list_cp_C_Y,
                         mu_q_nu_mu,
                         mu_q_nu_beta, Z,
                         mu_q_zeta, Sigma_q_zeta,
                         mu_q_nu_phi, Sigma_q_nu_phi,
                         mu_q_xi, mu_q_nu_psi,
                         mu_q_a, term_a,
                         mu_q_b_specific,
                         mu_q_recip_sigsq_eps,
                         S, n_s, p, d, L_f, L_s, M_f,
                         c_val = 1, n_cpus = 1) {

  # ---- Initialize output ----
  mu_q_zeta_new <- vector("list", S)
  Sigma_q_zeta_new <- vector("list", S)
  tr_qi_shared <- array(NA, dim = c(L_f, S, max(n_s)))
  mu_q_zeta_work <- mu_q_zeta
  trace_driver <- .driver_diagnostic_should_trace()

  for (s in 1:S) {
    L_ss <- .L_s_at(L_s, s)

    mu_q_zeta_new[[s]] <- vector("list", L_f)
    Sigma_q_zeta_new[[s]] <- vector("list", L_f)

    # Precompute C_i^T C_i * mu_q_nu_mu_j for all (i,j) in this study
    list_cp_C_nu_mu <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(j) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[j]])
    })

    for (l in seq_len(L_f)) {

      M_l <- M_f[l]

      # Prior precision (M_l x M_l identity)
      inv_Sigma_zeta <- diag(M_l)

      mu_q_zeta_new[[s]][[l]] <- matrix(0, nrow = n_s[s], ncol = M_l)
      Sigma_q_zeta_new[[s]][[l]] <- vector("list", n_s[s])
      if (trace_driver) {
        data_norm_squared <- linear_norm_squared <- mean_norm_squared <- 0
        data_envelope <- precision_trace <- precision_norm_squared <- 0
      }

      # Scalar precision factor for this (s,l)
      sum_sigma_a_l <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l])

      for (i in 1:n_s[s]) {

        C_si <- C[[s]][[i]]
        cp_C_si <- list_cp_C[[s]][[i]]

        # ---- H_phi matrix (M_l x M_l) ----
        # H = nu_phi^T * C^T C * nu_phi + diag(tr(C^T C * Sigma_nu_phi_m))
        H_mean <- crossprod(mu_q_nu_phi[[l]], cp_C_si) %*% mu_q_nu_phi[[l]]
        H_var <- diag(sapply(1:M_l, function(m) {
          tr(cp_C_si %*% Sigma_q_nu_phi[[l]][[m]])
        }), nrow = M_l, ncol = M_l)
        H_phi <- H_mean + H_var

        # ---- Posterior precision and covariance ----
        prec_zeta <- c_val * (sum_sigma_a_l * H_phi + inv_Sigma_zeta)
        Sigma_q_zeta_new[[s]][[l]][[i]] <- .inverse_spd(
          prec_zeta,
          context = sprintf("zeta[s=%d,l=%d,i=%d]", s, l, i))

        # ---- Linear term (residual excluding current factor l) ----
        sum_mu <- rep(0, ncol(C_si))  # length K_total
        if (trace_driver) {
          data_sum_mu <- rep(0, ncol(C_si))
          data_envelope_i <- 0
        }

        # Positive term: sum_j a_jl * sigma^{-2}_eps * (C^T y - C^T C nu_mu)
        for (jj in 1:p) {
          coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_a[jj, l]
          if (abs(coef) < 1e-15) next
          ct_res <- list_cp_C_Y[[s]][[i]][, jj] - list_cp_C_nu_mu[[i]][, jj]
          sum_mu <- sum_mu + coef * as.vector(ct_res)
          if (trace_driver) {
            data_delta <- coef * as.vector(ct_res)
            data_sum_mu <- data_sum_mu + data_delta
            data_envelope_i <- data_envelope_i +
              .driver_diagnostic_norm(data_delta)
          }
        }

        # Subtract beta contribution
        if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
          for (jj in 1:p) {
            coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_a[jj, l]
            if (abs(coef) < 1e-15) next
            for (r in 1:d) {
              ct_beta <- cp_C_si %*% mu_q_nu_beta[[jj]][[r]]
              sum_mu <- sum_mu - coef * Z[[s]][i, r] * as.vector(ct_beta)
            }
          }
        }

        # Subtract OTHER shared factor contributions (q != l)
        if (L_f > 1) {
          for (l_other in setdiff(1:L_f, l)) {
            f_other <- as.vector(
              C_si %*% mu_q_nu_phi[[l_other]] %*%
                mu_q_zeta_work[[s]][[l_other]][i, ]
            )
            ct_f <- crossprod(C_si, f_other)
            for (jj in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_a[jj, l] *
                      mu_q_a[jj, l_other]
              if (abs(coef) < 1e-15) next
              sum_mu <- sum_mu - coef * as.vector(ct_f)
            }
          }
        }

        # Subtract ALL specific factor contributions
        if (L_ss > 0L && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
          for (l_spec in seq_len(L_ss)) {
            g_other <- as.vector(
              C_si %*% mu_q_nu_psi[[s]][[l_spec]] %*% mu_q_xi[[s]][[l_spec]][i, ]
            )
            ct_g <- crossprod(C_si, g_other)
            for (jj in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_a[jj, l] *
                      mu_q_b_specific[[s]][jj, l_spec]
              if (abs(coef) < 1e-15) next
              sum_mu <- sum_mu - coef * as.vector(ct_g)
            }
          }
        }

        # ---- Posterior mean: project to M_l-dim via nu_phi^T ----
        score_driver <- crossprod(mu_q_nu_phi[[l]], sum_mu)
        mu_q_zeta_new[[s]][[l]][i, ] <- as.vector(
          c_val * Sigma_q_zeta_new[[s]][[l]][[i]] %*%
          score_driver
        )
        mu_q_zeta_work[[s]][[l]][i, ] <-
          mu_q_zeta_new[[s]][[l]][i, ]

        # ---- tr_qi for loadings (trace of H_phi * (Sigma + mu mu^T)) ----
        tr_qi_shared[l, s, i] <- tr(
          H_phi %*% (Sigma_q_zeta_new[[s]][[l]][[i]] +
                      tcrossprod(mu_q_zeta_new[[s]][[l]][i, ]))
        )
        if (trace_driver) {
          data_norm_squared <- data_norm_squared +
            .driver_diagnostic_norm(data_sum_mu)^2
          data_envelope <- data_envelope + data_envelope_i
          linear_norm_squared <- linear_norm_squared +
            .driver_diagnostic_norm(score_driver)^2
          mean_norm_squared <- mean_norm_squared +
            .driver_diagnostic_norm(mu_q_zeta_new[[s]][[l]][i, ])^2
          precision_trace <- precision_trace + sum(diag(prec_zeta))
          precision_norm_squared <- precision_norm_squared +
            .driver_diagnostic_norm(prec_zeta)^2
        }
      }
      if (trace_driver) {
        data_norm <- sqrt(data_norm_squared)
        precision_norm <- sqrt(precision_norm_squared)
        linear_norm <- sqrt(linear_norm_squared)
        .driver_diagnostic_record(
          block = "shared_score", study = s, factor = l,
          metrics = list(
            data_driver_norm = data_norm,
            data_envelope = data_envelope,
            data_coherence = data_norm / max(data_envelope, .Machine$double.eps),
            linear_driver_norm = linear_norm,
            precision_trace = precision_trace,
            precision_frobenius = precision_norm,
            prior_precision_trace = c_val * n_s[s] * M_l,
            driver_precision_ratio =
              linear_norm / max(precision_norm, .Machine$double.eps),
            posterior_mean_norm = sqrt(mean_norm_squared)
          )
        )
      }
    }
  }

  create_named_list(mu_q_zeta = mu_q_zeta_work,
                    Sigma_q_zeta = Sigma_q_zeta_new,
                    tr_qi_shared = tr_qi_shared)
}
