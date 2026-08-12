# =============================================================================
# update_xi.R — Study-specific standardised computational scores chi^{(sl)}_{si}
#
# Covers parameter block 4.4 (derivations.md §3.4).
# The function/object name `xi` is retained only for backward compatibility.
# Statistically, every xi object in this file represents chi with unit prior
# variance.  It becomes the theoretical FPCA score xi only after post-processing.
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

#' Update study-specific standardised computational scores chi
#'
#' The exported function name and argument names retain `xi` for backward
#' compatibility; these objects are working scores chi, not canonical FPCA
#' scores xi.
#'
#' @noRd
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

  if (!.has_specific(L_s) || is.null(mu_q_nu_psi)) return(NULL)

  mu_q_xi_new <- vector("list", S)
  Sigma_q_xi_new <- vector("list", S)
  tr_xi_specific <- array(NA, dim = c(S, .max_L_s(L_s), max(n_s)))
  mu_q_xi_work <- mu_q_xi
  trace_driver <- .driver_diagnostic_should_trace()

  for (s in 1:S) {
    L_ss <- .L_s_at(L_s, s)

    mu_q_xi_new[[s]] <- vector("list", L_ss)
    Sigma_q_xi_new[[s]] <- vector("list", L_ss)
    if (L_ss == 0L) next

    # Precompute C_i^T C_i * mu_q_nu_mu_j
    list_cp_C_nu_mu <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(j) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[j]])
    })

    for (l in seq_len(L_ss)) {

      M_sl <- M_s[[s]][l]

      # Prior precision: I_{M_sl}
      inv_Sigma_xi <- diag(M_sl)

      mu_q_xi_new[[s]][[l]] <- matrix(0, nrow = n_s[s], ncol = M_sl)
      Sigma_q_xi_new[[s]][[l]] <- vector("list", n_s[s])
      if (trace_driver) {
        data_norm_squared <- linear_norm_squared <- mean_norm_squared <- 0
        data_envelope <- precision_trace <- precision_norm_squared <- 0
      }

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
        Sigma_q_xi_new[[s]][[l]][[i]] <- .inverse_spd(
          prec_xi,
          context = sprintf("xi[s=%d,l=%d,i=%d]", s, l, i))

        # ---- Linear term ----
        sum_mu <- rep(0, ncol(C_si))
        if (trace_driver) {
          data_sum_mu <- rep(0, ncol(C_si))
          data_envelope_i <- 0
        }

        # Positive: sum_j b_{sjl} * sigma^{-2}_eps * (C^T y - C^T C nu_mu)
        for (jj in 1:p) {
          coef <- mu_q_recip_sigsq_eps[s, jj] * mu_q_b_specific[[s]][jj, l]
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
        if (L_ss > 1L) {
          for (l_other in setdiff(seq_len(L_ss), l)) {
            g_other <- as.vector(
              C_si %*% mu_q_nu_psi[[s]][[l_other]] %*%
                mu_q_xi_work[[s]][[l_other]][i, ]
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
        score_driver <- crossprod(mu_q_nu_psi[[s]][[l]], sum_mu)
        mu_q_xi_new[[s]][[l]][i, ] <- as.vector(
          c_val * Sigma_q_xi_new[[s]][[l]][[i]] %*%
          score_driver
        )
        mu_q_xi_work[[s]][[l]][i, ] <-
          mu_q_xi_new[[s]][[l]][i, ]

        # ---- Trace for loadings ----
        tr_xi_specific[s, l, i] <- tr(
          H_psi %*% (Sigma_q_xi_new[[s]][[l]][[i]] +
                      tcrossprod(mu_q_xi_new[[s]][[l]][i, ]))
        )
        if (trace_driver) {
          data_norm_squared <- data_norm_squared +
            .driver_diagnostic_norm(data_sum_mu)^2
          data_envelope <- data_envelope + data_envelope_i
          linear_norm_squared <- linear_norm_squared +
            .driver_diagnostic_norm(score_driver)^2
          mean_norm_squared <- mean_norm_squared +
            .driver_diagnostic_norm(mu_q_xi_new[[s]][[l]][i, ])^2
          precision_trace <- precision_trace + sum(diag(prec_xi))
          precision_norm_squared <- precision_norm_squared +
            .driver_diagnostic_norm(prec_xi)^2
        }
      }
      if (trace_driver) {
        data_norm <- sqrt(data_norm_squared)
        precision_norm <- sqrt(precision_norm_squared)
        linear_norm <- sqrt(linear_norm_squared)
        .driver_diagnostic_record(
          block = "specific_score", study = s, factor = l,
          metrics = list(
            data_driver_norm = data_norm,
            data_envelope = data_envelope,
            data_coherence = data_norm / max(data_envelope, .Machine$double.eps),
            linear_driver_norm = linear_norm,
            precision_trace = precision_trace,
            precision_frobenius = precision_norm,
            prior_precision_trace = c_val * n_s[s] * M_sl,
            driver_precision_ratio =
              linear_norm / max(precision_norm, .Machine$double.eps),
            posterior_mean_norm = sqrt(mean_norm_squared)
          )
        )
      }
    }
  }

  create_named_list(mu_q_xi = mu_q_xi_work,
                    Sigma_q_xi = Sigma_q_xi_new,
                    tr_xi_specific = tr_xi_specific)
}
