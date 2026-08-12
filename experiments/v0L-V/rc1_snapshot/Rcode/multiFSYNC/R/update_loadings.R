# =============================================================================
# update_loadings.R — Loadings a_{jl} and b_{sjl} with Spike-and-Slab PPI
#
# Covers parameter blocks 4.7 (shared a) and 4.8 (specific b).
#
# Under the fixed Lebesgue/counting-Dirac reference measure used by the
# tempered objective, the PPI logit is
#   c * (E[log ω] - E[log(1-ω)])
#   + 0.5 * (mu^2/sigma^2 + log(sigma^2))
#   + 0.5 * (1 - c) * log(2*pi).
# The last term is the continuous-slab normalising constant.  It vanishes at
# T = 1 but must not be dropped while T > 1.
#
# Shared loadings a: cross-study aggregation
# Specific loadings b: within study only
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

# Select E[log(omega)] - E[log(1-omega)] for one factor.  Factor-level
# omega is stored as a length-L vector, whereas variable-level omega is a
# p x L matrix.  Keeping this distinction here prevents R's vector recycling
# from silently using factor indices as variable indices.
.loading_log_prior_odds <- function(mu_log_omega, mu_log_1_omega,
                                    factor_index, p, label) {
  if (is.matrix(mu_log_omega) || is.matrix(mu_log_1_omega)) {
    if (!is.matrix(mu_log_omega) || !is.matrix(mu_log_1_omega) ||
        nrow(mu_log_omega) != p || nrow(mu_log_1_omega) != p ||
        ncol(mu_log_omega) < factor_index ||
        ncol(mu_log_1_omega) < factor_index) {
      stop(label, " variable-level omega moments must both be p x L matrices.")
    }
    return(mu_log_omega[, factor_index] -
             mu_log_1_omega[, factor_index])
  }

  if (length(mu_log_omega) < factor_index ||
      length(mu_log_1_omega) < factor_index) {
    stop(label, " factor-level omega moments must have length at least L.")
  }
  rep(mu_log_omega[factor_index] - mu_log_1_omega[factor_index], p)
}

# Exact inclusion logit for a standard-normal slab under
# L_T(q) = E_q[log p] - T E_q[log q], c = 1/T.  Conditional on inclusion,
# q(a | gamma = 1) is N(slab_mean, slab_variance).  The reference-measure
# term is required because the included state is continuous while the spike
# state is a point mass.
.tempered_spike_slab_logit <- function(log_prior_odds, slab_mean,
                                        slab_variance, c_val = 1) {
  if (length(c_val) != 1L || !is.finite(c_val) || c_val <= 0) {
    stop("c_val must be one finite positive number.")
  }
  if (any(!is.finite(slab_mean)) || any(!is.finite(slab_variance)) ||
      any(slab_variance <= 0) || any(!is.finite(log_prior_odds))) {
    stop("Spike-and-slab logit inputs must be finite and variances positive.")
  }

  c_val * log_prior_odds +
    0.5 * (slab_mean^2 / slab_variance + log(slab_variance)) +
    0.5 * (1 - c_val) * log(2 * pi)
}

#' Update shared loadings a_{jl} and PPI gamma_{jl}
#'
#' Cross-study aggregation: uses all studies' data.
#'
#' @param tr_qi_shared L_f x S x max(n_s) array of E[||f^{(l)}_{si}||^2]
#'
#' @noRd
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
  if (.has_specific(L_s) && !is.null(mu_q_nu_psi)) {
    for (s in 1:S) {
      L_ss <- .L_s_at(L_s, s)
      g_ct_contrib[[s]] <- lapply(1:n_s[s], function(i) {
        contrib <- matrix(0, nrow = K_total, ncol = p)
        for (l in seq_len(L_ss)) {
          g_ct <- crossprod(C[[s]][[i]], as.vector(
            C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l]] %*% mu_q_xi[[s]][[l]][i, ]))
          # Outer product is required: g_ct has length K_total while the
          # loading vector has length p.  Elementwise multiplication would
          # trigger R's recycling rule and corrupt the residual whenever
          # p != K_total.
          contrib <- contrib + tcrossprod(as.vector(g_ct),
                                           mu_q_b_specific[[s]][, l])
        }
        # contrib is K_total x p: column j = sum_l b_{sjl} * g^{(l)}_si contribution
        contrib
      })
    }
  }

  # ---- Update per factor l (accumulate into matrices) ----
  mu_q_normal_a_new <- matrix(0, nrow = p, ncol = L_f)
  Sigma_q_normal_a_new <- matrix(1, nrow = p, ncol = L_f)
  mu_q_gamma_a_new <- matrix(0.5, nrow = p, ncol = L_f)
  mu_q_a_new <- matrix(0, nrow = p, ncol = L_f)
  term_a_new <- matrix(0, nrow = p, ncol = L_f)
  mu_q_a_work <- mu_q_a
  term_a_work <- term_a
  trace_driver <- .driver_diagnostic_should_trace()

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
    if (trace_driver) data_envelope <- 0

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
              term <- term - mu_q_a_work[jj, l_other] *
                as.numeric(crossprod(tcp_phi_zeta_li, as.vector(ct_f_other)))
            }
          }

          # Subtract specific factors
          if (.L_s_at(L_s, s) > 0L && !is.null(g_ct_contrib[[s]])) {
            term <- term - as.numeric(
              crossprod(tcp_phi_zeta_li, g_ct_contrib[[s]][[i]][, jj]))
          }

          sum_mu_a[jj] <- sum_mu_a[jj] + eps_prec_j[jj] * term
          if (trace_driver) {
            data_envelope <- data_envelope + abs(eps_prec_j[jj] * term)
          }
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

    # ---- PPI from the joint tempered spike-and-slab normaliser ----
    log_prior_odds <- .loading_log_prior_odds(
      mu_q_log_omega_a, mu_q_log_1_omega_a, l, p, "Shared")
    logit_gamma <- .tempered_spike_slab_logit(
      log_prior_odds = log_prior_odds,
      slab_mean = mu_q_normal_a_new[, l],
      slab_variance = Sigma_q_normal_a_new[, l],
      c_val = c_val)

    unrestricted_gamma <- stats::plogis(logit_gamma)
    gate_active <- .driver_diagnostic_dense_gate_active()
    mu_q_gamma_a_new[, l] <- if (gate_active) 1 else unrestricted_gamma

    # ---- Spike-and-slab moments ----
    mu_q_a_new[, l] <- mu_q_gamma_a_new[, l] * mu_q_normal_a_new[, l]
    term_a_new[, l] <- (Sigma_q_normal_a_new[, l] + mu_q_normal_a_new[, l]^2) *
                       mu_q_gamma_a_new[, l]
    mu_q_a_work[, l] <- mu_q_a_new[, l]
    term_a_work[, l] <- term_a_new[, l]
    if (trace_driver) {
      driver_norm <- .driver_diagnostic_norm(sum_mu_a)
      precision_values <- 1 / Sigma_q_normal_a_new[, l]
      .driver_diagnostic_record(
        block = "shared_loading", factor = l,
        metrics = list(
          data_driver_norm = driver_norm,
          data_envelope = data_envelope,
          data_coherence =
            driver_norm / max(data_envelope, .Machine$double.eps),
          linear_driver_norm = driver_norm,
          precision_trace = sum(precision_values),
          precision_frobenius = .driver_diagnostic_norm(precision_values),
          prior_precision_trace = c_val * p,
          driver_precision_ratio =
            driver_norm /
            max(.driver_diagnostic_norm(precision_values),
                .Machine$double.eps),
          posterior_mean_norm = .driver_diagnostic_norm(mu_q_a_new[, l]),
          slab_mean_norm =
            .driver_diagnostic_norm(mu_q_normal_a_new[, l]),
          slab_variance_mean = mean(Sigma_q_normal_a_new[, l]),
          unrestricted_ppi_min = min(unrestricted_gamma),
          unrestricted_ppi_mean = mean(unrestricted_gamma),
          unrestricted_ppi_max = max(unrestricted_gamma),
          ppi_min = min(mu_q_gamma_a_new[, l]),
          ppi_mean = mean(mu_q_gamma_a_new[, l]),
          ppi_max = max(mu_q_gamma_a_new[, l]),
          log_prior_odds_mean = mean(log_prior_odds),
          gate_active = gate_active
        )
      )
    }
  }

  create_named_list(mu_q_normal_a = mu_q_normal_a_new,
                    Sigma_q_normal_a = Sigma_q_normal_a_new,
                    mu_q_gamma_a = mu_q_gamma_a_new,
                    mu_q_a = mu_q_a_work,
                    term_a = term_a_work)
}


#' Update study-specific loadings b_{sjl} and PPI
#'
#' Within-study only: each study's loadings are independent.
#'
#' @param tr_xi_specific S x L_s x max(n_s) array of E[||g^{(l)}_{si}||^2]
#'
#' @noRd
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

  if (!.has_specific(L_s) || is.null(mu_q_nu_psi)) return(NULL)

  mu_q_normal_b_new <- Sigma_q_normal_b_new <- vector("list", S)
  mu_q_gamma_b_new <- mu_q_b_spec_new <- term_b_spec_new <- vector("list", S)
  mu_q_b_specific_work <- mu_q_b_specific
  term_b_specific_work <- term_b_specific
  trace_driver <- .driver_diagnostic_should_trace()

  for (s in 1:S) {
    L_ss <- .L_s_at(L_s, s)

    mu_q_normal_b_new[[s]] <- matrix(0, p, L_ss)
    Sigma_q_normal_b_new[[s]] <- matrix(1, p, L_ss)
    mu_q_gamma_b_new[[s]] <- matrix(0.5, p, L_ss)
    if (L_ss == 0L) {
      mu_q_b_spec_new[[s]] <- mu_q_b_specific_work[[s]]
      term_b_spec_new[[s]] <- term_b_specific_work[[s]]
      next
    }

    # Precompute within study
    list_cp_C_nu_mu_s <- lapply(1:n_s[s], function(i) {
      sapply(1:p, function(jj) list_cp_C[[s]][[i]] %*% mu_q_nu_mu[[s]][[jj]])
    })

    # Precompute shared factor contributions for residual
    f_ct_contrib_s <- lapply(1:n_s[s], function(i) {
      contrib <- matrix(0, nrow = K_total, ncol = p)
      for (l in seq_len(L_f)) {
        f_ct <- crossprod(C[[s]][[i]], as.vector(
          C[[s]][[i]] %*% mu_q_nu_phi[[l]] %*% mu_q_zeta[[s]][[l]][i, ]))
        contrib <- contrib + tcrossprod(as.vector(f_ct), mu_q_a[, l])
      }
      contrib
    })

    for (l in seq_len(L_ss)) {

      # Aggregate E[||g||^2] within study s
      rs_tr_xi_sl <- 0
      for (i in 1:n_s[s]) {
        rs_tr_xi_sl <- rs_tr_xi_sl + tr_xi_specific[s, l, i]
      }

      tcp_psi_xi_sl <- tcrossprod(mu_q_nu_psi[[s]][[l]], mu_q_xi[[s]][[l]])

      sum_mu_b <- rep(0, p)
      if (trace_driver) data_envelope <- 0
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
          if (L_ss > 1L) {
            for (l_other in setdiff(seq_len(L_ss), l)) {
              ct_g_other <- list_cp_C[[s]][[i]] %*%
                tcrossprod(mu_q_nu_psi[[s]][[l_other]],
                           mu_q_xi[[s]][[l_other]])[, i]
              term <- term - mu_q_b_specific_work[[s]][jj, l_other] *
                as.numeric(crossprod(tcp_psi_xi_li, as.vector(ct_g_other)))
            }
          }

          sum_mu_b[jj] <- sum_mu_b[jj] + eps_prec_j[jj] * term
          if (trace_driver) {
            data_envelope <- data_envelope + abs(eps_prec_j[jj] * term)
          }
        }
      }

      # ---- Slab variance and mean ----
      Sigma_q_normal_b_new[[s]][, l] <-
        1 / (c_val * (eps_prec_j * rs_tr_xi_sl + 1))
      mu_q_normal_b_new[[s]][, l] <-
        c_val * Sigma_q_normal_b_new[[s]][, l] * sum_mu_b

      # ---- PPI from the joint tempered spike-and-slab normaliser ----
      log_prior_odds <- .loading_log_prior_odds(
        mu_q_log_omega_b[[s]], mu_q_log_1_omega_b[[s]], l, p,
        paste0("Study ", s, " specific"))
      logit_gamma <- .tempered_spike_slab_logit(
        log_prior_odds = log_prior_odds,
        slab_mean = mu_q_normal_b_new[[s]][, l],
        slab_variance = Sigma_q_normal_b_new[[s]][, l],
        c_val = c_val)

      unrestricted_gamma <- stats::plogis(logit_gamma)
      gate_active <- .driver_diagnostic_dense_gate_active()
      mu_q_gamma_b_new[[s]][, l] <-
        if (gate_active) 1 else unrestricted_gamma
      mu_q_b_specific_work[[s]][, l] <-
        mu_q_gamma_b_new[[s]][, l] * mu_q_normal_b_new[[s]][, l]
      term_b_specific_work[[s]][, l] <-
        (Sigma_q_normal_b_new[[s]][, l] +
         mu_q_normal_b_new[[s]][, l]^2) * mu_q_gamma_b_new[[s]][, l]
      if (trace_driver) {
        driver_norm <- .driver_diagnostic_norm(sum_mu_b)
        precision_values <- 1 / Sigma_q_normal_b_new[[s]][, l]
        .driver_diagnostic_record(
          block = "specific_loading", study = s, factor = l,
          metrics = list(
            data_driver_norm = driver_norm,
            data_envelope = data_envelope,
            data_coherence =
              driver_norm / max(data_envelope, .Machine$double.eps),
            linear_driver_norm = driver_norm,
            precision_trace = sum(precision_values),
            precision_frobenius =
              .driver_diagnostic_norm(precision_values),
            prior_precision_trace = c_val * p,
            driver_precision_ratio =
              driver_norm /
              max(.driver_diagnostic_norm(precision_values),
                  .Machine$double.eps),
            posterior_mean_norm =
              .driver_diagnostic_norm(mu_q_b_specific_work[[s]][, l]),
            slab_mean_norm =
              .driver_diagnostic_norm(mu_q_normal_b_new[[s]][, l]),
            slab_variance_mean =
              mean(Sigma_q_normal_b_new[[s]][, l]),
            unrestricted_ppi_min = min(unrestricted_gamma),
            unrestricted_ppi_mean = mean(unrestricted_gamma),
            unrestricted_ppi_max = max(unrestricted_gamma),
            ppi_min = min(mu_q_gamma_b_new[[s]][, l]),
            ppi_mean = mean(mu_q_gamma_b_new[[s]][, l]),
            ppi_max = max(mu_q_gamma_b_new[[s]][, l]),
            log_prior_odds_mean = mean(log_prior_odds),
            gate_active = gate_active
          )
        )
      }
    }

    mu_q_b_spec_new[[s]] <- mu_q_b_specific_work[[s]]
    term_b_spec_new[[s]] <- term_b_specific_work[[s]]
  }

  create_named_list(mu_q_normal_b = mu_q_normal_b_new,
                    Sigma_q_normal_b = Sigma_q_normal_b_new,
                    mu_q_gamma_b = mu_q_gamma_b_new,
                    mu_q_b_specific = mu_q_b_spec_new,
                    term_b_specific = term_b_spec_new)
}
