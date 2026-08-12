# =============================================================================
# update_eigenfunctions.R — Computational time-function spline coefficients
#
# Covers parameter blocks:
#   4.5  nu_{theta,ml}  — shared computational time-function coefficients
#   4.6  nu_{kappa,slm} — study-specific computational time-function coefficients
#
# Legacy API mapping:
#   nu_phi objects represent nu_theta; nu_psi objects represent nu_kappa.
# These functions are unconstrained during CAVI and must not be called FPCA
# eigenfunctions before orthonormalise_multi() recovers the canonical basis.
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

#' Update shared computational time-function coefficients nu_theta
#'
#' Cross-study aggregation: uses data from ALL studies to update
#' shared working time functions.  Each shared factor l can have a different
#' number of components M_f[l].
#' The function and object names retain `nu_phi` for backward compatibility.
#' @param lambda_orth Optional numerical separation penalty.  The default 0
#'   gives the model-derived CAVI update.  Positive values are experimental
#'   regularisation and are not an identification constraint.
#'
#' @noRd
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
                           c_val = 1, n_cpus = 1,
                           lambda_orth = 0) {

  # ---- Step 2: Dimension assertions ----
  stopifnot(length(mu_q_nu_phi) == L_f)
  if (length(lambda_orth) != 1L || !is.finite(lambda_orth) || lambda_orth < 0) {
    stop("lambda_orth must be one finite non-negative number.")
  }
  mu_q_nu_phi_work <- mu_q_nu_phi
  trace_driver <- .driver_diagnostic_should_trace()
  # per-factor component counts may differ, no longer assert uniform ncol

  mu_q_nu_phi_new <- vector("list", L_f)
  Sigma_q_nu_phi_new <- vector("list", L_f)
  inv_Sigma_q_nu_phi_new <- vector("list", L_f)

  for (l in seq_len(L_f)) {

    M_l <- M_f[l]

    mu_q_nu_phi_new[[l]] <- matrix(0, nrow = K_total, ncol = M_l)
    Sigma_q_nu_phi_new[[l]] <- vector("list", M_l)
    inv_Sigma_q_nu_phi_new[[l]] <- vector("list", M_l)

    component_update_order <- .private_phi_component_order(M_l)
    for (m in component_update_order) {

      trace_decomposition <-
        .phi_driver_decomposition_should_trace(l, m)
      if (trace_decomposition) {
        decomposition_sources <- c(
          "data_minus_mean", "covariate",
          "same_factor_other_component",
          "other_shared_factor", "specific_factor"
        )
        decomposition_details <- c(
          "data_minus_mean",
          if (d > 0L) paste0("covariate_", seq_len(d)) else character(),
          paste0(
            "same_factor_component_",
            setdiff(seq_len(M_l), m)
          ),
          paste0(
            "other_shared_factor_",
            setdiff(seq_len(L_f), l)
          ),
          if (.has_specific(L_s)) {
            paste0("specific_factor_", seq_len(.max_L_s(L_s)))
          } else {
            character()
          }
        )
        n_sources <- length(decomposition_sources)
        n_details <- length(decomposition_details)
        source_vectors <- matrix(
          0, nrow = n_sources, ncol = K_total,
          dimnames = list(decomposition_sources, NULL)
        )
        study_source_vectors <- array(
          0, dim = c(S, n_sources, K_total),
          dimnames = list(
            paste0("study_", seq_len(S)),
            decomposition_sources, NULL
          )
        )
        variable_source_vectors <- array(
          0, dim = c(p, n_sources, K_total),
          dimnames = list(
            paste0("variable_", seq_len(p)),
            decomposition_sources, NULL
          )
        )
        study_variable_source_vectors <- array(
          0, dim = c(S, p, n_sources, K_total),
          dimnames = list(
            paste0("study_", seq_len(S)),
            paste0("variable_", seq_len(p)),
            decomposition_sources, NULL
          )
        )
        subject_source_vectors <- array(
          0, dim = c(sum(n_s), n_sources, K_total),
          dimnames = list(
            paste0("subject_", seq_len(sum(n_s))),
            decomposition_sources, NULL
          )
        )
        detail_vectors <- matrix(
          0, nrow = n_details, ncol = K_total,
          dimnames = list(decomposition_details, NULL)
        )
        study_detail_vectors <- array(
          0, dim = c(S, n_details, K_total),
          dimnames = list(
            paste0("study_", seq_len(S)),
            decomposition_details, NULL
          )
        )
        variable_detail_vectors <- array(
          0, dim = c(p, n_details, K_total),
          dimnames = list(
            paste0("variable_", seq_len(p)),
            decomposition_details, NULL
          )
        )
        subject_detail_vectors <- array(
          0, dim = c(sum(n_s), n_details, K_total),
          dimnames = list(
            paste0("subject_", seq_len(sum(n_s))),
            decomposition_details, NULL
          )
        )
        subject_offsets <- c(0L, cumsum(n_s))[seq_len(S)]
        add_decomposition_delta <- function(source, study, subject,
                                            variable, delta,
                                            detail = source) {
          source_index <- match(source, decomposition_sources)
          detail_index <- match(detail, decomposition_details)
          if (is.na(detail_index)) {
            stop("Unknown phi-driver decomposition detail: ", detail)
          }
          global_subject <- subject_offsets[study] + subject
          source_vectors[source_index, ] <<-
            source_vectors[source_index, ] + delta
          study_source_vectors[study, source_index, ] <<-
            study_source_vectors[study, source_index, ] + delta
          variable_source_vectors[variable, source_index, ] <<-
            variable_source_vectors[variable, source_index, ] + delta
          study_variable_source_vectors[
            study, variable, source_index,
          ] <<- study_variable_source_vectors[
            study, variable, source_index,
          ] + delta
          subject_source_vectors[global_subject, source_index, ] <<-
            subject_source_vectors[global_subject, source_index, ] + delta
          detail_vectors[detail_index, ] <<-
            detail_vectors[detail_index, ] + delta
          study_detail_vectors[study, detail_index, ] <<-
            study_detail_vectors[study, detail_index, ] + delta
          variable_detail_vectors[variable, detail_index, ] <<-
            variable_detail_vectors[variable, detail_index, ] + delta
          subject_detail_vectors[global_subject, detail_index, ] <<-
            subject_detail_vectors[global_subject, detail_index, ] + delta
          invisible(NULL)
        }
        precision_by_study <- array(
          0, dim = c(S, K_total, K_total)
        )
        precision_by_study_variable <- array(
          0, dim = c(S, p, K_total, K_total)
        )
        phi_work_before_coordinate <- lapply(
          mu_q_nu_phi_work, function(value) value
        )
      }

      # --- Prior precision ---
      inv_prior <- blkdiag(inv_Sigma_beta,
                           mu_q_recip_sigsq_phi[[l]][m] * diag(K))

      # --- Per-study weighted sum: Σ_s w_s * Σ_i E[zeta^2] * C_i^T C_i ---
      # w_s = Σ_j sigma^{-2}_{eps,sj} * E[a²_{jl}] (study-specific scalar)
      prec_data <- matrix(0, K_total, K_total)
      for (s in 1:S) {
        w_s <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l])
        if (trace_decomposition) {
          score_weighted_cp <- matrix(0, K_total, K_total)
        }
        for (i in 1:n_s[s]) {
          Ez2 <- Sigma_q_zeta[[s]][[l]][[i]][m, m] +
                 mu_q_zeta[[s]][[l]][i, m]^2
          prec_data <- prec_data + w_s * Ez2 * list_cp_C[[s]][[i]]
          if (trace_decomposition) {
            score_weighted_cp <- score_weighted_cp +
              Ez2 * list_cp_C[[s]][[i]]
          }
        }
        if (trace_decomposition) {
          precision_by_study[s, , ] <- w_s * score_weighted_cp
          for (j in seq_len(p)) {
            precision_by_study_variable[s, j, , ] <-
              mu_q_recip_sigsq_eps[s, j] * term_a[j, l] *
              score_weighted_cp
          }
        }
      }

      # Optional numerical separation penalty.  It is assembled separately
      # from the likelihood precision so shared and specific updates use the
      # same scale.  lambda_orth = 0 (the default) adds nothing.
      orth_prec <- matrix(0, K_total, K_total)
      if (lambda_orth > 0 && L_f > 1) {
        for (l_other in setdiff(seq_len(L_f), l)) {
          for (s in 1:S) {
            w_s_other <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l_other])
            for (i in 1:n_s[s]) {
              f_other <- as.vector(C[[s]][[i]] %*%
                mu_q_nu_phi_work[[l_other]] %*%
                mu_q_zeta[[s]][[l_other]][i, ])
              f_ct <- crossprod(C[[s]][[i]], f_other)
              orth_prec <- orth_prec + w_s_other * tcrossprod(f_ct)
            }
          }
        }
      }

      if (lambda_orth > 0 && .has_specific(L_s) && !is.null(mu_q_nu_psi) &&
          !is.null(mu_q_b_specific)) {
        for (s in 1:S) {
          L_ss <- .L_s_at(L_s, s)
          for (l_spec in seq_len(L_ss)) {
            w_s_spec <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l_spec])
            for (i in 1:n_s[s]) {
              g_spec <- as.vector(C[[s]][[i]] %*% mu_q_nu_psi[[s]][[l_spec]] %*% mu_q_xi[[s]][[l_spec]][i, ])
              g_ct <- crossprod(C[[s]][[i]], g_spec)
              orth_prec <- orth_prec + w_s_spec * tcrossprod(g_ct)
            }
          }
        }
      }

      # --- Posterior precision ---
      prec <- c_val * (prec_data + lambda_orth * orth_prec + inv_prior)
      Sigma_q_nu_phi_new[[l]][[m]] <- .inverse_spd(
        prec, context = sprintf("nu_phi[l=%d,m=%d]", l, m))
      inv_Sigma_q_nu_phi_new[[l]][[m]] <- inv_prior

      # --- Linear term (Appendix C.4): cross-study aggregation ---
      sum_term <- rep(0, K_total)
      if (trace_driver) {
        data_sum_term <- rep(0, K_total)
        data_envelope <- 0
      }

      for (s in 1:S) {
        L_ss <- .L_s_at(L_s, s)
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
            if (trace_decomposition) {
              add_decomposition_delta(
                "data_minus_mean", s, i, j,
                coef * as.vector(ct_y - ct_mu)
              )
            }
            if (trace_driver) {
              data_delta <- coef * as.vector(ct_y - ct_mu)
              data_sum_term <- data_sum_term + data_delta
              data_envelope <- data_envelope +
                .driver_diagnostic_norm(data_delta)
            }
          }

          # ---- Subtract beta contribution ----
          if (d > 0 && !is.null(mu_q_nu_beta) && !is.null(Z)) {
            for (j in 1:p) {
              coef <- mu_q_recip_sigsq_eps[s, j] * mu_q_a[j, l] * zeta_silm
              if (abs(coef) < 1e-15) next
              for (r in 1:d) {
                ct_beta <- cp_C_si %*% mu_q_nu_beta[[j]][[r]]
                sum_term <- sum_term - coef * Z[[s]][i, r] * as.vector(ct_beta)
                if (trace_decomposition) {
                  add_decomposition_delta(
                    "covariate", s, i, j,
                    -coef * Z[[s]][i, r] * as.vector(ct_beta),
                    detail = paste0("covariate_", r)
                  )
                }
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

              ct_phi_mt <- cp_C_si %*%
                mu_q_nu_phi_work[[l]][, m_tilde]

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] * term_a[j, l] * Ez_prod
                sum_term <- sum_term - coef * as.vector(ct_phi_mt)
                if (trace_decomposition) {
                  add_decomposition_delta(
                    "same_factor_other_component", s, i, j,
                    -coef * as.vector(ct_phi_mt),
                    detail = paste0(
                      "same_factor_component_", m_tilde
                    )
                  )
                }
              }
            }
          }

          # ---- Subtract other shared factors (l_tilde != l) ----
          if (L_f > 1) {
            for (l_tilde in setdiff(1:L_f, l)) {
              f_other <- as.vector(
                C_si %*% mu_q_nu_phi_work[[l_tilde]] %*%
                mu_q_zeta[[s]][[l_tilde]][i, ]
              )
              ct_f_other <- crossprod(C_si, f_other)

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        mu_q_a[j, l] * mu_q_a[j, l_tilde] * zeta_silm
                if (abs(coef) < 1e-15) next
                sum_term <- sum_term - coef * as.vector(ct_f_other)
                if (trace_decomposition) {
                  add_decomposition_delta(
                    "other_shared_factor", s, i, j,
                    -coef * as.vector(ct_f_other),
                    detail = paste0(
                      "other_shared_factor_", l_tilde
                    )
                  )
                }
              }
            }
          }

          # ---- Subtract all specific factors ----
          if (L_ss > 0L && !is.null(mu_q_nu_psi) && !is.null(mu_q_b_specific)) {
            for (l_spec in seq_len(L_ss)) {
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
                if (trace_decomposition) {
                  add_decomposition_delta(
                    "specific_factor", s, i, j,
                    -coef * as.vector(ct_g_other),
                    detail = paste0("specific_factor_", l_spec)
                  )
                }
              }
            }
          }
        }
      }

      # --- Posterior mean ---
      mu_q_nu_phi_new[[l]][, m] <- as.vector(
        c_val * Sigma_q_nu_phi_new[[l]][[m]] %*% sum_term
      )
      mu_q_nu_phi_work[[l]][, m] <- mu_q_nu_phi_new[[l]][, m]
      if (trace_decomposition) {
        decomposition_environment <-
          .get_phi_driver_decomposition_environment()
        reconstructed_driver <- colSums(source_vectors)
        precision_from_studies <- apply(
          precision_by_study, c(2L, 3L), sum
        )
        precision_from_study_variables <- apply(
          precision_by_study_variable, c(3L, 4L), sum
        )
        score_mean <- lapply(seq_len(S), function(study) {
          mu_q_zeta[[study]][[l]][, m]
        })
        score_second <- lapply(seq_len(S), function(study) {
          vapply(seq_len(n_s[study]), function(subject) {
            Sigma_q_zeta[[study]][[l]][[subject]][m, m] +
              mu_q_zeta[[study]][[l]][subject, m]^2
          }, numeric(1))
        })
        loading_score_product <- lapply(seq_len(S), function(study) {
          outer(score_mean[[study]], mu_q_a[, l])
        })
        untempered_covariance <-
          c_val * Sigma_q_nu_phi_new[[l]][[m]]
        mean_from_untempered_system <- as.vector(
          untempered_covariance %*% sum_term
        )
        .phi_driver_decomposition_record(list(
          metadata = list(
            iteration = decomposition_environment$context$iteration,
            temperature = decomposition_environment$context$temperature,
            annealing = decomposition_environment$context$annealing,
            factor = l,
            component = m,
            inverse_temperature = c_val,
            S = S, n_s = n_s, p = p, d = d,
            L_f = L_f, L_s = L_s, M_f = M_f,
            K = K, K_total = K_total
          ),
          source_vectors = source_vectors,
          study_source_vectors = study_source_vectors,
          variable_source_vectors = variable_source_vectors,
          study_variable_source_vectors =
            study_variable_source_vectors,
          subject_source_vectors = subject_source_vectors,
          detail_vectors = detail_vectors,
          study_detail_vectors = study_detail_vectors,
          variable_detail_vectors = variable_detail_vectors,
          subject_detail_vectors = subject_detail_vectors,
          linear_driver = sum_term,
          reconstructed_driver = reconstructed_driver,
          linear_driver_additivity_max_abs = max(abs(
            sum_term - reconstructed_driver
          )),
          precision_data = prec_data,
          precision_by_study = precision_by_study,
          precision_by_study_variable =
            precision_by_study_variable,
          precision_from_studies = precision_from_studies,
          precision_from_study_variables =
            precision_from_study_variables,
          precision_study_additivity_max_abs = max(abs(
            prec_data - precision_from_studies
          )),
          precision_study_variable_additivity_max_abs = max(abs(
            prec_data - precision_from_study_variables
          )),
          prior_precision = inv_prior,
          orthogonal_precision = orth_prec,
          untempered_posterior_precision =
            prec_data + lambda_orth * orth_prec + inv_prior,
          tempered_posterior_precision = prec,
          tempered_posterior_covariance =
            Sigma_q_nu_phi_new[[l]][[m]],
          untempered_covariance_from_tempered =
            untempered_covariance,
          posterior_mean = mu_q_nu_phi_new[[l]][, m],
          mean_from_untempered_system =
            mean_from_untempered_system,
          temperature_mean_cancellation_max_abs = max(abs(
            mu_q_nu_phi_new[[l]][, m] -
              mean_from_untempered_system
          )),
          target_loading_mean = mu_q_a[, l],
          target_loading_second_moment = term_a[, l],
          target_score_mean = score_mean,
          target_score_second_moment = score_second,
          target_factor_score_mean = lapply(
            seq_len(S), function(study) {
              mu_q_zeta[[study]][[l]]
            }
          ),
          target_factor_score_covariance = lapply(
            seq_len(S), function(study) {
              Sigma_q_zeta[[study]][[l]]
            }
          ),
          loading_score_product = loading_score_product,
          noise_precision = mu_q_recip_sigsq_eps,
          phi_work_before_coordinate =
            phi_work_before_coordinate,
          phi_work_after_coordinate = lapply(
            mu_q_nu_phi_work, function(value) value
          )
        ))
      }
      if (trace_driver) {
        data_norm <- .driver_diagnostic_norm(data_sum_term)
        precision_norm <- .driver_diagnostic_norm(prec)
        .driver_diagnostic_record(
          block = "shared_time_function", factor = l, component = m,
          metrics = list(
            data_driver_norm = data_norm,
            data_envelope = data_envelope,
            data_coherence = data_norm / max(data_envelope, .Machine$double.eps),
            linear_driver_norm = .driver_diagnostic_norm(sum_term),
            precision_trace = sum(diag(prec)),
            precision_frobenius = precision_norm,
            prior_precision_trace = c_val * sum(diag(inv_prior)),
            driver_precision_ratio =
              .driver_diagnostic_norm(sum_term) /
              max(precision_norm, .Machine$double.eps),
            posterior_mean_norm =
              .driver_diagnostic_norm(mu_q_nu_phi_new[[l]][, m])
          )
        )
      }

    } # end m loop
  } # end l loop

  create_named_list(mu_q_nu_phi = mu_q_nu_phi_work,
                    Sigma_q_nu_phi = Sigma_q_nu_phi_new,
                    inv_Sigma_q_nu_phi = inv_Sigma_q_nu_phi_new)
}


#' Update study-specific computational time-function coefficients nu_kappa
#'
#' Within-study only: each study's working time functions are updated
#' independently using only that study's data.  The number of components
#' can vary by study and factor: M_s[[s]][l].
#' The function and object names retain `nu_psi` for backward compatibility.
#'
#' @noRd
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
                           c_val = 1, n_cpus = 1,
                           lambda_orth = 0) {

  if (is.null(mu_q_nu_psi) || !.has_specific(L_s)) return(NULL)
  if (length(lambda_orth) != 1L || !is.finite(lambda_orth) || lambda_orth < 0) {
    stop("lambda_orth must be one finite non-negative number.")
  }
  mu_q_nu_psi_work <- mu_q_nu_psi
  trace_driver <- .driver_diagnostic_should_trace()

  mu_q_nu_psi_new <- vector("list", S)
  Sigma_q_nu_psi_new <- vector("list", S)
  inv_Sigma_q_nu_psi_new <- vector("list", S)

  for (s in 1:S) {
    L_ss <- .L_s_at(L_s, s)

    mu_q_nu_psi_new[[s]] <- vector("list", L_ss)
    Sigma_q_nu_psi_new[[s]] <- vector("list", L_ss)
    inv_Sigma_q_nu_psi_new[[s]] <- vector("list", L_ss)
    if (L_ss == 0L) next

    for (l in seq_len(L_ss)) {

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

        # Optional numerical separation penalty, on the same precision scale
        # used by update_nu_phi().  It is absent when lambda_orth = 0.
        orth_prec <- matrix(0, K_total, K_total)

        # 1) Specific vs other specific factors in the same study
        if (lambda_orth > 0 && L_ss > 1L) {
          for (l_other in setdiff(seq_len(L_ss), l)) {
            w_s_other <- sum(mu_q_recip_sigsq_eps[s, ] * term_b_specific[[s]][, l_other])
            for (i in 1:n_s[s]) {
              g_other <- as.vector(C[[s]][[i]] %*%
                mu_q_nu_psi_work[[s]][[l_other]] %*%
                mu_q_xi[[s]][[l_other]][i, ])
              g_ct <- crossprod(C[[s]][[i]], g_other)
              orth_prec <- orth_prec + w_s_other * tcrossprod(g_ct)
            }
          }
        }
        
        # 2) Specific vs all shared factors (prevent shared from absorbing specific signal)
        if (lambda_orth > 0 && L_f > 0) {
          for (l_shared in seq_len(L_f)) {
            w_sh <- sum(mu_q_recip_sigsq_eps[s, ] * term_a[, l_shared])
            for (i in 1:n_s[s]) {
              f_sh <- as.vector(C[[s]][[i]] %*% mu_q_nu_phi[[l_shared]] %*% mu_q_zeta[[s]][[l_shared]][i, ])
              f_ct <- crossprod(C[[s]][[i]], f_sh)
              orth_prec <- orth_prec + w_sh * tcrossprod(f_ct)
            }
          }
        }

        # --- Posterior precision ---
        prec <- c_val * (sum_vec_b_sl * list_sum_psi_slm +
                         lambda_orth * orth_prec + inv_prior)
        Sigma_q_nu_psi_new[[s]][[l]][[m]] <- .inverse_spd(
          prec, context = sprintf("nu_psi[s=%d,l=%d,m=%d]", s, l, m))
        inv_Sigma_q_nu_psi_new[[s]][[l]][[m]] <- inv_prior

        # --- Linear term: within-study s only ---
        sum_term <- rep(0, K_total)
        if (trace_driver) {
          data_sum_term <- rep(0, K_total)
          data_envelope <- 0
        }

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
            if (trace_driver) {
              data_delta <- coef * as.vector(ct_y - ct_mu)
              data_sum_term <- data_sum_term + data_delta
              data_envelope <- data_envelope +
                .driver_diagnostic_norm(data_delta)
            }
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

              ct_psi_mt <- cp_C_si %*%
                mu_q_nu_psi_work[[s]][[l]][, m_tilde]

              for (j in 1:p) {
                coef <- mu_q_recip_sigsq_eps[s, j] *
                        term_b_specific[[s]][j, l] * Ex_prod
                sum_term <- sum_term - coef * as.vector(ct_psi_mt)
              }
            }
          }

          # ---- Subtract other specific factors (l_tilde != l) ----
          if (L_ss > 1L) {
            for (l_tilde in setdiff(seq_len(L_ss), l)) {
              g_other <- as.vector(
                C_si %*% mu_q_nu_psi_work[[s]][[l_tilde]] %*%
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
        mu_q_nu_psi_work[[s]][[l]][, m] <-
          mu_q_nu_psi_new[[s]][[l]][, m]
        if (trace_driver) {
          data_norm <- .driver_diagnostic_norm(data_sum_term)
          precision_norm <- .driver_diagnostic_norm(prec)
          .driver_diagnostic_record(
            block = "specific_time_function", study = s, factor = l,
            component = m,
            metrics = list(
              data_driver_norm = data_norm,
              data_envelope = data_envelope,
              data_coherence =
                data_norm / max(data_envelope, .Machine$double.eps),
              linear_driver_norm = .driver_diagnostic_norm(sum_term),
              precision_trace = sum(diag(prec)),
              precision_frobenius = precision_norm,
              prior_precision_trace = c_val * sum(diag(inv_prior)),
              driver_precision_ratio =
                .driver_diagnostic_norm(sum_term) /
                max(precision_norm, .Machine$double.eps),
              posterior_mean_norm =
                .driver_diagnostic_norm(mu_q_nu_psi_new[[s]][[l]][, m])
            )
          )
        }
      }
    }
  }

  create_named_list(mu_q_nu_psi = mu_q_nu_psi_work,
                    Sigma_q_nu_psi = Sigma_q_nu_psi_new,
                    inv_Sigma_q_nu_psi = inv_Sigma_q_nu_psi_new)
}
