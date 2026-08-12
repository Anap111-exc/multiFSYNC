# =============================================================================
# continuation_state.R -- Validated T=1 posterior continuation support
#
# These helpers are deliberately unexported.  They let a reduced candidate
# model inherit a complete working-scale variational state without changing
# any CAVI coordinate update, prior, or ELBO expression.
# =============================================================================

.continuation_parameter_names <- function() {
  c(
    "mu_q_nu_mu", "Sigma_q_nu_mu", "inv_Sigma_q_nu_mu",
    "mu_q_nu_beta", "Sigma_q_nu_beta",
    "mu_q_nu_phi", "Sigma_q_nu_phi", "inv_Sigma_q_nu_phi",
    "mu_q_nu_psi", "Sigma_q_nu_psi", "inv_Sigma_q_nu_psi",
    "mu_q_zeta", "Sigma_q_zeta",
    "mu_q_xi", "Sigma_q_xi",
    "mu_q_normal_a", "Sigma_q_normal_a", "mu_q_gamma_a",
    "mu_q_a", "term_a",
    "mu_q_normal_b", "Sigma_q_normal_b", "mu_q_gamma_b",
    "mu_q_b_specific", "term_b_specific",
    "mu_q_recip_sigsq_eps", "mu_q_recip_a_eps",
    "mu_q_recip_sigsq_mu", "mu_q_recip_a_mu",
    "mu_q_recip_sigsq_beta", "mu_q_recip_a_beta",
    "mu_q_recip_sigsq_phi", "mu_q_recip_a_phi",
    "mu_q_recip_sigsq_psi", "mu_q_recip_a_psi"
  )
}

.continuation_all_finite <- function(x) {
  if (is.null(x) || length(x) == 0L) return(TRUE)
  if (is.list(x)) {
    return(all(vapply(x, .continuation_all_finite, logical(1))))
  }
  is.numeric(x) && all(is.finite(x))
}

.continuation_same_shape <- function(x, template) {
  if (is.null(x) || is.null(template)) return(is.null(x) && is.null(template))
  if (is.list(x) || is.list(template)) {
    if (!is.list(x) || !is.list(template) || length(x) != length(template)) {
      return(FALSE)
    }
    if (!length(x)) return(TRUE)
    return(all(unlist(
      Map(.continuation_same_shape, x, template), use.names = FALSE
    )))
  }
  if (is.matrix(x) || is.matrix(template)) {
    return(is.matrix(x) && is.matrix(template) &&
             identical(dim(x), dim(template)))
  }
  is.numeric(x) && is.numeric(template) && length(x) == length(template)
}

.continuation_equal_integer_list <- function(x, y) {
  is.list(x) && is.list(y) && length(x) == length(y) &&
    all(unlist(
      Map(function(a, b) identical(as.integer(a), as.integer(b)), x, y),
      use.names = FALSE
    ))
}

.validate_continuation_state <- function(
    state, template_parameters,
    S, n_s, p, d, L_f, L_s, M_f, M_s, K,
    bool_var_spec_prob, response_center, response_scale,
    time_g) {
  if (!inherits(state, "multiFSYNC_continuation_state") ||
      !is.list(state) || !identical(state$version, 1L)) {
    stop("continuation_state must be created by the internal reduced-state helper.")
  }
  metadata_checks <- c(
    identical(state$S, as.integer(S)),
    identical(state$n_s, as.integer(n_s)),
    identical(state$p, as.integer(p)),
    identical(state$d, as.integer(d)),
    identical(state$L_f, as.integer(L_f)),
    identical(.normalize_L_s(state$L_s, S), as.integer(L_s)),
    identical(state$M_f, as.integer(M_f)),
    .continuation_equal_integer_list(state$M_s, M_s),
    identical(state$K, as.integer(K)),
    identical(state$bool_var_spec_prob, isTRUE(bool_var_spec_prob))
  )
  if (!all(metadata_checks)) {
    stop("continuation_state metadata does not match the requested model dimensions or omega branch.")
  }
  numeric_match <- function(x, y, tolerance = 1e-12) {
    is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      max(abs(as.numeric(x) - as.numeric(y)), 0) <= tolerance
  }
  if (!numeric_match(state$response_center, response_center) ||
      !numeric_match(state$response_scale, response_scale) ||
      !numeric_match(state$time_g, time_g)) {
    stop("continuation_state scaling or dense time grid differs from the requested fit.")
  }

  required <- .continuation_parameter_names()
  if (!is.list(state$parameters) ||
      !identical(names(state$parameters), required) ||
      !identical(names(template_parameters), required)) {
    stop("continuation_state has an incomplete or unexpected parameter set.")
  }
  shape_ok <- Map(
    .continuation_same_shape, state$parameters, template_parameters
  )
  if (!all(unlist(shape_ok, use.names = FALSE))) {
    failed <- required[!unlist(shape_ok, use.names = FALSE)]
    stop("continuation_state parameter shapes do not match: ",
         paste(failed, collapse = ", "), ".")
  }
  finite_ok <- vapply(
    state$parameters, .continuation_all_finite, logical(1)
  )
  if (!all(finite_ok)) {
    stop("continuation_state contains non-finite values in: ",
         paste(required[!finite_ok], collapse = ", "), ".")
  }

  pars <- state$parameters
  consistency_tolerance <- 1e-10
  close_enough <- function(x, y) {
    length(x) == length(y) &&
      max(abs(as.numeric(x) - as.numeric(y)), 0) <=
        consistency_tolerance * (1 + max(abs(as.numeric(y)), 0))
  }
  expected_a <- pars$mu_q_gamma_a * pars$mu_q_normal_a
  expected_term_a <- pars$mu_q_gamma_a *
    (pars$Sigma_q_normal_a + pars$mu_q_normal_a^2)
  if (!close_enough(pars$mu_q_a, expected_a) ||
      !close_enough(pars$term_a, expected_term_a)) {
    stop("continuation_state shared loading moments are internally inconsistent.")
  }
  if (.has_specific(L_s)) {
    for (study in seq_len(S)) {
      expected_b <- pars$mu_q_gamma_b[[study]] *
        pars$mu_q_normal_b[[study]]
      expected_term_b <- pars$mu_q_gamma_b[[study]] *
        (pars$Sigma_q_normal_b[[study]] +
           pars$mu_q_normal_b[[study]]^2)
      if (!close_enough(pars$mu_q_b_specific[[study]], expected_b) ||
          !close_enough(pars$term_b_specific[[study]], expected_term_b)) {
        stop("continuation_state specific loading moments are internally inconsistent in study ",
             study, ".")
      }
    }
  }
  invisible(TRUE)
}

# Map factor indices in the reported, norm-sorted representation back to the
# original working CAVI columns and retain every study-specific upper-bound
# candidate.  The returned object contains no data or simulation truth.
.make_reduced_continuation_state <- function(
    fit, shared_hat_indices = seq_len(fit$L_f)) {
  required_fit <- c(
    "S", "n_s", "p", "d", "L_f", "L_s", "M_f", "M_s", "K",
    "bool_var_spec_prob", "response_center", "response_scale", "time_g",
    "factor_order_shared", .continuation_parameter_names()
  )
  missing_fit <- setdiff(required_fit, names(fit))
  if (length(missing_fit)) {
    stop("fit lacks full variational fields required for continuation: ",
         paste(missing_fit, collapse = ", "), ".")
  }
  shared_hat_indices <- as.integer(shared_hat_indices)
  if (anyNA(shared_hat_indices) || anyDuplicated(shared_hat_indices) ||
      any(shared_hat_indices < 1L | shared_hat_indices > fit$L_f)) {
    if (length(shared_hat_indices)) {
      stop("shared_hat_indices must be unique reported-factor indices.")
    }
  }
  working_shared <- if (length(shared_hat_indices)) {
    as.integer(fit$factor_order_shared[shared_hat_indices])
  } else {
    integer()
  }
  if (anyNA(working_shared) || anyDuplicated(working_shared)) {
    stop("fit$factor_order_shared is not a valid permutation.")
  }

  common_names <- c(
    "mu_q_nu_mu", "Sigma_q_nu_mu", "inv_Sigma_q_nu_mu",
    "mu_q_nu_beta", "Sigma_q_nu_beta",
    "mu_q_nu_psi", "Sigma_q_nu_psi", "inv_Sigma_q_nu_psi",
    "mu_q_xi", "Sigma_q_xi",
    "mu_q_normal_b", "Sigma_q_normal_b", "mu_q_gamma_b",
    "mu_q_b_specific", "term_b_specific",
    "mu_q_recip_sigsq_eps", "mu_q_recip_a_eps",
    "mu_q_recip_sigsq_mu", "mu_q_recip_a_mu",
    "mu_q_recip_sigsq_beta", "mu_q_recip_a_beta",
    "mu_q_recip_sigsq_psi", "mu_q_recip_a_psi"
  )
  parameters <- fit[common_names]
  parameters$mu_q_nu_phi <- fit$mu_q_nu_phi[working_shared]
  parameters$Sigma_q_nu_phi <- fit$Sigma_q_nu_phi[working_shared]
  parameters$inv_Sigma_q_nu_phi <-
    fit$inv_Sigma_q_nu_phi[working_shared]
  parameters$mu_q_zeta <- lapply(fit$mu_q_zeta, `[`, working_shared)
  parameters$Sigma_q_zeta <- lapply(fit$Sigma_q_zeta, `[`, working_shared)
  parameters$mu_q_normal_a <-
    fit$mu_q_normal_a[, working_shared, drop = FALSE]
  parameters$Sigma_q_normal_a <-
    fit$Sigma_q_normal_a[, working_shared, drop = FALSE]
  parameters$mu_q_gamma_a <-
    fit$mu_q_gamma_a[, working_shared, drop = FALSE]
  parameters$mu_q_a <- fit$mu_q_a[, working_shared, drop = FALSE]
  parameters$term_a <- fit$term_a[, working_shared, drop = FALSE]
  parameters$mu_q_recip_sigsq_phi <-
    fit$mu_q_recip_sigsq_phi[working_shared]
  parameters$mu_q_recip_a_phi <- fit$mu_q_recip_a_phi[working_shared]
  parameters <- parameters[.continuation_parameter_names()]

  source_elbo <- if (length(fit$ELBO)) tail(fit$ELBO, 1L) else NA_real_
  L_s_by_study <- .normalize_L_s(fit$L_s, fit$S)
  diagnostics <- list(
    used = TRUE,
    source_fit_seed = fit$fit_seed,
    source_iterations = fit$i_iter,
    source_final_elbo = source_elbo,
    shared_reported_indices = shared_hat_indices,
    shared_working_indices = working_shared,
    inherited_specific_upper_bound = L_s_by_study,
    temperature = 1
  )
  structure(
    list(
      version = 1L,
      S = as.integer(fit$S), n_s = as.integer(fit$n_s),
      p = as.integer(fit$p), d = as.integer(fit$d),
      L_f = as.integer(length(working_shared)),
      L_s = .compact_L_s(L_s_by_study),
      M_f = as.integer(fit$M_f[working_shared]),
      M_s = lapply(fit$M_s, as.integer),
      K = as.integer(fit$K),
      bool_var_spec_prob = isTRUE(fit$bool_var_spec_prob),
      response_center = as.numeric(fit$response_center),
      response_scale = as.numeric(fit$response_scale),
      time_g = as.numeric(fit$time_g),
      parameters = parameters,
      diagnostics = diagnostics
    ),
    class = "multiFSYNC_continuation_state"
  )
}
