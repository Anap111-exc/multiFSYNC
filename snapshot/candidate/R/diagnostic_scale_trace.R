# =============================================================================
# diagnostic_scale_trace.R -- private initialization-scale path diagnostics
#
# This file deliberately exposes no public API.  The diagnostic wrapper uses a
# private option to pass tracing instructions through bayesSYNC_multi() without
# changing its public formals.  Unless that option is present, the main fit does
# not construct trace state or request extra fitted-value/RSS work.
# =============================================================================

.multiFSYNC_scale_trace_option <-
  "multiFSYNC.__private_scale_trace_8d72c46f__"

.default_scale_trace_control <- function() {
  list(
    include_initial = TRUE,
    include_annealing = TRUE,
    t1_sweeps = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
    block_annealing_sweeps = integer()
  )
}

.validate_scale_trace_control <- function(trace_control) {
  if (is.null(trace_control)) {
    trace_control <- .default_scale_trace_control()
  }
  if (!is.list(trace_control) || is.null(names(trace_control)) ||
      any(!nzchar(names(trace_control))) ||
      anyDuplicated(names(trace_control))) {
    stop("trace_control must be a uniquely named list.")
  }
  allowed <- c(
    "include_initial", "include_annealing", "t1_sweeps",
    "block_annealing_sweeps"
  )
  unknown <- setdiff(names(trace_control), allowed)
  if (length(unknown)) {
    stop(
      "Unknown trace_control field(s): ",
      paste(unknown, collapse = ", "), "."
    )
  }
  defaults <- .default_scale_trace_control()
  trace_control <- utils::modifyList(defaults, trace_control)

  for (field in c("include_initial", "include_annealing")) {
    value <- trace_control[[field]]
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("trace_control$", field, " must be TRUE or FALSE.")
    }
  }

  validate_sweeps <- function(sweeps, field) {
    if (!is.numeric(sweeps) ||
        any(!is.finite(sweeps)) ||
        any(sweeps < 1) ||
        any(abs(sweeps - round(sweeps)) > sqrt(.Machine$double.eps)) ||
        any(sweeps > .Machine$integer.max)) {
      stop("trace_control$", field, " must contain positive integers.")
    }
    sort(unique(as.integer(round(sweeps))))
  }
  trace_control$t1_sweeps <- validate_sweeps(
    trace_control$t1_sweeps, "t1_sweeps"
  )
  trace_control$block_annealing_sweeps <- validate_sweeps(
    trace_control$block_annealing_sweeps,
    "block_annealing_sweeps"
  )
  trace_control
}

#' Run a fit with private initialization-scale tracing
#'
#' This intentionally unexported helper leaves the public bayesSYNC_multi()
#' interface unchanged.  It restores both an existing private option value and
#' an absent option exactly, including when fitting exits with an error.
#'
#' @keywords internal
.bayesSYNC_multi_scale_trace <- function(..., trace_control = NULL) {
  trace_control <- .validate_scale_trace_control(trace_control)
  option_name <- .multiFSYNC_scale_trace_option
  old_options <- options()
  option_existed <- option_name %in% names(old_options)
  old_value <- old_options[[option_name]]

  on.exit({
    restore <- stats::setNames(
      list(if (option_existed) old_value else NULL),
      option_name
    )
    options(restore)
  }, add = TRUE)

  options(stats::setNames(list(trace_control), option_name))
  bayesSYNC_multi(...)
}

.scale_trace_numeric <- function(x) {
  if (is.null(x)) return(numeric())
  as.numeric(unlist(x, recursive = TRUE, use.names = FALSE))
}

.scale_trace_finite <- function(x) {
  x <- .scale_trace_numeric(x)
  x[is.finite(x)]
}

.scale_trace_rms <- function(x, absent = NA_real_) {
  x <- .scale_trace_finite(x)
  if (!length(x)) return(absent)
  sqrt(mean(x^2))
}

.scale_trace_mean <- function(x, absent = NA_real_) {
  x <- .scale_trace_finite(x)
  if (!length(x)) return(absent)
  mean(x)
}

.scale_trace_triplet <- function(x, absent = NA_real_) {
  x <- .scale_trace_finite(x)
  if (!length(x)) {
    return(c(min = absent, median = absent, max = absent))
  }
  c(min = min(x), median = stats::median(x), max = max(x))
}

.scale_trace_score_second <- function(means, covariances) {
  values <- numeric()
  if (!length(means)) return(values)
  for (s in seq_along(means)) {
    if (!length(means[[s]])) next
    for (l in seq_along(means[[s]])) {
      mean_sl <- means[[s]][[l]]
      covariance_sl <- covariances[[s]][[l]]
      if (!length(mean_sl)) next
      for (i in seq_len(nrow(mean_sl))) {
        values <- c(
          values,
          mean_sl[i, ]^2 + diag(covariance_sl[[i]])
        )
      }
    }
  }
  values
}

.scale_trace_vector_second <- function(means, covariances) {
  if (is.null(means) || is.null(covariances)) return(numeric())
  if (is.atomic(means) && is.null(dim(means)) &&
      is.matrix(covariances)) {
    return(as.numeric(means)^2 + diag(covariances))
  }
  values <- numeric()
  for (index in seq_along(means)) {
    values <- c(
      values,
      .scale_trace_vector_second(means[[index]], covariances[[index]])
    )
  }
  values
}

.scale_trace_dense_second <- function(C_g, coefficient_means,
                                      coefficient_covariances) {
  mean_values <- numeric()
  second_values <- numeric()
  if (!length(coefficient_means)) {
    return(list(mean = mean_values, second = second_values))
  }

  append_components <- function(mean_matrix, covariance_list) {
    local_mean <- numeric()
    local_second <- numeric()
    if (is.null(dim(mean_matrix))) {
      mean_matrix <- matrix(mean_matrix, ncol = 1L)
    }
    for (m in seq_len(ncol(mean_matrix))) {
      curve_mean <- drop(C_g %*% mean_matrix[, m])
      projected_covariance <- C_g %*% covariance_list[[m]]
      curve_variance <- rowSums(projected_covariance * C_g)
      local_mean <- c(local_mean, curve_mean)
      local_second <- c(
        local_second,
        curve_mean^2 + pmax(curve_variance, 0)
      )
    }
    list(mean = local_mean, second = local_second)
  }

  # Shared coefficients are factor -> coefficient matrix.
  if (is.matrix(coefficient_means[[1L]])) {
    for (l in seq_along(coefficient_means)) {
      values <- append_components(
        coefficient_means[[l]],
        coefficient_covariances[[l]]
      )
      mean_values <- c(mean_values, values$mean)
      second_values <- c(second_values, values$second)
    }
  } else {
    # Specific coefficients are study -> factor -> coefficient matrix.
    for (s in seq_along(coefficient_means)) {
      for (l in seq_along(coefficient_means[[s]])) {
        values <- append_components(
          coefficient_means[[s]][[l]],
          coefficient_covariances[[s]][[l]]
        )
        mean_values <- c(mean_values, values$mean)
        second_values <- c(second_values, values$second)
      }
    }
  }
  list(mean = mean_values, second = second_values)
}

.scale_trace_mean_contributions <- function(Y, C, Z, state) {
  observed <- numeric()
  mean_part <- numeric()
  beta_part <- numeric()
  shared_part <- numeric()
  specific_part <- numeric()

  for (s in seq_len(state$S)) {
    for (i in seq_len(state$n_s[s])) {
      C_si <- C[[s]][[i]]
      n_i <- nrow(C_si)
      for (j in seq_len(state$p)) {
        observed_one <- as.numeric(Y[[s]][[i]][[j]])
        mean_one <- drop(C_si %*% state$mu_q_nu_mu[[s]][[j]])

        beta_one <- numeric(n_i)
        if (state$d > 0L) {
          for (r in seq_len(state$d)) {
            beta_one <- beta_one +
              Z[[s]][i, r] *
              drop(C_si %*% state$mu_q_nu_beta[[j]][[r]])
          }
        }

        shared_one <- numeric(n_i)
        for (l in seq_len(state$L_f)) {
          factor_curve <- drop(
            (C_si %*% state$mu_q_nu_phi[[l]]) %*%
              state$mu_q_zeta[[s]][[l]][i, ]
          )
          shared_one <- shared_one +
            state$mu_q_a[j, l] * factor_curve
        }

        specific_one <- numeric(n_i)
        for (l in seq_len(.L_s_at(state$L_s, s))) {
          factor_curve <- drop(
            (C_si %*% state$mu_q_nu_psi[[s]][[l]]) %*%
              state$mu_q_xi[[s]][[l]][i, ]
          )
          specific_one <- specific_one +
            state$mu_q_b_specific[[s]][j, l] * factor_curve
        }

        observed <- c(observed, observed_one)
        mean_part <- c(mean_part, mean_one)
        beta_part <- c(beta_part, beta_one)
        shared_part <- c(shared_part, shared_one)
        specific_part <- c(specific_part, specific_one)
      }
    }
  }

  fitted <- mean_part + beta_part + shared_part + specific_part
  list(
    observed = observed,
    fitted = fitted,
    mean = mean_part,
    beta = beta_part,
    shared = shared_part,
    specific = specific_part
  )
}

.scale_trace_elbo_component_names <- c(
  "data_likelihood",
  "osullivan_mu", "osullivan_beta",
  "osullivan_phi", "osullivan_psi",
  "score_zeta", "score_xi",
  "spike_slab_shared", "spike_slab_specific",
  "omega_shared", "omega_specific",
  "half_cauchy_eps", "half_cauchy_mu",
  "half_cauchy_beta", "half_cauchy_phi",
  "half_cauchy_psi"
)

.make_scale_trace_row <- function(stage, iteration, t1_sweep, temperature,
                                  Y, C, C_g, Z, state,
                                  expected_rss, cached_fitted = NULL,
                                  elbo_result = NULL,
                                  elbo_reason = NA_character_,
                                  block = "end_of_sweep") {
  contributions <- .scale_trace_mean_contributions(Y, C, Z, state)
  raw_prediction_fields <- c(
    "fitted", "mean", "beta", "shared", "specific"
  )
  raw_prediction_nonfinite <- vapply(
    contributions[raw_prediction_fields],
    function(values) sum(!is.finite(values)),
    integer(1)
  )
  expected_rss_values <- .scale_trace_finite(expected_rss)
  total_observations <- length(contributions$observed)

  zeta_second <- .scale_trace_score_second(
    state$mu_q_zeta, state$Sigma_q_zeta
  )
  xi_second <- if (.has_specific(state$L_s)) {
    .scale_trace_score_second(state$mu_q_xi, state$Sigma_q_xi)
  } else {
    numeric()
  }
  score_second <- c(zeta_second, xi_second)
  mu_second <- .scale_trace_vector_second(
    state$mu_q_nu_mu, state$Sigma_q_nu_mu
  )
  beta_second <- .scale_trace_vector_second(
    state$mu_q_nu_beta, state$Sigma_q_nu_beta
  )

  slab_shared_second <- if (state$L_f > 0L) {
    state$Sigma_q_normal_a + state$mu_q_normal_a^2
  } else {
    numeric()
  }
  slab_specific_second <- if (.has_specific(state$L_s)) {
    unlist(Map(
      function(mu, variance) variance + mu^2,
      state$mu_q_normal_b, state$Sigma_q_normal_b
    ), use.names = FALSE)
  } else {
    numeric()
  }
  effective_shared_second <- if (state$L_f > 0L) {
    state$term_a
  } else {
    numeric()
  }
  effective_specific_second <- if (.has_specific(state$L_s)) {
    .scale_trace_numeric(state$term_b_specific)
  } else {
    numeric()
  }

  dense_phi <- .scale_trace_dense_second(
    C_g, state$mu_q_nu_phi, state$Sigma_q_nu_phi
  )
  dense_psi <- if (.has_specific(state$L_s)) {
    .scale_trace_dense_second(
      C_g, state$mu_q_nu_psi, state$Sigma_q_nu_psi
    )
  } else {
    list(mean = numeric(), second = numeric())
  }

  ppi_shared_raw <- if (state$L_f > 0L) {
    .scale_trace_numeric(state$mu_q_gamma_a)
  } else {
    numeric()
  }
  ppi_specific_raw <- if (.has_specific(state$L_s)) {
    .scale_trace_numeric(state$mu_q_gamma_b)
  } else {
    numeric()
  }
  ppi_shared <- ppi_shared_raw[is.finite(ppi_shared_raw)]
  ppi_specific <- ppi_specific_raw[is.finite(ppi_specific_raw)]
  ppi_raw <- c(ppi_shared_raw, ppi_specific_raw)
  ppi <- c(ppi_shared, ppi_specific)
  omega_shared <- if (state$L_f > 0L) {
    .scale_trace_finite(
      state$c_1_omega_a / (state$c_1_omega_a + state$d_1_omega_a)
    )
  } else {
    numeric()
  }
  omega_specific <- if (.has_specific(state$L_s)) {
    .scale_trace_finite(Map(
      function(shape1, shape2) shape1 / (shape1 + shape2),
      state$c_1_omega_b, state$d_1_omega_b
    ))
  } else {
    numeric()
  }
  omega <- c(omega_shared, omega_specific)
  omega_summary <- .scale_trace_triplet(omega)
  omega_shared_summary <- .scale_trace_triplet(omega_shared)
  omega_specific_summary <- .scale_trace_triplet(omega_specific)

  smooth_mu <- .scale_trace_finite(state$mu_q_recip_sigsq_mu)
  smooth_beta <- .scale_trace_finite(state$mu_q_recip_sigsq_beta)
  smooth_phi <- .scale_trace_finite(state$mu_q_recip_sigsq_phi)
  smooth_psi <- .scale_trace_finite(state$mu_q_recip_sigsq_psi)
  smooth_all <- c(smooth_mu, smooth_beta, smooth_phi, smooth_psi)
  smooth_summary <- .scale_trace_triplet(smooth_all)

  aux_eps <- .scale_trace_finite(state$mu_q_recip_a_eps)
  aux_mu <- .scale_trace_finite(state$mu_q_recip_a_mu)
  aux_beta <- .scale_trace_finite(state$mu_q_recip_a_beta)
  aux_phi <- .scale_trace_finite(state$mu_q_recip_a_phi)
  aux_psi <- .scale_trace_finite(state$mu_q_recip_a_psi)
  aux_all <- c(aux_eps, aux_mu, aux_beta, aux_phi, aux_psi)
  aux_summary <- .scale_trace_triplet(aux_all)

  noise_summary <- .scale_trace_triplet(state$mu_q_recip_sigsq_eps)
  cached_fitted <- .scale_trace_finite(cached_fitted)
  fitted_cache_difference <- if (length(cached_fitted)) {
    if (length(cached_fitted) != length(contributions$fitted)) {
      Inf
    } else {
      max(abs(cached_fitted - contributions$fitted))
    }
  } else {
    NA_real_
  }

  elbo_components <- stats::setNames(
    rep(NA_real_, length(.scale_trace_elbo_component_names)),
    paste0("elbo_", .scale_trace_elbo_component_names)
  )
  elbo <- NA_real_
  objective_scope <- NA_character_
  if (!is.null(elbo_result)) {
    elbo <- as.numeric(elbo_result$total)
    objective_scope <- if (is.null(elbo_result$scope)) {
      NA_character_
    } else {
      as.character(elbo_result$scope)
    }
    available <- intersect(
      names(elbo_result$components),
      .scale_trace_elbo_component_names
    )
    elbo_components[paste0("elbo_", available)] <-
      as.numeric(elbo_result$components[available])
  }

  row <- c(
    list(
      stage = as.character(stage),
      block = as.character(block),
      iteration = as.integer(iteration),
      t1_sweep = as.integer(t1_sweep),
      temperature = as.numeric(temperature),
      elbo = elbo,
      objective_scope = objective_scope,
      elbo_reason = as.character(elbo_reason),
      observed_rms = .scale_trace_rms(contributions$observed),
      fitted_rms = .scale_trace_rms(contributions$fitted),
      mean_contribution_rms = .scale_trace_rms(contributions$mean, 0),
      beta_contribution_rms = .scale_trace_rms(contributions$beta, 0),
      shared_contribution_rms = .scale_trace_rms(contributions$shared, 0),
      specific_contribution_rms = .scale_trace_rms(
        contributions$specific, 0
      ),
      raw_observation_count = length(contributions$observed),
      raw_observation_nonfinite_count =
        sum(!is.finite(contributions$observed)),
      raw_prediction_value_count = sum(vapply(
        contributions[raw_prediction_fields],
        length,
        integer(1)
      )),
      raw_fitted_nonfinite_count =
        unname(raw_prediction_nonfinite[["fitted"]]),
      raw_mean_contribution_nonfinite_count =
        unname(raw_prediction_nonfinite[["mean"]]),
      raw_beta_contribution_nonfinite_count =
        unname(raw_prediction_nonfinite[["beta"]]),
      raw_shared_contribution_nonfinite_count =
        unname(raw_prediction_nonfinite[["shared"]]),
      raw_specific_contribution_nonfinite_count =
        unname(raw_prediction_nonfinite[["specific"]]),
      raw_factor_contribution_nonfinite_count =
        unname(raw_prediction_nonfinite[["shared"]] +
                 raw_prediction_nonfinite[["specific"]]),
      raw_all_prediction_nonfinite_count =
        sum(raw_prediction_nonfinite),
      expected_rss_total = sum(expected_rss_values),
      expected_rss_per_obs = sum(expected_rss_values) / total_observations,
      mean_residual_sse = sum(
        (contributions$observed - contributions$fitted)^2
      ),
      fitted_cache_max_abs_diff = fitted_cache_difference,
      noise_precision_min = unname(noise_summary["min"]),
      noise_precision_median = unname(noise_summary["median"]),
      noise_precision_max = unname(noise_summary["max"]),
      mu_rms = .scale_trace_rms(state$mu_q_nu_mu),
      beta_rms = .scale_trace_rms(state$mu_q_nu_beta),
      phi_rms = .scale_trace_rms(state$mu_q_nu_phi),
      psi_rms = .scale_trace_rms(state$mu_q_nu_psi),
      zeta_rms = .scale_trace_rms(state$mu_q_zeta),
      xi_rms = .scale_trace_rms(state$mu_q_xi),
      a_rms = .scale_trace_rms(state$mu_q_a),
      b_rms = .scale_trace_rms(state$mu_q_b_specific),
      mu_q_second_moment = .scale_trace_mean(mu_second),
      beta_q_second_moment = .scale_trace_mean(beta_second),
      score_q_second_moment = .scale_trace_mean(score_second),
      zeta_score_q_second_moment = .scale_trace_mean(zeta_second),
      xi_score_q_second_moment = .scale_trace_mean(xi_second),
      slab_loading_q_second_moment = .scale_trace_mean(
        c(slab_shared_second, slab_specific_second)
      ),
      shared_slab_loading_q_second_moment =
        .scale_trace_mean(slab_shared_second),
      specific_slab_loading_q_second_moment =
        .scale_trace_mean(slab_specific_second),
      effective_loading_q_second_moment = .scale_trace_mean(
        c(effective_shared_second, effective_specific_second)
      ),
      shared_effective_loading_q_second_moment =
        .scale_trace_mean(effective_shared_second),
      specific_effective_loading_q_second_moment =
        .scale_trace_mean(effective_specific_second),
      dense_phi_mean_rms = .scale_trace_rms(dense_phi$mean),
      dense_phi_qsecond_rms = sqrt(
        .scale_trace_mean(dense_phi$second)
      ),
      dense_psi_mean_rms = .scale_trace_rms(dense_psi$mean),
      dense_psi_qsecond_rms = sqrt(
        .scale_trace_mean(dense_psi$second)
      ),
      ppi_min = unname(.scale_trace_triplet(ppi)["min"]),
      ppi_mean = .scale_trace_mean(ppi),
      ppi_max = unname(.scale_trace_triplet(ppi)["max"]),
      ppi_gt_half = sum(ppi > 0.5),
      raw_ppi_count = length(ppi_raw),
      raw_ppi_nonfinite_count = sum(!is.finite(ppi_raw)),
      raw_ppi_out_of_range_count = sum(
        is.finite(ppi_raw) & (ppi_raw < 0 | ppi_raw > 1)
      ),
      shared_ppi_mean = .scale_trace_mean(ppi_shared),
      specific_ppi_mean = .scale_trace_mean(ppi_specific),
      omega_min = unname(omega_summary["min"]),
      omega_mean = .scale_trace_mean(omega),
      omega_max = unname(omega_summary["max"]),
      shared_omega_min = unname(omega_shared_summary["min"]),
      shared_omega_mean = .scale_trace_mean(omega_shared),
      shared_omega_max = unname(omega_shared_summary["max"]),
      specific_omega_min = unname(omega_specific_summary["min"]),
      specific_omega_mean = .scale_trace_mean(omega_specific),
      specific_omega_max = unname(omega_specific_summary["max"]),
      smooth_precision_min = unname(smooth_summary["min"]),
      smooth_precision_median = unname(smooth_summary["median"]),
      smooth_precision_max = unname(smooth_summary["max"]),
      smooth_mu_precision_mean = .scale_trace_mean(smooth_mu),
      smooth_beta_precision_mean = .scale_trace_mean(smooth_beta),
      smooth_phi_precision_mean = .scale_trace_mean(smooth_phi),
      smooth_psi_precision_mean = .scale_trace_mean(smooth_psi),
      half_cauchy_aux_min = unname(aux_summary["min"]),
      half_cauchy_aux_median = unname(aux_summary["median"]),
      half_cauchy_aux_max = unname(aux_summary["max"]),
      half_cauchy_aux_eps_mean = .scale_trace_mean(aux_eps),
      half_cauchy_aux_mu_mean = .scale_trace_mean(aux_mu),
      half_cauchy_aux_beta_mean = .scale_trace_mean(aux_beta),
      half_cauchy_aux_phi_mean = .scale_trace_mean(aux_phi),
      half_cauchy_aux_psi_mean = .scale_trace_mean(aux_psi)
    ),
    as.list(elbo_components)
  )
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
}
