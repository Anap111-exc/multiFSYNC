# =============================================================================
# diagnostic_annealing_light_trace.R -- private cached sweep diagnostics
#
# The option is absent in ordinary fits.  The hook is called after the omega
# update and records only already available RSS summaries plus inexpensive
# posterior moment summaries.  It never recomputes fitted trajectories or an
# ELBO/tempered objective and draws no random numbers.
# =============================================================================

.multiFSYNC_annealing_light_trace_option <-
  "multiFSYNC.__private_annealing_light_trace_6e58d90b__"

.validate_annealing_light_trace_control <- function(control) {
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))) {
    stop("light trace control must be a uniquely named list.")
  }
  allowed <- "trace_iterations"
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop(
      "Unknown light trace control field(s): ",
      paste(unknown, collapse = ", "), "."
    )
  }
  iterations <- control$trace_iterations
  if (!is.numeric(iterations) || any(!is.finite(iterations)) ||
      any(iterations < 1) ||
      any(abs(iterations - round(iterations)) >
            sqrt(.Machine$double.eps)) ||
      any(iterations > .Machine$integer.max)) {
    stop("control$trace_iterations must contain positive integers.")
  }
  list(trace_iterations = sort(unique(as.integer(round(iterations)))))
}

.new_annealing_light_trace_environment <- function(control) {
  environment <- new.env(parent = emptyenv())
  environment$control <- .validate_annealing_light_trace_control(control)
  environment$global_rows <- list()
  environment$factor_rows <- list()
  environment$study_rows <- list()
  environment
}

.get_annealing_light_trace_environment <- function() {
  environment <- getOption(
    .multiFSYNC_annealing_light_trace_option, NULL
  )
  if (is.null(environment)) return(NULL)
  if (!is.environment(environment) ||
      is.null(environment$control) ||
      !is.list(environment$control)) {
    stop("The private light trace option is malformed.")
  }
  environment
}

.annealing_light_score_summary <- function(
    factor, mu_q_zeta, Sigma_q_zeta) {
  mean_values <- numeric()
  second_values <- numeric()
  for (study in seq_along(mu_q_zeta)) {
    means <- mu_q_zeta[[study]][[factor]]
    covariances <- Sigma_q_zeta[[study]][[factor]]
    mean_values <- c(mean_values, as.numeric(means))
    for (subject in seq_len(nrow(means))) {
      second_values <- c(
        second_values,
        means[subject, ]^2 + diag(covariances[[subject]])
      )
    }
  }
  list(mean = mean_values, second = second_values)
}

.annealing_light_trace_record <- function(
    iteration, temperature, annealing,
    C_g, expected_rss_sum, total_obs_sj,
    mu_q_nu_mu, Sigma_q_nu_mu,
    mu_q_nu_phi, Sigma_q_nu_phi,
    mu_q_zeta, Sigma_q_zeta,
    mu_q_normal_a, Sigma_q_normal_a,
    mu_q_gamma_a, mu_q_a, term_a,
    mu_q_recip_sigsq_eps,
    mu_q_recip_sigsq_mu, mu_q_recip_sigsq_phi,
    mu_q_recip_a_eps, mu_q_recip_a_mu, mu_q_recip_a_phi,
    c_1_omega_a, d_1_omega_a, L_f) {
  environment <- .get_annealing_light_trace_environment()
  if (is.null(environment)) return(invisible(NULL))
  if (!iteration %in% environment$control$trace_iterations) {
    return(invisible(NULL))
  }

  total_observations <- sum(total_obs_sj)
  dense_all <- .scale_trace_dense_second(
    C_g, mu_q_nu_phi, Sigma_q_nu_phi
  )
  zeta_second_all <- .scale_trace_score_second(
    mu_q_zeta, Sigma_q_zeta
  )
  slab_second_all <- Sigma_q_normal_a + mu_q_normal_a^2
  omega <- c_1_omega_a / (c_1_omega_a + d_1_omega_a)
  noise_summary <- .scale_trace_triplet(mu_q_recip_sigsq_eps)

  global_row <- data.frame(
    iteration = as.integer(iteration),
    temperature = as.numeric(temperature),
    annealing = isTRUE(annealing),
    expected_rss_total = sum(expected_rss_sum),
    expected_rss_per_obs = sum(expected_rss_sum) / total_observations,
    mu_rms = .scale_trace_rms(mu_q_nu_mu),
    mu_q_second_moment = .scale_trace_mean(
      .scale_trace_vector_second(mu_q_nu_mu, Sigma_q_nu_mu)
    ),
    phi_rms = .scale_trace_rms(mu_q_nu_phi),
    dense_phi_mean_rms = .scale_trace_rms(dense_all$mean),
    dense_phi_qsecond_rms = sqrt(
      .scale_trace_mean(dense_all$second)
    ),
    zeta_rms = .scale_trace_rms(mu_q_zeta),
    zeta_score_q_second_moment =
      .scale_trace_mean(zeta_second_all),
    slab_loading_rms = .scale_trace_rms(mu_q_normal_a),
    shared_slab_loading_q_second_moment =
      .scale_trace_mean(slab_second_all),
    a_rms = .scale_trace_rms(mu_q_a),
    shared_effective_loading_q_second_moment =
      .scale_trace_mean(term_a),
    shared_ppi_mean = .scale_trace_mean(mu_q_gamma_a),
    shared_ppi_min = min(mu_q_gamma_a),
    shared_ppi_max = max(mu_q_gamma_a),
    shared_omega_mean = .scale_trace_mean(omega),
    shared_omega_min = min(omega),
    shared_omega_max = max(omega),
    noise_precision_min = unname(noise_summary["min"]),
    noise_precision_median = unname(noise_summary["median"]),
    noise_precision_max = unname(noise_summary["max"]),
    smooth_mu_precision_mean =
      .scale_trace_mean(mu_q_recip_sigsq_mu),
    smooth_phi_precision_mean =
      .scale_trace_mean(mu_q_recip_sigsq_phi),
    half_cauchy_aux_eps_mean =
      .scale_trace_mean(mu_q_recip_a_eps),
    half_cauchy_aux_mu_mean =
      .scale_trace_mean(mu_q_recip_a_mu),
    half_cauchy_aux_phi_mean =
      .scale_trace_mean(mu_q_recip_a_phi),
    stringsAsFactors = FALSE
  )
  environment$global_rows[[length(environment$global_rows) + 1L]] <-
    global_row

  for (study in seq_len(nrow(expected_rss_sum))) {
    environment$study_rows[[length(environment$study_rows) + 1L]] <-
      data.frame(
      iteration = as.integer(iteration),
      temperature = as.numeric(temperature),
      annealing = isTRUE(annealing),
      study = as.integer(study),
      expected_rss_total = sum(expected_rss_sum[study, ]),
      observation_count = sum(total_obs_sj[study, ]),
      expected_rss_per_obs =
        sum(expected_rss_sum[study, ]) / sum(total_obs_sj[study, ]),
      noise_precision_mean =
        mean(mu_q_recip_sigsq_eps[study, ]),
      noise_precision_median =
        stats::median(mu_q_recip_sigsq_eps[study, ]),
      stringsAsFactors = FALSE
    )
  }

  for (factor in seq_len(L_f)) {
    dense_factor <- .scale_trace_dense_second(
      C_g,
      list(mu_q_nu_phi[[factor]]),
      list(Sigma_q_nu_phi[[factor]])
    )
    score_factor <- .annealing_light_score_summary(
      factor, mu_q_zeta, Sigma_q_zeta
    )
    factor_omega <- omega[factor]
    environment$factor_rows[[length(environment$factor_rows) + 1L]] <-
      data.frame(
      iteration = as.integer(iteration),
      temperature = as.numeric(temperature),
      annealing = isTRUE(annealing),
      factor = as.integer(factor),
      phi_coefficient_rms =
        .scale_trace_rms(mu_q_nu_phi[[factor]]),
      dense_phi_mean_rms =
        .scale_trace_rms(dense_factor$mean),
      dense_phi_qsecond_rms = sqrt(
        .scale_trace_mean(dense_factor$second)
      ),
      score_mean_rms = .scale_trace_rms(score_factor$mean),
      score_q_second_moment =
        .scale_trace_mean(score_factor$second),
      slab_loading_mean_rms =
        .scale_trace_rms(mu_q_normal_a[, factor]),
      slab_loading_q_second_moment =
        .scale_trace_mean(slab_second_all[, factor]),
      effective_loading_mean_rms =
        .scale_trace_rms(mu_q_a[, factor]),
      effective_loading_q_second_moment =
        .scale_trace_mean(term_a[, factor]),
      ppi_mean = mean(mu_q_gamma_a[, factor]),
      ppi_min = min(mu_q_gamma_a[, factor]),
      ppi_max = max(mu_q_gamma_a[, factor]),
      omega_mean = factor_omega,
      scale_product_proxy =
        .scale_trace_rms(dense_factor$mean) *
        .scale_trace_rms(score_factor$mean) *
        .scale_trace_rms(mu_q_a[, factor]),
      stringsAsFactors = FALSE
    )
  }
  invisible(NULL)
}
