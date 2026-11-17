# =============================================================================
# diagnostic_block_rss_probe.R -- private selective block RSS diagnostics
#
# Ordinary fits have no probe option. Only explicitly selected annealing
# sweeps and block names request compute_rss_cache(). The probe draws no random
# numbers and never changes a variational coordinate.
# =============================================================================

.multiFSYNC_block_rss_probe_option <-
  "multiFSYNC.__private_block_rss_probe_9f8426c1__"

.validate_block_rss_probe_control <- function(control) {
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))) {
    stop("block RSS probe control must be a uniquely named list.")
  }
  allowed <- c("trace_sweeps", "trace_blocks")
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop(
      "Unknown block RSS probe control field(s): ",
      paste(unknown, collapse = ", "), "."
    )
  }
  sweeps <- control$trace_sweeps
  if (!is.numeric(sweeps) || any(!is.finite(sweeps)) ||
      any(sweeps < 1) ||
      any(abs(sweeps - round(sweeps)) >
            sqrt(.Machine$double.eps)) ||
      any(sweeps > .Machine$integer.max)) {
    stop("control$trace_sweeps must contain positive integers.")
  }
  valid_blocks <- c(
    "before_updates", "after_mu", "after_beta", "after_phi",
    "after_psi", "after_zeta", "after_xi",
    "after_shared_loadings", "after_specific_loadings",
    "after_variances", "after_omega"
  )
  blocks <- control$trace_blocks
  if (!is.character(blocks) || !length(blocks) ||
      anyNA(blocks) || any(!nzchar(blocks)) ||
      any(!blocks %in% valid_blocks)) {
    stop("control$trace_blocks contains an invalid block name.")
  }
  list(
    trace_sweeps = sort(unique(as.integer(round(sweeps)))),
    trace_blocks = unique(blocks)
  )
}

.new_block_rss_probe_environment <- function(control) {
  environment <- new.env(parent = emptyenv())
  environment$control <- .validate_block_rss_probe_control(control)
  environment$global_rows <- list()
  environment$factor_rows <- list()
  environment
}

.get_block_rss_probe_environment <- function() {
  environment <- getOption(.multiFSYNC_block_rss_probe_option, NULL)
  if (is.null(environment)) return(NULL)
  if (!is.environment(environment) ||
      is.null(environment$control) ||
      !is.list(environment$control)) {
    stop("The private block RSS probe option is malformed.")
  }
  environment
}

.block_rss_probe_should_trace <- function(block, iteration, annealing) {
  environment <- .get_block_rss_probe_environment()
  !is.null(environment) &&
    isTRUE(annealing) &&
    iteration %in% environment$control$trace_sweeps &&
    block %in% environment$control$trace_blocks
}

.block_rss_probe_record <- function(
    block, iteration, temperature,
    Y, fitted_values, expected_rss_sum, total_obs_sj, C_g,
    mu_q_nu_phi, Sigma_q_nu_phi,
    mu_q_zeta, Sigma_q_zeta,
    mu_q_normal_a, Sigma_q_normal_a,
    mu_q_gamma_a, mu_q_a, term_a,
    mu_q_recip_sigsq_eps, L_f) {
  environment <- .get_block_rss_probe_environment()
  if (is.null(environment)) return(invisible(NULL))

  observed <- .scale_trace_numeric(Y)
  fitted <- .scale_trace_numeric(fitted_values)
  if (length(observed) != length(fitted)) {
    stop("Block RSS probe fitted-value length mismatch.")
  }
  dense_all <- .scale_trace_dense_second(
    C_g, mu_q_nu_phi, Sigma_q_nu_phi
  )
  global_row <- data.frame(
    block = as.character(block),
    iteration = as.integer(iteration),
    temperature = as.numeric(temperature),
    expected_rss_total = sum(expected_rss_sum),
    expected_rss_per_obs =
      sum(expected_rss_sum) / sum(total_obs_sj),
    mean_residual_sse = sum((observed - fitted)^2),
    mean_residual_rms = sqrt(mean((observed - fitted)^2)),
    phi_rms = .scale_trace_rms(mu_q_nu_phi),
    dense_phi_mean_rms = .scale_trace_rms(dense_all$mean),
    dense_phi_qsecond_rms = sqrt(
      .scale_trace_mean(dense_all$second)
    ),
    zeta_rms = .scale_trace_rms(mu_q_zeta),
    a_rms = .scale_trace_rms(mu_q_a),
    shared_ppi_mean = .scale_trace_mean(mu_q_gamma_a),
    noise_precision_median =
      stats::median(mu_q_recip_sigsq_eps),
    stringsAsFactors = FALSE
  )
  environment$global_rows[[length(environment$global_rows) + 1L]] <-
    global_row

  slab_second <- Sigma_q_normal_a + mu_q_normal_a^2
  for (factor in seq_len(L_f)) {
    dense_factor <- .scale_trace_dense_second(
      C_g,
      list(mu_q_nu_phi[[factor]]),
      list(Sigma_q_nu_phi[[factor]])
    )
    score_factor <- .annealing_light_score_summary(
      factor, mu_q_zeta, Sigma_q_zeta
    )
    environment$factor_rows[[length(environment$factor_rows) + 1L]] <-
      data.frame(
      block = as.character(block),
      iteration = as.integer(iteration),
      temperature = as.numeric(temperature),
      factor = as.integer(factor),
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
        .scale_trace_mean(slab_second[, factor]),
      effective_loading_mean_rms =
        .scale_trace_rms(mu_q_a[, factor]),
      effective_loading_q_second_moment =
        .scale_trace_mean(term_a[, factor]),
      ppi_mean = mean(mu_q_gamma_a[, factor]),
      scale_product_proxy =
        .scale_trace_rms(dense_factor$mean) *
        .scale_trace_rms(score_factor$mean) *
        .scale_trace_rms(mu_q_a[, factor]),
      stringsAsFactors = FALSE
    )
  }
  invisible(NULL)
}
