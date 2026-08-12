# =============================================================================
# diagnostic_driver_trace.R -- private first-sweep driver and gate diagnostics
#
# This file deliberately adds no public argument to bayesSYNC_multi().  A
# private option carries a mutable environment through the fit so internal
# coordinate updates can report the linear driver, precision scale, posterior
# mean, and loading-gate evidence.  The same environment can temporarily hold
# loading inclusion coordinates at gamma = 1 for a bounded number of sweeps.
# Once released, every update is the original spike-and-slab CAVI update for
# the unchanged T = 1 model.
# =============================================================================

.multiFSYNC_driver_diagnostic_option <-
  "multiFSYNC.__private_driver_diagnostic_4a61ef20__"

.default_driver_diagnostic_control <- function() {
  list(
    trace_sweeps = 1:5,
    dense_gate_sweeps = 0L,
    random_scale_calibration = "none",
    pre_score_sweeps = 0L
  )
}

.validate_driver_diagnostic_control <- function(control) {
  if (is.null(control)) control <- .default_driver_diagnostic_control()
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))) {
    stop("driver diagnostic control must be a uniquely named list.")
  }
  allowed <- c(
    "trace_sweeps", "dense_gate_sweeps",
    "random_scale_calibration", "pre_score_sweeps"
  )
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop(
      "Unknown driver diagnostic control field(s): ",
      paste(unknown, collapse = ", "), "."
    )
  }
  control <- utils::modifyList(.default_driver_diagnostic_control(), control)

  sweeps <- control$trace_sweeps
  if (!is.numeric(sweeps) || any(!is.finite(sweeps)) ||
      any(sweeps < 1) ||
      any(abs(sweeps - round(sweeps)) > sqrt(.Machine$double.eps)) ||
      any(sweeps > .Machine$integer.max)) {
    stop("control$trace_sweeps must contain positive integers.")
  }
  control$trace_sweeps <- sort(unique(as.integer(round(sweeps))))

  dense <- control$dense_gate_sweeps
  if (length(dense) != 1L || !is.numeric(dense) || !is.finite(dense) ||
      dense < 0 ||
      abs(dense - round(dense)) > sqrt(.Machine$double.eps) ||
      dense > .Machine$integer.max) {
    stop("control$dense_gate_sweeps must be one non-negative integer.")
  }
  control$dense_gate_sweeps <- as.integer(round(dense))

  calibration <- control$random_scale_calibration
  if (length(calibration) != 1L || !is.character(calibration) ||
      is.na(calibration) ||
      !calibration %in% c("none", "function", "all")) {
    stop(
      "control$random_scale_calibration must be one of ",
      "'none', 'function', or 'all'."
    )
  }

  pre_score <- control$pre_score_sweeps
  if (length(pre_score) != 1L || !is.numeric(pre_score) ||
      !is.finite(pre_score) || pre_score < 0 ||
      abs(pre_score - round(pre_score)) > sqrt(.Machine$double.eps) ||
      pre_score > .Machine$integer.max) {
    stop("control$pre_score_sweeps must be one non-negative integer.")
  }
  control$pre_score_sweeps <- as.integer(round(pre_score))
  control
}

.new_driver_diagnostic_environment <- function(control) {
  environment <- new.env(parent = emptyenv())
  environment$control <- .validate_driver_diagnostic_control(control)
  environment$rows <- list()
  environment$sequence <- 0L
  environment$iteration <- NA_integer_
  environment$temperature <- NA_real_
  environment$annealing <- NA
  environment$phase <- "main"
  environment$calibration_diagnostics <- data.frame()
  environment
}

.get_driver_diagnostic_environment <- function() {
  environment <- getOption(.multiFSYNC_driver_diagnostic_option, NULL)
  if (is.null(environment)) return(NULL)
  if (!is.environment(environment) ||
      is.null(environment$control) ||
      !is.list(environment$control)) {
    stop("The private driver diagnostic option is malformed.")
  }
  environment
}

.driver_diagnostic_set_context <- function(iteration, temperature, annealing) {
  environment <- .get_driver_diagnostic_environment()
  if (is.null(environment)) return(invisible(NULL))
  environment$iteration <- as.integer(iteration)
  environment$temperature <- as.numeric(temperature)
  environment$annealing <- isTRUE(annealing)
  environment$phase <- "main"
  invisible(NULL)
}

.driver_diagnostic_set_phase <- function(phase) {
  environment <- .get_driver_diagnostic_environment()
  if (is.null(environment)) return(invisible(NULL))
  if (length(phase) != 1L || !is.character(phase) || is.na(phase) ||
      !phase %in% c("main", "pre_score")) {
    stop("Private driver diagnostic phase must be 'main' or 'pre_score'.")
  }
  environment$phase <- phase
  invisible(NULL)
}

.driver_diagnostic_should_trace <- function() {
  environment <- .get_driver_diagnostic_environment()
  !is.null(environment) &&
    is.finite(environment$iteration) &&
    environment$iteration %in% environment$control$trace_sweeps
}

.driver_diagnostic_dense_gate_sweeps <- function() {
  environment <- .get_driver_diagnostic_environment()
  if (is.null(environment)) return(0L)
  environment$control$dense_gate_sweeps
}

.driver_diagnostic_initial_dense_gate <- function() {
  .driver_diagnostic_dense_gate_sweeps() > 0L
}

.driver_diagnostic_dense_gate_active <- function() {
  environment <- .get_driver_diagnostic_environment()
  !is.null(environment) &&
    is.finite(environment$iteration) &&
    environment$iteration <= environment$control$dense_gate_sweeps
}

.driver_diagnostic_random_scale_calibration <- function() {
  environment <- .get_driver_diagnostic_environment()
  if (is.null(environment)) return("none")
  environment$control$random_scale_calibration
}

.driver_diagnostic_pre_score_sweeps <- function() {
  environment <- .get_driver_diagnostic_environment()
  if (is.null(environment)) return(0L)
  environment$control$pre_score_sweeps
}

.driver_diagnostic_pre_score_active <- function() {
  environment <- .get_driver_diagnostic_environment()
  !is.null(environment) &&
    is.finite(environment$iteration) &&
    environment$iteration <= environment$control$pre_score_sweeps
}

.driver_diagnostic_norm <- function(x) {
  sqrt(sum(as.numeric(x)^2))
}

.driver_diagnostic_record <- function(
    block, study = NA_integer_, factor = NA_integer_,
    component = NA_integer_, metrics = list()) {
  if (!.driver_diagnostic_should_trace()) return(invisible(NULL))
  environment <- .get_driver_diagnostic_environment()

  defaults <- list(
    data_driver_norm = NA_real_,
    data_envelope = NA_real_,
    data_coherence = NA_real_,
    linear_driver_norm = NA_real_,
    precision_trace = NA_real_,
    precision_frobenius = NA_real_,
    prior_precision_trace = NA_real_,
    driver_precision_ratio = NA_real_,
    posterior_mean_norm = NA_real_,
    slab_mean_norm = NA_real_,
    slab_variance_mean = NA_real_,
    unrestricted_ppi_min = NA_real_,
    unrestricted_ppi_mean = NA_real_,
    unrestricted_ppi_max = NA_real_,
    ppi_min = NA_real_,
    ppi_mean = NA_real_,
    ppi_max = NA_real_,
    log_prior_odds_mean = NA_real_,
    gate_active = .driver_diagnostic_dense_gate_active()
  )
  unknown <- setdiff(names(metrics), names(defaults))
  if (length(unknown)) {
    stop("Unknown driver diagnostic metric(s): ",
         paste(unknown, collapse = ", "), ".")
  }
  values <- utils::modifyList(defaults, metrics)
  environment$sequence <- environment$sequence + 1L
  environment$rows[[length(environment$rows) + 1L]] <- data.frame(
    sequence = environment$sequence,
    iteration = environment$iteration,
    temperature = environment$temperature,
    annealing = environment$annealing,
    phase = environment$phase,
    block = as.character(block),
    study = as.integer(study),
    factor = as.integer(factor),
    component = as.integer(component),
    data_driver_norm = as.numeric(values$data_driver_norm),
    data_envelope = as.numeric(values$data_envelope),
    data_coherence = as.numeric(values$data_coherence),
    linear_driver_norm = as.numeric(values$linear_driver_norm),
    precision_trace = as.numeric(values$precision_trace),
    precision_frobenius = as.numeric(values$precision_frobenius),
    prior_precision_trace = as.numeric(values$prior_precision_trace),
    driver_precision_ratio = as.numeric(values$driver_precision_ratio),
    posterior_mean_norm = as.numeric(values$posterior_mean_norm),
    slab_mean_norm = as.numeric(values$slab_mean_norm),
    slab_variance_mean = as.numeric(values$slab_variance_mean),
    unrestricted_ppi_min = as.numeric(values$unrestricted_ppi_min),
    unrestricted_ppi_mean = as.numeric(values$unrestricted_ppi_mean),
    unrestricted_ppi_max = as.numeric(values$unrestricted_ppi_max),
    ppi_min = as.numeric(values$ppi_min),
    ppi_mean = as.numeric(values$ppi_mean),
    ppi_max = as.numeric(values$ppi_max),
    log_prior_odds_mean = as.numeric(values$log_prior_odds_mean),
    gate_active = isTRUE(values$gate_active),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

.driver_calibration_scale <- function(values) {
  scale <- sqrt(mean(as.numeric(values)^2))
  if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) NA_real_ else scale
}

.calibrate_random_factor_state <- function(
    C_g, time_g,
    mu_q_nu_phi, mu_q_zeta, mu_q_normal_a, Sigma_q_normal_a, mu_q_gamma_a,
    mu_q_nu_psi, mu_q_xi, mu_q_normal_b, Sigma_q_normal_b, mu_q_gamma_b,
    S, L_f, L_s, M_f, M_s, mode) {
  if (!mode %in% c("function", "all")) {
    stop("Private random-state calibration mode must be 'function' or 'all'.")
  }
  weights <- .trap_weights(time_g)
  diagnostics <- list()
  diagnostic_index <- 0L
  add_diagnostic <- function(block, study, factor, component,
                             before, after, scale, success) {
    diagnostic_index <<- diagnostic_index + 1L
    diagnostics[[diagnostic_index]] <<- data.frame(
      block = block,
      study = as.integer(study),
      factor = as.integer(factor),
      component = as.integer(component),
      before = as.numeric(before),
      after = as.numeric(after),
      scale = as.numeric(scale),
      success = isTRUE(success),
      stringsAsFactors = FALSE
    )
  }

  for (l in seq_len(L_f)) {
    for (m in seq_len(M_f[l])) {
      dense <- as.vector(C_g %*% mu_q_nu_phi[[l]][, m])
      scale <- sqrt(sum(weights * dense^2))
      success <- is.finite(scale) &&
        scale > sqrt(.Machine$double.eps)
      if (success) mu_q_nu_phi[[l]][, m] <- mu_q_nu_phi[[l]][, m] / scale
      dense_after <- as.vector(C_g %*% mu_q_nu_phi[[l]][, m])
      add_diagnostic(
        "shared_function_l2", NA_integer_, l, m,
        sqrt(sum(weights * dense^2)),
        sqrt(sum(weights * dense_after^2)), scale, success
      )
    }
    if (identical(mode, "all")) {
      for (m in seq_len(M_f[l])) {
        values <- unlist(lapply(
          seq_len(S), function(s) mu_q_zeta[[s]][[l]][, m]
        ), use.names = FALSE)
        scale <- .driver_calibration_scale(values)
        success <- is.finite(scale)
        if (success) {
          for (s in seq_len(S)) {
            mu_q_zeta[[s]][[l]][, m] <-
              mu_q_zeta[[s]][[l]][, m] / scale
          }
        }
        values_after <- unlist(lapply(
          seq_len(S), function(s) mu_q_zeta[[s]][[l]][, m]
        ), use.names = FALSE)
        add_diagnostic(
          "shared_score_rms", NA_integer_, l, m,
          .driver_calibration_scale(values),
          .driver_calibration_scale(values_after), scale, success
        )
      }
      scale <- .driver_calibration_scale(mu_q_normal_a[, l])
      success <- is.finite(scale)
      before <- .driver_calibration_scale(mu_q_normal_a[, l])
      if (success) mu_q_normal_a[, l] <- mu_q_normal_a[, l] / scale
      add_diagnostic(
        "shared_slab_rms", NA_integer_, l, NA_integer_,
        before, .driver_calibration_scale(mu_q_normal_a[, l]),
        scale, success
      )
    }
  }

  if (.has_specific(L_s)) {
    for (s in seq_len(S)) {
      for (l in seq_len(.L_s_at(L_s, s))) {
        for (m in seq_len(M_s[[s]][l])) {
          dense <- as.vector(C_g %*% mu_q_nu_psi[[s]][[l]][, m])
          scale <- sqrt(sum(weights * dense^2))
          success <- is.finite(scale) &&
            scale > sqrt(.Machine$double.eps)
          if (success) {
            mu_q_nu_psi[[s]][[l]][, m] <-
              mu_q_nu_psi[[s]][[l]][, m] / scale
          }
          dense_after <- as.vector(
            C_g %*% mu_q_nu_psi[[s]][[l]][, m]
          )
          add_diagnostic(
            "specific_function_l2", s, l, m,
            sqrt(sum(weights * dense^2)),
            sqrt(sum(weights * dense_after^2)), scale, success
          )
        }
        if (identical(mode, "all")) {
          for (m in seq_len(M_s[[s]][l])) {
            before <- .driver_calibration_scale(mu_q_xi[[s]][[l]][, m])
            scale <- before
            success <- is.finite(scale)
            if (success) {
              mu_q_xi[[s]][[l]][, m] <-
                mu_q_xi[[s]][[l]][, m] / scale
            }
            add_diagnostic(
              "specific_score_rms", s, l, m,
              before,
              .driver_calibration_scale(mu_q_xi[[s]][[l]][, m]),
              scale, success
            )
          }
          before <- .driver_calibration_scale(mu_q_normal_b[[s]][, l])
          scale <- before
          success <- is.finite(scale)
          if (success) {
            mu_q_normal_b[[s]][, l] <-
              mu_q_normal_b[[s]][, l] / scale
          }
          add_diagnostic(
            "specific_slab_rms", s, l, NA_integer_,
            before, .driver_calibration_scale(mu_q_normal_b[[s]][, l]),
            scale, success
          )
        }
      }
    }
  }

  mu_q_a <- mu_q_gamma_a * mu_q_normal_a
  term_a <- (Sigma_q_normal_a + mu_q_normal_a^2) * mu_q_gamma_a
  if (.has_specific(L_s)) {
    mu_q_b_specific <- lapply(seq_len(S), function(s) {
      mu_q_gamma_b[[s]] * mu_q_normal_b[[s]]
    })
    term_b_specific <- lapply(seq_len(S), function(s) {
      (Sigma_q_normal_b[[s]] + mu_q_normal_b[[s]]^2) *
        mu_q_gamma_b[[s]]
    })
  } else {
    mu_q_b_specific <- term_b_specific <- NULL
  }
  diagnostic_table <- if (length(diagnostics)) {
    result <- do.call(rbind, diagnostics)
    rownames(result) <- NULL
    result
  } else {
    data.frame()
  }
  environment <- .get_driver_diagnostic_environment()
  if (!is.null(environment)) {
    environment$calibration_diagnostics <- diagnostic_table
  }
  create_named_list(
    mu_q_nu_phi, mu_q_zeta, mu_q_normal_a, mu_q_a, term_a,
    mu_q_nu_psi, mu_q_xi, mu_q_normal_b,
    mu_q_b_specific, term_b_specific,
    diagnostics = diagnostic_table
  )
}

#' Run a fit with private first-sweep driver diagnostics
#'
#' This helper is intentionally unexported.  `dense_gate_sweeps = r` holds the
#' inclusion coordinates at gamma = 1 during total sweeps 1,...,r, suppresses
#' stopping during that window, and then releases the original gamma update.
#'
#' @keywords internal
.bayesSYNC_multi_driver_diagnostic <- function(..., control = NULL) {
  environment <- .new_driver_diagnostic_environment(control)
  option_name <- .multiFSYNC_driver_diagnostic_option
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

  options(stats::setNames(list(environment), option_name))
  fit <- bayesSYNC_multi(...)
  trace <- if (length(environment$rows)) {
    result <- do.call(rbind, environment$rows)
    rownames(result) <- NULL
    result
  } else {
    data.frame()
  }
  dense <- environment$control$dense_gate_sweeps
  attr(trace, "diagnostic_control") <- environment$control
  attr(trace, "dense_gate_released") <- dense == 0L || fit$i_iter > dense
  attr(trace, "post_release_sweeps") <- max(0L, fit$i_iter - dense)
  fit$driver_trace <- trace
  fit$driver_diagnostic_control <- environment$control
  fit$random_scale_calibration_diagnostics <-
    environment$calibration_diagnostics
  fit
}
