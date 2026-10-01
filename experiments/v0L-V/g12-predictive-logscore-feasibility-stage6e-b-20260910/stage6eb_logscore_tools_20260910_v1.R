# Stage 6E-B: truth-free conditional predictive log-score feasibility tools.
#
# This file is deliberately experiment-local.  It does not change the
# preregistered Stage 6B RMSE/NRMSE interface or the multiFSYNC package API.

s6eb_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

s6eb_as_symmetric <- function(value, name) {
  value <- as.matrix(value)
  s6eb_assert(nrow(value) == ncol(value),
               paste0(name, " must be square."))
  s6eb_assert(all(is.finite(value)),
               paste0(name, " contains a non-finite value."))
  (value + t(value)) / 2
}

s6eb_checked_nonnegative <- function(value, name,
                                     tolerance_multiplier = 1e-9,
                                     reference_scale = NULL) {
  value <- as.numeric(value)
  s6eb_assert(all(is.finite(value)),
               paste0(name, " contains a non-finite value."))
  scale <- abs(value)
  if (!is.null(reference_scale)) {
    reference_scale <- as.numeric(reference_scale)
    s6eb_assert(length(reference_scale) %in% c(1L, length(value)) &&
                 all(is.finite(reference_scale)),
                 paste0(name, " has an invalid numerical reference scale."))
    scale <- pmax(scale, abs(reference_scale))
  }
  tolerance <- tolerance_multiplier * pmax(1, scale)
  s6eb_assert(all(value >= -tolerance),
               paste0(name, " is materially negative."))
  list(value = pmax(value, 0), clamped = sum(value < 0))
}

s6eb_design_at_time <- function(fit, time) {
  s6eb_assert(is.numeric(fit$time_g) && length(fit$time_g) >= 2L,
               "Fit lacks a valid dense time grid.")
  s6eb_assert(is.matrix(fit$C_g) && nrow(fit$C_g) == length(fit$time_g),
               "Fit lacks a spline design aligned with the dense time grid.")
  s6eb_assert(length(time) == 1L && is.finite(time),
               "Prediction time must be one finite number.")
  time_tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(fit$time_g)))
  s6eb_assert(time >= min(fit$time_g) - time_tolerance &&
               time <= max(fit$time_g) + time_tolerance,
               "Prediction time lies outside the fitted dense grid.")
  vapply(seq_len(ncol(fit$C_g)), function(column) {
    stats::approx(
      x = fit$time_g,
      y = fit$C_g[, column],
      xout = time,
      method = "linear",
      rule = 2,
      ties = "ordered"
    )$y
  }, numeric(1L))
}

s6eb_quadratic_variance <- function(design, covariance, name) {
  covariance <- s6eb_as_symmetric(covariance, name)
  s6eb_assert(nrow(covariance) == length(design),
               paste0(name, " has the wrong dimension."))
  checked <- s6eb_checked_nonnegative(
    drop(crossprod(design, covariance %*% design)), name
  )
  checked
}

s6eb_factor_process_moments <- function(score_mean, score_covariance,
                                         coefficient_mean,
                                         coefficient_covariances,
                                         design) {
  score_mean <- as.numeric(score_mean)
  coefficient_mean <- as.matrix(coefficient_mean)
  component_count <- length(score_mean)
  s6eb_assert(component_count >= 1L,
               "A factor process must contain at least one FPCA component.")
  s6eb_assert(ncol(coefficient_mean) == component_count &&
               nrow(coefficient_mean) == length(design),
               "Factor coefficient means have incompatible dimensions.")
  s6eb_assert(is.list(coefficient_covariances) &&
               length(coefficient_covariances) == component_count,
               "Factor coefficient covariance list has the wrong length.")
  score_covariance <- s6eb_as_symmetric(
    score_covariance, "score covariance"
  )
  s6eb_assert(all(dim(score_covariance) == c(component_count,
                                               component_count)),
               "Score covariance has the wrong dimension.")

  function_mean <- as.numeric(design %*% coefficient_mean)
  function_variance_checks <- lapply(seq_len(component_count), function(component) {
    s6eb_quadratic_variance(
      design,
      coefficient_covariances[[component]],
      paste0("coefficient variance for component ", component)
    )
  })
  function_variance <- vapply(
    function_variance_checks, `[[`, numeric(1L), "value"
  )
  function_variance_clamps <- sum(vapply(
    function_variance_checks, `[[`, integer(1L), "clamped"
  ))

  second_function <- tcrossprod(function_mean)
  diag(second_function) <- diag(second_function) + function_variance
  second_score <- score_covariance + tcrossprod(score_mean)
  process_mean <- sum(score_mean * function_mean)
  process_second_moment <- sum(second_score * second_function)
  checked_variance <- s6eb_checked_nonnegative(
    process_second_moment - process_mean^2,
    "factor process posterior variance",
    reference_scale = max(process_second_moment, process_mean^2)
  )

  list(
    mean = process_mean,
    second_moment = process_second_moment,
    variance = checked_variance$value,
    function_mean = function_mean,
    function_variance = function_variance,
    numerical_clamp_count =
      function_variance_clamps + checked_variance$clamped
  )
}

s6eb_factor_contribution <- function(loading_mean, loading_second_moment,
                                      process_moments, label) {
  loading_mean <- as.numeric(loading_mean)
  loading_second_moment <- as.numeric(loading_second_moment)
  s6eb_assert(length(loading_mean) == length(loading_second_moment),
               paste0(label, " loading moments have incompatible dimensions."))
  loading_variance <- s6eb_checked_nonnegative(
    loading_second_moment - loading_mean^2,
    paste0(label, " loading posterior variance"),
    reference_scale = pmax(loading_second_moment, loading_mean^2)
  )
  process_variance <- as.numeric(process_moments$variance)
  low_rank <- if (process_variance > 0) {
    sqrt(process_variance) * loading_mean
  } else {
    numeric(length(loading_mean))
  }
  list(
    mean = loading_mean * process_moments$mean,
    diagonal = loading_variance$value * process_moments$second_moment,
    low_rank = low_rank,
    loading_variance_clamp_count = loading_variance$clamped,
    process_numerical_clamp_count =
      process_moments$numerical_clamp_count
  )
}

s6eb_predictive_moments <- function(fit, study, subject, time) {
  required <- c(
    "S", "n_s", "p", "d", "L_f", "L_s_by_study", "time_g", "C_g",
    "mu_q_nu_mu_original", "Sigma_q_nu_mu_original",
    "mu_q_nu_phi", "Sigma_q_nu_phi", "mu_q_zeta", "Sigma_q_zeta",
    "mu_q_a_original", "term_a_original", "sigsq_eps_hat"
  )
  s6eb_assert(all(required %in% names(fit)),
               "Fit lacks posterior objects required for Stage 6E-B.")
  s6eb_assert(as.integer(fit$d) == 0L,
               "Stage 6E-B v1 supports the registered d=0 route only.")
  study <- as.integer(study)
  subject <- as.integer(subject)
  s6eb_assert(study >= 1L && study <= as.integer(fit$S),
               "Study index is outside the fit.")
  s6eb_assert(subject >= 1L && subject <= as.integer(fit$n_s[[study]]),
               "Subject index is outside the fit.")
  variable_count <- as.integer(fit$p)
  design <- s6eb_design_at_time(fit, time)

  predictive_mean <- numeric(variable_count)
  diagonal_variance <- as.numeric(fit$sigsq_eps_hat[study, ])
  s6eb_assert(length(diagonal_variance) == variable_count &&
               all(is.finite(diagonal_variance)) &&
               all(diagonal_variance > 0),
               "Residual posterior mean variances must be positive and finite.")
  mean_variance_clamp_count <- 0L
  for (variable in seq_len(variable_count)) {
    coefficient_mean <-
      as.numeric(fit$mu_q_nu_mu_original[[study]][[variable]])
    s6eb_assert(length(coefficient_mean) == length(design),
                 "Mean coefficient vector has the wrong dimension.")
    predictive_mean[[variable]] <- sum(design * coefficient_mean)
    mean_variance <- s6eb_quadratic_variance(
      design,
      fit$Sigma_q_nu_mu_original[[study]][[variable]],
      paste0("mean-function variance for variable ", variable)
    )
    diagonal_variance[[variable]] <-
      diagonal_variance[[variable]] + mean_variance$value
    mean_variance_clamp_count <-
      mean_variance_clamp_count + mean_variance$clamped
  }

  low_rank_columns <- list()
  loading_variance_clamp_count <- 0L
  process_numerical_clamp_count <- 0L
  factor_index <- 0L
  add_factor <- function(process_moments, loading_mean,
                         loading_second_moment, label) {
    contribution <- s6eb_factor_contribution(
      loading_mean, loading_second_moment, process_moments, label
    )
    predictive_mean <<- predictive_mean + contribution$mean
    diagonal_variance <<- diagonal_variance + contribution$diagonal
    factor_index <<- factor_index + 1L
    low_rank_columns[[factor_index]] <<- contribution$low_rank
    loading_variance_clamp_count <<-
      loading_variance_clamp_count +
      contribution$loading_variance_clamp_count
    process_numerical_clamp_count <<-
      process_numerical_clamp_count +
      contribution$process_numerical_clamp_count
    invisible(NULL)
  }

  shared_factor_count <- as.integer(fit$L_f)
  if (shared_factor_count > 0L) {
    for (factor in seq_len(shared_factor_count)) {
      process <- s6eb_factor_process_moments(
        score_mean = fit$mu_q_zeta[[study]][[factor]][subject, ],
        score_covariance = fit$Sigma_q_zeta[[study]][[factor]][[subject]],
        coefficient_mean = fit$mu_q_nu_phi[[factor]],
        coefficient_covariances = fit$Sigma_q_nu_phi[[factor]],
        design = design
      )
      add_factor(
        process,
        fit$mu_q_a_original[, factor],
        fit$term_a_original[, factor],
        paste0("shared factor ", factor)
      )
    }
  }

  specific_factor_count <- as.integer(fit$L_s_by_study[[study]])
  if (specific_factor_count > 0L) {
    specific_required <- c(
      "mu_q_nu_psi", "Sigma_q_nu_psi", "mu_q_xi", "Sigma_q_xi",
      "mu_q_b_specific_original", "term_b_specific_original"
    )
    s6eb_assert(all(specific_required %in% names(fit)),
                 "Fit lacks study-specific posterior moments.")
    for (factor in seq_len(specific_factor_count)) {
      process <- s6eb_factor_process_moments(
        score_mean = fit$mu_q_xi[[study]][[factor]][subject, ],
        score_covariance = fit$Sigma_q_xi[[study]][[factor]][[subject]],
        coefficient_mean = fit$mu_q_nu_psi[[study]][[factor]],
        coefficient_covariances =
          fit$Sigma_q_nu_psi[[study]][[factor]],
        design = design
      )
      add_factor(
        process,
        fit$mu_q_b_specific_original[[study]][, factor],
        fit$term_b_specific_original[[study]][, factor],
        paste0("study ", study, " specific factor ", factor)
      )
    }
  }

  low_rank <- if (length(low_rank_columns)) {
    do.call(cbind, low_rank_columns)
  } else {
    matrix(numeric(), nrow = variable_count, ncol = 0L)
  }
  s6eb_assert(all(is.finite(predictive_mean)) &&
               all(is.finite(diagonal_variance)) &&
               all(diagonal_variance > 0) &&
               all(is.finite(low_rank)),
               "Predictive moments are non-finite or non-positive.")

  list(
    mean = predictive_mean,
    diagonal_variance = diagonal_variance,
    low_rank = low_rank,
    marginal_variance = diagonal_variance + rowSums(low_rank^2),
    variance_scope = "full_mean_field_moments_gaussian_approximation",
    exact_posterior_predictive_density = FALSE,
    factor_mode = "all",
    mean_variance_clamp_count = mean_variance_clamp_count,
    loading_variance_clamp_count = loading_variance_clamp_count,
    process_numerical_clamp_count = process_numerical_clamp_count
  )
}

s6eb_lowrank_gaussian_nlpd <- function(observed, mean,
                                        diagonal_variance, low_rank,
                                        variance_floor = 1e-10) {
  observed <- as.numeric(observed)
  mean <- as.numeric(mean)
  diagonal_variance <- as.numeric(diagonal_variance)
  variable_count <- length(observed)
  s6eb_assert(length(mean) == variable_count &&
               length(diagonal_variance) == variable_count,
               "Predictive Gaussian inputs have incompatible dimensions.")
  if (is.null(dim(low_rank))) {
    low_rank <- matrix(as.numeric(low_rank), nrow = variable_count)
  } else {
    low_rank <- as.matrix(low_rank)
  }
  s6eb_assert(nrow(low_rank) == variable_count,
               "Low-rank covariance factor has the wrong row count.")
  s6eb_assert(all(is.finite(c(observed, mean, diagonal_variance, low_rank))),
               "Predictive Gaussian inputs must be finite.")
  s6eb_assert(is.finite(variance_floor) && variance_floor > 0,
               "variance_floor must be positive and finite.")
  s6eb_assert(all(diagonal_variance >= 0),
               "Diagonal predictive variance cannot be negative.")
  floor_count <- sum(diagonal_variance < variance_floor)
  diagonal_used <- pmax(diagonal_variance, variance_floor)
  residual <- observed - mean

  if (ncol(low_rank) == 0L) {
    log_determinant <- sum(log(diagonal_used))
    quadratic <- sum(residual^2 / diagonal_used)
    reduced_log_determinant <- 0
  } else {
    weighted_low_rank <- low_rank / sqrt(diagonal_used)
    reduced <- diag(ncol(low_rank)) + crossprod(weighted_low_rank)
    reduced_cholesky <- chol((reduced + t(reduced)) / 2)
    reduced_log_determinant <- 2 * sum(log(diag(reduced_cholesky)))
    projection <- as.numeric(crossprod(
      low_rank, residual / diagonal_used
    ))
    reduced_solution <- backsolve(
      reduced_cholesky,
      forwardsolve(t(reduced_cholesky), projection)
    )
    log_determinant <- sum(log(diagonal_used)) +
      reduced_log_determinant
    quadratic <- sum(residual^2 / diagonal_used) -
      sum(projection * reduced_solution)
    quadratic_tolerance <- 1e-9 * max(1, sum(residual^2 / diagonal_used))
    s6eb_assert(quadratic >= -quadratic_tolerance,
                 "Woodbury quadratic form is materially negative.")
    quadratic <- max(quadratic, 0)
  }
  total_negative_log_density <- 0.5 * (
    variable_count * log(2 * pi) + log_determinant + quadratic
  )
  list(
    total_negative_log_density = total_negative_log_density,
    nlpd_per_variable = total_negative_log_density / variable_count,
    log_determinant = log_determinant,
    quadratic = quadratic,
    reduced_log_determinant = reduced_log_determinant,
    diagonal_floor_count = floor_count,
    covariance_rank_update = ncol(low_rank)
  )
}

s6eb_dense_gaussian_nlpd <- function(observed, mean,
                                      diagonal_variance, low_rank) {
  observed <- as.numeric(observed)
  mean <- as.numeric(mean)
  diagonal_variance <- as.numeric(diagonal_variance)
  low_rank <- as.matrix(low_rank)
  covariance <- diag(diagonal_variance, nrow = length(diagonal_variance)) +
    tcrossprod(low_rank)
  covariance <- (covariance + t(covariance)) / 2
  cholesky <- chol(covariance)
  residual <- observed - mean
  solution <- backsolve(cholesky, forwardsolve(t(cholesky), residual))
  total_negative_log_density <- 0.5 * (
    length(observed) * log(2 * pi) +
      2 * sum(log(diag(cholesky))) +
      sum(residual * solution)
  )
  list(
    total_negative_log_density = total_negative_log_density,
    nlpd_per_variable = total_negative_log_density / length(observed),
    covariance = covariance
  )
}

s6eb_score_holdout <- function(fit, holdout_plan,
                                variance_floor = 1e-10,
                                interval_level = 0.95) {
  s6eb_assert(is.list(holdout_plan) &&
               identical(holdout_plan$truth_used, FALSE),
               "Invalid or truth-contaminated holdout plan.")
  s6eb_assert(is.list(holdout_plan$records) &&
               length(holdout_plan$records) >= 1L,
               "Holdout plan has no records.")
  s6eb_assert(interval_level > 0 && interval_level < 1,
               "interval_level must lie strictly between zero and one.")
  normal_quantile <- stats::qnorm((1 + interval_level) / 2)

  rows <- lapply(holdout_plan$records, function(record) {
    moments <- s6eb_predictive_moments(
      fit, record$study, record$subject, record$time
    )
    observed <- as.numeric(record$y)
    s6eb_assert(length(observed) == length(moments$mean) &&
                 all(is.finite(observed)),
                 "Held-out response vector is invalid.")
    log_score <- s6eb_lowrank_gaussian_nlpd(
      observed = observed,
      mean = moments$mean,
      diagonal_variance = moments$diagonal_variance,
      low_rank = moments$low_rank,
      variance_floor = variance_floor
    )
    residual <- observed - moments$mean
    marginal_sd <- sqrt(pmax(moments$marginal_variance, variance_floor))
    marginal_nlpd <- mean(
      0.5 * log(2 * pi * marginal_sd^2) +
        0.5 * residual^2 / marginal_sd^2
    )
    covered <- abs(residual) <= normal_quantile * marginal_sd
    data.frame(
      study = as.integer(record$study),
      subject = as.integer(record$subject),
      holdout_time = as.numeric(record$time),
      original_index = as.integer(record$original_index),
      train_time_count = as.integer(record$train_time_count),
      variables = length(observed),
      mse = mean(residual^2),
      mae = mean(abs(residual)),
      observed_second_moment = mean(observed^2),
      joint_nlpd_per_variable = log_score$nlpd_per_variable,
      marginal_nlpd = marginal_nlpd,
      marginal_interval_coverage = mean(covered),
      mean_standardized_squared_error = mean((residual / marginal_sd)^2),
      mean_predictive_sd = mean(marginal_sd),
      minimum_predictive_variance = min(moments$marginal_variance),
      maximum_predictive_variance = max(moments$marginal_variance),
      covariance_rank_update = log_score$covariance_rank_update,
      diagonal_floor_count = log_score$diagonal_floor_count,
      mean_variance_clamp_count = moments$mean_variance_clamp_count,
      loading_variance_clamp_count =
        moments$loading_variance_clamp_count,
      process_numerical_clamp_count =
        moments$process_numerical_clamp_count,
      stringsAsFactors = FALSE
    )
  })
  subject_time <- do.call(rbind, rows)
  mse <- mean(subject_time$mse)
  denominator <- sqrt(mean(subject_time$observed_second_moment))
  aggregate <- data.frame(
    factor_mode = "all",
    subject_time_count = nrow(subject_time),
    variables_per_time_min = min(subject_time$variables),
    variables_per_time_max = max(subject_time$variables),
    heldout_rmse = sqrt(mse),
    heldout_nrmse = if (denominator > 0) sqrt(mse) / denominator else NA_real_,
    heldout_mae = mean(subject_time$mae),
    joint_nlpd_per_variable = mean(subject_time$joint_nlpd_per_variable),
    marginal_nlpd = mean(subject_time$marginal_nlpd),
    marginal_interval_level = interval_level,
    marginal_interval_coverage =
      mean(subject_time$marginal_interval_coverage),
    mean_standardized_squared_error =
      mean(subject_time$mean_standardized_squared_error),
    mean_predictive_sd = mean(subject_time$mean_predictive_sd),
    minimum_predictive_variance =
      min(subject_time$minimum_predictive_variance),
    maximum_predictive_variance =
      max(subject_time$maximum_predictive_variance),
    maximum_covariance_rank_update =
      max(subject_time$covariance_rank_update),
    total_diagonal_floor_count = sum(subject_time$diagonal_floor_count),
    total_numerical_clamp_count = sum(
      subject_time$mean_variance_clamp_count +
        subject_time$loading_variance_clamp_count +
        subject_time$process_numerical_clamp_count
    ),
    predictive_variance_scope =
      "full_mean_field_moments_gaussian_approximation",
    exact_posterior_predictive_density = FALSE,
    feasibility_status = "experimental_not_preregistered_primary",
    replaces_preregistered_rmse = FALSE,
    primary_dimension_score_ready = FALSE,
    truth_used = FALSE,
    stringsAsFactors = FALSE
  )
  list(subject_time = subject_time, aggregate = aggregate)
}
