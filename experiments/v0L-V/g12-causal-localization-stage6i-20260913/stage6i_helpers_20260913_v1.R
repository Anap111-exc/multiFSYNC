s6i_swap_two <- function(value) {
  if (is.null(value)) return(value)
  if (is.list(value) && length(value) == 2L) return(value[c(2L, 1L)])
  if (is.atomic(value) && length(value) == 2L) return(value[c(2L, 1L)])
  stop("Expected a two-study object.", call. = FALSE)
}

s6i_swap_observation <- function(observation) {
  stopifnot(is.list(observation), length(observation$Y) == 2L,
            length(observation$time_obs) == 2L)
  result <- observation
  for (name in c("Y", "time_obs", "Z")) {
    if (!is.null(result[[name]])) result[[name]] <- s6i_swap_two(result[[name]])
  }
  if (!is.null(result$dimensions$n_s)) {
    result$dimensions$n_s <- rev(result$dimensions$n_s)
  }
  result$data_id <- paste0(observation$data_id, "__studies_swapped")
  result$study_permutation <- c(2L, 1L)
  result
}

s6i_swap_truth_bundle <- function(bundle) {
  stopifnot(is.list(bundle), is.list(bundle$data),
            is.list(bundle$data$true_params))
  result <- bundle
  data <- result$data
  for (name in c("Y", "Z", "time_obs", "C")) {
    if (!is.null(data[[name]])) data[[name]] <- s6i_swap_two(data[[name]])
  }
  truth <- data$true_params
  study_list_names <- c(
    "b_true", "gamma_b_true", "omega_b_true", "zeta_true", "xi_true",
    "eta_true", "chi_true", "zeta_fpca_true", "xi_fpca_true",
    "psi_dense_list", "kappa_dense_list", "lambda_psi_true",
    "psi_basis_meta", "M_s", "mu_true_values", "f_true_values",
    "g_true_values", "beta_true_values", "signal_true_values",
    "noise_true_values", "n_obs_by_subject"
  )
  for (name in intersect(study_list_names, names(truth))) {
    value <- truth[[name]]
    if (!is.null(value)) truth[[name]] <- s6i_swap_two(value)
  }
  for (name in intersect(c("n_s", "L_s_by_study"), names(truth))) {
    truth[[name]] <- rev(truth[[name]])
  }
  if (!is.null(truth$phase_mu_mat) && is.matrix(truth$phase_mu_mat)) {
    if (ncol(truth$phase_mu_mat) == 2L) {
      truth$phase_mu_mat <- truth$phase_mu_mat[, c(2L, 1L), drop = FALSE]
    } else if (nrow(truth$phase_mu_mat) == 2L) {
      truth$phase_mu_mat <- truth$phase_mu_mat[c(2L, 1L), , drop = FALSE]
    }
  }
  data$true_params <- truth
  result$data <- data
  result$data_id <- paste0(bundle$data_id, "__studies_swapped")
  result$study_permutation <- c(2L, 1L)
  result
}

s6i_trapezoid_weights <- function(time) {
  time <- as.numeric(time)
  stopifnot(length(time) >= 2L, all(is.finite(time)), all(diff(time) > 0))
  delta <- diff(time)
  c(delta[[1L]] / 2,
    (delta[-1L] + delta[-length(delta)]) / 2,
    delta[[length(delta)]] / 2)
}

s6i_abs_cosine <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  abs(sum(x * y) / denominator)
}

s6i_project_dense_block <- function(block, truth_time, fit_time, C_g) {
  block <- as.matrix(block)
  target <- vapply(seq_len(ncol(block)), function(component) {
    stats::approx(truth_time, block[, component], xout = fit_time,
                  rule = 2)$y
  }, numeric(length(fit_time)))
  if (!is.matrix(target)) target <- matrix(target, ncol = ncol(block))
  coefficients <- qr.solve(C_g, target)
  fitted <- C_g %*% coefficients
  denominator <- sqrt(sum(target^2))
  relative_rmse <- if (denominator > 0) {
    sqrt(sum((fitted - target)^2)) / denominator
  } else 0
  list(coefficients = coefficients, target = target, fitted = fitted,
       relative_rmse = relative_rmse)
}

s6i_project_mean <- function(data) {
  truth <- data$true_params
  result <- vector("list", truth$S)
  errors <- numeric(truth$S)
  for (study in seq_len(truth$S)) {
    design <- do.call(rbind, data$C[[study]])
    result[[study]] <- vector("list", truth$p)
    squared_error <- 0
    squared_truth <- 0
    for (variable in seq_len(truth$p)) {
      target <- unlist(lapply(
        truth$mu_true_values[[study]], function(value) value[, variable]
      ), use.names = FALSE)
      coefficient <- as.numeric(qr.solve(design, target))
      result[[study]][[variable]] <- coefficient
      squared_error <- squared_error + sum((design %*% coefficient - target)^2)
      squared_truth <- squared_truth + sum(target^2)
    }
    errors[[study]] <- sqrt(
      squared_error / max(squared_truth, .Machine$double.eps)
    )
  }
  list(coefficients = result, relative_rmse = errors)
}

s6i_near_loading <- function(truth_loading, target_cosine = 0.995) {
  truth_loading <- as.numeric(truth_loading)
  if (target_cosine >= 1) return(truth_loading)
  active <- which(abs(truth_loading) > 0)
  result <- numeric(length(truth_loading))
  if (length(active) < 2L) return(truth_loading)
  unit_truth <- truth_loading[active] / sqrt(sum(truth_loading[active]^2))
  orthogonal <- stats::rnorm(length(active))
  orthogonal <- orthogonal - sum(orthogonal * unit_truth) * unit_truth
  orthogonal_norm <- sqrt(sum(orthogonal^2))
  if (!is.finite(orthogonal_norm) || orthogonal_norm <= 1e-12) {
    orthogonal <- rev(unit_truth)
    orthogonal <- orthogonal - sum(orthogonal * unit_truth) * unit_truth
    orthogonal_norm <- sqrt(sum(orthogonal^2))
  }
  orthogonal <- orthogonal / orthogonal_norm
  truth_norm <- sqrt(sum(truth_loading[active]^2))
  result[active] <- truth_norm * (
    target_cosine * unit_truth + sqrt(1 - target_cosine^2) * orthogonal
  )
  result
}

s6i_set_loading_moments <- function(parameters, shared_loading,
                                    specific_loadings,
                                    active_probability = 0.99,
                                    inactive_probability = 0.01,
                                    slab_variance = 1e-4) {
  set_one <- function(target) {
    target <- as.numeric(target)
    probability <- ifelse(abs(target) > 0, active_probability,
                          inactive_probability)
    normal_mean <- target / probability
    list(
      normal_mean = normal_mean,
      normal_variance = rep(slab_variance, length(target)),
      probability = probability,
      mixture_mean = target,
      second_moment = probability * (slab_variance + normal_mean^2)
    )
  }
  shared <- set_one(shared_loading)
  parameters$mu_q_normal_a[, 1L] <- shared$normal_mean
  parameters$Sigma_q_normal_a[, 1L] <- shared$normal_variance
  parameters$mu_q_gamma_a[, 1L] <- shared$probability
  parameters$mu_q_a[, 1L] <- shared$mixture_mean
  parameters$term_a[, 1L] <- shared$second_moment
  for (study in seq_along(specific_loadings)) {
    specific <- set_one(specific_loadings[[study]])
    parameters$mu_q_normal_b[[study]][, 1L] <- specific$normal_mean
    parameters$Sigma_q_normal_b[[study]][, 1L] <- specific$normal_variance
    parameters$mu_q_gamma_b[[study]][, 1L] <- specific$probability
    parameters$mu_q_b_specific[[study]][, 1L] <- specific$mixture_mean
    parameters$term_b_specific[[study]][, 1L] <- specific$second_moment
  }
  parameters
}

s6i_small_covariance_list <- function(existing, dimension, variance) {
  lapply(existing, function(value) {
    result <- diag(variance, dimension)
    diagnostic <- attr(value, "solver_diagnostic", exact = TRUE)
    if (!is.null(diagnostic)) attr(result, "solver_diagnostic") <- diagnostic
    result
  })
}

s6i_map_truth_state <- function(state, truth_bundle, template_fit,
                                mode = c("exact", "near"), seed) {
  mode <- match.arg(mode)
  set.seed(as.integer(seed))
  data <- truth_bundle$data
  truth <- data$true_params
  parameters <- state$parameters
  mean_projection <- s6i_project_mean(data)
  shared_projection <- s6i_project_dense_block(
    truth$theta_dense_list[[1L]], truth$t_grid_dense,
    template_fit$time_g, template_fit$C_g
  )
  specific_projection <- lapply(seq_len(truth$S), function(study) {
    s6i_project_dense_block(
      truth$kappa_dense_list[[study]][[1L]], truth$t_grid_dense,
      template_fit$time_g, template_fit$C_g
    )
  })
  parameters$mu_q_nu_mu <- mean_projection$coefficients
  target_cosine <- if (mode == "exact") 1 else 0.995
  shared_loading <- s6i_near_loading(truth$a_true[, 1L], target_cosine)
  specific_loadings <- lapply(truth$b_true, function(block) {
    s6i_near_loading(block[, 1L], target_cosine)
  })
  parameters$mu_q_nu_phi[[1L]] <- shared_projection$coefficients
  parameters$mu_q_zeta <- lapply(seq_len(truth$S), function(study) {
    list(as.matrix(truth$eta_true[[study]][[1L]]))
  })
  for (study in seq_len(truth$S)) {
    parameters$mu_q_nu_psi[[study]][[1L]] <-
      specific_projection[[study]]$coefficients
    parameters$mu_q_xi[[study]][[1L]] <-
      as.matrix(truth$chi_true[[study]][[1L]])
  }
  K_total <- nrow(parameters$mu_q_nu_phi[[1L]])
  M_shared <- ncol(parameters$mu_q_nu_phi[[1L]])
  parameters$Sigma_q_nu_phi[[1L]] <- s6i_small_covariance_list(
    parameters$Sigma_q_nu_phi[[1L]], K_total, 1e-3
  )
  for (study in seq_len(truth$S)) {
    parameters$Sigma_q_zeta[[study]][[1L]] <- lapply(
      parameters$Sigma_q_zeta[[study]][[1L]],
      function(value) diag(0.05, M_shared)
    )
    M_specific <- ncol(parameters$mu_q_nu_psi[[study]][[1L]])
    parameters$Sigma_q_nu_psi[[study]][[1L]] <-
      s6i_small_covariance_list(
        parameters$Sigma_q_nu_psi[[study]][[1L]], K_total, 1e-3
      )
    parameters$Sigma_q_xi[[study]][[1L]] <- lapply(
      parameters$Sigma_q_xi[[study]][[1L]],
      function(value) diag(0.05, M_specific)
    )
  }
  parameters <- s6i_set_loading_moments(
    parameters, shared_loading, specific_loadings
  )
  parameters$mu_q_recip_sigsq_eps <- 1 / truth$sigma2_eps_true
  state$parameters <- parameters
  state$diagnostics <- c(state$diagnostics, list(
    stage6i_path = paste0(mode, "_projected_truth"),
    truth_available_to_fit = TRUE,
    loading_target_cosine = target_cosine,
    projected_mean_relative_rmse = mean_projection$relative_rmse,
    projected_shared_function_relative_rmse =
      shared_projection$relative_rmse,
    projected_specific_function_relative_rmse = vapply(
      specific_projection, `[[`, numeric(1L), "relative_rmse"
    )
  ))
  state
}

s6i_validate_state <- function(state, template_state) {
  validator <- getFromNamespace(".validate_continuation_state", "multiFSYNC")
  L_s <- if (length(state$L_s) == 1L) rep(state$L_s, state$S) else state$L_s
  validator(
    state = state, template_parameters = template_state$parameters,
    S = state$S, n_s = state$n_s, p = state$p, d = state$d,
    L_f = state$L_f, L_s = as.integer(L_s), M_f = state$M_f,
    M_s = state$M_s, K = state$K,
    bool_var_spec_prob = state$bool_var_spec_prob,
    response_center = state$response_center,
    response_scale = state$response_scale, time_g = state$time_g
  )
  invisible(TRUE)
}

s6i_second_moment_energy <- function(mu, Sigma, gram) {
  as.numeric(crossprod(mu, gram %*% mu) + sum(gram * t(Sigma)))
}

s6i_generalized_eigenvalues <- function(numerator, denominator,
                                        ridge_fraction = 1e-10) {
  dimension <- nrow(denominator)
  scale <- sum(diag(denominator)) / dimension
  ridge <- max(.Machine$double.eps, ridge_fraction * scale)
  regularized <- (denominator + t(denominator)) / 2 + diag(ridge, dimension)
  root <- chol(regularized)
  root_inverse <- backsolve(root, diag(dimension))
  transformed <- crossprod(root_inverse, numerator %*% root_inverse)
  sort(eigen((transformed + t(transformed)) / 2,
             symmetric = TRUE, only.values = TRUE)$values,
       decreasing = TRUE)
}
