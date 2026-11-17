v0i_abs_cosine <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  abs(sum(x * y) / denominator)
}

v0i_trapezoid_weights <- function(time) {
  time <- as.numeric(time)
  if (length(time) < 2L || any(!is.finite(time)) || any(diff(time) <= 0)) {
    stop("time must be a finite, strictly increasing grid.")
  }
  delta <- diff(time)
  c(delta[[1L]] / 2,
    (delta[-1L] + delta[-length(delta)]) / 2,
    delta[[length(delta)]] / 2)
}

v0i_kernel_hs_inner <- function(functions_1, functions_2, weights) {
  functions_1 <- as.matrix(functions_1)
  functions_2 <- as.matrix(functions_2)
  if (nrow(functions_1) != nrow(functions_2) ||
      nrow(functions_1) != length(weights)) {
    stop("Function blocks and integration weights have incompatible shapes.")
  }
  cross_gram <- crossprod(functions_1, weights * functions_2)
  sum(cross_gram^2)
}

v0i_kernel_hs_cosine <- function(functions_1, functions_2, weights) {
  numerator <- v0i_kernel_hs_inner(functions_1, functions_2, weights)
  denominator <- sqrt(
    v0i_kernel_hs_inner(functions_1, functions_1, weights) *
      v0i_kernel_hs_inner(functions_2, functions_2, weights)
  )
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}

v0i_matrix_rank <- function(matrix, tolerance = 1e-8) {
  if (!ncol(as.matrix(matrix))) return(0L)
  as.integer(qr(as.matrix(matrix), tol = tolerance)$rank)
}

v0i_truth_audit <- function(data, group = NA_character_, strict = TRUE) {
  if (!is.list(data) || is.null(data$true_params)) {
    stop("data must contain true_params.")
  }
  truth <- data$true_params
  required <- c(
    "a_true", "b_true", "theta_dense_list", "kappa_dense_list",
    "lambda_phi_true", "lambda_psi_true", "t_grid_dense",
    "f_true_values", "g_true_values", "mu_true_values",
    "signal_true_values", "noise_true_values", "sigma2_eps_true",
    "L_f", "L_s", "M_f", "M_s", "S", "n_s", "p", "K"
  )
  missing <- setdiff(required, names(truth))
  if (length(missing)) {
    stop("truth is missing: ", paste(missing, collapse = ", "), ".")
  }
  S <- as.integer(truth$S)
  L_s <- if (length(truth$L_s) == 1L) {
    rep(as.integer(truth$L_s), S)
  } else {
    as.integer(truth$L_s)
  }
  dimensions_ok <-
    S == 2L && identical(as.integer(truth$L_f), 1L) &&
    identical(L_s, c(1L, 1L)) &&
    identical(as.integer(truth$M_f), 2L) &&
    length(truth$M_s) == 2L &&
    all(vapply(truth$M_s, function(x) identical(as.integer(x), 2L), logical(1))) &&
    identical(as.integer(truth$K), 5L) &&
    identical(as.integer(truth$p), 500L) &&
    identical(as.integer(truth$n_s), c(30L, 30L))

  A <- as.matrix(truth$a_true)
  B <- lapply(truth$b_true, as.matrix)
  loading_geometry <- data.frame(
    group = group,
    pair = c("A_vs_B1", "A_vs_B2", "B1_vs_B2"),
    absolute_cosine = c(
      v0i_abs_cosine(A[, 1L], B[[1L]][, 1L]),
      v0i_abs_cosine(A[, 1L], B[[2L]][, 1L]),
      v0i_abs_cosine(B[[1L]][, 1L], B[[2L]][, 1L])
    ),
    stringsAsFactors = FALSE
  )

  study_spaces <- lapply(B, function(block) cbind(A, block))
  intersection_dimension <-
    v0i_matrix_rank(study_spaces[[1L]]) +
    v0i_matrix_rank(study_spaces[[2L]]) -
    v0i_matrix_rank(cbind(study_spaces[[1L]], study_spaces[[2L]]))

  weights <- v0i_trapezoid_weights(truth$t_grid_dense)
  theta <- as.matrix(truth$theta_dense_list[[1L]])
  kappa <- lapply(truth$kappa_dense_list, function(study) {
    as.matrix(study[[1L]])
  })
  kernel_geometry <- data.frame(
    group = group,
    pair = c("shared_vs_specific_1", "shared_vs_specific_2",
             "specific_1_vs_specific_2"),
    hs_cosine = c(
      v0i_kernel_hs_cosine(theta, kappa[[1L]], weights),
      v0i_kernel_hs_cosine(theta, kappa[[2L]], weights),
      v0i_kernel_hs_cosine(kappa[[1L]], kappa[[2L]], weights)
    ),
    stringsAsFactors = FALSE
  )
  operator_geometry <- data.frame(
    group = group,
    study = seq_len(S),
    loading_cosine_squared = loading_geometry$absolute_cosine[1:2]^2,
    kernel_hs_cosine = kernel_geometry$hs_cosine[1:2],
    full_atom_hs_cosine = loading_geometry$absolute_cosine[1:2]^2 *
      kernel_geometry$hs_cosine[1:2],
    stringsAsFactors = FALSE
  )

  eigen_rows <- list(data.frame(
    group = group, block = "shared", study = NA_integer_, factor = 1L,
    component = seq_along(truth$lambda_phi_true[[1L]]),
    eigenvalue = as.numeric(truth$lambda_phi_true[[1L]]),
    stringsAsFactors = FALSE
  ))
  for (study in seq_len(S)) {
    eigen_rows[[length(eigen_rows) + 1L]] <- data.frame(
      group = group, block = "specific", study = study, factor = 1L,
      component = seq_along(truth$lambda_psi_true[[study]][[1L]]),
      eigenvalue = as.numeric(truth$lambda_psi_true[[study]][[1L]]),
      stringsAsFactors = FALSE
    )
  }
  eigenvalues <- do.call(rbind, eigen_rows)
  eigen_gap <- min(vapply(split(eigenvalues$eigenvalue,
                                interaction(eigenvalues$block,
                                            eigenvalues$study,
                                            drop = TRUE)),
                           function(x) if (length(x) > 1L) {
                             min(abs(diff(sort(x, decreasing = TRUE))))
                           } else Inf, numeric(1)))

  contribution_rows <- list()
  for (study in seq_len(S)) {
    shared_process <- unlist(lapply(
      truth$f_true_values[[study]], function(value) value[, 1L]
    ), use.names = FALSE)
    specific_process <- unlist(lapply(
      truth$g_true_values[[study]], function(value) value[, 1L]
    ), use.names = FALSE)
    contribution_rows[[length(contribution_rows) + 1L]] <- data.frame(
      group = group, study = study,
      role = c("shared", "specific"),
      process_rms = c(
        sqrt(mean(shared_process^2)), sqrt(mean(specific_process^2))
      ),
      loading_norm = c(
        sqrt(sum(A[, 1L]^2)), sqrt(sum(B[[study]][, 1L]^2))
      ),
      complete_contribution_rms = c(
        sqrt(mean(shared_process^2) * sum(A[, 1L]^2) / truth$p),
        sqrt(mean(specific_process^2) * sum(B[[study]][, 1L]^2) / truth$p)
      ),
      stringsAsFactors = FALSE
    )
  }
  contribution_energy <- do.call(rbind, contribution_rows)

  reconstruction_error <- 0
  decomposition_error <- 0
  for (study in seq_len(S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      observed <- do.call(cbind, data$Y[[study]][[subject]])
      signal <- truth$signal_true_values[[study]][[subject]]
      noise <- truth$noise_true_values[[study]][[subject]]
      reconstruction_error <- max(
        reconstruction_error, abs(observed - signal - noise)
      )
      reconstructed_signal <- truth$mu_true_values[[study]][[subject]] +
        truth$f_true_values[[study]][[subject]][, 1L, drop = FALSE] %*%
          t(A[, 1L, drop = FALSE]) +
        truth$g_true_values[[study]][[subject]][, 1L, drop = FALSE] %*%
          t(B[[study]][, 1L, drop = FALSE])
      decomposition_error <- max(
        decomposition_error, abs(signal - reconstructed_signal)
      )
    }
  }

  checks <- data.frame(
    group = group,
    check = c(
      "frozen_dimensions", "A_orthogonal_B1", "A_orthogonal_B2",
      "specific_directions_distinct", "common_intersection_is_A",
      "positive_eigenvalues", "eigenvalue_gap",
      "positive_realized_contribution", "observation_reconstruction",
      "signal_decomposition"
    ),
    value = c(
      as.numeric(dimensions_ok),
      loading_geometry$absolute_cosine[[1L]],
      loading_geometry$absolute_cosine[[2L]],
      loading_geometry$absolute_cosine[[3L]],
      intersection_dimension,
      min(eigenvalues$eigenvalue),
      eigen_gap,
      min(contribution_energy$complete_contribution_rms),
      reconstruction_error,
      decomposition_error
    ),
    criterion = c(
      "equal_to_1", "at_most_1e-10", "at_most_1e-10",
      "at_most_0.10", "equal_to_1", "greater_than_0",
      "greater_than_0.05", "greater_than_0.01",
      "at_most_1e-10", "at_most_1e-10"
    ),
    passed = c(
      dimensions_ok,
      loading_geometry$absolute_cosine[[1L]] <= 1e-10,
      loading_geometry$absolute_cosine[[2L]] <= 1e-10,
      loading_geometry$absolute_cosine[[3L]] <= 0.10,
      intersection_dimension == 1L,
      min(eigenvalues$eigenvalue) > 0,
      eigen_gap > 0.05,
      min(contribution_energy$complete_contribution_rms) > 0.01,
      reconstruction_error <= 1e-10,
      decomposition_error <= 1e-10
    ),
    stringsAsFactors = FALSE
  )
  if (strict && !all(checks$passed)) {
    failed <- checks$check[!checks$passed]
    stop("Frozen truth audit failed: ", paste(failed, collapse = ", "), ".")
  }
  list(
    checks = checks,
    loading_geometry = loading_geometry,
    kernel_geometry = kernel_geometry,
    operator_geometry = operator_geometry,
    eigenvalues = eigenvalues,
    contribution_energy = contribution_energy
  )
}

v0i_project_dense_block <- function(block, truth_time, fit_time, C_g) {
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
  } else {
    0
  }
  list(coefficients = coefficients, target = target, fitted = fitted,
       relative_rmse = relative_rmse)
}

v0i_project_mean <- function(data) {
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
    errors[[study]] <- sqrt(squared_error / max(squared_truth, .Machine$double.eps))
  }
  list(coefficients = result, relative_rmse = errors)
}

v0i_random_on_support <- function(truth_loading) {
  truth_loading <- as.numeric(truth_loading)
  active <- which(abs(truth_loading) > 0)
  result <- numeric(length(truth_loading))
  if (!length(active)) return(result)
  draw <- stats::rnorm(length(active))
  draw_norm <- sqrt(sum(draw^2))
  if (!is.finite(draw_norm) || draw_norm <= 0) stop("Degenerate random loading draw.")
  result[active] <- draw / draw_norm * sqrt(sum(truth_loading^2))
  result
}

v0i_near_loading <- function(truth_loading, target_cosine = 0.995) {
  truth_loading <- as.numeric(truth_loading)
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

v0i_set_loading_moments <- function(parameters, shared_loading,
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
    parameters$Sigma_q_normal_b[[study]][, 1L] <-
      specific$normal_variance
    parameters$mu_q_gamma_b[[study]][, 1L] <- specific$probability
    parameters$mu_q_b_specific[[study]][, 1L] <- specific$mixture_mean
    parameters$term_b_specific[[study]][, 1L] <- specific$second_moment
  }
  parameters
}

v0i_small_covariance_list <- function(existing, dimension, variance) {
  lapply(existing, function(value) {
    result <- diag(variance, dimension)
    diagnostic <- attr(value, "solver_diagnostic", exact = TRUE)
    if (!is.null(diagnostic)) attr(result, "solver_diagnostic") <- diagnostic
    result
  })
}

v0i_map_complete_truth <- function(state, data, template_fit,
                                   mode = c("near", "swap"), seed) {
  mode <- match.arg(mode)
  set.seed(as.integer(seed))
  truth <- data$true_params
  parameters <- state$parameters
  mean_projection <- v0i_project_mean(data)
  shared_projection <- v0i_project_dense_block(
    truth$theta_dense_list[[1L]], truth$t_grid_dense,
    template_fit$time_g, template_fit$C_g
  )
  specific_projection <- lapply(seq_len(truth$S), function(study) {
    v0i_project_dense_block(
      truth$kappa_dense_list[[study]][[1L]], truth$t_grid_dense,
      template_fit$time_g, template_fit$C_g
    )
  })
  parameters$mu_q_nu_mu <- mean_projection$coefficients

  if (mode == "near") {
    shared_loading <- v0i_near_loading(truth$a_true[, 1L])
    specific_loadings <- lapply(truth$b_true, function(block) {
      v0i_near_loading(block[, 1L])
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
  } else {
    shared_loading <- truth$b_true[[1L]][, 1L]
    specific_loadings <- list(
      truth$a_true[, 1L], truth$b_true[[2L]][, 1L]
    )
    parameters$mu_q_nu_phi[[1L]] <-
      specific_projection[[1L]]$coefficients
    parameters$mu_q_zeta[[1L]][[1L]] <-
      as.matrix(truth$chi_true[[1L]][[1L]])
    parameters$mu_q_zeta[[2L]][[1L]][] <- 0
    parameters$mu_q_nu_psi[[1L]][[1L]] <-
      shared_projection$coefficients
    parameters$mu_q_xi[[1L]][[1L]] <-
      as.matrix(truth$eta_true[[1L]][[1L]])
    parameters$mu_q_nu_psi[[2L]][[1L]] <-
      specific_projection[[2L]]$coefficients
    parameters$mu_q_xi[[2L]][[1L]] <-
      as.matrix(truth$chi_true[[2L]][[1L]])
  }

  K_total <- nrow(parameters$mu_q_nu_phi[[1L]])
  M_shared <- ncol(parameters$mu_q_nu_phi[[1L]])
  parameters$Sigma_q_nu_phi[[1L]] <- v0i_small_covariance_list(
    parameters$Sigma_q_nu_phi[[1L]], K_total, 1e-3
  )
  for (study in seq_len(truth$S)) {
    parameters$Sigma_q_zeta[[study]][[1L]] <- lapply(
      parameters$Sigma_q_zeta[[study]][[1L]],
      function(value) diag(0.05, M_shared)
    )
    M_specific <- ncol(parameters$mu_q_nu_psi[[study]][[1L]])
    parameters$Sigma_q_nu_psi[[study]][[1L]] <-
      v0i_small_covariance_list(
        parameters$Sigma_q_nu_psi[[study]][[1L]], K_total, 1e-3
      )
    parameters$Sigma_q_xi[[study]][[1L]] <- lapply(
      parameters$Sigma_q_xi[[study]][[1L]],
      function(value) diag(0.05, M_specific)
    )
  }
  parameters <- v0i_set_loading_moments(
    parameters, shared_loading, specific_loadings
  )
  parameters$mu_q_recip_sigsq_eps <- 1 / truth$sigma2_eps_true
  state$parameters <- parameters
  state$diagnostics <- c(state$diagnostics, list(
    v0i_path = if (mode == "near") "near_projected_truth" else
      "reference_study_swap",
    truth_available_to_fit = TRUE,
    projected_mean_relative_rmse = mean_projection$relative_rmse,
    projected_shared_function_relative_rmse =
      shared_projection$relative_rmse,
    projected_specific_function_relative_rmse = vapply(
      specific_projection, `[[`, numeric(1), "relative_rmse"
    )
  ))
  state
}

v0i_make_correct_block_random_state <- function(state, data, seed) {
  set.seed(as.integer(seed))
  truth <- data$true_params
  shared_loading <- v0i_random_on_support(truth$a_true[, 1L])
  specific_loadings <- lapply(truth$b_true, function(block) {
    v0i_random_on_support(block[, 1L])
  })
  state$parameters <- v0i_set_loading_moments(
    state$parameters, shared_loading, specific_loadings
  )
  state$diagnostics <- c(state$diagnostics, list(
    v0i_path = "correct_block_random_numeric",
    truth_available_to_fit = TRUE,
    information_injected = "parameter_block_and_loading_support_only"
  ))
  state
}

v0i_validate_state <- function(multi_env, state, template_state) {
  L_s_by_study <- if (length(state$L_s) == 1L) {
    rep(as.integer(state$L_s), state$S)
  } else {
    as.integer(state$L_s)
  }
  multi_env$.validate_continuation_state(
    state = state,
    template_parameters = template_state$parameters,
    S = state$S, n_s = state$n_s, p = state$p, d = state$d,
    L_f = state$L_f, L_s = L_s_by_study,
    M_f = state$M_f, M_s = state$M_s, K = state$K,
    bool_var_spec_prob = state$bool_var_spec_prob,
    response_center = state$response_center,
    response_scale = state$response_scale,
    time_g = state$time_g
  )
  invisible(TRUE)
}

v0i_initial_state_metrics <- function(state, data, path, group) {
  truth <- data$true_params
  parameters <- state$parameters
  comparisons <- data.frame(
    group = group, path = path,
    estimated_block = c("shared", "shared", "shared",
                        "specific_1", "specific_1", "specific_1",
                        "specific_2", "specific_2", "specific_2"),
    truth_direction = rep(c("A", "B1", "B2"), 3L),
    absolute_cosine = c(
      v0i_abs_cosine(parameters$mu_q_a[, 1L], truth$a_true[, 1L]),
      v0i_abs_cosine(parameters$mu_q_a[, 1L], truth$b_true[[1L]][, 1L]),
      v0i_abs_cosine(parameters$mu_q_a[, 1L], truth$b_true[[2L]][, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[1L]][, 1L], truth$a_true[, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[1L]][, 1L], truth$b_true[[1L]][, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[1L]][, 1L], truth$b_true[[2L]][, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[2L]][, 1L], truth$a_true[, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[2L]][, 1L], truth$b_true[[1L]][, 1L]),
      v0i_abs_cosine(parameters$mu_q_b_specific[[2L]][, 1L], truth$b_true[[2L]][, 1L])
    ),
    stringsAsFactors = FALSE
  )
  ppi <- data.frame(
    group = group, path = path,
    block = c("shared", "specific_1", "specific_2"),
    selected_at_half = c(
      sum(parameters$mu_q_gamma_a[, 1L] >= 0.5),
      sum(parameters$mu_q_gamma_b[[1L]][, 1L] >= 0.5),
      sum(parameters$mu_q_gamma_b[[2L]][, 1L] >= 0.5)
    ),
    expected_active = c(
      sum(parameters$mu_q_gamma_a[, 1L]),
      sum(parameters$mu_q_gamma_b[[1L]][, 1L]),
      sum(parameters$mu_q_gamma_b[[2L]][, 1L])
    ),
    stringsAsFactors = FALSE
  )
  list(comparisons = comparisons, ppi = ppi)
}

v0i_complete_operator_error <- function(fit, data) {
  truth <- data$true_params
  weights <- v0i_trapezoid_weights(fit$time_g)
  truth_shared <- vapply(seq_len(ncol(truth$theta_dense_list[[1L]])),
                         function(component) {
    stats::approx(truth$t_grid_dense,
                  truth$theta_dense_list[[1L]][, component],
                  xout = fit$time_g, rule = 2)$y
  }, numeric(length(fit$time_g)))
  truth_specific <- lapply(seq_len(truth$S), function(study) {
    block <- truth$kappa_dense_list[[study]][[1L]]
    vapply(seq_len(ncol(block)), function(component) {
      stats::approx(truth$t_grid_dense, block[, component],
                    xout = fit$time_g, rule = 2)$y
    }, numeric(length(fit$time_g)))
  })
  estimated_shared <- fit$list_Phi_hat[[1L]] %*%
    diag(sqrt(fit$list_eigenvalues[[1L]]),
         nrow = length(fit$list_eigenvalues[[1L]]))
  estimated_specific <- lapply(seq_len(truth$S), function(study) {
    fit$list_Phi_hat_spec[[study]][[1L]] %*%
      diag(sqrt(fit$list_eigenvalues_spec[[study]][[1L]]),
           nrow = length(fit$list_eigenvalues_spec[[study]][[1L]]))
  })

  atom_inner <- function(load_1, fun_1, load_2, fun_2) {
    sum(load_1 * load_2)^2 * v0i_kernel_hs_inner(fun_1, fun_2, weights)
  }
  rows <- lapply(seq_len(truth$S), function(study) {
    truth_load <- list(truth$a_true[, 1L], truth$b_true[[study]][, 1L])
    truth_fun <- list(truth_shared, truth_specific[[study]])
    estimated_load <- list(fit$mu_q_a_hat[, 1L],
                           fit$mu_q_b_specific_hat[[study]][, 1L])
    estimated_fun <- list(estimated_shared, estimated_specific[[study]])
    truth_norm <- sum(vapply(seq_len(2L), function(i) {
      sum(vapply(seq_len(2L), function(j) {
        atom_inner(truth_load[[i]], truth_fun[[i]],
                   truth_load[[j]], truth_fun[[j]])
      }, numeric(1)))
    }, numeric(1)))
    estimated_norm <- sum(vapply(seq_len(2L), function(i) {
      sum(vapply(seq_len(2L), function(j) {
        atom_inner(estimated_load[[i]], estimated_fun[[i]],
                   estimated_load[[j]], estimated_fun[[j]])
      }, numeric(1)))
    }, numeric(1)))
    cross_inner <- sum(vapply(seq_len(2L), function(i) {
      sum(vapply(seq_len(2L), function(j) {
        atom_inner(truth_load[[i]], truth_fun[[i]],
                   estimated_load[[j]], estimated_fun[[j]])
      }, numeric(1)))
    }, numeric(1)))
    squared_error <- max(0, truth_norm + estimated_norm - 2 * cross_inner)
    data.frame(
      study = study,
      truth_operator_hs_norm = sqrt(truth_norm),
      estimated_operator_hs_norm = sqrt(estimated_norm),
      covariance_operator_relative_error =
        sqrt(squared_error / max(truth_norm, .Machine$double.eps)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
