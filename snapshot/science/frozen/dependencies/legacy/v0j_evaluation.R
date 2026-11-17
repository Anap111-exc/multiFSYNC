v0j_bind_rows <- function(rows, empty = data.frame()) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

v0j_weighted_norm <- function(value, weights) {
  value <- as.numeric(value)
  weights <- as.numeric(weights)
  sqrt(sum(weights * value^2))
}

v0j_normalize_function <- function(value, weights) {
  value <- as.numeric(value)
  norm <- v0j_weighted_norm(value, weights)
  if (!is.finite(norm) || norm <= .Machine$double.eps) {
    return(rep(0, length(value)))
  }
  value / norm
}

v0j_weighted_ise <- function(estimate, truth, weights) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
  weights <- as.numeric(weights)
  if (length(estimate) != length(truth) || length(truth) != length(weights) ||
      !all(is.finite(c(estimate, truth, weights)))) return(NA_real_)
  sum(weights * (estimate - truth)^2)
}

v0j_project_function <- function(
    truth_value, truth_grid, fit_grid, fit_design) {
  truth_value <- as.numeric(truth_value)
  truth_grid <- as.numeric(truth_grid)
  fit_grid <- as.numeric(fit_grid)
  fit_design <- as.matrix(fit_design)
  target_fit <- as.numeric(stats::approx(
    truth_grid, truth_value, xout = fit_grid, rule = 2, ties = mean
  )$y)
  fit_weights <- v0h_trapezoid_weights(fit_grid)
  weighted_design <- sqrt(fit_weights) * fit_design
  coefficient <- qr.solve(weighted_design, sqrt(fit_weights) * target_fit)
  projected_fit <- as.numeric(fit_design %*% coefficient)
  projected_truth <- as.numeric(stats::approx(
    fit_grid, projected_fit, xout = truth_grid, rule = 2, ties = mean
  )$y)
  list(
    coefficient = coefficient,
    fit_grid_value = projected_fit,
    truth_grid_value = projected_truth
  )
}

v0j_factor_eigenvalues <- function(fit, study, block, factor, effective) {
  value <- if (block == "shared") {
    fit$list_eigenvalues[[factor]]
  } else {
    fit$list_eigenvalues_spec[[study]][[factor]]
  }
  value <- as.numeric(value)
  effective <- max(0L, min(as.integer(effective), length(value)))
  if (effective) value[seq_len(effective)] else numeric()
}

v0j_true_eigenvalues <- function(data, study, block, factor) {
  if (block == "shared") {
    as.numeric(data$true_params$lambda_phi_true[[factor]])
  } else {
    as.numeric(data$true_params$lambda_psi_true[[study]][[factor]])
  }
}

v0j_feature_component_errors <- function(fit, data, evaluation) {
  factors <- evaluation$functional_recovery
  factors <- factors[factors$scope == "ppi_selected", , drop = FALSE]
  truth_grid <- as.numeric(data$true_params$t_grid_dense)
  truth_weights <- v0h_trapezoid_weights(truth_grid)
  rows <- list()
  row_index <- 0L
  for (index in seq_len(nrow(factors))) {
    factor_row <- factors[index, , drop = FALSE]
    study <- as.integer(factor_row$study)
    true_block <- as.character(factor_row$true_role)
    true_factor <- as.integer(factor_row$true_factor)
    truth_object <- v0h_true_function_and_score(
      data, study, true_block, true_factor
    )
    truth_functions <- as.matrix(truth_object$functions)
    matched_factor <- isTRUE(factor_row$matched_by_loading[[1L]])
    estimated_object <- NULL
    component_match <- NULL
    if (matched_factor) {
      estimated_object <- v0h_estimated_function_and_score(
        fit, study, as.character(factor_row$estimated_block),
        as.integer(factor_row$estimated_factor)
      )
      interpolated <- v0h_interpolate_matrix(
        fit$time_g, estimated_object$functions, truth_grid
      )
      similarity <- matrix(
        0, nrow = ncol(truth_functions), ncol = ncol(interpolated)
      )
      if (length(similarity)) {
        for (truth_component in seq_len(nrow(similarity))) {
          for (estimated_component in seq_len(ncol(similarity))) {
            similarity[truth_component, estimated_component] <-
              v0h_weighted_cosine(
                truth_functions[, truth_component],
                interpolated[, estimated_component], truth_weights
              )
          }
        }
      }
      component_match <- v0h_one_to_one_match(similarity, threshold = 0)
    }
    for (truth_component in seq_len(ncol(truth_functions))) {
      row_index <- row_index + 1L
      truth_function <- v0j_normalize_function(
        truth_functions[, truth_component], truth_weights
      )
      projection <- v0j_project_function(
        truth_function, truth_grid, fit$time_g, fit$C_g
      )
      projected <- v0j_normalize_function(
        projection$truth_grid_value, truth_weights
      )
      projection_sign <- sign(sum(truth_weights * truth_function * projected))
      if (!is.finite(projection_sign) || projection_sign == 0) projection_sign <- 1
      projected <- projection_sign * projected
      projection_floor <- v0j_weighted_ise(
        projected, truth_function, truth_weights
      )
      estimated_component <- if (matched_factor) {
        component_match$assignment[[truth_component]]
      } else 0L
      component_matched <- matched_factor && estimated_component > 0L
      if (component_matched) {
        estimated_fit <- estimated_object$functions[, estimated_component]
        estimated_truth <- as.numeric(stats::approx(
          fit$time_g, estimated_fit, xout = truth_grid,
          rule = 2, ties = mean
        )$y)
        estimated_truth <- v0j_normalize_function(
          estimated_truth, truth_weights
        )
        estimated_sign <- sign(sum(
          truth_weights * truth_function * estimated_truth
        ))
        if (!is.finite(estimated_sign) || estimated_sign == 0) estimated_sign <- 1
        estimated_truth <- estimated_sign * estimated_truth
        total_ise <- v0j_weighted_ise(
          estimated_truth, truth_function, truth_weights
        )
        estimate_to_projection_ise <- v0j_weighted_ise(
          estimated_truth, projected, truth_weights
        )
        absolute_inner <- abs(sum(
          truth_weights * truth_function * estimated_truth
        ))
      } else {
        estimated_component <- NA_integer_
        estimated_sign <- NA_real_
        total_ise <- NA_real_
        estimate_to_projection_ise <- NA_real_
        absolute_inner <- 0
      }
      row_index_now <- row_index
      rows[[row_index_now]] <- data.frame(
        scope = "ppi_selected", study = study,
        true_id = as.character(factor_row$true_id),
        true_role = true_block, true_factor = true_factor,
        loading_structure_status =
          as.character(factor_row$loading_structure_status),
        candidate_id = as.character(factor_row$candidate_id),
        estimated_block = as.character(factor_row$estimated_block),
        estimated_factor = suppressWarnings(as.integer(
          factor_row$estimated_factor
        )),
        true_component = truth_component,
        estimated_component = estimated_component,
        factor_matched = matched_factor,
        component_matched = component_matched,
        absolute_inner_product = absolute_inner,
        function_sign = estimated_sign,
        total_sign_aligned_ise = total_ise,
        projection_floor_ise = projection_floor,
        estimate_to_projection_ise = estimate_to_projection_ise,
        stringsAsFactors = FALSE
      )
    }
  }
  v0j_bind_rows(rows)
}

v0j_kernel_hs_norm_squared <- function(kernel, weights) {
  kernel <- as.matrix(kernel)
  weights <- as.numeric(weights)
  if (nrow(kernel) != length(weights) || ncol(kernel) != length(weights)) {
    stop("kernel and weights have incompatible dimensions")
  }
  sum((sqrt(weights) * kernel * rep(sqrt(weights), each = nrow(kernel)))^2)
}

v0j_factor_kernel_errors <- function(fit, data, evaluation) {
  factors <- evaluation$functional_recovery
  factors <- factors[factors$scope == "ppi_selected", , drop = FALSE]
  truth_grid <- as.numeric(data$true_params$t_grid_dense)
  weights <- v0h_trapezoid_weights(truth_grid)
  rows <- lapply(seq_len(nrow(factors)), function(index) {
    factor_row <- factors[index, , drop = FALSE]
    study <- as.integer(factor_row$study)
    true_block <- as.character(factor_row$true_role)
    true_factor <- as.integer(factor_row$true_factor)
    truth_object <- v0h_true_function_and_score(
      data, study, true_block, true_factor
    )
    truth_functions <- as.matrix(truth_object$functions)
    truth_eigenvalues <- v0j_true_eigenvalues(
      data, study, true_block, true_factor
    )
    truth_count <- min(ncol(truth_functions), length(truth_eigenvalues))
    truth_functions <- truth_functions[, seq_len(truth_count), drop = FALSE]
    truth_eigenvalues <- truth_eigenvalues[seq_len(truth_count)]
    truth_kernel <- truth_functions %*% diag(truth_eigenvalues, truth_count) %*%
      t(truth_functions)
    projected_functions <- vapply(seq_len(truth_count), function(component) {
      projection <- v0j_project_function(
        truth_functions[, component], truth_grid, fit$time_g, fit$C_g
      )
      v0j_normalize_function(projection$truth_grid_value, weights)
    }, numeric(length(truth_grid)))
    if (!is.matrix(projected_functions)) {
      projected_functions <- matrix(projected_functions, ncol = truth_count)
    }
    projected_kernel <- projected_functions %*%
      diag(truth_eigenvalues, truth_count) %*% t(projected_functions)
    truth_norm_sq <- v0j_kernel_hs_norm_squared(truth_kernel, weights)
    projection_error_sq <- v0j_kernel_hs_norm_squared(
      projected_kernel - truth_kernel, weights
    )
    matched <- isTRUE(factor_row$matched_by_loading[[1L]])
    if (matched) {
      estimated_block <- as.character(factor_row$estimated_block)
      estimated_factor <- as.integer(factor_row$estimated_factor)
      estimated <- v0h_estimated_function_and_score(
        fit, study, estimated_block, estimated_factor
      )
      estimated_functions <- v0h_interpolate_matrix(
        fit$time_g, estimated$functions, truth_grid
      )
      estimated_eigenvalues <- v0j_factor_eigenvalues(
        fit, study, estimated_block, estimated_factor,
        estimated$effective_M
      )
      estimated_count <- min(
        ncol(estimated_functions), length(estimated_eigenvalues)
      )
      if (estimated_count) {
        estimated_functions <- estimated_functions[
          , seq_len(estimated_count), drop = FALSE
        ]
        estimated_kernel <- estimated_functions %*%
          diag(estimated_eigenvalues[seq_len(estimated_count)], estimated_count) %*%
          t(estimated_functions)
        total_error_sq <- v0j_kernel_hs_norm_squared(
          estimated_kernel - truth_kernel, weights
        )
        projection_distance_sq <- v0j_kernel_hs_norm_squared(
          estimated_kernel - projected_kernel, weights
        )
      } else {
        total_error_sq <- projection_distance_sq <- NA_real_
      }
    } else {
      estimated_count <- 0L
      total_error_sq <- projection_distance_sq <- NA_real_
    }
    data.frame(
      study = study, true_id = as.character(factor_row$true_id),
      true_role = true_block, true_factor = true_factor,
      loading_structure_status =
        as.character(factor_row$loading_structure_status),
      candidate_id = as.character(factor_row$candidate_id),
      factor_matched = matched, true_M = truth_count,
      estimated_effective_M = estimated_count,
      truth_kernel_hs_norm = sqrt(truth_norm_sq),
      kernel_relative_ise = if (is.finite(total_error_sq)) {
        total_error_sq / max(truth_norm_sq, .Machine$double.eps)
      } else NA_real_,
      kernel_relative_hs_error = if (is.finite(total_error_sq)) {
        sqrt(total_error_sq / max(truth_norm_sq, .Machine$double.eps))
      } else NA_real_,
      kernel_projection_floor_relative_ise =
        projection_error_sq / max(truth_norm_sq, .Machine$double.eps),
      kernel_estimate_to_projection_relative_ise =
        if (is.finite(projection_distance_sq)) {
          projection_distance_sq / max(truth_norm_sq, .Machine$double.eps)
        } else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  v0j_bind_rows(rows)
}

v0j_dense_signal_metrics <- function(fit, data) {
  truth <- data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  true_L_s <- if (!is.null(truth$L_s_by_study)) {
    as.integer(unlist(truth$L_s_by_study, use.names = FALSE))
  } else if (length(truth$L_s) == truth$S) {
    as.integer(unlist(truth$L_s, use.names = FALSE))
  } else {
    rep.int(as.integer(truth$L_s), truth$S)
  }
  if (length(true_L_s) != truth$S || any(!is.finite(true_L_s))) {
    stop("Could not recover study-specific true factor counts.")
  }
  true_mean_amplitude <- lapply(seq_len(truth$S), function(study) {
    vapply(seq_len(truth$p), function(variable) {
      numerator <- 0
      denominator <- 0
      phase <- truth$phase_mu_mat[study, variable]
      for (subject in seq_len(truth$n_s[[study]])) {
        time <- data$time_obs[[study]][[subject]]
        base <- sin(2 * pi * time + phase) +
          0.5 * sin(4 * pi * time + 2 * phase)
        saved <- truth$mu_true_values[[study]][[subject]][, variable]
        numerator <- numerator + sum(base * saved)
        denominator <- denominator + sum(base^2)
      }
      if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
        stop("Could not recover frozen mean-function amplitude.")
      }
      numerator / denominator
    }, numeric(1))
  })
  shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
  specific_mask <- lapply(
    fit$factor_ppi_specific, function(value) as.numeric(value >= 0.5)
  )
  accumulators <- list(
    all = c(sse = 0, truth = 0, n = 0),
    selected = c(sse = 0, truth = 0, n = 0),
    mean_only = c(sse = 0, truth = 0, n = 0)
  )
  subject_offset <- 0L
  for (study in seq_len(truth$S)) {
    estimated_mean <- vapply(seq_len(truth$p), function(variable) {
      as.numeric(stats::approx(
        fit$time_g, fit$list_mu_hat[[study]][[variable]],
        xout = grid, rule = 2, ties = mean
      )$y)
    }, numeric(length(grid)))
    true_mean <- vapply(seq_len(truth$p), function(variable) {
      phase <- truth$phase_mu_mat[study, variable]
      true_mean_amplitude[[study]][[variable]] * (
        sin(2 * pi * grid + phase) +
          0.5 * sin(4 * pi * grid + 2 * phase)
      )
    }, numeric(length(grid)))
    saved_mean_error <- max(vapply(seq_len(truth$n_s[[study]]), function(subject) {
      time <- data$time_obs[[study]][[subject]]
      expected <- vapply(seq_len(truth$p), function(variable) {
        phase <- truth$phase_mu_mat[study, variable]
        true_mean_amplitude[[study]][[variable]] * (
          sin(2 * pi * time + phase) +
            0.5 * sin(4 * pi * time + 2 * phase)
        )
      }, numeric(length(time)))
      max(abs(expected - truth$mu_true_values[[study]][[subject]]))
    }, numeric(1)))
    if (!is.finite(saved_mean_error) || saved_mean_error > 1e-10) {
      stop("Frozen dense mean reconstruction assumption failed.")
    }
    estimated_shared_functions <- lapply(seq_len(fit$L_f), function(factor) {
      v0h_interpolate_matrix(
        fit$time_g, fit$list_Phi_hat[[factor]], grid
      )
    })
    estimated_specific_functions <- lapply(
      seq_along(fit$list_Phi_hat_spec[[study]]), function(factor) {
        v0h_interpolate_matrix(
          fit$time_g, fit$list_Phi_hat_spec[[study]][[factor]], grid
        )
      }
    )
    for (subject in seq_len(truth$n_s[[study]])) {
      global_subject <- subject_offset + subject
      estimated_shared_path <- matrix(
        0, nrow = length(grid), ncol = fit$L_f
      )
      for (factor in seq_len(fit$L_f)) {
        estimated_shared_path[, factor] <-
          estimated_shared_functions[[factor]] %*%
          fit$list_Zeta_hat[[factor]][global_subject, ]
      }
      estimated_specific_path <- matrix(
        0, nrow = length(grid),
        ncol = length(estimated_specific_functions)
      )
      if (ncol(estimated_specific_path)) {
        for (factor in seq_len(ncol(estimated_specific_path))) {
          estimated_specific_path[, factor] <-
            estimated_specific_functions[[factor]] %*%
            fit$list_Zeta_hat_spec[[study]][[factor]][subject, ]
        }
      }
      true_shared_path <- matrix(0, length(grid), truth$L_f)
      for (factor in seq_len(truth$L_f)) {
        true_shared_path[, factor] <-
          truth$theta_dense_list[[factor]] %*%
          truth$eta_true[[study]][[factor]][subject, ]
      }
      true_specific_path <- matrix(0, length(grid), true_L_s[[study]])
      if (ncol(true_specific_path)) {
        for (factor in seq_len(ncol(true_specific_path))) {
          true_specific_path[, factor] <-
            truth$kappa_dense_list[[study]][[factor]] %*%
            truth$chi_true[[study]][[factor]][subject, ]
        }
      }
      true_signal <- true_mean +
        true_shared_path %*% t(truth$a_true) +
        true_specific_path %*% t(truth$b_true[[study]])
      estimate_all <- estimated_mean +
        estimated_shared_path %*% t(fit$mu_q_a_hat) +
        estimated_specific_path %*% t(fit$mu_q_b_specific_hat[[study]])
      estimate_selected <- estimated_mean +
        estimated_shared_path %*%
          t(sweep(fit$mu_q_a_hat, 2L, shared_mask, "*")) +
        estimated_specific_path %*%
          t(sweep(
            fit$mu_q_b_specific_hat[[study]], 2L,
            specific_mask[[study]], "*"
          ))
      estimates <- list(
        all = estimate_all,
        selected = estimate_selected,
        mean_only = estimated_mean
      )
      for (name in names(estimates)) {
        error <- estimates[[name]] - true_signal
        accumulators[[name]][["sse"]] <-
          accumulators[[name]][["sse"]] + sum(error^2)
        accumulators[[name]][["truth"]] <-
          accumulators[[name]][["truth"]] + sum(true_signal^2)
        accumulators[[name]][["n"]] <-
          accumulators[[name]][["n"]] + length(true_signal)
      }
    }
    subject_offset <- subject_offset + truth$n_s[[study]]
  }
  rows <- lapply(names(accumulators), function(name) {
    value <- accumulators[[name]]
    data.frame(
      scope = name,
      grid_size = length(grid),
      squared_error = value[["sse"]],
      truth_squared_energy = value[["truth"]],
      scalar_count = value[["n"]],
      rmse = sqrt(value[["sse"]] / value[["n"]]),
      truth_rms = sqrt(value[["truth"]] / value[["n"]]),
      nrmse = sqrt(value[["sse"]] / max(value[["truth"]], 1e-12)),
      relative_mise = value[["sse"]] / max(value[["truth"]], 1e-12),
      stringsAsFactors = FALSE
    )
  })
  v0j_bind_rows(rows)
}

v0j_auc <- function(label, score) {
  label <- as.integer(label)
  score <- as.numeric(score)
  keep <- is.finite(score) & label %in% c(0L, 1L)
  label <- label[keep]
  score <- score[keep]
  positive <- sum(label == 1L)
  negative <- sum(label == 0L)
  if (!positive || !negative) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[label == 1L]) - positive * (positive + 1) / 2) /
    (positive * negative)
}

v0j_loading_recovery <- function(fit, data, evaluation) {
  matches <- evaluation$global_direction_matches
  matches <- matches[matches$scope == "ppi_selected", , drop = FALSE]
  rows <- lapply(seq_len(nrow(matches)), function(index) {
    row <- matches[index, , drop = FALSE]
    truth_loading <- if (row$true_block == "shared") {
      data$true_params$a_true[, row$true_factor]
    } else {
      data$true_params$b_true[[row$true_study]][, row$true_factor]
    }
    truth_support <- if (row$true_block == "shared") {
      data$true_params$gamma_a_true[, row$true_factor]
    } else {
      data$true_params$gamma_b_true[[row$true_study]][, row$true_factor]
    }
    matched <- isTRUE(row$matched[[1L]])
    if (matched && row$candidate_block == "shared") {
      estimate <- fit$mu_q_a_hat[, row$candidate_factor]
      ppi <- fit$mu_q_gamma_a_hat[, row$candidate_factor]
    } else if (matched) {
      estimate <- fit$mu_q_b_specific_hat[[row$candidate_study]][,
        row$candidate_factor]
      ppi <- fit$mu_q_gamma_b_hat[[row$candidate_study]][,
        row$candidate_factor]
    } else {
      estimate <- ppi <- rep(NA_real_, length(truth_loading))
    }
    if (matched) {
      sign_value <- v0h_loading_sign(estimate, truth_loading)
      aligned <- sign_value * estimate
      relative_error <- sqrt(sum((aligned - truth_loading)^2)) /
        max(sqrt(sum(truth_loading^2)), 1e-12)
      norm_ratio <- sqrt(sum(estimate^2)) /
        max(sqrt(sum(truth_loading^2)), 1e-12)
      support_auc <- v0j_auc(truth_support, ppi)
    } else {
      sign_value <- relative_error <- norm_ratio <- support_auc <- NA_real_
    }
    data.frame(
      true_id = row$true_id, true_block = row$true_block,
      true_study = row$true_study, candidate_id = row$candidate_id,
      candidate_block = row$candidate_block,
      candidate_study = row$candidate_study,
      status = row$status, matched = matched,
      loading_abs_cosine = row$loading_abs_cosine,
      loading_sign = sign_value,
      loading_relative_l2_error = relative_error,
      loading_norm_ratio = norm_ratio,
      loading_support_ppi_auc = support_auc,
      stringsAsFactors = FALSE
    )
  })
  v0j_bind_rows(rows)
}

v0j_threshold_sensitivity <- function(fit, data,
                                      thresholds = c(0.6, 0.7, 0.8, 0.9)) {
  truth_global <- v0h_truth_global_catalog(data)
  rows <- list()
  index <- 0L
  for (scope in c("ppi_selected", "all_candidates")) {
    candidate <- v0h_candidate_global_catalog(fit, scope, 0.5)
    for (threshold in thresholds) {
      global <- v0h_classify_catalog_match(
        truth_global, candidate, scope, "global_direction", 0L,
        threshold, threshold
      )
      index <- index + 1L
      rows[[index]] <- data.frame(
        scope = scope, threshold = threshold, unit = "global_direction",
        study = 0L,
        correct = sum(global$truth$status == "correct"),
        misplaced = sum(global$truth$status == "misplaced"),
        missing = sum(global$truth$status == "missing"),
        extra = sum(global$candidate$status == "extra"),
        duplicate = sum(global$candidate$status == "duplicate"),
        stringsAsFactors = FALSE
      )
      for (study in seq_len(fit$S)) {
        catalogs <- v0h_study_catalogs(truth_global, candidate, study)
        local <- v0h_classify_catalog_match(
          catalogs$truth, catalogs$candidate, scope, "study_role", study,
          threshold, threshold
        )
        index <- index + 1L
        rows[[index]] <- data.frame(
          scope = scope, threshold = threshold, unit = "study_role",
          study = study,
          correct = sum(local$truth$status == "correct"),
          misplaced = sum(local$truth$status == "misplaced"),
          missing = sum(local$truth$status == "missing"),
          extra = sum(local$candidate$status == "extra"),
          duplicate = sum(local$candidate$status == "duplicate"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  v0j_bind_rows(rows)
}

v0j_unthresholded_assignment <- function(fit, data) {
  truth_global <- v0h_truth_global_catalog(data)
  rows <- list()
  index <- 0L
  for (scope in c("ppi_selected", "all_candidates")) {
    candidate <- v0h_candidate_global_catalog(fit, scope, 0.5)
    units <- c(list(list(
      unit = "global_direction", study = 0L,
      truth = truth_global, candidate = candidate
    )), lapply(seq_len(fit$S), function(study) {
      catalogs <- v0h_study_catalogs(truth_global, candidate, study)
      list(unit = "study_role", study = study,
           truth = catalogs$truth, candidate = catalogs$candidate)
    }))
    for (unit in units) {
      similarity <- v0h_catalog_similarity(unit$truth, unit$candidate)
      assignment <- v0h_one_to_one_match(similarity, threshold = 0)
      for (truth_index in seq_len(nrow(unit$truth$table))) {
        candidate_index <- assignment$assignment[[truth_index]]
        matched <- candidate_index > 0L
        index <- index + 1L
        rows[[index]] <- data.frame(
          scope = scope, unit = unit$unit, study = unit$study,
          true_id = unit$truth$table$id[[truth_index]],
          true_block = unit$truth$table$block[[truth_index]],
          true_study = unit$truth$table$study[[truth_index]],
          candidate_id = if (matched) {
            unit$candidate$table$id[[candidate_index]]
          } else NA_character_,
          candidate_block = if (matched) {
            unit$candidate$table$block[[candidate_index]]
          } else NA_character_,
          candidate_study = if (matched) {
            unit$candidate$table$study[[candidate_index]]
          } else NA_integer_,
          loading_abs_cosine = if (matched) {
            similarity[truth_index, candidate_index]
          } else 0,
          correct_parameter_block = matched && v0h_correct_block(
            unit$truth$table[truth_index, , drop = FALSE],
            unit$candidate$table[candidate_index, , drop = FALSE],
            global = unit$unit == "global_direction"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  v0j_bind_rows(rows)
}

v0j_ppi_diagnostics <- function(fit) {
  v0j_bind_rows(c(
    list(v0g_ppi_table(
      fit$mu_q_gamma_a_hat, fit$factor_ppi_shared, "shared", 0L
    )),
    lapply(seq_len(fit$S), function(study) {
      v0g_ppi_table(
        fit$mu_q_gamma_b_hat[[study]],
        fit$factor_ppi_specific[[study]], "specific", study
      )
    })
  ))
}

v0j_evaluate_fit <- function(fit, data) {
  required <- c(
    "v0h_evaluate_fit", "v0i_complete_operator_error",
    "v0g_reported_paths_by_study", "v0g_shared_candidate_metrics"
  )
  if (!all(vapply(required, exists, logical(1), mode = "function"))) {
    stop("Source v0G, v0H and v0I evaluation helpers first.")
  }
  primary <- v0h_evaluate_fit(fit, data, v0h_default_thresholds())
  feature_components <- v0j_feature_component_errors(fit, data, primary)
  kernel_errors <- v0j_factor_kernel_errors(fit, data, primary)
  dense_signal <- v0j_dense_signal_metrics(fit, data)
  loading_recovery <- v0j_loading_recovery(fit, data, primary)
  threshold_sensitivity <- v0j_threshold_sensitivity(fit, data)
  unthresholded <- v0j_unthresholded_assignment(fit, data)
  operator <- v0i_complete_operator_error(fit, data)
  paths <- v0g_reported_paths_by_study(fit, data)
  ppi <- v0j_ppi_diagnostics(fit)
  shared_candidate <- v0g_shared_candidate_metrics(fit, paths, 0.10)
  block <- primary$block_subspace
  get_metric <- function(metric, field) {
    v0h_extract_metric(block, "ppi_selected", metric, field)
  }
  dense_selected <- dense_signal[dense_signal$scope == "selected", , drop = FALSE]
  dense_all <- dense_signal[dense_signal$scope == "all", , drop = FALSE]
  matched_component <- feature_components$component_matched
  matched_kernel <- kernel_errors$factor_matched
  summary_extra <- data.frame(
    P_A = get_metric("R_A", "subspace_purity"),
    P_B1 = get_metric("R_B1", "subspace_purity"),
    P_B2 = get_metric("R_B2", "subspace_purity"),
    min_correct_block_R = min(c(
      get_metric("R_A", "subspace_recall"),
      get_metric("R_B1", "subspace_recall"),
      get_metric("R_B2", "subspace_recall")
    )),
    max_cross_block_L = max(c(
      get_metric("L_B1_to_A", "subspace_recall"),
      get_metric("L_B2_to_A", "subspace_recall"),
      get_metric("L_A_to_B1", "subspace_recall"),
      get_metric("L_A_to_B2", "subspace_recall")
    )),
    dense_signal_nrmse_ppi_selected = dense_selected$nrmse,
    dense_signal_relative_mise_ppi_selected = dense_selected$relative_mise,
    dense_signal_nrmse_all_candidates = dense_all$nrmse,
    feature_components_total = nrow(feature_components),
    feature_components_matched = sum(matched_component),
    feature_total_ise_mean = if (any(matched_component)) {
      mean(feature_components$total_sign_aligned_ise[matched_component])
    } else NA_real_,
    feature_projection_floor_ise_mean =
      mean(feature_components$projection_floor_ise),
    feature_estimate_to_projection_ise_mean = if (any(matched_component)) {
      mean(feature_components$estimate_to_projection_ise[matched_component])
    } else NA_real_,
    factor_kernels_total = nrow(kernel_errors),
    factor_kernels_matched = sum(matched_kernel),
    kernel_relative_ise_mean = if (any(matched_kernel)) {
      mean(kernel_errors$kernel_relative_ise[matched_kernel])
    } else NA_real_,
    covariance_operator_relative_error_mean =
      mean(operator$covariance_operator_relative_error),
    loading_relative_l2_error_mean = if (any(loading_recovery$matched)) {
      mean(loading_recovery$loading_relative_l2_error[
        loading_recovery$matched
      ])
    } else NA_real_,
    loading_support_ppi_auc_mean = if (any(loading_recovery$matched)) {
      mean(loading_recovery$loading_support_ppi_auc[
        loading_recovery$matched
      ], na.rm = TRUE)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  list(
    primary = primary,
    feature_component_errors = feature_components,
    factor_kernel_errors = kernel_errors,
    dense_signal_metrics = dense_signal,
    loading_recovery = loading_recovery,
    threshold_sensitivity = threshold_sensitivity,
    unthresholded_assignment = unthresholded,
    covariance_operator_error = operator,
    ppi_diagnostics = ppi,
    shared_candidate_metrics = shared_candidate,
    summary = cbind(primary$summary, summary_extra)
  )
}
