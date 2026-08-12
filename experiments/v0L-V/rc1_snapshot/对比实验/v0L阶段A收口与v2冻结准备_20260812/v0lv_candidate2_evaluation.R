# Candidate2 evaluation layer. Truth-bearing objects may enter only after the
# truth-free selection marker has been frozen and unseal has been authorized.

V0LV_CANDIDATE2_EVALUATOR_VERSION <-
  "v0lv_candidate2_evaluator_1.0.0-20260812"
V0LV_CANDIDATE2_POOLED_MATCH_TARGET <- "span[A0,B10,B20]"
V0LV_CANDIDATE2_FACTOR_PPI_THRESHOLD <- 0.5
V0LV_CANDIDATE2_BOOTSTRAP_REPLICATES <- 9999L
V0LV_CANDIDATE2_BOOTSTRAP_SEED <- 82640001L

candidate2_evaluation_source_paths <- function(project_root) {
  file.path(project_root, "对比实验", c(
    file.path("v0G中性缩减重拟合_20260803", "v0g_evaluation.R"),
    file.path("v0H载荷优先结构评估_20260803", "v0h_evaluation.R"),
    file.path("v0I真值初始化短诊断_20260805", "v0i_truth_diagnostic.R"),
    file.path("v0J无真值可达性析因_20260808", "v0j_evaluation.R")
  ))
}

candidate2_load_frozen_multi_evaluator <- function(project_root) {
  paths <- candidate2_evaluation_source_paths(project_root)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Frozen multiFSYNC evaluator source is missing: ",
         paste(missing, collapse = ", "))
  }
  environment <- new.env(parent = globalenv())
  for (path in paths) sys.source(path, envir = environment, keep.source = TRUE)
  # Frozen v0J uses vapply(required, exists, ...) without an explicit envir.
  # Under its historical global sourcing this found sibling helpers; preserve
  # that lookup contract locally without exporting evaluator symbols globally.
  evaluator_environment <- environment
  environment$exists <- function(x, where = -1, envir = NULL,
                                 mode = "any", inherits = TRUE) {
    base::exists(x, envir = evaluator_environment,
                 mode = mode, inherits = inherits)
  }
  if (!is.function(environment$v0j_evaluate_fit)) {
    stop("Frozen v0J evaluator did not load v0j_evaluate_fit().")
  }
  attr(environment, "source_paths") <- normalizePath(
    paths, winslash = "/", mustWork = TRUE
  )
  environment
}

candidate2_truth_data <- function(unsealed_truth) {
  if (!is.list(unsealed_truth) ||
      !identical(unsealed_truth$bundle_class,
                 "v0lv_sealed_truth_candidate2") ||
      !is.list(unsealed_truth$data) || is.null(unsealed_truth$data$true_params)) {
    stop("Evaluation requires an authorized, unsealed candidate2 truth bundle.")
  }
  unsealed_truth$data
}

candidate2_evaluate_multi <- function(
    fit, unsealed_truth, project_root, evaluator_environment = NULL,
    mode = c("formal", "smoke"), binding_dir = NULL) {
  mode <- match.arg(mode)
  if (mode == "formal" && !is.null(evaluator_environment)) {
    stop("Formal evaluation forbids evaluator injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (is.null(evaluator_environment)) {
    evaluator_environment <- candidate2_load_frozen_multi_evaluator(project_root)
  }
  data <- candidate2_truth_data(unsealed_truth)
  result <- evaluator_environment$v0j_evaluate_fit(fit, data)
  result$evaluation_contract <- data.frame(
    evaluator_version = V0LV_CANDIDATE2_EVALUATOR_VERSION,
    evaluator = "frozen_v0G_v0H_v0I_v0J",
    factor_ppi_rule = ">=0.5",
    oracle_L_M = TRUE,
    truth_used_for_fit_or_selection = FALSE,
    stringsAsFactors = FALSE
  )
  result
}

candidate2_safe_cosine <- function(x, y, absolute = TRUE) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || !length(x) ||
      any(!is.finite(x)) || any(!is.finite(y))) return(NA_real_)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator)) return(NA_real_)
  # A finite zero loading/function candidate remains one of the three pooled
  # oracle candidates. Its directional similarity is zero, not missing.
  if (denominator <= 0) return(0)
  value <- sum(x * y) / denominator
  if (absolute) abs(value) else value
}

candidate2_safe_cor <- function(x, y, absolute = FALSE) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L ||
      any(!is.finite(x)) || any(!is.finite(y)) ||
      stats::sd(x) <= 0 || stats::sd(y) <= 0) return(NA_real_)
  value <- stats::cor(x, y)
  if (absolute) abs(value) else value
}

candidate2_nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate); truth <- as.numeric(truth)
  if (length(estimate) != length(truth) || !length(truth) ||
      any(!is.finite(estimate)) || any(!is.finite(truth))) return(NA_real_)
  denominator <- sqrt(mean(truth^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sqrt(mean((estimate - truth)^2)) / denominator
}

candidate2_trapezoid <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    stop("Evaluation grid must be finite and strictly increasing.")
  }
  delta <- diff(grid)
  c(delta[[1L]] / 2,
    (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}

candidate2_interpolate_columns <- function(grid, value, target) {
  value <- as.matrix(value)
  if (!ncol(value)) return(matrix(numeric(), nrow = length(target), ncol = 0L))
  result <- vapply(seq_len(ncol(value)), function(column) {
    stats::approx(grid, value[, column], xout = target,
                  rule = 2, ties = mean)$y
  }, numeric(length(target)))
  matrix(result, nrow = length(target), ncol = ncol(value))
}

candidate2_permutations <- function(values) {
  values <- as.integer(values)
  if (length(values) <= 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(index) {
    rest <- candidate2_permutations(values[-index])
    cbind(values[[index]], rest)
  }))
}

candidate2_pooled_truth_catalog <- function(data) {
  truth <- data$true_params
  truth_L_s <- as.integer(unlist(truth$L_s, use.names = FALSE))
  if (length(truth_L_s) == 1L) truth_L_s <- rep(truth_L_s, truth$S)
  if (truth$S != 2L || truth$L_f != 1L ||
      !identical(truth_L_s, c(1L, 1L))) {
    stop("Candidate2 pooled oracle requires S=2 and L=(1,1,1).")
  }
  vectors <- list(
    A0 = as.numeric(truth$a_true[, 1L]),
    B10 = as.numeric(truth$b_true[[1L]][, 1L]),
    B20 = as.numeric(truth$b_true[[2L]][, 1L])
  )
  table <- data.frame(
    oracle_id = names(vectors), oracle_role = c("shared", "specific", "specific"),
    active_studies = c("1;2", "1", "2"), inactive_studies = c("", "2", "1"),
    stringsAsFactors = FALSE
  )
  list(table = table, vectors = vectors)
}

candidate2_match_pooled_directions <- function(fit, data) {
  if (!is.list(fit) || !is.matrix(fit$B_hat) || ncol(fit$B_hat) != 3L ||
      length(fit$factor_ppi) != 3L) {
    stop("Pooled candidate2 fit must contain exactly Q=3 loading directions.")
  }
  truth <- candidate2_pooled_truth_catalog(data)
  similarity <- matrix(NA_real_, nrow = 3L, ncol = 3L,
                       dimnames = list(truth$table$oracle_id,
                                       paste0("pooled_", seq_len(3L))))
  signed <- similarity
  for (row in seq_len(3L)) for (column in seq_len(3L)) {
    signed[row, column] <- candidate2_safe_cosine(
      truth$vectors[[row]], fit$B_hat[, column], absolute = FALSE
    )
    similarity[row, column] <- abs(signed[row, column])
  }
  permutations <- candidate2_permutations(seq_len(3L))
  objectives <- apply(permutations, 1L, function(permutation) {
    sum(similarity[cbind(seq_len(3L), permutation)])
  })
  best <- which(objectives == max(objectives))[1L]
  assignment <- as.integer(permutations[best, ])
  rows <- truth$table
  rows$pooled_factor <- assignment
  rows$loading_abs_cosine <- similarity[cbind(seq_len(3L), assignment)]
  rows$loading_sign <- ifelse(
    signed[cbind(seq_len(3L), assignment)] < 0, -1, 1
  )
  rows$factor_ppi <- as.numeric(fit$factor_ppi[assignment])
  rows$auxiliary_factor_retained <-
    rows$factor_ppi >= V0LV_CANDIDATE2_FACTOR_PPI_THRESHOLD
  rows$evaluated_in_primary <- TRUE
  rows$matching_target <- V0LV_CANDIDATE2_POOLED_MATCH_TARGET
  rows$matching_objective <- objectives[[best]]
  rows$ppi_used_for_matching <- FALSE
  rows
}

candidate2_truth_path <- function(data, study, oracle_id, subject) {
  if (oracle_id == "A0") {
    as.numeric(data$true_params$f_true_values[[study]][[subject]][, 1L])
  } else if (oracle_id == "B10" && study == 1L) {
    as.numeric(data$true_params$g_true_values[[1L]][[subject]][, 1L])
  } else if (oracle_id == "B20" && study == 2L) {
    as.numeric(data$true_params$g_true_values[[2L]][[subject]][, 1L])
  } else rep(0, length(data$time_obs[[study]][[subject]]))
}

candidate2_truth_loading <- function(data, oracle_id) {
  if (oracle_id == "A0") data$true_params$a_true[, 1L]
  else if (oracle_id == "B10") data$true_params$b_true[[1L]][, 1L]
  else if (oracle_id == "B20") data$true_params$b_true[[2L]][, 1L]
  else stop("Unknown oracle direction.")
}

candidate2_global_subject <- function(data, study, subject) {
  as.integer(sum(data$true_params$n_s[seq_len(study - 1L)]) + subject)
}

candidate2_estimated_pooled_path <- function(fit, data, study, subject, factor) {
  global_subject <- candidate2_global_subject(data, study, subject)
  value <- as.numeric(fit$list_h_hat[[global_subject]][[factor]])
  target <- data$time_obs[[study]][[subject]]
  expected <- length(target)
  if (length(value) == length(fit$time_g)) {
    value <- as.numeric(stats::approx(
      fit$time_g, value, xout = target, rule = 2, ties = mean
    )$y)
  }
  if (length(value) != expected) {
    stop("Pooled factor path is not aligned to the subject observation grid.")
  }
  value
}

candidate2_component_assignment <- function(truth_functions, estimate_functions,
                                            weights) {
  truth_functions <- as.matrix(truth_functions)
  estimate_functions <- as.matrix(estimate_functions)
  similarity <- matrix(NA_real_, ncol(truth_functions), ncol(estimate_functions))
  for (row in seq_len(nrow(similarity))) for (column in seq_len(ncol(similarity))) {
    similarity[row, column] <- candidate2_safe_cosine(
      sqrt(weights) * truth_functions[, row],
      sqrt(weights) * estimate_functions[, column]
    )
  }
  permutations <- candidate2_permutations(seq_len(ncol(estimate_functions)))
  if (ncol(truth_functions) != ncol(estimate_functions)) {
    stop("Candidate2 pooled oracle fixes both true and fitted FPCA M at two.")
  }
  objectives <- apply(permutations, 1L, function(permutation) {
    sum(similarity[cbind(seq_len(ncol(truth_functions)), permutation)])
  })
  assignment <- as.integer(permutations[which(objectives == max(objectives))[1L], ])
  list(assignment = assignment, similarity = similarity)
}

candidate2_weighted_ise <- function(estimate, truth, weights) {
  sum(weights * (as.numeric(estimate) - as.numeric(truth))^2)
}

candidate2_projection_floor <- function(truth_function, design, weights) {
  design <- as.matrix(design)
  weighted <- sqrt(weights) * design
  coefficient <- qr.solve(weighted, sqrt(weights) * truth_function, tol = 1e-10)
  projection <- as.numeric(design %*% coefficient)
  list(projection = projection,
       ise = candidate2_weighted_ise(projection, truth_function, weights))
}

candidate2_pooled_spline_design <- function(fit, data, bayes_environment) {
  if (is.null(bayes_environment) ||
      !is.function(bayes_environment$get_grid_objects)) {
    stop("The loaded bayesSYNC reference environment is required for projection floors.")
  }
  pooled_time <- unname(unlist(data$time_obs, recursive = FALSE))
  p <- data$true_params$p
  multivariate_time <- lapply(pooled_time, function(value) rep(list(value), p))
  grid <- bayes_environment$get_grid_objects(
    multivariate_time, K = fit$K, n_g = NULL, time_g = fit$time_g
  )
  as.matrix(grid$C_g[[1L]])
}

candidate2_pooled_function_metrics <- function(
    fit, data, matches, bayes_environment) {
  design <- candidate2_pooled_spline_design(fit, data, bayes_environment)
  grid <- as.numeric(fit$time_g)
  weights <- candidate2_trapezoid(grid)
  rows <- list(); index <- 0L
  for (row in seq_len(nrow(matches))) {
    oracle_id <- matches$oracle_id[[row]]
    factor <- matches$pooled_factor[[row]]
    truth_dense <- if (oracle_id == "A0") {
      data$true_params$phi_dense_list[[1L]]
    } else if (oracle_id == "B10") {
      data$true_params$psi_dense_list[[1L]][[1L]]
    } else data$true_params$psi_dense_list[[2L]][[1L]]
    truth <- candidate2_interpolate_columns(
      data$true_params$t_grid_dense, truth_dense, grid
    )
    estimate <- as.matrix(fit$list_list_Phi_hat[[factor]])
    component_match <- candidate2_component_assignment(truth, estimate, weights)
    active <- as.integer(strsplit(matches$active_studies[[row]], ";", fixed = TRUE)[[1L]])
    truth_scores <- if (oracle_id == "B10") {
      data$true_params$xi_fpca_true[[1L]][[1L]]
    } else if (oracle_id == "B20") {
      data$true_params$xi_fpca_true[[2L]][[1L]]
    } else NULL
    if (oracle_id == "A0") {
      truth_scores <- do.call(rbind, lapply(active, function(study) {
        data$true_params$zeta_fpca_true[[study]][[1L]]
      }))
    }
    active_global <- unlist(lapply(active, function(study) {
      offset <- sum(data$true_params$n_s[seq_len(study - 1L)])
      offset + seq_len(data$true_params$n_s[[study]])
    }), use.names = FALSE)
    estimate_scores <- as.matrix(fit$list_Zeta_hat[[factor]])[active_global, , drop = FALSE]
    for (component in seq_len(ncol(truth))) {
      estimated_component <- component_match$assignment[[component]]
      raw_inner <- sum(weights * truth[, component] * estimate[, estimated_component])
      function_sign <- if (!is.finite(raw_inner) || raw_inner >= 0) 1 else -1
      aligned_function <- function_sign * estimate[, estimated_component]
      aligned_score <- matches$loading_sign[[row]] * function_sign *
        estimate_scores[, estimated_component]
      floor <- candidate2_projection_floor(truth[, component], design, weights)
      total_ise <- candidate2_weighted_ise(aligned_function,
                                           truth[, component], weights)
      index <- index + 1L
      rows[[index]] <- data.frame(
        oracle_id = oracle_id, pooled_factor = factor,
        true_component = component, estimated_component = estimated_component,
        function_abs_inner_product =
          component_match$similarity[component, estimated_component],
        feature_total_ise = total_ise,
        feature_projection_floor_ise = floor$ise,
        feature_excess_ise = max(0, total_ise - floor$ise),
        fpca_score_signed_cor = candidate2_safe_cor(
          aligned_score, truth_scores[, component]
        ),
        fpca_score_abs_cor = candidate2_safe_cor(
          aligned_score, truth_scores[, component], absolute = TRUE
        ),
        evaluated_studies = matches$active_studies[[row]],
        inactive_studies_excluded_from_score_metric = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

candidate2_pooled_direction_metrics <- function(fit, data, matches) {
  rows <- list()
  for (row in seq_len(nrow(matches))) {
    oracle_id <- matches$oracle_id[[row]]
    factor <- matches$pooled_factor[[row]]
    sign <- matches$loading_sign[[row]]
    active <- as.integer(strsplit(matches$active_studies[[row]], ";", fixed = TRUE)[[1L]])
    inactive_text <- matches$inactive_studies[[row]]
    inactive <- if (nzchar(inactive_text)) as.integer(inactive_text) else integer()
    estimated_active_path <- true_active_path <- numeric()
    estimated_active_contribution <- true_active_contribution <- numeric()
    inactive_contribution <- numeric()
    truth_loading <- candidate2_truth_loading(data, oracle_id)
    for (study in active) for (subject in seq_len(data$true_params$n_s[[study]])) {
      estimated_path <- candidate2_estimated_pooled_path(
        fit, data, study, subject, factor
      )
      true_path <- candidate2_truth_path(data, study, oracle_id, subject)
      estimated_active_path <- c(estimated_active_path, sign * estimated_path)
      true_active_path <- c(true_active_path, true_path)
      estimated_active_contribution <- c(
        estimated_active_contribution,
        as.numeric(outer(estimated_path, fit$B_hat[, factor]))
      )
      true_active_contribution <- c(
        true_active_contribution,
        as.numeric(outer(true_path, truth_loading))
      )
    }
    if (length(inactive)) for (study in inactive) {
      for (subject in seq_len(data$true_params$n_s[[study]])) {
        estimated_path <- candidate2_estimated_pooled_path(
          fit, data, study, subject, factor
        )
        inactive_contribution <- c(
          inactive_contribution,
          as.numeric(outer(estimated_path, fit$B_hat[, factor]))
        )
      }
    }
    active_energy <- sum(estimated_active_contribution^2)
    inactive_energy <- sum(inactive_contribution^2)
    leakage_fraction <- if (length(inactive_contribution) &&
                            active_energy + inactive_energy > 0) {
      inactive_energy / (active_energy + inactive_energy)
    } else NA_real_
    rows[[row]] <- data.frame(
      oracle_id = oracle_id, pooled_factor = factor,
      loading_abs_cosine = matches$loading_abs_cosine[[row]],
      loading_relative_l2_error = candidate2_nrmse(
        sign * fit$B_hat[, factor], truth_loading
      ),
      factor_process_signed_cor = candidate2_safe_cor(
        estimated_active_path, true_active_path
      ),
      factor_process_abs_cor = candidate2_safe_cor(
        estimated_active_path, true_active_path, absolute = TRUE
      ),
      factor_process_nrmse = candidate2_nrmse(
        estimated_active_path, true_active_path
      ),
      complete_contribution_nrmse = candidate2_nrmse(
        estimated_active_contribution, true_active_contribution
      ),
      active_studies = matches$active_studies[[row]],
      inactive_studies = inactive_text,
      inactive_study_leakage_rms = if (length(inactive_contribution)) {
        sqrt(mean(inactive_contribution^2))
      } else NA_real_,
      inactive_study_leakage_energy_fraction = leakage_fraction,
      inactive_leakage_not_applicable_reason = if (!length(inactive)) {
        "shared_direction_is_active_in_both_studies"
      } else NA_character_,
      misplaced_status = NA_character_,
      misplaced_not_applicable_reason =
        "pooled_bayesSYNC_has_no_shared_specific_parameter_blocks",
      block_purity = NA_real_,
      block_purity_not_applicable_reason =
        "pooled_bayesSYNC_has_no_shared_specific_parameter_blocks",
      factor_ppi = matches$factor_ppi[[row]],
      auxiliary_factor_retained = matches$auxiliary_factor_retained[[row]],
      evaluated_in_primary = TRUE,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

candidate2_pooled_reconstruction_metrics <- function(fit, data) {
  estimated_observed <- truth_observed <- numeric()
  estimated_dense <- truth_dense <- numeric()
  mean_amplitude <- 0.6
  global_subject <- 0L
  for (study in seq_len(data$true_params$S)) {
    for (subject in seq_len(data$true_params$n_s[[study]])) {
      global_subject <- global_subject + 1L
      estimated_dense_subject_reported <- do.call(
        cbind, fit$list_Y_hat[[global_subject]]
      )
      target <- data$time_obs[[study]][[subject]]
      estimated_subject <- candidate2_interpolate_columns(
        fit$time_g, estimated_dense_subject_reported, target
      )
      truth_subject <- data$true_params$signal_true_values[[study]][[subject]]
      estimated_observed <- c(estimated_observed, as.numeric(estimated_subject))
      truth_observed <- c(truth_observed, as.numeric(truth_subject))
      grid <- fit$time_g
      mean_dense <- vapply(seq_len(data$true_params$p), function(variable) {
        mean_amplitude * (
          sin(2 * pi * grid + data$true_params$phase_mu_mat[study, variable]) +
            0.5 * sin(4 * pi * grid +
                        2 * data$true_params$phase_mu_mat[study, variable])
        )
      }, numeric(length(grid)))
      shared_path <- candidate2_interpolate_columns(
        data$true_params$t_grid_dense,
        data$true_params$phi_dense_list[[1L]], grid
      ) %*% data$true_params$zeta_fpca_true[[study]][[1L]][subject, ]
      specific_path <- candidate2_interpolate_columns(
        data$true_params$t_grid_dense,
        data$true_params$psi_dense_list[[study]][[1L]], grid
      ) %*% data$true_params$xi_fpca_true[[study]][[1L]][subject, ]
      true_dense_subject <- mean_dense +
        outer(as.numeric(shared_path), data$true_params$a_true[, 1L]) +
        outer(as.numeric(specific_path),
              data$true_params$b_true[[study]][, 1L])
      estimated_mean <- do.call(cbind, fit$list_mu_hat)
      estimated_dense_subject <- estimated_mean
      for (factor in seq_len(fit$Q)) {
        path <- fit$list_list_Phi_hat[[factor]] %*%
          fit$list_Zeta_hat[[factor]][global_subject, ]
        estimated_dense_subject <- estimated_dense_subject +
          outer(as.numeric(path), fit$B_hat[, factor])
      }
      estimated_dense <- c(estimated_dense, as.numeric(estimated_dense_subject))
      truth_dense <- c(truth_dense, as.numeric(true_dense_subject))
    }
  }
  data.frame(
    observed_signal_nrmse = candidate2_nrmse(
      estimated_observed, truth_observed
    ),
    dense_grid_signal_nrmse = candidate2_nrmse(
      estimated_dense, truth_dense
    ),
    observed_points = length(truth_observed),
    dense_grid_points = length(truth_dense),
    stringsAsFactors = FALSE
  )
}

candidate2_evaluate_pooled <- function(
    fit, unsealed_truth, project_root, bayes_environment = NULL,
    mode = c("formal", "smoke"), binding_dir = NULL) {
  mode <- match.arg(mode)
  if (mode == "formal" && !is.null(bayes_environment)) {
    stop("Formal pooled evaluation forbids reference-environment injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (is.null(bayes_environment)) {
    bayes_environment <- candidate2_load_bayes_reference(project_root)
  }
  data <- candidate2_truth_data(unsealed_truth)
  matches <- candidate2_match_pooled_directions(fit, data)
  direction <- candidate2_pooled_direction_metrics(fit, data, matches)
  functions <- candidate2_pooled_function_metrics(
    fit, data, matches, bayes_environment
  )
  reconstruction <- candidate2_pooled_reconstruction_metrics(fit, data)
  list(
    direction_matches = matches,
    direction_metrics = direction,
    feature_and_score_metrics = functions,
    reconstruction_metrics = reconstruction,
    contract = data.frame(
      evaluator_version = V0LV_CANDIDATE2_EVALUATOR_VERSION,
      studies_combined_and_fit_once = TRUE,
      study_labels_passed_to_model = FALSE,
      Q = 3L, M = 2L,
      all_three_oracle_candidates_evaluated = TRUE,
      factor_ppi_rule = ">=0.5_auxiliary_only",
      cross_model_elbo_comparison = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

candidate2_cluster_bootstrap <- function(
    values, clusters, B = V0LV_CANDIDATE2_BOOTSTRAP_REPLICATES,
    seed = V0LV_CANDIDATE2_BOOTSTRAP_SEED,
    statistic = stats::median) {
  values <- as.numeric(values); clusters <- as.character(clusters)
  B <- c2_integer(B, "B", 1L); seed <- c2_integer(seed, "seed")
  if (length(values) != length(clusters) || !length(values) ||
      any(!is.finite(values)) || any(!nzchar(clusters))) {
    stop("Cluster bootstrap inputs are malformed.")
  }
  ids <- sort(unique(clusters))
  set.seed(seed)
  replicates <- replicate(B, {
    sampled <- sample(ids, length(ids), replace = TRUE)
    sampled_values <- unlist(lapply(sampled, function(id) values[clusters == id]),
                             use.names = FALSE)
    statistic(sampled_values)
  })
  c(lower = unname(stats::quantile(replicates, 0.025, type = 8)),
    estimate = statistic(values),
    upper = unname(stats::quantile(replicates, 0.975, type = 8)))
}

candidate2_exact_sign_test <- function(difference) {
  difference <- as.numeric(difference)
  difference <- difference[is.finite(difference)]
  nonzero <- difference[difference != 0]
  if (!length(nonzero)) {
    return(data.frame(nonzero_n = 0L, improved_n = 0L, worsened_n = 0L,
                      tied_n = length(difference), two_sided_p = 1,
                      stringsAsFactors = FALSE))
  }
  test <- stats::binom.test(sum(nonzero < 0), length(nonzero), p = 0.5,
                            alternative = "two.sided")
  data.frame(
    nonzero_n = length(nonzero), improved_n = sum(nonzero < 0),
    worsened_n = sum(nonzero > 0), tied_n = sum(difference == 0),
    two_sided_p = unname(test$p.value), stringsAsFactors = FALSE
  )
}

candidate2_summarize_paired_differences <- function(paired_metrics) {
  required <- c("data_id", "comparison", "metric", "difference")
  if (!is.data.frame(paired_metrics) ||
      !all(required %in% names(paired_metrics))) {
    stop("paired_metrics lacks candidate2 paired-difference columns.")
  }
  groups <- split(paired_metrics,
                  interaction(paired_metrics$comparison, paired_metrics$metric,
                              drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(groups, function(group) {
    difference <- as.numeric(group$difference)
    if (any(!is.finite(difference)) || anyDuplicated(group$data_id)) {
      stop("Each paired metric requires one finite difference per data cluster.")
    }
    bootstrap <- candidate2_cluster_bootstrap(
      difference, group$data_id,
      B = V0LV_CANDIDATE2_BOOTSTRAP_REPLICATES,
      seed = V0LV_CANDIDATE2_BOOTSTRAP_SEED
    )
    sign <- candidate2_exact_sign_test(difference)
    data.frame(
      comparison = group$comparison[[1L]], metric = group$metric[[1L]],
      data_n = length(difference), median_difference = stats::median(difference),
      q1_difference = unname(stats::quantile(difference, 0.25, type = 8)),
      q3_difference = unname(stats::quantile(difference, 0.75, type = 8)),
      cluster_bootstrap_lower = bootstrap[["lower"]],
      cluster_bootstrap_upper = bootstrap[["upper"]],
      exact_sign_two_sided_p = sign$two_sided_p,
      improved_count = sign$improved_n, worsened_count = sign$worsened_n,
      tied_count = sign$tied_n,
      median_no_worsening_or_advantage_gate = stats::median(difference) <= 0,
      worsening_count_is_descriptive_only = TRUE,
      invented_noninferiority_margin_used = FALSE,
      success_gate_name = "中位无恶化/优势门槛",
      bootstrap_replicates = V0LV_CANDIDATE2_BOOTSTRAP_REPLICATES,
      bootstrap_seed = V0LV_CANDIDATE2_BOOTSTRAP_SEED,
      stringsAsFactors = FALSE
    )
  }))
}
