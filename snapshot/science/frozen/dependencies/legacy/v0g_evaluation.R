v0g_safe_abs_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(0)
  value <- suppressWarnings(stats::cor(x, y))
  if (is.finite(value)) abs(value) else 0
}

v0g_safe_abs_cosine <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(0)
  min(1, abs(sum(x * y) / denominator))
}

v0g_interpolate_dense <- function(grid, value, target) {
  as.numeric(stats::approx(
    grid, as.numeric(value), xout = target,
    rule = 2, ties = mean
  )$y)
}

v0g_reported_paths_by_study <- function(fit, data) {
  n_study <- as.integer(data$true_params$n_s)
  shared <- vector("list", length(n_study))
  specific <- vector("list", length(n_study))
  subject_offset <- 0L
  for (study in seq_along(n_study)) {
    L_ss <- length(fit$list_Phi_hat_spec[[study]])
    n_rows <- sum(lengths(data$time_obs[[study]]))
    shared[[study]] <- matrix(0, nrow = n_rows, ncol = fit$L_f)
    specific[[study]] <- matrix(0, nrow = n_rows, ncol = L_ss)
    row_offset <- 0L
    for (subject in seq_len(n_study[[study]])) {
      target <- data$time_obs[[study]][[subject]]
      rows <- row_offset + seq_along(target)
      global_subject <- subject_offset + subject
      for (factor in seq_len(fit$L_f)) {
        dense <- fit$list_Phi_hat[[factor]] %*%
          fit$list_Zeta_hat[[factor]][global_subject, ]
        shared[[study]][rows, factor] <- v0g_interpolate_dense(
          fit$time_g, dense, target
        )
      }
      for (factor in seq_len(L_ss)) {
        dense <- fit$list_Phi_hat_spec[[study]][[factor]] %*%
          fit$list_Zeta_hat_spec[[study]][[factor]][subject, ]
        specific[[study]][rows, factor] <- v0g_interpolate_dense(
          fit$time_g, dense, target
        )
      }
      row_offset <- row_offset + length(target)
    }
    subject_offset <- subject_offset + n_study[[study]]
  }
  list(shared = shared, specific = specific)
}

v0g_truth_role_paths <- function(data) {
  list(
    shared = lapply(seq_len(2L), function(study) {
      do.call(rbind, lapply(
        data$true_params$f_true_values[[study]],
        function(value) as.matrix(value)[, 1L, drop = FALSE]
      ))
    }),
    specific = lapply(seq_len(2L), function(study) {
      do.call(rbind, lapply(
        data$true_params$g_true_values[[study]],
        function(value) as.matrix(value)[, 1L, drop = FALSE]
      ))
    })
  )
}

v0g_role_coverage <- function(fit, data, paths = NULL, thresholds = NULL) {
  if (is.null(paths)) paths <- v0g_reported_paths_by_study(fit, data)
  if (is.null(thresholds)) {
    thresholds <- list(
      trajectory_abs_cor = 0.8,
      loading_abs_cosine = 0.8,
      factor_ppi = 0.5
    )
  }
  truth_paths <- v0g_truth_role_paths(data)
  rows <- list()
  index <- 0L
  for (study in seq_len(2L)) {
    for (block in c("shared", "specific")) {
      estimated_paths <- paths[[block]][[study]]
      loading <- if (block == "shared") {
        fit$mu_q_a_hat
      } else {
        fit$mu_q_b_specific_hat[[study]]
      }
      ppi <- if (block == "shared") {
        fit$factor_ppi_shared
      } else {
        fit$factor_ppi_specific[[study]]
      }
      if (!ncol(estimated_paths)) next
      for (factor in seq_len(ncol(estimated_paths))) {
        for (true_role in c("shared", "specific")) {
          true_loading <- if (true_role == "shared") {
            data$true_params$a_true[, 1L]
          } else {
            data$true_params$b_true[[study]][, 1L]
          }
          path_cor <- v0g_safe_abs_cor(
            estimated_paths[, factor], truth_paths[[true_role]][[study]][, 1L]
          )
          loading_cosine <- v0g_safe_abs_cosine(
            loading[, factor], true_loading
          )
          index <- index + 1L
          rows[[index]] <- data.frame(
            study = study,
            estimated_block = block,
            candidate_factor = factor,
            factor_ppi = as.numeric(ppi[[factor]]),
            retained = as.numeric(ppi[[factor]]) >= thresholds$factor_ppi,
            true_role = true_role,
            trajectory_abs_cor = path_cor,
            loading_abs_cosine = loading_cosine,
            joint_score = (path_cor + loading_cosine) / 2,
            strict_aligned =
              path_cor >= thresholds$trajectory_abs_cor &&
              loading_cosine >= thresholds$loading_abs_cosine,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  alignment <- if (length(rows)) do.call(rbind, rows) else data.frame()
  coverage_rows <- list()
  index <- 0L
  for (study in seq_len(2L)) {
    for (true_role in c("shared", "specific")) {
      available <- alignment[
        alignment$study == study & alignment$true_role == true_role &
          alignment$retained,
        , drop = FALSE
      ]
      index <- index + 1L
      if (!nrow(available)) {
        coverage_rows[[index]] <- data.frame(
          study = study,
          true_role = true_role,
          best_estimated_block = NA_character_,
          best_candidate_factor = NA_integer_,
          factor_ppi = 0,
          trajectory_abs_cor = 0,
          loading_abs_cosine = 0,
          joint_score = 0,
          strict_aligned = FALSE,
          correct_type = FALSE,
          misplaced = FALSE,
          missing = TRUE,
          stringsAsFactors = FALSE
        )
      } else {
        eligible <- if (any(available$strict_aligned)) {
          available[available$strict_aligned, , drop = FALSE]
        } else {
          available
        }
        best <- eligible[which.max(eligible$joint_score), , drop = FALSE]
        aligned <- isTRUE(best$strict_aligned[[1L]])
        coverage_rows[[index]] <- data.frame(
          study = study,
          true_role = true_role,
          best_estimated_block = best$estimated_block,
          best_candidate_factor = best$candidate_factor,
          factor_ppi = best$factor_ppi,
          trajectory_abs_cor = best$trajectory_abs_cor,
          loading_abs_cosine = best$loading_abs_cosine,
          joint_score = best$joint_score,
          strict_aligned = aligned,
          correct_type = aligned && best$estimated_block == true_role,
          misplaced = aligned && best$estimated_block != true_role,
          missing = !aligned,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(alignment = alignment, coverage = do.call(rbind, coverage_rows))
}

v0g_prediction_diagnostics <- function(fit, data, paths) {
  flatten <- function(value) unlist(value, recursive = TRUE, use.names = FALSE)
  rms <- function(value) sqrt(mean(as.numeric(value)^2))
  signal_nrmse <- function(estimate, truth) {
    estimate <- as.numeric(estimate)
    truth <- as.numeric(truth)
    sqrt(mean((estimate - truth)^2)) / max(rms(truth), 1e-12)
  }
  shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
  specific_mask <- lapply(
    fit$factor_ppi_specific,
    function(value) as.numeric(value >= 0.5)
  )
  all_prediction <- vector("list", fit$S)
  selected_prediction <- vector("list", fit$S)
  mean_prediction <- vector("list", fit$S)
  for (study in seq_len(fit$S)) {
    all_prediction[[study]] <- vector("list", fit$n_s[study])
    selected_prediction[[study]] <- vector("list", fit$n_s[study])
    mean_prediction[[study]] <- vector("list", fit$n_s[study])
    local_offset <- 0L
    for (subject in seq_len(fit$n_s[study])) {
      target_time <- data$time_obs[[study]][[subject]]
      n_time <- length(target_time)
      rows <- local_offset + seq_len(n_time)
      mean_matrix <- vapply(seq_len(fit$p), function(variable) {
        v0g_interpolate_dense(
          fit$time_g, fit$list_mu_hat[[study]][[variable]], target_time
        )
      }, numeric(n_time))
      shared_all <- paths$shared[[study]][rows, , drop = FALSE] %*%
        t(fit$mu_q_a_hat)
      specific_all <- paths$specific[[study]][rows, , drop = FALSE] %*%
        t(fit$mu_q_b_specific_hat[[study]])
      shared_selected <- paths$shared[[study]][rows, , drop = FALSE] %*%
        t(sweep(fit$mu_q_a_hat, 2L, shared_mask, "*"))
      specific_selected <- paths$specific[[study]][rows, , drop = FALSE] %*%
        t(sweep(
          fit$mu_q_b_specific_hat[[study]], 2L,
          specific_mask[[study]], "*"
        ))
      mean_prediction[[study]][[subject]] <- mean_matrix
      all_prediction[[study]][[subject]] <-
        mean_matrix + shared_all + specific_all
      selected_prediction[[study]][[subject]] <-
        mean_matrix + shared_selected + specific_selected
      local_offset <- local_offset + n_time
    }
  }
  truth_signal <- flatten(data$true_params$signal_true_values)
  c(
    signal_nrmse_all_candidates =
      signal_nrmse(flatten(all_prediction), truth_signal),
    signal_nrmse_ppi_selected =
      signal_nrmse(flatten(selected_prediction), truth_signal),
    signal_nrmse_mean_only =
      signal_nrmse(flatten(mean_prediction), truth_signal)
  )
}

v0g_ppi_table <- function(variable_ppi, factor_ppi, block, study) {
  variable_ppi <- as.matrix(variable_ppi)
  union_ppi <- 1 - apply(1 - pmin(pmax(variable_ppi, 0), 1), 2L, prod)
  data.frame(
    block = block,
    study = study,
    candidate_factor = seq_len(ncol(variable_ppi)),
    factor_ppi = as.numeric(factor_ppi),
    recomputed_union_ppi = union_ppi,
    union_abs_error = abs(as.numeric(factor_ppi) - union_ppi),
    expected_active_loading_count = colSums(variable_ppi),
    retained = as.numeric(factor_ppi) >= 0.5,
    stringsAsFactors = FALSE
  )
}

v0g_shared_candidate_metrics <- function(
    fit, paths, cross_study_threshold = 0.10) {
  rows <- lapply(seq_len(fit$L_f), function(factor) {
    loading <- fit$mu_q_a_hat[, factor]
    contribution <- vapply(seq_len(2L), function(study) {
      sqrt(mean(paths$shared[[study]][, factor]^2)) *
        sqrt(mean(loading^2))
    }, numeric(1))
    pooled_path <- unlist(lapply(
      paths$shared, function(value) value[, factor]
    ), use.names = FALSE)
    pooled <- sqrt(mean(pooled_path^2)) * sqrt(mean(loading^2))
    data.frame(
      candidate_factor = factor,
      factor_ppi = as.numeric(fit$factor_ppi_shared[[factor]]),
      retained = as.numeric(fit$factor_ppi_shared[[factor]]) >= 0.5,
      loading_norm = sqrt(sum(loading^2)),
      contribution_rms_pooled = pooled,
      contribution_rms_study_1 = contribution[[1L]],
      contribution_rms_study_2 = contribution[[2L]],
      study_min_max_ratio = if (max(contribution) > 0) {
        min(contribution) / max(contribution)
      } else 0,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  retained_pool_max <- max(result$contribution_rms_pooled[result$retained], 0)
  result$standardized_min_contribution <- if (retained_pool_max > 0) {
    pmin(
      result$contribution_rms_study_1,
      result$contribution_rms_study_2
    ) / retained_pool_max
  } else 0
  result$cross_study_supported_0_1 <-
    result$retained &
    result$standardized_min_contribution >= cross_study_threshold
  result
}

v0g_evaluate_fit <- function(fit, data, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- list(
      trajectory_abs_cor = 0.8,
      loading_abs_cosine = 0.8,
      factor_ppi = 0.5,
      cross_study_standardized_min_contribution = 0.10
    )
  }
  paths <- v0g_reported_paths_by_study(fit, data)
  role <- v0g_role_coverage(fit, data, paths, thresholds)
  prediction <- v0g_prediction_diagnostics(fit, data, paths)
  selected_shared <- sum(fit$factor_ppi_shared >= thresholds$factor_ppi)
  selected_specific <- vapply(
    fit$factor_ppi_specific,
    function(value) sum(value >= thresholds$factor_ppi), integer(1)
  )
  ppi <- do.call(rbind, c(
    list(v0g_ppi_table(
      fit$mu_q_gamma_a_hat, fit$factor_ppi_shared, "shared", 0L
    )),
    lapply(seq_len(2L), function(study) {
      v0g_ppi_table(
        fit$mu_q_gamma_b_hat[[study]],
        fit$factor_ppi_specific[[study]], "specific", study
      )
    })
  ))
  shared_metrics <- v0g_shared_candidate_metrics(
    fit, paths,
    thresholds$cross_study_standardized_min_contribution
  )
  list(
    paths = paths,
    role = role,
    prediction = prediction,
    selected_shared = selected_shared,
    selected_specific = selected_specific,
    count_vector = c(selected_shared, selected_specific),
    ppi = ppi,
    shared_metrics = shared_metrics
  )
}
