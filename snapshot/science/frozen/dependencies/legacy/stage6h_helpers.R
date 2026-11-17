s6h_finite_mean <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}

# This is the Stage 6G additional-metric calculation carried forward without
# changing its scientific definitions.
s6h_multi_additional_metrics <- function(
    fit, truth, evaluation, evaluator_environment) {
  data <- truth$data
  reported_paths <- evaluator_environment$v0g_reported_paths_by_study(
    fit, data
  )
  true_paths <- evaluator_environment$v0h_truth_paths_by_study(data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected", , drop = FALSE]
  expected_roles <- data$true_params$S * data$true_params$L_f +
    sum(as.integer(unlist(data$true_params$L_s_by_study,
                          use.names = FALSE)))
  pp_assert(nrow(functional) == expected_roles,
            "Unexpected Stage 6H study-role truth-row count.")
  process_rows <- lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_role <- as.character(item$true_role[[1L]])
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_role == "shared") {
      data$true_params$a_true[, true_factor]
    } else {
      data$true_params$b_true[[study]][, true_factor]
    }
    true_path <- true_paths[[true_role]][[study]][, true_factor]
    true_contribution <- as.numeric(outer(true_path, true_loading))
    matched <- isTRUE(item$matched_by_loading[[1L]])
    if (matched) {
      estimated_role <- as.character(item$estimated_block[[1L]])
      estimated_factor <- as.integer(item$estimated_factor[[1L]])
      estimated_loading <- if (estimated_role == "shared") {
        fit$mu_q_a_hat[, estimated_factor]
      } else {
        fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      }
      estimated_path <- reported_paths[[estimated_role]][[study]][,
                                                                 estimated_factor]
      aligned_path <- as.numeric(item$loading_sign[[1L]]) * estimated_path
      process_error <- pp_nrmse(aligned_path, true_path)
      contribution_error <- pp_nrmse(
        as.numeric(outer(estimated_path, estimated_loading)),
        true_contribution
      )
      loading_norm_ratio <- pp_vector_l2(estimated_loading) /
        pp_vector_l2(true_loading)
      loading_cosine <- pp_abs_cosine(estimated_loading, true_loading)
      estimated_path_rms <- pp_vector_rms(estimated_path)
    } else {
      estimated_role <- NA_character_
      estimated_factor <- NA_integer_
      process_error <- contribution_error <- loading_norm_ratio <-
        loading_cosine <- estimated_path_rms <- NA_real_
    }
    data.frame(
      study = study, true_role = true_role, true_factor = true_factor,
      true_id = item$true_id[[1L]], estimated_block = estimated_role,
      estimated_factor = estimated_factor,
      loading_structure_status = item$loading_structure_status[[1L]],
      matched_by_loading = matched,
      loading_norm_ratio = loading_norm_ratio,
      loading_abs_cosine = loading_cosine,
      true_path_rms = pp_vector_rms(true_path),
      estimated_path_rms = estimated_path_rms,
      factor_process_nrmse_matched = process_error,
      factor_process_nrmse_missing_as_zero = if (matched) process_error else 1,
      complete_contribution_nrmse_matched = contribution_error,
      complete_contribution_nrmse_missing_as_zero =
        if (matched) contribution_error else 1,
      stringsAsFactors = FALSE
    )
  })
  process <- do.call(rbind, process_rows)
  process_summary <- data.frame(
    study_role_matched = sum(process$matched_by_loading),
    study_role_truth_total = nrow(process),
    factor_process_nrmse_mean_matched =
      s6h_finite_mean(process$factor_process_nrmse_matched),
    factor_process_nrmse_mean_missing_as_zero =
      mean(process$factor_process_nrmse_missing_as_zero),
    complete_contribution_nrmse_mean_matched =
      s6h_finite_mean(process$complete_contribution_nrmse_matched),
    complete_contribution_nrmse_mean_missing_as_zero =
      mean(process$complete_contribution_nrmse_missing_as_zero),
    stringsAsFactors = FALSE
  )

  loading_scale_rows <- lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_role <- as.character(item$true_role[[1L]])
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_role == "shared") {
      data$true_params$a_true[, true_factor]
    } else {
      data$true_params$b_true[[study]][, true_factor]
    }
    matched <- isTRUE(item$matched_by_loading[[1L]])
    estimated_block <- if (matched) {
      as.character(item$estimated_block[[1L]])
    } else NA_character_
    estimated_factor <- if (matched) {
      as.integer(item$estimated_factor[[1L]])
    } else NA_integer_
    if (!matched) {
      raw <- canonical <- rep(NA_real_, length(true_loading))
      factor_scale <- factor_ppi <- NA_real_
    } else if (estimated_block == "shared") {
      original_factor <- as.integer(fit$factor_order_shared[[estimated_factor]])
      raw <- fit$mu_q_a_original[, original_factor]
      canonical <- fit$mu_q_a_hat[, estimated_factor]
      factor_scale <- fit$factor_scale_shared[[estimated_factor]]
      factor_ppi <- fit$factor_ppi_shared[[estimated_factor]]
    } else {
      original_factor <- as.integer(
        fit$factor_order_specific[[study]][[estimated_factor]]
      )
      raw <- fit$mu_q_b_specific_original[[study]][, original_factor]
      canonical <- fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      factor_scale <- fit$factor_scale_specific[[study]][[estimated_factor]]
      factor_ppi <- fit$factor_ppi_specific[[study]][[estimated_factor]]
    }
    data.frame(
      study = study, true_role = true_role, true_factor = true_factor,
      estimated_block = estimated_block,
      estimated_factor = estimated_factor,
      matched_by_loading = matched, factor_ppi = factor_ppi,
      factor_scale = factor_scale,
      raw_loading_norm_ratio = pp_vector_l2(raw) / pp_vector_l2(true_loading),
      canonical_loading_norm_ratio =
        pp_vector_l2(canonical) / pp_vector_l2(true_loading),
      raw_loading_abs_cosine = pp_abs_cosine(raw, true_loading),
      canonical_loading_abs_cosine = pp_abs_cosine(canonical, true_loading),
      stringsAsFactors = FALSE
    )
  })
  loading_scale <- do.call(rbind, loading_scale_rows)
  loading_scale_summary <- data.frame(
    raw_loading_norm_ratio_mean =
      s6h_finite_mean(loading_scale$raw_loading_norm_ratio),
    canonical_loading_norm_ratio_mean =
      s6h_finite_mean(loading_scale$canonical_loading_norm_ratio),
    raw_loading_abs_cosine_mean =
      s6h_finite_mean(loading_scale$raw_loading_abs_cosine),
    canonical_loading_abs_cosine_mean =
      s6h_finite_mean(loading_scale$canonical_loading_abs_cosine),
    factor_scale_mean = s6h_finite_mean(loading_scale$factor_scale),
    stringsAsFactors = FALSE
  )
  list(
    process_contribution_by_role = process,
    process_contribution_summary = process_summary,
    loading_scale_by_block = loading_scale,
    loading_scale_summary = loading_scale_summary
  )
}

s6h_evaluate_fit <- function(fit, truth, evaluator_environment) {
  adapter <- truth
  adapter$bundle_class <- "v0lv_sealed_truth_candidate2"
  result <- candidate2_evaluate_multi(
    fit, adapter, PP_SNAPSHOT,
    evaluator_environment = evaluator_environment, mode = "smoke"
  )
  result$evaluation_contract$stage6h_original_truth_bundle_class <-
    truth$bundle_class
  result$evaluation_contract$post_unseal_oracle_diagnostic <- TRUE
  result$evaluation_contract$truth_used_for_fit_stopping_or_elbo_selection <-
    FALSE
  result$evaluation_contract$all_endpoints_frozen_before_truth_evaluation <-
    TRUE
  result$evaluation_contract$cross_dimension_selection <- FALSE
  result$evaluation_contract$continuation <- FALSE
  additional <- s6h_multi_additional_metrics(
    fit, truth, result, evaluator_environment
  )
  dimension_activity <- pp_dimension_activity_metrics(
    fit, adapter, result, evaluator_environment
  )
  c(result, additional, dimension_activity)
}
