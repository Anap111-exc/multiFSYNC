#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else "."
evaluation_root <- file.path(root, "truth_evaluation_20260912_v1")
selection_root <- file.path(root, "truth_free_selection")
analysis_root <- file.path(root, "stage6g_scientific_analysis_20260912_v1")
if (dir.exists(analysis_root)) {
  stop("Stage 6G analysis directory already exists; refusing overwrite.",
       call. = FALSE)
}
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(analysis_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
safe_mean <- function(value) {
  value <- as.numeric(value); value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
safe_median <- function(value) {
  value <- as.numeric(value); value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}
safe_min <- function(value) {
  value <- as.numeric(value); value <- value[is.finite(value)]
  if (length(value)) min(value) else NA_real_
}
safe_max <- function(value) {
  value <- as.numeric(value); value <- value[is.finite(value)]
  if (length(value)) max(value) else NA_real_
}
safe_cor <- function(x, y, method = "spearman") {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  stats::cor(x[keep], y[keep], method = method)
}
parse_rank <- function(value) {
  if (is.na(value) || !nzchar(as.character(value))) return(integer())
  as.integer(strsplit(as.character(value), ";", fixed = TRUE)[[1L]])
}

all <- read_csv(file.path(evaluation_root,
                          "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv"))
winners <- read_csv(file.path(
  evaluation_root, "CONFIG_SPLIT_WINNER_SCIENTIFIC_RESULTS.csv"
))
one_se_endpoints <- read_csv(file.path(
  evaluation_root, "ONE_SE_SELECTED_ENDPOINT_SCIENTIFIC_RESULTS.csv"
))
minimum_endpoints <- read_csv(file.path(
  evaluation_root, "MINIMUM_LOSS_ENDPOINT_SCIENTIFIC_RESULTS.csv"
))
selections <- read_csv(file.path(
  selection_root, "FINAL_TRUTH_FREE_DIMENSION_SELECTIONS.csv"
))
status <- read_csv(file.path(evaluation_root, "EVALUATION_STATUS.csv"))
metric_qc <- read_csv(file.path(evaluation_root,
                                "METRIC_COMPLETENESS_QC.csv"))
runtime <- read_csv(file.path(root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))

stopifnot(
  nrow(all) == 120L, nrow(winners) == 30L,
  nrow(one_se_endpoints) == 6L, nrow(minimum_endpoints) == 6L,
  nrow(selections) == 3L, nrow(status) == 120L, all(status$ok),
  nrow(metric_qc) == 120L, all(metric_qc$metric_contract_passed),
  length(unique(all$data_id)) == 3L,
  all(table(all$data_id, all$split_id, all$fit_config_id) == 4L),
  !any(all$continuation_used)
)

all$representationally_adequate <- vapply(seq_len(nrow(all)), function(index) {
  row <- all[index, , drop = FALSE]
  factor_capacity <- row$fit_L_f >= row$truth_L_f &&
    row$fit_L_s_1 >= row$truth_L_s_1 &&
    row$fit_L_s_2 >= row$truth_L_s_2
  ranks <- c(parse_rank(row$fit_M_f), parse_rank(row$fit_M_s_1),
             parse_rank(row$fit_M_s_2))
  factor_capacity && length(ranks) > 0L &&
    all(ranks >= row$truth_M_common)
}, logical(1L))
all$ppi_factor_counts_exact <- with(
  all,
  selected_shared == truth_L_f &
    selected_specific_study_1 == truth_L_s_1 &
    selected_specific_study_2 == truth_L_s_2
)
all$core_truth_recovered <- with(
  all,
  global_correct == truth_global_direction_count &
    global_misplaced == 0L & global_missing == 0L &
    study_correct == truth_study_role_count &
    study_misplaced == 0L & study_missing == 0L
)
all$exact_thresholded_factor_structure <- with(
  all,
  core_truth_recovered & ppi_factor_counts_exact &
    global_extra == 0L & global_duplicate == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L
)
all$full_true_feature_component_coverage <- with(
  all, feature_components_matched == feature_components_total
)
all$core_with_full_feature_coverage <- with(
  all, core_truth_recovered & full_true_feature_component_coverage
)
all$generating_effective_M_exact <- vapply(
  seq_len(nrow(all)), function(index) {
    row <- all[index, , drop = FALSE]
    if (!isTRUE(row$generating_dimension_configuration)) return(FALSE)
    shared <- parse_rank(row$shared_effective_M)
    specific_1 <- parse_rank(row$specific_effective_M_study_1)
    specific_2 <- parse_rank(row$specific_effective_M_study_2)
    length(shared) == row$truth_L_f &&
      length(specific_1) == row$truth_L_s_1 &&
      length(specific_2) == row$truth_L_s_2 &&
      all(c(shared, specific_1, specific_2) == row$truth_M_common)
  }, logical(1L)
)

derived_columns <- c(
  "fit_id", "data_id", "split_id", "fit_config_id",
  "expected_truth_fit_config_id", "representationally_adequate",
  "generating_dimension_configuration", "seed_index", "fit_seed",
  "final_elbo", "practical_converged", "convergence_status",
  "within_configuration_split_winner", "selected_by_one_se",
  "selected_by_minimum_loss", "selected_counts", "ppi_factor_counts_exact",
  "core_truth_recovered", "exact_thresholded_factor_structure",
  "full_true_feature_component_coverage",
  "core_with_full_feature_coverage", "generating_effective_M_exact",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "matched_function_recall_mean", "matched_trajectory_abs_cor_mean",
  "matched_score_abs_cor_mean", "feature_components_matched",
  "feature_components_total", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean", "loading_support_ppi_auc_mean",
  "factor_process_nrmse_mean_matched",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "global_correct", "global_misplaced", "global_missing", "global_extra",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate",
  "activity_ppi_retained_duplicate_or_extra",
  "activity_duplicate_or_extra_contribution_rms_max",
  "activity_effective_M_mean", "activity_effective_M_at_rank_cap",
  "fit_elapsed_seconds", "fit_warning_count"
)
start_level <- all[derived_columns]
start_level <- start_level[order(start_level$data_id, start_level$split_id,
                                 start_level$fit_config_id,
                                 start_level$seed_index), , drop = FALSE]

split_config_summary <- do.call(rbind, lapply(
  split(all, interaction(all$data_id, all$split_id, all$fit_config_id,
                         drop = TRUE, lex.order = TRUE)), function(rows) {
    winner <- rows[rows$within_configuration_split_winner, , drop = FALSE]
    stopifnot(nrow(rows) == 4L, nrow(winner) == 1L)
    data.frame(
      data_id = rows$data_id[[1L]], split_id = rows$split_id[[1L]],
      fit_config_id = rows$fit_config_id[[1L]],
      expected_truth_fit_config_id = rows$expected_truth_fit_config_id[[1L]],
      representationally_adequate = rows$representationally_adequate[[1L]],
      generating_dimension_configuration =
        rows$generating_dimension_configuration[[1L]],
      starts = nrow(rows), practical_converged = sum(rows$practical_converged),
      core_hits = sum(rows$core_truth_recovered),
      exact_structure_hits = sum(rows$exact_thresholded_factor_structure),
      full_feature_coverage_hits =
        sum(rows$full_true_feature_component_coverage),
      core_full_feature_hits = sum(rows$core_with_full_feature_coverage),
      any_core_hit = any(rows$core_truth_recovered),
      any_exact_structure_hit = any(rows$exact_thresholded_factor_structure),
      any_core_full_feature_hit = any(rows$core_with_full_feature_coverage),
      winner_fit_id = winner$fit_id,
      winner_seed_index = winner$seed_index,
      winner_core_hit = winner$core_truth_recovered,
      winner_exact_structure_hit = winner$exact_thresholded_factor_structure,
      winner_full_feature_coverage =
        winner$full_true_feature_component_coverage,
      winner_core_full_feature = winner$core_with_full_feature_coverage,
      elbo_missed_available_core =
        any(rows$core_truth_recovered) && !winner$core_truth_recovered,
      elbo_missed_available_exact =
        any(rows$exact_thresholded_factor_structure) &&
          !winner$exact_thresholded_factor_structure,
      winner_selected_counts = winner$selected_counts,
      winner_heldout_nrmse = winner$heldout_nrmse,
      winner_observed_signal_nrmse = winner$signal_nrmse_ppi_selected,
      winner_dense_signal_nrmse = winner$dense_signal_nrmse_ppi_selected,
      winner_feature_ise = winner$feature_total_ise_mean,
      winner_loading_relative_l2 = winner$loading_relative_l2_error_mean,
      winner_process_nrmse = winner$factor_process_nrmse_mean_matched,
      winner_contribution_nrmse =
        winner$complete_contribution_nrmse_mean_matched,
      elbo_vs_dense_spearman = safe_cor(
        rows$final_elbo, rows$dense_signal_nrmse_ppi_selected
      ),
      stringsAsFactors = FALSE
    )
  }
))
split_config_summary <- split_config_summary[order(
  split_config_summary$data_id, split_config_summary$split_id,
  split_config_summary$fit_config_id
), , drop = FALSE]

scenario_config_summary <- do.call(rbind, lapply(
  split(all, interaction(all$data_id, all$fit_config_id, drop = TRUE,
                         lex.order = TRUE)), function(rows) {
    row_winners <- rows[rows$within_configuration_split_winner, , drop = FALSE]
    data.frame(
      data_id = rows$data_id[[1L]], fit_config_id = rows$fit_config_id[[1L]],
      expected_truth_fit_config_id = rows$expected_truth_fit_config_id[[1L]],
      representationally_adequate = rows$representationally_adequate[[1L]],
      generating_dimension_configuration =
        rows$generating_dimension_configuration[[1L]],
      endpoints = nrow(rows), split_winners = nrow(row_winners),
      practical_converged = sum(rows$practical_converged),
      core_hits = sum(rows$core_truth_recovered),
      exact_structure_hits = sum(rows$exact_thresholded_factor_structure),
      full_feature_coverage_hits =
        sum(rows$full_true_feature_component_coverage),
      core_full_feature_hits = sum(rows$core_with_full_feature_coverage),
      split_winner_core_hits = sum(row_winners$core_truth_recovered),
      split_winner_exact_hits =
        sum(row_winners$exact_thresholded_factor_structure),
      split_winner_core_full_feature_hits =
        sum(row_winners$core_with_full_feature_coverage),
      ppi_factor_counts_exact = sum(rows$ppi_factor_counts_exact),
      generating_effective_M_exact = sum(rows$generating_effective_M_exact),
      observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
      dense_signal_nrmse_mean =
        safe_mean(rows$dense_signal_nrmse_ppi_selected),
      feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
      loading_relative_l2_mean =
        safe_mean(rows$loading_relative_l2_error_mean),
      process_nrmse_missing_as_zero_mean =
        safe_mean(rows$factor_process_nrmse_mean_missing_as_zero),
      contribution_nrmse_missing_as_zero_mean =
        safe_mean(rows$complete_contribution_nrmse_mean_missing_as_zero),
      extra_or_duplicate_retained_mean = safe_mean(
        rows$activity_ppi_retained_duplicate_or_extra
      ), stringsAsFactors = FALSE
    )
  }
))
scenario_config_summary <- scenario_config_summary[order(
  scenario_config_summary$data_id, scenario_config_summary$fit_config_id
), , drop = FALSE]

reachability <- do.call(rbind, lapply(split(all, all$data_id), function(rows) {
  adequate <- rows[rows$representationally_adequate, , drop = FALSE]
  generating <- rows[rows$generating_dimension_configuration, , drop = FALSE]
  generating_winners <- generating[
    generating$within_configuration_split_winner, , drop = FALSE
  ]
  strata <- split_config_summary[
    split_config_summary$data_id == rows$data_id[[1L]], , drop = FALSE
  ]
  gen_strata <- strata[strata$generating_dimension_configuration,
                       , drop = FALSE]
  data.frame(
    data_id = rows$data_id[[1L]],
    expected_truth_fit_config_id = rows$expected_truth_fit_config_id[[1L]],
    endpoints = nrow(rows), adequate_endpoints = nrow(adequate),
    adequate_core_hits = sum(adequate$core_truth_recovered),
    adequate_exact_hits = sum(adequate$exact_thresholded_factor_structure),
    generating_config_endpoints = nrow(generating),
    generating_config_core_hits = sum(generating$core_truth_recovered),
    generating_config_exact_hits =
      sum(generating$exact_thresholded_factor_structure),
    generating_config_core_full_feature_hits =
      sum(generating$core_with_full_feature_coverage),
    generating_config_effective_M_exact =
      sum(generating$generating_effective_M_exact),
    generating_config_split_winner_core_hits =
      sum(generating_winners$core_truth_recovered),
    generating_config_split_winner_exact_hits =
      sum(generating_winners$exact_thresholded_factor_structure),
    generating_strata_with_available_core = sum(gen_strata$any_core_hit),
    generating_strata_where_elbo_missed_core =
      sum(gen_strata$elbo_missed_available_core),
    generating_strata_with_available_exact =
      sum(gen_strata$any_exact_structure_hit),
    generating_strata_where_elbo_missed_exact =
      sum(gen_strata$elbo_missed_available_exact),
    truth_basin_reachable_in_generating_config =
      any(generating$core_truth_recovered),
    exact_structure_reachable_in_generating_config =
      any(generating$exact_thresholded_factor_structure),
    stringsAsFactors = FALSE
  )
}))

candidate_scores <- do.call(rbind, lapply(selections$data_id, function(data_id) {
  value <- read_csv(file.path(selection_root, data_id,
                              "MSE_ONE_SE_CANDIDATES.csv"))
  value$data_id <- data_id
  value
}))
candidate_scores <- candidate_scores[order(candidate_scores$data_id,
                                           candidate_scores$complexity_proxy),
                                     , drop = FALSE]

endpoint_metric_summary <- function(rows) {
  data.frame(
    endpoints = nrow(rows), core_hits = sum(rows$core_truth_recovered),
    exact_structure_hits = sum(rows$exact_thresholded_factor_structure),
    full_feature_coverage_hits =
      sum(rows$full_true_feature_component_coverage),
    ppi_factor_count_exact_hits = sum(rows$ppi_factor_counts_exact),
    observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean =
      safe_mean(rows$dense_signal_nrmse_ppi_selected),
    matched_function_recall_mean =
      safe_mean(rows$matched_function_recall_mean),
    feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
    kernel_relative_ise_mean = safe_mean(rows$kernel_relative_ise_mean),
    covariance_operator_error_mean =
      safe_mean(rows$covariance_operator_relative_error_mean),
    loading_relative_l2_mean =
      safe_mean(rows$loading_relative_l2_error_mean),
    loading_support_auc_mean = safe_mean(rows$loading_support_ppi_auc_mean),
    score_abs_cor_mean = safe_mean(rows$matched_score_abs_cor_mean),
    process_nrmse_missing_as_zero_mean =
      safe_mean(rows$factor_process_nrmse_mean_missing_as_zero),
    contribution_nrmse_missing_as_zero_mean =
      safe_mean(rows$complete_contribution_nrmse_mean_missing_as_zero),
    stringsAsFactors = FALSE
  )
}

selection_assessment <- do.call(rbind, lapply(
  selections$data_id, function(data_id) {
    selection <- selections[selections$data_id == data_id, , drop = FALSE]
    scores <- candidate_scores[candidate_scores$data_id == data_id,
                               , drop = FALSE]
    expected <- unique(all$expected_truth_fit_config_id[all$data_id == data_id])
    selected_score <- scores[scores$fit_config_id ==
                               selection$selected_config_id, , drop = FALSE]
    true_score <- scores[scores$fit_config_id == expected, , drop = FALSE]
    minimum_score <- scores[scores$fit_config_id ==
                              selection$minimum_loss_config_id,
                            , drop = FALSE]
    selected_rows <- all[all$data_id == data_id &
                           all$selected_by_one_se, , drop = FALSE]
    minimum_rows <- all[all$data_id == data_id &
                          all$selected_by_minimum_loss, , drop = FALSE]
    truth_winners <- all[
      all$data_id == data_id &
        all$fit_config_id == expected &
        all$within_configuration_split_winner, , drop = FALSE
    ]
    stopifnot(length(expected) == 1L, nrow(selected_score) == 1L,
              nrow(true_score) == 1L, nrow(minimum_score) == 1L,
              nrow(selected_rows) == 2L, nrow(minimum_rows) == 2L,
              nrow(truth_winners) == 2L)
    data.frame(
      data_id = data_id, truth_config_id = expected,
      one_se_selected_config_id = selection$selected_config_id,
      minimum_loss_config_id = selection$minimum_loss_config_id,
      one_se_dimension_correct = selection$selected_config_id == expected,
      minimum_loss_dimension_correct =
        selection$minimum_loss_config_id == expected,
      one_se_selected_mean_mse = selected_score$mean_loss,
      truth_config_mean_mse = true_score$mean_loss,
      minimum_mean_mse = minimum_score$mean_loss,
      one_se_minus_truth_mse =
        selected_score$mean_loss - true_score$mean_loss,
      truth_minus_minimum_mse =
        true_score$mean_loss - minimum_score$mean_loss,
      truth_config_within_one_se = true_score$within_one_se,
      split_rank_spearman = selection$split_rank_spearman,
      split_minimum_winner_agreement =
        selection$split_minimum_winner_agreement,
      one_se_endpoint_core_hits = sum(selected_rows$core_truth_recovered),
      one_se_endpoint_exact_hits =
        sum(selected_rows$exact_thresholded_factor_structure),
      one_se_endpoint_full_feature_hits =
        sum(selected_rows$full_true_feature_component_coverage),
      minimum_endpoint_core_hits = sum(minimum_rows$core_truth_recovered),
      minimum_endpoint_exact_hits =
        sum(minimum_rows$exact_thresholded_factor_structure),
      truth_config_winner_core_hits =
        sum(truth_winners$core_truth_recovered),
      truth_config_winner_exact_hits =
        sum(truth_winners$exact_thresholded_factor_structure),
      one_se_dense_nrmse_mean =
        safe_mean(selected_rows$dense_signal_nrmse_ppi_selected),
      truth_config_winner_dense_nrmse_mean =
        safe_mean(truth_winners$dense_signal_nrmse_ppi_selected),
      one_se_feature_ise_mean =
        safe_mean(selected_rows$feature_total_ise_mean),
      truth_config_winner_feature_ise_mean =
        safe_mean(truth_winners$feature_total_ise_mean),
      one_se_loading_relative_l2_mean =
        safe_mean(selected_rows$loading_relative_l2_error_mean),
      truth_config_winner_loading_relative_l2_mean =
        safe_mean(truth_winners$loading_relative_l2_error_mean),
      stringsAsFactors = FALSE
    )
  }
))

selected_metric_summary <- do.call(rbind, lapply(
  selections$data_id, function(data_id) {
    expected <- unique(all$expected_truth_fit_config_id[all$data_id == data_id])
    values <- list(
      one_se = all[all$data_id == data_id & all$selected_by_one_se,
                   , drop = FALSE],
      minimum_loss = all[all$data_id == data_id &
                           all$selected_by_minimum_loss, , drop = FALSE],
      truth_configuration_winner = all[
        all$data_id == data_id & all$fit_config_id == expected &
          all$within_configuration_split_winner, , drop = FALSE
      ]
    )
    do.call(rbind, lapply(names(values), function(label) {
      cbind(data.frame(data_id = data_id, endpoint_set = label,
                       stringsAsFactors = FALSE),
            endpoint_metric_summary(values[[label]]))
    }))
  }
))

convergence_summary <- do.call(rbind, lapply(
  split(all, interaction(all$practical_converged,
                         all$generating_dimension_configuration,
                         drop = TRUE, lex.order = TRUE)), function(rows) {
    data.frame(
      scope = "all_scenarios",
      practical_converged = rows$practical_converged[[1L]],
      generating_dimension_configuration =
        rows$generating_dimension_configuration[[1L]],
      endpoints = nrow(rows), core_recovery_rate =
        mean(rows$core_truth_recovered),
      exact_structure_rate = mean(rows$exact_thresholded_factor_structure),
      observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
      dense_signal_nrmse_mean =
        safe_mean(rows$dense_signal_nrmse_ppi_selected),
      feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
      stringsAsFactors = FALSE
    )
  }
))
overall_convergence <- do.call(rbind, lapply(
  split(all, all$practical_converged), function(rows) data.frame(
    scope = "all_endpoints", practical_converged =
      rows$practical_converged[[1L]],
    generating_dimension_configuration = NA,
    endpoints = nrow(rows), core_recovery_rate =
      mean(rows$core_truth_recovered),
    exact_structure_rate = mean(rows$exact_thresholded_factor_structure),
    observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean =
      safe_mean(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
    stringsAsFactors = FALSE
  ))
)
convergence_summary <- rbind(overall_convergence, convergence_summary)

associations <- do.call(rbind, lapply(
  split(all, interaction(all$data_id, all$split_id, all$fit_config_id,
                         drop = TRUE, lex.order = TRUE)), function(rows) {
    data.frame(
      data_id = rows$data_id[[1L]], split_id = rows$split_id[[1L]],
      fit_config_id = rows$fit_config_id[[1L]], endpoints = nrow(rows),
      elbo_vs_observed = safe_cor(rows$final_elbo,
                                  rows$signal_nrmse_ppi_selected),
      elbo_vs_dense = safe_cor(rows$final_elbo,
                               rows$dense_signal_nrmse_ppi_selected),
      elbo_vs_feature_ise = safe_cor(rows$final_elbo,
                                     rows$feature_total_ise_mean),
      stringsAsFactors = FALSE
    )
  }
))

runtime$total_fit_hours <- runtime$endpoints *
  runtime$elapsed_median_seconds / 3600
runtime$median_minutes_per_fit <- runtime$elapsed_median_seconds / 60
runtime$max_minutes_per_fit <- runtime$elapsed_max_seconds / 60
runtime$peak_memory_mib <- runtime$peak_memory_max_bytes / 1024^2

stage_summary <- data.frame(
  endpoints = nrow(all), data_scenarios = length(unique(all$data_id)),
  split_configuration_strata = nrow(split_config_summary),
  objective_eligible = sum(is.finite(all$final_elbo)),
  practical_converged = sum(all$practical_converged),
  core_truth_hits = sum(all$core_truth_recovered),
  exact_thresholded_structure_hits =
    sum(all$exact_thresholded_factor_structure),
  full_feature_coverage_hits =
    sum(all$full_true_feature_component_coverage),
  generating_config_endpoints = sum(all$generating_dimension_configuration),
  generating_config_core_hits = sum(
    all$core_truth_recovered & all$generating_dimension_configuration
  ),
  generating_config_exact_hits = sum(
    all$exact_thresholded_factor_structure &
      all$generating_dimension_configuration
  ),
  generating_config_split_winner_core_hits = sum(
    all$core_truth_recovered & all$generating_dimension_configuration &
      all$within_configuration_split_winner
  ),
  generating_config_split_winner_exact_hits = sum(
    all$exact_thresholded_factor_structure &
      all$generating_dimension_configuration &
      all$within_configuration_split_winner
  ),
  strata_with_available_core = sum(split_config_summary$any_core_hit),
  strata_where_elbo_missed_available_core =
    sum(split_config_summary$elbo_missed_available_core),
  strata_with_available_exact =
    sum(split_config_summary$any_exact_structure_hit),
  strata_where_elbo_missed_available_exact =
    sum(split_config_summary$elbo_missed_available_exact),
  one_se_correct_dimensions = sum(
    selection_assessment$one_se_dimension_correct
  ),
  minimum_loss_correct_dimensions = sum(
    selection_assessment$minimum_loss_dimension_correct
  ),
  evaluation_errors = sum(!status$ok), fit_warnings = sum(all$fit_warning_count),
  new_fits_after_unseal = 0L, continuations_after_unseal = 0L,
  stringsAsFactors = FALSE
)

qc <- data.frame(
  check = c(
    "one_hundred_twenty_endpoints_evaluated", "zero_evaluation_errors",
    "all_metric_contracts_passed", "thirty_split_config_winners",
    "six_one_se_endpoints", "six_minimum_loss_endpoints",
    "three_truth_scenarios", "truth_free_selection_frozen_before_unseal",
    "no_new_fit_after_unseal", "no_continuation_after_unseal"
  ),
  passed = c(
    nrow(all) == 120L, all(status$ok),
    all(metric_qc$metric_contract_passed), nrow(winners) == 30L,
    nrow(one_se_endpoints) == 6L, nrow(minimum_endpoints) == 6L,
    length(unique(all$data_id)) == 3L,
    file.exists(file.path(selection_root,
                          "DIMENSION_SELECTION_FROZEN.txt")),
    length(list.files(file.path(root, "fits"), "^FIT_COMPLETE[.]txt$",
                      recursive = TRUE)) == 120L,
    !any(all$continuation_used)
  ), stringsAsFactors = FALSE
)

write_csv(start_level, "START_LEVEL_DERIVED_SCIENCE.csv")
write_csv(split_config_summary, "DATA_SPLIT_CONFIGURATION_SUMMARY.csv")
write_csv(scenario_config_summary, "SCENARIO_CONFIGURATION_SUMMARY.csv")
write_csv(reachability, "SCENARIO_REACHABILITY_ATTRIBUTION.csv")
write_csv(candidate_scores, "TRUTH_FREE_CANDIDATE_RISKS.csv")
write_csv(selection_assessment,
          "DIMENSION_CRITERION_SCIENTIFIC_ASSESSMENT.csv")
write_csv(selected_metric_summary, "SELECTED_ENDPOINT_METRIC_SUMMARY.csv")
write_csv(convergence_summary, "PRACTICAL_CONVERGENCE_SCIENCE_SUMMARY.csv")
write_csv(associations, "ELBO_SCIENCE_ASSOCIATIONS.csv")
write_csv(runtime, "RUNTIME_INTERPRETATION.csv")
write_csv(stage_summary, "STAGE6G_KEY_RESULTS.csv")
write_csv(qc, "ANALYSIS_QC.csv")
writeLines(c(
  "status=STAGE6G_SCIENTIFIC_ANALYSIS_COMPLETE",
  paste0("completed_utc=",
         format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("endpoints=", stage_summary$endpoints),
  paste0("data_scenarios=", stage_summary$data_scenarios),
  paste0("core_truth_hits=", stage_summary$core_truth_hits),
  paste0("exact_thresholded_structure_hits=",
         stage_summary$exact_thresholded_structure_hits),
  paste0("one_se_correct_dimensions=",
         stage_summary$one_se_correct_dimensions, "/3"),
  paste0("minimum_loss_correct_dimensions=",
         stage_summary$minimum_loss_correct_dimensions, "/3"),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(analysis_root, "STAGE6G_ANALYSIS_COMPLETE.txt"),
useBytes = TRUE)
cat(
  "STAGE6G_ANALYSIS_COMPLETE endpoints=", nrow(all),
  " core_hits=", stage_summary$core_truth_hits,
  " exact_hits=", stage_summary$exact_thresholded_structure_hits,
  " one_se_correct=", stage_summary$one_se_correct_dimensions, "/3",
  " minimum_correct=", stage_summary$minimum_loss_correct_dimensions, "/3\n",
  sep = ""
)
