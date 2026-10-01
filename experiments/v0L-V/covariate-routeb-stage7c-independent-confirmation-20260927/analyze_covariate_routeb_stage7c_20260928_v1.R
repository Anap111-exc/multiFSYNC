#!/usr/bin/env Rscript

# Read-only post hoc decomposition of the completed Stage 7C evaluation.
# This script does not fit, continue, select, or rewrite any endpoint.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_argument) == 1L)
script_path <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_path)
evaluation_root <- file.path(
  script_root, "scientific_evaluation_20260928_v1")
summary_root <- file.path(evaluation_root,
                          "descriptive_summary_20260928_v1")
output_root <- file.path(evaluation_root,
                         "posthoc_read_only_analysis_20260928_v2")
stopifnot(dir.exists(evaluation_root), dir.exists(summary_root),
          file.exists(file.path(evaluation_root, "EVALUATION_COMPLETE.txt")),
          !dir.exists(output_root))
dir.create(output_root, recursive = TRUE)

read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(output_root, name), row.names = FALSE, quote = TRUE, na = "")
finite_stat <- function(value, fun) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) fun(value) else NA_real_
}

science <- read_csv(file.path(evaluation_root, "ALL_START_SCIENCE.csv"))
winners <- read_csv(file.path(evaluation_root, "WINNER_SCIENCE.csv"))
pool <- read_csv(file.path(summary_root, "POOL_BASIN_AND_CONVERGENCE.csv"))
stopifnot(nrow(science) == 108L, nrow(winners) == 9L, nrow(pool) == 9L,
          all(table(science$data_id) == 12L), all(winners$truth_free_winner))

purity_by_pool <- do.call(rbind, lapply(split(science, science$data_id), function(x) {
  summarize_metric <- function(value, prefix) {
    stats <- c(min = min(value), median = stats::median(value), max = max(value),
               above_half = sum(value >= 0.5))
    stats::setNames(as.list(stats), paste0(prefix, "_", names(stats)))
  }
  data.frame(
    logical_data_id = x$logical_data_id[[1L]],
    data_id = x$data_id[[1L]], beta_arm = x$beta_arm[[1L]],
    summarize_metric(x$P_A, "P_A_subspace_purity"),
    summarize_metric(x$P_B1, "P_B1_subspace_purity"),
    summarize_metric(x$P_B2, "P_B2_subspace_purity"),
    stringsAsFactors = FALSE, check.names = FALSE)
}))
write_csv(purity_by_pool, "SUBSPACE_PURITY_BY_POOL.csv")

pattern_key <- paste(science$selected_shared,
                     science$selected_specific_study_1,
                     science$selected_specific_study_2, sep = "/")
pattern_split <- split(seq_len(nrow(science)),
                       paste(science$data_id, pattern_key, sep = "|||"))
patterns <- do.call(rbind, lapply(pattern_split, function(index) {
  x <- science[index, , drop = FALSE]
  data.frame(
    logical_data_id = x$logical_data_id[[1L]],
    data_id = x$data_id[[1L]], beta_arm = x$beta_arm[[1L]],
    selected_pattern = paste(x$selected_shared[[1L]],
                             x$selected_specific_study_1[[1L]],
                             x$selected_specific_study_2[[1L]], sep = "/"),
    starts = nrow(x), core_hits = sum(x$core_truth_recovered),
    strict_hits = sum(x$strict_practical_converged),
    max_elbo = max(x$final_ordinary_elbo),
    stringsAsFactors = FALSE)
}))
patterns <- patterns[order(patterns$data_id, -patterns$max_elbo), ]
write_csv(patterns, "SELECTED_STRUCTURE_PATTERNS_BY_POOL.csv")

strict_core <- do.call(rbind, lapply(split(science, science$data_id), function(x) {
  data.frame(
    logical_data_id = x$logical_data_id[[1L]],
    data_id = x$data_id[[1L]], beta_arm = x$beta_arm[[1L]],
    strict_and_core = sum(x$strict_practical_converged & x$core_truth_recovered),
    strict_not_core = sum(x$strict_practical_converged & !x$core_truth_recovered),
    nonstrict_but_core = sum(!x$strict_practical_converged & x$core_truth_recovered),
    neither = sum(!x$strict_practical_converged & !x$core_truth_recovered),
    stringsAsFactors = FALSE)
}))
write_csv(strict_core, "STRICT_CONVERGENCE_VS_CORE_BY_POOL.csv")

winner_columns <- intersect(c(
  "logical_data_id", "data_id", "beta_arm", "fit_id", "start_index",
  "strict_practical_converged", "core_truth_recovered",
  "exact_thresholded_factor_structure",
  "full_true_feature_component_coverage", "selected_shared",
  "selected_specific_study_1", "selected_specific_study_2",
  "feature_components_matched", "feature_components_total",
  "beta_estimated_energy_total", "beta_true_energy_total",
  "beta_relative_ise_when_identifiable", "beta_function_l2_cosine_mean",
  "beta_contribution_obs_nrmse", "beta_contribution_dense_nrmse",
  "obs_full_signal_nrmse_selected", "dense_full_signal_nrmse_selected",
  "mean_dense_nrmse", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "loading_relative_l2_error_mean", "matched_score_abs_cor_mean",
  "factor_process_nrmse_mean_matched",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean",
  "fit_elapsed_seconds"), names(winners))
write_csv(winners[, winner_columns, drop = FALSE], "WINNER_KEY_METRICS.csv")

paired <- read_csv(file.path(
  summary_root, "WINNER_C1_C2_CAUSAL_COMPARISON_BY_SEED.csv"))
paired$obs_full_signal_nrmse_difference_C2_minus_C1 <-
  paired$obs_full_signal_nrmse_selected_C2 -
  paired$obs_full_signal_nrmse_selected_C1
paired$dense_full_signal_nrmse_difference_C2_minus_C1 <-
  paired$dense_full_signal_nrmse_selected_C2 -
  paired$dense_full_signal_nrmse_selected_C1
paired$mean_dense_nrmse_difference_C2_minus_C1 <-
  paired$mean_dense_nrmse_C2 - paired$mean_dense_nrmse_C1
paired$complete_contribution_missing_zero_difference_C2_minus_C1 <-
  paired$complete_contribution_nrmse_mean_missing_as_zero_C2 -
  paired$complete_contribution_nrmse_mean_missing_as_zero_C1
write_csv(paired, "WINNER_C1_C2_DECOMPOSITION.csv")

arm_winner <- do.call(rbind, lapply(split(winners, winners$beta_arm), function(x) {
  data.frame(
    beta_arm = x$beta_arm[[1L]], winners = nrow(x),
    core_winners = sum(x$core_truth_recovered),
    exact_winners = sum(x$exact_thresholded_factor_structure),
    full_feature_winners = sum(x$full_true_feature_component_coverage),
    median_obs_full_signal_nrmse =
      finite_stat(x$obs_full_signal_nrmse_selected, stats::median),
    median_dense_full_signal_nrmse =
      finite_stat(x$dense_full_signal_nrmse_selected, stats::median),
    median_mean_dense_nrmse = finite_stat(x$mean_dense_nrmse, stats::median),
    median_covariance_operator_relative_error =
      finite_stat(x$covariance_operator_relative_error_mean, stats::median),
    median_complete_contribution_missing_zero_nrmse = finite_stat(
      x$complete_contribution_nrmse_mean_missing_as_zero, stats::median),
    stringsAsFactors = FALSE)
}))
write_csv(arm_winner, "ARM_WINNER_SCIENCE_SUMMARY.csv")

writeLines(c(
  "status=STAGE7C_READ_ONLY_POSTHOC_ANALYSIS_V2_COMPLETE",
  "fits_read=108", "winners_read=9", "new_fit_started=FALSE",
  "continuation_started=FALSE", "selection_changed=FALSE",
  "truth_free_winners_unchanged=TRUE"),
  file.path(output_root, "POSTHOC_READ_ONLY_ANALYSIS_COMPLETE.txt"))
cat("STAGE7C_POSTHOC_READ_ONLY_ANALYSIS_PASS fits=108 winners=9\n")
