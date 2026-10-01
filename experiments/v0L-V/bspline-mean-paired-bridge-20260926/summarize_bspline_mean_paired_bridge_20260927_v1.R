#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(value, message) if (!isTRUE(value)) abort(message)

arguments <- commandArgs(trailingOnly = TRUE)
root_argument <- grep("^--output-root=", arguments, value = TRUE)
assert(length(root_argument) == 1L,
       "Pass exactly one --output-root=<completed mean bridge result>.")
output_root <- sub("^--output-root=", "", root_argument[[1L]])
evaluation_root <- file.path(output_root,
                             "scientific_evaluation_20260927_v1")
assert(file.exists(file.path(
  evaluation_root, "MEAN_BRIDGE_EVALUATION_COMPLETE.txt")),
  "The completed mean-bridge scientific evaluation is missing.")

read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
finite_mean <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
finite_median <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}
finite_range <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) paste(format(range(value), digits = 8), collapse = ";")
  else ""
}

science <- read_csv(file.path(
  evaluation_root, "ALL_MEAN_BRIDGE_SCIENTIFIC_RESULTS.csv"))
winners <- read_csv(file.path(
  evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))
endpoints <- read_csv(file.path(
  output_root, "truth_free_selection", "ALL_TRUTH_FREE_ENDPOINTS.csv"))
assert(nrow(science) == 48L && nrow(winners) == 4L &&
         nrow(endpoints) == 48L && setequal(science$fit_id, endpoints$fit_id),
       "Endpoint/evaluation accounting changed.")

summary_root <- file.path(evaluation_root,
                          "descriptive_summary_20260927_v1")
assert(!dir.exists(summary_root),
       "Descriptive summary exists; refusing overwrite.")
dir.create(summary_root, recursive = TRUE)

strata <- split(science, science$selection_stratum_id)
stratum_summary <- do.call(rbind, lapply(strata, function(group) {
  group <- group[order(group$start_index), , drop = FALSE]
  ordered <- order(-group$final_ordinary_elbo, group$fit_id)
  winner <- group[group$truth_free_winner, , drop = FALSE]
  correct <- group$core_truth_recovered
  data.frame(
    selection_stratum_id = group$selection_stratum_id[[1L]],
    logical_data_id = group$logical_data_id[[1L]],
    data_id = group$data_id[[1L]], mean_arm = group$mean_arm[[1L]],
    starts = nrow(group), objective_eligible = sum(group$objective_eligible),
    strict_practical = sum(group$strict_practical_converged),
    core_hits = sum(correct),
    core_hit_starts = paste(group$start_index[correct], collapse = ";"),
    ppi_count_exact_hits = sum(group$ppi_factor_counts_exact),
    exact_structure_hits = sum(group$exact_thresholded_factor_structure),
    full_feature_coverage_hits =
      sum(group$full_true_feature_component_coverage),
    winner_fit_id = winner$fit_id[[1L]],
    winner_start = winner$start_index[[1L]],
    winner_core = winner$core_truth_recovered[[1L]],
    winner_exact_structure =
      winner$exact_thresholded_factor_structure[[1L]],
    winner_selected_counts = winner$selected_counts[[1L]],
    winner_feature_coverage = paste0(
      winner$feature_components_matched[[1L]], "/",
      winner$feature_components_total[[1L]]),
    winner_margin_to_second =
      group$final_ordinary_elbo[ordered[[1L]]] -
      group$final_ordinary_elbo[ordered[[2L]]],
    best_correct_minus_best_incorrect = if (any(correct) && any(!correct)) {
      max(group$final_ordinary_elbo[correct]) -
        max(group$final_ordinary_elbo[!correct])
    } else NA_real_,
    winner_mean_observed_nrmse = winner$mean_observed_nrmse[[1L]],
    winner_mean_dense_nrmse = winner$mean_dense_nrmse[[1L]],
    winner_mean_projection_floor_ise =
      winner$mean_projection_floor_ise_mean[[1L]],
    winner_mean_excess_ise = winner$mean_excess_ise_mean[[1L]],
    winner_observed_signal_nrmse =
      winner$signal_nrmse_ppi_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_signal_nrmse_ppi_selected[[1L]],
    winner_feature_total_ise = winner$feature_total_ise_mean[[1L]],
    winner_kernel_relative_ise = winner$kernel_relative_ise_mean[[1L]],
    winner_loading_relative_l2 =
      winner$loading_relative_l2_error_mean[[1L]],
    winner_loading_abs_cosine =
      winner$canonical_loading_abs_cosine_mean[[1L]],
    winner_score_abs_cor =
      winner$correctly_placed_score_abs_cor_mean[[1L]],
    winner_process_nrmse_missing_zero =
      winner$factor_process_nrmse_mean_missing_as_zero[[1L]],
    winner_contribution_nrmse_missing_zero =
      winner$complete_contribution_nrmse_mean_missing_as_zero[[1L]],
    stringsAsFactors = FALSE)
}))
stratum_summary <- stratum_summary[order(
  stratum_summary$logical_data_id, stratum_summary$mean_arm), , drop = FALSE]
write_csv(stratum_summary,
          file.path(summary_root, "STRATUM_BASIN_AND_WINNER_SUMMARY.csv"))

metric_columns <- c(
  mean_observed_nrmse = "mean_observed_nrmse",
  mean_dense_nrmse = "mean_dense_nrmse",
  mean_dense_ise = "mean_dense_ise_mean",
  mean_projection_floor_ise = "mean_projection_floor_ise_mean",
  mean_excess_ise = "mean_excess_ise_mean",
  observed_signal_nrmse = "signal_nrmse_ppi_selected",
  dense_signal_nrmse = "dense_signal_nrmse_ppi_selected",
  feature_total_ise = "feature_total_ise_mean",
  feature_projection_floor_ise = "feature_projection_floor_ise_mean",
  kernel_relative_ise = "kernel_relative_ise_mean",
  covariance_operator_relative_error =
    "covariance_operator_relative_error_mean",
  loading_relative_l2 = "loading_relative_l2_error_mean",
  loading_abs_cosine = "canonical_loading_abs_cosine_mean",
  score_abs_cor = "correctly_placed_score_abs_cor_mean",
  trajectory_abs_cor = "correctly_placed_trajectory_abs_cor_mean",
  process_nrmse_missing_zero =
    "factor_process_nrmse_mean_missing_as_zero",
  contribution_nrmse_missing_zero =
    "complete_contribution_nrmse_mean_missing_as_zero")

arm_summary <- do.call(rbind, lapply(split(science, science$mean_arm), function(x) {
  correct <- x[x$core_truth_recovered, , drop = FALSE]
  winner <- x[x$truth_free_winner, , drop = FALSE]
  per_data_core <- tapply(x$core_truth_recovered, x$logical_data_id, sum)
  base <- data.frame(
    mean_arm = x$mean_arm[[1L]], starts = nrow(x),
    logical_datasets = length(unique(x$logical_data_id)),
    core_hits = nrow(correct), core_hit_rate = nrow(correct) / nrow(x),
    ppi_count_exact_hits = sum(x$ppi_factor_counts_exact),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    strict_practical_hits = sum(x$strict_practical_converged),
    datasets_with_core_endpoint = sum(per_data_core > 0),
    winner_core_hits = sum(winner$core_truth_recovered),
    winner_exact_structure_hits =
      sum(winner$exact_thresholded_factor_structure),
    mean_fit_minutes = mean(x$fit_elapsed_seconds) / 60,
    fit_minutes_range = finite_range(x$fit_elapsed_seconds / 60),
    stringsAsFactors = FALSE)
  for (name in names(metric_columns)) {
    column <- metric_columns[[name]]
    base[[paste0("core_median_", name)]] <- finite_median(correct[[column]])
    base[[paste0("winner_mean_", name)]] <- finite_mean(winner[[column]])
    base[[paste0("winner_range_", name)]] <- finite_range(winner[[column]])
  }
  base
}))
write_csv(arm_summary, file.path(summary_root, "MEAN_ARM_SCIENCE_SUMMARY.csv"))

winner_columns <- c(
  "logical_data_id", "mean_arm", "data_id", "fit_id", "start_index",
  "fit_seed", "core_truth_recovered", "exact_thresholded_factor_structure",
  "selected_counts", "global_correct", "global_misplaced", "global_missing",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "feature_components_matched", "feature_components_total",
  "mean_observed_rmse", "mean_observed_nrmse", "mean_dense_rmse",
  "mean_dense_nrmse", "mean_dense_ise_mean",
  "mean_projection_floor_ise_mean", "mean_excess_ise_mean",
  "mean_estimate_to_projection_ise_mean",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_total_ise_mean", "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "canonical_loading_abs_cosine_mean",
  "correctly_placed_score_abs_cor_mean",
  "correctly_placed_trajectory_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "strict_practical_converged", "fit_elapsed_seconds")
write_csv(winners[, winner_columns, drop = FALSE],
          file.path(summary_root, "WINNER_SCIENTIFIC_PRESENTATION.csv"))

periodic <- science[science$mean_arm == "periodic_phase", , drop = FALSE]
splinet <- science[science$mean_arm == "independent_splinet", , drop = FALSE]
paired <- merge(
  periodic, splinet,
  by = c("logical_data_id", "paired_start_id", "start_index", "fit_seed"),
  suffixes = c("_P", "_B"), all = FALSE, sort = TRUE)
assert(nrow(paired) == 24L,
       "Expected exactly 24 paired start endpoints across mean arms.")
paired$transition <- ifelse(
  paired$core_truth_recovered_P & paired$core_truth_recovered_B,
  "P_core_B_core", ifelse(
    paired$core_truth_recovered_P & !paired$core_truth_recovered_B,
    "P_core_B_not", ifelse(
      !paired$core_truth_recovered_P & paired$core_truth_recovered_B,
      "P_not_B_core", "P_not_B_not")))
paired_metrics <- c(
  "mean_observed_nrmse", "mean_dense_nrmse",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_total_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean",
  "canonical_loading_abs_cosine_mean",
  "correctly_placed_score_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "fit_elapsed_seconds")
for (metric in paired_metrics) {
  paired[[paste0(metric, "_B_minus_P")]] <-
    paired[[paste0(metric, "_B")]] - paired[[paste0(metric, "_P")]]
}
paired_keep <- c(
  "logical_data_id", "paired_start_id", "start_index", "fit_seed",
  "fit_id_P", "fit_id_B", "core_truth_recovered_P",
  "core_truth_recovered_B", "transition",
  unlist(lapply(paired_metrics, function(metric) {
    c(paste0(metric, "_P"), paste0(metric, "_B"),
      paste0(metric, "_B_minus_P"))
  })))
write_csv(paired[, paired_keep, drop = FALSE],
          file.path(summary_root, "PAIRED_START_ENDPOINT_TRANSITIONS.csv"))

transition_levels <- c(
  "P_core_B_core", "P_core_B_not", "P_not_B_core", "P_not_B_not")
transition_summary <- do.call(rbind, lapply(
  c("ALL", sort(unique(paired$logical_data_id))), function(data_id) {
    value <- if (data_id == "ALL") paired else
      paired[paired$logical_data_id == data_id, , drop = FALSE]
    counts <- table(factor(value$transition, levels = transition_levels))
    data.frame(
      logical_data_id = data_id, paired_starts = nrow(value),
      P_core_B_core = unname(counts[[1L]]),
      P_core_B_not = unname(counts[[2L]]),
      P_not_B_core = unname(counts[[3L]]),
      P_not_B_not = unname(counts[[4L]]),
      net_B_minus_P_core_hits = unname(counts[[3L]] - counts[[2L]]),
      stringsAsFactors = FALSE)
  }))
write_csv(transition_summary,
          file.path(summary_root, "PAIRED_CORE_TRANSITION_SUMMARY.csv"))

paired_metric_summary <- do.call(rbind, lapply(paired_metrics, function(metric) {
  difference <- paired[[paste0(metric, "_B_minus_P")]]
  data.frame(
    metric = metric, paired_n = sum(is.finite(difference)),
    B_minus_P_mean = finite_mean(difference),
    B_minus_P_median = finite_median(difference),
    B_minus_P_range = finite_range(difference),
    B_lower_than_P = sum(difference < 0, na.rm = TRUE),
    equal = sum(difference == 0, na.rm = TRUE),
    B_higher_than_P = sum(difference > 0, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
write_csv(paired_metric_summary,
          file.path(summary_root, "PAIRED_METRIC_DIFFERENCE_SUMMARY.csv"))

convergence <- as.data.frame(table(
  mean_arm = science$mean_arm,
  strict_practical_converged = science$strict_practical_converged,
  core_truth_recovered = science$core_truth_recovered),
  stringsAsFactors = FALSE)
names(convergence)[[4L]] <- "fits"
write_csv(convergence,
          file.path(summary_root, "CONVERGENCE_BY_ARM_AND_CORE.csv"))

warning_paths <- file.path(output_root, "fits", science$fit_id, "WARNINGS.csv")
warnings <- do.call(rbind, lapply(seq_along(warning_paths), function(index) {
  value <- read_csv(warning_paths[[index]])
  if (!nrow(value)) return(NULL)
  data.frame(
    fit_id = science$fit_id[[index]],
    logical_data_id = science$logical_data_id[[index]],
    mean_arm = science$mean_arm[[index]],
    start_index = science$start_index[[index]],
    message = value$message, stringsAsFactors = FALSE)
}))
if (is.null(warnings)) warnings <- data.frame(
  fit_id = character(), logical_data_id = character(),
  mean_arm = character(), start_index = integer(), message = character())
write_csv(warnings, file.path(summary_root, "FIT_WARNINGS.csv"))

projection <- do.call(rbind, lapply(strata, function(group) {
  value <- group$mean_projection_floor_ise_mean
  assert(max(value) - min(value) < 1e-12,
         "Mean projection floor changed across starts within a stratum.")
  data.frame(
    logical_data_id = group$logical_data_id[[1L]],
    data_id = group$data_id[[1L]], mean_arm = group$mean_arm[[1L]],
    truth_mean_projection_floor_ise = value[[1L]],
    truth_mean_projection_relative_ise =
      group$mean_projection_relative_ise_mean[[1L]],
    stringsAsFactors = FALSE)
}))
projection <- projection[order(
  projection$logical_data_id, projection$mean_arm), , drop = FALSE]
write_csv(projection,
          file.path(summary_root, "MEAN_PROJECTION_FLOOR_BY_DATA_ARM.csv"))

per_data_arm <- do.call(rbind, lapply(strata, function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(
    logical_data_id = x$logical_data_id[[1L]],
    mean_arm = x$mean_arm[[1L]], core_hits = sum(x$core_truth_recovered),
    core_hit_rate = mean(x$core_truth_recovered),
    winner_core = winner$core_truth_recovered[[1L]],
    winner_mean_dense_nrmse = winner$mean_dense_nrmse[[1L]],
    winner_observed_signal_nrmse =
      winner$signal_nrmse_ppi_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_signal_nrmse_ppi_selected[[1L]],
    stringsAsFactors = FALSE)
}))
per_data_arm <- per_data_arm[order(
  per_data_arm$logical_data_id, per_data_arm$mean_arm), , drop = FALSE]
write_csv(per_data_arm,
          file.path(summary_root, "DATASET_ARM_DECISION_EVIDENCE.csv"))

writeLines(c(
  "status=BSPLINE_MEAN_BRIDGE_DESCRIPTIVE_SUMMARY_COMPLETE",
  paste0("fits=", nrow(science)), paste0("winners=", nrow(winners)),
  paste0("paired_starts=", nrow(paired)),
  paste0("fit_warnings=", nrow(warnings)),
  "new_fit_or_continuation_started=FALSE",
  "truth_used_only_after_frozen_selection=TRUE",
  "cross_mean_arm_ELBO_comparison=FALSE",
  "broad_generalization_claim=FALSE"
), file.path(summary_root, "DESCRIPTIVE_SUMMARY_COMPLETE.txt"),
useBytes = TRUE)

cat("MEAN_BRIDGE_DESCRIPTIVE_SUMMARY_COMPLETE\n")
