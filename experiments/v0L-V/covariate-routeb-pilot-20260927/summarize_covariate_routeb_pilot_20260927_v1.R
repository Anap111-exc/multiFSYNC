#!/usr/bin/env Rscript

# Descriptive post-evaluation summaries for the covariate Route-B pilot.
# No truth is opened here; this script reads only completed evaluation tables.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
`%||%` <- function(left, right) if (is.null(left)) right else left

parse_cli <- function(arguments) {
  result <- list()
  for (argument in arguments) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE))
      stop("Arguments must use --key=value syntax.", call. = FALSE)
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    result[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
  }
  result
}
assert <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
write_lines_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
}
finite_summary <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (!length(value)) return(c(n = 0, mean = NA, sd = NA, median = NA,
                               q10 = NA, q90 = NA, min = NA, max = NA))
  c(n = length(value), mean = mean(value),
    sd = if (length(value) > 1L) stats::sd(value) else NA_real_,
    median = stats::median(value),
    q10 = as.numeric(stats::quantile(value, 0.1, names = FALSE)),
    q90 = as.numeric(stats::quantile(value, 0.9, names = FALSE)),
    min = min(value), max = max(value))
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
evaluation_name <- cli$`evaluation-name` %||%
  "scientific_evaluation_20260927_v1"
evaluation_root <- file.path(output_root, evaluation_name)
summary_root <- file.path(evaluation_root,
                          "descriptive_summary_20260927_v1")
assert(nzchar(output_root) && dir.exists(evaluation_root),
       "A completed evaluation root is required.")
assert(file.exists(file.path(evaluation_root, "EVALUATION_COMPLETE.txt")),
       "Scientific evaluation is incomplete.")
assert(!dir.exists(summary_root), "Summary output exists; refusing overwrite.")
dir.create(summary_root, recursive = TRUE)

science <- read_csv(file.path(evaluation_root, "ALL_START_SCIENCE.csv"))
winners <- read_csv(file.path(evaluation_root, "WINNER_SCIENCE.csv"))
assert(nrow(science) %in% c(3L, 36L) && nrow(winners) == 3L &&
         setequal(science$beta_arm, c("C0", "C1", "C2")) &&
         all(table(science$beta_arm) == nrow(science) / 3L) &&
         all(winners$truth_free_winner),
       "Science or winner table dimensions are invalid.")

basin <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  core_elbo <- x$final_ordinary_elbo[x$core_truth_recovered]
  noncore_elbo <- x$final_ordinary_elbo[!x$core_truth_recovered]
  data.frame(
    beta_arm = x$beta_arm[[1L]], starts = nrow(x),
    objective_eligible = sum(x$objective_eligible),
    strict_converged = sum(x$strict_practical_converged),
    core_hits = sum(x$core_truth_recovered),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    full_feature_coverage_hits = sum(x$full_true_feature_component_coverage),
    strict_and_core = sum(x$strict_practical_converged & x$core_truth_recovered),
    strict_not_core = sum(x$strict_practical_converged & !x$core_truth_recovered),
    nonstrict_but_core = sum(!x$strict_practical_converged & x$core_truth_recovered),
    neither_strict_nor_core = sum(!x$strict_practical_converged &
                                   !x$core_truth_recovered),
    max_elbo_core = if (length(core_elbo)) max(core_elbo) else NA_real_,
    max_elbo_noncore = if (length(noncore_elbo)) max(noncore_elbo) else NA_real_,
    max_core_minus_max_noncore_elbo = if (length(core_elbo) &&
                                             length(noncore_elbo))
      max(core_elbo) - max(noncore_elbo) else NA_real_,
    winner_fit_id = winner$fit_id[[1L]],
    winner_start_index = winner$start_index[[1L]],
    winner_strict = winner$strict_practical_converged[[1L]],
    winner_core = winner$core_truth_recovered[[1L]],
    winner_exact = winner$exact_thresholded_factor_structure[[1L]],
    winner_full_feature_coverage =
      winner$full_true_feature_component_coverage[[1L]],
    stringsAsFactors = FALSE)
}))
write_csv_new(basin, file.path(summary_root,
                                "ARM_BASIN_AND_CONVERGENCE.csv"))

metrics <- intersect(c(
  "beta_estimated_energy_total", "beta_absolute_ise_total",
  "beta_relative_ise_when_identifiable", "beta_zero_imputed_relative_ise",
  "beta_function_l2_cosine_mean", "beta_contribution_obs_rmse",
  "beta_contribution_obs_nrmse", "beta_contribution_dense_rmse",
  "beta_contribution_dense_nrmse", "obs_full_signal_nrmse_selected",
  "dense_full_signal_nrmse_selected", "obs_residual_signal_nrmse_selected",
  "dense_residual_signal_nrmse_selected",
  "obs_factor_contribution_nrmse_selected",
  "dense_factor_contribution_nrmse_selected", "obs_mean_nrmse",
  "dense_mean_nrmse", "dense_minus_observed_signal_nrmse",
  "factor_beta_absorption_projection_obs",
  "factor_beta_absorption_projection_dense",
  "residual_beta_absorption_projection_obs",
  "residual_beta_absorption_projection_dense",
  "matched_function_recall_mean", "matched_trajectory_abs_cor_mean",
  "matched_score_abs_cor_mean", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean", "loading_support_ppi_auc_mean",
  "factor_process_nrmse_mean_matched",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "mean_dense_nrmse", "mean_projection_relative_ise_mean",
  "mean_excess_ise_mean", "fit_elapsed_seconds"
), names(science))

descriptive <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  do.call(rbind, lapply(metrics, function(metric) {
    value <- finite_summary(x[[metric]])
    data.frame(beta_arm = x$beta_arm[[1L]], metric = metric,
               n = value[["n"]], mean = value[["mean"]],
               sd = value[["sd"]], median = value[["median"]],
               q10 = value[["q10"]], q90 = value[["q90"]],
               min = value[["min"]], max = value[["max"]],
               stringsAsFactors = FALSE)
  }))
}))
write_csv_new(descriptive,
              file.path(summary_root, "ALL_START_METRIC_SUMMARY_BY_ARM.csv"))

winner_columns <- intersect(c(
  "beta_arm", "fit_id", "start_index", "fit_seed", "final_ordinary_elbo",
  "strict_practical_converged", "core_truth_recovered",
  "exact_thresholded_factor_structure", "full_true_feature_component_coverage",
  "selected_shared", "selected_specific_study_1",
  "selected_specific_study_2", "global_correct", "global_misplaced",
  "global_missing", "global_extra", "global_duplicate", "study_correct",
  "study_misplaced", "study_missing", "study_candidate_extra",
  "study_candidate_duplicate", "beta_estimate_available",
  "beta_true_energy_total", "beta_estimated_energy_total",
  "beta_absolute_ise_total", "beta_relative_ise_when_identifiable",
  "beta_zero_imputed_relative_ise", "beta_function_l2_cosine_mean",
  "beta_contribution_obs_rmse", "beta_contribution_obs_nrmse",
  "beta_contribution_dense_nrmse", "obs_full_signal_nrmse_selected",
  "dense_full_signal_nrmse_selected", "obs_residual_signal_nrmse_selected",
  "dense_residual_signal_nrmse_selected",
  "obs_factor_contribution_nrmse_selected",
  "dense_factor_contribution_nrmse_selected",
  "factor_beta_absorption_projection_obs",
  "factor_beta_absorption_projection_dense", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "matched_score_abs_cor_mean", "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched", "mean_dense_nrmse",
  "fit_elapsed_seconds"
), names(winners))
write_csv_new(winners[, winner_columns, drop = FALSE],
              file.path(summary_root, "WINNER_SCIENCE_PRESENTATION.csv"))

pair_metrics <- intersect(c(
  "obs_full_signal_nrmse_selected", "dense_full_signal_nrmse_selected",
  "obs_residual_signal_nrmse_selected",
  "dense_residual_signal_nrmse_selected",
  "obs_factor_contribution_nrmse_selected",
  "dense_factor_contribution_nrmse_selected", "feature_total_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "matched_score_abs_cor_mean", "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched", "mean_dense_nrmse"
), names(science))
c1 <- science[science$beta_arm == "C1", ]
c2 <- science[science$beta_arm == "C2", ]
paired <- merge(c1, c2, by = c("logical_data_id", "start_index", "fit_seed"),
                suffixes = c("_C1", "_C2"), sort = TRUE)
assert(nrow(paired) == nrow(science) / 3L, "Paired C1/C2 join failed.")
pair_difference <- paired[, c("logical_data_id", "start_index", "fit_seed")]
for (metric in pair_metrics) {
  pair_difference[[paste0(metric, "_C2_minus_C1")]] <-
    paired[[paste0(metric, "_C2")]] - paired[[paste0(metric, "_C1")]]
}
pair_difference$core_C1 <- paired$core_truth_recovered_C1
pair_difference$core_C2 <- paired$core_truth_recovered_C2
pair_difference$exact_C1 <- paired$exact_thresholded_factor_structure_C1
pair_difference$exact_C2 <- paired$exact_thresholded_factor_structure_C2
write_csv_new(pair_difference,
              file.path(summary_root, "PAIRED_START_C1_C2_DIFFERENCES.csv"))

paired_summary <- do.call(rbind, lapply(pair_metrics, function(metric) {
  value <- pair_difference[[paste0(metric, "_C2_minus_C1")]]
  valid <- is.finite(value)
  data.frame(
    metric = metric, paired_n = sum(valid),
    mean_C2_minus_C1 = if (any(valid)) mean(value[valid]) else NA_real_,
    median_C2_minus_C1 = if (any(valid)) stats::median(value[valid]) else NA_real_,
    min_C2_minus_C1 = if (any(valid)) min(value[valid]) else NA_real_,
    max_C2_minus_C1 = if (any(valid)) max(value[valid]) else NA_real_,
    C2_worse_count_for_lower_is_better_metric = if (any(valid))
      sum(value[valid] > 0) else NA_integer_,
    stringsAsFactors = FALSE)
}))
write_csv_new(paired_summary,
              file.path(summary_root, "PAIRED_START_C1_C2_SUMMARY.csv"))

rank_rows <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  do.call(rbind, lapply(intersect(c(
    "obs_full_signal_nrmse_selected", "dense_full_signal_nrmse_selected",
    "beta_relative_ise_when_identifiable", "feature_total_ise_mean",
    "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
    "factor_process_nrmse_mean_matched",
    "complete_contribution_nrmse_mean_matched"), names(x)), function(metric) {
      valid <- is.finite(x[[metric]]) & is.finite(x$final_ordinary_elbo)
      data.frame(
        beta_arm = x$beta_arm[[1L]], metric = metric,
        valid_n = sum(valid),
        spearman_elbo_vs_metric = if (sum(valid) >= 3L)
          suppressWarnings(stats::cor(x$final_ordinary_elbo[valid],
                                      x[[metric]][valid], method = "spearman"))
          else NA_real_,
        winner_metric_rank_lower_is_better = {
          winner_value <- x[[metric]][x$truth_free_winner]
          if (length(winner_value) == 1L && is.finite(winner_value))
            rank(x[[metric]], na.last = "keep", ties.method = "min")[x$truth_free_winner]
          else NA_real_
        },
        stringsAsFactors = FALSE)
    }))
  }))
write_csv_new(rank_rows,
              file.path(summary_root, "ELBO_SCIENCE_ASSOCIATION.csv"))

c1_truth_energy <- unique(science$beta_true_energy_total[science$beta_arm == "C1"])
assert(length(c1_truth_energy) == 1L && is.finite(c1_truth_energy),
       "C1 truth beta energy is not constant.")
beta_summary <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(
    beta_arm = x$beta_arm[[1L]], starts = nrow(x),
    estimate_available = unique(x$beta_estimate_available),
    median_estimated_beta_energy = stats::median(
      x$beta_estimated_energy_total, na.rm = TRUE),
    median_estimated_energy_over_C1_truth_energy = stats::median(
      x$beta_estimated_energy_total / c1_truth_energy, na.rm = TRUE),
    median_beta_relative_ise = stats::median(
      x$beta_relative_ise_when_identifiable, na.rm = TRUE),
    median_beta_l2_cosine = stats::median(
      x$beta_function_l2_cosine_mean, na.rm = TRUE),
    winner_estimated_beta_energy = winner$beta_estimated_energy_total[[1L]],
    winner_estimated_energy_over_C1_truth_energy =
      winner$beta_estimated_energy_total[[1L]] / c1_truth_energy,
    winner_beta_relative_ise =
      winner$beta_relative_ise_when_identifiable[[1L]],
    winner_beta_l2_cosine = winner$beta_function_l2_cosine_mean[[1L]],
    stringsAsFactors = FALSE)
}))
# Avoid NaN for the deliberately unavailable C2 beta estimate.
numeric_columns <- vapply(beta_summary, is.numeric, logical(1L))
for (name in names(beta_summary)[numeric_columns]) {
  beta_summary[[name]][is.nan(beta_summary[[name]])] <- NA_real_
}
write_csv_new(beta_summary,
              file.path(summary_root, "BETA_ARM_SUMMARY.csv"))

runtime <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  data.frame(beta_arm = x$beta_arm[[1L]], fits = nrow(x),
             total_fit_hours = sum(x$fit_elapsed_seconds) / 3600,
             median_fit_minutes = stats::median(x$fit_elapsed_seconds) / 60,
             min_fit_minutes = min(x$fit_elapsed_seconds) / 60,
             max_fit_minutes = max(x$fit_elapsed_seconds) / 60,
             stringsAsFactors = FALSE)
}))
write_csv_new(runtime, file.path(summary_root, "RUNTIME_SUMMARY.csv"))

write_lines_new(c(
  "status=COVARIATE_ROUTEB_DESCRIPTIVE_SUMMARY_COMPLETE",
  paste0("fits_summarized=", nrow(science)),
  paste0("winners_summarized=", nrow(winners)),
  "truth_free_winners_unchanged=TRUE",
  "new_fit_or_continuation_started=FALSE",
  "broad_generalization_claim=FALSE"
), file.path(summary_root, "DESCRIPTIVE_SUMMARY_COMPLETE.txt"))
cat("COVARIATE_ROUTEB_SUMMARY_COMPLETE fits=", nrow(science),
    " winners=", nrow(winners), "\n", sep = "")
