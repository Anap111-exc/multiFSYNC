#!/usr/bin/env Rscript

# Descriptive, post-truth Stage 7A summaries. Never used for endpoint selection.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[[1L]]))
root <- args[[1L]]
read <- function(name) utils::read.csv(file.path(root, name),
                                       stringsAsFactors = FALSE)
all <- read("ALL_START_SCIENCE.csv")
winner <- read("WINNER_SCIENCE.csv")
status <- read("EVALUATION_STATUS.csv")
stopifnot(nrow(all) == 144L, nrow(winner) == 12L,
          nrow(status) == 144L, all(status$ok),
          all(winner$truth_free_winner),
          all(table(all$data_id) == 12L),
          all(table(winner$logical_data_id) == 2L))

out <- file.path(root, "descriptive_analysis_20260920_v1")
stopifnot(!dir.exists(out))
dir.create(out)
write_new <- function(value, name) {
  path <- file.path(out, name)
  stopifnot(!file.exists(path))
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
summary_numeric <- function(x) c(mean = mean(x, na.rm = TRUE),
                                  median = median(x, na.rm = TRUE),
                                  min = min(x, na.rm = TRUE),
                                  max = max(x, na.rm = TRUE))

dataset <- do.call(rbind, lapply(split(all, all$data_id), function(x) {
  w <- x[x$truth_free_winner, , drop = FALSE]
  stopifnot(nrow(w) == 1L)
  correct <- x[x$core_truth_recovered, , drop = FALSE]
  incorrect <- x[!x$core_truth_recovered, , drop = FALSE]
  descending <- sort(x$final_ordinary_elbo, decreasing = TRUE)
  data.frame(
    data_id = x$data_id[[1L]], logical_data_id = x$logical_data_id[[1L]],
    beta_arm = x$beta_arm[[1L]],
    correct_basin_starts = nrow(correct),
    exact_structure_starts = sum(x$exact_thresholded_factor_structure),
    strict_practical_starts = sum(x$strict_practical_converged),
    winner_fit_id = w$fit_id[[1L]], winner_start = w$start_index[[1L]],
    winner_correct_basin = w$core_truth_recovered[[1L]],
    winner_strict = w$strict_practical_converged[[1L]],
    winner_beta_relative_ise = w$beta_relative_ise[[1L]],
    winner_observed_signal_nrmse =
      w$obs_full_signal_nrmse_selected[[1L]],
    winner_dense_signal_nrmse =
      w$dense_full_signal_nrmse_selected[[1L]],
    winner_loading_relative_l2 = w$loading_relative_l2_error_mean[[1L]],
    winner_covariance_operator_relative_error =
      w$covariance_operator_relative_error_mean[[1L]],
    winner_beta_ise_rank_among_starts =
      rank(x$beta_relative_ise, ties.method = "min")[x$truth_free_winner],
    winner_dense_nrmse_rank_among_starts =
      rank(x$dense_full_signal_nrmse_selected,
           ties.method = "min")[x$truth_free_winner],
    elbo_margin_to_second = descending[[1L]] - descending[[2L]],
    correct_minus_best_incorrect_elbo =
      if (nrow(incorrect)) max(correct$final_ordinary_elbo) -
        max(incorrect$final_ordinary_elbo) else NA_real_)
}))
dataset <- dataset[order(dataset$logical_data_id, dataset$beta_arm), ]
write_new(dataset, "DATASET_ARM_SUMMARY.csv")

metrics <- c(
  "beta_relative_ise", "beta_contribution_obs_nrmse",
  "obs_full_signal_nrmse_selected", "dense_full_signal_nrmse_selected",
  "dense_minus_observed_signal_nrmse",
  "loading_relative_l2_error_mean", "correctly_placed_function_recall_mean",
  "feature_total_ise_mean", "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "covariance_operator_relative_error_mean",
  "correctly_placed_score_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "fit_elapsed_seconds")
stopifnot(all(metrics %in% names(winner)))
paired <- merge(
  winner[winner$beta_arm == "sine", c("logical_data_id", metrics)],
  winner[winner$beta_arm == "bspline", c("logical_data_id", metrics)],
  by = "logical_data_id", suffixes = c("_sine", "_bspline"), sort = TRUE)
stopifnot(nrow(paired) == 6L)
paired_summary <- do.call(rbind, lapply(metrics, function(metric) {
  sine <- paired[[paste0(metric, "_sine")]]
  bspline <- paired[[paste0(metric, "_bspline")]]
  delta <- bspline - sine
  data.frame(metric = metric,
             sine_mean = mean(sine), sine_median = median(sine),
             sine_min = min(sine), sine_max = max(sine),
             bspline_mean = mean(bspline), bspline_median = median(bspline),
             bspline_min = min(bspline), bspline_max = max(bspline),
             paired_delta_mean = mean(delta),
             paired_delta_median = median(delta),
             paired_delta_min = min(delta),
             paired_delta_max = max(delta),
             bspline_greater_pairs = sum(delta > 0),
             bspline_less_pairs = sum(delta < 0))
}))
write_new(paired_summary, "PAIRED_WINNER_METRIC_SUMMARY.csv")

arm <- do.call(rbind, lapply(split(all, all$beta_arm), function(x) {
  w <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(beta_arm = x$beta_arm[[1L]], starts = nrow(x),
             datasets = length(unique(x$data_id)),
             core_hits = sum(x$core_truth_recovered),
             exact_structure_hits = sum(x$exact_thresholded_factor_structure),
             full_feature_coverage_hits =
               sum(x$full_true_feature_component_coverage),
             strict_practical_hits = sum(x$strict_practical_converged),
             winner_core_hits = sum(w$core_truth_recovered),
             winner_exact_structure_hits =
               sum(w$exact_thresholded_factor_structure),
             winner_strict_hits = sum(w$strict_practical_converged),
             winner_dense_nrmse_above_one =
               sum(w$dense_full_signal_nrmse_selected > 1),
             winner_beta_ise_median = median(w$beta_relative_ise),
             winner_obs_signal_nrmse_median =
               median(w$obs_full_signal_nrmse_selected),
             winner_dense_signal_nrmse_median =
               median(w$dense_full_signal_nrmse_selected))
}))
write_new(arm, "ARM_SUMMARY.csv")

convergence <- do.call(rbind, lapply(split(all, all$strict_practical_converged),
                                    function(x) {
  data.frame(strict_practical_converged =
               x$strict_practical_converged[[1L]], starts = nrow(x),
             core_hits = sum(x$core_truth_recovered),
             exact_structure_hits = sum(x$exact_thresholded_factor_structure),
             obs_nrmse_median = median(x$obs_full_signal_nrmse_selected),
             dense_nrmse_median =
               median(x$dense_full_signal_nrmse_selected))
}))
write_new(convergence, "STRICT_CONVERGENCE_SCIENCE.csv")

same_starts <- merge(
  all[all$beta_arm == "sine",
      c("logical_data_id", "start_index", "core_truth_recovered")],
  all[all$beta_arm == "bspline",
      c("logical_data_id", "start_index", "core_truth_recovered")],
  by = c("logical_data_id", "start_index"),
  suffixes = c("_sine", "_bspline"))
stopifnot(nrow(same_starts) == 72L)
pairing <- data.frame(
  paired_starts = nrow(same_starts),
  both_correct = sum(same_starts$core_truth_recovered_sine &
                       same_starts$core_truth_recovered_bspline),
  sine_only_correct = sum(same_starts$core_truth_recovered_sine &
                            !same_starts$core_truth_recovered_bspline),
  bspline_only_correct = sum(!same_starts$core_truth_recovered_sine &
                               same_starts$core_truth_recovered_bspline),
  neither_correct = sum(!same_starts$core_truth_recovered_sine &
                          !same_starts$core_truth_recovered_bspline))
write_new(pairing, "PAIRED_START_BASIN_CONCORDANCE.csv")

writeLines(c("status=STAGE7A_DESCRIPTIVE_ANALYSIS_COMPLETE",
             "post_truth_only=TRUE", "used_for_fit_or_selection=FALSE",
             "paired_logical_datasets=6", "starts_evaluated=144",
             "winner_count=12", "no_new_fits_or_continuations=TRUE"),
           file.path(out, "ANALYSIS_COMPLETE.txt"))
cat("STAGE7A_DESCRIPTIVE_ANALYSIS_COMPLETE\n")
