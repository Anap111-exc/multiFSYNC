#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
evaluation_root <- file.path(root, "evaluation")
selection_root <- file.path(root, "truth_free_selection")
out_root <- file.path(root, "analysis")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(out_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}

all <- read_csv(file.path(evaluation_root,
                          "ALL_STAGE6B_SCIENTIFIC_RESULTS.csv"))
config_winners <- read_csv(file.path(
  evaluation_root, "CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"
))
final <- read_csv(file.path(
  evaluation_root, "FINAL_SELECTED_SCIENTIFIC_RESULT.csv"
))
factor_candidates <- read_csv(file.path(
  selection_root, "FACTOR_DIMENSION_CANDIDATES.csv"
))
fpca_candidates <- read_csv(file.path(
  selection_root, "FPCA_DIMENSION_CANDIDATES.csv"
))
runtime <- read_csv(file.path(root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
status <- read_csv(file.path(evaluation_root, "EVALUATION_STATUS.csv"))
metric_qc <- read_csv(file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))

stopifnot(nrow(all) == 20L, nrow(config_winners) == 5L, nrow(final) == 1L,
          nrow(status) == 20L, all(status$ok), nrow(metric_qc) == 20L,
          all(metric_qc$metric_contract_passed))

all$core_truth_recovered <- with(
  all,
  global_correct == 3L & global_misplaced == 0L & global_missing == 0L &
    study_correct == 4L & study_misplaced == 0L & study_missing == 0L
)
all$exact_factor_structure_recovered <- with(
  all,
  core_truth_recovered & global_extra == 0L & global_duplicate == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L &
    selected_counts == "1,1,1"
)
all$true_effective_M <- with(
  all,
  shared_effective_M == 2L & specific_effective_M_study_1 == 2L &
    specific_effective_M_study_2 == 2L
)

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  stats::cor(x[keep], y[keep], method = method)
}

by_config <- lapply(split(all, all$fit_config_id), function(rows) {
  elbo_row <- rows[which.max(rows$final_elbo), , drop = FALSE]
  holdout_row <- rows[which.min(rows$heldout_nrmse_all), , drop = FALSE]
  dense_row <- rows[which.min(rows$dense_signal_nrmse_ppi_selected), , drop = FALSE]
  data.frame(
    fit_config_id = rows$fit_config_id[[1L]],
    starts = nrow(rows),
    core_truth_recovered = sum(rows$core_truth_recovered),
    exact_factor_structure_recovered =
      sum(rows$exact_factor_structure_recovered),
    true_effective_M = sum(rows$true_effective_M),
    practical_converged = sum(rows$practical_converged),
    elbo_winner_seed = elbo_row$seed_index,
    holdout_best_seed = holdout_row$seed_index,
    dense_truth_best_seed = dense_row$seed_index,
    elbo_and_holdout_same_start = elbo_row$fit_id == holdout_row$fit_id,
    elbo_and_dense_same_start = elbo_row$fit_id == dense_row$fit_id,
    heldout_nrmse_mean = mean(rows$heldout_nrmse_all),
    heldout_nrmse_min = min(rows$heldout_nrmse_all),
    heldout_nrmse_max = max(rows$heldout_nrmse_all),
    observed_signal_nrmse_mean = mean(rows$signal_nrmse_ppi_selected),
    observed_signal_nrmse_min = min(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean = mean(rows$dense_signal_nrmse_ppi_selected),
    dense_signal_nrmse_min = min(rows$dense_signal_nrmse_ppi_selected),
    dense_signal_nrmse_max = max(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = mean(rows$feature_total_ise_mean),
    feature_ise_min = min(rows$feature_total_ise_mean),
    covariance_operator_error_mean =
      mean(rows$covariance_operator_relative_error_mean),
    covariance_operator_error_min =
      min(rows$covariance_operator_relative_error_mean),
    elbo_holdout_spearman = safe_cor(rows$final_elbo,
                                     rows$heldout_nrmse_all),
    elbo_dense_truth_spearman = safe_cor(
      rows$final_elbo, rows$dense_signal_nrmse_ppi_selected
    ),
    stringsAsFactors = FALSE
  )
})
by_config <- do.call(rbind, by_config)
by_config <- by_config[match(unique(all$fit_config_id),
                             by_config$fit_config_id), , drop = FALSE]

derived <- all[, c(
  "fit_id", "fit_config_id", "seed_index", "fit_seed", "final_elbo",
  "heldout_nrmse_all", "heldout_rmse_all", "practical_converged",
  "within_configuration_winner", "final_selected_endpoint",
  "selected_counts", "core_truth_recovered",
  "exact_factor_structure_recovered",
  "true_effective_M", "signal_nrmse_ppi_selected",
  "dense_signal_nrmse_ppi_selected", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "feature_total_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "complete_contribution_nrmse_mean_matched"
)]
derived <- derived[order(derived$fit_config_id, derived$seed_index), , drop = FALSE]

config_winners$core_truth_recovered <- all$core_truth_recovered[
  match(config_winners$fit_id, all$fit_id)
]
config_winners$exact_factor_structure_recovered <-
  all$exact_factor_structure_recovered[
  match(config_winners$fit_id, all$fit_id)
]
config_winners$true_effective_M <- all$true_effective_M[
  match(config_winners$fit_id, all$fit_id)
]
reference <- final[1L, , drop = FALSE]
config_winners$heldout_nrmse_ratio_to_final <-
  config_winners$heldout_nrmse_all / reference$heldout_nrmse_all
config_winners$observed_signal_nrmse_ratio_to_final <-
  config_winners$signal_nrmse_ppi_selected /
  reference$signal_nrmse_ppi_selected
config_winners$dense_signal_nrmse_ratio_to_final <-
  config_winners$dense_signal_nrmse_ppi_selected /
  reference$dense_signal_nrmse_ppi_selected
config_winners$feature_ise_ratio_to_final <-
  config_winners$feature_total_ise_mean / reference$feature_total_ise_mean
config_winners$covariance_error_ratio_to_final <-
  config_winners$covariance_operator_relative_error_mean /
  reference$covariance_operator_relative_error_mean

association <- data.frame(
  scope = c("all_20_fits", "all_20_fits", "all_20_fits", "all_20_fits",
            "five_config_winners", "five_config_winners",
            "five_config_winners"),
  x = c("final_elbo", "final_elbo", "heldout_nrmse", "heldout_nrmse",
        "heldout_nrmse", "heldout_nrmse", "heldout_nrmse"),
  y = c("heldout_nrmse", "dense_signal_nrmse", "dense_signal_nrmse",
        "feature_ise", "dense_signal_nrmse", "feature_ise",
        "covariance_operator_error"),
  spearman = c(
    safe_cor(all$final_elbo, all$heldout_nrmse_all),
    safe_cor(all$final_elbo, all$dense_signal_nrmse_ppi_selected),
    safe_cor(all$heldout_nrmse_all,
             all$dense_signal_nrmse_ppi_selected),
    safe_cor(all$heldout_nrmse_all, all$feature_total_ise_mean),
    safe_cor(config_winners$heldout_nrmse_all,
             config_winners$dense_signal_nrmse_ppi_selected),
    safe_cor(config_winners$heldout_nrmse_all,
             config_winners$feature_total_ise_mean),
    safe_cor(config_winners$heldout_nrmse_all,
             config_winners$covariance_operator_relative_error_mean)
  ),
  stringsAsFactors = FALSE
)

factor_candidates$screen <- "factor_count"
fpca_candidates$screen <- "fpca_cap"
candidate_scores <- rbind(factor_candidates, fpca_candidates)
candidate_scores$relative_mse_vs_screen_min <- ave(
  candidate_scores$mean_subject_time_mse, candidate_scores$screen,
  FUN = function(value) value / min(value) - 1
)

runtime$total_fit_hours <- runtime$fits * runtime$elapsed_median_seconds / 3600
runtime$median_minutes_per_fit <- runtime$elapsed_median_seconds / 60
runtime$max_minutes_per_fit <- runtime$elapsed_max_seconds / 60
runtime$peak_memory_mib <- runtime$peak_memory_max_bytes / 1024^2

qc <- data.frame(
  check = c(
    "twenty_fits_evaluated", "zero_evaluation_errors",
    "all_metric_contracts_passed", "one_final_endpoint",
    "truth_free_selection_chose_true_dimensions",
    "no_continuation", "no_final_full_data_refit"
  ),
  passed = c(
    nrow(all) == 20L, all(status$ok), all(metric_qc$metric_contract_passed),
    nrow(final) == 1L,
    final$fit_L_f == 1L && final$fit_L_s_1 == 1L &&
      final$fit_L_s_2 == 1L && final$fit_M_f == "2" &&
      final$fit_M_s_1 == "2" && final$fit_M_s_2 == "2",
    !any(all$continuation_used),
    !file.exists(file.path(root, "FINAL_FULL_DATA_REFIT_COMPLETE.txt"))
  ),
  stringsAsFactors = FALSE
)

write_csv(derived, "START_LEVEL_SCIENTIFIC_RESULTS.csv")
write_csv(by_config, "CONFIGURATION_ALL_START_SUMMARY.csv")
write_csv(config_winners, "CONFIGURATION_WINNER_COMPARISON.csv")
write_csv(candidate_scores, "TRUTH_FREE_DIMENSION_SCORE_COMPARISON.csv")
write_csv(association, "ELBO_HOLDOUT_TRUTH_ASSOCIATIONS.csv")
write_csv(runtime, "RUNTIME_INTERPRETATION.csv")
write_csv(qc, "ANALYSIS_QC.csv")

cat("STAGE6B_V2_ANALYSIS_COMPLETE fits=", nrow(all),
    " exact_factor_structure=",
    sum(all$exact_factor_structure_recovered),
    " final_fit=", final$fit_id, "\n", sep = "")
