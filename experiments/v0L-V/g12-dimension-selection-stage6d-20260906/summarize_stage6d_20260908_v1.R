#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else "."
evaluation_root <- file.path(root, "evaluation")
selection_root <- file.path(root, "truth_free_selection")
out_root <- file.path(root, "analysis")
if (dir.exists(out_root)) {
  stop("Stage 6D analysis directory already exists; refusing overwrite.",
       call. = FALSE)
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(out_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
safe_mean <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
safe_median <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}
safe_min <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) min(value) else NA_real_
}
safe_cor <- function(x, y, method = "spearman") {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  stats::cor(x[keep], y[keep], method = method)
}

all <- read_csv(file.path(evaluation_root,
                          "ALL_STAGE6D_SCIENTIFIC_RESULTS.csv"))
winners <- read_csv(file.path(
  evaluation_root, "CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"
))
axis_selected <- read_csv(file.path(
  evaluation_root, "AXIS_SELECTED_SCIENTIFIC_RESULTS.csv"
))
selections <- read_csv(file.path(selection_root,
                                 "DIMENSION_SELECTIONS.csv"))
runtime <- read_csv(file.path(root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
status <- read_csv(file.path(evaluation_root, "EVALUATION_STATUS.csv"))
metric_qc <- read_csv(file.path(evaluation_root,
                                "METRIC_COMPLETENESS_QC.csv"))

stopifnot(nrow(all) == 180L, nrow(winners) == 15L,
          nrow(axis_selected) == 6L, nrow(selections) == 3L,
          nrow(status) == 180L, all(status$ok),
          nrow(metric_qc) == 180L, all(metric_qc$metric_contract_passed),
          length(unique(all$data_id)) == 3L,
          all(table(all$data_id, all$fit_config_id) == 12L))

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
all$oracle_dimension_configuration <- with(
  all,
  fit_L_f == 1L & fit_L_s_1 == 1L & fit_L_s_2 == 1L &
    fit_M_f == "2" & fit_M_s_1 == "2" & fit_M_s_2 == "2"
)

winner_flags <- all[, c(
  "fit_id", "core_truth_recovered", "exact_factor_structure_recovered",
  "true_effective_M"
)]
winners <- merge(winners, winner_flags, by = "fit_id", all.x = TRUE,
                 suffixes = c("", "_derived"), sort = FALSE)
axis_selected <- merge(axis_selected, winner_flags, by = "fit_id",
                       all.x = TRUE, suffixes = c("", "_derived"),
                       sort = FALSE)
axis_selected$dimension_selected_correctly <- with(
  axis_selected,
  (selection_axis == "factor_count" & selected_config_id == "factor_L1_M2") |
    (selection_axis == "fpca_cap" & selected_config_id == "factor_L1_M2")
)

group_summary <- function(rows) {
  winner <- rows[rows$within_configuration_winner, , drop = FALSE]
  stopifnot(nrow(winner) == 1L)
  dense_best <- rows[which.min(rows$dense_signal_nrmse_ppi_selected),
                     , drop = FALSE]
  holdout_best <- rows[which.min(rows$heldout_nrmse_all), , drop = FALSE]
  data.frame(
    data_id = rows$data_id[[1L]], fit_config_id = rows$fit_config_id[[1L]],
    starts = nrow(rows), practical_converged = sum(rows$practical_converged),
    core_truth_recovered = sum(rows$core_truth_recovered),
    exact_factor_structure_recovered =
      sum(rows$exact_factor_structure_recovered),
    true_effective_M = sum(rows$true_effective_M),
    configurations_with_selected_1_1_1 = sum(rows$selected_counts == "1,1,1"),
    selected_shared_mean = safe_mean(rows$selected_shared),
    selected_specific_1_mean = safe_mean(rows$selected_specific_study_1),
    selected_specific_2_mean = safe_mean(rows$selected_specific_study_2),
    candidate_extra_mean = safe_mean(rows$global_extra +
                                       rows$study_candidate_extra),
    winner_fit_id = winner$fit_id,
    winner_seed_index = winner$seed_index,
    winner_practical_converged = winner$practical_converged,
    winner_core_truth_recovered = winner$core_truth_recovered,
    winner_exact_factor_structure = winner$exact_factor_structure_recovered,
    winner_true_effective_M = winner$true_effective_M,
    winner_heldout_nrmse = winner$heldout_nrmse_all,
    winner_observed_signal_nrmse = winner$signal_nrmse_ppi_selected,
    winner_dense_signal_nrmse = winner$dense_signal_nrmse_ppi_selected,
    winner_feature_ise = winner$feature_total_ise_mean,
    winner_kernel_relative_ise = winner$kernel_relative_ise_mean,
    winner_covariance_operator_error =
      winner$covariance_operator_relative_error_mean,
    winner_loading_relative_l2 = winner$loading_relative_l2_error_mean,
    winner_factor_process_nrmse = winner$factor_process_nrmse_mean_matched,
    winner_complete_contribution_nrmse =
      winner$complete_contribution_nrmse_mean_matched,
    holdout_best_seed_index = holdout_best$seed_index,
    dense_best_seed_index = dense_best$seed_index,
    elbo_and_holdout_same_start = winner$fit_id == holdout_best$fit_id,
    elbo_and_dense_same_start = winner$fit_id == dense_best$fit_id,
    winner_dense_regret = winner$dense_signal_nrmse_ppi_selected -
      dense_best$dense_signal_nrmse_ppi_selected,
    heldout_nrmse_mean = safe_mean(rows$heldout_nrmse_all),
    heldout_nrmse_min = safe_min(rows$heldout_nrmse_all),
    observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean = safe_mean(rows$dense_signal_nrmse_ppi_selected),
    dense_signal_nrmse_min = safe_min(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
    feature_ise_min = safe_min(rows$feature_total_ise_mean),
    covariance_operator_error_mean =
      safe_mean(rows$covariance_operator_relative_error_mean),
    covariance_operator_error_min =
      safe_min(rows$covariance_operator_relative_error_mean),
    elbo_holdout_spearman = safe_cor(rows$final_elbo,
                                     rows$heldout_nrmse_all),
    elbo_dense_spearman = safe_cor(
      rows$final_elbo, rows$dense_signal_nrmse_ppi_selected
    ),
    holdout_dense_spearman = safe_cor(
      rows$heldout_nrmse_all, rows$dense_signal_nrmse_ppi_selected
    ), stringsAsFactors = FALSE
  )
}

by_data_config <- do.call(rbind, lapply(
  split(all, interaction(all$data_id, all$fit_config_id, drop = TRUE)),
  group_summary
))
by_data_config <- by_data_config[order(by_data_config$data_id,
                                       by_data_config$fit_config_id),
                                 , drop = FALSE]

availability <- do.call(rbind, lapply(split(all, all$data_id), function(rows) {
  oracle <- rows[rows$fit_config_id == "factor_L1_M2", , drop = FALSE]
  factor_over <- rows[rows$fit_config_id %in%
                        c("factor_L2_M2", "factor_L3_M2"), , drop = FALSE]
  fpca_over <- rows[rows$fit_config_id %in%
                      c("fpca_L1_M3", "fpca_L1_M4"), , drop = FALSE]
  row_winners <- rows[rows$within_configuration_winner, , drop = FALSE]
  data.frame(
    data_id = rows$data_id[[1L]], starts = nrow(rows),
    any_core_truth_recovered = any(rows$core_truth_recovered),
    any_exact_factor_structure = any(rows$exact_factor_structure_recovered),
    total_core_hits = sum(rows$core_truth_recovered),
    total_exact_hits = sum(rows$exact_factor_structure_recovered),
    oracle_config_core_hits = sum(oracle$core_truth_recovered),
    oracle_config_exact_hits = sum(oracle$exact_factor_structure_recovered),
    factor_overspec_core_hits = sum(factor_over$core_truth_recovered),
    factor_overspec_exact_hits =
      sum(factor_over$exact_factor_structure_recovered),
    fpca_overspec_core_hits = sum(fpca_over$core_truth_recovered),
    fpca_overspec_exact_hits =
      sum(fpca_over$exact_factor_structure_recovered),
    config_winners_core_hits = sum(row_winners$core_truth_recovered),
    config_winners_exact_hits =
      sum(row_winners$exact_factor_structure_recovered),
    configurations_with_any_core = length(unique(
      rows$fit_config_id[rows$core_truth_recovered]
    )),
    configurations_with_any_exact = length(unique(
      rows$fit_config_id[rows$exact_factor_structure_recovered]
    )), stringsAsFactors = FALSE
  )
}))

config_aggregate <- do.call(rbind, lapply(split(all, all$fit_config_id),
                                          function(rows) data.frame(
  fit_config_id = rows$fit_config_id[[1L]], starts = nrow(rows),
  data_sets = length(unique(rows$data_id)),
  core_truth_recovered = sum(rows$core_truth_recovered),
  core_truth_recovery_rate = mean(rows$core_truth_recovered),
  exact_factor_structure_recovered =
    sum(rows$exact_factor_structure_recovered),
  exact_factor_structure_rate = mean(rows$exact_factor_structure_recovered),
  true_effective_M = sum(rows$true_effective_M),
  practical_converged = sum(rows$practical_converged),
  practical_convergence_rate = mean(rows$practical_converged),
  configuration_winner_core_hits = sum(
    rows$core_truth_recovered & rows$within_configuration_winner
  ),
  configuration_winner_exact_hits = sum(
    rows$exact_factor_structure_recovered & rows$within_configuration_winner
  ),
  heldout_nrmse_mean = safe_mean(rows$heldout_nrmse_all),
  observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
  dense_signal_nrmse_mean = safe_mean(rows$dense_signal_nrmse_ppi_selected),
  dense_signal_nrmse_median = safe_median(rows$dense_signal_nrmse_ppi_selected),
  feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
  covariance_operator_error_mean =
    safe_mean(rows$covariance_operator_relative_error_mean),
  loading_relative_l2_mean = safe_mean(rows$loading_relative_l2_error_mean),
  stringsAsFactors = FALSE
)))

convergence_summary <- do.call(rbind, lapply(
  split(all, interaction(all$fit_config_id, all$practical_converged,
                         drop = TRUE)), function(rows) data.frame(
    scope = "by_configuration", fit_config_id = rows$fit_config_id[[1L]],
    practical_converged = rows$practical_converged[[1L]],
    fits = nrow(rows), core_truth_recovery_rate = mean(rows$core_truth_recovered),
    exact_factor_structure_rate = mean(rows$exact_factor_structure_recovered),
    heldout_nrmse_mean = safe_mean(rows$heldout_nrmse_all),
    observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean = safe_mean(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
    covariance_operator_error_mean =
      safe_mean(rows$covariance_operator_relative_error_mean),
    stringsAsFactors = FALSE
  )
))
overall_convergence <- do.call(rbind, lapply(
  split(all, all$practical_converged), function(rows) data.frame(
    scope = "all_fits", fit_config_id = "ALL",
    practical_converged = rows$practical_converged[[1L]],
    fits = nrow(rows), core_truth_recovery_rate = mean(rows$core_truth_recovered),
    exact_factor_structure_rate = mean(rows$exact_factor_structure_recovered),
    heldout_nrmse_mean = safe_mean(rows$heldout_nrmse_all),
    observed_signal_nrmse_mean = safe_mean(rows$signal_nrmse_ppi_selected),
    dense_signal_nrmse_mean = safe_mean(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = safe_mean(rows$feature_total_ise_mean),
    covariance_operator_error_mean =
      safe_mean(rows$covariance_operator_relative_error_mean),
    stringsAsFactors = FALSE
  )
))
convergence_summary <- rbind(overall_convergence, convergence_summary)

associations <- do.call(rbind, lapply(
  split(all, interaction(all$data_id, all$fit_config_id, drop = TRUE)),
  function(rows) data.frame(
    data_id = rows$data_id[[1L]], fit_config_id = rows$fit_config_id[[1L]],
    starts = nrow(rows),
    elbo_vs_holdout = safe_cor(rows$final_elbo, rows$heldout_nrmse_all),
    elbo_vs_observed = safe_cor(rows$final_elbo,
                                rows$signal_nrmse_ppi_selected),
    elbo_vs_dense = safe_cor(rows$final_elbo,
                             rows$dense_signal_nrmse_ppi_selected),
    elbo_vs_feature_ise = safe_cor(rows$final_elbo,
                                   rows$feature_total_ise_mean),
    holdout_vs_dense = safe_cor(rows$heldout_nrmse_all,
                                rows$dense_signal_nrmse_ppi_selected),
    holdout_vs_feature_ise = safe_cor(rows$heldout_nrmse_all,
                                      rows$feature_total_ise_mean),
    stringsAsFactors = FALSE
  )
))

candidate_scores <- do.call(rbind, lapply(selections$data_id, function(data_id) {
  factor <- read_csv(file.path(selection_root, data_id,
                               "FACTOR_CANDIDATES.csv"))
  factor$data_id <- data_id; factor$selection_axis <- "factor_count"
  fpca <- read_csv(file.path(selection_root, data_id,
                             "FPCA_CANDIDATES.csv"))
  fpca$data_id <- data_id; fpca$selection_axis <- "fpca_cap"
  rbind(factor, fpca)
}))
candidate_scores$relative_mse_vs_axis_min <- ave(
  candidate_scores$mean_subject_time_mse,
  interaction(candidate_scores$data_id, candidate_scores$selection_axis),
  FUN = function(value) value / min(value) - 1
)

axis_selected$dimension_selected_correctly <- as.logical(
  axis_selected$dimension_selected_correctly
)
joint_selection <- selections
joint_selection$factor_axis_correct <-
  joint_selection$selected_factor_count == 1L
joint_selection$fpca_axis_correct <- joint_selection$selected_fpca_cap == 2L
joint_selection$both_axes_correct <- with(
  joint_selection, factor_axis_correct & fpca_axis_correct
)

runtime$total_fit_hours_using_median <-
  runtime$fits * runtime$elapsed_median_seconds / 3600
runtime$median_minutes_per_fit <- runtime$elapsed_median_seconds / 60
runtime$max_minutes_per_fit <- runtime$elapsed_max_seconds / 60
runtime$peak_memory_mib <- runtime$peak_memory_max_bytes / 1024^2

derived_columns <- c(
  "fit_id", "data_id", "fit_config_id", "seed_index", "fit_seed",
  "final_elbo", "heldout_nrmse_all", "practical_converged",
  "within_configuration_winner", "factor_axis_selected_endpoint",
  "fpca_axis_selected_endpoint", "selected_counts",
  "core_truth_recovered", "exact_factor_structure_recovered",
  "true_effective_M", "signal_nrmse_ppi_selected",
  "dense_signal_nrmse_ppi_selected", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "feature_components_matched", "feature_components_total",
  "feature_total_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched", "global_correct",
  "global_misplaced", "global_missing", "global_extra", "global_duplicate",
  "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate"
)
derived <- all[derived_columns]
derived <- derived[order(derived$data_id, derived$fit_config_id,
                         derived$seed_index), , drop = FALSE]

qc <- data.frame(
  check = c(
    "one_hundred_eighty_fits_evaluated", "zero_evaluation_errors",
    "all_metric_contracts_passed", "fifteen_config_winners",
    "six_axis_selections", "three_data_sets", "no_continuation",
    "selection_frozen_before_truth", "no_new_fit_after_unseal"
  ),
  passed = c(
    nrow(all) == 180L, all(status$ok),
    all(metric_qc$metric_contract_passed), nrow(winners) == 15L,
    nrow(axis_selected) == 6L, length(unique(all$data_id)) == 3L,
    !any(all$continuation_used),
    file.exists(file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt")),
    length(list.files(file.path(root, "fits"), "^FIT_COMPLETE[.]txt$",
                      recursive = TRUE)) == 180L
  ), stringsAsFactors = FALSE
)

write_csv(derived, "START_LEVEL_SCIENTIFIC_RESULTS.csv")
write_csv(by_data_config, "DATA_CONFIGURATION_SUMMARY.csv")
write_csv(availability, "DATA_CANDIDATE_AVAILABILITY.csv")
write_csv(config_aggregate, "CONFIGURATION_AGGREGATE_SUMMARY.csv")
write_csv(winners, "CONFIGURATION_WINNER_COMPARISON.csv")
write_csv(axis_selected, "AXIS_SELECTED_ENDPOINTS.csv")
write_csv(joint_selection, "TRUTH_FREE_SELECTION_ACCURACY.csv")
write_csv(candidate_scores, "TRUTH_FREE_DIMENSION_SCORE_COMPARISON.csv")
write_csv(convergence_summary, "PRACTICAL_CONVERGENCE_SCIENCE_SUMMARY.csv")
write_csv(associations, "ELBO_HOLDOUT_TRUTH_ASSOCIATIONS.csv")
write_csv(runtime, "RUNTIME_INTERPRETATION.csv")
write_csv(qc, "ANALYSIS_QC.csv")

summary <- data.frame(
  fits = nrow(all), data_sets = length(unique(all$data_id)),
  core_truth_hits = sum(all$core_truth_recovered),
  exact_factor_structure_hits = sum(all$exact_factor_structure_recovered),
  configuration_winner_core_hits = sum(
    all$core_truth_recovered & all$within_configuration_winner
  ),
  configuration_winner_exact_hits = sum(
    all$exact_factor_structure_recovered & all$within_configuration_winner
  ),
  factor_axis_correct_data_sets = sum(joint_selection$factor_axis_correct),
  fpca_axis_correct_data_sets = sum(joint_selection$fpca_axis_correct),
  both_axes_correct_data_sets = sum(joint_selection$both_axes_correct),
  practical_converged = sum(all$practical_converged),
  evaluation_errors = sum(!status$ok), stringsAsFactors = FALSE
)
write_csv(summary, "STAGE6D_KEY_RESULTS.csv")
writeLines(c(
  "status=STAGE6D_ANALYSIS_COMPLETE",
  paste0("completed_utc=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("fits=", summary$fits), paste0("data_sets=", summary$data_sets),
  paste0("core_truth_hits=", summary$core_truth_hits),
  paste0("exact_factor_structure_hits=",
         summary$exact_factor_structure_hits),
  paste0("both_axes_correct_data_sets=",
         summary$both_axes_correct_data_sets),
  "new_fits_after_unseal=0", "continuation_started=FALSE",
  "formal_paper_mc_result=FALSE"
), file.path(out_root, "STAGE6D_ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE6D_ANALYSIS_COMPLETE fits=", nrow(all),
    " core_hits=", sum(all$core_truth_recovered),
    " exact_hits=", sum(all$exact_factor_structure_recovered),
    " both_axes_correct=", sum(joint_selection$both_axes_correct), "/3\n",
    sep = "")
