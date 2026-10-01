#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else "."
evaluation_root <- file.path(root, "truth_evaluation_20260912_v1")
analysis_root <- file.path(root, "stage6g_scientific_analysis_20260912_v1")
all_path <- file.path(evaluation_root, "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv")
global_path <- file.path(analysis_root,
                         "GLOBAL_TRUTH_DIRECTION_RECOVERY_SUMMARY.csv")
role_path <- file.path(analysis_root, "STUDY_ROLE_RECOVERY_SUMMARY.csv")
if (!file.exists(all_path) || !dir.exists(analysis_root)) {
  stop("Completed Stage 6G evaluation and analysis are required.", call. = FALSE)
}
if (file.exists(global_path) || file.exists(role_path)) {
  stop("Role-recovery summary exists; refusing overwrite.", call. = FALSE)
}
all <- utils::read.csv(all_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
global_rows <- vector("list", nrow(all))
role_rows <- vector("list", nrow(all))
for (index in seq_len(nrow(all))) {
  row <- all[index, , drop = FALSE]
  path <- file.path(evaluation_root, row$fit_id, "evaluation.rds")
  if (!file.exists(path)) stop("Missing evaluation: ", row$fit_id,
                               call. = FALSE)
  evaluation <- readRDS(path)
  global <- evaluation$primary$global_direction_matches
  global <- global[global$scope == "ppi_selected", , drop = FALSE]
  global$data_id <- row$data_id
  global$fit_config_id <- row$fit_config_id
  global$generating_dimension_configuration <-
    row$generating_dimension_configuration
  global_rows[[index]] <- global
  role <- evaluation$primary$functional_recovery
  role <- role[role$scope == "ppi_selected", , drop = FALSE]
  role$data_id <- row$data_id
  role$fit_config_id <- row$fit_config_id
  role$generating_dimension_configuration <-
    row$generating_dimension_configuration
  role_rows[[index]] <- role
}
global_detail <- do.call(rbind, global_rows)
role_detail <- do.call(rbind, role_rows)
safe_mean <- function(value) {
  value <- as.numeric(value); value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}

global_groups <- split(global_detail, interaction(
  global_detail$data_id, global_detail$fit_config_id,
  global_detail$true_id, drop = TRUE, lex.order = TRUE
))
global_summary <- do.call(rbind, lapply(global_groups, function(value) {
  data.frame(
    data_id = value$data_id[[1L]],
    fit_config_id = value$fit_config_id[[1L]],
    generating_dimension_configuration =
      value$generating_dimension_configuration[[1L]],
    true_id = value$true_id[[1L]], true_role = value$true_block[[1L]],
    true_study = value$true_study[[1L]],
    true_factor = value$true_factor[[1L]], endpoints = nrow(value),
    correct = sum(value$status == "correct"),
    misplaced = sum(value$status == "misplaced"),
    missing = sum(value$status == "missing"),
    matched_loading_abs_cosine_mean = safe_mean(
      value$loading_abs_cosine[value$matched]
    ),
    loading_abs_cosine_missing_as_zero_mean =
      mean(value$loading_abs_cosine), stringsAsFactors = FALSE
  )
}))
global_summary <- global_summary[order(
  global_summary$data_id, global_summary$fit_config_id,
  global_summary$true_id
), , drop = FALSE]

role_groups <- split(role_detail, interaction(
  role_detail$data_id, role_detail$fit_config_id, role_detail$study,
  role_detail$true_id, drop = TRUE, lex.order = TRUE
))
role_summary <- do.call(rbind, lapply(role_groups, function(value) {
  matched <- value$matched_by_loading
  data.frame(
    data_id = value$data_id[[1L]],
    fit_config_id = value$fit_config_id[[1L]],
    generating_dimension_configuration =
      value$generating_dimension_configuration[[1L]],
    study = value$study[[1L]], true_id = value$true_id[[1L]],
    true_role = value$true_role[[1L]],
    true_factor = value$true_factor[[1L]], endpoints = nrow(value),
    correct = sum(value$loading_structure_status == "correct"),
    misplaced = sum(value$loading_structure_status == "misplaced"),
    missing = sum(value$loading_structure_status == "missing"),
    matched = sum(matched),
    loading_abs_cosine_mean_matched = safe_mean(
      value$loading_abs_cosine[matched]
    ),
    function_subspace_recall_mean_matched = safe_mean(
      value$function_subspace_recall[matched]
    ),
    trajectory_abs_cor_mean_matched = safe_mean(
      value$trajectory_abs_cor[matched]
    ),
    score_abs_cor_mean_matched = safe_mean(value$score_abs_cor_mean[matched]),
    true_M = value$true_M[[1L]],
    estimated_effective_M_mean_matched = safe_mean(
      value$estimated_effective_M[matched]
    ), stringsAsFactors = FALSE
  )
}))
role_summary <- role_summary[order(
  role_summary$data_id, role_summary$fit_config_id, role_summary$study,
  role_summary$true_id
), , drop = FALSE]

utils::write.csv(global_summary, global_path, row.names = FALSE,
                 quote = TRUE, na = "")
utils::write.csv(role_summary, role_path, row.names = FALSE,
                 quote = TRUE, na = "")
cat("STAGE6G_ROLE_RECOVERY_COMPLETE global_rows=", nrow(global_summary),
    " study_role_rows=", nrow(role_summary), "\n", sep = "")
