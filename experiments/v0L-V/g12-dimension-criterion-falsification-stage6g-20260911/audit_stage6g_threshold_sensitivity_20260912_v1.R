#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
root <- if (length(arguments)) arguments[[1L]] else "."
evaluation_root <- file.path(root, "truth_evaluation_20260912_v1")
analysis_root <- file.path(root, "stage6g_scientific_analysis_20260912_v1")
all_path <- file.path(evaluation_root, "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv")
output_path <- file.path(analysis_root,
                         "LOADING_THRESHOLD_SENSITIVITY_SUMMARY.csv")
if (!file.exists(all_path) || !dir.exists(analysis_root)) {
  stop("Completed Stage 6G evaluation and analysis are required.", call. = FALSE)
}
if (file.exists(output_path)) {
  stop("Threshold-sensitivity summary exists; refusing overwrite.",
       call. = FALSE)
}
all <- utils::read.csv(all_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
rows <- vector("list", nrow(all))
for (index in seq_len(nrow(all))) {
  fit_id <- all$fit_id[[index]]
  path <- file.path(evaluation_root, fit_id, "threshold_sensitivity.csv")
  if (!file.exists(path)) stop("Missing threshold table: ", fit_id,
                               call. = FALSE)
  value <- utils::read.csv(path, stringsAsFactors = FALSE,
                           check.names = FALSE)
  value <- value[value$scope == "ppi_selected", , drop = FALSE]
  thresholds <- sort(unique(value$threshold))
  rows[[index]] <- do.call(rbind, lapply(thresholds, function(threshold) {
    at_threshold <- value[value$threshold == threshold, , drop = FALSE]
    global <- at_threshold[at_threshold$unit == "global_direction",
                           , drop = FALSE]
    study <- at_threshold[at_threshold$unit == "study_role", , drop = FALSE]
    if (nrow(global) != 1L || nrow(study) != 2L) {
      stop("Malformed threshold rows: ", fit_id, call. = FALSE)
    }
    data.frame(
      fit_id = fit_id, data_id = all$data_id[[index]],
      fit_config_id = all$fit_config_id[[index]],
      generating_dimension_configuration =
        all$generating_dimension_configuration[[index]],
      threshold = threshold,
      threshold_core_recovered =
        global$correct == all$truth_global_direction_count[[index]] &&
        global$misplaced == 0L && global$missing == 0L &&
        sum(study$correct) == all$truth_study_role_count[[index]] &&
        sum(study$misplaced) == 0L && sum(study$missing) == 0L,
      stringsAsFactors = FALSE
    )
  }))
}
detail <- do.call(rbind, rows)
groups <- split(detail, interaction(
  detail$data_id, detail$generating_dimension_configuration,
  detail$threshold, drop = TRUE, lex.order = TRUE
))
summary <- do.call(rbind, lapply(groups, function(value) data.frame(
  data_id = value$data_id[[1L]],
  generating_dimension_configuration =
    value$generating_dimension_configuration[[1L]],
  loading_abs_cosine_threshold = value$threshold[[1L]],
  endpoints = nrow(value),
  core_recovered = sum(value$threshold_core_recovered),
  core_recovery_rate = mean(value$threshold_core_recovered),
  stringsAsFactors = FALSE
)))
summary <- summary[order(summary$data_id,
                         summary$generating_dimension_configuration,
                         summary$loading_abs_cosine_threshold), , drop = FALSE]
utils::write.csv(summary, output_path, row.names = FALSE, quote = TRUE,
                 na = "")
cat("STAGE6G_THRESHOLD_SENSITIVITY_COMPLETE rows=", nrow(summary),
    " core_hits=", sum(detail$threshold_core_recovered), "\n", sep = "")
