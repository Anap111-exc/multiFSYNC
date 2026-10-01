#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(value, message) if (!isTRUE(value)) abort(message)
`%||%` <- function(left, right) if (is.null(left)) right else left

arguments <- commandArgs(trailingOnly = TRUE)
root_argument <- grep("^--output-root=", arguments, value = TRUE)
assert(length(root_argument) == 1L,
       "Pass exactly one --output-root=<completed Route-B confirmation result>.")
output_root <- sub("^--output-root=", "", root_argument[[1L]])
evaluation_root <- file.path(output_root, "scientific_evaluation_20260926_v1")
assert(file.exists(file.path(
  evaluation_root, "ROUTE_B_CONFIRMATION_EVALUATION_COMPLETE.txt")),
  "The completed Route-B confirmation scientific evaluation is missing.")

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
  evaluation_root, "ALL_ROUTE_B_SCIENTIFIC_RESULTS.csv"))
winners <- read_csv(file.path(
  evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))
endpoints <- read_csv(file.path(
  output_root, "truth_free_selection", "ALL_TRUTH_FREE_ENDPOINTS.csv"))
assert(nrow(science) == 48L && nrow(winners) == 6L && nrow(endpoints) == 48L &&
         setequal(science$fit_id, endpoints$fit_id),
       "Endpoint/evaluation accounting changed.")

summary_root <- file.path(evaluation_root,
                          "descriptive_summary_20260926_v1")
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
    data_id = group$data_id[[1L]], fit_K = group$fit_K[[1L]],
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
    winner_observed_signal_nrmse =
      winner$signal_nrmse_ppi_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_signal_nrmse_ppi_selected[[1L]],
    winner_feature_total_ise = winner$feature_total_ise_mean[[1L]],
    winner_projection_floor_ise =
      winner$feature_projection_floor_ise_mean[[1L]],
    winner_kernel_relative_ise = winner$kernel_relative_ise_mean[[1L]],
    winner_covariance_operator_relative_error =
      winner$covariance_operator_relative_error_mean[[1L]],
    winner_loading_relative_l2 = winner$loading_relative_l2_error_mean[[1L]],
    winner_loading_abs_cosine =
      winner$canonical_loading_abs_cosine_mean[[1L]],
    winner_score_abs_cor =
      winner$correctly_placed_score_abs_cor_mean[[1L]],
    winner_trajectory_abs_cor =
      winner$correctly_placed_trajectory_abs_cor_mean[[1L]],
    winner_process_nrmse_missing_zero =
      winner$factor_process_nrmse_mean_missing_as_zero[[1L]],
    winner_contribution_nrmse_missing_zero =
      winner$complete_contribution_nrmse_mean_missing_as_zero[[1L]],
    stringsAsFactors = FALSE)
}))
stratum_summary <- stratum_summary[order(
  stratum_summary$data_id, stratum_summary$fit_K), , drop = FALSE]
write_csv(stratum_summary,
          file.path(summary_root, "STRATUM_BASIN_AND_WINNER_SUMMARY.csv"))

metric_columns <- c(
  observed_signal_nrmse = "signal_nrmse_ppi_selected",
  dense_signal_nrmse = "dense_signal_nrmse_ppi_selected",
  feature_total_ise = "feature_total_ise_mean",
  projection_floor_ise = "feature_projection_floor_ise_mean",
  estimate_to_projection_ise = "feature_estimate_to_projection_ise_mean",
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

k_summary <- do.call(rbind, lapply(split(science, science$fit_K), function(x) {
  correct <- x[x$core_truth_recovered, , drop = FALSE]
  winner <- x[x$truth_free_winner, , drop = FALSE]
  base <- data.frame(
    fit_K = x$fit_K[[1L]], starts = nrow(x), datasets = length(unique(x$data_id)),
    core_hits = nrow(correct), core_hit_rate = nrow(correct) / nrow(x),
    ppi_count_exact_hits = sum(x$ppi_factor_counts_exact),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    strict_practical_hits = sum(x$strict_practical_converged),
    datasets_with_core_endpoint = length(unique(correct$data_id)),
    winner_core_hits = sum(winner$core_truth_recovered),
    mean_fit_minutes = mean(x$fit_elapsed_seconds) / 60,
    fit_minutes_range = finite_range(x$fit_elapsed_seconds / 60),
    stringsAsFactors = FALSE)
  for (name in names(metric_columns)) {
    base[[paste0("correct_median_", name)]] <-
      finite_median(correct[[metric_columns[[name]]]])
    base[[paste0("correct_range_", name)]] <-
      finite_range(correct[[metric_columns[[name]]]])
  }
  base
}))
write_csv(k_summary, file.path(summary_root, "K_CORRECT_ENDPOINT_SUMMARY.csv"))

winner_columns <- c(
  "data_id", "fit_K", "fit_id", "start_index", "fit_seed",
  "core_truth_recovered", "exact_thresholded_factor_structure",
  "selected_counts", "global_correct", "global_misplaced", "global_missing",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "feature_components_matched", "feature_components_total",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_total_ise_mean", "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "canonical_loading_abs_cosine_mean",
  "correctly_placed_score_abs_cor_mean",
  "correctly_placed_trajectory_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero",
  "fit_elapsed_seconds")
write_csv(winners[, winner_columns, drop = FALSE],
          file.path(summary_root, "WINNER_SCIENTIFIC_PRESENTATION.csv"))

projection <- do.call(rbind, lapply(strata, function(group) {
  values <- group$feature_projection_floor_ise_mean
  assert(max(values) - min(values) < 1e-12,
         "Projection floor changed across starts within a data-by-K stratum.")
  data.frame(data_id = group$data_id[[1L]], fit_K = group$fit_K[[1L]],
             truth_projection_floor_ise_mean = values[[1L]])
}))
projection <- projection[order(projection$data_id, projection$fit_K), ]
write_csv(projection,
          file.path(summary_root, "TRUTH_PROJECTION_FLOOR_BY_DATA_AND_K.csv"))

core_long <- science[, c("data_id", "start_index", "fit_K",
                         "core_truth_recovered"), drop = FALSE]
core_wide <- reshape(core_long, idvar = c("data_id", "start_index"),
                     timevar = "fit_K", direction = "wide")
names(core_wide) <- sub("core_truth_recovered[.]", "core_K", names(core_wide))
core_wide <- core_wide[order(core_wide$data_id, core_wide$start_index), ]
write_csv(core_wide, file.path(summary_root, "PAIRED_START_CORE_PATTERN.csv"))

convergence <- as.data.frame.matrix(table(
  strict_practical_converged = science$strict_practical_converged,
  core_truth_recovered = science$core_truth_recovered))
convergence$strict_practical_converged <- rownames(convergence)
rownames(convergence) <- NULL
convergence <- convergence[, c("strict_practical_converged", "FALSE", "TRUE")]
names(convergence)[2:3] <- c("core_false", "core_true")
write_csv(convergence,
          file.path(summary_root, "CONVERGENCE_BY_CORE_RECOVERY.csv"))

warning_paths <- file.path(output_root, "fits", science$fit_id, "WARNINGS.csv")
warnings <- do.call(rbind, lapply(seq_along(warning_paths), function(index) {
  value <- read_csv(warning_paths[[index]])
  if (!nrow(value)) return(NULL)
  data.frame(fit_id = science$fit_id[[index]],
             data_id = science$data_id[[index]],
             fit_K = science$fit_K[[index]],
             start_index = science$start_index[[index]],
             message = value$message, stringsAsFactors = FALSE)
}))
if (is.null(warnings)) warnings <- data.frame(
  fit_id = character(), data_id = character(), fit_K = integer(),
  start_index = integer(), message = character())
write_csv(warnings, file.path(summary_root, "FIT_WARNINGS.csv"))

geometry <- do.call(rbind, lapply(unique(science$data_id), function(data_id) {
  qc <- read_csv(file.path(output_root, "data", data_id, "sealed_truth",
                           "TRUTH_GENERATION_QC.csv"))
  detail <- qc$detail[qc$check_id == "no_exact_cross_block_duplicate"]
  data.frame(
    data_id = data_id,
    max_cross_block_principal_correlation =
      as.numeric(sub("^max_principal_correlation=", "", detail)),
    continuous_orthonormality_passed =
      qc$passed[qc$check_id == "continuous_orthonormality"],
    strong_loading_identification_passed =
      qc$passed[qc$check_id == "strong_loading_identification"],
    stringsAsFactors = FALSE)
}))
write_csv(geometry, file.path(summary_root, "TRUTH_GEOMETRY_CONTEXT.csv"))

selection_gate <- data.frame(
  fit_K = k_summary$fit_K,
  every_dataset_produced_core_endpoint =
    k_summary$datasets_with_core_endpoint == k_summary$datasets,
  every_truth_free_winner_core =
    k_summary$winner_core_hits == k_summary$datasets,
  core_hit_rate = k_summary$core_hit_rate,
  mean_fit_minutes = k_summary$mean_fit_minutes,
  broad_generalization_claim = FALSE)
write_csv(selection_gate,
          file.path(summary_root, "PROTOCOL_DECISION_EVIDENCE.csv"))

writeLines(c(
  "status=ROUTE_B_CONFIRMATION_DESCRIPTIVE_SUMMARY_COMPLETE",
  paste0("fits=", nrow(science)), paste0("winners=", nrow(winners)),
  paste0("fit_warnings=", nrow(warnings)),
  "new_fit_or_continuation_started=FALSE",
  "truth_used_only_after_frozen_selection=TRUE",
  "cross_K_ELBO_comparison=FALSE",
  "broad_generalization_claim=FALSE"
), file.path(summary_root, "DESCRIPTIVE_SUMMARY_COMPLETE.txt"),
useBytes = TRUE)

cat("ROUTE_B_CONFIRMATION_DESCRIPTIVE_SUMMARY_COMPLETE\n")


