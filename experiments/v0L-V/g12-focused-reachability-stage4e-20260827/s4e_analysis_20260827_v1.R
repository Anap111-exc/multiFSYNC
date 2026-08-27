#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- dirname(script_path)
evaluation_root <- file.path(root, "truth_evaluation_20260827_v1")
output_root <- file.path(root, "analysis_20260827_v1")
if (dir.exists(output_root)) {
  stop("Refusing to overwrite Stage-4E analysis output.", call. = FALSE)
}
dir.create(output_root, recursive = TRUE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(output_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
safe_median <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}

actions <- read_csv(file.path(evaluation_root, "EVALUATION_ACTIONS.csv"))
winners <- read_csv(file.path(root, "TRUTH_FREE_ELBO_WINNERS_8_STARTS.csv"))
inventory <- read_csv(file.path(root, "ALL_24_TRUTH_FREE_ENDPOINTS.csv"))
assert(nrow(actions) == 24L && all(actions$status == "COMPLETE"),
       "Stage-4E endpoint evaluation is incomplete.")
assert(nrow(winners) == 3L && nrow(inventory) == 24L,
       "Stage-4E selection inputs are incomplete.")

summary_paths <- file.path(
  evaluation_root, "endpoints", actions$fit_id, "summary.csv"
)
assert(all(file.exists(summary_paths)), "An endpoint summary is missing.")
summaries <- do.call(rbind, lapply(seq_along(summary_paths), function(index) {
  value <- read_csv(summary_paths[[index]])
  assert(nrow(value) == 1L, "Every endpoint summary must have one row.")
  value$fit_id <- actions$fit_id[[index]]
  value
}))
rownames(summaries) <- NULL

summaries$evaluated_final_elbo <- summaries$final_elbo
summaries$final_elbo <- NULL
endpoint <- merge(actions, summaries, by = "fit_id", sort = FALSE)
endpoint <- endpoint[match(actions$fit_id, endpoint$fit_id), , drop = FALSE]
required <- c(
  "selected_counts", "global_correct", "truth_global_direction_count",
  "global_misplaced", "global_missing", "global_extra", "global_duplicate",
  "study_correct", "truth_study_role_count", "study_misplaced",
  "study_missing", "study_candidate_extra", "study_candidate_duplicate",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "loading_relative_l2_error_mean", "feature_components_matched",
  "feature_components_total", "factor_kernels_matched",
  "factor_kernels_total"
)
assert(all(required %in% names(endpoint)),
       "Frozen evaluator summary lacks a required Stage-4E field.")

endpoint$structure_exact <- with(
  endpoint,
  selected_counts == "1,1,1" &
    global_correct == truth_global_direction_count &
    global_misplaced == 0L & global_missing == 0L &
    global_extra == 0L & global_duplicate == 0L &
    study_correct == truth_study_role_count &
    study_misplaced == 0L & study_missing == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L
)
endpoint$structure_class <- ifelse(
  endpoint$structure_exact, "exact_1_1_1_identity",
  ifelse(endpoint$selected_counts == "1,1,1",
         "correct_counts_wrong_identity", "underselected_or_missing")
)
endpoint$feature_component_coverage <- with(
  endpoint, feature_components_matched / feature_components_total
)
endpoint$factor_kernel_coverage <- with(
  endpoint, factor_kernels_matched / factor_kernels_total
)
endpoint$elbo_identity_preserved <- abs(
  endpoint$frozen_final_elbo - endpoint$evaluated_final_elbo
) < 1e-8
endpoint$elbo_rank_within_data <- ave(
  -endpoint$frozen_final_elbo, endpoint$data_id,
  FUN = function(value) rank(value, ties.method = "first")
)

key_columns <- intersect(c(
  "fit_id", "data_id", "scenario_id", "seed_index", "source_stage",
  "endpoint_role", "elbo_rank_within_data", "frozen_final_elbo",
  "objective_eligible", "practical_converged", "slow_case",
  "selected_counts", "structure_exact", "structure_class",
  "global_correct", "global_misplaced", "global_missing", "global_extra",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "loading_relative_l2_error_mean", "feature_component_coverage",
  "factor_kernel_coverage", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "feature_total_ise_mean", "feature_estimate_to_projection_ise_mean",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean",
  "P_A", "P_B1", "P_B2", "R_A", "R_B1", "R_B2", "max_cross_block_L",
  "elapsed_seconds", "warning_count"
), names(endpoint))
write_csv(endpoint[, key_columns, drop = FALSE],
          "PER_ENDPOINT_REACHABILITY.csv")

data_rows <- lapply(unique(endpoint$data_id), function(data_id) {
  group <- endpoint[endpoint$data_id == data_id, , drop = FALSE]
  group <- group[order(group$seed_index), , drop = FALSE]
  winner <- group[group$endpoint_role == "winner", , drop = FALSE]
  exact <- group[group$structure_exact, , drop = FALSE]
  incorrect <- group[!group$structure_exact, , drop = FALSE]
  assert(nrow(group) == 8L && nrow(winner) == 1L,
         "Every target must have eight endpoints and one winner.")
  exact_best <- if (nrow(exact)) max(exact$frozen_final_elbo) else NA_real_
  wrong_best <- if (nrow(incorrect)) max(incorrect$frozen_final_elbo)
    else NA_real_
  decision <- if (!nrow(exact)) {
    "not_reached_in_eight_starts"
  } else if (winner$structure_exact[[1L]]) {
    "reachable_and_elbo_supported"
  } else {
    "reachable_but_elbo_prefers_incorrect"
  }
  data.frame(
    data_id = data_id, scenario_id = group$scenario_id[[1L]],
    total_starts = nrow(group), inherited_starts = sum(group$seed_index <= 2L),
    new_starts = sum(group$seed_index >= 3L),
    exact_hits_first_two = sum(group$structure_exact[group$seed_index <= 2L]),
    exact_hits_new_six = sum(group$structure_exact[group$seed_index >= 3L]),
    exact_hits_all_eight = nrow(exact),
    exact_hit_rate_all_eight = nrow(exact) / nrow(group),
    winner_fit_id = winner$fit_id[[1L]],
    winner_seed_index = winner$seed_index[[1L]],
    winner_final_elbo = winner$frozen_final_elbo[[1L]],
    winner_structure_exact = winner$structure_exact[[1L]],
    winner_structure_class = winner$structure_class[[1L]],
    best_exact_fit_id = if (nrow(exact))
      exact$fit_id[[which.max(exact$frozen_final_elbo)]] else "",
    best_exact_elbo = exact_best, best_incorrect_elbo = wrong_best,
    best_exact_minus_best_incorrect_elbo = exact_best - wrong_best,
    decision_class = decision,
    winner_observed_signal_nrmse = winner$signal_nrmse_ppi_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_signal_nrmse_ppi_selected[[1L]],
    winner_loading_relative_l2 =
      winner$loading_relative_l2_error_mean[[1L]],
    winner_feature_component_coverage =
      winner$feature_component_coverage[[1L]],
    winner_factor_kernel_coverage = winner$factor_kernel_coverage[[1L]],
    stringsAsFactors = FALSE
  )
})
data_summary <- do.call(rbind, data_rows)
rownames(data_summary) <- NULL
write_csv(data_summary, "PER_DATA_REACHABILITY_DECISION.csv")

group_rows <- lapply(c(FALSE, TRUE), function(exact) {
  group <- endpoint[endpoint$structure_exact == exact, , drop = FALSE]
  data.frame(
    structure_exact = exact, endpoints = nrow(group),
    median_final_elbo = safe_median(group$frozen_final_elbo),
    median_observed_signal_nrmse =
      safe_median(group$signal_nrmse_ppi_selected),
    median_dense_signal_nrmse =
      safe_median(group$dense_signal_nrmse_ppi_selected),
    median_loading_relative_l2 =
      safe_median(group$loading_relative_l2_error_mean),
    median_feature_component_coverage =
      safe_median(group$feature_component_coverage),
    median_factor_kernel_coverage =
      safe_median(group$factor_kernel_coverage),
    stringsAsFactors = FALSE
  )
})
write_csv(do.call(rbind, group_rows), "STRUCTURE_GROUP_SUMMARY.csv")

qc <- data.frame(
  check_id = c(
    "twenty_four_evaluations", "three_targets_eight_starts",
    "three_truth_free_winners_preserved", "no_truth_based_selection",
    "final_elbo_identity_preserved", "structure_class_complete",
    "three_reachability_decisions", "development_not_formal"
  ),
  passed = c(
    nrow(endpoint) == 24L && all(actions$status == "COMPLETE"),
    all(table(endpoint$data_id) == 8L),
    sum(endpoint$endpoint_role == "winner") == 3L &&
      setequal(endpoint$fit_id[endpoint$endpoint_role == "winner"],
               winners$selected_fit_id),
    !any(actions$truth_used_for_selection) &&
      all(actions$selection_frozen_before_truth_read),
    all(endpoint$elbo_identity_preserved),
    !any(is.na(endpoint$structure_class) | endpoint$structure_class == ""),
    nrow(data_summary) == 3L,
    !any(actions$formal_v0lv_result)
  ), stringsAsFactors = FALSE
)
write_csv(qc, "ANALYSIS_QC.csv")
assert(all(qc$passed), "Stage-4E analysis QC failed.")
writeLines(c(
  "status=PASS",
  paste0("completed_utc=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ",
                                  tz = "UTC")),
  "endpoints=24", "target_data_sets=3", "starts_per_data=8",
  "truth_free_selection_preserved=TRUE", "formal_v0lv_result=FALSE"
), file.path(output_root, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE4E_ANALYSIS_PASS endpoints=24 targets=3 starts=8\n")
