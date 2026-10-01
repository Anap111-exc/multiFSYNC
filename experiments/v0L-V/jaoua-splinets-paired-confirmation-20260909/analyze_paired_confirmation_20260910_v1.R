#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
}
output_root <- get_arg("output-root")
if (!nzchar(output_root) || !dir.exists(output_root)) {
  stop("--output-root must identify the paired experiment.", call. = FALSE)
}
evaluation_root <- file.path(output_root, "truth_evaluation_20260910_v1")
analysis_root <- file.path(output_root, "paired_analysis_20260910_v1")
if (!file.exists(file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))) {
  stop("Paired truth evaluation is incomplete.", call. = FALSE)
}
if (dir.exists(analysis_root)) {
  stop("Paired analysis already exists; refusing overwrite.", call. = FALSE)
}
dir.create(analysis_root, recursive = TRUE, mode = "0700")

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(analysis_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

science <- read_csv(file.path(
  evaluation_root, "ALL_PAIRED_SCIENTIFIC_RESULTS.csv"))
winners <- read_csv(file.path(
  output_root, "truth_free_selection", "G12_TRUTH_FREE_WINNERS.csv"))
endpoints <- read_csv(file.path(
  output_root, "truth_free_selection", "ALL_TRUTH_FREE_ENDPOINTS.csv"))
status <- read_csv(file.path(evaluation_root, "EVALUATION_STATUS.csv"))

required <- c(
  "fit_id", "data_id", "logical_data_id", "basis_arm",
  "truth_basis_method", "paired_start_id", "seed_index", "final_elbo",
  "practical_converged", "truth_free_winner", "selected_counts",
  "truth_global_direction_count", "truth_study_role_count",
  "global_correct", "global_misplaced", "global_missing", "global_extra",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "loading_relative_l2_error_mean", "canonical_loading_abs_cosine_mean",
  "feature_components_total", "feature_components_matched",
  "feature_total_ise_mean", "feature_estimate_to_projection_ise_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched", "factor_kernels_total",
  "factor_kernels_matched", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "elapsed_seconds",
  "warning_count"
)
assert(nrow(science) == 72L && all(required %in% names(science)),
       "Paired compact scientific table has the wrong contract.")
assert(nrow(winners) == 6L && nrow(endpoints) == 72L &&
         nrow(status) == 72L && all(status$ok),
       "Paired endpoint/selection/evaluation counts are incomplete.")
assert(all(table(science$basis_arm) == 36L) &&
         all(table(science$logical_data_id, science$basis_arm) == 12L) &&
         all(table(science$paired_start_id) == 2L) &&
         sum(science$truth_free_winner) == 6L,
       "Paired balance or winner mapping failed.")

science$structure_exact <- with(
  science,
  selected_counts == "1,1,1" &
    global_correct == truth_global_direction_count &
    global_misplaced == 0L & global_missing == 0L &
    global_extra == 0L & global_duplicate == 0L &
    study_correct == truth_study_role_count &
    study_misplaced == 0L & study_missing == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L
)
science$structure_class <- ifelse(
  science$structure_exact, "exact_1_1_1_identity",
  ifelse(science$selected_counts == "1,1,1",
         "correct_counts_wrong_identity", "underselected_or_missing")
)
science$feature_component_coverage <- with(
  science, feature_components_matched / feature_components_total)
science$factor_kernel_coverage <- with(
  science, factor_kernels_matched / factor_kernels_total)
science$elbo_rank_within_arm <- ave(
  -science$final_elbo, science$data_id,
  FUN = function(value) rank(value, ties.method = "first"))

write_csv(science, "PER_ENDPOINT_PAIRED_RESULTS.csv")

strata <- do.call(rbind, lapply(split(science, science$data_id), function(x) {
  x <- x[order(x$seed_index), , drop = FALSE]
  winner <- x[x$truth_free_winner, , drop = FALSE]
  exact <- x[x$structure_exact, , drop = FALSE]
  wrong <- x[!x$structure_exact, , drop = FALSE]
  assert(nrow(x) == 12L && nrow(winner) == 1L,
         "Each basis/data stratum must have 12 starts and one winner.")
  exact_best <- if (nrow(exact)) max(exact$final_elbo) else NA_real_
  wrong_best <- if (nrow(wrong)) max(wrong$final_elbo) else NA_real_
  decision <- if (!nrow(exact)) {
    "not_reached_in_12_starts"
  } else if (winner$structure_exact[[1L]]) {
    "reachable_and_elbo_supported"
  } else {
    "reachable_but_elbo_prefers_incorrect"
  }
  data.frame(
    data_id = x$data_id[[1L]], logical_data_id = x$logical_data_id[[1L]],
    basis_arm = x$basis_arm[[1L]],
    truth_basis_method = x$truth_basis_method[[1L]], total_starts = nrow(x),
    exact_hits = nrow(exact), exact_hit_rate = nrow(exact) / nrow(x),
    correct_counts_wrong_identity =
      sum(x$structure_class == "correct_counts_wrong_identity"),
    underselected_or_missing =
      sum(x$structure_class == "underselected_or_missing"),
    practical_converged = sum(x$practical_converged),
    winner_fit_id = winner$fit_id[[1L]],
    winner_seed_index = winner$seed_index[[1L]],
    winner_structure_exact = winner$structure_exact[[1L]],
    winner_structure_class = winner$structure_class[[1L]],
    best_exact_elbo = exact_best, best_incorrect_elbo = wrong_best,
    best_exact_minus_best_incorrect_elbo = exact_best - wrong_best,
    decision_class = decision,
    winner_observed_signal_nrmse =
      winner$signal_nrmse_ppi_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_signal_nrmse_ppi_selected[[1L]],
    winner_loading_relative_l2 =
      winner$loading_relative_l2_error_mean[[1L]],
    winner_loading_abs_cosine =
      winner$canonical_loading_abs_cosine_mean[[1L]],
    winner_feature_component_coverage =
      winner$feature_component_coverage[[1L]],
    winner_feature_total_ise = winner$feature_total_ise_mean[[1L]],
    winner_feature_excess_ise =
      winner$feature_estimate_to_projection_ise_mean[[1L]],
    winner_factor_process_nrmse =
      winner$factor_process_nrmse_mean_matched[[1L]],
    winner_complete_contribution_nrmse =
      winner$complete_contribution_nrmse_mean_matched[[1L]],
    winner_kernel_relative_ise = winner$kernel_relative_ise_mean[[1L]],
    winner_covariance_operator_relative_error =
      winner$covariance_operator_relative_error_mean[[1L]],
    stringsAsFactors = FALSE
  )
}))
rownames(strata) <- NULL
write_csv(strata, "PER_DATA_ARM_REACHABILITY.csv")

metric_names <- c(
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "loading_relative_l2_error_mean", "canonical_loading_abs_cosine_mean",
  "feature_component_coverage", "feature_total_ise_mean",
  "feature_estimate_to_projection_ise_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched", "factor_kernel_coverage",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean"
)

summarise_group <- function(x, endpoint_group) {
  base <- data.frame(
    basis_arm = x$basis_arm[[1L]], endpoint_group = endpoint_group,
    endpoints = nrow(x), exact_hits = sum(x$structure_exact),
    exact_hit_rate = mean(x$structure_exact),
    practical_converged = sum(x$practical_converged),
    practical_convergence_rate = mean(x$practical_converged),
    warning_total = sum(x$warning_count),
    median_elapsed_seconds = safe_median(x$elapsed_seconds),
    total_elapsed_hours = sum(x$elapsed_seconds) / 3600,
    stringsAsFactors = FALSE
  )
  for (metric in metric_names) {
    base[[paste0("median_", metric)]] <- safe_median(x[[metric]])
    base[[paste0("min_", metric)]] <- safe_min(x[[metric]])
    base[[paste0("max_", metric)]] <- safe_max(x[[metric]])
  }
  base
}
arm_summary <- do.call(rbind, lapply(c("legacy", "jaoua"), function(arm) {
  x <- science[science$basis_arm == arm, , drop = FALSE]
  rbind(summarise_group(x, "all_starts"),
        summarise_group(x[x$truth_free_winner, , drop = FALSE], "winners"),
        if (any(x$structure_exact)) {
          summarise_group(x[x$structure_exact, , drop = FALSE], "exact_starts")
        } else NULL)
}))
rownames(arm_summary) <- NULL
write_csv(arm_summary, "ARM_SCIENTIFIC_SUMMARY.csv")

paired <- split(science, science$paired_start_id)
paired_rows <- do.call(rbind, lapply(paired, function(x) {
  x <- x[match(c("legacy", "jaoua"), x$basis_arm), , drop = FALSE]
  assert(nrow(x) == 2L && !any(is.na(x$basis_arm)), "Malformed paired start.")
  out <- data.frame(
    paired_start_id = x$paired_start_id[[1L]],
    logical_data_id = x$logical_data_id[[1L]],
    seed_index = x$seed_index[[1L]],
    legacy_exact = x$structure_exact[[1L]],
    jaoua_exact = x$structure_exact[[2L]],
    transition = if (x$structure_exact[[1L]] && x$structure_exact[[2L]]) {
      "both_exact"
    } else if (!x$structure_exact[[1L]] && x$structure_exact[[2L]]) {
      "jaoua_only_exact"
    } else if (x$structure_exact[[1L]] && !x$structure_exact[[2L]]) {
      "legacy_only_exact"
    } else "neither_exact",
    stringsAsFactors = FALSE
  )
  for (metric in metric_names) {
    out[[paste0("legacy_", metric)]] <- x[[metric]][[1L]]
    out[[paste0("jaoua_", metric)]] <- x[[metric]][[2L]]
    out[[paste0("jaoua_minus_legacy_", metric)]] <-
      x[[metric]][[2L]] - x[[metric]][[1L]]
  }
  out
}))
rownames(paired_rows) <- NULL
write_csv(paired_rows, "PAIRED_START_COMPARISONS.csv")

transition_summary <- as.data.frame(table(paired_rows$transition),
                                    stringsAsFactors = FALSE)
names(transition_summary) <- c("transition", "paired_starts")
transition_summary$proportion <- transition_summary$paired_starts /
  nrow(paired_rows)
write_csv(transition_summary, "PAIRED_STRUCTURE_TRANSITIONS.csv")

paired_metric_summary <- do.call(rbind, lapply(metric_names, function(metric) {
  difference <- paired_rows[[paste0("jaoua_minus_legacy_", metric)]]
  finite <- is.finite(difference)
  data.frame(
    metric = metric, paired_starts = sum(finite),
    mean_jaoua_minus_legacy = safe_mean(difference),
    median_jaoua_minus_legacy = safe_median(difference),
    min_jaoua_minus_legacy = safe_min(difference),
    max_jaoua_minus_legacy = safe_max(difference),
    jaoua_lower = sum(difference[finite] < 0),
    ties = sum(difference[finite] == 0),
    jaoua_higher = sum(difference[finite] > 0),
    stringsAsFactors = FALSE
  )
}))
write_csv(paired_metric_summary, "PAIRED_METRIC_DIFFERENCES.csv")

convergence_summary <- do.call(rbind, lapply(
  split(science, interaction(science$basis_arm,
                             science$practical_converged, drop = TRUE)),
  function(x) data.frame(
    basis_arm = x$basis_arm[[1L]],
    practical_converged = x$practical_converged[[1L]],
    endpoints = nrow(x), exact_hits = sum(x$structure_exact),
    exact_hit_rate = mean(x$structure_exact),
    median_observed_signal_nrmse =
      safe_median(x$signal_nrmse_ppi_selected),
    median_dense_signal_nrmse =
      safe_median(x$dense_signal_nrmse_ppi_selected),
    median_loading_relative_l2 =
      safe_median(x$loading_relative_l2_error_mean),
    median_feature_coverage = safe_median(x$feature_component_coverage),
    stringsAsFactors = FALSE
  )
))
rownames(convergence_summary) <- NULL
write_csv(convergence_summary, "CONVERGENCE_DESCRIPTIVE_SUMMARY.csv")

qc <- data.frame(
  check_id = c(
    "seventy_two_evaluations", "balanced_36_per_arm",
    "three_logical_data_pairs", "twelve_starts_per_data_arm",
    "thirty_six_paired_starts", "six_truth_free_winners_preserved",
    "no_truth_based_fit_or_selection", "pairing_qc_all_passed",
    "development_not_formal"
  ),
  passed = c(
    nrow(science) == 72L && all(status$ok),
    all(table(science$basis_arm) == 36L),
    length(unique(science$logical_data_id)) == 3L,
    all(table(science$data_id) == 12L),
    length(unique(science$paired_start_id)) == 36L &&
      all(table(science$paired_start_id) == 2L),
    sum(science$truth_free_winner) == 6L &&
      setequal(science$fit_id[science$truth_free_winner], winners$fit_id),
    !any(endpoints$truth_used_for_fit_stopping_or_selection),
    all(vapply(list.files(file.path(output_root, "pairing_qc"),
                          "PAIRING_QC[.]csv$", full.names = TRUE),
               function(path) all(read_csv(path)$passed), logical(1L))) &&
      length(list.files(file.path(output_root, "pairing_qc"),
                        "PAIRING_QC[.]csv$")) == 3L,
    !any(endpoints$formal_v0lv_result) &&
      !any(endpoints$formal_paper_mc_result)
  ), stringsAsFactors = FALSE
)
write_csv(qc, "ANALYSIS_QC.csv")
assert(all(qc$passed), "Paired analysis QC failed.")

writeLines(c(
  "status=PASS",
  paste0("completed_utc=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ",
                                  tz = "UTC")),
  "evaluated_fits=72", "basis_arms=2", "logical_data_sets=3",
  "paired_starts=36", "truth_free_winners=6",
  "new_fits_after_truth_unseal=0", "continuation_started=FALSE",
  "development_only=TRUE", "formal_v0lv_result=FALSE",
  "formal_paper_mc_result=FALSE"
), file.path(analysis_root, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
cat("PAIRED_ANALYSIS_PASS fits=72 paired_starts=36 winners=6\n")
