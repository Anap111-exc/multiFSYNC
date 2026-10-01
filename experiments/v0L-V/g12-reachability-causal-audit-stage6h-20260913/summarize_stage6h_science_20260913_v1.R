#!/usr/bin/env Rscript

## Stage 6H post-unseal descriptive scientific summary.
##
## This script is deliberately read-only with respect to fits and evaluations.
## It combines the frozen Stage 6H endpoint table with already-created per-fit
## evaluation tables.  No fit, continuation, endpoint selection, or metric
## recomputation is performed.

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}

stage6h_root <- arg_value(
  "--stage6h-root",
  "/root/multifsync-results/stage6h_20260913_v1"
)
stage6g_eval_root <- arg_value(
  "--stage6g-evaluation-root",
  "/root/multifsync-results/stage6g_20260911_v1/truth_evaluation_20260912_v1"
)
output_dir <- arg_value(
  "--output-dir",
  file.path(stage6h_root, "truth_evaluation_20260913_v1",
            "scientific_analysis_20260913_v2")
)

stage6h_eval_root <- file.path(stage6h_root, "truth_evaluation_20260913_v1")
all_path <- file.path(stage6h_eval_root, "ALL_STAGE6H_SCIENTIFIC_RESULTS.csv")

if (!file.exists(all_path)) stop("Missing frozen Stage 6H endpoint table: ", all_path)
if (!dir.exists(stage6g_eval_root)) stop("Missing Stage 6G evaluation root: ", stage6g_eval_root)
if (dir.exists(output_dir) || file.exists(output_dir)) {
  stop("Refusing to overwrite existing analysis output: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (!file.exists(path)) stop("Missing required evaluation table: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_table <- function(x, name) {
  write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_cor <- function(x, y, method = "spearman") {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  unname(cor(x[keep], y[keep], method = method))
}

all_results <- read_table(all_path)
stopifnot(nrow(all_results) == 36L)

per_fit_tables <- list()
for (i in seq_len(nrow(all_results))) {
  row <- all_results[i, , drop = FALSE]
  eval_root <- if (isTRUE(row$stage6h_new_fit)) stage6h_eval_root else stage6g_eval_root
  eval_dir <- file.path(eval_root, row$fit_id)
  if (!dir.exists(eval_dir)) stop("Missing evaluation directory: ", eval_dir)

  for (nm in c("loading_recovery", "threshold_sensitivity",
               "process_contribution_by_role", "ppi_diagnostics")) {
    tab <- read_table(file.path(eval_dir, paste0(nm, ".csv")))
    tab$fit_id <- row$fit_id
    tab$condition_id <- row$condition_id
    tab$seed_index <- row$seed_index
    tab$stage6h_new_fit <- row$stage6h_new_fit
    per_fit_tables[[nm]][[length(per_fit_tables[[nm]]) + 1L]] <- tab
  }
}

loading <- do.call(rbind, per_fit_tables$loading_recovery)
threshold <- do.call(rbind, per_fit_tables$threshold_sensitivity)
process <- do.call(rbind, per_fit_tables$process_contribution_by_role)
ppi <- do.call(rbind, per_fit_tables$ppi_diagnostics)

loading_groups <- split(loading, interaction(loading$condition_id,
                                             loading$true_id, drop = TRUE))
role_summary <- do.call(rbind, lapply(loading_groups, function(z) {
  data.frame(
    condition_id = z$condition_id[[1L]],
    true_id = z$true_id[[1L]],
    true_block = z$true_block[[1L]],
    true_study = z$true_study[[1L]],
    endpoints = nrow(z),
    correct = sum(z$status == "correct"),
    misplaced = sum(z$status == "misplaced"),
    missing = sum(z$status == "missing"),
    mean_loading_abs_cosine_when_matched = safe_mean(z$loading_abs_cosine[z$matched]),
    stringsAsFactors = FALSE
  )
}))
role_summary <- role_summary[order(role_summary$condition_id,
                                   role_summary$true_study,
                                   role_summary$true_id), ]
write_table(role_summary, "ROLE_DIRECTION_RECOVERY_SUMMARY.csv")

selected_key <- paste(all_results$condition_id, all_results$selected_counts,
                      sep = "\r")
selected_distribution <- as.data.frame(table(selected_key), stringsAsFactors = FALSE)
selected_parts <- strsplit(selected_distribution$selected_key, "\r", fixed = TRUE)
selected_distribution$condition_id <- vapply(selected_parts, `[[`, character(1), 1L)
selected_distribution$selected_counts <- vapply(selected_parts, `[[`, character(1), 2L)
selected_distribution <- selected_distribution[, c("condition_id", "selected_counts", "Freq")]
names(selected_distribution)[3L] <- "endpoints"
selected_distribution <- selected_distribution[selected_distribution$endpoints > 0L, ]
selected_distribution <- selected_distribution[order(selected_distribution$condition_id,
                                                       -selected_distribution$endpoints,
                                                       selected_distribution$selected_counts), ]
write_table(selected_distribution, "SELECTED_COUNT_DISTRIBUTION.csv")

near_key <- paste(all_results$condition_id,
                  all_results$global_correct,
                  all_results$global_misplaced,
                  all_results$global_missing,
                  all_results$study_correct,
                  all_results$study_misplaced,
                  all_results$study_missing,
                  sep = "\r")
near_distribution <- as.data.frame(table(near_key), stringsAsFactors = FALSE)
near_parts <- strsplit(near_distribution$near_key, "\r", fixed = TRUE)
near_distribution <- data.frame(
  condition_id = vapply(near_parts, `[[`, character(1), 1L),
  global_correct = as.integer(vapply(near_parts, `[[`, character(1), 2L)),
  global_misplaced = as.integer(vapply(near_parts, `[[`, character(1), 3L)),
  global_missing = as.integer(vapply(near_parts, `[[`, character(1), 4L)),
  study_correct = as.integer(vapply(near_parts, `[[`, character(1), 5L)),
  study_misplaced = as.integer(vapply(near_parts, `[[`, character(1), 6L)),
  study_missing = as.integer(vapply(near_parts, `[[`, character(1), 7L)),
  endpoints = near_distribution$Freq,
  stringsAsFactors = FALSE
)
near_distribution <- near_distribution[near_distribution$endpoints > 0L, ]
near_distribution <- near_distribution[order(near_distribution$condition_id,
                                               -near_distribution$endpoints), ]
write_table(near_distribution, "IDENTITY_PATTERN_DISTRIBUTION.csv")

threshold_groups <- split(threshold,
                          interaction(threshold$condition_id, threshold$scope,
                                      threshold$threshold, drop = TRUE))
threshold_summary <- do.call(rbind, lapply(threshold_groups, function(z) {
  per_fit <- split(z, z$fit_id)
  core <- vapply(per_fit, function(q) {
    global <- q[q$unit == "global_direction", , drop = FALSE]
    study <- q[q$unit == "study_role", , drop = FALSE]
    nrow(global) == 1L && global$correct == 3L && global$misplaced == 0L &&
      global$missing == 0L && global$extra == 0L && global$duplicate == 0L &&
      nrow(study) == 2L && all(study$correct == 2L) &&
      all(study$misplaced == 0L) && all(study$missing == 0L) &&
      all(study$extra == 0L) && all(study$duplicate == 0L)
  }, logical(1))
  data.frame(
    condition_id = z$condition_id[[1L]],
    scope = z$scope[[1L]],
    threshold = z$threshold[[1L]],
    endpoints = length(core),
    core_hits = sum(core),
    core_hit_rate = mean(core),
    stringsAsFactors = FALSE
  )
}))
threshold_summary <- threshold_summary[order(threshold_summary$condition_id,
                                             threshold_summary$scope,
                                             threshold_summary$threshold), ]
write_table(threshold_summary, "THRESHOLD_CORE_SENSITIVITY_SUMMARY.csv")

metric_names <- c(
  "final_elbo", "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_components_matched", "feature_components_total", "feature_total_ise_mean",
  "factor_kernels_matched", "factor_kernels_total", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "matched_function_recall_mean", "matched_trajectory_abs_cor_mean",
  "matched_score_abs_cor_mean", "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero"
)
metric_rows <- list()
for (condition in unique(all_results$condition_id)) {
  z <- all_results[all_results$condition_id == condition, , drop = FALSE]
  for (metric in metric_names) {
    x <- suppressWarnings(as.numeric(z[[metric]]))
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      condition_id = condition,
      metric = metric,
      n_finite = sum(is.finite(x)),
      mean = if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_,
      sd = if (sum(is.finite(x)) > 1L) sd(x[is.finite(x)]) else NA_real_,
      median = if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_,
      q25 = if (any(is.finite(x))) unname(quantile(x[is.finite(x)], 0.25)) else NA_real_,
      q75 = if (any(is.finite(x))) unname(quantile(x[is.finite(x)], 0.75)) else NA_real_,
      min = if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_,
      max = if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
metric_summary <- do.call(rbind, metric_rows)
write_table(metric_summary, "CONDITION_METRIC_DISTRIBUTIONS.csv")

convergence_groups <- split(all_results,
                            interaction(all_results$condition_id,
                                        all_results$practical_converged,
                                        drop = TRUE))
convergence_summary <- do.call(rbind, lapply(convergence_groups, function(z) {
  data.frame(
    condition_id = z$condition_id[[1L]],
    practical_converged = z$practical_converged[[1L]],
    endpoints = nrow(z),
    core_hits = sum(z$core_truth_recovered),
    ppi_factor_count_exact_hits = sum(z$ppi_factor_counts_exact),
    mean_observed_signal_nrmse = safe_mean(z$signal_nrmse_ppi_selected),
    mean_dense_signal_nrmse = safe_mean(z$dense_signal_nrmse_ppi_selected),
    mean_global_correct = safe_mean(z$global_correct),
    mean_study_correct = safe_mean(z$study_correct),
    stringsAsFactors = FALSE
  )
}))
convergence_summary <- convergence_summary[order(convergence_summary$condition_id,
                                                 -convergence_summary$practical_converged), ]
write_table(convergence_summary, "CONVERGENCE_SCIENCE_SUMMARY.csv")

elbo_rows <- lapply(split(all_results, all_results$condition_id), function(z) {
  data.frame(
    condition_id = z$condition_id[[1L]],
    endpoints = nrow(z),
    spearman_elbo_vs_global_correct = safe_cor(z$final_elbo, z$global_correct),
    spearman_elbo_vs_study_correct = safe_cor(z$final_elbo, z$study_correct),
    spearman_elbo_vs_negative_observed_nrmse = safe_cor(z$final_elbo,
                                                        -z$signal_nrmse_ppi_selected),
    spearman_elbo_vs_negative_dense_nrmse = safe_cor(z$final_elbo,
                                                     -z$dense_signal_nrmse_ppi_selected),
    stringsAsFactors = FALSE
  )
})
write_table(do.call(rbind, elbo_rows), "ELBO_SCIENCE_ASSOCIATIONS.csv")

ppi_groups <- split(ppi, interaction(ppi$condition_id, ppi$block,
                                    ppi$study, drop = TRUE))
ppi_summary <- do.call(rbind, lapply(ppi_groups, function(z) {
  data.frame(
    condition_id = z$condition_id[[1L]],
    block = z$block[[1L]],
    study = z$study[[1L]],
    endpoints = nrow(z),
    retained = sum(z$retained),
    median_factor_ppi = safe_median(z$factor_ppi),
    median_expected_active_loading_count = safe_median(z$expected_active_loading_count),
    stringsAsFactors = FALSE
  )
}))
ppi_summary <- ppi_summary[order(ppi_summary$condition_id,
                                 ppi_summary$block, ppi_summary$study), ]
write_table(ppi_summary, "PPI_BLOCK_SUMMARY.csv")

process_groups <- split(process, interaction(process$condition_id,
                                             process$true_id,
                                             process$study, drop = TRUE))
process_summary <- do.call(rbind, lapply(process_groups, function(z) {
  data.frame(
    condition_id = z$condition_id[[1L]],
    true_id = z$true_id[[1L]],
    study = z$study[[1L]],
    endpoints = nrow(z),
    loading_matched = sum(z$matched_by_loading),
    mean_factor_process_nrmse_missing_as_zero =
      safe_mean(z$factor_process_nrmse_missing_as_zero),
    mean_complete_contribution_nrmse_missing_as_zero =
      safe_mean(z$complete_contribution_nrmse_missing_as_zero),
    stringsAsFactors = FALSE
  )
}))
process_summary <- process_summary[order(process_summary$condition_id,
                                         process_summary$study,
                                         process_summary$true_id), ]
write_table(process_summary, "PROCESS_CONTRIBUTION_BY_ROLE_SUMMARY.csv")

closest <- all_results[order(all_results$condition_id,
                             -all_results$global_correct,
                             -all_results$study_correct,
                             all_results$global_missing,
                             all_results$study_missing,
                             -all_results$final_elbo),
                       c("condition_id", "fit_id", "seed_index", "final_elbo",
                         "practical_converged", "selected_counts", "global_correct",
                         "global_misplaced", "global_missing", "study_correct",
                         "study_misplaced", "study_missing", "core_truth_recovered",
                         "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
                         "feature_components_matched", "feature_components_total",
                         "covariance_operator_relative_error_mean"), drop = FALSE]
closest$descriptive_post_unseal_not_selection <- TRUE
closest <- do.call(rbind, lapply(split(closest, closest$condition_id), head, 3L))
write_table(closest, "CLOSEST_ENDPOINTS_POST_UNSEAL.csv")

qc <- data.frame(
  check_id = c(
    "frozen_endpoint_rows", "new_endpoint_rows", "imported_endpoint_rows",
    "loading_rows", "threshold_rows", "process_rows", "ppi_rows",
    "all_three_truth_directions_per_endpoint", "all_threshold_scopes_present",
    "no_endpoint_selection_recomputed"
  ),
  observed = c(
    nrow(all_results), sum(all_results$stage6h_new_fit),
    sum(!all_results$stage6h_new_fit), nrow(loading), nrow(threshold),
    nrow(process), nrow(ppi),
    all(table(loading$fit_id) == 3L),
    all(c("ppi_selected", "all_candidates") %in% threshold$scope), TRUE
  ),
  expected = c("36", "28", "8", "108", "864", "144", "108",
               "TRUE", "TRUE", "TRUE"),
  passed = c(
    nrow(all_results) == 36L, sum(all_results$stage6h_new_fit) == 28L,
    sum(!all_results$stage6h_new_fit) == 8L, nrow(loading) == 108L,
    nrow(threshold) == 864L, nrow(process) == 144L, nrow(ppi) == 108L,
    all(table(loading$fit_id) == 3L),
    all(c("ppi_selected", "all_candidates") %in% threshold$scope), TRUE
  ),
  stringsAsFactors = FALSE
)
write_table(qc, "ANALYSIS_QC.csv")
if (!all(qc$passed)) stop("Stage 6H descriptive analysis QC failed")

writeLines(c(
  "STAGE6H_POST_UNSEAL_DESCRIPTIVE_ANALYSIS_COMPLETE",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "fits_or_continuations_started=FALSE",
  "endpoint_selection_recomputed=FALSE",
  "truth_used=TRUE",
  "scope=post_unseal_oracle_diagnostic"
), file.path(output_dir, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)

cat("Stage 6H post-unseal descriptive analysis complete: ", output_dir, "\n", sep = "")
