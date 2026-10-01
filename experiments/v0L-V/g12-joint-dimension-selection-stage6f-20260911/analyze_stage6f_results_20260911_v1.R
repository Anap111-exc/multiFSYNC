#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
get_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  value <- arguments[startsWith(arguments, prefix)]
  if (length(value) != 1L) {
    stop(paste0("Exactly one ", prefix, " argument is required."),
         call. = FALSE)
  }
  sub(prefix, "", value[[1L]], fixed = TRUE)
}
output_root <- get_argument("output-root")
analysis_root <- get_argument("analysis-root")
if (!dir.exists(output_root)) stop("Stage 6F output root is missing.",
                                   call. = FALSE)
if (dir.exists(analysis_root) || file.exists(analysis_root)) {
  stop("Analysis root already exists; refusing overwrite.", call. = FALSE)
}
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop(paste0("Missing file: ", path), call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(analysis_root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
assert <- function(value, message) if (!isTRUE(value)) {
  stop(message, call. = FALSE)
}
marker_value <- function(path, key) {
  lines <- readLines(path, warn = FALSE)
  prefix <- paste0(key, "=")
  value <- lines[startsWith(lines, prefix)]
  if (!length(value)) return(NA_character_)
  sub(prefix, "", value[[1L]], fixed = TRUE)
}

required_markers <- c(
  "STAGE6F_COMPLETE.txt", "STAGE6F_TRUTH_FREE_COMPLETE.txt",
  file.path("truth_free_selection", "START_SELECTION_FROZEN.txt"),
  file.path("truth_free_selection", "DIMENSION_SELECTION_FROZEN.txt"),
  file.path("truth_evaluation", "EVALUATION_COMPLETE.txt")
)
assert(all(file.exists(file.path(output_root, required_markers))),
       "Stage 6F completion/freeze markers are incomplete.")
fit_errors <- list.files(output_root, "^FIT_ERROR[.]txt$", recursive = TRUE)
score_errors <- list.files(output_root, "^SCORE_ERROR[.]txt$", recursive = TRUE)
evaluation_errors <- list.files(output_root, "^EVALUATION_ERROR[.]txt$",
                                recursive = TRUE)
assert(!length(c(fit_errors, score_errors, evaluation_errors)),
       "Stage 6F contains retained errors.")

manifest <- read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
endpoints <- read_csv(file.path(output_root, "truth_free_selection",
                                "ALL_ENDPOINTS.csv"))
winners <- read_csv(file.path(output_root, "truth_free_selection",
                              "CONFIG_SPLIT_WINNERS.csv"))
selection <- read_csv(file.path(output_root, "truth_free_selection",
                                "FINAL_TRUTH_FREE_DIMENSION_SELECTION.csv"))
mse_candidates <- read_csv(file.path(output_root, "truth_free_selection",
                                     "MSE_ONE_SE_CANDIDATES.csv"))
logscore_candidates <- read_csv(file.path(
  output_root, "truth_free_selection", "LOGSCORE_ONE_SE_CANDIDATES.csv"
))
calibration <- read_csv(file.path(output_root, "truth_free_selection",
                                  "PREDICTIVE_CALIBRATION.csv"))
split_risks <- read_csv(file.path(output_root, "truth_free_selection",
                                  "SPLIT_SPECIFIC_RISKS.csv"))
rank_stability <- read_csv(file.path(
  output_root, "truth_free_selection", "CROSS_SPLIT_RANK_STABILITY.csv"
))
science <- read_csv(file.path(
  output_root, "truth_evaluation",
  "SPLIT_CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"
))
science_by_config <- read_csv(file.path(
  output_root, "truth_evaluation",
  "CONFIG_SCIENTIFIC_RESULTS_ACROSS_SPLITS.csv"
))

assert(nrow(manifest) == 32L && sum(manifest$reuse_stage6b) == 12L &&
         sum(manifest$execute_new_fit) == 20L,
       "Manifest endpoint accounting differs from 32=12+20.")
assert(nrow(endpoints) == 32L && all(endpoints$objective_eligible) &&
         nrow(winners) == 8L && nrow(selection) == 2L &&
         nrow(science) == 8L,
       "Stage 6F endpoint, winner, selection, or evaluation counts differ.")
assert(!any(endpoints$truth_used_for_fit_stopping_or_selection) &&
         !any(selection$truth_used),
       "Truth-free selection contract failed.")

endpoint_groups <- split(
  endpoints, interaction(endpoints$split_id, endpoints$fit_config_id,
                         drop = TRUE, lex.order = TRUE)
)
endpoint_summary <- do.call(rbind, lapply(endpoint_groups, function(rows) {
  winner <- winners[
    winners$selection_stratum_id == rows$selection_stratum_id[[1L]],
    , drop = FALSE
  ]
  finite_elapsed <- as.numeric(rows$elapsed_seconds[
    is.finite(as.numeric(rows$elapsed_seconds))
  ])
  data.frame(
    split_id = rows$split_id[[1L]],
    fit_config_id = rows$fit_config_id[[1L]], endpoints = nrow(rows),
    objective_eligible = sum(rows$objective_eligible),
    practical_converged = sum(rows$practical_converged),
    slow_case = sum(rows$convergence_status == "slow_case"),
    warnings = sum(as.integer(rows$warning_count)),
    final_elbo_min = min(rows$final_elbo),
    final_elbo_max = max(rows$final_elbo),
    final_elbo_spread = max(rows$final_elbo) - min(rows$final_elbo),
    winner_fit_seed = winner$fit_seed[[1L]],
    winner_final_elbo = winner$final_elbo[[1L]],
    winner_practical_converged = winner$practical_converged[[1L]],
    new_runtime_count = length(finite_elapsed),
    new_runtime_median_seconds = if (length(finite_elapsed)) {
      stats::median(finite_elapsed)
    } else NA_real_,
    new_runtime_max_seconds = if (length(finite_elapsed)) {
      max(finite_elapsed)
    } else NA_real_, stringsAsFactors = FALSE
  )
}))
rownames(endpoint_summary) <- NULL
endpoint_summary <- endpoint_summary[order(endpoint_summary$split_id,
                                           endpoint_summary$fit_config_id), ]
write_csv(endpoint_summary, "ENDPOINT_OPERATIONAL_SUMMARY.csv")

ppi_rows <- lapply(seq_len(nrow(manifest)), function(index) {
  row <- manifest[index, , drop = FALSE]
  fit <- readRDS(file.path(output_root, "fits", row$fit_id[[1L]], "fit.rds"))
  shared <- as.numeric(unlist(fit$factor_ppi_shared, use.names = FALSE))
  specific_1 <- as.numeric(unlist(fit$factor_ppi_specific[[1L]],
                                  use.names = FALSE))
  specific_2 <- as.numeric(unlist(fit$factor_ppi_specific[[2L]],
                                  use.names = FALSE))
  result <- data.frame(
    fit_id = row$fit_id[[1L]], split_id = row$split_id[[1L]],
    fit_config_id = row$fit_config_id[[1L]], fit_seed = row$fit_seed[[1L]],
    selected_shared = sum(shared >= 0.5),
    selected_specific_1 = sum(specific_1 >= 0.5),
    selected_specific_2 = sum(specific_2 >= 0.5),
    shared_ppi = paste(format(shared, digits = 10), collapse = ";"),
    specific_1_ppi = paste(format(specific_1, digits = 10), collapse = ";"),
    specific_2_ppi = paste(format(specific_2, digits = 10), collapse = ";"),
    stringsAsFactors = FALSE
  )
  rm(fit)
  gc(FALSE)
  result
})
ppi_rows <- do.call(rbind, ppi_rows)
ppi_rows$selected_counts <- paste(
  ppi_rows$selected_shared, ppi_rows$selected_specific_1,
  ppi_rows$selected_specific_2, sep = "/"
)
ppi_rows$nominal_count_correct <-
  ppi_rows$selected_shared == 1L & ppi_rows$selected_specific_1 == 1L &
  ppi_rows$selected_specific_2 == 1L
write_csv(ppi_rows, "ALL_ENDPOINT_PPI_SELECTION.csv")
ppi_summary <- do.call(rbind, lapply(split(
  ppi_rows, interaction(ppi_rows$split_id, ppi_rows$fit_config_id,
                        drop = TRUE, lex.order = TRUE)
), function(rows) {
  data.frame(
    split_id = rows$split_id[[1L]],
    fit_config_id = rows$fit_config_id[[1L]], endpoints = nrow(rows),
    nominal_1_1_1_count = sum(rows$nominal_count_correct),
    selected_counts_observed = paste(sort(unique(rows$selected_counts)),
                                     collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
rownames(ppi_summary) <- NULL
write_csv(ppi_summary, "PPI_SELECTION_SUMMARY.csv")

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
logscore_tool <- file.path(
  dirname(script_root),
  "g12-predictive-logscore-feasibility-stage6e-b-20260910",
  "stage6eb_logscore_tools_20260910_v1.R"
)
source(logscore_tool, local = TRUE)

decomposition <- lapply(seq_len(nrow(winners)), function(index) {
  winner <- winners[index, , drop = FALSE]
  fit <- readRDS(file.path(output_root, "fits", winner$fit_id[[1L]],
                           "fit.rds"))
  holdout <- readRDS(file.path(output_root, "splits",
                               winner$split_id[[1L]], "holdout_plan.rds"))
  rows <- lapply(holdout$records, function(record) {
    moments <- s6eb_predictive_moments(
      fit, record$study, record$subject, record$time
    )
    observed <- as.numeric(record$y)
    score <- s6eb_lowrank_gaussian_nlpd(
      observed, moments$mean, moments$diagonal_variance, moments$low_rank
    )
    diagonal_used <- pmax(moments$diagonal_variance, 1e-10)
    residual <- observed - moments$mean
    diagonal_quadratic <- sum(residual^2 / diagonal_used)
    lowrank_variance <- rowSums(moments$low_rank^2)
    data.frame(
      study = record$study, subject = record$subject,
      joint_nlpd_per_variable = score$nlpd_per_variable,
      gaussian_constant_component = 0.5 * log(2 * pi),
      logdet_component = 0.5 * score$log_determinant / length(observed),
      quadratic_component = 0.5 * score$quadratic / length(observed),
      lowrank_logdet_penalty =
        0.5 * score$reduced_log_determinant / length(observed),
      lowrank_quadratic_benefit =
        0.5 * (diagonal_quadratic - score$quadratic) / length(observed),
      diagonal_only_nlpd = 0.5 * (
        log(2 * pi) + mean(log(diagonal_used)) +
          diagonal_quadratic / length(observed)
      ),
      covariance_rank_update = score$covariance_rank_update,
      mean_diagonal_variance = mean(diagonal_used),
      mean_lowrank_variance = mean(lowrank_variance),
      lowrank_fraction_of_marginal_variance = mean(
        lowrank_variance / (diagonal_used + lowrank_variance)
      ), stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  numeric_names <- names(rows)[vapply(rows, is.numeric, logical(1L))]
  numeric_names <- setdiff(numeric_names, c("study", "subject"))
  result <- data.frame(
    split_id = winner$split_id[[1L]],
    fit_config_id = winner$fit_config_id[[1L]],
    fit_seed = winner$fit_seed[[1L]], subject_time_vectors = nrow(rows),
    stringsAsFactors = FALSE
  )
  for (name in numeric_names) result[[name]] <- mean(rows[[name]])
  rm(fit, holdout, rows)
  gc(FALSE)
  result
})
decomposition <- do.call(rbind, decomposition)
write_csv(decomposition, "JOINT_LOGSCORE_DECOMPOSITION_BY_SPLIT.csv")
decomposition_config <- do.call(rbind, lapply(
  split(decomposition, decomposition$fit_config_id), function(rows) {
    result <- data.frame(fit_config_id = rows$fit_config_id[[1L]],
                         split_winners = nrow(rows), stringsAsFactors = FALSE)
    numeric_names <- names(rows)[vapply(rows, is.numeric, logical(1L))]
    numeric_names <- setdiff(numeric_names,
                             c("fit_seed", "subject_time_vectors"))
    for (name in numeric_names) result[[name]] <- mean(rows[[name]])
    result
  }
))
rownames(decomposition_config) <- NULL
write_csv(decomposition_config, "JOINT_LOGSCORE_DECOMPOSITION_BY_CONFIG.csv")

science$exact_structure_recovery <- with(science,
  global_correct == 3L & global_misplaced == 0L & global_missing == 0L &
    global_extra == 0L & global_duplicate == 0L & study_correct == 4L &
    study_misplaced == 0L & study_missing == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L
)
science$true_core_complete <- with(science,
  global_correct == 3L & global_misplaced == 0L & global_missing == 0L &
    study_correct == 4L & study_misplaced == 0L & study_missing == 0L
)
science$structural_error_total <- with(science,
  global_misplaced + global_missing + global_extra + global_duplicate +
    study_misplaced + study_missing + study_candidate_extra +
    study_candidate_duplicate
)
selected_science_columns <- c(
  "split_id", "fit_config_id", "fit_seed", "practical_converged",
  "exact_structure_recovery", "true_core_complete",
  "structural_error_total", "global_correct", "global_misplaced",
  "global_missing", "global_extra", "global_duplicate", "study_correct",
  "study_misplaced", "study_missing", "study_candidate_extra",
  "study_candidate_duplicate", "signal_nrmse_all_candidates",
  "dense_signal_nrmse_all_candidates", "feature_components_total",
  "feature_components_matched", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "loading_support_ppi_auc_mean", "factor_process_nrmse_mean_matched",
  "complete_contribution_nrmse_mean_matched"
)
write_csv(science[selected_science_columns],
          "SCIENTIFIC_RESULTS_BY_SPLIT_CONFIG.csv")

criterion_audit <- merge(selection, calibration,
                         by.x = "selected_config_id",
                         by.y = "fit_config_id", all.x = TRUE, sort = FALSE)
criterion_audit <- merge(criterion_audit, science_by_config,
                         by.x = "selected_config_id",
                         by.y = "fit_config_id", all.x = TRUE, sort = FALSE)
criterion_audit$selected_generating_dimensions <-
  criterion_audit$selected_config_id == "joint_L1_M2"
write_csv(criterion_audit, "CRITERION_DECISION_AUDIT.csv")

prepared_time <- marker_value(file.path(output_root, "STAGE6F_PREPARED.txt"),
                              "prepared_utc")
completed_time <- marker_value(file.path(output_root, "truth_evaluation",
                                         "EVALUATION_COMPLETE.txt"),
                               "completed_utc")
elapsed_new <- as.numeric(endpoints$elapsed_seconds[
  !endpoints$imported_from_stage6b &
    is.finite(as.numeric(endpoints$elapsed_seconds))
])
run_summary <- data.frame(
  registered_endpoints = nrow(manifest), reused_endpoints =
    sum(manifest$reuse_stage6b), new_fits = sum(manifest$execute_new_fit),
  fit_errors = length(fit_errors), score_errors = length(score_errors),
  evaluation_errors = length(evaluation_errors),
  objective_eligible = sum(endpoints$objective_eligible),
  practical_converged = sum(endpoints$practical_converged),
  warning_total = sum(as.integer(endpoints$warning_count)),
  new_fit_cpu_hours = sum(elapsed_new) / 3600,
  new_fit_median_minutes = stats::median(elapsed_new) / 60,
  new_fit_max_minutes = max(elapsed_new) / 60,
  prepared_utc = prepared_time, evaluation_completed_utc = completed_time,
  wall_clock_hours = as.numeric(difftime(
    as.POSIXct(completed_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    as.POSIXct(prepared_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    units = "hours"
  )), stringsAsFactors = FALSE
)
write_csv(run_summary, "RUN_SUMMARY.csv")
write_csv(mse_candidates, "MSE_ONE_SE_CANDIDATES_COPY.csv")
write_csv(logscore_candidates, "LOGSCORE_ONE_SE_CANDIDATES_COPY.csv")
write_csv(rank_stability, "CROSS_SPLIT_RANK_STABILITY_COPY.csv")
write_csv(split_risks, "SPLIT_SPECIFIC_RISKS_COPY.csv")

integrity <- data.frame(
  check_id = c(
    "endpoint_accounting_32_12_20", "all_objective_eligible",
    "no_fit_score_evaluation_errors", "eight_winners_evaluated",
    "truth_free_selection_precedes_evaluation", "no_continuation",
    "joint_logscore_decomposition_identity"
  ),
  passed = c(
    nrow(manifest) == 32L && sum(manifest$reuse_stage6b) == 12L &&
      sum(manifest$execute_new_fit) == 20L,
    all(endpoints$objective_eligible),
    !length(c(fit_errors, score_errors, evaluation_errors)),
    nrow(winners) == 8L && nrow(science) == 8L,
    !any(selection$truth_used),
    !any(endpoints$continuation_used),
    max(abs(decomposition$joint_nlpd_per_variable -
              (decomposition$gaussian_constant_component +
                 decomposition$logdet_component +
                 decomposition$quadratic_component))) < 1e-10
  ), stringsAsFactors = FALSE
)
write_csv(integrity, "ANALYSIS_INTEGRITY.csv")
writeLines(c(
  paste0("status=", if (all(integrity$passed)) "PASS" else "FAIL"),
  paste0("checks_passed=", sum(integrity$passed)),
  paste0("checks_total=", nrow(integrity)),
  "new_fits_started_by_analysis=0", "continuations_started_by_analysis=0",
  "selection_recomputed=FALSE"
), file.path(analysis_root, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
if (!all(integrity$passed)) stop("Stage 6F analysis integrity failed.",
                                 call. = FALSE)
cat("STAGE6F_ANALYSIS_COMPLETE checks=", nrow(integrity), "\n", sep = "")
