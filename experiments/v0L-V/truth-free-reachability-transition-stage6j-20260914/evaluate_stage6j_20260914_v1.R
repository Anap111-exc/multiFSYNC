#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
experiments_root <- dirname(script_root)
source(file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)
source(file.path(
  experiments_root, "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R"
), local = FALSE)

PP_EXPERIMENT_ID <- "STAGE6J_TRUTH_FREE_REACHABILITY_TRANSITION_V1_20260914"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6J_OUTPUT_ROOT")
workers <- as.integer(cli$workers %||% pp_env("STAGE6J_EVAL_WORKERS", "7"))
stage6g_output_root <- pp_env(
  "STAGE6J_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6j_library <- pp_env("STAGE6J_R_LIB")
if (nzchar(stage6j_library)) Sys.setenv(STAGE6A_R_LIB = stage6j_library)
pp_assert(nzchar(output_root) && dir.exists(output_root) &&
            .Platform$OS.type == "unix" && workers %in% 1:8,
          "Stage 6J evaluation requires a prepared Unix output and 1--8 workers.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

authorization_path <- file.path(
  script_root, "STAGE6J_TRUTH_UNSEAL_AUTHORIZATION_20260914.txt"
)
authorization <- readLines(authorization_path, warn = FALSE)
pp_assert(all(c(
  "status=AUTHORIZED_POST_FREEZE_STAGE6J_TRUTH_UNSEAL",
  "truth_unseal_authorized=TRUE",
  "scientific_evaluation_authorized=TRUE",
  "new_fit_authorized=FALSE",
  "continuation_authorized=FALSE",
  "cross_dimension_selection_authorized=FALSE"
) %in% authorization), "Stage 6J truth-unseal authorization is incomplete.")
pp_assert(file.exists(file.path(output_root,
                                "STAGE6J_TRUTH_FREE_COMPLETE.txt")),
          "Stage 6J truth-free endpoints must be frozen before evaluation.")

manifest <- pp_read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
terminals <- pp_read_csv(file.path(output_root, "TRUTH_FREE_ENDPOINTS.csv"))
selection <- pp_read_csv(file.path(output_root, "TRUTH_FREE_SELECTION.csv"))
pp_assert(nrow(manifest) == 8L && nrow(terminals) == 8L &&
            nrow(selection) == 1L &&
            setequal(manifest$fit_id, terminals$fit_id) &&
            all(terminals$objective_eligible) &&
            !any(manifest$truth_available_to_fit) &&
            !any(manifest$truth_available_to_stopping) &&
            !any(manifest$truth_available_to_selection),
          "Stage 6J frozen endpoint accounting failed.")

source_binding <- pp_read_csv(file.path(output_root,
                                        "SOURCE_ARTIFACT_BINDING.csv"))
source_info <- file.info(source_binding$path)
source_match <- file.exists(source_binding$path) &
  as.numeric(source_info$size) == as.numeric(source_binding$bytes) &
  vapply(source_binding$path, pp_sha256, character(1L)) == source_binding$sha256
pp_assert(all(source_match), "A read-only Stage 6J source artifact changed.")

truth_path <- file.path(
  stage6g_output_root, "data", "s6g_minimal_L111_M2", "sealed_truth",
  "truth_bundle.rds"
)
pp_assert(file.exists(truth_path), "The sealed Stage 6G truth bundle is missing.")
truth <- readRDS(truth_path)
pp_assert(identical(truth$data_id, "s6g_minimal_L111_M2") &&
            identical(as.integer(truth$data$true_params$L_f), 1L) &&
            all(as.integer(unlist(truth$data$true_params$L_s_by_study,
                                  use.names = FALSE)) == 1L) &&
            all(as.integer(truth$data$true_params$M_f) == 2L),
          "Stage 6J truth identity failed.")

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
evaluation_root <- file.path(output_root,
                             "scientific_evaluation_20260914_v1")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6J scientific evaluation output already exists.")
dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")
pp_write_csv(data.frame(
  role = "unsealed_truth", path = truth_path,
  bytes = as.numeric(file.info(truth_path)$size),
  sha256 = pp_sha256(truth_path),
  endpoints_frozen_before_unseal = TRUE,
  truth_used_for_fit_stopping_or_selection = FALSE,
  authorization_path = authorization_path,
  authorization_sha256 = pp_sha256(authorization_path),
  stringsAsFactors = FALSE
), file.path(evaluation_root, "UNSEALED_TRUTH_BINDING.csv"))

s6j_add_recovery_flags <- function(value) {
  value$ppi_factor_counts_exact <- with(
    value, selected_shared == truth_L_f &
      selected_specific_study_1 == truth_L_s_1 &
      selected_specific_study_2 == truth_L_s_2
  )
  value$core_truth_recovered <- with(
    value, global_correct == truth_global_direction_count &
      global_misplaced == 0L & global_missing == 0L &
      study_correct == truth_study_role_count &
      study_misplaced == 0L & study_missing == 0L
  )
  value$exact_thresholded_factor_structure <- with(
    value, core_truth_recovered & ppi_factor_counts_exact &
      global_extra == 0L & global_duplicate == 0L &
      study_candidate_extra == 0L & study_candidate_duplicate == 0L
  )
  value$full_true_feature_component_coverage <- with(
    value, feature_components_matched == feature_components_total
  )
  value$core_with_full_feature_coverage <- with(
    value, core_truth_recovered & full_true_feature_component_coverage
  )
  value
}

evaluate_one <- function(index) {
  row <- manifest[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  fit_path <- file.path(output_root, "fits", fit_id, "fit.rds")
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  warning_rows <- list()
  started <- Sys.time()
  result <- tryCatch(withCallingHandlers({
    s6h_evaluate_fit(readRDS(fit_path), truth, evaluator_environment)
  }, warning = function(condition) {
    warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
      class = paste(class(condition), collapse = ";"),
      message = conditionMessage(condition), stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  }), error = function(condition) condition)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(class = character(), message = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  if (inherits(result, "condition")) {
    pp_write_lines(c("status=ERROR", paste0("message=", conditionMessage(result))),
                   file.path(result_dir, "EVALUATION_ERROR.txt"))
    return(list(status = data.frame(
      fit_id = fit_id, ok = FALSE, elapsed_seconds = elapsed,
      error = conditionMessage(result), stringsAsFactors = FALSE
    ), compact = NULL, trace = NULL, provenance = NULL))
  }
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(result_dir, paste0(name, ".csv")))
  }
  terminal <- terminals[terminals$fit_id == fit_id, , drop = FALSE]
  fit <- readRDS(fit_path)
  base <- data.frame(
    fit_id = fit_id,
    initialization_arm = row$initialization_arm[[1L]],
    seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
    final_elbo = terminal$final_ordinary_elbo[[1L]],
    practical_converged = terminal$practical_converged[[1L]],
    convergence_status = terminal$convergence_status[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    fit_warning_count = terminal$warning_count[[1L]],
    truth_free_winner = fit_id == selection$fit_id[[1L]],
    truth_used_for_fit = FALSE,
    truth_used_for_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE
  )
  activity <- result$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  compact <- cbind(base, result$summary,
                   result$process_contribution_summary,
                   result$loading_scale_summary, activity)
  trace <- fit$driver_trace
  trace <- trace[trace$iteration %in%
                   c(1L, 2L, 5L, 10L, 25L, 50L, 75L, 99L, 100L, 101L,
                     105L, 110L, 120L, 129L, 149L, 199L, 299L, 399L, 499L),
                 , drop = FALSE]
  trace$fit_id <- fit_id
  trace$initialization_arm <- row$initialization_arm[[1L]]
  trace$seed_index <- row$seed_index[[1L]]
  trace <- trace[, c("fit_id", "initialization_arm", "seed_index",
                     setdiff(names(trace),
                             c("fit_id", "initialization_arm", "seed_index"))),
                 drop = FALSE]
  init <- fit$initialization_diagnostics
  provenance <- data.frame(
    fit_id = fit_id,
    initialization_arm = row$initialization_arm[[1L]],
    seed_index = row$seed_index[[1L]],
    requested_method = init$requested_method %||% NA_character_,
    used_method = init$used_method %||% NA_character_,
    initialization_complete = isTRUE(init$complete),
    fallback_count = as.integer(init$fallback_count %||% NA_integer_),
    explained_fraction = as.numeric(init$explained_fraction %||% NA_real_),
    initialization_elapsed_seconds = as.numeric(init$elapsed_seconds %||% NA_real_),
    independence_id = fit$initialization_independence_id %||% NA_character_,
    fit_sha256 = pp_sha256(fit_path), stringsAsFactors = FALSE
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "cross_dimension_selection=FALSE", "continuation=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(status = data.frame(
    fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed, error = "",
    stringsAsFactors = FALSE
  ), compact = compact, trace = trace, provenance = provenance)
}

results <- parallel::mclapply(
  seq_len(nrow(manifest)), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok), "Stage 6J evaluation failed; no selective rerun.")
science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
science <- science[match(manifest$fit_id, science$fit_id), , drop = FALSE]
science <- s6j_add_recovery_flags(science)
pp_write_csv(science, file.path(evaluation_root,
                                "ALL_STAGE6J_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(pp_rbind_fill(lapply(results, `[[`, "trace")),
             file.path(evaluation_root, "TRACE_KEY_CHECKPOINTS.csv"))
provenance <- pp_rbind_fill(lapply(results, `[[`, "provenance"))
pp_write_csv(provenance,
             file.path(evaluation_root, "INITIALIZATION_PROVENANCE.csv"))

finite_mean <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
finite_median <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) median(value) else NA_real_
}
summarize_arm <- function(arm) {
  value <- science[science$initialization_arm == arm, , drop = FALSE]
  data.frame(
    initialization_arm = arm, endpoints = nrow(value),
    core_hits = sum(value$core_truth_recovered),
    exact_structure_hits = sum(value$exact_thresholded_factor_structure),
    full_feature_hits = sum(value$full_true_feature_component_coverage),
    core_with_full_feature_hits = sum(value$core_with_full_feature_coverage),
    truth_free_winner_hits = sum(value$truth_free_winner &
                                   value$core_truth_recovered),
    max_elbo = max(value$final_elbo),
    distinct_exact_elbo_count = length(unique(value$final_elbo)),
    mean_observed_signal_nrmse = finite_mean(value$signal_nrmse_ppi_selected),
    mean_dense_signal_nrmse = finite_mean(value$dense_signal_nrmse_ppi_selected),
    mean_covariance_error = finite_mean(
      value$covariance_operator_relative_error_mean
    ),
    mean_loading_error = finite_mean(value$loading_relative_l2_error_mean),
    mean_matched_feature_recall = finite_mean(value$matched_function_recall_mean),
    median_elapsed_seconds = finite_median(value$fit_elapsed_seconds),
    stringsAsFactors = FALSE
  )
}
arm_summary <- do.call(rbind, lapply(unique(science$initialization_arm),
                                    summarize_arm))
pp_write_csv(arm_summary,
             file.path(evaluation_root, "ARM_SCIENTIFIC_SUMMARY.csv"))

pair_columns <- c(
  "seed_index", "fit_id", "final_elbo", "core_truth_recovered",
  "exact_thresholded_factor_structure", "full_true_feature_component_coverage",
  "global_correct", "global_misplaced", "global_missing", "study_correct",
  "study_misplaced", "study_missing", "signal_nrmse_ppi_selected",
  "dense_signal_nrmse_ppi_selected", "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean", "matched_function_recall_mean",
  "correctly_placed_trajectory_abs_cor_mean",
  "correctly_placed_score_abs_cor_mean",
  "complete_contribution_nrmse_mean_missing_as_zero"
)
pair_columns <- intersect(pair_columns, names(science))
g12 <- science[science$initialization_arm == "random_gram_unit_energy",
               pair_columns, drop = FALSE]
residual <- science[science$initialization_arm == "residual_fpca",
                    pair_columns, drop = FALSE]
paired <- merge(g12, residual, by = "seed_index", all = TRUE,
                suffixes = c("_g12", "_residual_fpca"), sort = TRUE)
pp_write_csv(paired,
             file.path(evaluation_root, "PAIRED_SCIENTIFIC_RESULTS.csv"))

elbo_science <- science[, c(
  "fit_id", "initialization_arm", "seed_index", "final_elbo",
  "truth_free_winner", "core_truth_recovered",
  "exact_thresholded_factor_structure", "full_true_feature_component_coverage",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean"
), drop = FALSE]
elbo_science$elbo_rank_all <- rank(-elbo_science$final_elbo,
                                   ties.method = "min")
elbo_science <- elbo_science[order(-elbo_science$final_elbo), , drop = FALSE]
pp_write_csv(elbo_science,
             file.path(evaluation_root, "ELBO_SCIENCE_ASSOCIATION.csv"))
pp_write_csv(science[science$truth_free_winner, , drop = FALSE],
             file.path(evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))

g12_hits <- sum(science$core_truth_recovered[
  science$initialization_arm == "random_gram_unit_energy"
])
residual_hits <- sum(science$core_truth_recovered[
  science$initialization_arm == "residual_fpca"
])
winner_core <- science$core_truth_recovered[science$truth_free_winner][[1L]]
max_correct_elbo <- if (any(science$core_truth_recovered)) {
  max(science$final_elbo[science$core_truth_recovered])
} else NA_real_
max_wrong_elbo <- if (any(!science$core_truth_recovered)) {
  max(science$final_elbo[!science$core_truth_recovered])
} else NA_real_
decision <- if (residual_hits >= 2L &&
                (is.na(max_wrong_elbo) || max_correct_elbo >= max_wrong_elbo)) {
  "RESIDUAL_FPCA_CANDIDATE_FOR_NEW_DATA_CONFIRMATION"
} else if (g12_hits > 0L && residual_hits > 0L && !winner_core) {
  "CORRECT_ENDPOINTS_EXIST_BUT_ELBO_SELECTED_WRONG"
} else if (g12_hits == 0L && residual_hits == 0L) {
  "BOTH_TRUTH_FREE_INITIALIZATION_ARMS_FAILED_CORE_REACHABILITY"
} else if (residual_hits > 0L && residual_hits < 2L) {
  "RESIDUAL_FPCA_SIGNAL_INSUFFICIENT_FOR_CONFIRMATION"
} else {
  "MIXED_REACHABILITY_EVIDENCE_REQUIRES_TRAJECTORY_INTERPRETATION"
}
pp_write_csv(data.frame(
  g12_core_hits = g12_hits, residual_fpca_core_hits = residual_hits,
  truth_free_winner_fit_id = selection$fit_id[[1L]],
  truth_free_winner_core = winner_core,
  max_correct_endpoint_elbo = max_correct_elbo,
  max_wrong_endpoint_elbo = max_wrong_elbo,
  decision = decision,
  cross_dimension_criterion_tested = FALSE,
  formal_paper_mc = FALSE, stringsAsFactors = FALSE
), file.path(evaluation_root, "STAGE6J_DECISION.csv"))

pp_write_lines(c(
  "status=STAGE6J_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("completed_utc=", pp_iso_time()),
  "fits_evaluated=8", "evaluation_errors=0",
  "truth_unsealed_after_endpoint_freeze=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fits_started_by_evaluation=FALSE",
  "continuation_started=FALSE", "cross_dimension_selection=FALSE",
  paste0("g12_core_hits=", g12_hits),
  paste0("residual_fpca_core_hits=", residual_hits),
  paste0("truth_free_winner_core=", winner_core),
  paste0("decision=", decision)
), file.path(evaluation_root, "STAGE6J_EVALUATION_COMPLETE.txt"))
cat("STAGE6J_EVALUATION_COMPLETE g12_core=", g12_hits,
    " residual_core=", residual_hits, " winner_core=", winner_core,
    " decision=", decision, "\n", sep = "")
