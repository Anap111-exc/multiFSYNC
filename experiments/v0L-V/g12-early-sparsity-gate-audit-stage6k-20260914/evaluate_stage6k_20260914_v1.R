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

PP_EXPERIMENT_ID <- "STAGE6K_EARLY_SPARSITY_GATE_CAUSAL_AUDIT_V1_20260914"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6K_OUTPUT_ROOT")
workers <- as.integer(cli$workers %||% pp_env("STAGE6K_EVAL_WORKERS", "4"))
stage6g_output_root <- pp_env(
  "STAGE6K_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6j_output_root <- pp_env(
  "STAGE6K_STAGE6J_OUTPUT_ROOT",
  "/root/multifsync-results/stage6j_20260914_v1"
)
stage6k_library <- pp_env("STAGE6K_R_LIB")
if (nzchar(stage6k_library)) Sys.setenv(STAGE6A_R_LIB = stage6k_library)
pp_assert(nzchar(output_root) && dir.exists(output_root) &&
            .Platform$OS.type == "unix" && workers %in% 1:4,
          "Stage 6K evaluation requires a prepared Unix output and 1--4 workers.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

authorization <- readLines(file.path(
  script_root, "STAGE6K_EXECUTION_AUTHORIZATION_20260914.txt"
), warn = FALSE)
pp_assert(all(c(
  "status=AUTHORIZED_TRUTH_FREE_STAGE6K_EXECUTION",
  "post_freeze_scientific_evaluation_authorized=TRUE",
  "new_data_authorized=FALSE", "continuation_authorized=FALSE",
  "additional_gate_lengths_authorized=FALSE"
) %in% authorization), "Stage 6K evaluation authorization is incomplete.")
pp_assert(file.exists(file.path(output_root,
                                "STAGE6K_TRUTH_FREE_COMPLETE.txt")),
          "Stage 6K truth-free endpoints must be frozen before evaluation.")

manifest <- pp_read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
terminals <- pp_read_csv(file.path(output_root, "TRUTH_FREE_ENDPOINTS.csv"))
selection <- pp_read_csv(file.path(output_root, "TRUTH_FREE_SELECTION.csv"))
pp_assert(nrow(manifest) == 4L && nrow(terminals) == 4L &&
            nrow(selection) == 1L &&
            setequal(manifest$fit_id, terminals$fit_id) &&
            all(terminals$objective_eligible) &&
            all(terminals$dense_gate_sweeps == 10L) &&
            all(terminals$dense_gate_released) &&
            !any(manifest$truth_available_to_fit) &&
            !any(manifest$truth_available_to_stopping) &&
            !any(manifest$truth_available_to_selection),
          "Stage 6K frozen endpoint accounting failed.")

source_binding <- pp_read_csv(file.path(output_root,
                                        "SOURCE_ARTIFACT_BINDING.csv"))
source_info <- file.info(source_binding$path)
source_match <- file.exists(source_binding$path) &
  as.numeric(source_info$size) == as.numeric(source_binding$bytes) &
  vapply(source_binding$path, pp_sha256, character(1L)) == source_binding$sha256
pp_assert(all(source_match), "A read-only Stage 6K source artifact changed.")

truth_path <- file.path(
  stage6g_output_root, "data", "s6g_minimal_L111_M2", "sealed_truth",
  "truth_bundle.rds"
)
control_science_path <- file.path(
  stage6j_output_root, "scientific_evaluation_20260914_v1",
  "ALL_STAGE6J_SCIENTIFIC_RESULTS.csv"
)
pp_assert(file.exists(truth_path) && file.exists(control_science_path),
          "Truth or Stage 6J control science is missing.")
truth <- readRDS(truth_path)
pp_assert(identical(truth$data_id, "s6g_minimal_L111_M2") &&
            identical(as.integer(truth$data$true_params$L_f), 1L) &&
            all(as.integer(unlist(truth$data$true_params$L_s_by_study,
                                  use.names = FALSE)) == 1L) &&
            all(as.integer(truth$data$true_params$M_f) == 2L),
          "Stage 6K truth identity failed.")

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
evaluation_root <- file.path(output_root, "scientific_evaluation_20260914_v1")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6K scientific evaluation output already exists.")
dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")
pp_write_csv(data.frame(
  role = c("unsealed_truth", "stage6j_control_science"),
  path = c(truth_path, control_science_path),
  bytes = as.numeric(file.info(c(truth_path, control_science_path))$size),
  sha256 = vapply(c(truth_path, control_science_path), pp_sha256, character(1L)),
  endpoints_frozen_before_truth_read = TRUE,
  truth_used_for_fit_stopping_or_selection = FALSE,
  stringsAsFactors = FALSE
), file.path(evaluation_root, "EVALUATION_SOURCE_BINDING.csv"))

s6k_add_recovery_flags <- function(value) {
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
    ), compact = NULL, trace = NULL))
  }
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(result_dir, paste0(name, ".csv")))
  }
  terminal <- terminals[terminals$fit_id == fit_id, , drop = FALSE]
  fit <- readRDS(fit_path)
  base <- data.frame(
    fit_id = fit_id,
    paired_control_fit_id = row$paired_control_fit_id[[1L]],
    initialization_arm = row$initialization_arm[[1L]],
    seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
    dense_gate_sweeps = 10L,
    final_elbo = terminal$final_ordinary_elbo[[1L]],
    practical_converged = terminal$practical_converged[[1L]],
    convergence_status = terminal$convergence_status[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    fit_warning_count = terminal$warning_count[[1L]],
    truth_free_winner = fit_id == selection$fit_id[[1L]],
    truth_used_for_fit_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE
  )
  activity <- result$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  compact <- cbind(base, result$summary,
                   result$process_contribution_summary,
                   result$loading_scale_summary, activity)
  trace <- fit$driver_trace
  trace <- trace[trace$iteration %in% c(
    1L, 2L, 5L, 10L, 11L, 12L, 25L, 50L, 75L, 99L,
    100L, 120L, 199L, 299L, 399L, 499L
  ), , drop = FALSE]
  trace$fit_id <- fit_id
  trace$arm <- "gate10"
  trace$seed_index <- row$seed_index[[1L]]
  trace <- trace[, c("fit_id", "arm", "seed_index",
                     setdiff(names(trace), c("fit_id", "arm", "seed_index"))),
                 drop = FALSE]
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "new_fit_started_by_evaluation=FALSE", "continuation=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(status = data.frame(
    fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed, error = "",
    stringsAsFactors = FALSE
  ), compact = compact, trace = trace)
}

results <- parallel::mclapply(
  seq_len(nrow(manifest)), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok), "Stage 6K evaluation failed; no selective rerun.")
gate_science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
gate_science <- gate_science[match(manifest$fit_id, gate_science$fit_id),
                             , drop = FALSE]
gate_science <- s6k_add_recovery_flags(gate_science)
pp_write_csv(gate_science, file.path(
  evaluation_root, "ALL_STAGE6K_GATE10_SCIENTIFIC_RESULTS.csv"
))
pp_write_csv(pp_rbind_fill(lapply(results, `[[`, "trace")), file.path(
  evaluation_root, "GATE10_TRACE_KEY_CHECKPOINTS.csv"
))

control <- pp_read_csv(control_science_path)
control <- control[
  control$initialization_arm == "random_gram_unit_energy" &
    control$seed_index %in% 1:4, , drop = FALSE
]
pp_assert(nrow(control) == 4L && !anyDuplicated(control$seed_index) &&
            sum(control$core_truth_recovered) == 0L,
          "Stage 6J paired G12 control identity failed.")

pair_columns <- c(
  "seed_index", "fit_id", "final_elbo", "core_truth_recovered",
  "exact_thresholded_factor_structure", "ppi_factor_counts_exact",
  "global_correct", "global_misplaced", "global_missing",
  "study_correct", "study_misplaced", "study_missing",
  "selected_shared", "selected_specific_study_1",
  "selected_specific_study_2", "signal_nrmse_ppi_selected",
  "dense_signal_nrmse_ppi_selected",
  "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean", "matched_function_recall_mean",
  "correctly_placed_trajectory_abs_cor_mean",
  "correctly_placed_score_abs_cor_mean",
  "complete_contribution_nrmse_mean_missing_as_zero"
)
pair_columns <- Reduce(intersect, list(pair_columns,
                                       names(control), names(gate_science)))
paired <- merge(
  control[, pair_columns, drop = FALSE],
  gate_science[, pair_columns, drop = FALSE],
  by = "seed_index", suffixes = c("_control", "_gate10"), sort = TRUE
)
paired$delta_elbo_gate_minus_control <-
  paired$final_elbo_gate10 - paired$final_elbo_control
for (metric in intersect(c(
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "complete_contribution_nrmse_mean_missing_as_zero"
), pair_columns)) {
  paired[[paste0("delta_", metric, "_gate_minus_control")]] <-
    paired[[paste0(metric, "_gate10")]] -
    paired[[paste0(metric, "_control")]]
}
pp_write_csv(paired, file.path(
  evaluation_root, "PAIRED_GATE10_VS_STAGE6J_CONTROL.csv"
))

gate_hits <- sum(gate_science$core_truth_recovered)
correct_rows <- gate_science$core_truth_recovered
correct_elbo_not_lower <- gate_hits > 0L && all(
  paired$delta_elbo_gate_minus_control[
    match(gate_science$seed_index[correct_rows], paired$seed_index)
  ] >= 0
)
winner_core <- gate_science$core_truth_recovered[
  gate_science$truth_free_winner
][[1L]]
decision <- if (gate_hits >= 2L && correct_elbo_not_lower) {
  "CAUSAL_SUPPORT_EARLY_SPARSITY_FEEDBACK_MATERIAL_BARRIER_NO_DEFAULT_CHANGE"
} else if (gate_hits >= 2L) {
  "REACHABILITY_RESCUE_WITH_OBJECTIVE_CONFLICT_NO_DEFAULT_CHANGE"
} else if (gate_hits == 1L) {
  "WEAK_UNSTABLE_GATE_SIGNAL_NO_DEFAULT_CHANGE"
} else {
  "ZERO_OF_FOUR_CLOSE_DENSE_GATE_ROUTE"
}
decision_table <- data.frame(
  paired_controls = 4L, control_core_hits = 0L,
  gate10_endpoints = 4L, gate10_core_hits = gate_hits,
  correct_gate_endpoint_elbo_not_lower_than_paired_control =
    correct_elbo_not_lower,
  truth_free_gate10_winner_fit_id = selection$fit_id[[1L]],
  truth_free_gate10_winner_core = winner_core,
  additional_gate_lengths_authorized = FALSE,
  decision = decision, formal_paper_mc = FALSE,
  stringsAsFactors = FALSE
)
pp_write_csv(decision_table,
             file.path(evaluation_root, "STAGE6K_DECISION.csv"))
pp_write_csv(gate_science[gate_science$truth_free_winner, , drop = FALSE],
             file.path(evaluation_root, "TRUTH_FREE_GATE10_WINNER_SCIENCE.csv"))
pp_write_lines(c(
  "status=STAGE6K_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("completed_utc=", pp_iso_time()),
  "new_fits_evaluated=4", "paired_controls_readonly=4",
  "evaluation_errors=0", "truth_read_after_endpoint_freeze=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fits_started_by_evaluation=FALSE", "continuation_started=FALSE",
  "additional_gate_lengths_started=FALSE",
  paste0("gate10_core_hits=", gate_hits),
  paste0("truth_free_gate10_winner_core=", winner_core),
  paste0("decision=", decision)
), file.path(evaluation_root, "STAGE6K_EVALUATION_COMPLETE.txt"))
cat("STAGE6K_EVALUATION_COMPLETE gate10_core=", gate_hits,
    " winner_core=", winner_core, " decision=", decision, "\n", sep = "")
