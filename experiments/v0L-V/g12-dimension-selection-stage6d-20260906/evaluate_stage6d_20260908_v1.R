#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(dirname(script_root),
                 "g12-dimension-overspec-stage6a-20260831",
                 "stage6a_common_20260831_v1.R"), local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6D_INDEPENDENT_DIMENSION_SELECTION_V1_20260906"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% pp_env("STAGE6D_EVAL_WORKERS", "7"))
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify the completed Stage 6D output.")
pp_assert(.Platform$OS.type == "unix" && length(workers) == 1L &&
            is.finite(workers) && workers %in% 1:7,
          "Stage 6D evaluation requires Unix and 1--7 workers.")
stage6d_library <- pp_env("STAGE6D_R_LIB")
if (nzchar(stage6d_library)) Sys.setenv(STAGE6A_R_LIB = stage6d_library)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

selection_root <- file.path(output_root, "truth_free_selection")
selection_path <- file.path(selection_root, "DIMENSION_SELECTIONS.csv")
winner_path <- file.path(selection_root, "CONFIG_WINNERS.csv")
selection_marker <- file.path(selection_root,
                              "TRUTH_FREE_SELECTION_FROZEN.txt")
authorization_path <- file.path(output_root,
                                "TRUTH_UNSEAL_AUTHORIZATION.txt")
experiment_complete <- file.path(output_root, "STAGE6D_COMPLETE.txt")
evaluation_root <- file.path(output_root, "evaluation")

pp_assert(file.exists(selection_path) && file.exists(winner_path) &&
            file.exists(selection_marker) && file.exists(experiment_complete),
          "The complete frozen truth-free Stage 6D selection is missing.")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6D evaluation already exists; refusing overwrite.")
selection_hash <- pp_sha256(selection_path)
winner_hash <- pp_sha256(winner_path)
pp_assert(file.exists(authorization_path),
          "Explicit Stage 6D truth-unseal authorization is missing.")
authorization <- readLines(authorization_path, warn = FALSE)
required_authorization <- c(
  "status=TRUTH_UNSEAL_AND_EVALUATION_AUTHORIZED",
  paste0("dimension_selection_sha256=", selection_hash),
  paste0("configuration_winner_sha256=", winner_hash),
  "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
  "new_fits_authorized=FALSE", "continuation_authorized=FALSE",
  "full_data_refit_authorized=FALSE",
  "authorization_scope=stage6d_existing_180_fits_only"
)
pp_assert(all(required_authorization %in% authorization),
          "Truth-unseal authorization does not bind both frozen selections.")

fit_root <- file.path(output_root, "fits")
fit_dirs <- list.dirs(fit_root, recursive = FALSE, full.names = TRUE)
fit_dirs <- fit_dirs[file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))]
pp_assert(length(fit_dirs) == 180L,
          "Stage 6D truth evaluation requires exactly 180 complete fits.")
pp_assert(length(list.files(fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE,
                            full.names = TRUE)) == 0L,
          "A Stage 6D fit error is present.")

records <- lapply(fit_dirs, function(path) {
  record <- readRDS(file.path(path, "terminal_record.rds"))
  pp_assert(identical(record$terminal_schema, "G12_STAGE6D_TERMINAL_V1") &&
              isTRUE(record$objective_eligible) &&
              !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
              !isTRUE(record$truth_used_for_holdout_scoring) &&
              identical(record$terminal_status, "fixed_400_complete"),
            paste0("Ineligible Stage 6D endpoint: ", basename(path)))
  record
})
names(records) <- vapply(records, `[[`, character(1L), "fit_id")
records <- records[order(names(records))]
pp_assert(!anyDuplicated(names(records)), "Duplicate Stage 6D fit_id.")

winners <- pp_read_csv(winner_path)
selections <- pp_read_csv(selection_path)
pp_assert(nrow(winners) == 15L && all(table(winners$data_id) == 5L) &&
            nrow(selections) == 3L &&
            !any(selections$selection_used_truth),
          "Malformed frozen Stage 6D selection tables.")
winner_ids <- winners$fit_id

factor_map <- merge(
  selections[c("data_id", "selected_factor_config_id")],
  winners[c("data_id", "fit_config_id", "fit_id")],
  by.x = c("data_id", "selected_factor_config_id"),
  by.y = c("data_id", "fit_config_id"), sort = FALSE
)
fpca_map <- merge(
  selections[c("data_id", "selected_fpca_config_id")],
  winners[c("data_id", "fit_config_id", "fit_id")],
  by.x = c("data_id", "selected_fpca_config_id"),
  by.y = c("data_id", "fit_config_id"), sort = FALSE
)
pp_assert(nrow(factor_map) == 3L && nrow(fpca_map) == 3L,
          "Axis selections do not map to six frozen configuration winners.")
factor_selected_ids <- factor_map$fit_id
fpca_selected_ids <- fpca_map$fit_id

pp_load_package()
pp_source_runtime(TRUE)
truth_paths <- file.path(output_root, "data", selections$data_id,
                         "sealed_truth", "truth_bundle.rds")
pp_assert(all(file.exists(truth_paths)), "A Stage 6D truth bundle is missing.")
truths <- lapply(truth_paths, readRDS)
names(truths) <- selections$data_id
pp_assert(all(vapply(seq_along(truths), function(index) {
  truth <- truths[[index]]
  identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
    identical(truth$data_id, names(truths)[[index]]) &&
    identical(truth$experiment_id, PP_EXPERIMENT_ID)
}, logical(1L))), "Stage 6D truth identity check failed.")
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

dir.create(evaluation_root, recursive = TRUE, mode = "0700")
pp_write_lines(c(
  "status=EVALUATION_STARTED", paste0("started_utc=", pp_iso_time()),
  paste0("workers=", workers), "registered_fits=180",
  paste0("dimension_selection_sha256=", selection_hash),
  paste0("configuration_winner_sha256=", winner_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_STARTED.txt"))

parse_rank <- function(value) {
  as.integer(strsplit(as.character(value), ";", fixed = TRUE)[[1L]])
}
expected_tables <- c(
  "feature_component_errors", "factor_kernel_errors",
  "dense_signal_metrics", "loading_recovery", "threshold_sensitivity",
  "unthresholded_assignment", "covariance_operator_error",
  "ppi_diagnostics", "shared_candidate_metrics", "summary",
  "evaluation_contract", "process_contribution_by_role",
  "process_contribution_summary", "loading_scale_by_block",
  "loading_scale_summary", "candidate_activity",
  "candidate_activity_summary", "fpca_rank_activity",
  "fpca_component_activity"
)

evaluate_one <- function(index) {
  record <- records[[index]]
  fit_id <- record$fit_id
  fit_dir <- file.path(fit_root, fit_id)
  output_dir <- file.path(evaluation_root, fit_id)
  dir.create(output_dir, recursive = TRUE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    value <- pp_evaluate_fit(
      readRDS(file.path(fit_dir, "fit.rds")), record,
      truths[[record$data_id]], evaluator_environment
    )
    oracle_dimensions <- record$fit_L_f == 1L &&
      record$fit_L_s_1 == 1L && record$fit_L_s_2 == 1L &&
      all(parse_rank(record$fit_M_f) == 2L) &&
      all(parse_rank(record$fit_M_s_1) == 2L) &&
      all(parse_rank(record$fit_M_s_2) == 2L)
    value$evaluation_contract$oracle_L_M <- oracle_dimensions
    value$evaluation_contract$stage6d_truth_free_selection_frozen <- TRUE
    value$evaluation_contract$heldout_observations_excluded_from_fit <- TRUE
    value
  }, warning = function(condition) {
    warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
      class = paste(class(condition), collapse = ";"),
      message = conditionMessage(condition),
      call = paste(deparse(conditionCall(condition)), collapse = " "),
      stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  }), error = function(condition) condition)
  ended <- Sys.time()
  warning_table <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(class = character(), message = character(), call = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warning_table,
               file.path(output_dir, "EVALUATION_WARNINGS.csv"))
  elapsed <- as.numeric(difftime(ended, started, units = "secs"))
  if (inherits(result, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(result))
    ), file.path(output_dir, "EVALUATION_ERROR.txt"))
    return(list(status = data.frame(
      fit_id = fit_id, ok = FALSE, elapsed_seconds = elapsed,
      warning_count = nrow(warning_table), error = conditionMessage(result),
      stringsAsFactors = FALSE
    ), compact = NULL, metric = NULL))
  }

  pp_save_rds(result, file.path(output_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(output_dir, paste0(name, ".csv")))
  }
  science <- cbind(result$summary, result$process_contribution_summary,
                   result$loading_scale_summary)
  base <- data.frame(
    fit_id = fit_id, data_id = record$data_id,
    scenario_id = record$scenario_id, method_id = record$method_id,
    fit_config_id = record$fit_config_id,
    selection_stratum_id = record$selection_stratum_id,
    seed_index = record$seed_index, fit_seed = record$fit_seed,
    fit_L_f = record$fit_L_f, fit_L_s_1 = record$fit_L_s_1,
    fit_L_s_2 = record$fit_L_s_2, fit_M_f = record$fit_M_f,
    fit_M_s_1 = record$fit_M_s_1, fit_M_s_2 = record$fit_M_s_2,
    final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    convergence_status = record$convergence_status,
    heldout_rmse_all = record$heldout_rmse_all,
    heldout_nrmse_all = record$heldout_nrmse_all,
    within_configuration_winner = fit_id %in% winner_ids,
    factor_axis_selected_endpoint = fit_id %in% factor_selected_ids,
    fpca_axis_selected_endpoint = fit_id %in% fpca_selected_ids,
    selected_by_any_axis = fit_id %in%
      unique(c(factor_selected_ids, fpca_selected_ids)),
    stringsAsFactors = FALSE
  )
  compact <- cbind(base, science)
  tables_present <- all(expected_tables %in% names(result)) &&
    all(vapply(result[expected_tables], function(value) {
      is.data.frame(value) && nrow(value) > 0L
    }, logical(1L)))
  checkpoint_complete <- identical(
    names(readRDS(file.path(fit_dir, "fit.rds"))$practical_checkpoints),
    c("380", "400")
  )
  metric <- data.frame(
    fit_id = fit_id, expected_table_count = length(expected_tables),
    expected_tables_present = tables_present,
    checkpoint_380_400_complete = checkpoint_complete,
    metric_contract_passed = tables_present && checkpoint_complete,
    stringsAsFactors = FALSE
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("dimension_selection_sha256=", selection_hash),
    paste0("configuration_winner_sha256=", winner_hash),
    "selection_frozen_before_truth_read=TRUE",
    "fit_or_selection_used_truth=FALSE", "new_fits_after_unseal=0",
    "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(output_dir, "EVALUATION_COMPLETE.txt"))
  list(
    status = data.frame(
      fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed,
      warning_count = nrow(warning_table), error = "",
      stringsAsFactors = FALSE
    ), compact = compact, metric = metric
  )
}

results <- parallel::mclapply(
  seq_along(records), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
if (!all(status$ok)) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("failed_fit_ids=", paste(status$fit_id[!status$ok], collapse = ";"))
  ), file.path(evaluation_root, "EVALUATION_ERROR.txt"))
  pp_abort("Stage 6D evaluation errors retained; no selective rerun.")
}

metric_qc <- do.call(rbind, lapply(results, `[[`, "metric"))
pp_assert(nrow(metric_qc) == 180L && all(metric_qc$metric_contract_passed),
          "One or more Stage 6D evaluation metric contracts failed.")
all_science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
config_winners <- all_science[
  all_science$within_configuration_winner, , drop = FALSE
]
pp_assert(nrow(all_science) == 180L && nrow(config_winners) == 15L,
          "Stage 6D evaluated winner mapping failed.")

axis_map <- rbind(
  data.frame(data_id = factor_map$data_id, selection_axis = "factor_count",
             selected_config_id = factor_map$selected_factor_config_id,
             fit_id = factor_map$fit_id, stringsAsFactors = FALSE),
  data.frame(data_id = fpca_map$data_id, selection_axis = "fpca_cap",
             selected_config_id = fpca_map$selected_fpca_config_id,
             fit_id = fpca_map$fit_id, stringsAsFactors = FALSE)
)
axis_science <- merge(axis_map, all_science, by = c("data_id", "fit_id"),
                      sort = FALSE)
axis_science <- axis_science[order(axis_science$data_id,
                                   axis_science$selection_axis), , drop = FALSE]
pp_assert(nrow(axis_science) == 6L,
          "Stage 6D axis-selected endpoint mapping failed.")

pp_write_csv(metric_qc,
             file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
pp_write_csv(all_science,
             file.path(evaluation_root, "ALL_STAGE6D_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(config_winners,
             file.path(evaluation_root, "CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(axis_science,
             file.path(evaluation_root, "AXIS_SELECTED_SCIENTIFIC_RESULTS.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "fits_evaluated=180", "evaluation_errors=0",
  "configuration_winners_evaluated=15", "axis_selections_evaluated=6",
  paste0("dimension_selection_sha256=", selection_hash),
  paste0("configuration_winner_sha256=", winner_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE6D_EVALUATION_COMPLETE fits=180 errors=0 axis_selections=6\n")
