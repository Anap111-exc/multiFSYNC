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

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
truth_path <- cli$`truth-path` %||% ""
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify the completed Stage 6B v2 output.")
pp_assert(nzchar(truth_path) && file.exists(truth_path),
          "--truth-path must identify the sealed Stage 6A truth bundle.")

selection_path <- file.path(
  output_root, "truth_free_selection",
  "FINAL_TRUTH_FREE_DIMENSION_SELECTION.csv"
)
selection_marker <- file.path(
  output_root, "truth_free_selection",
  "FINAL_TRUTH_FREE_DIMENSION_SELECTION_FROZEN.txt"
)
authorization_path <- file.path(output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt")
pilot_complete <- file.path(output_root, "STAGE6B_PILOT_COMPLETE.txt")
evaluation_root <- file.path(output_root, "evaluation")

pp_assert(file.exists(selection_path) && file.exists(selection_marker) &&
            file.exists(pilot_complete),
          "The complete frozen truth-free Stage 6B selection is missing.")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6B evaluation already exists; refusing overwrite.")
selection_hash <- pp_sha256(selection_path)
authorization <- readLines(authorization_path, warn = FALSE)
required_authorization <- c(
  "status=TRUTH_UNSEAL_AND_EVALUATION_AUTHORIZED",
  paste0("selection_sha256=", selection_hash),
  "truth_unseal_authorized=TRUE",
  "evaluation_authorized=TRUE",
  "new_fits_authorized=FALSE",
  "continuation_authorized=FALSE",
  "final_full_data_refit_authorized=FALSE",
  "authorization_scope=stage6b_v2_existing_20_fits_only"
)
pp_assert(all(required_authorization %in% authorization),
          "Truth-unseal authorization is missing or does not bind selection.")

fit_root <- file.path(output_root, "fits")
fit_dirs <- list.dirs(fit_root, recursive = FALSE, full.names = TRUE)
fit_dirs <- fit_dirs[file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))]
pp_assert(length(fit_dirs) == 20L,
          "Stage 6B truth evaluation requires exactly 20 complete fits.")
pp_assert(length(list.files(fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE,
                            full.names = TRUE)) == 0L,
          "A Stage 6B fit error is present.")

records <- lapply(fit_dirs, function(path) {
  record <- readRDS(file.path(path, "terminal_record.rds"))
  pp_assert(isTRUE(record$objective_eligible) &&
              !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
              !isTRUE(record$truth_used_for_holdout_scoring) &&
              identical(record$terminal_status, "fixed_400_complete"),
            paste0("Ineligible Stage 6B endpoint: ", basename(path)))
  record
})
names(records) <- vapply(records, `[[`, character(1L), "fit_id")
records <- records[order(names(records))]
pp_assert(!anyDuplicated(names(records)), "Duplicate Stage 6B fit_id.")

factor_winners <- pp_read_csv(file.path(
  output_root, "truth_free_selection", "factor_screen_start_selection",
  "CONFIG_WINNERS.csv"
))
fpca_winners <- pp_read_csv(file.path(
  output_root, "truth_free_selection", "fpca_screen_start_selection",
  "CONFIG_WINNERS.csv"
))
winner_ids <- unique(c(factor_winners$fit_id, fpca_winners$fit_id))
final_selection <- pp_read_csv(selection_path)
pp_assert(nrow(final_selection) == 1L &&
            !isTRUE(final_selection$selection_used_truth[[1L]]),
          "Malformed final truth-free selection.")
final_config <- final_selection$selected_fit_config_id[[1L]]
winner_records <- records[names(records) %in% winner_ids]
final_endpoint_ids <- names(winner_records)[vapply(
  winner_records, function(record) identical(record$fit_config_id, final_config),
  logical(1L)
)]
pp_assert(length(final_endpoint_ids) == 1L,
          "Final selected configuration does not map to one frozen start winner.")

pp_load_package()
pp_source_runtime(TRUE)
truth <- readRDS(truth_path)
pp_assert(identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
            identical(truth$data_id, "g12s6a_01"),
          "Stage 6B truth identity check failed.")
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

dir.create(evaluation_root, recursive = TRUE, mode = "0700")
statuses <- vector("list", length(records))
compact <- vector("list", length(records))
metric_checks <- vector("list", length(records))

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

for (index in seq_along(records)) {
  record <- records[[index]]
  fit_id <- record$fit_id
  fit_dir <- file.path(fit_root, fit_id)
  output_dir <- file.path(evaluation_root, fit_id)
  dir.create(output_dir, recursive = TRUE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    value <- pp_evaluate_fit(
      readRDS(file.path(fit_dir, "fit.rds")), record, truth,
      evaluator_environment
    )
    oracle_dimensions <- record$fit_L_f == 1L &&
      record$fit_L_s_1 == 1L && record$fit_L_s_2 == 1L &&
      all(parse_rank(record$fit_M_f) == 2L) &&
      all(parse_rank(record$fit_M_s_1) == 2L) &&
      all(parse_rank(record$fit_M_s_2) == 2L)
    value$evaluation_contract$oracle_L_M <- oracle_dimensions
    value$evaluation_contract$stage6b_truth_free_selection_frozen <- TRUE
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
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(class = character(), message = character(), call = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(output_dir, "EVALUATION_WARNINGS.csv"))
  if (inherits(result, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
      paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(result))
    ), file.path(output_dir, "EVALUATION_ERROR.txt"))
    statuses[[index]] <- data.frame(
      fit_id = fit_id, ok = FALSE,
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      warning_count = nrow(warnings), error = conditionMessage(result),
      stringsAsFactors = FALSE
    )
    next
  }

  pp_save_rds(result, file.path(output_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(output_dir, paste0(name, ".csv")))
  }
  science <- cbind(
    result$summary, result$process_contribution_summary,
    result$loading_scale_summary
  )
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
    final_dimension_configuration = identical(record$fit_config_id, final_config),
    final_selected_endpoint = fit_id %in% final_endpoint_ids,
    stringsAsFactors = FALSE
  )
  compact[[index]] <- cbind(base, science)
  tables_present <- all(expected_tables %in% names(result)) &&
    all(vapply(result[expected_tables], function(value) {
      is.data.frame(value) && nrow(value) > 0L
    }, logical(1L)))
  checkpoint_complete <- identical(
    names(readRDS(file.path(fit_dir, "fit.rds"))$practical_checkpoints),
    c("380", "400")
  )
  metric_checks[[index]] <- data.frame(
    fit_id = fit_id, expected_table_count = length(expected_tables),
    expected_tables_present = tables_present,
    checkpoint_380_400_complete = checkpoint_complete,
    metric_contract_passed = tables_present && checkpoint_complete,
    stringsAsFactors = FALSE
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("selection_sha256=", selection_hash),
    "selection_frozen_before_truth_read=TRUE",
    "fit_or_selection_used_truth=FALSE", "new_fits_after_unseal=0",
    "continuation_started=FALSE", "final_full_data_refit_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_dir, "EVALUATION_COMPLETE.txt"))
  statuses[[index]] <- data.frame(
    fit_id = fit_id, ok = TRUE,
    elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
    warning_count = nrow(warnings), error = "", stringsAsFactors = FALSE
  )
}

status <- do.call(rbind, statuses)
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
if (!all(status$ok)) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("failed_fit_ids=", paste(status$fit_id[!status$ok], collapse = ";"))
  ), file.path(evaluation_root, "EVALUATION_ERROR.txt"))
  pp_abort("Stage 6B evaluation errors retained; no selective rerun.")
}

metric_qc <- do.call(rbind, metric_checks)
pp_assert(all(metric_qc$metric_contract_passed),
          "One or more Stage 6B evaluation metric contracts failed.")
all_science <- pp_rbind_fill(compact)
config_winners <- all_science[
  all_science$within_configuration_winner, , drop = FALSE
]
final_result <- all_science[
  all_science$final_selected_endpoint, , drop = FALSE
]
pp_assert(nrow(config_winners) == length(winner_ids) && nrow(final_result) == 1L,
          "Stage 6B evaluated winner mapping failed.")
pp_write_csv(metric_qc,
             file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
pp_write_csv(all_science,
             file.path(evaluation_root, "ALL_STAGE6B_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(config_winners,
             file.path(evaluation_root, "CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(final_result,
             file.path(evaluation_root, "FINAL_SELECTED_SCIENTIFIC_RESULT.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "fits_evaluated=20", "evaluation_errors=0",
  paste0("configuration_winners_evaluated=", nrow(config_winners)),
  paste0("final_selected_fit_id=", final_endpoint_ids),
  paste0("selection_sha256=", selection_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "final_full_data_refit_started=FALSE",
  "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE6B_EVALUATION_COMPLETE fits=20 errors=0 final_fit=",
    final_endpoint_ids, "\n", sep = "")
