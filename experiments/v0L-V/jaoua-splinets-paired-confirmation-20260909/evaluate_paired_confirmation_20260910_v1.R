#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(
  script_root, "..", "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)

PP_ROOT <- script_root
PP_EXPERIMENT_ID <- "JAOUA_SPLINETS_PAIRED_CONFIRMATION_V1_20260909"
PP_EXPECTED_FITS <- 72L
PP_EXPECTED_WINNERS <- 6L

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% pp_env("JSPC_EVAL_WORKERS", "7"))
library_path <- cli$library %||% pp_env("JSPC_R_LIB")
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify the completed paired experiment.")
pp_assert(.Platform$OS.type == "unix" && length(workers) == 1L &&
            is.finite(workers) && workers %in% 1:8,
          "Paired evaluation requires Unix and 1--8 workers.")
if (nzchar(library_path)) Sys.setenv(STAGE6A_R_LIB = library_path)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

selection_root <- file.path(output_root, "truth_free_selection")
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "G12_TRUTH_FREE_WINNERS.csv")
selection_marker <- file.path(selection_root,
                              "TRUTH_FREE_SELECTION_FROZEN.txt")
experiment_complete <- file.path(output_root,
                                 "FITS_AND_SELECTION_COMPLETE.txt")
authorization_path <- file.path(output_root,
                                "TRUTH_UNSEAL_AUTHORIZATION.txt")
evaluation_root <- file.path(output_root, "truth_evaluation_20260910_v1")

pp_assert(all(file.exists(c(endpoints_path, winners_path, selection_marker,
                            experiment_complete))),
          "Complete frozen truth-free inputs are missing.")
pp_assert(!dir.exists(evaluation_root),
          "Paired evaluation directory already exists; refusing overwrite.")
endpoints_hash <- pp_sha256(endpoints_path)
selection_hash <- pp_sha256(winners_path)
pp_assert(file.exists(authorization_path),
          "Explicit paired truth-unseal authorization is missing.")
authorization <- readLines(authorization_path, warn = FALSE)
required_authorization <- c(
  "status=TRUTH_UNSEAL_AND_EVALUATION_AUTHORIZED",
  paste0("all_endpoints_sha256=", endpoints_hash),
  paste0("winner_selection_sha256=", selection_hash),
  "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
  "new_fits_authorized=FALSE", "continuation_authorized=FALSE",
  "authorization_scope=existing_72_paired_fits_only"
)
pp_assert(all(required_authorization %in% authorization),
          "Authorization does not bind the frozen paired endpoints/winners.")

fit_root <- file.path(output_root, "fits")
fit_dirs <- list.dirs(fit_root, recursive = FALSE, full.names = TRUE)
fit_dirs <- fit_dirs[file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))]
pp_assert(length(fit_dirs) == PP_EXPECTED_FITS,
          "Paired truth evaluation requires exactly 72 complete fits.")
pp_assert(length(list.files(fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE,
                            full.names = TRUE)) == 0L,
          "A paired fit error is present.")

records <- lapply(fit_dirs, function(path) {
  record <- readRDS(file.path(path, "terminal_record.rds"))
  pp_assert(
    identical(record$terminal_schema,
              "JAOUA_SPLINETS_PAIRED_CONFIRMATION_TERMINAL_V1") &&
      isTRUE(record$objective_eligible) &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$continuation_used) &&
      identical(record$terminal_status, "fixed_400_complete"),
    paste0("Ineligible paired endpoint: ", basename(path))
  )
  record
})
names(records) <- vapply(records, `[[`, character(1L), "fit_id")
records <- records[order(names(records))]
pp_assert(!anyDuplicated(names(records)), "Duplicate paired fit_id.")

endpoints <- pp_read_csv(endpoints_path)
winners <- pp_read_csv(winners_path)
pp_assert(nrow(endpoints) == PP_EXPECTED_FITS &&
            nrow(winners) == PP_EXPECTED_WINNERS &&
            all(table(endpoints$selection_stratum_id) == 12L) &&
            all(table(winners$selection_stratum_id) == 1L) &&
            !any(endpoints$truth_used_for_fit_stopping_or_selection),
          "Malformed frozen paired endpoint/winner tables.")
winner_ids <- winners$fit_id
pp_assert(all(winner_ids %in% names(records)),
          "A frozen winner does not map to a terminal record.")

pp_load_package()
pp_source_runtime(TRUE)
data_ids <- unique(endpoints$data_id)
truth_paths <- file.path(output_root, "data", data_ids,
                         "sealed_truth", "truth_bundle.rds")
pp_assert(length(data_ids) == 6L && all(file.exists(truth_paths)),
          "A paired truth bundle is missing.")
truths <- lapply(truth_paths, readRDS)
names(truths) <- data_ids
pp_assert(all(vapply(seq_along(truths), function(index) {
  truth <- truths[[index]]
  identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
    identical(truth$data_id, names(truths)[[index]]) &&
    identical(truth$experiment_id, PP_EXPERIMENT_ID)
}, logical(1L))), "Paired truth identity check failed.")
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

dir.create(evaluation_root, recursive = TRUE, mode = "0700")
pp_write_lines(c(
  "status=EVALUATION_STARTED", paste0("started_utc=", pp_iso_time()),
  paste0("workers=", workers), "registered_fits=72",
  paste0("all_endpoints_sha256=", endpoints_hash),
  paste0("winner_selection_sha256=", selection_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "development_only=TRUE",
  "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_STARTED.txt"))

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
  output_dir <- file.path(evaluation_root, "endpoints", fit_id)
  dir.create(output_dir, recursive = TRUE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    value <- pp_evaluate_fit(
      readRDS(file.path(fit_dir, "fit.rds")), record,
      truths[[record$data_id]], evaluator_environment
    )
    value$evaluation_contract$oracle_L_M <- TRUE
    value$evaluation_contract$paired_generator_sensitivity <- TRUE
    value$evaluation_contract$selection_frozen_before_truth_read <- TRUE
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
    logical_data_id = record$logical_data_id,
    basis_arm = record$basis_arm,
    truth_basis_method = record$truth_basis_method,
    paired_start_id = record$paired_start_id,
    scenario_id = record$scenario_id, method_id = record$method_id,
    fit_config_id = record$fit_config_id,
    selection_stratum_id = record$selection_stratum_id,
    seed_index = record$seed_index, fit_seed = record$fit_seed,
    final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    convergence_status = record$convergence_status,
    elapsed_seconds = record$elapsed_seconds,
    warning_count = if (is.data.frame(record$warnings)) {
      nrow(record$warnings)
    } else {
      length(record$warnings %||% character())
    },
    truth_free_winner = fit_id %in% winner_ids,
    stringsAsFactors = FALSE
  )
  compact <- cbind(base, science)
  tables_present <- all(expected_tables %in% names(result)) &&
    all(vapply(result[expected_tables], function(value) {
      is.data.frame(value) && nrow(value) > 0L
    }, logical(1L)))
  metric <- data.frame(
    fit_id = fit_id, expected_table_count = length(expected_tables),
    expected_tables_present = tables_present,
    metric_contract_passed = tables_present,
    stringsAsFactors = FALSE
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("all_endpoints_sha256=", endpoints_hash),
    paste0("winner_selection_sha256=", selection_hash),
    "selection_frozen_before_truth_read=TRUE",
    "fit_or_selection_used_truth=FALSE", "new_fits_after_unseal=0",
    "continuation_started=FALSE", "development_only=TRUE",
    "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"
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
  pp_abort("Paired evaluation errors retained; no selective rerun.")
}

metric_qc <- do.call(rbind, lapply(results, `[[`, "metric"))
pp_assert(nrow(metric_qc) == PP_EXPECTED_FITS &&
            all(metric_qc$metric_contract_passed),
          "One or more paired evaluation metric contracts failed.")
all_science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
winner_science <- all_science[all_science$truth_free_winner, , drop = FALSE]
pp_assert(nrow(all_science) == PP_EXPECTED_FITS &&
            nrow(winner_science) == PP_EXPECTED_WINNERS,
          "Paired evaluated winner mapping failed.")

pp_write_csv(metric_qc,
             file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
pp_write_csv(all_science,
             file.path(evaluation_root, "ALL_PAIRED_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(winner_science,
             file.path(evaluation_root, "WINNER_SCIENTIFIC_RESULTS.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "fits_evaluated=72", "evaluation_errors=0", "winners_evaluated=6",
  paste0("all_endpoints_sha256=", endpoints_hash),
  paste0("winner_selection_sha256=", selection_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "development_only=TRUE",
  "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("PAIRED_EVALUATION_COMPLETE fits=72 errors=0 winners=6\n")
