#!/usr/bin/env Rscript

options(warn = 1)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4e_common_20260827_v1.R"),
       local = FALSE)

output_root <- file.path(UC_ROOT, "truth_evaluation_20260827_v1")
uc_assert(!dir.exists(output_root),
          "Refusing to overwrite Stage-4E truth evaluation output.")
uc_assert(file.exists(file.path(UC_ROOT,
                                "TRUTH_FREE_SELECTION_COMPLETE.txt")),
          "Truth-free Stage-4E selection has not been frozen.")
uc_assert(file.exists(file.path(
  UC_ROOT, "STAGE4E_TRUTH_EVALUATION_SCOPE_20260827.txt"
)), "Stage-4E truth-evaluation scope record is missing.")

inputs <- s4e_verify_inputs(load_runtime = TRUE)
inventory <- uc_read_csv(file.path(
  UC_ROOT, "ALL_24_TRUTH_FREE_ENDPOINTS.csv"
))
winners <- uc_read_csv(file.path(
  UC_ROOT, "TRUTH_FREE_ELBO_WINNERS_8_STARTS.csv"
))
uc_assert(
  nrow(inventory) == 24L && all(table(inventory$data_id) == 8L) &&
    all(inventory$objective_eligible) && !any(inventory$truth_used) &&
    nrow(winners) == 3L && !any(winners$selection_used_truth) &&
    !any(winners$selection_used_structure) &&
    !any(winners$selection_used_NRMSE_ISE_RPL),
  "Frozen truth-free Stage-4E selection is invalid."
)

fit_paths <- ifelse(
  inventory$source_stage == "stage4c_inherited",
  file.path(S4E_PARENT_STAGE4C_ROOT, "fits", inventory$fit_id, "fit.rds"),
  file.path(UC_ROOT, "fits", inventory$fit_id, "fit.rds")
)
truth_paths <- file.path(
  S4E_PARENT_STAGE4C_ROOT, "data", inventory$data_id,
  "sealed_truth", "truth_bundle.rds"
)
uc_assert(all(file.exists(fit_paths)) && all(file.exists(truth_paths)),
          "A fit or previously unsealed Stage-4C truth bundle is missing.")

wrapper <- uc_find_one(UC_SNAPSHOT, "^v0lv_candidate2_evaluation[.]R$")
sys.source(wrapper, envir = globalenv(), keep.source = TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(UC_SNAPSHOT)

dir.create(output_root, recursive = TRUE, mode = "0700")
action_rows <- vector("list", nrow(inventory))
for (index in seq_len(nrow(inventory))) {
  row <- inventory[index, , drop = FALSE]
  started <- Sys.time()
  warning_messages <- character()
  error_message <- ""
  evaluation <- NULL
  tryCatch({
    evaluation <- withCallingHandlers({
      fit <- readRDS(fit_paths[[index]])
      truth <- readRDS(truth_paths[[index]])
      uc_assert(identical(truth$bundle_class,
                          "v0lv_sealed_truth_candidate2") &&
                  identical(truth$data_id, row$data_id[[1L]]),
                "Truth bundle identity mismatch.")
      candidate2_evaluate_multi(
        fit, truth, UC_SNAPSHOT,
        evaluator_environment = evaluator_environment,
        mode = "smoke", binding_dir = NULL
      )
    }, warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  }, error = function(condition) {
    error_message <<- conditionMessage(condition)
  })
  ended <- Sys.time()
  status <- if (nzchar(error_message)) "ERROR" else "COMPLETE"
  if (identical(status, "COMPLETE")) {
    endpoint_dir <- file.path(output_root, "endpoints", row$fit_id[[1L]])
    dir.create(endpoint_dir, recursive = TRUE, mode = "0700")
    saveRDS(evaluation, file.path(endpoint_dir, "evaluation.rds"), version = 3)
    for (name in names(evaluation)) {
      if (is.data.frame(evaluation[[name]])) {
        utils::write.csv(
          evaluation[[name]], file.path(endpoint_dir, paste0(name, ".csv")),
          row.names = FALSE, quote = TRUE, na = ""
        )
      }
    }
  }
  action_rows[[index]] <- data.frame(
    fit_id = row$fit_id[[1L]], data_id = row$data_id[[1L]],
    scenario_id = row$scenario_id[[1L]],
    seed_index = row$seed_index[[1L]], source_stage = row$source_stage[[1L]],
    endpoint_role = if (row$fit_id[[1L]] %in% winners$selected_fit_id)
      "winner" else "not_selected",
    frozen_final_elbo = row$final_elbo[[1L]], status = status,
    warning_count = length(warning_messages),
    warnings = paste(unique(warning_messages), collapse = " | "),
    error = error_message, started_utc = uc_iso_time(started),
    ended_utc = uc_iso_time(ended),
    elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
    selection_frozen_before_truth_read = TRUE,
    truth_used_for_fit = FALSE, truth_used_for_stopping = FALSE,
    truth_used_for_selection = FALSE, truth_used_for_evaluation = TRUE,
    continuation_run = FALSE, automatic_800_run = FALSE,
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )
  utils::write.csv(
    do.call(rbind, action_rows[seq_len(index)]),
    file.path(output_root, "EVALUATION_ACTIONS.csv"),
    row.names = FALSE, quote = TRUE, na = ""
  )
}

actions <- do.call(rbind, action_rows)
if (any(actions$status != "COMPLETE")) {
  writeLines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("completed_utc=", uc_iso_time()),
    paste0("failed_endpoints=", paste(
      actions$fit_id[actions$status != "COMPLETE"], collapse = ";"
    ))
  ), file.path(output_root, "EVALUATION_ERROR.txt"), useBytes = TRUE)
  uc_abort("Stage-4E evaluation retained one or more errors.")
}

writeLines(c(
  "status=PASS", paste0("completed_utc=", uc_iso_time()),
  "evaluated_endpoints=24", "target_truth_bundles=3",
  "truth_free_winners_preserved=3",
  "selection_frozen_before_truth_read=TRUE",
  "truth_used_for_selection=FALSE", "truth_used_for_evaluation=TRUE",
  "continuation_run=FALSE", "automatic_800_run=FALSE",
  "formal_v0lv_result=FALSE"
), file.path(output_root, "EVALUATION_COMPLETE.txt"), useBytes = TRUE)
utils::capture.output(
  sessionInfo(), file = file.path(output_root, "SESSION_INFO.txt")
)
cat("STAGE4E_TRUTH_EVALUATION_PASS endpoints=24 targets=3\n")
