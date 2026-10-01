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
source(file.path(experiments_root,
                 "g12-dimension-overspec-stage6a-20260831",
                 "stage6a_common_20260831_v1.R"), local = FALSE)
source(file.path(script_root,
                 "stage6f_joint_selection_tools_20260911_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6F_JOINT_DIMENSION_SELECTION_V1_20260911"
stage6f_library <- pp_env("STAGE6F_R_LIB")
if (nzchar(stage6f_library)) Sys.setenv(STAGE6A_R_LIB = stage6f_library)
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6F_OUTPUT_ROOT")
truth_path <- cli$`truth-path` %||% pp_env("STAGE6F_TRUTH_BUNDLE")
authorization_token <- cli$`authorization-token` %||% ""

pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify the completed Stage 6F output.")
pp_assert(nzchar(truth_path) && file.exists(truth_path),
          "--truth-path must identify the sealed Stage 6A truth bundle.")
pp_assert(identical(authorization_token, "USER_AUTHORIZED_STAGE6F_20260911"),
          "Explicit Stage 6F truth-evaluation authorization token is missing.")

selection_root <- file.path(output_root, "truth_free_selection")
winner_path <- file.path(selection_root, "CONFIG_SPLIT_WINNERS.csv")
dimension_path <- file.path(
  selection_root, "FINAL_TRUTH_FREE_DIMENSION_SELECTION.csv"
)
dimension_marker <- file.path(selection_root,
                              "DIMENSION_SELECTION_FROZEN.txt")
truth_free_marker <- file.path(output_root, "STAGE6F_TRUTH_FREE_COMPLETE.txt")
authorization_path <- file.path(output_root,
                                "TRUTH_EVALUATION_AUTHORIZED.txt")
evaluation_root <- file.path(output_root, "truth_evaluation")
pp_assert(all(file.exists(c(winner_path, dimension_path, dimension_marker,
                            truth_free_marker))),
          "Frozen Stage 6F truth-free selection is incomplete.")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6F truth evaluation already exists; refusing overwrite.")
pp_assert(length(list.files(file.path(output_root, "fits"),
                            "^FIT_COMPLETE[.]txt$", recursive = TRUE,
                            full.names = TRUE)) == 32L,
          "Stage 6F truth evaluation requires 32 complete endpoints.")
pp_assert(length(list.files(output_root, "^(FIT|SCORE)_ERROR[.]txt$",
                            recursive = TRUE, full.names = TRUE)) == 0L,
          "Stage 6F contains a retained fit or score error.")

winner_hash <- pp_sha256(winner_path)
dimension_hash <- pp_sha256(dimension_path)
if (!file.exists(authorization_path)) {
  pp_write_lines(c(
    "status=TRUTH_EVALUATION_AUTHORIZED",
    "authorization_source=user_request_20260911_add_20_fits_and_assess_criterion",
    "authorization_token=USER_AUTHORIZED_STAGE6F_20260911",
    paste0("start_winner_table_sha256=", winner_hash),
    paste0("dimension_selection_sha256=", dimension_hash),
    "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
    "new_fits_authorized_after_unseal=FALSE",
    "continuation_authorized=FALSE", "full_data_refit_authorized=FALSE",
    "scope=eight_stage6f_split_configuration_winners_only"
  ), authorization_path)
}
authorization <- readLines(authorization_path, warn = FALSE)
pp_assert(all(c(
  "status=TRUTH_EVALUATION_AUTHORIZED",
  paste0("start_winner_table_sha256=", winner_hash),
  paste0("dimension_selection_sha256=", dimension_hash),
  "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
  "new_fits_authorized_after_unseal=FALSE",
  "continuation_authorized=FALSE"
) %in% authorization),
"Stage 6F authorization marker is incomplete or selection-mismatched.")

winners <- pp_read_csv(winner_path)
selections <- pp_read_csv(dimension_path)
pp_assert(nrow(winners) == 8L && nrow(selections) == 2L &&
            !any(selections$truth_used),
          "Stage 6F truth-free winner/selection tables are malformed.")
records <- lapply(winners$fit_id, function(fit_id) {
  record <- readRDS(file.path(output_root, "fits", fit_id,
                              "terminal_record.rds"))
  pp_assert(isTRUE(record$objective_eligible) &&
              !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
              identical(record$terminal_status, "fixed_400_complete"),
            paste0("Ineligible Stage 6F winner: ", fit_id))
  record
})
names(records) <- winners$fit_id

pp_load_package()
pp_source_runtime(TRUE)
truth <- readRDS(truth_path)
pp_assert(identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
            identical(truth$data_id, "g12s6a_01"),
          "Stage 6F truth identity check failed.")
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")
statuses <- vector("list", nrow(winners))
compact <- vector("list", nrow(winners))
for (index in seq_len(nrow(winners))) {
  winner <- winners[index, , drop = FALSE]
  fit_id <- winner$fit_id[[1L]]
  record <- records[[fit_id]]
  fit_dir <- file.path(output_root, "fits", fit_id)
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    value <- pp_evaluate_fit(
      readRDS(file.path(fit_dir, "fit.rds")), record, truth,
      evaluator_environment
    )
    value$evaluation_contract$oracle_L_M <-
      identical(winner$fit_config_id[[1L]], "joint_L1_M2")
    value$evaluation_contract$stage6f_truth_free_selection_frozen <- TRUE
    value$evaluation_contract$heldout_observations_excluded_from_fit <- TRUE
    value$evaluation_contract$selection_hash <- dimension_hash
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
  pp_write_csv(warnings, file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  if (inherits(result, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(result))
    ), file.path(result_dir, "EVALUATION_ERROR.txt"))
    statuses[[index]] <- data.frame(
      fit_id = fit_id, ok = FALSE,
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      warning_count = nrow(warnings), error = conditionMessage(result),
      stringsAsFactors = FALSE
    )
    next
  }
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]],
                 file.path(result_dir, paste0(name, ".csv")))
  }
  base <- data.frame(
    fit_id = fit_id, split_id = winner$split_id[[1L]],
    fit_config_id = winner$fit_config_id[[1L]],
    seed_index = winner$seed_index[[1L]], fit_seed = winner$fit_seed[[1L]],
    final_elbo = winner$final_elbo[[1L]],
    practical_converged = winner$practical_converged[[1L]],
    imported_from_stage6b = winner$imported_from_stage6b[[1L]],
    generating_dimensions = identical(winner$fit_config_id[[1L]],
                                      "joint_L1_M2"),
    selected_by_mse_one_se = winner$fit_config_id[[1L]] %in%
      selections$selected_config_id[
        selections$criterion_id == "paired_holdout_mse"
      ],
    selected_by_joint_nlpd_one_se = winner$fit_config_id[[1L]] %in%
      selections$selected_config_id[
        selections$criterion_id == "experimental_joint_nlpd"
      ],
    stringsAsFactors = FALSE
  )
  compact[[index]] <- cbind(
    base, result$summary, result$process_contribution_summary,
    result$loading_scale_summary
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("dimension_selection_sha256=", dimension_hash),
    "selection_frozen_before_truth_read=TRUE",
    "fit_start_and_dimension_selection_used_truth=FALSE",
    "new_fits_after_unseal=0", "continuation_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
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
  pp_abort("Stage 6F evaluation errors retained; no selective rerun.")
}

all_science <- pp_rbind_fill(compact)
pp_assert(nrow(all_science) == 8L,
          "Stage 6F compact scientific result must contain eight winners.")
pp_write_csv(all_science,
             file.path(evaluation_root,
                       "SPLIT_CONFIG_WINNER_SCIENTIFIC_RESULTS.csv"))

numeric_columns <- names(all_science)[vapply(all_science, is.numeric,
                                             logical(1L))]
excluded_numeric <- c("seed_index", "fit_seed", "final_elbo")
metric_columns <- setdiff(numeric_columns, excluded_numeric)
config_science <- do.call(rbind, lapply(
  sort(unique(all_science$fit_config_id)), function(config_id) {
    rows <- all_science[all_science$fit_config_id == config_id,
                        , drop = FALSE]
    base <- data.frame(
      fit_config_id = config_id, split_winners = nrow(rows),
      generating_dimensions = identical(config_id, "joint_L1_M2"),
      stringsAsFactors = FALSE
    )
    for (column in metric_columns) {
      value <- as.numeric(rows[[column]])
      base[[paste0(column, "_mean")]] <-
        if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
    }
    base
  }
))
pp_write_csv(config_science,
             file.path(evaluation_root,
                       "CONFIG_SCIENTIFIC_RESULTS_ACROSS_SPLITS.csv"))

assessment <- merge(
  selections,
  config_science,
  by.x = "selected_config_id", by.y = "fit_config_id",
  all.x = TRUE, sort = FALSE
)
assessment$selected_generating_dimensions <-
  assessment$selected_config_id == "joint_L1_M2"
assessment$pilot_scope <- "one_data_seed_two_paired_splits_development_only"
pp_write_csv(assessment,
             file.path(evaluation_root, "CRITERION_SCIENTIFIC_ASSESSMENT.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "split_configuration_winners_evaluated=8", "evaluation_errors=0",
  paste0("mse_selected=", selections$selected_config_id[
    selections$criterion_id == "paired_holdout_mse"
  ]),
  paste0("joint_nlpd_selected=", selections$selected_config_id[
    selections$criterion_id == "experimental_joint_nlpd"
  ]),
  paste0("dimension_selection_sha256=", dimension_hash),
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE6F_EVALUATION_COMPLETE winners=8 errors=0\n")
