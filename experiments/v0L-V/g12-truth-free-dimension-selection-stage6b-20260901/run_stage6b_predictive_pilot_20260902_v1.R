#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}

stage6a_root <- file.path(dirname(script_root),
                          "g12-dimension-overspec-stage6a-20260831")
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)
source(file.path(script_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = FALSE)
source(file.path(script_root,
                 "stage6b_predictive_pilot_tools_20260902_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6B_PREDICTIVE_PILOT_V1_20260902"
PP_WORKERS <- as.integer(pp_env("STAGE6B_WORKERS", "3"))
pp_assert(length(PP_WORKERS) == 1L && is.finite(PP_WORKERS) &&
            PP_WORKERS %in% 1:3,
          "STAGE6B_WORKERS must be an integer from 1 to 3.")
stage6b_library <- pp_env("STAGE6B_R_LIB")
if (nzchar(stage6b_library)) Sys.setenv(STAGE6A_R_LIB = stage6b_library)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6B_OUTPUT_ROOT")
observation_path <- cli$`observation-bundle` %||%
  pp_env("STAGE6B_OBSERVATION_BUNDLE")

config_path <- file.path(script_root, "FIT_CONFIGS_STAGE6B_PILOT_V1.csv")
seed_path <- file.path(script_root, "START_SEEDS_STAGE6B_PILOT_V1.csv")

require_output_root <- function(create = TRUE) {
  pp_assert(nzchar(output_root), "--output-root is required.")
  if (create && !dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  output_root
}

load_registration <- function() {
  configs <- s6bp_read_csv(config_path)
  seeds <- s6bp_read_csv(seed_path)
  s6bp_validate_registration(configs, seeds)
  list(
    configs = configs,
    seeds = seeds,
    fits = s6bp_build_fit_manifest(configs, seeds)
  )
}

verify_environment <- function() {
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  registration <- load_registration()
  pp_load_package()
  pp_source_runtime(FALSE)
  pp_assert(nrow(registration$fits) == 36L &&
              nrow(s6bp_phase_manifest(registration$fits,
                                       "factor_screen")) == 12L,
            "Stage 6B registration integrity failed.")
  registration
}

prepared_paths <- function() list(
  train = file.path(output_root, "split", "train_observation.rds"),
  holdout = file.path(output_root, "split", "holdout_plan.rds"),
  manifest = file.path(output_root, "FIT_MANIFEST_CONDITIONAL.csv"),
  marker = file.path(output_root, "STAGE6B_PREPARED.txt")
)

prepare_pilot <- function() {
  require_output_root(TRUE)
  registration <- verify_environment()
  pp_assert(nzchar(observation_path) && file.exists(observation_path),
            "The Stage 6A observation-only bundle is missing.")
  paths <- prepared_paths()
  source_sha <- pp_sha256(observation_path)
  if (file.exists(paths$marker)) {
    pp_assert(all(file.exists(unlist(paths[c("train", "holdout", "manifest")]))) &&
                any(readLines(paths$marker, warn = FALSE) ==
                      paste0("source_observation_sha256=", source_sha)),
              "Existing Stage 6B preparation is incomplete or source-mismatched.")
    cat("STAGE6B_ALREADY_PREPARED output_root=", output_root, "\n", sep = "")
    return(invisible(paths))
  }
  pp_assert(!dir.exists(file.path(output_root, "fits")) &&
              !dir.exists(file.path(output_root, "truth_free_selection")),
            "Refusing to prepare over an existing Stage 6B fit/selection tree.")

  observation <- readRDS(observation_path)
  pp_assert(identical(observation$bundle_class,
                      "v0lv_observation_only_candidate2") &&
              identical(observation$data_id, "g12s6a_01") &&
              !length(s6b_forbidden_observation_names(observation)),
            "Source observation identity or isolation failed.")
  split <- s6b_make_middle_time_holdout(observation, min_train_times = 5L)
  pp_assert(split$holdout_plan$subject_time_count == 60L &&
              all(vapply(split$holdout_plan$records, `[[`, integer(1L),
                         "p") == 500L) &&
              identical(split$holdout_plan$truth_used, FALSE),
            "Registered full-size Stage 6B split contract failed.")
  pp_save_rds(split$train_observation, paths$train)
  pp_save_rds(split$holdout_plan, paths$holdout)
  pp_write_csv(registration$fits, paths$manifest)

  records <- split$holdout_plan$records
  split_summary <- data.frame(
    split_id = split$holdout_plan$split_id,
    study = vapply(records, `[[`, integer(1L), "study"),
    subject = vapply(records, `[[`, integer(1L), "subject"),
    original_index = vapply(records, `[[`, integer(1L), "original_index"),
    holdout_time = vapply(records, `[[`, numeric(1L), "time"),
    train_time_count = vapply(records, `[[`, integer(1L), "train_time_count"),
    variables = vapply(records, `[[`, integer(1L), "p"),
    truth_used = FALSE,
    stringsAsFactors = FALSE
  )
  pp_write_csv(split_summary,
               file.path(output_root, "split", "HOLDOUT_SPLIT_SUMMARY.csv"))
  pp_write_lines(c(
    "status=STAGE6B_PREPARED",
    paste0("prepared_utc=", pp_iso_time()),
    paste0("source_observation_sha256=", source_sha),
    "source_data_id=g12s6a_01", "new_data_generated=FALSE",
    "heldout_subject_time_vectors=60", "variables_per_vector=500",
    "fit_manifest_rows_registered_conditionally=36",
    "maximum_fit_rows_to_execute=20", "truth_read=FALSE",
    "truth_access_authorized=FALSE", "continuation_authorized=FALSE",
    "formal_paper_mc_result=FALSE"
  ), paths$marker)
  cat("STAGE6B_PREPARED heldout=60 conditional_fits=36 max_executed=20\n")
  invisible(paths)
}

load_prepared <- function() {
  require_output_root(FALSE)
  paths <- prepared_paths()
  pp_assert(all(file.exists(unlist(paths))),
            "Stage 6B preparation files are incomplete.")
  train <- readRDS(paths$train)
  holdout <- readRDS(paths$holdout)
  pp_assert(!length(s6b_forbidden_observation_names(train)) &&
              identical(holdout$truth_used, FALSE) &&
              holdout$subject_time_count == 60L,
            "Prepared Stage 6B split is invalid or truth-bearing.")
  list(train = train, holdout = holdout)
}

factor_selection_path <- function() file.path(
  output_root, "truth_free_selection", "FACTOR_DIMENSION_SELECTION.csv"
)

selected_factor_count <- function() {
  path <- factor_selection_path()
  pp_assert(file.exists(path), "The factor-count selection is not frozen.")
  selection <- pp_read_csv(path)
  pp_assert(nrow(selection) == 1L &&
              as.integer(selection$selected_factor_count[[1L]]) %in% 1:3 &&
              identical(selection$selection_used_truth[[1L]], FALSE),
            "The frozen factor-count selection is malformed.")
  as.integer(selection$selected_factor_count[[1L]])
}

fit_allowed <- function(row) {
  if (identical(row$phase[[1L]], "factor_screen")) return(TRUE)
  selected <- selected_factor_count()
  identical(as.integer(row$fit_L_f[[1L]]), selected)
}

write_fit_error <- function(fit_dir, fit_id, message) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("fit_id=", fit_id), paste0("message=", message),
    paste0("ended_utc=", pp_iso_time()), "truth_read=FALSE"
  ), file.path(fit_dir, "FIT_ERROR.txt"))
}

run_fit <- function(fit_id) {
  require_output_root(FALSE)
  registration <- verify_environment()
  split <- load_prepared()
  pp_assert(!file.exists(file.path(output_root,
                                   "TRUTH_UNSEAL_AUTHORIZATION.txt")),
            "Refusing new Stage 6B fits after truth-access authorization.")
  row <- registration$fits[registration$fits$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L && fit_allowed(row),
            "fit_id is absent or not authorized by the sequential selection.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")) &&
                file.exists(file.path(fit_dir, "fit.rds")) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Pre-existing fit is incomplete/error and cannot be rerun: ",
                     fit_id))
    cat("FIT_ALREADY_COMPLETE fit_id=", fit_id, "\n", sep = "")
    return(invisible(0L))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- pp_fit_g12(split$train, row, smoke = FALSE)
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  result$captured$warnings <- s6bp_set_warning_phase(
    result$captured$warnings,
    "stage6b_predictive_pilot_g12_fixed400"
  )
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6B_PREDICTIVE_PILOT_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$holdout_split_id <- split$holdout$split_id
  record$heldout_subject_time_vectors <- split$holdout$subject_time_count
  record$primary_dimension_score <- "heldout_nrmse"
  record$holdout_scoring_complete <- FALSE
  record$truth_used_for_fit_stopping_or_selection <- FALSE
  record$truth_used_for_holdout_scoring <- FALSE
  record$final_full_data_refit <- FALSE

  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(), call = character(),
               phase = character(), stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
    write_fit_error(fit_dir, fit_id, record$error$message %||% "unknown")
    pp_abort("Stage 6B fit failed and was retained: ", fit_id)
  }

  pp_write_csv(data.frame(
    t1_sweep = seq_along(result$fit$ELBO),
    ordinary_elbo = as.numeric(result$fit$ELBO)
  ), file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  if (is.data.frame(result$fit$practical_diagnostics)) {
    pp_write_csv(result$fit$practical_diagnostics,
                 file.path(fit_dir, "PRACTICAL_DIAGNOSTICS.csv"))
  }
  if (is.data.frame(result$fit$random_scale_calibration_diagnostics)) {
    pp_write_csv(result$fit$random_scale_calibration_diagnostics,
                 file.path(fit_dir, "INITIALIZATION_DIAGNOSTICS.csv"))
  }

  scores <- tryCatch(list(
    all = s6b_score_holdout(result$fit, split$holdout, "all"),
    ppi = s6b_score_holdout(result$fit, split$holdout, "ppi_0.5")
  ), error = function(condition) condition)
  if (inherits(scores, "condition")) {
    record$holdout_scoring_error <- conditionMessage(scores)
    pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
    write_fit_error(fit_dir, fit_id,
                    paste0("holdout scoring: ", conditionMessage(scores)))
    pp_abort("Stage 6B holdout scoring failed and was retained: ", fit_id)
  }
  pp_write_csv(scores$all$subject_time,
               file.path(fit_dir, "HOLDOUT_SUBJECT_TIME_ALL.csv"))
  pp_write_csv(scores$all$aggregate,
               file.path(fit_dir, "HOLDOUT_SUMMARY_ALL.csv"))
  pp_write_csv(scores$ppi$subject_time,
               file.path(fit_dir, "HOLDOUT_SUBJECT_TIME_PPI_0_5.csv"))
  pp_write_csv(scores$ppi$aggregate,
               file.path(fit_dir, "HOLDOUT_SUMMARY_PPI_0_5.csv"))
  record$holdout_scoring_complete <- TRUE
  record$heldout_rmse_all <- scores$all$aggregate$heldout_rmse[[1L]]
  record$heldout_nrmse_all <- scores$all$aggregate$heldout_nrmse[[1L]]
  record$heldout_rmse_ppi_0_5 <- scores$ppi$aggregate$heldout_rmse[[1L]]
  record$heldout_nrmse_ppi_0_5 <- scores$ppi$aggregate$heldout_nrmse[[1L]]
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "terminal_status=fixed_400_complete",
    paste0("objective_eligible=", record$objective_eligible),
    "holdout_scoring_complete=TRUE", "primary_score=heldout_nrmse",
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "truth_used_for_holdout_scoring=FALSE", "n_cpus=1",
    "continuation_used=FALSE", "final_full_data_refit=FALSE",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), " heldout_nrmse=",
      format(record$heldout_nrmse_all, digits = 8), "\n", sep = "")
  invisible(0L)
}

collect_terminals <- function(manifest) {
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  complete <- file.path(output_root, "fits", manifest$fit_id,
                        "FIT_COMPLETE.txt")
  score_paths <- file.path(output_root, "fits", manifest$fit_id,
                           "HOLDOUT_SUBJECT_TIME_ALL.csv")
  pp_assert(all(file.exists(paths)) && all(file.exists(complete)) &&
              all(file.exists(score_paths)),
            "Stage 6B phase has missing terminal, completion or score files.")
  records <- lapply(paths, readRDS)
  pp_assert(all(vapply(records, function(record) {
    identical(record$terminal_schema,
              "G12_STAGE6B_PREDICTIVE_PILOT_TERMINAL_V1") &&
      isTRUE(record$holdout_scoring_complete) &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$truth_used_for_holdout_scoring)
  }, logical(1L))), "Stage 6B terminal contract failed.")
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  pp_assert(identical(rows$fit_id, manifest$fit_id) &&
              all(rows$terminal_status == "fixed_400_complete") &&
              all(rows$objective_eligible) &&
              !any(rows$truth_used_for_fit_stopping_or_selection),
            "Stage 6B phase contains an ineligible or truth-bearing endpoint.")
  list(records = records, rows = rows)
}

selection_phase_dir <- function(phase) file.path(
  output_root, "truth_free_selection", paste0(phase, "_start_selection")
)

freeze_start_winners <- function(manifest, phase) {
  phase_dir <- selection_phase_dir(phase)
  marker <- file.path(phase_dir, "START_SELECTION_FROZEN.txt")
  if (file.exists(marker)) {
    winners <- pp_read_csv(file.path(phase_dir, "CONFIG_WINNERS.csv"))
    pp_assert(nrow(winners) == length(unique(manifest$fit_config_id)),
              "Existing Stage 6B start-winner table is malformed.")
    return(winners)
  }
  pp_assert(!dir.exists(phase_dir),
            "Incomplete Stage 6B start-selection directory exists.")
  terminal <- collect_terminals(manifest)$rows
  winners <- pp_select_g12_winners(terminal)
  pp_assert(nrow(winners) == length(unique(manifest$fit_config_id)),
            "Stage 6B start selection did not yield one winner per config.")
  dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  pp_write_csv(terminal, file.path(phase_dir, "ALL_ENDPOINTS.csv"))
  pp_write_csv(winners, file.path(phase_dir, "CONFIG_WINNERS.csv"))
  hash <- pp_sha256(file.path(phase_dir, "CONFIG_WINNERS.csv"))
  pp_write_lines(c(
    "status=START_SELECTION_FROZEN", paste0("phase=", phase),
    paste0("completed_utc=", pp_iso_time()),
    paste0("endpoints=", nrow(terminal)),
    paste0("config_winners=", nrow(winners)),
    paste0("winner_table_sha256=", hash),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_configuration",
    "heldout_score_read_for_start_selection=FALSE",
    "selection_used_truth=FALSE"
  ), marker)
  winners
}

read_winner_losses <- function(winners) {
  values <- lapply(seq_len(nrow(winners)), function(index) {
    fit_id <- winners$fit_id[[index]]
    path <- file.path(output_root, "fits", fit_id,
                      "HOLDOUT_SUBJECT_TIME_ALL.csv")
    rows <- pp_read_csv(path)
    rows$fit_config_id <- winners$fit_config_id[[index]]
    rows$winner_fit_id <- fit_id
    rows
  })
  do.call(rbind, values)
}

freeze_factor_dimension <- function(registration, factor_winners) {
  selection_root <- file.path(output_root, "truth_free_selection")
  path <- factor_selection_path()
  if (file.exists(path)) return(pp_read_csv(path))
  metadata <- registration$configs[
    registration$configs$phase == "factor_screen", , drop = FALSE
  ]
  result <- s6bp_select_one_se(
    read_winner_losses(factor_winners), metadata, "factor_screen"
  )
  pp_write_csv(result$summary,
               file.path(selection_root, "FACTOR_DIMENSION_CANDIDATES.csv"))
  selected_metadata <- metadata[
    metadata$fit_config_id == result$selected_config_id, , drop = FALSE
  ]
  selection <- data.frame(
    phase = "factor_screen",
    selected_fit_config_id = result$selected_config_id,
    minimum_loss_fit_config_id = result$minimum_loss_config_id,
    selected_factor_count = as.integer(selected_metadata$fit_L_f[[1L]]),
    primary_score = "heldout_nrmse",
    one_se_scale = "mean_subject_time_mse",
    one_se_threshold_mse = result$threshold_mse,
    selection_used_truth = FALSE,
    heldout_score_used_for_start_selection = FALSE,
    stringsAsFactors = FALSE
  )
  pp_write_csv(selection, path)
  hash <- pp_sha256(path)
  pp_write_lines(c(
    "status=FACTOR_DIMENSION_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("selected_factor_count=", selection$selected_factor_count),
    paste0("selection_sha256=", hash),
    "selection_rule=heldout_nrmse_one_standard_error_simplest",
    "selection_used_truth=FALSE"
  ), file.path(selection_root, "FACTOR_DIMENSION_SELECTION_FROZEN.txt"))
  selection
}

freeze_fpca_dimension <- function(registration, selected_L,
                                  factor_winners, fpca_winners) {
  selection_root <- file.path(output_root, "truth_free_selection")
  final_path <- file.path(selection_root,
                          "FINAL_TRUTH_FREE_DIMENSION_SELECTION.csv")
  if (file.exists(final_path)) return(pp_read_csv(final_path))
  base_id <- paste0("factor_L", selected_L, "_M2")
  base_winner <- factor_winners[
    factor_winners$fit_config_id == base_id, , drop = FALSE
  ]
  pp_assert(nrow(base_winner) == 1L && nrow(fpca_winners) == 2L,
            "The Stage 6B FPCA winner set must contain M=2,3,4.")
  winners <- rbind(base_winner, fpca_winners)
  base_metadata <- registration$configs[
    registration$configs$fit_config_id == base_id, , drop = FALSE
  ]
  base_metadata$phase <- "fpca_screen"
  base_metadata$complexity_value <- 2L
  conditional_metadata <- registration$configs[
    registration$configs$phase == "fpca_screen" &
      as.integer(registration$configs$fit_L_f) == selected_L,
    , drop = FALSE
  ]
  metadata <- rbind(base_metadata, conditional_metadata)
  result <- s6bp_select_one_se(
    read_winner_losses(winners), metadata, "fpca_screen"
  )
  pp_write_csv(result$summary,
               file.path(selection_root, "FPCA_DIMENSION_CANDIDATES.csv"))
  selected_metadata <- metadata[
    metadata$fit_config_id == result$selected_config_id, , drop = FALSE
  ]
  selected_M <- as.integer(selected_metadata$complexity_value[[1L]])
  selection <- data.frame(
    selected_factor_count = selected_L,
    selected_fpca_cap = selected_M,
    selected_fit_config_id = result$selected_config_id,
    minimum_loss_fit_config_id = result$minimum_loss_config_id,
    primary_score = "heldout_nrmse",
    one_se_scale = "mean_subject_time_mse",
    one_se_threshold_mse = result$threshold_mse,
    within_configuration_start_rule =
      "maximum_eligible_ordinary_T1_ELBO",
    cross_dimension_rule =
      "heldout_NRMSE_equivalent_MSE_one_standard_error_simplest",
    selection_used_truth = FALSE,
    new_data_generated = FALSE,
    final_full_data_refit_started = FALSE,
    continuation_started = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
  pp_write_csv(selection, final_path)
  hash <- pp_sha256(final_path)
  pp_write_lines(c(
    "status=FINAL_TRUTH_FREE_DIMENSION_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("selected_factor_count=", selected_L),
    paste0("selected_fpca_cap=", selected_M),
    paste0("selection_sha256=", hash), "selection_used_truth=FALSE",
    "truth_access_authorized=FALSE", "final_full_data_refit_started=FALSE",
    "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(selection_root,
               "FINAL_TRUTH_FREE_DIMENSION_SELECTION_FROZEN.txt"))
  selection
}

write_status <- function(state, phase, total, detail = "") {
  complete <- if (dir.exists(file.path(output_root, "fits"))) {
    length(list.files(file.path(output_root, "fits"),
                      pattern = "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  } else 0L
  errors <- if (dir.exists(file.path(output_root, "fits"))) {
    length(list.files(file.path(output_root, "fits"),
                      pattern = "^FIT_ERROR[.]txt$", recursive = TRUE))
  } else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, phase = phase,
    detail = detail, phase_registered_fits = total,
    all_complete_fits = complete, all_error_fits = errors,
    outer_workers = PP_WORKERS, per_fit_n_cpus = 1L,
    truth_read_by_fit_or_selection = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

run_phase <- function(manifest, phase) {
  marker <- file.path(selection_phase_dir(phase),
                      "START_SELECTION_FROZEN.txt")
  if (file.exists(marker)) return(freeze_start_winners(manifest, phase))
  pp_assert(.Platform$OS.type == "unix",
            "The Stage 6B full pilot supervisor requires Unix.")
  log_dir <- file.path(output_root, "supervisor_logs", phase)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  write_status("RUNNING", phase, nrow(manifest),
               "registered fits; no selective rerun")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    started <- Sys.time()
    status <- system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=fit",
                   paste0("--fit-id=", fit_id),
                   paste0("--output-root=", output_root)),
      stdout = stdout, stderr = stderr
    )
    ended <- Sys.time()
    pp_write_csv(data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      exit_status = as.integer(status), stringsAsFactors = FALSE
    ), file.path(log_dir, paste0(key, ".exit.csv")))
    list(fit_id = fit_id, status = as.integer(status))
  }
  results <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = PP_WORKERS,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    write_status("ERROR", phase, nrow(manifest),
                 paste0("retained_no_rerun=", paste(failed, collapse = ";")))
    pp_abort("Stage 6B fit errors retained without selective rerun: ",
             paste(failed, collapse = ", "))
  }
  winners <- freeze_start_winners(manifest, phase)
  write_status("PHASE_COMPLETE", phase, nrow(manifest),
               paste0("config_winners=", nrow(winners)))
  winners
}

summarize_pilot <- function(registration, selected_L) {
  factor_manifest <- s6bp_phase_manifest(registration$fits, "factor_screen")
  fpca_manifest <- s6bp_phase_manifest(registration$fits, "fpca_screen",
                                       selected_L)
  manifest <- rbind(factor_manifest, fpca_manifest)
  terminal <- collect_terminals(manifest)$rows
  runtime <- do.call(rbind, lapply(
    split(terminal, terminal$fit_config_id), function(rows) data.frame(
      fit_config_id = rows$fit_config_id[[1L]], fits = nrow(rows),
      elapsed_median_seconds = stats::median(rows$elapsed_seconds),
      elapsed_max_seconds = max(rows$elapsed_seconds),
      peak_memory_max_bytes = max(rows$peak_memory_bytes),
      warning_total = sum(rows$warning_count),
      objective_eligible = sum(rows$objective_eligible),
      practical_converged = sum(rows$practical_converged),
      stringsAsFactors = FALSE
    )
  ))
  pp_write_csv(runtime, file.path(output_root,
                                  "RUNTIME_AND_MEMORY_SUMMARY.csv"))
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  sizes <- file.info(files)$size
  pp_write_csv(data.frame(
    generated_at_utc = pp_iso_time(), files = length(files),
    disk_bytes = sum(sizes[is.finite(sizes)], na.rm = TRUE),
    registered_conditional_fit_rows = nrow(registration$fits),
    executed_fits = nrow(manifest), complete_fits = nrow(terminal),
    objective_eligible = sum(terminal$objective_eligible),
    worker_count = PP_WORKERS, n_cpus_per_fit = 1L,
    new_data_generated = FALSE, truth_read = FALSE,
    continuation_started = FALSE, final_full_data_refit_started = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "RESOURCE_SUMMARY.csv"))
  invisible(runtime)
}

run_pilot <- function() {
  require_output_root(FALSE)
  registration <- verify_environment()
  load_prepared()
  pp_assert(!file.exists(file.path(output_root,
                                   "TRUTH_UNSEAL_AUTHORIZATION.txt")),
            "Refusing Stage 6B execution after truth-access authorization.")
  factor_manifest <- s6bp_phase_manifest(registration$fits, "factor_screen")
  factor_winners <- run_phase(factor_manifest, "factor_screen")
  factor_selection <- freeze_factor_dimension(registration, factor_winners)
  selected_L <- as.integer(factor_selection$selected_factor_count[[1L]])
  fpca_manifest <- s6bp_phase_manifest(registration$fits, "fpca_screen",
                                       selected_L)
  fpca_winners <- run_phase(fpca_manifest, "fpca_screen")
  final_selection <- freeze_fpca_dimension(
    registration, selected_L, factor_winners, fpca_winners
  )
  summarize_pilot(registration, selected_L)
  pp_write_lines(c(
    "status=STAGE6B_PREDICTIVE_PILOT_COMPLETE",
    paste0("completed_utc=", pp_iso_time()), "executed_fits=20",
    paste0("selected_factor_count=",
           final_selection$selected_factor_count[[1L]]),
    paste0("selected_fpca_cap=",
           final_selection$selected_fpca_cap[[1L]]),
    "selection_used_truth=FALSE", "truth_access_authorized=FALSE",
    "new_data_generated=FALSE", "continuation_started=FALSE",
    "final_full_data_refit_started=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "STAGE6B_PILOT_COMPLETE.txt"))
  write_status("COMPLETE_TRUTH_REMAINS_SEALED", "all", 20L,
               "factor and FPCA dimension selections frozen")
  cat("STAGE6B_PILOT_COMPLETE selected_L=",
      final_selection$selected_factor_count[[1L]], " selected_M=",
      final_selection$selected_fpca_cap[[1L]], "\n", sep = "")
  invisible(final_selection)
}

if (action == "check") {
  registration <- verify_environment()
  cat("STAGE6B_CHECK_PASS conditional_fits=", nrow(registration$fits),
      " factor_fits=12 max_executed=20 workers=", PP_WORKERS, "\n", sep = "")
} else if (action == "prepare") {
  prepare_pilot()
} else if (action == "fit") {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required for action=fit.")
  run_fit(fit_id)
} else if (action == "run") {
  run_pilot()
} else if (action == "summarize") {
  registration <- verify_environment()
  summarize_pilot(registration, selected_factor_count())
} else {
  pp_abort("Unknown Stage 6B action: ", action)
}
