#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}

experiment_parent <- dirname(script_root)
stage6a_root <- file.path(experiment_parent,
                          "g12-dimension-overspec-stage6a-20260831")
stage6b_root <- file.path(experiment_parent,
                          "g12-truth-free-dimension-selection-stage6b-20260901")
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)
source(file.path(stage6b_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = FALSE)
source(file.path(stage6b_root,
                 "stage6b_predictive_pilot_tools_20260902_v1.R"),
       local = FALSE)
source(file.path(script_root, "stage6d_tools_20260906_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6D_INDEPENDENT_DIMENSION_SELECTION_V1_20260906"
S6D_WORKERS <- as.integer(pp_env("STAGE6D_WORKERS", "7"))
pp_assert(length(S6D_WORKERS) == 1L && is.finite(S6D_WORKERS) &&
            S6D_WORKERS %in% 1:7,
          "STAGE6D_WORKERS must be an integer from 1 to 7.")
stage6d_library <- pp_env("STAGE6D_R_LIB")
if (nzchar(stage6d_library)) Sys.setenv(STAGE6A_R_LIB = stage6d_library)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6D_OUTPUT_ROOT")

data_path <- file.path(script_root, "DATA_SEEDS_STAGE6D_V1.csv")
config_path <- file.path(script_root, "FIT_CONFIGS_STAGE6D_V1.csv")
seed_path <- file.path(script_root, "START_SEEDS_STAGE6D_V1.csv")

require_output_root <- function(create = TRUE) {
  pp_assert(nzchar(output_root), "--output-root is required for this action.")
  if (create && !dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  output_root
}

load_registration <- function() {
  data <- s6d_read_csv(data_path)
  configs <- s6d_read_csv(config_path)
  seeds <- s6d_read_csv(seed_path)
  s6d_validate_registration(data, configs, seeds)
  list(data = data, configs = configs, seeds = seeds,
       fits = s6d_build_fit_manifest(data, configs, seeds))
}

verify_environment <- function() {
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  registration <- load_registration()
  pp_load_package()
  pp_source_runtime(FALSE)
  pp_assert(nrow(registration$fits) == 180L &&
              length(unique(registration$fits$selection_stratum_id)) == 15L &&
              all(table(registration$fits$selection_stratum_id) == 12L),
            "Stage 6D registration integrity failed.")
  registration
}

data_paths <- function(data_id) list(
  observation = file.path(output_root, "data", data_id,
                          "observation_bundle.rds"),
  truth = file.path(output_root, "data", data_id, "sealed_truth",
                    "truth_bundle.rds"),
  seal = file.path(output_root, "data", data_id,
                   "DATA_SEAL_COMPLETE.txt"),
  train = file.path(output_root, "data", data_id, "split",
                    "train_observation.rds"),
  holdout = file.path(output_root, "data", data_id, "split",
                      "holdout_plan.rds"),
  split_summary = file.path(output_root, "data", data_id, "split",
                            "HOLDOUT_SPLIT_SUMMARY.csv")
)

truth_authorization_path <- function() file.path(
  output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"
)

assert_truth_locked <- function() {
  pp_assert(!file.exists(truth_authorization_path()),
            "Refusing Stage 6D fitting/selection after truth authorization.")
  invisible(TRUE)
}

write_source_binding <- function() {
  relative <- c(
    "DATA_SEEDS_STAGE6D_V1.csv", "FIT_CONFIGS_STAGE6D_V1.csv",
    "START_SEEDS_STAGE6D_V1.csv",
    "PROTOCOL_G12_STAGE6D_INDEPENDENT_DIMENSION_SELECTION_V1_20260906.md",
    "stage6d_tools_20260906_v1.R", "run_stage6d_20260906_v1.R",
    file.path("..", "g12-dimension-overspec-stage6a-20260831",
              "stage6a_common_20260831_v1.R"),
    file.path("..", "g12-truth-free-dimension-selection-stage6b-20260901",
              "stage6b_truth_free_tools_20260901_v1.R"),
    file.path("..", "g12-truth-free-dimension-selection-stage6b-20260901",
              "stage6b_predictive_pilot_tools_20260902_v1.R")
  )
  absolute <- file.path(script_root, relative)
  pp_assert(all(file.exists(absolute)), "A Stage 6D bound source is missing.")
  info <- file.info(absolute)
  pp_write_csv(data.frame(
    experiment_id = PP_EXPERIMENT_ID,
    relative_to_stage6d_script_root = gsub("\\\\", "/", relative),
    bytes = as.numeric(info$size),
    sha256 = vapply(absolute, pp_sha256, character(1L)),
    stringsAsFactors = FALSE
  ), file.path(output_root, "SOURCE_BINDING_RUNTIME.csv"))
  git_commit <- tryCatch(
    system2("git", c("-C", shQuote(PP_REPO_ROOT), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(condition) character()
  )
  if (length(git_commit) == 1L && grepl("^[0-9a-f]{40}$", git_commit)) {
    pp_write_lines(paste0("git_commit=", git_commit),
                   file.path(output_root, "SOURCE_GIT_COMMIT.txt"))
  }
  invisible(TRUE)
}

prepare_experiment <- function(smoke = FALSE) {
  require_output_root(TRUE)
  assert_truth_locked()
  registration <- verify_environment()
  marker <- file.path(output_root, if (smoke) {
    "STAGE6D_MICRO_PREPARED.txt"
  } else "STAGE6D_PREPARED.txt")
  if (file.exists(marker)) {
    expected_ids <- if (smoke) registration$data$data_id[[1L]] else
      registration$data$data_id
    expected <- unlist(lapply(expected_ids, data_paths),
                       use.names = FALSE)
    pp_assert(all(file.exists(expected)),
              "Existing Stage 6D preparation is incomplete.")
    return(invisible(registration))
  }
  pp_assert(!dir.exists(file.path(output_root, "fits")) &&
              !dir.exists(file.path(output_root, "truth_free_selection")),
            "Refusing to prepare over an existing Stage 6D fit/selection tree.")
  pp_write_csv(registration$data,
               file.path(output_root, "DATA_MANIFEST.csv"))
  pp_write_csv(registration$configs,
               file.path(output_root, "FIT_CONFIGS.csv"))
  pp_write_csv(registration$seeds,
               file.path(output_root, "START_SEEDS.csv"))
  pp_write_csv(registration$fits,
               file.path(output_root, "FIT_MANIFEST.csv"))
  write_source_binding()

  data_rows <- if (smoke) registration$data[1L, , drop = FALSE] else
    registration$data
  for (index in seq_len(nrow(data_rows))) {
    row <- data_rows[index, , drop = FALSE]
    pp_generate_one(row, output_root, smoke = smoke)
    observation <- pp_load_observation(row$data_id[[1L]], output_root)
    split <- s6b_make_middle_time_holdout(
      observation, min_train_times = if (smoke) 3L else 5L
    )
    expected_units <- if (smoke) 8L else 60L
    pp_assert(split$holdout_plan$subject_time_count == expected_units &&
                identical(split$holdout_plan$truth_used, FALSE) &&
                !length(s6b_forbidden_observation_names(
                  split$train_observation
                )),
              "Stage 6D holdout split contract failed.")
    paths <- data_paths(row$data_id[[1L]])
    pp_save_rds(split$train_observation, paths$train)
    pp_save_rds(split$holdout_plan, paths$holdout)
    records <- split$holdout_plan$records
    pp_write_csv(data.frame(
      data_id = row$data_id[[1L]], split_id = split$holdout_plan$split_id,
      study = vapply(records, `[[`, integer(1L), "study"),
      subject = vapply(records, `[[`, integer(1L), "subject"),
      original_index = vapply(records, `[[`, integer(1L), "original_index"),
      holdout_time = vapply(records, `[[`, numeric(1L), "time"),
      train_time_count = vapply(records, `[[`, integer(1L),
                                "train_time_count"),
      variables = vapply(records, `[[`, integer(1L), "p"),
      truth_used = FALSE, stringsAsFactors = FALSE
    ), paths$split_summary)
  }
  pp_write_lines(c(
    paste0("status=", if (smoke) "STAGE6D_MICRO_PREPARED" else
             "STAGE6D_PREPARED"),
    paste0("prepared_utc=", pp_iso_time()),
    paste0("data_sets=", nrow(data_rows)),
    paste0("registered_full_fit_rows=", nrow(registration$fits)),
    paste0("planned_executed_fit_rows=", if (smoke) 5L else 180L),
    "truth_read_by_fit_or_selection=FALSE",
    "truth_access_authorized=FALSE", "continuation_authorized=FALSE",
    "formal_paper_mc_result=FALSE"
  ), marker)
  invisible(registration)
}

load_split <- function(data_id) {
  paths <- data_paths(data_id)
  pp_assert(file.exists(paths$train) && file.exists(paths$holdout),
            paste0("Prepared Stage 6D split is missing: ", data_id))
  train <- readRDS(paths$train)
  holdout <- readRDS(paths$holdout)
  pp_assert(!length(s6b_forbidden_observation_names(train)) &&
              identical(holdout$truth_used, FALSE),
            paste0("Stage 6D split is invalid or truth-bearing: ", data_id))
  list(train = train, holdout = holdout)
}

write_fit_error <- function(fit_dir, fit_id, message) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
    paste0("message=", message), paste0("ended_utc=", pp_iso_time()),
    "truth_read=FALSE"
  ), file.path(fit_dir, "FIT_ERROR.txt"))
}

run_fit <- function(fit_id) {
  require_output_root(FALSE)
  assert_truth_locked()
  registration <- verify_environment()
  row <- registration$fits[registration$fits$fit_id == fit_id,
                           , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Stage 6D fit_id is not registered.")
  split <- load_split(row$data_id[[1L]])
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")) &&
                file.exists(file.path(fit_dir, "fit.rds")) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Pre-existing fit is incomplete/error and cannot be ",
                     "selectively rerun: ", fit_id))
    cat("FIT_ALREADY_COMPLETE fit_id=", fit_id, "\n", sep = "")
    return(invisible(0L))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- pp_fit_g12(split$train, row, smoke = FALSE)
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  result$captured$warnings <- s6d_set_warning_phase(
    result$captured$warnings, "stage6d_g12_fixed400"
  )
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6D_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$holdout_split_id <- split$holdout$split_id
  record$heldout_subject_time_vectors <- split$holdout$subject_time_count
  record$primary_dimension_score <- "heldout_nrmse"
  record$holdout_scoring_complete <- FALSE
  record$truth_used_for_holdout_scoring <- FALSE
  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(), call = character(),
               phase = character(), stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
    write_fit_error(fit_dir, fit_id, record$error$message %||% "unknown")
    pp_abort("Stage 6D fit failed and was retained: ", fit_id)
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
    pp_abort("Stage 6D holdout scoring failed and was retained: ", fit_id)
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
    "continuation_used=FALSE", paste0("ended_utc=", pp_iso_time())
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
            "Stage 6D has missing terminal, completion or score files.")
  records <- lapply(paths, readRDS)
  pp_assert(all(vapply(records, function(record) {
    identical(record$terminal_schema, "G12_STAGE6D_TERMINAL_V1") &&
      isTRUE(record$holdout_scoring_complete) &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$truth_used_for_holdout_scoring)
  }, logical(1L))), "Stage 6D terminal contract failed.")
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  pp_assert(identical(rows$fit_id, manifest$fit_id) &&
              all(rows$terminal_status == "fixed_400_complete") &&
              all(rows$objective_eligible) &&
              !any(rows$truth_used_for_fit_stopping_or_selection),
            "Stage 6D contains an ineligible or truth-bearing endpoint.")
  list(records = records, rows = rows)
}

write_status <- function(state, detail = "") {
  complete <- if (dir.exists(file.path(output_root, "fits"))) {
    length(list.files(file.path(output_root, "fits"),
                      pattern = "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  } else 0L
  errors <- if (dir.exists(file.path(output_root, "fits"))) {
    length(list.files(file.path(output_root, "fits"),
                      pattern = "^FIT_ERROR[.]txt$", recursive = TRUE))
  } else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_fits = 180L, complete_fits = complete, error_fits = errors,
    outer_workers = S6D_WORKERS, per_fit_n_cpus = 1L,
    truth_read_by_fit_or_selection = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

read_winner_losses <- function(winners) {
  values <- lapply(seq_len(nrow(winners)), function(index) {
    fit_id <- winners$fit_id[[index]]
    rows <- pp_read_csv(file.path(output_root, "fits", fit_id,
                                  "HOLDOUT_SUBJECT_TIME_ALL.csv"))
    rows$data_id <- winners$data_id[[index]]
    rows$fit_config_id <- winners$fit_config_id[[index]]
    rows$winner_fit_id <- fit_id
    rows
  })
  do.call(rbind, values)
}

freeze_truth_free_selection <- function(registration, terminals) {
  selection_root <- file.path(output_root, "truth_free_selection")
  marker <- file.path(selection_root,
                      "TRUTH_FREE_SELECTION_FROZEN.txt")
  if (file.exists(marker)) return(invisible(TRUE))
  pp_assert(!dir.exists(selection_root),
            "Incomplete Stage 6D truth-free selection directory exists.")
  winners <- pp_select_g12_winners(terminals)
  pp_assert(nrow(winners) == 15L &&
              all(table(winners$data_id) == 5L),
            "Stage 6D did not produce 15 configuration winners.")
  dir.create(selection_root, recursive = TRUE, showWarnings = FALSE,
             mode = "0700")
  pp_write_csv(terminals, file.path(selection_root, "ALL_ENDPOINTS.csv"))
  pp_write_csv(winners, file.path(selection_root, "CONFIG_WINNERS.csv"))
  winner_hash <- pp_sha256(file.path(selection_root, "CONFIG_WINNERS.csv"))
  pp_write_lines(c(
    "status=START_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()), "endpoints=180",
    "configuration_winners=15",
    paste0("winner_table_sha256=", winner_hash),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_configuration",
    "heldout_score_read_for_start_selection=FALSE",
    "selection_used_truth=FALSE"
  ), file.path(selection_root, "START_SELECTION_FROZEN.txt"))

  losses <- read_winner_losses(winners)
  selections <- lapply(registration$data$data_id, function(data_id) {
    result <- s6d_select_axes(losses, registration$configs, data_id)
    data_dir <- file.path(selection_root, data_id)
    pp_write_csv(result$factor$summary,
                 file.path(data_dir, "FACTOR_CANDIDATES.csv"))
    pp_write_csv(result$fpca$summary,
                 file.path(data_dir, "FPCA_CANDIDATES.csv"))
    selected <- data.frame(
      data_id = data_id,
      selected_factor_config_id = result$factor$selected_config_id,
      minimum_factor_loss_config_id = result$factor$minimum_loss_config_id,
      selected_factor_count = as.integer(sub(".*_L([0-9]+)_M.*", "\\1",
                                               result$factor$selected_config_id)),
      selected_fpca_config_id = result$fpca$selected_config_id,
      minimum_fpca_loss_config_id = result$fpca$minimum_loss_config_id,
      selected_fpca_cap = as.integer(sub(".*_M([0-9]+)$", "\\1",
                                          result$fpca$selected_config_id)),
      primary_score = "heldout_nrmse",
      one_se_scale = "mean_subject_time_mse",
      within_configuration_start_rule =
        "maximum_eligible_ordinary_T1_ELBO",
      cross_dimension_rule =
        "heldout_MSE_one_standard_error_simplest_separately_by_axis",
      selection_used_truth = FALSE,
      continuation_started = FALSE,
      stringsAsFactors = FALSE
    )
    pp_write_csv(selected, file.path(data_dir, "DIMENSION_SELECTION.csv"))
    selected
  })
  selections <- do.call(rbind, selections)
  pp_write_csv(selections,
               file.path(selection_root, "DIMENSION_SELECTIONS.csv"))
  pp_write_csv(data.frame(
    selected_factor_count = 1:3,
    data_sets = vapply(1:3, function(value) {
      sum(selections$selected_factor_count == value)
    }, integer(1L)), stringsAsFactors = FALSE
  ), file.path(selection_root, "FACTOR_SELECTION_COUNTS.csv"))
  pp_write_csv(data.frame(
    selected_fpca_cap = 2:4,
    data_sets = vapply(2:4, function(value) {
      sum(selections$selected_fpca_cap == value)
    }, integer(1L)), stringsAsFactors = FALSE
  ), file.path(selection_root, "FPCA_SELECTION_COUNTS.csv"))
  selection_hash <- pp_sha256(
    file.path(selection_root, "DIMENSION_SELECTIONS.csv")
  )
  pp_write_lines(c(
    "status=TRUTH_FREE_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("configuration_winner_sha256=", winner_hash),
    paste0("dimension_selection_sha256=", selection_hash),
    "data_sets=3", "factor_axis_candidates_per_data=3",
    "fpca_axis_candidates_per_data=3", "selection_used_truth=FALSE",
    "truth_access_authorized=FALSE", "continuation_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), marker)
  invisible(selections)
}

summarize_resources <- function(registration, terminals) {
  groups <- split(terminals, terminals$fit_config_id)
  runtime <- do.call(rbind, lapply(groups, function(rows) data.frame(
    fit_config_id = rows$fit_config_id[[1L]], fits = nrow(rows),
    elapsed_median_seconds = stats::median(rows$elapsed_seconds),
    elapsed_max_seconds = max(rows$elapsed_seconds),
    peak_memory_max_bytes = max(rows$peak_memory_bytes),
    warning_total = sum(rows$warning_count),
    objective_eligible = sum(rows$objective_eligible),
    practical_converged = sum(rows$practical_converged),
    stringsAsFactors = FALSE
  )))
  pp_write_csv(runtime,
               file.path(output_root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  sizes <- file.info(files)$size
  pp_write_csv(data.frame(
    generated_at_utc = pp_iso_time(), files = length(files),
    disk_bytes = sum(sizes[is.finite(sizes)], na.rm = TRUE),
    data_sets = nrow(registration$data), executed_fits = nrow(terminals),
    objective_eligible = sum(terminals$objective_eligible),
    worker_count = S6D_WORKERS, n_cpus_per_fit = 1L,
    truth_read = FALSE, continuation_started = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "RESOURCE_SUMMARY.csv"))
  invisible(runtime)
}

run_experiment <- function() {
  require_output_root(FALSE)
  assert_truth_locked()
  registration <- verify_environment()
  pp_assert(file.exists(file.path(output_root, "STAGE6D_PREPARED.txt")),
            "Prepare Stage 6D before running fits.")
  pp_assert(.Platform$OS.type == "unix",
            "The Stage 6D full supervisor requires Unix.")
  log_dir <- file.path(output_root, "supervisor_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  write_status("RUNNING", "180 registered fits; no selective rerun")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    started <- Sys.time()
    status <- tryCatch(system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=fit",
                   paste0("--fit-id=", fit_id),
                   paste0("--output-root=", output_root)),
      stdout = stdout, stderr = stderr
    ), error = function(condition) 999L)
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
    registration$fits$fit_id, run_one, mc.cores = S6D_WORKERS,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    write_status("ERROR",
                 paste0("retained_no_rerun=", paste(failed, collapse = ";")))
    pp_abort("Stage 6D fit errors retained without selective rerun: ",
             paste(failed, collapse = ", "))
  }
  terminal <- collect_terminals(registration$fits)$rows
  freeze_truth_free_selection(registration, terminal)
  summarize_resources(registration, terminal)
  pp_write_lines(c(
    "status=STAGE6D_COMPLETE_TRUTH_REMAINS_SEALED",
    paste0("completed_utc=", pp_iso_time()), "data_sets=3",
    "executed_fits=180", "configuration_winners=15",
    "selection_used_truth=FALSE", "truth_access_authorized=FALSE",
    "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "STAGE6D_COMPLETE.txt"))
  write_status("COMPLETE_TRUTH_REMAINS_SEALED",
               "15 start winners and six axis selections frozen")
  cat("STAGE6D_COMPLETE fits=180 truth_remains_sealed=TRUE\n")
  invisible(TRUE)
}

run_micro_smoke <- function() {
  require_output_root(TRUE)
  registration <- prepare_experiment(smoke = TRUE)
  data_id <- registration$data$data_id[[1L]]
  split <- load_split(data_id)
  manifest <- registration$fits[
    registration$fits$data_id == data_id &
      registration$fits$seed_index == 1L, , drop = FALSE
  ]
  pp_assert(nrow(manifest) == 5L, "Micro smoke must contain five fits.")
  records <- vector("list", nrow(manifest))
  losses <- vector("list", nrow(manifest))
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    result <- pp_fit_g12(split$train, row, smoke = TRUE)
    pp_assert(!is.null(result$fit) && isTRUE(result$objective$eligible),
              paste0("Micro fit failed: ", row$fit_id[[1L]]))
    record <- pp_make_terminal_record(result, row, smoke = TRUE)
    record$holdout_scoring_complete <- TRUE
    record$truth_used_for_holdout_scoring <- FALSE
    records[[index]] <- record
    score <- s6b_score_holdout(result$fit, split$holdout, "all")
    rows <- score$subject_time
    rows$data_id <- data_id
    rows$fit_config_id <- row$fit_config_id[[1L]]
    rows$winner_fit_id <- row$fit_id[[1L]]
    losses[[index]] <- rows
  }
  terminals <- do.call(rbind, lapply(records, pp_terminal_row))
  winners <- pp_select_g12_winners(terminals)
  pp_assert(nrow(winners) == 5L,
            "Micro smoke did not produce five start winners.")
  axes <- s6d_select_axes(do.call(rbind, losses), registration$configs,
                          data_id)
  pp_write_csv(terminals, file.path(output_root, "MICRO_TERMINALS.csv"))
  pp_write_csv(axes$factor$summary,
               file.path(output_root, "MICRO_FACTOR_SELECTION.csv"))
  pp_write_csv(axes$fpca$summary,
               file.path(output_root, "MICRO_FPCA_SELECTION.csv"))
  pp_write_lines(c(
    "status=PASS", paste0("completed_utc=", pp_iso_time()),
    "micro_data_sets=1", "micro_fits=5", "configuration_winners=5",
    "factor_axis_selection=PASS", "fpca_axis_selection=PASS",
    "truth_read_by_fit_or_selection=FALSE"
  ), file.path(output_root, "STAGE6D_MICRO_SMOKE_PASS.txt"))
  cat("STAGE6D_MICRO_SMOKE_PASS fits=5 truth_read=FALSE\n")
  invisible(TRUE)
}

if (action == "check") {
  registration <- verify_environment()
  cat("STAGE6D_CHECK_PASS data=", nrow(registration$data),
      " configs=", nrow(registration$configs),
      " starts=", nrow(registration$seeds),
      " fits=", nrow(registration$fits),
      " workers=", S6D_WORKERS, " n_cpus_per_fit=1\n", sep = "")
} else if (action == "prepare") {
  prepare_experiment(FALSE)
  cat("STAGE6D_PREPARED data=3 fits=180 truth_sealed=TRUE\n")
} else if (action == "fit") {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required for action=fit.")
  run_fit(fit_id)
} else if (action == "run") {
  run_experiment()
} else if (action == "micro-smoke") {
  run_micro_smoke()
} else {
  pp_abort("Unknown Stage 6D action: ", action)
}
