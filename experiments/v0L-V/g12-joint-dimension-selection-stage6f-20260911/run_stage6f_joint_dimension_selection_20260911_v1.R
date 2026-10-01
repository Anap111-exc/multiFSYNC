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
stage6a_root <- file.path(experiments_root,
                          "g12-dimension-overspec-stage6a-20260831")
stage6b_code_root <- file.path(
  experiments_root, "g12-truth-free-dimension-selection-stage6b-20260901"
)
stage6eb_root <- file.path(
  experiments_root, "g12-predictive-logscore-feasibility-stage6e-b-20260910"
)
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)
source(file.path(stage6b_code_root,
                 "stage6b_truth_free_tools_20260901_v1.R"),
       local = FALSE)
source(file.path(stage6b_code_root,
                 "stage6b_predictive_pilot_tools_20260902_v1.R"),
       local = FALSE)
source(file.path(stage6eb_root,
                 "stage6eb_logscore_tools_20260910_v1.R"),
       local = FALSE)
source(file.path(script_root,
                 "stage6f_joint_selection_tools_20260911_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6F_JOINT_DIMENSION_SELECTION_V1_20260911"
stage6f_library <- pp_env("STAGE6F_R_LIB")
if (nzchar(stage6f_library)) Sys.setenv(STAGE6A_R_LIB = stage6f_library)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6F_OUTPUT_ROOT")
stage6b_output_root <- cli$`stage6b-root` %||%
  pp_env("STAGE6F_STAGE6B_ROOT")
observation_path <- cli$`observation-bundle` %||%
  pp_env("STAGE6F_OBSERVATION_BUNDLE")

config_path <- file.path(script_root, "FIT_CONFIGS_STAGE6F_V1.csv")
seed_path <- file.path(script_root, "START_SEEDS_STAGE6F_V1.csv")
split_path <- file.path(script_root, "HOLDOUT_SPLITS_STAGE6F_V1.csv")

load_registration <- function() {
  configs <- s6f_read_csv(config_path)
  seeds <- s6f_read_csv(seed_path)
  splits <- s6f_read_csv(split_path)
  s6f_validate_registration(configs, seeds, splits)
  list(
    configs = configs, seeds = seeds, splits = splits,
    fits = s6f_build_manifest(configs, seeds, splits)
  )
}

verify_environment <- function(load_runtime = TRUE) {
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  registration <- load_registration()
  if (load_runtime) {
    pp_load_package()
    pp_source_runtime(FALSE)
  }
  pp_assert(nrow(registration$fits) == 32L &&
              sum(registration$fits$reuse_stage6b) == 12L &&
              sum(registration$fits$execute_new_fit) == 20L,
            "Stage 6F registration integrity failed.")
  registration
}

require_roots <- function(create_output = FALSE) {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pp_assert(nzchar(stage6b_output_root) && dir.exists(stage6b_output_root),
            "--stage6b-root must identify the corrected Stage 6B v2 output.")
  pp_assert(nzchar(observation_path) && file.exists(observation_path),
            "--observation-bundle must identify the Stage 6A observation bundle.")
  if (create_output && !dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  invisible(TRUE)
}

prepared_paths <- function() list(
  manifest = file.path(output_root, "FIT_MANIFEST.csv"),
  middle_train = file.path(output_root, "splits", "split1_middle",
                           "train_observation.rds"),
  middle_holdout = file.path(output_root, "splits", "split1_middle",
                             "holdout_plan.rds"),
  upper_train = file.path(output_root, "splits", "split2_upper",
                          "train_observation.rds"),
  upper_holdout = file.path(output_root, "splits", "split2_upper",
                            "holdout_plan.rds"),
  marker = file.path(output_root, "STAGE6F_PREPARED.txt")
)

s6f_link_or_copy <- function(source, destination) {
  pp_assert(file.exists(source), paste0("Missing source file: ", source))
  pp_assert(!file.exists(destination), paste0(
    "Refusing to overwrite existing destination: ", destination
  ))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE,
             mode = "0700")
  linked <- isTRUE(file.link(source, destination))
  if (!linked) {
    copied <- file.copy(source, destination, overwrite = FALSE,
                        copy.mode = TRUE, copy.date = TRUE)
    pp_assert(isTRUE(copied), paste0("Could not copy source: ", source))
  }
  list(method = if (linked) "hardlink" else "copy",
       sha256 = pp_sha256(destination))
}

s6f_parse_marker_number <- function(path, prefix) {
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  hit <- lines[startsWith(lines, prefix)]
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(sub(prefix, "", hit[[1L]], fixed = TRUE)))
}

s6f_validate_fit_dimensions <- function(fit, row) {
  dimensions <- pp_fit_dimensions(row)
  returned_L_s <- as.integer(fit$L_s)
  if (length(returned_L_s) == 1L) returned_L_s <- rep(returned_L_s, 2L)
  pp_assert(
    identical(as.integer(fit$L_f), dimensions$L_f) &&
      identical(returned_L_s, dimensions$L_s) &&
      identical(as.integer(fit$M_f), dimensions$M_f) &&
      identical(lapply(fit$M_s, as.integer), dimensions$M_s),
    paste0("Fitted dimensions differ from manifest: ", row$fit_id[[1L]])
  )
  invisible(TRUE)
}

s6f_write_endpoint_artifacts <- function(result, row, fit_dir,
                                          imported = FALSE,
                                          source_fit_id = "",
                                          source_fit_sha256 = "") {
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6F_JOINT_SELECTION_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$split_id <- as.character(row$split_id[[1L]])
  record$holdout_scoring_complete <- FALSE
  record$holdout_available_to_start_selection <- FALSE
  record$truth_used_for_fit_stopping_or_selection <- FALSE
  record$truth_used_for_holdout_scoring <- FALSE
  record$imported_from_stage6b <- imported
  record$source_stage6b_fit_id <- source_fit_id
  record$source_fit_sha256 <- source_fit_sha256
  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(),
               call = character(), phase = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  record
}

import_stage6b_fit <- function(row) {
  pp_assert(isTRUE(row$reuse_stage6b[[1L]]) &&
              !isTRUE(row$execute_new_fit[[1L]]),
            "Only registered Stage 6B endpoints may be imported.")
  fit_dir <- file.path(output_root, "fits", row$fit_id[[1L]])
  complete_path <- file.path(fit_dir, "FIT_COMPLETE.txt")
  if (dir.exists(fit_dir)) {
    pp_assert(file.exists(complete_path) &&
                file.exists(file.path(fit_dir, "fit.rds")) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Incomplete imported endpoint exists: ", row$fit_id[[1L]]))
    return(invisible(TRUE))
  }
  source_dir <- file.path(
    stage6b_output_root, "fits", row$source_stage6b_fit_id[[1L]]
  )
  source_fit_path <- file.path(source_dir, "fit.rds")
  pp_assert(file.exists(file.path(source_dir, "FIT_COMPLETE.txt")) &&
              file.exists(source_fit_path),
            paste0("Stage 6B source endpoint is incomplete: ",
                   row$source_stage6b_fit_id[[1L]]))
  fit <- readRDS(source_fit_path)
  s6f_validate_fit_dimensions(fit, row)
  pp_assert(identical(fit$initialization, "random") &&
              as.integer(fit$n_cpus_used) == 1L &&
              identical(fit$pre_score_interface$function_initialization,
                        "gram_unit_energy") &&
              as.integer(fit$pre_score_interface$pre_score_sweeps) == 1L &&
              as.integer(fit$t1_sweeps) == 400L,
            "Imported Stage 6B fit route differs from Stage 6F registration.")
  objective <- r_route_v3_candidate2_objective_status(fit)
  pp_assert(isTRUE(objective$eligible),
            "Imported Stage 6B endpoint is not objective eligible.")
  source_warnings <- file.path(source_dir, "WARNINGS.csv")
  warnings <- if (file.exists(source_warnings)) pp_read_csv(source_warnings) else
    data.frame(class = character(), message = character(),
               call = character(), phase = character(),
               stringsAsFactors = FALSE)
  warnings <- s6bp_set_warning_phase(warnings, "stage6f_imported_stage6b")
  captured <- list(
    started_at_utc = NA_character_, ended_at_utc = NA_character_,
    elapsed_seconds = s6f_parse_marker_number(
      file.path(source_dir, "FIT_COMPLETE.txt"), "elapsed_seconds="
    ),
    peak_memory_bytes = NA_real_, warnings = warnings, error = NULL
  )
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  link <- s6f_link_or_copy(source_fit_path, file.path(fit_dir, "fit.rds"))
  result <- list(fit = fit, captured = captured, objective = objective)
  record <- s6f_write_endpoint_artifacts(
    result, row, fit_dir, imported = TRUE,
    source_fit_id = row$source_stage6b_fit_id[[1L]],
    source_fit_sha256 = link$sha256
  )
  pp_write_csv(data.frame(
    t1_sweep = seq_along(fit$ELBO), ordinary_elbo = as.numeric(fit$ELBO)
  ), file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", row$fit_id[[1L]]),
    "endpoint_origin=reused_stage6b_v2",
    paste0("source_fit_id=", row$source_stage6b_fit_id[[1L]]),
    paste0("source_fit_sha256=", link$sha256),
    paste0("link_method=", link$method),
    paste0("objective_eligible=", record$objective_eligible),
    "ordinary_T1_sweeps=400", "n_cpus=1", "holdout_scored=FALSE",
    "truth_read=FALSE"
  ), complete_path)
  invisible(TRUE)
}

prepare_stage6f <- function() {
  require_roots(TRUE)
  registration <- verify_environment(TRUE)
  paths <- prepared_paths()
  source_observation_sha <- pp_sha256(observation_path)
  if (file.exists(paths$marker)) {
    pp_assert(all(file.exists(unlist(paths[c(
      "manifest", "middle_train", "middle_holdout", "upper_train",
      "upper_holdout"
    )]))) && any(readLines(paths$marker, warn = FALSE) ==
                    paste0("source_observation_sha256=",
                           source_observation_sha)),
              "Existing Stage 6F preparation is incomplete or mismatched.")
    cat("STAGE6F_ALREADY_PREPARED output_root=", output_root, "\n", sep = "")
    return(invisible(paths))
  }
  pp_assert(!dir.exists(file.path(output_root, "truth_free_selection")) &&
              !file.exists(file.path(output_root,
                                     "TRUTH_EVALUATION_AUTHORIZED.txt")),
            "Refusing to prepare over Stage 6F selection/evaluation artifacts.")
  observation <- readRDS(observation_path)
  pp_assert(identical(observation$bundle_class,
                      "v0lv_observation_only_candidate2") &&
              identical(observation$data_id, "g12s6a_01") &&
              !length(s6b_forbidden_observation_names(observation)),
            "Stage 6F source observation identity/isolation failed.")
  middle <- s6f_make_time_holdout(observation, "middle")
  upper <- s6f_make_time_holdout(observation, "upper")
  source_middle_train <- file.path(stage6b_output_root, "split",
                                   "train_observation.rds")
  source_middle_holdout <- file.path(stage6b_output_root, "split",
                                     "holdout_plan.rds")
  pp_assert(file.exists(source_middle_train) &&
              file.exists(source_middle_holdout) &&
              identical(middle$train_observation,
                        readRDS(source_middle_train)) &&
              identical(middle$holdout_plan,
                        readRDS(source_middle_holdout)),
            "Stage 6F middle split does not exactly reproduce Stage 6B v2.")
  middle_link_train <- s6f_link_or_copy(source_middle_train,
                                        paths$middle_train)
  middle_link_holdout <- s6f_link_or_copy(source_middle_holdout,
                                          paths$middle_holdout)
  pp_save_rds(upper$train_observation, paths$upper_train)
  pp_save_rds(upper$holdout_plan, paths$upper_holdout)
  pp_write_csv(registration$fits, paths$manifest)

  middle_records <- middle$holdout_plan$records
  upper_records <- upper$holdout_plan$records
  pp_assert(length(middle_records) == 60L && length(upper_records) == 60L &&
              all(vapply(middle_records, `[[`, integer(1L), "p") == 500L) &&
              all(vapply(upper_records, `[[`, integer(1L), "p") == 500L),
            "Stage 6F full-size split contract failed.")
  middle_index <- vapply(middle_records, `[[`, integer(1L), "original_index")
  upper_index <- vapply(upper_records, `[[`, integer(1L), "original_index")
  pp_assert(all(upper_index == middle_index + 1L),
            "Stage 6F paired holdout indices are not adjacent and distinct.")
  split_summary <- rbind(
    data.frame(
      split_id = "split1_middle",
      study = vapply(middle_records, `[[`, integer(1L), "study"),
      subject = vapply(middle_records, `[[`, integer(1L), "subject"),
      original_index = middle_index,
      holdout_time = vapply(middle_records, `[[`, numeric(1L), "time"),
      train_time_count = vapply(middle_records, `[[`, integer(1L),
                                "train_time_count"), truth_used = FALSE
    ),
    data.frame(
      split_id = "split2_upper",
      study = vapply(upper_records, `[[`, integer(1L), "study"),
      subject = vapply(upper_records, `[[`, integer(1L), "subject"),
      original_index = upper_index,
      holdout_time = vapply(upper_records, `[[`, numeric(1L), "time"),
      train_time_count = vapply(upper_records, `[[`, integer(1L),
                                "train_time_count"), truth_used = FALSE
    )
  )
  pp_write_csv(split_summary,
               file.path(output_root, "splits", "HOLDOUT_SPLIT_SUMMARY.csv"))
  for (index in which(registration$fits$reuse_stage6b)) {
    import_stage6b_fit(registration$fits[index, , drop = FALSE])
  }
  pp_write_lines(c(
    "status=STAGE6F_PREPARED", paste0("prepared_utc=", pp_iso_time()),
    paste0("source_observation_sha256=", source_observation_sha),
    paste0("middle_train_sha256=", middle_link_train$sha256),
    paste0("middle_holdout_sha256=", middle_link_holdout$sha256),
    paste0("upper_train_sha256=", pp_sha256(paths$upper_train)),
    paste0("upper_holdout_sha256=", pp_sha256(paths$upper_holdout)),
    "data_id=g12s6a_01", "data_generated=FALSE",
    "registered_endpoints=32", "reused_stage6b_endpoints=12",
    "new_fit_endpoints_authorized=20", "truth_read=FALSE",
    "holdout_used_for_start_selection=FALSE", "continuation_authorized=FALSE",
    "formal_paper_mc_result=FALSE"
  ), paths$marker)
  cat("STAGE6F_PREPARED endpoints=32 reused=12 new=20\n")
  invisible(paths)
}

load_prepared_split <- function(split_id) {
  pp_assert(split_id %in% c("split1_middle", "split2_upper"),
            "Unknown Stage 6F split_id.")
  root <- file.path(output_root, "splits", split_id)
  train <- readRDS(file.path(root, "train_observation.rds"))
  holdout <- readRDS(file.path(root, "holdout_plan.rds"))
  pp_assert(!length(s6b_forbidden_observation_names(train)) &&
              identical(holdout$truth_used, FALSE) &&
              holdout$subject_time_count == 60L,
            "Prepared Stage 6F split is invalid or truth-bearing.")
  list(train = train, holdout = holdout)
}

write_fit_error <- function(fit_dir, fit_id, message) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
    paste0("message=", message), paste0("ended_utc=", pp_iso_time()),
    "truth_read=FALSE"
  ), file.path(fit_dir, "FIT_ERROR.txt"))
}

run_new_fit <- function(fit_id) {
  require_roots(FALSE)
  registration <- verify_environment(TRUE)
  paths <- prepared_paths()
  pp_assert(file.exists(paths$marker) && file.exists(paths$manifest),
            "Stage 6F preparation is incomplete.")
  pp_assert(!file.exists(file.path(output_root,
                                   "TRUTH_EVALUATION_AUTHORIZED.txt")),
            "Refusing a Stage 6F fit after truth-evaluation authorization.")
  row <- registration$fits[registration$fits$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L && isTRUE(row$execute_new_fit[[1L]]) &&
              !isTRUE(row$reuse_stage6b[[1L]]),
            "fit_id is absent or is not one of the 20 authorized new fits.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")) &&
                file.exists(file.path(fit_dir, "fit.rds")) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Existing fit is incomplete/error; no selective rerun: ",
                     fit_id))
    cat("FIT_ALREADY_COMPLETE fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  split <- load_prepared_split(row$split_id[[1L]])
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- pp_fit_g12(split$train, row, smoke = FALSE)
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  result$captured$warnings <- s6bp_set_warning_phase(
    result$captured$warnings, "stage6f_joint_dimension_g12_fixed400"
  )
  record <- s6f_write_endpoint_artifacts(result, row, fit_dir)
  if (is.null(result$fit)) {
    write_fit_error(fit_dir, fit_id, record$error$message %||% "unknown")
    pp_abort("Stage 6F fit failed and was retained: ", fit_id)
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
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "endpoint_origin=new_stage6f_fit", "terminal_status=fixed_400_complete",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    "ordinary_T1_sweeps=400", "n_cpus=1", "holdout_scored=FALSE",
    "truth_read=FALSE", "continuation_used=FALSE",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

collect_endpoints <- function(manifest) {
  record_paths <- file.path(output_root, "fits", manifest$fit_id,
                            "terminal_record.rds")
  complete_paths <- file.path(output_root, "fits", manifest$fit_id,
                              "FIT_COMPLETE.txt")
  pp_assert(all(file.exists(record_paths)) && all(file.exists(complete_paths)),
            "Stage 6F has missing endpoint records or completion markers.")
  records <- lapply(record_paths, readRDS)
  pp_assert(all(vapply(records, function(record) {
    identical(record$terminal_schema,
              "G12_STAGE6F_JOINT_SELECTION_TERMINAL_V1") &&
      isTRUE(record$objective_eligible) &&
      identical(record$terminal_status, "fixed_400_complete") &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$holdout_scoring_complete)
  }, logical(1L))), "Stage 6F endpoint contract failed.")
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  rows$split_id <- manifest$split_id
  rows$complexity_proxy <- manifest$complexity_proxy
  rows$imported_from_stage6b <- manifest$reuse_stage6b
  pp_assert(identical(rows$fit_id, manifest$fit_id),
            "Stage 6F endpoint-to-manifest join failed.")
  list(records = records, rows = rows)
}

freeze_start_winners <- function() {
  require_roots(FALSE)
  registration <- verify_environment(TRUE)
  pp_assert(file.exists(prepared_paths()$marker),
            "Stage 6F preparation is incomplete.")
  selection_root <- file.path(output_root, "truth_free_selection")
  marker <- file.path(selection_root, "START_SELECTION_FROZEN.txt")
  winner_path <- file.path(selection_root, "CONFIG_SPLIT_WINNERS.csv")
  if (file.exists(marker)) {
    winners <- pp_read_csv(winner_path)
    pp_assert(nrow(winners) == 8L,
              "Existing Stage 6F winner table is malformed.")
    return(invisible(winners))
  }
  pp_assert(!dir.exists(selection_root),
            "Incomplete Stage 6F truth-free selection directory exists.")
  endpoints <- collect_endpoints(registration$fits)$rows
  winners <- pp_select_g12_winners(endpoints)
  pp_assert(nrow(winners) == 8L &&
              length(unique(winners$selection_stratum_id)) == 8L,
            "Stage 6F start selection did not yield eight winners.")
  dir.create(selection_root, recursive = TRUE, showWarnings = FALSE,
             mode = "0700")
  pp_write_csv(endpoints, file.path(selection_root, "ALL_ENDPOINTS.csv"))
  pp_write_csv(winners, winner_path)
  hash <- pp_sha256(winner_path)
  pp_write_lines(c(
    "status=START_SELECTION_FROZEN", paste0("completed_utc=", pp_iso_time()),
    "endpoints=32", "split_configuration_winners=8",
    paste0("winner_table_sha256=", hash),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_split_and_configuration",
    "heldout_response_read_for_start_selection=FALSE",
    "cross_dimension_elbo_comparison=FALSE", "selection_used_truth=FALSE"
  ), marker)
  cat("STAGE6F_START_SELECTION_FROZEN winners=8\n")
  invisible(winners)
}

score_winner <- function(winner) {
  fit_id <- winner$fit_id[[1L]]
  score_dir <- file.path(output_root, "truth_free_scores", fit_id)
  complete <- file.path(score_dir, "SCORE_COMPLETE.txt")
  if (dir.exists(score_dir)) {
    pp_assert(file.exists(complete),
              paste0("Incomplete Stage 6F score exists: ", fit_id))
    return(invisible(TRUE))
  }
  split <- load_prepared_split(winner$split_id[[1L]])
  fit <- readRDS(file.path(output_root, "fits", fit_id, "fit.rds"))
  dir.create(score_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  score <- tryCatch(list(
    rmse = s6b_score_holdout(fit, split$holdout, "all"),
    logscore = s6eb_score_holdout(fit, split$holdout)
  ), error = function(condition) condition)
  if (inherits(score, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(score)), "truth_read=FALSE"
    ), file.path(score_dir, "SCORE_ERROR.txt"))
    pp_abort("Stage 6F scoring failed and was retained: ", fit_id)
  }
  decorate <- function(value) {
    value$fit_id <- fit_id
    value$fit_config_id <- winner$fit_config_id[[1L]]
    value$split_id <- winner$split_id[[1L]]
    value
  }
  rmse_subject <- decorate(score$rmse$subject_time)
  log_subject <- decorate(score$logscore$subject_time)
  rmse_aggregate <- decorate(score$rmse$aggregate)
  log_aggregate <- decorate(score$logscore$aggregate)
  pp_write_csv(rmse_subject,
               file.path(score_dir, "HOLDOUT_SUBJECT_TIME_RMSE.csv"))
  pp_write_csv(rmse_aggregate,
               file.path(score_dir, "HOLDOUT_SUMMARY_RMSE.csv"))
  pp_write_csv(log_subject,
               file.path(score_dir, "HOLDOUT_SUBJECT_TIME_LOGSCORE.csv"))
  pp_write_csv(log_aggregate,
               file.path(score_dir, "HOLDOUT_SUMMARY_LOGSCORE.csv"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("fit_config_id=", winner$fit_config_id[[1L]]),
    paste0("split_id=", winner$split_id[[1L]]),
    "fit_selected_before_holdout_score=TRUE", "truth_read=FALSE",
    "rmse_interface=established_stage6b_all_factor",
    "joint_logscore=experimental_moment_matched_gaussian",
    paste0("completed_utc=", pp_iso_time())
  ), complete)
  invisible(TRUE)
}

score_winners <- function() {
  require_roots(FALSE)
  verify_environment(TRUE)
  selection_root <- file.path(output_root, "truth_free_selection")
  pp_assert(file.exists(file.path(selection_root,
                                  "START_SELECTION_FROZEN.txt")),
            "Start winners must be frozen before holdout scoring.")
  pp_assert(!file.exists(file.path(output_root,
                                   "TRUTH_EVALUATION_AUTHORIZED.txt")),
            "Refusing truth-free scoring after truth-evaluation authorization.")
  winners <- pp_read_csv(file.path(selection_root,
                                   "CONFIG_SPLIT_WINNERS.csv"))
  pp_assert(nrow(winners) == 8L,
            "Stage 6F requires eight start-selected winners.")
  for (index in seq_len(nrow(winners))) {
    score_winner(winners[index, , drop = FALSE])
  }
  marker <- file.path(output_root, "truth_free_scores",
                      "WINNER_SCORING_COMPLETE.txt")
  pp_write_lines(c(
    "status=WINNER_SCORING_COMPLETE", "winners_scored=8",
    "score_errors=0", "truth_read=FALSE",
    paste0("completed_utc=", pp_iso_time())
  ), marker)
  cat("STAGE6F_WINNER_SCORING_COMPLETE winners=8\n")
  invisible(TRUE)
}

read_score_rows <- function(winners, filename) {
  values <- lapply(seq_len(nrow(winners)), function(index) {
    path <- file.path(output_root, "truth_free_scores",
                      winners$fit_id[[index]], filename)
    pp_assert(file.exists(path), paste0("Missing Stage 6F score: ", path))
    pp_read_csv(path)
  })
  do.call(rbind, values)
}

freeze_dimensions <- function() {
  require_roots(FALSE)
  registration <- verify_environment(FALSE)
  selection_root <- file.path(output_root, "truth_free_selection")
  score_marker <- file.path(output_root, "truth_free_scores",
                            "WINNER_SCORING_COMPLETE.txt")
  final_marker <- file.path(selection_root,
                            "DIMENSION_SELECTION_FROZEN.txt")
  if (file.exists(final_marker)) {
    cat("STAGE6F_DIMENSION_SELECTION_ALREADY_FROZEN\n")
    return(invisible(TRUE))
  }
  pp_assert(file.exists(score_marker) &&
              file.exists(file.path(selection_root,
                                    "START_SELECTION_FROZEN.txt")),
            "Stage 6F winner scoring/start selection is incomplete.")
  winners <- pp_read_csv(file.path(selection_root,
                                   "CONFIG_SPLIT_WINNERS.csv"))
  rmse_rows <- read_score_rows(winners, "HOLDOUT_SUBJECT_TIME_RMSE.csv")
  log_rows <- read_score_rows(winners, "HOLDOUT_SUBJECT_TIME_LOGSCORE.csv")
  pp_assert(nrow(rmse_rows) == 480L && nrow(log_rows) == 480L,
            "Stage 6F must contain 8 x 60 winner score rows.")
  configs <- registration$configs
  mse_subject <- s6f_subject_cluster_losses(
    rmse_rows, c("mse", "mae", "observed_second_moment")
  )
  log_subject <- s6f_subject_cluster_losses(
    log_rows, c(
      "joint_nlpd_per_variable", "marginal_nlpd",
      "marginal_interval_coverage", "mean_standardized_squared_error",
      "mean_predictive_sd", "diagonal_floor_count",
      "mean_variance_clamp_count", "loading_variance_clamp_count",
      "process_numerical_clamp_count"
    )
  )
  mse_selection <- s6f_select_one_se(
    mse_subject, configs, "mse", "paired_holdout_mse"
  )
  observed_rms <- sqrt(mean(mse_subject$observed_second_moment))
  mse_selection$summary$heldout_rmse <-
    sqrt(mse_selection$summary$mean_loss)
  mse_selection$summary$heldout_nrmse <-
    mse_selection$summary$heldout_rmse / observed_rms
  log_selection <- s6f_select_one_se(
    log_subject, configs, "joint_nlpd_per_variable",
    "experimental_joint_nlpd"
  )
  mse_stability <- s6f_split_ranking_stability(
    rmse_rows, "mse", "paired_holdout_mse"
  )
  log_stability <- s6f_split_ranking_stability(
    log_rows, "joint_nlpd_per_variable", "experimental_joint_nlpd"
  )
  calibration <- do.call(rbind, lapply(configs$fit_config_id,
                                       function(config_id) {
    rows <- log_rows[log_rows$fit_config_id == config_id, , drop = FALSE]
    data.frame(
      fit_config_id = config_id, subject_time_vectors = nrow(rows),
      marginal_interval_coverage = mean(rows$marginal_interval_coverage),
      mean_standardized_squared_error =
        mean(rows$mean_standardized_squared_error),
      marginal_nlpd = mean(rows$marginal_nlpd),
      joint_nlpd_per_variable = mean(rows$joint_nlpd_per_variable),
      mean_predictive_sd = mean(rows$mean_predictive_sd),
      total_diagonal_floor_count = sum(rows$diagonal_floor_count),
      total_numerical_clamp_count = sum(
        rows$mean_variance_clamp_count +
          rows$loading_variance_clamp_count +
          rows$process_numerical_clamp_count
      ), stringsAsFactors = FALSE
    )
  }))
  selections <- data.frame(
    criterion_id = c("paired_holdout_mse", "experimental_joint_nlpd"),
    selected_config_id = c(mse_selection$selected_config_id,
                           log_selection$selected_config_id),
    minimum_loss_config_id = c(mse_selection$minimum_loss_config_id,
                               log_selection$minimum_loss_config_id),
    one_se_threshold = c(mse_selection$threshold,
                         log_selection$threshold),
    truth_used = FALSE,
    criterion_status = c("preregistered_primary_interface",
                         "experimental_feasibility_contender"),
    stringsAsFactors = FALSE
  )
  stability_summary <- rbind(mse_stability$summary,
                             log_stability$summary)
  split_risks <- rbind(mse_stability$risks, log_stability$risks)
  split_pairs <- rbind(mse_stability$pairwise, log_stability$pairwise)
  pp_write_csv(mse_subject,
               file.path(selection_root, "SUBJECT_CLUSTERED_MSE.csv"))
  pp_write_csv(log_subject,
               file.path(selection_root, "SUBJECT_CLUSTERED_LOGSCORE.csv"))
  pp_write_csv(mse_selection$summary,
               file.path(selection_root, "MSE_ONE_SE_CANDIDATES.csv"))
  pp_write_csv(log_selection$summary,
               file.path(selection_root, "LOGSCORE_ONE_SE_CANDIDATES.csv"))
  pp_write_csv(calibration,
               file.path(selection_root, "PREDICTIVE_CALIBRATION.csv"))
  pp_write_csv(split_risks,
               file.path(selection_root, "SPLIT_SPECIFIC_RISKS.csv"))
  pp_write_csv(stability_summary,
               file.path(selection_root, "CROSS_SPLIT_RANK_STABILITY.csv"))
  pp_write_csv(split_pairs,
               file.path(selection_root, "CROSS_SPLIT_PAIRWISE_ORDERING.csv"))
  final_path <- file.path(selection_root,
                          "FINAL_TRUTH_FREE_DIMENSION_SELECTION.csv")
  pp_write_csv(selections, final_path)
  final_hash <- pp_sha256(final_path)
  winners_hash <- pp_sha256(file.path(selection_root,
                                      "CONFIG_SPLIT_WINNERS.csv"))
  pp_write_lines(c(
    "status=DIMENSION_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("start_winner_table_sha256=", winners_hash),
    paste0("dimension_selection_sha256=", final_hash),
    paste0("mse_one_se_selected=", mse_selection$selected_config_id),
    paste0("joint_nlpd_one_se_selected=", log_selection$selected_config_id),
    "start_selection_used_ordinary_elbo_within_configuration_only=TRUE",
    "dimension_selection_used_truth=FALSE", "truth_read=FALSE",
    "joint_logscore_exact_posterior_predictive_density=FALSE",
    "new_fits=20", "continuations=0", "formal_paper_mc_result=FALSE"
  ), final_marker)
  pp_write_lines(c(
    "status=TRUTH_FREE_STAGE_COMPLETE", "registered_endpoints=32",
    "reused_endpoints=12", "new_fits=20", "split_winners=8",
    "truth_read=FALSE", "continuations=0",
    paste0("completed_utc=", pp_iso_time())
  ), file.path(output_root, "STAGE6F_TRUTH_FREE_COMPLETE.txt"))
  cat("STAGE6F_DIMENSION_SELECTION_FROZEN mse=",
      mse_selection$selected_config_id, " joint_nlpd=",
      log_selection$selected_config_id, "\n", sep = "")
  invisible(TRUE)
}

status_stage6f <- function() {
  registration <- load_registration()
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) {
    length(list.files(fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE,
                      full.names = TRUE))
  } else 0L
  errors <- if (dir.exists(fit_root)) {
    length(list.files(fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE,
                      full.names = TRUE))
  } else 0L
  cat("registered_endpoints=", nrow(registration$fits), "\n", sep = "")
  cat("reused_registered=", sum(registration$fits$reuse_stage6b), "\n",
      sep = "")
  cat("new_registered=", sum(registration$fits$execute_new_fit), "\n",
      sep = "")
  cat("fit_complete=", complete, "\n", sep = "")
  cat("fit_errors=", errors, "\n", sep = "")
  cat("start_selection_frozen=", file.exists(file.path(
    output_root, "truth_free_selection", "START_SELECTION_FROZEN.txt"
  )), "\n", sep = "")
  cat("dimension_selection_frozen=", file.exists(file.path(
    output_root, "truth_free_selection", "DIMENSION_SELECTION_FROZEN.txt"
  )), "\n", sep = "")
  invisible(TRUE)
}

if (identical(action, "check")) {
  verify_environment(TRUE)
  cat("STAGE6F_CHECK_PASS registered=32 reused=12 new=20\n")
} else if (identical(action, "prepare")) {
  prepare_stage6f()
} else if (identical(action, "run-fit")) {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required for run-fit.")
  run_new_fit(fit_id)
} else if (identical(action, "freeze-starts")) {
  freeze_start_winners()
} else if (identical(action, "score-winners")) {
  score_winners()
} else if (identical(action, "freeze-dimensions")) {
  freeze_dimensions()
} else if (identical(action, "status")) {
  pp_assert(nzchar(output_root), "--output-root is required for status.")
  status_stage6f()
} else {
  pp_abort("Unknown --action: ", action)
}
