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
stage6a_root <- file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831"
)
stage6g_root <- file.path(
  experiments_root, "g12-dimension-criterion-falsification-stage6g-20260911"
)
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6H_REACHABILITY_CAUSAL_AUDIT_V1_20260913"
stage6h_library <- pp_env("STAGE6H_R_LIB")
if (nzchar(stage6h_library)) Sys.setenv(STAGE6A_R_LIB = stage6h_library)
STAGE6H_WORKERS <- as.integer(pp_env("STAGE6H_WORKERS", "7"))
pp_assert(length(STAGE6H_WORKERS) == 1L && is.finite(STAGE6H_WORKERS) &&
            STAGE6H_WORKERS %in% 1:7,
          "STAGE6H_WORKERS must be an integer from 1 to 7.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6H_OUTPUT_ROOT")
stage6g_output_root <- pp_env(
  "STAGE6H_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6g_evaluation_root <- pp_env(
  "STAGE6H_STAGE6G_EVALUATION_ROOT",
  file.path(stage6g_output_root, "truth_evaluation_20260912_v1")
)
start_path <- file.path(script_root, "START_SEEDS_STAGE6H_V1.csv")
authorization_path <- file.path(
  script_root, "STAGE6H_EXECUTION_AUTHORIZATION_20260913.txt"
)
package_binding_path <- file.path(
  stage6g_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
)

s6h_assert_authorized <- function() {
  pp_assert(file.exists(authorization_path),
            "Stage 6H execution authorization is missing.")
  lines <- readLines(authorization_path, warn = FALSE)
  required <- c(
    "status=AUTHORIZED_POST_UNSEAL_REACHABILITY_CAUSAL_AUDIT",
    "new_stage6h_fits_authorized=TRUE", "new_fit_count=28",
    "existing_stage6g_fits_reused_read_only=8",
    "scientific_evaluation_after_endpoint_freeze_authorized=TRUE",
    "truth_available_to_fit=FALSE", "truth_available_to_stopping=FALSE",
    "truth_available_to_elbo_selection=FALSE",
    "continuation_authorized=FALSE",
    "cross_dimension_selection_authorized=FALSE",
    "formal_paper_mc=FALSE", "stage6g_outputs_mutable=FALSE"
  )
  pp_assert(all(required %in% lines),
            "Stage 6H authorization has the wrong scope.")
  invisible(TRUE)
}

s6h_read_starts <- function() {
  starts <- pp_read_csv(start_path)
  pp_assert(
    identical(names(starts),
              c("seed_index", "fit_seed", "paired_across_conditions")) &&
      nrow(starts) == 12L &&
      identical(as.integer(starts$seed_index), 1:12) &&
      identical(as.integer(starts$fit_seed), 84611101:84611112) &&
      all(starts$paired_across_conditions) &&
      !anyDuplicated(starts$fit_seed),
    "Stage 6H start registration is malformed."
  )
  starts
}

s6h_build_manifests <- function() {
  starts <- s6h_read_starts()
  conditions <- data.frame(
    condition_index = 1:3,
    condition_id = c(
      "holdout_split1_middle", "holdout_split2_upper", "full_data_6_9"
    ),
    split_id = c("split1_middle", "split2_upper", "full_data"),
    input_kind = c("holdout_train", "holdout_train", "full_observation"),
    training_time_count_min = c(5L, 5L, 6L),
    training_time_count_max = c(8L, 8L, 9L),
    stringsAsFactors = FALSE
  )
  rows <- list()
  index <- 0L
  for (condition_row in split(conditions, seq_len(nrow(conditions)))) {
    for (start_row in split(starts, seq_len(nrow(starts)))) {
      index <- index + 1L
      imported <- condition_row$input_kind == "holdout_train" &&
        start_row$seed_index <= 4L
      source_fit_id <- if (imported) paste(
        "s6g_minimal_L111_M2", condition_row$split_id,
        "cand_L111_M2", sprintf("G12_%02d", start_row$seed_index),
        sep = "__"
      ) else ""
      fit_id <- if (imported) source_fit_id else paste(
        "s6h_minimal_L111_M2", condition_row$condition_id,
        sprintf("G12_%02d", start_row$seed_index), sep = "__"
      )
      rows[[index]] <- data.frame(
        experiment_id = PP_EXPERIMENT_ID,
        fit_id = fit_id,
        data_id = "s6g_minimal_L111_M2",
        scenario_id = "s6g_minimal_L111_M2",
        condition_id = condition_row$condition_id,
        condition_index = condition_row$condition_index,
        split_id = condition_row$split_id,
        input_kind = condition_row$input_kind,
        training_time_count_min = condition_row$training_time_count_min,
        training_time_count_max = condition_row$training_time_count_max,
        method_id = "G12", fit_config_id = "cand_L111_M2",
        selection_stratum_id = paste0(
          "s6g_minimal_L111_M2__", condition_row$condition_id,
          "__cand_L111_M2"
        ),
        route_id = "random__pre1__jaoua_anneal__G12__fixed400",
        seed_index = start_row$seed_index,
        fit_seed = start_row$fit_seed,
        fit_L_f = 1L, fit_L_s_1 = 1L, fit_L_s_2 = 1L,
        fit_M_f = "2", fit_M_s_1 = "2", fit_M_s_2 = "2",
        fit_M_common = 2L, complexity_proxy = 1782L,
        n_cpus = 1L, planned_annealing_sweeps = 99L,
        planned_ordinary_t1_sweeps = 400L,
        fit_source = if (imported) "stage6g_read_only" else "stage6h_new",
        source_fit_id = source_fit_id,
        truth_available_to_fit = FALSE,
        truth_available_to_stopping = FALSE,
        truth_available_to_elbo_selection = FALSE,
        continuation = FALSE, actual_racing = FALSE,
        engineering_diagnostic = TRUE, formal_paper_mc_result = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  all <- do.call(rbind, rows)
  new <- all[all$fit_source == "stage6h_new", , drop = FALSE]
  list(starts = starts, conditions = conditions, all = all, new = new)
}

s6h_validate_manifests <- function(registration) {
  all <- registration$all
  new <- registration$new
  forbidden <- c(
    "truth_basis_seed", "truth_L_f", "truth_L_s_1", "truth_L_s_2",
    "truth_M_common", "expected_truth_fit_config_id"
  )
  pp_assert(
    nrow(all) == 36L && nrow(new) == 28L &&
      sum(all$fit_source == "stage6g_read_only") == 8L &&
      !anyDuplicated(all$fit_id) && !anyDuplicated(new$fit_id) &&
      all(table(all$condition_id) == 12L) &&
      all(table(new$condition_id) == c(8L, 12L, 8L)[
        match(names(table(new$condition_id)),
              c("holdout_split1_middle", "full_data_6_9",
                "holdout_split2_upper"))
      ]) &&
      all(all$fit_L_f == 1L & all$fit_L_s_1 == 1L &
            all$fit_L_s_2 == 1L & all$fit_M_common == 2L) &&
      all(all$n_cpus == 1L) &&
      !any(all$truth_available_to_fit) &&
      !any(all$truth_available_to_stopping) &&
      !any(all$truth_available_to_elbo_selection) &&
      !any(all$continuation) && !any(all$actual_racing) &&
      !any(all$formal_paper_mc_result) &&
      !length(intersect(names(all), forbidden)),
    "Stage 6H endpoint registration failed."
  )
  for (condition_id in unique(all$condition_id)) {
    rows <- all[all$condition_id == condition_id, , drop = FALSE]
    pp_assert(identical(as.integer(rows$fit_seed), 84611101:84611112),
              paste0("Paired starts failed for ", condition_id, "."))
  }
  invisible(TRUE)
}

s6h_require_output_root <- function(create = FALSE) {
  pp_assert(nzchar(output_root), "--output-root is required.")
  if (create && !dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  output_root
}

s6h_runtime_paths <- function() list(
  all_manifest = file.path(output_root, "ANALYSIS_ENDPOINT_MANIFEST.csv"),
  new_manifest = file.path(output_root, "NEW_FIT_MANIFEST.csv"),
  starts = file.path(output_root, "START_SEEDS.csv"),
  conditions = file.path(output_root, "CONDITIONS.csv"),
  source_binding = file.path(output_root, "SOURCE_ARTIFACT_BINDING.csv"),
  prepared = file.path(output_root, "STAGE6H_PREPARED.txt")
)

s6h_load_runtime_registration <- function() {
  paths <- s6h_runtime_paths()
  pp_assert(all(file.exists(unlist(paths[c(
    "all_manifest", "new_manifest", "starts", "conditions", "prepared"
  )], use.names = FALSE))), "Stage 6H prepared manifests are incomplete.")
  value <- list(
    all = pp_read_csv(paths$all_manifest),
    new = pp_read_csv(paths$new_manifest),
    starts = pp_read_csv(paths$starts),
    conditions = pp_read_csv(paths$conditions)
  )
  s6h_validate_manifests(value)
  value
}

s6h_stage6g_paths <- function() list(
  observation = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2",
    "observation_bundle.rds"
  ),
  split1 = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2", "splits",
    "split1_middle", "train_observation.rds"
  ),
  split2 = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2", "splits",
    "split2_upper", "train_observation.rds"
  ),
  truth = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2",
    "sealed_truth", "truth_bundle.rds"
  ),
  manifest = file.path(stage6g_output_root, "FIT_MANIFEST.csv"),
  complete = file.path(stage6g_output_root,
                       "STAGE6G_TRUTH_FREE_COMPLETE.txt"),
  evaluation = file.path(stage6g_evaluation_root,
                         "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv"),
  evaluation_complete = file.path(stage6g_evaluation_root,
                                  "EVALUATION_COMPLETE.txt")
)

s6h_forbidden_observation_names <- function(value) {
  found <- character()
  walk <- function(object) {
    object_names <- names(object)
    if (!is.null(object_names)) {
      found <<- c(found, object_names[
        grepl("(^|_)(truth|true)(_|$)|signal_true|noise_true",
              object_names, ignore.case = TRUE)
      ])
    }
    if (is.list(object)) for (item in object) walk(item)
  }
  walk(value)
  unique(found)
}

s6h_verify_package_binding <- function() {
  pp_assert(file.exists(package_binding_path),
            "Stage 6G package-source binding is missing.")
  source_root <- pp_env("STAGE6H_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
  pp_assert(dir.exists(source_root),
            "STAGE6H_PACKAGE_SOURCE_ROOT does not exist.")
  binding <- pp_read_csv(package_binding_path)
  pp_assert(nrow(binding) == 29L &&
              identical(names(binding), c("relative_path", "sha256", "bytes")) &&
              !anyDuplicated(binding$relative_path),
            "Stage 6G package-source binding is malformed.")
  paths <- file.path(source_root, binding$relative_path)
  info <- file.info(paths)
  matches <- file.exists(paths) & is.finite(info$size) &
    as.numeric(info$size) == as.numeric(binding$bytes) &
    vapply(paths, pp_sha256, character(1L)) == binding$sha256
  pp_assert(all(matches), paste0(
    "Bound package source mismatch: ",
    paste(binding$relative_path[!matches], collapse = ";")
  ))
  invisible(TRUE)
}

s6h_verify_stage6g_source <- function(registration) {
  paths <- s6h_stage6g_paths()
  pp_assert(all(file.exists(unlist(paths, use.names = FALSE))),
            "Required read-only Stage 6G artifact is missing.")
  source_manifest <- pp_read_csv(paths$manifest)
  imported <- registration$all[
    registration$all$fit_source == "stage6g_read_only", , drop = FALSE
  ]
  source_rows <- source_manifest[
    source_manifest$fit_id %in% imported$source_fit_id, , drop = FALSE
  ]
  pp_assert(
    nrow(source_manifest) == 120L && nrow(source_rows) == 8L &&
      setequal(source_rows$fit_id, imported$source_fit_id) &&
      all(source_rows$scenario_id == "s6g_minimal_L111_M2") &&
      all(source_rows$fit_config_id == "cand_L111_M2") &&
      setequal(as.integer(source_rows$seed_index), 1:4),
    "Stage 6G imported endpoint identity failed."
  )
  fit_dirs <- file.path(stage6g_output_root, "fits", imported$source_fit_id)
  pp_assert(
    all(file.exists(file.path(fit_dirs, "fit.rds"))) &&
      all(file.exists(file.path(fit_dirs, "terminal_record.rds"))) &&
      all(file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))),
    "A read-only Stage 6G endpoint is incomplete."
  )
  invisible(paths)
}

s6h_verify_environment <- function(registration) {
  s6h_assert_authorized()
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  pp_load_package()
  pp_source_runtime(FALSE)
  s6h_verify_package_binding()
  s6h_verify_stage6g_source(registration)
  invisible(TRUE)
}

s6h_prepare <- function() {
  s6h_assert_authorized()
  registration <- s6h_build_manifests()
  s6h_validate_manifests(registration)
  pp_assert(!dir.exists(output_root),
            "Stage 6H output already exists; refusing overwrite.")
  s6h_verify_environment(registration)
  paths <- s6h_stage6g_paths()
  s6h_require_output_root(TRUE)
  runtime <- s6h_runtime_paths()
  pp_write_csv(registration$all, runtime$all_manifest)
  pp_write_csv(registration$new, runtime$new_manifest)
  pp_write_csv(registration$starts, runtime$starts)
  pp_write_csv(registration$conditions, runtime$conditions)
  bound_names <- c("observation", "split1", "split2", "truth",
                   "manifest", "complete", "evaluation",
                   "evaluation_complete")
  bound_paths <- unlist(paths[bound_names], use.names = FALSE)
  pp_write_csv(data.frame(
    artifact_role = bound_names,
    path = bound_paths,
    bytes = as.numeric(file.info(bound_paths)$size),
    sha256 = vapply(bound_paths, pp_sha256, character(1L)),
    mutable_by_stage6h = FALSE,
    available_to_fit = bound_names %in% c("observation", "split1", "split2"),
    stringsAsFactors = FALSE
  ), runtime$source_binding)
  pp_write_lines(c(
    "status=PREPARED", paste0("prepared_utc=", pp_iso_time()),
    "analysis_endpoints=36", "new_fits=28", "stage6g_read_only_fits=8",
    "conditions=3", "starts_per_condition=12", "n_cpus_per_fit=1",
    paste0("outer_workers=", STAGE6H_WORKERS),
    "truth_read_during_prepare=FALSE",
    "truth_available_to_fit_stopping_or_elbo_selection=FALSE",
    "post_unseal_oracle_configuration_diagnostic=TRUE",
    "continuation=FALSE", "cross_dimension_selection=FALSE",
    "formal_paper_mc=FALSE"
  ), runtime$prepared)
  cat("STAGE6H_PREPARED analysis=36 new_fits=28 imported=8\n")
  invisible(TRUE)
}

s6h_load_fit_input <- function(row) {
  paths <- s6h_stage6g_paths()
  path <- if (row$condition_id[[1L]] == "full_data_6_9") {
    paths$observation
  } else if (row$condition_id[[1L]] == "holdout_split1_middle") {
    paths$split1
  } else if (row$condition_id[[1L]] == "holdout_split2_upper") {
    paths$split2
  } else pp_abort("Unknown Stage 6H condition.")
  observation <- readRDS(path)
  pp_assert(all(c("Y", "time_obs", "Z", "dimensions") %in%
                  names(observation)) &&
              !length(s6h_forbidden_observation_names(observation)),
            "Stage 6H fit input is incomplete or truth-bearing.")
  observation
}

s6h_set_warning_phase <- function(warnings, phase) {
  if (!is.data.frame(warnings) || !nrow(warnings)) {
    return(data.frame(
      class = character(), message = character(), call = character(),
      phase = character(), stringsAsFactors = FALSE
    ))
  }
  warnings$phase <- phase
  warnings
}

s6h_write_fit_error <- function(fit_dir, fit_id, message) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
    paste0("message=", message), paste0("ended_utc=", pp_iso_time()),
    "truth_read=FALSE", "continuation=FALSE"
  ), file.path(fit_dir, "FIT_ERROR.txt"))
}

s6h_run_fit <- function(fit_id) {
  s6h_require_output_root(FALSE)
  registration <- s6h_load_runtime_registration()
  row <- registration$new[registration$new$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Stage 6H fit_id is not a registered new fit.")
  s6h_verify_environment(registration)
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    complete <- file.exists(file.path(fit_dir, "FIT_COMPLETE.txt"))
    error <- file.exists(file.path(fit_dir, "FIT_ERROR.txt"))
    pp_assert(xor(complete, error) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Incomplete existing Stage 6H fit directory: ", fit_id))
    cat(if (complete) "FIT_ALREADY_COMPLETE fit_id=" else
          "FIT_ERROR_RETAINED fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  observation <- s6h_load_fit_input(row)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- pp_fit_g12(observation, row, smoke = FALSE)
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  result$captured$warnings <- s6h_set_warning_phase(
    result$captured$warnings, "stage6h_g12_fixed400"
  )
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6H_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$condition_id <- as.character(row$condition_id[[1L]])
  record$input_kind <- as.character(row$input_kind[[1L]])
  record$fit_source <- "stage6h_new"
  record$truth_used_for_fit_stopping_or_selection <- FALSE
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  pp_write_csv(result$captured$warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    s6h_write_fit_error(
      fit_dir, fit_id, record$error$message %||% "unknown endpoint error"
    )
    cat("FIT_ERROR_RETAINED fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
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
    paste0("condition_id=", row$condition_id[[1L]]),
    "terminal_status=fixed_400_complete",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    "ordinary_T1_sweeps=400", "n_cpus=1", "truth_read=FALSE",
    "continuation=FALSE", paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

s6h_record_path <- function(row) {
  if (row$fit_source[[1L]] == "stage6g_read_only") {
    file.path(stage6g_output_root, "fits", row$source_fit_id[[1L]],
              "terminal_record.rds")
  } else file.path(output_root, "fits", row$fit_id[[1L]],
                   "terminal_record.rds")
}

s6h_fit_path <- function(row) {
  if (row$fit_source[[1L]] == "stage6g_read_only") {
    file.path(stage6g_output_root, "fits", row$source_fit_id[[1L]],
              "fit.rds")
  } else file.path(output_root, "fits", row$fit_id[[1L]], "fit.rds")
}

s6h_collect_all <- function(registration) {
  records <- lapply(seq_len(nrow(registration$all)), function(index) {
    row <- registration$all[index, , drop = FALSE]
    record_path <- s6h_record_path(row)
    pp_assert(file.exists(record_path),
              paste0("Missing terminal record: ", row$fit_id[[1L]]))
    record <- readRDS(record_path)
    pp_assert(
      identical(record$fit_id, row$fit_id[[1L]]) &&
        isTRUE(record$objective_eligible) &&
        !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
        !isTRUE(record$continuation_used) &&
        identical(as.integer(record$n_cpus), 1L),
      paste0("Ineligible endpoint: ", row$fit_id[[1L]])
    )
    record
  })
  rows <- pp_rbind_fill(lapply(records, pp_terminal_row))
  rows <- rows[match(registration$all$fit_id, rows$fit_id), , drop = FALSE]
  rows$condition_id <- registration$all$condition_id
  rows$split_id <- registration$all$split_id
  rows$input_kind <- registration$all$input_kind
  rows$fit_source <- registration$all$fit_source
  rows$source_fit_id <- registration$all$source_fit_id
  rows$training_time_count_min <-
    registration$all$training_time_count_min
  rows$training_time_count_max <-
    registration$all$training_time_count_max
  pp_assert(identical(as.character(rows$fit_id),
                      as.character(registration$all$fit_id)),
            "Stage 6H terminal-to-manifest join failed.")
  list(records = records, rows = rows)
}

s6h_freeze_endpoints <- function(registration) {
  frozen <- file.path(output_root,
                      "STAGE6H_ENDPOINTS_AND_ELBO_SELECTION_FROZEN.txt")
  pp_assert(!file.exists(frozen),
            "Stage 6H endpoint freeze already exists; refusing overwrite.")
  collected <- s6h_collect_all(registration)
  rows <- collected$rows
  winners <- do.call(rbind, lapply(split(rows, rows$condition_id), function(x) {
    eligible <- x$objective_eligible & is.finite(x$final_elbo) &
      x$terminal_status == "fixed_400_complete"
    pp_assert(any(eligible), paste0(
      "No eligible endpoint for condition ", x$condition_id[[1L]], "."
    ))
    candidates <- x[eligible, , drop = FALSE]
    candidates <- candidates[order(-candidates$final_elbo,
                                   candidates$fit_id), , drop = FALSE]
    winner <- candidates[1L, , drop = FALSE]
    winner$selection_rule <-
      "maximum_eligible_ordinary_T1_ELBO_within_fixed_dimension_condition"
    winner
  }))
  rownames(winners) <- NULL
  pp_assert(nrow(rows) == 36L && nrow(winners) == 3L &&
              !any(winners$truth_used_for_fit_stopping_or_selection),
            "Stage 6H endpoint freeze accounting failed.")
  pp_write_csv(rows, file.path(output_root, "ALL_ENDPOINT_TERMINALS.csv"))
  pp_write_csv(winners, file.path(output_root,
                                  "CONDITION_ELBO_WINNERS.csv"))
  pp_write_lines(c(
    "status=ENDPOINTS_AND_ELBO_SELECTION_FROZEN",
    paste0("frozen_utc=", pp_iso_time()), "analysis_endpoints=36",
    "new_stage6h_endpoints=28", "read_only_stage6g_endpoints=8",
    "condition_elbo_winners=3", "truth_read=FALSE",
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "cross_dimension_selection=FALSE", "continuation=FALSE"
  ), frozen)
  invisible(winners)
}

s6h_write_status <- function(state, detail = "") {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_new_fits = 28L, complete_new_fits = complete,
    error_new_fits = errors, analysis_endpoints = 36L,
    outer_workers = STAGE6H_WORKERS, per_fit_n_cpus = 1L,
    truth_read_by_fit_or_elbo_selection = FALSE,
    continuation = FALSE, formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

s6h_run_supervisor <- function() {
  s6h_require_output_root(FALSE)
  registration <- s6h_load_runtime_registration()
  s6h_verify_environment(registration)
  pp_assert(.Platform$OS.type == "unix",
            "The Stage 6H supervisor requires Unix.")
  log_dir <- file.path(output_root, "supervisor_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  s6h_write_status("RUNNING", "28 new fits; no selective rerun")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    started <- Sys.time()
    status <- tryCatch(system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=run-fit",
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
    registration$new$fit_id, run_one, mc.cores = STAGE6H_WORKERS,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    s6h_write_status("INCOMPLETE_PROCESS_FAILURE", paste(failed, collapse = ";"))
    pp_write_lines(c(
      "status=INCOMPLETE_PROCESS_FAILURE",
      paste0("failed_fit_processes=", paste(failed, collapse = ";")),
      "selective_rerun=FALSE", paste0("ended_utc=", pp_iso_time())
    ), file.path(output_root, "STAGE6H_INCOMPLETE.txt"))
    pp_abort("Stage 6H child-process failures retained: ",
             paste(failed, collapse = ", "))
  }
  s6h_freeze_endpoints(registration)
  s6h_write_status("FITS_COMPLETE_ENDPOINTS_FROZEN",
                   "28 new fits; 36 endpoints; 3 ELBO winners")
  pp_write_lines(c(
    "status=FITS_COMPLETE_ENDPOINTS_FROZEN",
    paste0("completed_utc=", pp_iso_time()), "new_fits=28",
    "analysis_endpoints=36", "condition_elbo_winners=3",
    "truth_read=FALSE", "scientific_evaluation_not_yet_started=TRUE",
    "continuation=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6H_FITS_COMPLETE.txt"))
  cat("STAGE6H_FITS_COMPLETE new=28 analysis=36 winners=3\n")
  invisible(TRUE)
}

s6h_status <- function() {
  s6h_require_output_root(FALSE)
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  frozen <- file.exists(file.path(
    output_root, "STAGE6H_ENDPOINTS_AND_ELBO_SELECTION_FROZEN.txt"
  ))
  evaluated <- file.exists(file.path(
    output_root, "truth_evaluation_20260913_v1", "EVALUATION_COMPLETE.txt"
  ))
  cat("STAGE6H_STATUS complete_new=", complete, "/28 errors=", errors,
      " endpoints_frozen=", frozen, " evaluated=", evaluated, "\n", sep = "")
  invisible(TRUE)
}

registration <- s6h_build_manifests()
s6h_validate_manifests(registration)

if (action == "check") {
  s6h_assert_authorized()
  cat("STAGE6H_CHECK_PASS analysis=36 new=28 imported=8 conditions=3 starts=12\n")
} else if (action == "preflight") {
  s6h_verify_environment(registration)
  cat("STAGE6H_PREFLIGHT_PASS package_binding=29 stage6g_imports=8 inputs=3\n")
} else if (action == "prepare") {
  s6h_prepare()
} else if (action == "run-fit") {
  fit_id <- cli$`fit-id` %||% pp_abort("--fit-id is required.")
  s6h_run_fit(fit_id)
} else if (action == "run") {
  s6h_run_supervisor()
} else if (action == "status") {
  s6h_status()
} else {
  pp_abort("Unknown Stage 6H action: ", action)
}
