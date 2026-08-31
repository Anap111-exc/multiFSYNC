#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(script_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
smoke <- pp_bool(cli$smoke, FALSE)
output_root <- cli$`output-root` %||% ""

require_output_root <- function() {
  pp_assert(nzchar(output_root), "--output-root is required for this action.")
  if (!dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  output_root
}

verify_environment <- function(include_evaluation = FALSE) {
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  pp_verify_source_binding()
  pp_load_package()
  pp_source_runtime(include_evaluation)
  manifests <- pp_load_manifests()
  qc <- pp_read_csv(file.path(PP_ROOT, "MANIFEST_QC.csv"))
  pp_assert(nrow(qc) == 14L && all(qc$passed), "Manifest QC is not 14/14.")
  manifests
}

generate_all <- function(smoke = FALSE) {
  require_output_root()
  manifests <- verify_environment(FALSE)
  pp_assert(!file.exists(file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt")),
            "Refusing generation after fit/selection completion.")
  statuses <- lapply(seq_len(nrow(manifests$data)), function(index) {
    row <- manifests$data[index, , drop = FALSE]
    created <- pp_generate_one(row, output_root, smoke)
    data.frame(
      data_id = row$data_id[[1L]], data_seed = row$data_seed[[1L]],
      scenario_id = row$scenario_id[[1L]],
      observation_created = TRUE, sealed_truth_created = TRUE,
      newly_created = isTRUE(created), truth_unseal_authorized = FALSE,
      truth_available_to_fit_stopping_or_selection = FALSE,
      smoke = smoke, checked_at_utc = pp_iso_time(),
      stringsAsFactors = FALSE
    )
  })
  status <- do.call(rbind, statuses)
  pp_write_csv(status, file.path(output_root, "DATA_GENERATION_STATUS.csv"))
  pp_write_lines(c(
    "status=DATA_GENERATION_AND_SEAL_COMPLETE",
    paste0("completed_utc=", pp_iso_time()),
    "data_sets=1", "observation_bundles=1", "sealed_truth_bundles=1",
    "truth_unseal_authorized=FALSE",
    "truth_available_to_fit_stopping_or_selection=FALSE",
    paste0("smoke=", toupper(as.character(smoke))),
    "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("DATA_GENERATION_AND_SEAL_COMPLETE data=1 smoke=", smoke, "\n", sep = "")
  invisible(status)
}

run_fit <- function(fit_id, smoke = FALSE) {
  require_output_root()
  manifests <- verify_environment(FALSE)
  pp_assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
            "Data generation/sealing is incomplete.")
  row <- manifests$fits[manifests$fits$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "fit_id must join to exactly one manifest row.")
  observation <- pp_load_observation(row$data_id[[1L]], output_root)
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Incomplete/error pre-existing fit directory: ", fit_id))
    cat("FIT_ALREADY_COMPLETE fit_id=", fit_id, "\n", sep = "")
    return(invisible(0L))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- if (row$method_id[[1L]] == "G12") {
    pp_fit_g12(observation, row, smoke)
  } else {
    pp_fit_pooled(observation, row, smoke)
  }
  record <- pp_make_terminal_record(result, row, smoke)
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
    if (row$method_id[[1L]] == "G12") {
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
      if (length(result$fit$practical_checkpoints)) {
        pp_save_rds(result$fit$practical_checkpoints,
                    file.path(fit_dir, "PRACTICAL_CHECKPOINTS.rds"))
      }
    } else {
      pp_write_csv(data.frame(
        total_iteration = seq_along(result$fit$ELBO_iter),
        elbo = as.numeric(result$fit$ELBO_iter)
      ), file.path(fit_dir, "REFERENCE_ELBO.csv"))
    }
  }
  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(), call = character(),
               phase = character(), stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (record$terminal_status == "error") {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN",
      paste0("fit_id=", fit_id),
      paste0("message=", record$error$message %||% "unknown"),
      paste0("ended_utc=", pp_iso_time())
    ), file.path(fit_dir, "FIT_ERROR.txt"))
    pp_abort("Fit failed and was retained: ", fit_id, ": ",
             record$error$message %||% "unknown")
  }
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("method_id=", record$method_id),
    paste0("terminal_status=", record$terminal_status),
    paste0("objective_eligible=", record$objective_eligible),
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "n_cpus=1", "continuation_used=FALSE", "automatic_800_used=FALSE",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " method=", record$method_id,
      " elapsed_seconds=", format(record$elapsed_seconds, digits = 8), "\n",
      sep = "")
  invisible(0L)
}

collect_terminals <- function(manifest) {
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  pp_assert(all(file.exists(paths)), "One or more terminal records are missing.")
  records <- lapply(paths, readRDS)
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  pp_assert(nrow(rows) == nrow(manifest) && !anyDuplicated(rows$fit_id) &&
              identical(rows$fit_id, manifest$fit_id),
            "Terminal records do not exactly join to the fit manifest.")
  list(records = records, rows = rows)
}

stage6a_probe_manifest <- function(fits) {
  result <- fits[fits$capacity_probe, , drop = FALSE]
  pp_assert(nrow(result) == 3L &&
              length(unique(result$fit_config_id)) == 3L &&
              length(unique(result$data_id)) == 1L &&
              all(result$seed_index == 1L),
            "Stage 6A probe/smoke manifest must contain three registered fits.")
  result
}

freeze_selection <- function(expected_fits = PP_EXPECTED_FITS, smoke = FALSE) {
  require_output_root()
  manifests <- verify_environment(FALSE)
  manifest <- if (smoke) {
    stage6a_probe_manifest(manifests$fits)
  } else manifests$fits
  pp_assert(nrow(manifest) == expected_fits, "Unexpected selection fit count.")
  terminal <- collect_terminals(manifest)
  rows <- terminal$rows
  pp_assert(all(rows$terminal_status != "error") &&
              all(rows$objective_eligible) &&
              !any(rows$truth_used_for_fit_stopping_or_selection) &&
              !any(rows$continuation_used),
            "Terminal pool contains an error, ineligible, truth-bearing, or continued fit.")
  winners <- pp_select_g12_winners(rows)
  expected_winners <- length(unique(manifest$selection_stratum_id))
  pp_assert(nrow(winners) == expected_winners &&
              all(table(winners$selection_stratum_id) == 1L),
            "Truth-free G12 winner table is invalid.")
  selection_dir <- file.path(output_root, "truth_free_selection")
  pp_assert(!dir.exists(selection_dir),
            "Truth-free selection directory already exists; refusing overwrite.")
  dir.create(selection_dir, recursive = TRUE, mode = "0700")
  endpoints_path <- file.path(selection_dir, "ALL_TRUTH_FREE_ENDPOINTS.csv")
  winners_path <- file.path(selection_dir, "G12_TRUTH_FREE_WINNERS.csv")
  pp_write_csv(rows, endpoints_path)
  pp_write_csv(winners, winners_path)
  selection_hash <- pp_sha256(winners_path)
  pp_write_lines(c(
    "status=TRUTH_FREE_G12_SELECTIONS_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("endpoints=", nrow(rows)),
    paste0("g12_winners=", nrow(winners)),
    paste0("selection_sha256=", selection_hash),
    "selection_used_truth=FALSE",
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_fit_config",
    "tie_break_rule=fit_id_lexicographic",
    "cross_model_ELBO_comparison=FALSE",
    "cross_dimension_ELBO_comparison=FALSE",
    "truth_unseal_authorized=FALSE", "continuation_started=FALSE",
    paste0("smoke=", toupper(as.character(smoke))),
    "formal_paper_mc_result=FALSE"
  ), file.path(selection_dir, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  pp_write_lines(c(
    "status=FITS_COMPLETE_AND_TRUTH_FREE_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("fits=", nrow(rows)), paste0("objective_eligible=", sum(rows$objective_eligible)),
    paste0("g12_winners=", nrow(winners)), "outer_workers=3",
    "per_fit_n_cpus=1", "truth_unseal_authorized=FALSE",
    "evaluation_started=FALSE", "continuation_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  pp_write_lines(c(
    "status=WAITING_FOR_EXPLICIT_TRUTH_UNSEAL_AUTHORIZATION",
    paste0("selection_sha256=", selection_hash),
    "truth_unseal_authorized=FALSE", "evaluation_started=FALSE",
    "new_fits_authorized=FALSE", "continuation_authorized=FALSE"
  ), file.path(output_root, "WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt"))
  cat("FITS_COMPLETE_AND_SELECTION_FROZEN fits=", nrow(rows),
      " winners=", nrow(winners), " selection_sha256=", selection_hash,
      "\n", sep = "")
  invisible(list(rows = rows, winners = winners, hash = selection_hash))
}

run_capacity_probe <- function() {
  require_output_root()
  manifests <- verify_environment(FALSE)
  pp_assert(.Platform$OS.type == "unix",
            "The concurrent capacity probe requires a Linux/Unix host.")
  pp_assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
            "Generate and seal data before the capacity probe.")
  pp_assert(!file.exists(file.path(output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt")) &&
              !dir.exists(file.path(output_root, "truth_free_selection")),
            "Capacity probing is forbidden after selection or truth unseal.")
  probe <- stage6a_probe_manifest(manifests$fits)
  log_dir <- file.path(output_root, "capacity_probe_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  time_binary <- Sys.which("time")
  if (!nzchar(time_binary) && file.exists("/usr/bin/time")) {
    time_binary <- "/usr/bin/time"
  }
  pp_assert(nzchar(time_binary), "GNU time is required for the capacity probe.")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    time_log <- file.path(log_dir, paste0(key, ".time.txt"))
    started <- Sys.time()
    status <- system2(
      time_binary,
      c("-v", "-o", shQuote(time_log), "Rscript", shQuote(PP_SCRIPT_FILE),
        "--action=fit", paste0("--fit-id=", fit_id),
        paste0("--output-root=", shQuote(output_root))),
      stdout = stdout, stderr = stderr
    )
    ended <- Sys.time()
    time_lines <- if (file.exists(time_log)) readLines(time_log, warn = FALSE) else
      character()
    rss_line <- grep("Maximum resident set size", time_lines, value = TRUE)
    max_rss_kb <- if (length(rss_line) == 1L) {
      suppressWarnings(as.numeric(sub("^.*:[[:space:]]*", "", rss_line)))
    } else NA_real_
    data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      max_rss_kb = max_rss_kb, exit_status = as.integer(status),
      stringsAsFactors = FALSE
    )
  }
  probe_results <- parallel::mclapply(
    probe$fit_id, run_one, mc.cores = 3L, mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  process <- do.call(rbind, probe_results)
  terminal <- collect_terminals(probe)$rows
  result <- merge(probe[, c("fit_id", "fit_config_id")], process,
                  by = "fit_id", sort = FALSE)
  result <- merge(result, terminal[, c(
    "fit_id", "terminal_status", "objective_eligible", "warning_count",
    "elapsed_seconds", "peak_memory_bytes"
  )], by = "fit_id", sort = FALSE, suffixes = c("_process", "_fit"))
  result <- result[match(probe$fit_id, result$fit_id), , drop = FALSE]
  df_lines <- system2("df", c("-Pk", shQuote(output_root)), stdout = TRUE)
  df_fields <- strsplit(trimws(tail(df_lines, 1L)), "[[:space:]]+")[[1L]]
  available_disk_kb <- suppressWarnings(as.numeric(df_fields[[4L]]))
  concurrent_rss_bytes <- sum(result$max_rss_kb, na.rm = FALSE) * 1024
  projected_wall_hours <- 1.10 *
    (12 * length(PP_DATA_IDS) / PP_WORKERS) *
    sum(result$elapsed_seconds_fit) / 3600
  checks <- c(
    three_processes_exit_zero = all(result$exit_status == 0L),
    three_fit_endpoints_complete = all(result$terminal_status == "fixed_400_complete"),
    three_objectives_eligible = all(result$objective_eligible),
    rss_observed_for_all = all(is.finite(result$max_rss_kb)),
    concurrent_rss_below_six_gib = is.finite(concurrent_rss_bytes) &&
      concurrent_rss_bytes < 6 * 1024^3,
    disk_free_above_ten_gib = is.finite(available_disk_kb) &&
      available_disk_kb > 10 * 1024^2
  )
  pp_write_csv(result, file.path(output_root, "CAPACITY_PROBE_RESULTS.csv"))
  pp_write_csv(data.frame(
    check_id = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE
  ), file.path(output_root, "CAPACITY_PROBE_QC.csv"))
  pp_write_lines(c(
    paste0("status=", if (all(checks)) "PASS" else "FAIL"),
    paste0("completed_utc=", pp_iso_time()),
    "registered_probe_fits=3", "probe_fits_are_part_of_36=TRUE",
    paste0("concurrent_rss_bytes_upper_sum=", format(concurrent_rss_bytes,
                                                     scientific = FALSE)),
    paste0("available_disk_kb=", format(available_disk_kb, scientific = FALSE)),
    paste0("projected_wall_hours_three_workers=",
           format(projected_wall_hours, digits = 6)),
    "truth_read=FALSE", "selection_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_root, if (all(checks))
    "CAPACITY_PROBE_COMPLETE.txt" else "CAPACITY_PROBE_ERROR.txt"))
  pp_assert(all(checks), paste0(
    "Capacity probe failed: ", paste(names(checks)[!checks], collapse = ";")
  ))
  cat("CAPACITY_PROBE_PASS fits=3 projected_wall_hours=",
      format(projected_wall_hours, digits = 6), "\n", sep = "")
  invisible(result)
}

run_supervisor <- function() {
  require_output_root()
  manifests <- verify_environment(FALSE)
  pp_assert(.Platform$OS.type == "unix",
            "The full supervisor requires Unix mclapply; use smoke on Windows.")
  pp_assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
            "Generate and seal data before starting fits.")
  pp_assert(file.exists(file.path(output_root, "CAPACITY_PROBE_COMPLETE.txt")),
            "The registered three-fit capacity probe has not passed.")
  pp_assert(!file.exists(file.path(output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt")),
            "Refusing to start fits after truth-unseal authorization.")
  manifest <- manifests$fits
  log_dir <- file.path(output_root, "supervisor_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  write_status <- function(state, detail = "") {
    complete <- length(list.files(file.path(output_root, "fits"),
                                  pattern = "^FIT_COMPLETE[.]txt$",
                                  recursive = TRUE))
    errors <- length(list.files(file.path(output_root, "fits"),
                                pattern = "^FIT_ERROR[.]txt$",
                                recursive = TRUE))
    status <- data.frame(
      timestamp_utc = pp_iso_time(), state = state, detail = detail,
      total_fits = nrow(manifest), complete_fits = complete,
      error_fits = errors, outer_workers = PP_WORKERS,
      per_fit_n_cpus = 1L, truth_read_by_fit_or_selection = FALSE,
      formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
    )
    pp_write_csv(status, file.path(output_root, "CURRENT_STATUS.csv"))
  }
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
    record <- data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      exit_status = as.integer(status), stringsAsFactors = FALSE
    )
    pp_write_csv(record, file.path(log_dir, paste0(key, ".exit.csv")))
    list(fit_id = fit_id, status = as.integer(status))
  }
  write_status("RUNNING", "36 registered fits; three completed probes reused; no selective rerun")
  results <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = PP_WORKERS,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    write_status("ERROR", paste0("retained=", paste(failed, collapse = ";")))
    pp_abort("Fit errors retained without selective rerun: ",
             paste(failed, collapse = ", "))
  }
  freeze_selection(PP_EXPECTED_FITS, FALSE)
  write_status("WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION",
               "36/36 complete; three within-configuration G12 winners frozen")
  invisible(TRUE)
}

validate_unseal_authorization <- function(smoke = FALSE) {
  selection_marker <- file.path(output_root, "truth_free_selection",
                                "TRUTH_FREE_SELECTION_FROZEN.txt")
  winners_path <- file.path(output_root, "truth_free_selection",
                            "G12_TRUTH_FREE_WINNERS.csv")
  pp_assert(file.exists(selection_marker) && file.exists(winners_path),
            "Truth-free selection freeze is missing.")
  hash <- pp_sha256(winners_path)
  path <- file.path(output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt")
  pp_assert(file.exists(path), "Explicit truth-unseal authorization is missing.")
  lines <- readLines(path, warn = FALSE)
  required <- c(
    "status=TRUTH_UNSEAL_AND_EVALUATION_AUTHORIZED",
    paste0("selection_sha256=", hash),
    "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
    "new_fits_authorized=FALSE", "continuation_authorized=FALSE",
    paste0("smoke=", toupper(as.character(smoke)))
  )
  pp_assert(all(required %in% lines),
            "Truth-unseal authorization does not bind the frozen selection.")
  hash
}

evaluate_all <- function(smoke = FALSE) {
  require_output_root()
  manifests <- verify_environment(TRUE)
  selection_hash <- validate_unseal_authorization(smoke)
  manifest <- if (smoke) {
    stage6a_probe_manifest(manifests$fits)
  } else manifests$fits
  pp_assert(!dir.exists(file.path(output_root, "evaluation")),
            "Evaluation directory already exists; refusing overwrite.")
  terminal <- collect_terminals(manifest)
  pp_assert(all(terminal$rows$objective_eligible) &&
              !any(terminal$rows$truth_used_for_fit_stopping_or_selection),
            "Evaluation terminal pool is ineligible or truth-bearing.")
  truths <- setNames(lapply(unique(manifest$data_id), function(data_id) {
    path <- file.path(output_root, "data", data_id, "sealed_truth",
                      "truth_bundle.rds")
    pp_assert(file.exists(path), paste0("Missing sealed truth: ", data_id))
    truth <- readRDS(path)
    pp_assert(identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
                identical(truth$data_id, data_id),
              paste0("Truth bundle identity failed: ", data_id))
    truth
  }), unique(manifest$data_id))
  evaluation_root <- file.path(output_root, "evaluation")
  dir.create(evaluation_root, recursive = TRUE, mode = "0700")
  evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
  statuses <- vector("list", nrow(manifest))
  compact <- vector("list", nrow(manifest))
  metric_checks <- vector("list", nrow(manifest))
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    fit_id <- row$fit_id[[1L]]
    fit <- readRDS(file.path(output_root, "fits", fit_id, "fit.rds"))
    record <- readRDS(file.path(output_root, "fits", fit_id,
                               "terminal_record.rds"))
    started <- Sys.time()
    warning_rows <- list()
    result <- tryCatch(withCallingHandlers({
      pp_evaluate_fit(
        fit, record, truths[[row$data_id[[1L]]]],
        evaluator_environment
      )
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
    output_dir <- file.path(evaluation_root, fit_id)
    dir.create(output_dir, recursive = TRUE, mode = "0700")
    warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
      data.frame(class = character(), message = character(), call = character(),
                 stringsAsFactors = FALSE)
    pp_write_csv(warnings, file.path(output_dir, "EVALUATION_WARNINGS.csv"))
    if (inherits(result, "condition")) {
      pp_write_lines(c(
        "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
        paste0("fit_id=", fit_id), paste0("message=", conditionMessage(result))
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
    base <- data.frame(
      fit_id = fit_id, data_id = row$data_id[[1L]],
      scenario_id = row$scenario_id[[1L]], method_id = row$method_id[[1L]],
      fit_config_id = row$fit_config_id[[1L]],
      selection_stratum_id = row$selection_stratum_id[[1L]],
      seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
      fit_L_f = row$fit_L_f[[1L]], fit_L_s_1 = row$fit_L_s_1[[1L]],
      fit_L_s_2 = row$fit_L_s_2[[1L]], fit_M_f = row$fit_M_f[[1L]],
      fit_M_s_1 = row$fit_M_s_1[[1L]], fit_M_s_2 = row$fit_M_s_2[[1L]],
      selected_g12_winner = fit_id %in% pp_read_csv(file.path(
        output_root, "truth_free_selection", "G12_TRUTH_FREE_WINNERS.csv"
      ))$fit_id,
      stringsAsFactors = FALSE
    )
    pp_assert(is.data.frame(result$summary) && nrow(result$summary) == 1L,
              paste0("Malformed multi evaluator summary: ", fit_id))
    science <- cbind(
      result$summary,
      result$process_contribution_summary,
      result$loading_scale_summary
    )
    compact[[index]] <- cbind(base, science)
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
    tables_present <- all(expected_tables %in% names(result)) &&
      all(vapply(result[expected_tables], function(value) {
        is.data.frame(value) && nrow(value) > 0L
      }, logical(1L)))
    checkpoint_complete <- if (smoke) TRUE else
      identical(names(fit$practical_checkpoints), c("380", "400"))
    metric_checks[[index]] <- data.frame(
      fit_id = fit_id, method_id = row$method_id[[1L]],
      expected_table_count = length(expected_tables),
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
      "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
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
    pp_abort("Evaluation errors retained; no selective rerun.")
  }
  all_science <- pp_rbind_fill(compact)
  metric_qc <- do.call(rbind, metric_checks)
  pp_assert(nrow(metric_qc) == nrow(manifest) &&
              all(metric_qc$metric_contract_passed),
            "One or more evaluation metric contracts are incomplete.")
  pp_write_csv(metric_qc,
               file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
  pp_write_csv(all_science,
               file.path(evaluation_root, "ALL_STAGE6A_SCIENTIFIC_RESULTS.csv"))
  winner_ids <- pp_read_csv(file.path(
    output_root, "truth_free_selection", "G12_TRUTH_FREE_WINNERS.csv"
  ))$fit_id
  paper_endpoints <- all_science[all_science$fit_id %in% winner_ids, , drop = FALSE]
  pp_write_csv(paper_endpoints,
               file.path(evaluation_root, "STAGE6A_WINNER_RESULTS.csv"))
  pp_write_lines(c(
    "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
    paste0("fits_evaluated=", nrow(manifest)), "evaluation_errors=0",
    paste0("g12_winners=", length(winner_ids)), "pooled_endpoints=0",
    paste0("selection_sha256=", selection_hash),
    "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
    "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
  cat("EVALUATION_COMPLETE fits=", nrow(manifest), " errors=0\n", sep = "")
  invisible(all_science)
}

summarize_run <- function(smoke = FALSE) {
  require_output_root()
  manifests <- verify_environment(FALSE)
  manifest <- if (smoke) {
    stage6a_probe_manifest(manifests$fits)
  } else manifests$fits
  terminal <- collect_terminals(manifest)$rows
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  sizes <- file.info(files)$size
  runtime <- do.call(rbind, lapply(split(terminal, terminal$fit_config_id), function(x) {
    data.frame(
      method_id = x$method_id[[1L]], fits = nrow(x),
      fit_config_id = x$fit_config_id[[1L]],
      elapsed_median_seconds = stats::median(x$elapsed_seconds, na.rm = TRUE),
      elapsed_max_seconds = max(x$elapsed_seconds, na.rm = TRUE),
      peak_memory_max_bytes = max(x$peak_memory_bytes, na.rm = TRUE),
      warning_total = sum(x$warning_count),
      objective_eligible = sum(x$objective_eligible),
      stringsAsFactors = FALSE
    )
  }))
  pp_write_csv(runtime, file.path(output_root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
  resource <- data.frame(
    generated_at_utc = pp_iso_time(), files = length(files),
    disk_bytes = sum(sizes[is.finite(sizes)], na.rm = TRUE),
    registered_fits = nrow(terminal), complete_fits = sum(terminal$terminal_status != "error"),
    objective_eligible = sum(terminal$objective_eligible),
    worker_count = if (smoke) 1L else PP_WORKERS,
    n_cpus_per_fit = 1L, smoke = smoke, formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
  pp_write_csv(resource, file.path(output_root, "RESOURCE_SUMMARY.csv"))
  cat("SUMMARY_COMPLETE fits=", nrow(terminal), " disk_bytes=",
      resource$disk_bytes, "\n", sep = "")
  invisible(list(runtime = runtime, resource = resource))
}

run_smoke <- function() {
  verify_environment(TRUE)
  if (!nzchar(output_root)) {
    output_root <<- file.path(tempdir(), paste0("g12_stage6a_smoke_",
                                               format(Sys.time(), "%Y%m%d%H%M%S")))
  }
  dir.create(output_root, recursive = TRUE, mode = "0700")
  generate_all(TRUE)
  manifests <- pp_load_manifests()
  smoke_manifest <- stage6a_probe_manifest(manifests$fits)
  for (fit_id in smoke_manifest$fit_id) run_fit(fit_id, TRUE)
  selection <- freeze_selection(nrow(smoke_manifest), TRUE)
  pp_write_lines(c(
    "status=TRUTH_UNSEAL_AND_EVALUATION_AUTHORIZED",
    paste0("selection_sha256=", selection$hash),
    "truth_unseal_authorized=TRUE", "evaluation_authorized=TRUE",
    "new_fits_authorized=FALSE", "continuation_authorized=FALSE",
    "smoke=TRUE", "authorization_scope=nonformal_micro_smoke_only"
  ), file.path(output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))
  evaluate_all(TRUE)
  summarize_run(TRUE)
  pp_write_lines(c(
    "status=PASS", paste0("completed_utc=", pp_iso_time()),
    "data_sets_generated=1", "data_sets_fitted=1", "fits=3",
    "g12_winners=3", "pooled_endpoints=0",
    "truth_free_selection_preceded_evaluation=TRUE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "SMOKE_COMPLETE.txt"))
  cat("STAGE6A_MICRO_SMOKE_PASS output_root=", output_root, "\n", sep = "")
  invisible(output_root)
}

if (action == "check") {
  manifests <- verify_environment(FALSE)
  cat("CHECK_PASS data=", nrow(manifests$data), " fits=", nrow(manifests$fits),
      " workers=", PP_WORKERS, "\n", sep = "")
} else if (action == "generate") {
  generate_all(smoke)
} else if (action == "fit") {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required for action=fit.")
  run_fit(fit_id, smoke)
} else if (action == "run") {
  pp_assert(!smoke, "Use action=smoke for the nonformal smoke.")
  run_supervisor()
} else if (action == "capacity") {
  pp_assert(!smoke, "Capacity probing uses the registered full-size probe fits.")
  run_capacity_probe()
} else if (action == "select") {
  freeze_selection(if (smoke) 3L else PP_EXPECTED_FITS, smoke)
} else if (action == "evaluate") {
  evaluate_all(smoke)
} else if (action == "summarize") {
  summarize_run(smoke)
} else if (action == "smoke") {
  run_smoke()
} else {
  pp_abort("Unknown action: ", action)
}
