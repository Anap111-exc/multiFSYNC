#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

cli <- uc_parse_cli(commandArgs(trailingOnly = TRUE))
fit_id <- cli$`fit-id`
check_only <- identical(tolower(cli$`check-only` %||% "false"), "true")
uc_assert(is.character(fit_id) && length(fit_id) == 1L && nzchar(fit_id),
          "--fit-id is required.")

s4c_verify_environment(load_runtime = TRUE, package_role = "development")
manifest <- uc_read_csv(file.path(UC_ROOT, "FIT_MANIFEST.csv"))
s4c_validate_fit_manifest(manifest)
row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
uc_assert(nrow(row) == 1L, "fit_id does not join to exactly one manifest row.")
observation <- uc_load_observation(row$data_id[[1L]])

if (check_only) {
  cat("WORKER_CHECK_ONLY_PASS fit_id=", fit_id,
      " fixed_t1=400 truth=0\n", sep = "")
  quit(save = "no", status = 0L)
}

output_dir <- file.path(UC_ROOT, "fits", fit_id)
if (dir.exists(output_dir)) {
  complete <- file.exists(file.path(output_dir, "FIT_COMPLETE.txt"))
  terminal <- file.exists(file.path(output_dir, "terminal_record.rds"))
  uc_assert(complete && terminal,
            paste0("Incomplete/error pre-existing fit directory: ", fit_id))
  record <- readRDS(file.path(output_dir, "terminal_record.rds"))
  uc_assert(identical(record$terminal_status, "fixed_400_complete"),
            "Pre-existing terminal record is not reusable fixed-400 output.")
  cat("FIT_ALREADY_COMPLETE fit_id=", fit_id, "\n", sep = "")
  quit(save = "no", status = 0L)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")

spec <- s4c_fit_control()
fit_function <- getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score")
p <- length(observation$Y[[1L]][[1L]])
call_args <- list(
  Y = observation$Y, Z = observation$Z,
  time_obs = observation$time_obs,
  L_f = observation$dimensions$L_f,
  L_s = observation$dimensions$L_s,
  M_f = observation$dimensions$M_f,
  M_s = observation$dimensions$M_s,
  K = observation$dimensions$K,
  anneal = spec$anneal, list_hyper = NULL,
  n_g = spec$n_g, time_g = NULL,
  tol_abs = spec$tol_abs, tol_rel = spec$tol_rel,
  maxit = spec$maxit, n_cpus = 1L, verbose = FALSE,
  seed = as.integer(row$fit_seed[[1L]]),
  bool_scale = FALSE, bool_var_spec_prob = FALSE,
  d_0 = as.integer(p), convergence_rule = spec$convergence_rule,
  lambda_orth = 0, practical_control = spec$practical_control,
  initialization = "random",
  initialization_control = spec$initialization_control,
  continuation_state = NULL,
  pre_score_sweeps = 1L,
  trace_sweeps = 1L,
  function_initialization = "gram_unit_energy"
)

captured <- c2_capture_conditions(
  s4c_with_private_traces(
    s4c_trace_control(), do.call(fit_function, call_args)
  ), phase = "fixed_400_main_fit"
)
trace_result <- captured$value
fit <- if (is.null(trace_result)) NULL else trace_result$fit
direction_snapshots <- if (is.null(trace_result)) list() else
  trace_result$direction_snapshots
if (!is.null(fit) && exists("c2_annotate_fit_warnings", mode = "function")) {
  captured$warnings <- c2_annotate_fit_warnings(captured$warnings, fit)
}

validation_error <- NULL
objective <- list(
  eligible = FALSE, checks = logical(), invalid_reasons = "terminal_error"
)
if (is.null(captured$error)) {
  validation_error <- tryCatch({
    uc_assert(is.list(fit) && is.list(direction_snapshots),
              "Fixed-400 run did not return fit and direction snapshots.")
    control <- fit$driver_diagnostic_control
    uc_assert(is.list(control) && control$pre_score_sweeps == 1L &&
                control$dense_gate_sweeps == 0L &&
                identical(control$random_scale_calibration, "function"),
              "Returned diagnostic control differs from G12 manifest.")
    interface <- fit$pre_score_interface
    uc_assert(is.list(interface) &&
                identical(interface$interface, "bayesSYNC_multi_pre_score") &&
                identical(interface$function_initialization,
                          "gram_unit_energy") &&
                identical(interface$random_scale_calibration, "function") &&
                interface$pre_score_sweeps == 1L,
              "Returned public interface provenance differs from G12 route.")
    calibration <- fit$random_scale_calibration_diagnostics
    uc_assert(is.data.frame(calibration) && nrow(calibration) == 6L &&
                all(calibration$success) &&
                max(abs(calibration$after - 1)) < 1e-10,
              "G12 unit Gram-energy calibration is incomplete.")
    uc_assert(
      identical(fit$initialization, "random") &&
        fit$annealing_sweeps == 99L && fit$n_g == 51L && fit$d_0 == p &&
        fit$n_cpus_used == 1L && identical(as.numeric(fit$lambda_orth), 0) &&
        identical(fit$convergence_rule, "practical") &&
        identical(fit$practical_control, spec$practical_control) &&
        fit$t1_sweeps == S4C_FIXED_T1 &&
        fit$practical_control$min_t1 == S4C_STOPPING_MIN_T1 &&
        fit$practical_control$consecutive == S4C_STOPPING_CONSECUTIVE &&
        fit$practical_control$max_t1 == S4C_FIXED_T1,
      "Returned fit violates fixed-400 G12 route."
    )
    pre <- fit$driver_trace$phase == "pre_score"
    uc_assert(any(pre) && all(fit$driver_trace$iteration[pre] == 1L),
              "Bounded first-sweep pre-score was not traced.")
    uc_assert(
      is.data.frame(fit$practical_diagnostics) &&
        nrow(fit$practical_diagnostics) == S4C_FIXED_T1 &&
        identical(as.integer(fit$practical_diagnostics$t1_sweep),
                  seq_len(S4C_FIXED_T1)) &&
        identical(names(fit$practical_checkpoints),
                  as.character(S4C_CHECKPOINTS)),
      "Fixed-400 practical trajectory is incomplete."
    )
    t1_trace <- fit$scale_trace[fit$scale_trace$stage == "t1", , drop = FALSE]
    uc_assert(
      sum(fit$scale_trace$stage == "initial") == 1L &&
        identical(as.integer(t1_trace$t1_sweep), S4C_CHECKPOINTS),
      "Fixed-400 scale trace is incomplete."
    )
    expected_snapshots <- as.character(c(0L, S4C_CHECKPOINTS))
    uc_assert(
      identical(names(direction_snapshots), expected_snapshots) &&
        all(vapply(direction_snapshots, function(item) {
          values <- c(
            item$shared_loading,
            unlist(item$specific_loading, use.names = FALSE),
            unlist(item$shared_feature, use.names = FALSE),
            unlist(item$specific_feature, use.names = FALSE),
            unlist(item$shared_score_mean, use.names = FALSE),
            unlist(item$specific_score_mean, use.names = FALSE),
            unlist(item$shared_score_kernel, use.names = FALSE),
            unlist(item$specific_score_kernel, use.names = FALSE),
            item$shared_contribution, item$specific_contribution
          )
          is.list(item) && item$cumulative_sweep %in% c(0L, S4C_CHECKPOINTS) &&
            all(is.finite(values))
        }, logical(1L))),
      "Direction snapshot sequence or finiteness is incomplete."
    )
    objective <- r_route_v3_candidate2_objective_status(fit)
    uc_assert(objective$eligible, paste0(
      "Fixed-400 endpoint is not objective eligible: ",
      paste(objective$invalid_reasons, collapse = ";")
    ))
    NULL
  }, error = function(condition) condition)
}

if (!is.null(validation_error)) {
  captured$error <- list(
    class = paste(class(validation_error), collapse = ";"),
    message = conditionMessage(validation_error),
    call = "stage4C fixed-400 return validation",
    trace_summary = "s4c_worker_20260826_v1"
  )
  fit <- NULL
  direction_snapshots <- list()
}

terminal_status <- if (is.null(captured$error)) "fixed_400_complete" else "error"
binding <- uc_read_csv(file.path(UC_ROOT, "SOURCE_BINDING.csv"))
record <- list(
  experiment_id = UC_EXPERIMENT_ID,
  fit_id = row$fit_id[[1L]], data_id = row$data_id[[1L]],
  scenario_id = row$scenario_id[[1L]],
  data_seed = as.integer(row$data_seed[[1L]]),
  method_id = "G12", function_initialization = "gram_unit_energy",
  calibration_mode = "function",
  seed_index = as.integer(row$seed_index[[1L]]),
  fit_seed = as.integer(row$fit_seed[[1L]]),
  dev_commit = binding$dev_commit[[1L]],
  parent_frozen_commit = UC_PARENT_COMMIT,
  parent_frozen_tag = UC_PARENT_TAG,
  started_at_utc = captured$started_at_utc,
  ended_at_utc = captured$ended_at_utc,
  elapsed_seconds = captured$elapsed_seconds,
  peak_memory_bytes = captured$peak_memory_bytes,
  terminal_status = terminal_status,
  error = captured$error %||% list(
    class = "", message = "", call = "", trace_summary = ""
  ),
  warnings = captured$warnings,
  actual_annealing_sweeps = if (is.null(fit)) NA_integer_ else
    as.integer(fit$annealing_sweeps),
  actual_T1_sweeps = if (is.null(fit)) NA_integer_ else
    as.integer(fit$t1_sweeps),
  objective_eligible = isTRUE(objective$eligible),
  objective_checks = objective$checks,
  objective_invalid_reasons = objective$invalid_reasons,
  final_elbo = if (is.null(fit)) NA_real_ else
    as.numeric(tail(fit$ELBO, 1L)),
  practical_converged = if (is.null(fit)) FALSE else
    isTRUE(fit$practical_converged),
  slow_case = if (is.null(fit)) TRUE else isTRUE(fit$slow_case),
  convergence_status = if (is.null(fit)) "error" else
    as.character(fit$convergence_status),
  convergence_reason = if (is.null(fit)) "terminal_error" else
    as.character(fit$convergence_reason),
  direction_snapshot_count = length(direction_snapshots),
  fixed_400_endpoint = TRUE,
  public_g12_stopping_profile_applied = TRUE,
  truth_used = FALSE, continuation_used = FALSE,
  automatic_800_used = FALSE, formal_v0lv_result = FALSE
)

saveRDS(record, file.path(output_dir, "terminal_record.rds"), version = 3)
if (!is.null(fit)) {
  saveRDS(fit, file.path(output_dir, "fit.rds"), version = 3)
  saveRDS(direction_snapshots,
          file.path(output_dir, "DIRECTION_SNAPSHOTS.rds"), version = 3)
  utils::write.csv(data.frame(
    t1_sweep = seq_along(fit$ELBO), elbo = as.numeric(fit$ELBO)
  ), file.path(output_dir, "ELBO_1_TO_400.csv"),
  row.names = FALSE, quote = TRUE)
  utils::write.csv(
    fit$practical_diagnostics,
    file.path(output_dir, "PRACTICAL_DIAGNOSTICS_1_TO_400.csv"),
    row.names = FALSE, quote = TRUE, na = ""
  )
  utils::write.csv(
    fit$scale_trace, file.path(output_dir, "SCALE_TRACE_0_TO_400.csv"),
    row.names = FALSE, quote = TRUE, na = ""
  )
  utils::write.csv(
    fit$random_scale_calibration_diagnostics,
    file.path(output_dir, "INITIALIZATION_DIAGNOSTICS.csv"),
    row.names = FALSE, quote = TRUE, na = ""
  )
}
if (is.data.frame(record$warnings)) {
  utils::write.csv(record$warnings, file.path(output_dir, "WARNINGS.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
}

if (terminal_status == "error") {
  writeLines(c(
    "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN",
    paste0("fit_id=", fit_id),
    paste0("message=", record$error$message),
    paste0("ended_utc=", uc_iso_time())
  ), file.path(output_dir, "FIT_ERROR.txt"), useBytes = TRUE)
  quit(save = "no", status = 2L)
}
writeLines(c(
  "status=COMPLETE", paste0("fit_id=", fit_id),
  "annealing_sweeps=99", "ordinary_T1_sweeps=400",
  "direction_snapshots=21", "public_g12_stopping_profile_applied=TRUE",
  paste0("practical_converged=", record$practical_converged),
  paste0("slow_case=", record$slow_case),
  "truth_used=FALSE", "continuation_used=FALSE",
  "automatic_800_used=FALSE",
  paste0("objective_eligible=", record$objective_eligible),
  paste0("ended_utc=", uc_iso_time())
), file.path(output_dir, "FIT_COMPLETE.txt"), useBytes = TRUE)
cat("FIT_COMPLETE fit_id=", fit_id,
    " fixed_t1=400 elapsed_seconds=", record$elapsed_seconds, "\n", sep = "")
