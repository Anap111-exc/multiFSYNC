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
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"), local = FALSE)

PP_EXPERIMENT_ID <- "STAGE6K_EARLY_SPARSITY_GATE_CAUSAL_AUDIT_V1_20260914"
stage6k_library <- pp_env("STAGE6K_R_LIB")
if (nzchar(stage6k_library)) Sys.setenv(STAGE6A_R_LIB = stage6k_library)
workers <- as.integer(pp_env("STAGE6K_WORKERS", "4"))
pp_assert(length(workers) == 1L && is.finite(workers) && workers %in% 1:4,
          "STAGE6K_WORKERS must be an integer from 1 to 4.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6K_OUTPUT_ROOT")
stage6g_output_root <- pp_env(
  "STAGE6K_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6j_output_root <- pp_env(
  "STAGE6K_STAGE6J_OUTPUT_ROOT",
  "/root/multifsync-results/stage6j_20260914_v1"
)
authorization_path <- file.path(
  script_root, "STAGE6K_EXECUTION_AUTHORIZATION_20260914.txt"
)
package_binding_path <- file.path(
  stage6g_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
)
trace_sweeps <- sort(unique(c(
  1:120, 129L, 149L, 199L, 299L, 399L, 499L
)))
report_trace_sweeps <- c(
  1L, 2L, 5L, 10L, 11L, 12L, 25L, 50L, 75L, 99L,
  100L, 120L, 199L, 299L, 399L, 499L
)

s6k_assert_authorized <- function() {
  pp_assert(file.exists(authorization_path),
            "Stage 6K execution authorization is missing.")
  lines <- readLines(authorization_path, warn = FALSE)
  required <- c(
    "status=AUTHORIZED_TRUTH_FREE_STAGE6K_EXECUTION",
    "new_fit_count=4", "paired_control_fit_count=4",
    "paired_fit_seed_count=4", "dense_gate_sweeps=10",
    "additional_gate_lengths_authorized=FALSE",
    "truth_available_to_fit=FALSE", "truth_available_to_stopping=FALSE",
    "truth_available_to_selection=FALSE",
    "post_freeze_scientific_evaluation_authorized=TRUE",
    "new_data_authorized=FALSE", "continuation_authorized=FALSE",
    "model_change_authorized=FALSE", "prior_change_authorized=FALSE",
    "cavi_change_authorized=FALSE", "elbo_change_authorized=FALSE",
    "public_interface_change_authorized=FALSE",
    "formal_paper_mc=FALSE", "selective_rerun_authorized=FALSE"
  )
  pp_assert(all(required %in% lines),
            "Stage 6K authorization has the wrong scope.")
  invisible(TRUE)
}

s6k_build_manifest <- function() {
  value <- data.frame(
    experiment_id = PP_EXPERIMENT_ID,
    fit_id = sprintf("s6k__G12_gate10__%02d", 1:4),
    data_id = "s6g_minimal_L111_M2",
    scenario_id = "s6g_minimal_L111_M2",
    method_id = "RANDOM_GRAM_UNIT_ENERGY_DENSE_GATE10",
    initialization_arm = "random_gram_unit_energy_dense_gate10",
    fit_config_id = "fixed_truth_L111_M2",
    selection_stratum_id = "stage6k_gate10_fixed_dimension",
    route_id = "random__gram_unit_energy__pre1__gate10__jaoua_anneal__fixed400",
    seed_index = 1:4,
    fit_seed = 84611101:84611104,
    paired_control_fit_id = sprintf("s6j__G12__%02d", 1:4),
    dense_gate_sweeps = 10L,
    fit_L_f = 1L, fit_L_s_1 = 1L, fit_L_s_2 = 1L,
    fit_M_f = "2", fit_M_s_1 = "2", fit_M_s_2 = "2",
    n_cpus = 1L, planned_annealing_sweeps = 99L,
    planned_ordinary_t1_sweeps = 400L,
    trace_total_sweeps = paste(trace_sweeps, collapse = ";"),
    truth_available_to_fit = FALSE,
    truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    scientific_evaluation_authorized_after_freeze = TRUE,
    cross_dimension_selection = FALSE,
    continuation = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
  value
}

s6k_validate_manifest <- function(manifest) {
  pp_assert(
    nrow(manifest) == 4L && !anyDuplicated(manifest$fit_id) &&
      identical(as.integer(manifest$fit_seed), 84611101:84611104) &&
      identical(as.integer(manifest$seed_index), 1:4) &&
      identical(as.character(manifest$paired_control_fit_id),
                sprintf("s6j__G12__%02d", 1:4)) &&
      all(manifest$dense_gate_sweeps == 10L) &&
      all(manifest$n_cpus == 1L) &&
      all(manifest$planned_annealing_sweeps == 99L) &&
      all(manifest$planned_ordinary_t1_sweeps == 400L) &&
      !any(manifest$truth_available_to_fit) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_selection) &&
      all(manifest$scientific_evaluation_authorized_after_freeze) &&
      !any(manifest$cross_dimension_selection) &&
      !any(manifest$continuation) &&
      !any(manifest$formal_paper_mc_result),
    "Stage 6K manifest validation failed."
  )
  invisible(TRUE)
}

s6k_source_paths <- function() {
  control_fit_ids <- sprintf("s6j__G12__%02d", 1:4)
  c(
    observation_bundle = file.path(
      stage6g_output_root, "data", "s6g_minimal_L111_M2",
      "observation_bundle.rds"
    ),
    stage6g_truth_free_complete = file.path(
      stage6g_output_root, "STAGE6G_TRUTH_FREE_COMPLETE.txt"
    ),
    stage6j_truth_free_complete = file.path(
      stage6j_output_root, "STAGE6J_TRUTH_FREE_COMPLETE.txt"
    ),
    stage6j_truth_free_endpoints = file.path(
      stage6j_output_root, "TRUTH_FREE_ENDPOINTS.csv"
    ),
    stats::setNames(file.path(
      stage6j_output_root, "fits", control_fit_ids, "fit.rds"
    ), paste0("paired_control_fit_", 1:4))
  )
}

s6k_verify_package_binding <- function() {
  pp_assert(file.exists(package_binding_path),
            "Stage 6G package source binding is missing.")
  source_root <- pp_env("STAGE6K_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
  binding <- pp_read_csv(package_binding_path)
  pp_assert(nrow(binding) == 29L && !anyDuplicated(binding$relative_path),
            "Stage 6G package binding is malformed.")
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

s6k_verify_environment <- function() {
  s6k_assert_authorized()
  paths <- s6k_source_paths()
  pp_assert(all(file.exists(paths)), paste0(
    "Required source artifact missing: ",
    paste(names(paths)[!file.exists(paths)], collapse = ";")
  ))
  pp_load_package()
  pp_source_runtime(FALSE)
  s6k_verify_package_binding()
  diagnostic <- getFromNamespace(
    ".bayesSYNC_multi_driver_diagnostic", "multiFSYNC"
  )
  pp_assert(is.function(diagnostic),
            "Private diagnostic runner is unavailable.")
  invisible(TRUE)
}

s6k_fit <- function(observation, row, smoke = FALSE) {
  if (smoke) {
    anneal <- c(1, 1.2, 3); maxit <- 14L; n_g <- 31L
    tol_abs <- 0; tol_rel <- 0; convergence_rule <- "parameters"
    practical_control <- NULL; requested_trace <- 1:14
  } else {
    anneal <- c(1, 1.9, 100); maxit <- 499L; n_g <- 51L
    tol_abs <- 1e-3; tol_rel <- 1e-5; convergence_rule <- "practical"
    practical_control <- multiFSYNC::g12_stopping_control()
    practical_control$checkpoints <- c(380L, 400L)
    requested_trace <- trace_sweeps
  }
  args <- list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = if (smoke) 3L else 5L,
    anneal = anneal, list_hyper = NULL, n_g = n_g, time_g = NULL,
    tol_abs = tol_abs, tol_rel = tol_rel, maxit = maxit,
    n_cpus = 1L, verbose = FALSE,
    seed = as.integer(row$fit_seed[[1L]]), bool_scale = FALSE,
    bool_var_spec_prob = FALSE, d_0 = 500L,
    convergence_rule = convergence_rule,
    practical_control = practical_control, lambda_orth = 0,
    initialization = "random",
    initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ),
    continuation_state = NULL,
    control = list(
      trace_sweeps = requested_trace,
      dense_gate_sweeps = 10L,
      random_scale_calibration = "function",
      pre_score_sweeps = 1L
    )
  )
  fit_function <- getFromNamespace(
    ".bayesSYNC_multi_driver_diagnostic", "multiFSYNC"
  )
  captured <- c2_capture_conditions(
    do.call(fit_function, args),
    phase = if (smoke) "stage6k_gate10_smoke" else "stage6k_gate10_fixed400"
  )
  fit <- captured$value
  objective <- list(eligible = FALSE, checks = logical(),
                    invalid_reasons = "terminal_error")
  validation_error <- NULL
  if (is.null(captured$error)) {
    validation_error <- tryCatch({
      control <- fit$driver_diagnostic_control
      calibration <- fit$random_scale_calibration_diagnostics
      trace <- fit$driver_trace
      checks <- c(
        returned_list = is.list(fit),
        random_initializer = identical(fit$initialization, "random"),
        one_cpu = identical(as.integer(fit$n_cpus_used), 1L),
        ordinary_elbo = identical(as.numeric(fit$lambda_orth), 0),
        diagnostic_control = is.list(control) &&
          identical(as.integer(control$dense_gate_sweeps), 10L) &&
          identical(as.integer(control$pre_score_sweeps), 1L) &&
          identical(control$random_scale_calibration, "function"),
        six_successful_g12_calibrations = is.data.frame(calibration) &&
          nrow(calibration) == 6L && all(calibration$success) &&
          max(abs(calibration$after - 1)) < 1e-10,
        trace_present = is.data.frame(trace) && nrow(trace) > 0L,
        gate_released = isTRUE(attr(trace, "dense_gate_released")),
        post_release_sweeps = as.integer(
          attr(trace, "post_release_sweeps")
        ) >= 1L,
        finite_elbo = length(fit$ELBO) > 0L && all(is.finite(fit$ELBO)),
        finite_loadings = all(is.finite(fit$mu_q_a)) &&
          all(is.finite(unlist(fit$mu_q_b_specific)))
      )
      if (!smoke) {
        objective <- r_route_v3_candidate2_objective_status(fit)
        checks <- c(checks,
          annealing_99 = identical(as.integer(fit$annealing_sweeps), 99L),
          t1_endpoint_valid = as.integer(fit$t1_sweeps) == 400L ||
            isTRUE(fit$practical_converged),
          endpoint_grid_51 = identical(as.integer(fit$n_g), 51L),
          requested_trace_complete = identical(
            sort(unique(trace$iteration)), trace_sweeps
          ),
          gate_active_through_10 = all(trace$gate_active[
            trace$phase == "main" & trace$iteration <= 10L
          ]),
          gate_inactive_after_10 = !any(trace$gate_active[
            trace$phase == "main" & trace$iteration > 10L
          ]),
          objective_eligible = objective$eligible
        )
      } else {
        objective <- list(
          eligible = all(checks), checks = checks,
          invalid_reasons = names(checks)[!checks]
        )
      }
      pp_assert(all(checks), paste(names(checks)[!checks], collapse = ";"))
      NULL
    }, error = function(condition) condition)
  }
  if (!is.null(validation_error)) {
    captured$error <- list(
      class = paste(class(validation_error), collapse = ";"),
      message = conditionMessage(validation_error),
      call = "s6k_fit", trace_summary = "route validation"
    )
    fit <- NULL
    objective$eligible <- FALSE
    objective$invalid_reasons <- conditionMessage(validation_error)
  }
  list(fit = fit, captured = captured, objective = objective)
}

s6k_prepare <- function() {
  pp_assert(nzchar(output_root) && !dir.exists(output_root),
            "Stage 6K output must be a new path.")
  s6k_verify_environment()
  manifest <- s6k_build_manifest()
  s6k_validate_manifest(manifest)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  pp_write_csv(manifest, file.path(output_root, "FIT_MANIFEST.csv"))
  sources <- s6k_source_paths()
  pp_write_csv(data.frame(
    role = names(sources), path = unname(sources),
    bytes = as.numeric(file.info(sources)$size),
    sha256 = vapply(sources, pp_sha256, character(1L)),
    mutable_by_stage6k = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "SOURCE_ARTIFACT_BINDING.csv"))
  pp_write_lines(c(
    "status=PREPARED", paste0("prepared_utc=", pp_iso_time()),
    "new_fits=4", "paired_readonly_controls=4", "dense_gate_sweeps=10",
    "truth_read_by_prepare=FALSE",
    "truth_available_to_fit_stopping_selection=FALSE",
    "scientific_evaluation_only_after_endpoint_freeze=TRUE",
    "model_prior_cavi_elbo_or_public_interface_changed=FALSE",
    "new_data=FALSE", "continuation=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6K_PREPARED.txt"))
  cat("STAGE6K_PREPARED fits=4 gate=10 truth_read=FALSE\n")
  invisible(TRUE)
}

s6k_load_manifest <- function() {
  value <- pp_read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
  s6k_validate_manifest(value)
  value
}

s6k_smoke <- function() {
  s6k_verify_environment()
  manifest <- if (nzchar(output_root) &&
                  file.exists(file.path(output_root, "FIT_MANIFEST.csv"))) {
    s6k_load_manifest()
  } else s6k_build_manifest()
  observation <- readRDS(s6k_source_paths()[["observation_bundle"]])
  result <- s6k_fit(observation, manifest[1L, , drop = FALSE], smoke = TRUE)
  status <- data.frame(
    fit_id = manifest$fit_id[[1L]],
    dense_gate_sweeps = 10L,
    gate_released = !is.null(result$fit) &&
      isTRUE(attr(result$fit$driver_trace, "dense_gate_released")),
    passed = !is.null(result$fit) && isTRUE(result$objective$eligible),
    warning_count = nrow(result$captured$warnings),
    error_message = result$captured$error$message %||% "",
    truth_read = FALSE, stringsAsFactors = FALSE
  )
  if (nzchar(output_root)) {
    pp_write_csv(status, file.path(output_root, "SMOKE_STATUS.csv"))
    pp_write_csv(result$captured$warnings,
                 file.path(output_root, "SMOKE_WARNINGS.csv"))
  }
  pp_assert(status$passed, "Stage 6K micro-smoke failed.")
  cat("STAGE6K_SMOKE_PASS gate_released=TRUE truth_read=FALSE\n")
  invisible(TRUE)
}

s6k_run_fit <- function(fit_id) {
  s6k_verify_environment()
  manifest <- s6k_load_manifest()
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Unknown Stage 6K fit_id.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(xor(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")),
                  file.exists(file.path(fit_dir, "FIT_ERROR.txt"))),
              "Existing Stage 6K fit directory is incomplete.")
    cat("FIT_ALREADY_TERMINAL fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  observation <- readRDS(s6k_source_paths()[["observation_bundle"]])
  result <- s6k_fit(observation, row, smoke = FALSE)
  if (!is.null(result$fit)) pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "STAGE6K_EARLY_SPARSITY_GATE_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$initialization_arm <- row$initialization_arm[[1L]]
  record$paired_control_fit_id <- row$paired_control_fit_id[[1L]]
  record$dense_gate_sweeps <- 10L
  record$truth_used_for_fit <- FALSE
  record$truth_used_for_stopping <- FALSE
  record$truth_used_for_selection <- FALSE
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  pp_write_csv(result$captured$warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", record$error$message %||% "unknown"),
      paste0("ended_utc=", pp_iso_time())
    ), file.path(fit_dir, "FIT_ERROR.txt"))
    cat("FIT_ERROR_RETAINED fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  fit <- result$fit
  pp_write_csv(data.frame(
    t1_sweep = seq_along(fit$ELBO), ordinary_elbo = as.numeric(fit$ELBO)
  ), file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  pp_write_csv(fit$driver_trace, file.path(fit_dir, "DENSE_DRIVER_TRACE.csv"))
  pp_write_csv(fit$driver_trace[
    fit$driver_trace$iteration %in% report_trace_sweeps, , drop = FALSE
  ], file.path(fit_dir, "REPORT_DRIVER_TRACE.csv"))
  pp_write_csv(fit$random_scale_calibration_diagnostics,
               file.path(fit_dir, "G12_CALIBRATION.csv"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("paired_control_fit_id=", row$paired_control_fit_id[[1L]]),
    "initialization=G12", "pre_score_sweeps=1", "dense_gate_sweeps=10",
    "dense_gate_released=TRUE", "truth_used_for_fit_stopping_selection=FALSE",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    paste0("annealing_sweeps=", fit$annealing_sweeps),
    paste0("ordinary_T1_sweeps=", fit$t1_sweeps), "n_cpus=1",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

s6k_status <- function(state, detail = "") {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_fits = 4L, complete_fits = complete, error_fits = errors,
    outer_workers = workers, per_fit_n_cpus = 1L,
    truth_unsealed = FALSE, formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

s6k_freeze_truth_free <- function() {
  s6k_verify_environment()
  manifest <- s6k_load_manifest()
  fit_dirs <- file.path(output_root, "fits", manifest$fit_id)
  pp_assert(all(file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))),
            "Not all Stage 6K fits are complete.")
  pp_assert(!any(file.exists(file.path(fit_dirs, "FIT_ERROR.txt"))),
            "Stage 6K errors are retained; truth-free freeze refused.")
  records <- lapply(fit_dirs, function(path) {
    readRDS(file.path(path, "terminal_record.rds"))
  })
  pp_assert(all(vapply(records, function(x) isTRUE(x$objective_eligible),
                       logical(1L))),
            "An ineligible endpoint cannot enter truth-free selection.")
  endpoints <- do.call(rbind, lapply(seq_along(records), function(index) {
    record <- records[[index]]
    fit <- readRDS(file.path(fit_dirs[[index]], "fit.rds"))
    data.frame(
      fit_id = record$fit_id,
      paired_control_fit_id = record$paired_control_fit_id,
      seed_index = record$seed_index, fit_seed = record$fit_seed,
      dense_gate_sweeps = 10L,
      final_ordinary_elbo = record$final_elbo,
      objective_eligible = record$objective_eligible,
      practical_converged = record$practical_converged,
      convergence_status = record$convergence_status,
      elapsed_seconds = record$elapsed_seconds,
      warning_count = nrow(record$warnings),
      annealing_sweeps = fit$annealing_sweeps,
      ordinary_t1_sweeps = fit$t1_sweeps,
      dense_gate_released = isTRUE(attr(
        fit$driver_trace, "dense_gate_released"
      )),
      post_release_sweeps = as.integer(attr(
        fit$driver_trace, "post_release_sweeps"
      )),
      truth_used_for_fit_stopping_or_selection = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  endpoints <- endpoints[order(endpoints$seed_index), , drop = FALSE]
  pp_write_csv(endpoints, file.path(output_root, "TRUTH_FREE_ENDPOINTS.csv"))
  winner <- endpoints[which.max(endpoints$final_ordinary_elbo), , drop = FALSE]
  winner$selection_rule <- "max_effective_ordinary_T1_ELBO_gate10_arm"
  winner$truth_used_for_selection <- FALSE
  pp_write_csv(winner, file.path(output_root, "TRUTH_FREE_SELECTION.csv"))
  pp_write_lines(c(
    "status=TRUTH_FREE_ENDPOINTS_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    "fits_complete=4", "fits_error=0", "objective_eligible=4",
    paste0("winner_fit_id=", winner$fit_id[[1L]]),
    "selection_rule=max_effective_ordinary_T1_ELBO_gate10_arm",
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "truth_unsealed=FALSE", "scientific_evaluation_run=FALSE",
    "new_data=FALSE", "continuation=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6K_TRUTH_FREE_COMPLETE.txt"))
  s6k_status("TRUTH_FREE_COMPLETE", "4 eligible gate10 endpoints; no truth read")
  cat("STAGE6K_TRUTH_FREE_COMPLETE winner=", winner$fit_id[[1L]],
      "\n", sep = "")
  invisible(TRUE)
}

s6k_run_supervisor <- function() {
  s6k_verify_environment()
  manifest <- s6k_load_manifest()
  pp_assert(.Platform$OS.type == "unix", "Supervisor requires Unix.")
  log_root <- file.path(output_root, "supervisor_logs")
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  s6k_status("RUNNING", "4 paired gate10 fits")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    started <- Sys.time()
    status <- tryCatch(system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=run-fit",
                   paste0("--fit-id=", fit_id),
                   paste0("--output-root=", output_root)),
      stdout = file.path(log_root, paste0(key, ".stdout.log")),
      stderr = file.path(log_root, paste0(key, ".stderr.log"))
    ), error = function(condition) 999L)
    ended <- Sys.time()
    data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      exit_status = as.integer(status), stringsAsFactors = FALSE
    )
  }
  rows <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = workers, mc.preschedule = FALSE
  )
  status <- do.call(rbind, rows)
  pp_write_csv(status, file.path(output_root, "SUPERVISOR_TASK_STATUS.csv"))
  failed <- status$fit_id[status$exit_status != 0L]
  if (length(failed)) {
    s6k_status("ERROR_RETAINED", paste(failed, collapse = ";"))
    pp_abort("Stage 6K child failures retained: ", paste(failed, collapse = ","))
  }
  s6k_freeze_truth_free()
  invisible(TRUE)
}

if (action == "check") {
  s6k_assert_authorized()
  manifest <- s6k_build_manifest()
  s6k_validate_manifest(manifest)
  cat("STAGE6K_CHECK_PASS fits=4 pairs=4 gate=10 truth_available=FALSE\n")
} else if (action == "preflight") {
  s6k_verify_environment()
  cat("STAGE6K_PREFLIGHT_PASS package_binding=29/29 paired_controls=4 truth_read=FALSE\n")
} else if (action == "prepare") {
  s6k_prepare()
} else if (action == "smoke") {
  s6k_smoke()
} else if (action == "run-fit") {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required.")
  s6k_run_fit(fit_id)
} else if (action == "run") {
  s6k_run_supervisor()
} else if (action == "freeze-truth-free") {
  s6k_freeze_truth_free()
} else {
  pp_abort("Unknown Stage 6K action: ", action)
}
