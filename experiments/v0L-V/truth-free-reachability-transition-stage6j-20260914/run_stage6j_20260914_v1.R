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

PP_EXPERIMENT_ID <- "STAGE6J_TRUTH_FREE_REACHABILITY_TRANSITION_V1_20260914"
stage6j_library <- pp_env("STAGE6J_R_LIB")
if (nzchar(stage6j_library)) Sys.setenv(STAGE6A_R_LIB = stage6j_library)
workers <- as.integer(pp_env("STAGE6J_WORKERS", "8"))
pp_assert(length(workers) == 1L && is.finite(workers) && workers %in% 1:8,
          "STAGE6J_WORKERS must be an integer from 1 to 8.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6J_OUTPUT_ROOT")
stage6g_output_root <- pp_env(
  "STAGE6J_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
authorization_path <- file.path(
  script_root, "STAGE6J_EXECUTION_AUTHORIZATION_20260914.txt"
)
package_binding_path <- file.path(
  stage6g_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
)
trace_sweeps <- sort(unique(c(1:120, 129L, 149L, 199L, 299L, 399L, 499L)))

s6j_assert_authorized <- function() {
  pp_assert(file.exists(authorization_path),
            "Stage 6J execution authorization is missing.")
  lines <- readLines(authorization_path, warn = FALSE)
  required <- c(
    "status=AUTHORIZED_TRUTH_FREE_STAGE6J_EXECUTION",
    "new_fit_count=8", "paired_fit_seed_count=4",
    "random_gram_unit_energy_fits=4", "residual_fpca_fits=4",
    "truth_available_to_fit=FALSE", "truth_available_to_stopping=FALSE",
    "truth_available_to_selection=FALSE", "truth_unseal_authorized=FALSE",
    "scientific_evaluation_authorized=FALSE",
    "model_change_authorized=FALSE", "prior_change_authorized=FALSE",
    "cavi_change_authorized=FALSE", "elbo_change_authorized=FALSE",
    "dimension_criterion_test=FALSE", "formal_paper_mc=FALSE",
    "source_outputs_mutable=FALSE"
  )
  pp_assert(all(required %in% lines),
            "Stage 6J authorization has the wrong scope.")
  invisible(TRUE)
}

s6j_build_manifest <- function() {
  base <- expand.grid(
    initialization_arm = c("random_gram_unit_energy", "residual_fpca"),
    seed_index = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  base <- base[order(base$initialization_arm, base$seed_index), , drop = FALSE]
  base$fit_id <- sprintf(
    "s6j__%s__%02d",
    ifelse(base$initialization_arm == "random_gram_unit_energy",
           "G12", "residualFPCA"),
    base$seed_index
  )
  base$fit_seed <- 84611100L + base$seed_index
  base$experiment_id <- PP_EXPERIMENT_ID
  base$data_id <- "s6g_minimal_L111_M2"
  base$scenario_id <- "s6g_minimal_L111_M2"
  base$method_id <- ifelse(
    base$initialization_arm == "random_gram_unit_energy",
    "RANDOM_GRAM_UNIT_ENERGY", "RESIDUAL_FPCA"
  )
  base$fit_config_id <- "fixed_truth_L111_M2"
  base$selection_stratum_id <- "stage6j_fixed_dimension_all_arms"
  base$route_id <- ifelse(
    base$initialization_arm == "random_gram_unit_energy",
    "random__gram_unit_energy__pre1__jaoua_anneal__fixed400",
    "residual_fpca__pre1__jaoua_anneal__fixed400"
  )
  base$fit_L_f <- 1L
  base$fit_L_s_1 <- 1L
  base$fit_L_s_2 <- 1L
  base$fit_M_f <- "2"
  base$fit_M_s_1 <- "2"
  base$fit_M_s_2 <- "2"
  base$n_cpus <- 1L
  base$planned_annealing_sweeps <- 99L
  base$planned_ordinary_t1_sweeps <- 400L
  base$trace_total_sweeps <- paste(trace_sweeps, collapse = ";")
  base$truth_available_to_fit <- FALSE
  base$truth_available_to_stopping <- FALSE
  base$truth_available_to_selection <- FALSE
  base$truth_unseal_authorized <- FALSE
  base$scientific_evaluation_authorized <- FALSE
  base$cross_dimension_selection <- FALSE
  base$formal_paper_mc_result <- FALSE
  base[, c(
    "experiment_id", "fit_id", "data_id", "scenario_id", "method_id",
    "initialization_arm", "fit_config_id", "selection_stratum_id",
    "route_id", "seed_index", "fit_seed", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2", "n_cpus",
    "planned_annealing_sweeps", "planned_ordinary_t1_sweeps",
    "trace_total_sweeps", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "truth_unseal_authorized", "scientific_evaluation_authorized",
    "cross_dimension_selection", "formal_paper_mc_result"
  )]
}

s6j_validate_manifest <- function(manifest) {
  pp_assert(
    nrow(manifest) == 8L && !anyDuplicated(manifest$fit_id) &&
      all(table(manifest$initialization_arm) == 4L) &&
      identical(sort(unique(manifest$fit_seed)), 84611101:84611104) &&
      all(table(manifest$fit_seed) == 2L) &&
      all(manifest$n_cpus == 1L) &&
      all(manifest$planned_annealing_sweeps == 99L) &&
      all(manifest$planned_ordinary_t1_sweeps == 400L) &&
      !any(manifest$truth_available_to_fit) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_selection) &&
      !any(manifest$truth_unseal_authorized) &&
      !any(manifest$scientific_evaluation_authorized) &&
      !any(manifest$cross_dimension_selection) &&
      !any(manifest$formal_paper_mc_result),
    "Stage 6J manifest validation failed."
  )
  invisible(TRUE)
}

s6j_source_paths <- function() list(
  observation = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2",
    "observation_bundle.rds"
  ),
  stage6g_complete = file.path(
    stage6g_output_root, "STAGE6G_TRUTH_FREE_COMPLETE.txt"
  )
)

s6j_verify_package_binding <- function() {
  pp_assert(file.exists(package_binding_path),
            "Stage 6G package source binding is missing.")
  source_root <- pp_env("STAGE6J_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
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

s6j_verify_environment <- function() {
  s6j_assert_authorized()
  paths <- s6j_source_paths()
  pp_assert(all(file.exists(unlist(paths))),
            "A required Stage 6G observation-only artifact is missing.")
  pp_load_package()
  pp_source_runtime(FALSE)
  s6j_verify_package_binding()
  invisible(TRUE)
}

s6j_fit <- function(observation, row, smoke = FALSE) {
  arm <- row$initialization_arm[[1L]]
  pp_assert(arm %in% c("random_gram_unit_energy", "residual_fpca"),
            "Unknown Stage 6J initialization arm.")
  if (smoke) {
    anneal <- c(1, 1.2, 3); maxit <- 6L; n_g <- 31L
    tol_abs <- 0; tol_rel <- 0; convergence_rule <- "parameters"
    practical_control <- NULL; requested_trace <- 1:6
  } else {
    anneal <- c(1, 1.9, 100); maxit <- 499L; n_g <- 51L
    tol_abs <- 1e-3; tol_rel <- 1e-5; convergence_rule <- "practical"
    practical_control <- multiFSYNC::g12_stopping_control()
    practical_control$checkpoints <- c(380L, 400L)
    requested_trace <- trace_sweeps
  }
  initialization <- if (arm == "residual_fpca") "residual_fpca" else "random"
  function_initialization <- if (arm == "residual_fpca") {
    "current_1_over_m"
  } else {
    "gram_unit_energy"
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
    initialization = initialization,
    initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ),
    continuation_state = NULL, pre_score_sweeps = 1L,
    trace_sweeps = requested_trace,
    function_initialization = function_initialization
  )
  captured <- c2_capture_conditions(do.call(
    getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score"), args
  ), phase = paste0("stage6j_", arm, if (smoke) "_smoke" else "_fixed400"))
  fit <- captured$value
  objective <- list(eligible = FALSE, checks = logical(),
                    invalid_reasons = "terminal_error")
  validation_error <- NULL
  if (is.null(captured$error)) {
    validation_error <- tryCatch({
      pp_assert(is.list(fit), "Stage 6J fit did not return a list.")
      checks <- c(
        correct_initializer = identical(fit$initialization, initialization),
        one_cpu = identical(as.integer(fit$n_cpus_used), 1L),
        ordinary_elbo = identical(as.numeric(fit$lambda_orth), 0),
        one_pre_score = identical(
          as.integer(fit$pre_score_interface$pre_score_sweeps), 1L
        ),
        finite_elbo = length(fit$ELBO) > 0L && all(is.finite(fit$ELBO)),
        finite_loadings = all(is.finite(fit$mu_q_a)) &&
          all(is.finite(unlist(fit$mu_q_b_specific)))
      )
      if (arm == "random_gram_unit_energy") {
        calibration <- fit$random_scale_calibration_diagnostics
        checks <- c(checks,
          g12_provenance = identical(
            fit$pre_score_interface$function_initialization,
            "gram_unit_energy"
          ),
          six_successful_calibrations = is.data.frame(calibration) &&
            nrow(calibration) == 6L && all(calibration$success) &&
            max(abs(calibration$after - 1)) < 1e-10
        )
      } else {
        init <- fit$initialization_diagnostics
        checks <- c(checks,
          residual_fpca_provenance = is.list(init) &&
            identical(init$requested_method, "residual_fpca") &&
            init$used_method %in% c("residual_fpca", "residual_fpca_partial") &&
            length(init$fallback_count) == 1L &&
            is.finite(init$fallback_count) && init$fallback_count >= 0L,
          no_random_calibration = is.data.frame(
            fit$random_scale_calibration_diagnostics
          ) && nrow(fit$random_scale_calibration_diagnostics) == 0L
        )
      }
      if (!smoke) {
        checks <- c(checks,
          annealing_99 = identical(as.integer(fit$annealing_sweeps), 99L),
          ordinary_t1_400 = identical(as.integer(fit$t1_sweeps), 400L),
          endpoint_grid_51 = identical(as.integer(fit$n_g), 51L),
          trace_complete = is.data.frame(fit$driver_trace) &&
            identical(sort(unique(fit$driver_trace$iteration)), trace_sweeps)
        )
        objective <- r_route_v3_candidate2_objective_status(fit)
        checks <- c(checks, objective_eligible = objective$eligible)
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
      call = "s6j_fit", trace_summary = "route validation"
    )
    fit <- NULL
    objective$eligible <- FALSE
    objective$invalid_reasons <- conditionMessage(validation_error)
  }
  list(fit = fit, captured = captured, objective = objective)
}

s6j_prepare <- function() {
  s6j_assert_authorized()
  pp_assert(nzchar(output_root) && !dir.exists(output_root),
            "Stage 6J output must be a new path.")
  s6j_verify_environment()
  manifest <- s6j_build_manifest()
  s6j_validate_manifest(manifest)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  pp_write_csv(manifest, file.path(output_root, "FIT_MANIFEST.csv"))
  sources <- s6j_source_paths()
  pp_write_csv(data.frame(
    role = names(sources), path = unlist(sources, use.names = FALSE),
    bytes = as.numeric(file.info(unlist(sources, use.names = FALSE))$size),
    sha256 = vapply(unlist(sources, use.names = FALSE), pp_sha256, character(1L)),
    mutable_by_stage6j = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "SOURCE_ARTIFACT_BINDING.csv"))
  pp_write_lines(c(
    "status=PREPARED", paste0("prepared_utc=", pp_iso_time()),
    "new_fits=8", "paired_fit_seeds=4", "truth_read_by_prepare=FALSE",
    "truth_available_to_fit_stopping_selection=FALSE",
    "scientific_evaluation_authorized=FALSE", "model_prior_cavi_elbo_changed=FALSE",
    "cross_dimension_selection=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6J_PREPARED.txt"))
  cat("STAGE6J_PREPARED fits=8 truth_read=FALSE\n")
  invisible(TRUE)
}

s6j_load_manifest <- function() {
  path <- file.path(output_root, "FIT_MANIFEST.csv")
  pp_assert(file.exists(path), "Stage 6J is not prepared.")
  value <- pp_read_csv(path)
  s6j_validate_manifest(value)
  value
}

s6j_smoke <- function() {
  s6j_verify_environment()
  manifest <- if (nzchar(output_root) &&
                  file.exists(file.path(output_root, "FIT_MANIFEST.csv"))) {
    s6j_load_manifest()
  } else s6j_build_manifest()
  observation <- readRDS(s6j_source_paths()$observation)
  rows <- manifest[match(c("random_gram_unit_energy", "residual_fpca"),
                         manifest$initialization_arm), , drop = FALSE]
  warning_rows <- list()
  results <- lapply(seq_len(nrow(rows)), function(i) {
    result <- s6j_fit(observation, rows[i, , drop = FALSE], smoke = TRUE)
    current_warnings <- result$captured$warnings
    if (nrow(current_warnings)) {
      current_warnings$initialization_arm <- rows$initialization_arm[[i]]
      warning_rows[[length(warning_rows) + 1L]] <<- current_warnings
    }
    data.frame(
      initialization_arm = rows$initialization_arm[[i]],
      passed = !is.null(result$fit) && isTRUE(result$objective$eligible),
      warning_count = nrow(result$captured$warnings),
      error_message = result$captured$error$message %||% "",
      stringsAsFactors = FALSE
    )
  })
  status <- do.call(rbind, results)
  if (nzchar(output_root)) {
    pp_write_csv(status, file.path(output_root, "SMOKE_STATUS.csv"))
    smoke_warnings <- if (length(warning_rows)) {
      do.call(rbind, warning_rows)
    } else {
      empty <- c2_empty_warnings()
      empty$initialization_arm <- character()
      empty
    }
    pp_write_csv(smoke_warnings,
                 file.path(output_root, "SMOKE_WARNINGS.csv"))
  }
  pp_assert(all(status$passed), "Stage 6J paired micro-smoke failed.")
  cat("STAGE6J_SMOKE_PASS arms=2 truth_read=FALSE\n")
  invisible(TRUE)
}

s6j_run_fit <- function(fit_id) {
  s6j_verify_environment()
  manifest <- s6j_load_manifest()
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Unknown Stage 6J fit_id.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(xor(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")),
                  file.exists(file.path(fit_dir, "FIT_ERROR.txt"))),
              "Existing Stage 6J fit directory is incomplete.")
    cat("FIT_ALREADY_TERMINAL fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  observation <- readRDS(s6j_source_paths()$observation)
  result <- s6j_fit(observation, row, smoke = FALSE)
  if (!is.null(result$fit)) pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "STAGE6J_TRUTH_FREE_REACHABILITY_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$initialization_arm <- row$initialization_arm[[1L]]
  record$truth_used_for_fit_stopping_or_selection <- FALSE
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
  init <- fit$initialization_diagnostics
  if (is.data.frame(init$factor_diagnostics)) {
    pp_write_csv(init$factor_diagnostics,
                 file.path(fit_dir, "INITIALIZATION_FACTOR_DIAGNOSTICS.csv"))
  }
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("initialization_arm=", row$initialization_arm[[1L]]),
    "truth_used_for_fit_stopping_or_selection=FALSE",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    "annealing_sweeps=99", "ordinary_T1_sweeps=400", "n_cpus=1",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

s6j_write_status <- function(state, detail = "") {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_fits = 8L, complete_fits = complete, error_fits = errors,
    outer_workers = workers, per_fit_n_cpus = 1L,
    truth_unsealed = FALSE, formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

s6j_freeze_truth_free <- function() {
  s6j_verify_environment()
  manifest <- s6j_load_manifest()
  fit_dirs <- file.path(output_root, "fits", manifest$fit_id)
  pp_assert(all(file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))),
            "Not all Stage 6J fits are complete.")
  pp_assert(!any(file.exists(file.path(fit_dirs, "FIT_ERROR.txt"))),
            "Stage 6J errors are retained; truth-free freeze refused.")
  records <- lapply(fit_dirs, function(path) readRDS(file.path(path, "terminal_record.rds")))
  pp_assert(all(vapply(records, function(x) isTRUE(x$objective_eligible), logical(1L))),
            "An ineligible endpoint cannot enter truth-free selection.")
  endpoints <- do.call(rbind, lapply(seq_along(records), function(i) {
    record <- records[[i]]
    fit <- readRDS(file.path(fit_dirs[[i]], "fit.rds"))
    init <- fit$initialization_diagnostics
    data.frame(
      fit_id = record$fit_id,
      initialization_arm = record$initialization_arm,
      seed_index = record$seed_index, fit_seed = record$fit_seed,
      final_ordinary_elbo = record$final_elbo,
      objective_eligible = record$objective_eligible,
      practical_converged = record$practical_converged,
      convergence_status = record$convergence_status,
      elapsed_seconds = record$elapsed_seconds,
      warning_count = nrow(record$warnings),
      initialization_used_method = init$used_method %||% "random",
      initialization_complete = isTRUE(init$complete),
      initialization_fallback_count = as.integer(init$fallback_count %||% 0L),
      trace_row_count = nrow(fit$driver_trace),
      trace_iteration_count = length(unique(fit$driver_trace$iteration)),
      truth_used_for_fit_stopping_or_selection = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  endpoints <- endpoints[order(endpoints$initialization_arm,
                               endpoints$seed_index), , drop = FALSE]
  pp_write_csv(endpoints, file.path(output_root, "TRUTH_FREE_ENDPOINTS.csv"))
  winner <- endpoints[which.max(endpoints$final_ordinary_elbo), , drop = FALSE]
  winner$selection_rule <- "max_effective_ordinary_T1_ELBO_fixed_dimension"
  winner$truth_used_for_selection <- FALSE
  pp_write_csv(winner, file.path(output_root, "TRUTH_FREE_SELECTION.csv"))
  pp_write_lines(c(
    "status=TRUTH_FREE_ENDPOINTS_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    "fits_complete=8", "fits_error=0", "objective_eligible=8",
    paste0("winner_fit_id=", winner$fit_id[[1L]]),
    paste0("winner_initialization_arm=", winner$initialization_arm[[1L]]),
    "selection_rule=max_effective_ordinary_T1_ELBO_fixed_dimension",
    "truth_used_for_fit_stopping_or_selection=FALSE",
    "truth_unsealed=FALSE", "scientific_evaluation_run=FALSE",
    "cross_dimension_selection=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6J_TRUTH_FREE_COMPLETE.txt"))
  s6j_write_status("TRUTH_FREE_COMPLETE", "8 eligible endpoints; no truth read")
  cat("STAGE6J_TRUTH_FREE_COMPLETE winner=", winner$fit_id[[1L]], "\n", sep = "")
  invisible(TRUE)
}

s6j_run_supervisor <- function() {
  s6j_verify_environment()
  manifest <- s6j_load_manifest()
  pp_assert(.Platform$OS.type == "unix", "Supervisor requires Unix.")
  log_root <- file.path(output_root, "supervisor_logs")
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  s6j_write_status("RUNNING", "8 paired truth-free fits")
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
    s6j_write_status("ERROR_RETAINED", paste(failed, collapse = ";"))
    pp_abort("Stage 6J child failures retained: ", paste(failed, collapse = ","))
  }
  s6j_freeze_truth_free()
  invisible(TRUE)
}

if (action == "check") {
  s6j_assert_authorized()
  manifest <- s6j_build_manifest()
  s6j_validate_manifest(manifest)
  cat("STAGE6J_CHECK_PASS fits=8 pairs=4 truth_available=FALSE\n")
} else if (action == "preflight") {
  s6j_verify_environment()
  cat("STAGE6J_PREFLIGHT_PASS package_binding=29/29 truth_read=FALSE\n")
} else if (action == "prepare") {
  s6j_prepare()
} else if (action == "smoke") {
  s6j_smoke()
} else if (action == "run-fit") {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required.")
  s6j_run_fit(fit_id)
} else if (action == "run") {
  s6j_run_supervisor()
} else if (action == "freeze-truth-free") {
  s6j_freeze_truth_free()
} else {
  pp_abort("Unknown Stage 6J action: ", action)
}
