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
source(file.path(script_root, "stage6i_helpers_20260913_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6I_CAUSAL_LOCALIZATION_V1_20260913"
stage6i_library <- pp_env("STAGE6I_R_LIB")
if (nzchar(stage6i_library)) Sys.setenv(STAGE6A_R_LIB = stage6i_library)
workers <- as.integer(pp_env("STAGE6I_WORKERS", "7"))
pp_assert(length(workers) == 1L && is.finite(workers) && workers %in% 1:7,
          "STAGE6I_WORKERS must be an integer from 1 to 7.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6I_OUTPUT_ROOT")
stage6g_output_root <- pp_env(
  "STAGE6I_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6h_output_root <- pp_env(
  "STAGE6I_STAGE6H_OUTPUT_ROOT",
  "/root/multifsync-results/stage6h_20260913_v1"
)
authorization_path <- file.path(
  script_root, "STAGE6I_EXECUTION_AUTHORIZATION_20260913.txt"
)
package_binding_path <- file.path(
  stage6g_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
)

s6i_assert_authorized <- function() {
  pp_assert(file.exists(authorization_path),
            "Stage 6I authorization is missing.")
  lines <- readLines(authorization_path, warn = FALSE)
  required <- c(
    "status=AUTHORIZED_POST_UNSEAL_CAUSAL_LOCALIZATION",
    "study_swap_g12_fits_authorized=TRUE",
    "study_swap_g12_fit_count=4",
    "truth_informed_oracle_fits_authorized=TRUE",
    "truth_informed_oracle_fit_count=6",
    "one_sweep_state_template_count=2",
    "existing_stage6h_full_data_fits_reused_read_only=4",
    "truth_available_to_oracle_initialization=TRUE",
    "truth_available_to_study_swap_g12_fit=FALSE",
    "cross_dimension_selection_authorized=FALSE",
    "new_dimension_criterion_authorized=FALSE",
    "model_or_prior_change_authorized=FALSE",
    "cavi_or_elbo_change_authorized=FALSE",
    "formal_paper_mc=FALSE", "source_outputs_mutable=FALSE"
  )
  pp_assert(all(required %in% lines),
            "Stage 6I authorization has the wrong scope.")
  invisible(TRUE)
}

s6i_build_manifest <- function() {
  swap <- data.frame(
    fit_id = sprintf("s6i__study_swap__G12_%02d", 1:4),
    fit_kind = "study_swap_g12", study_order = "swapped",
    oracle_mode = "none", replicate = 1:4,
    seed_index = 1:4, fit_seed = 84611101:84611104,
    state_seed = NA_integer_, source_fit_id = sprintf(
      "s6h_minimal_L111_M2__full_data_6_9__G12_%02d", 1:4
    ),
    route_id = "random__pre1__jaoua_anneal__G12__fixed400",
    truth_available_to_fit = FALSE,
    continuation_state_used = FALSE,
    stringsAsFactors = FALSE
  )
  oracle_base <- data.frame(
    study_order = rep(c("original", "swapped"), each = 3L),
    oracle_mode = rep(c("exact", "near", "near"), 2L),
    replicate = rep(c(1L, 1L, 2L), 2L),
    seed_index = rep(1:3, 2L),
    fit_seed = rep(84612101:84612103, 2L),
    state_seed = rep(84612001:84612003, 2L),
    stringsAsFactors = FALSE
  )
  oracle_base$fit_id <- with(oracle_base, sprintf(
    "s6i__oracle_%s__%s_%02d", oracle_mode, study_order, replicate
  ))
  oracle <- data.frame(
    fit_id = oracle_base$fit_id,
    fit_kind = "truth_neighborhood",
    study_order = oracle_base$study_order,
    oracle_mode = oracle_base$oracle_mode,
    replicate = oracle_base$replicate,
    seed_index = oracle_base$seed_index,
    fit_seed = oracle_base$fit_seed,
    state_seed = oracle_base$state_seed,
    source_fit_id = "",
    route_id = paste0(
      oracle_base$oracle_mode,
      "_projected_truth__ordinary_T1__fixed400"
    ),
    truth_available_to_fit = TRUE,
    continuation_state_used = TRUE,
    stringsAsFactors = FALSE
  )
  result <- rbind(swap, oracle)
  result$experiment_id <- PP_EXPERIMENT_ID
  result$data_id <- "s6g_minimal_L111_M2"
  result$scenario_id <- "s6g_minimal_L111_M2"
  result$method_id <- ifelse(result$fit_kind == "study_swap_g12",
                             "G12_STUDY_SWAP", "ORACLE_DIAGNOSTIC")
  result$fit_config_id <- "cand_L111_M2"
  result$selection_stratum_id <- paste0(
    "stage6i__", result$fit_kind, "__", result$study_order
  )
  result$fit_L_f <- 1L
  result$fit_L_s_1 <- 1L
  result$fit_L_s_2 <- 1L
  result$fit_M_f <- "2"
  result$fit_M_s_1 <- "2"
  result$fit_M_s_2 <- "2"
  result$n_cpus <- 1L
  result$planned_ordinary_t1_sweeps <- 400L
  result$planned_annealing_sweeps <- ifelse(
    result$fit_kind == "study_swap_g12", 99L, 0L
  )
  result$truth_available_to_stopping <- FALSE
  result$truth_available_to_selection <- FALSE
  result$cross_dimension_selection <- FALSE
  result$formal_paper_mc_result <- FALSE
  result[, c(
    "experiment_id", "fit_id", "fit_kind", "data_id", "scenario_id",
    "study_order", "oracle_mode", "replicate", "method_id",
    "fit_config_id", "selection_stratum_id", "route_id", "seed_index",
    "fit_seed", "state_seed", "source_fit_id", "fit_L_f",
    "fit_L_s_1", "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2",
    "n_cpus", "planned_annealing_sweeps",
    "planned_ordinary_t1_sweeps", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "continuation_state_used", "cross_dimension_selection",
    "formal_paper_mc_result"
  )]
}

s6i_validate_manifest <- function(manifest) {
  pp_assert(
    nrow(manifest) == 10L && !anyDuplicated(manifest$fit_id) &&
      sum(manifest$fit_kind == "study_swap_g12") == 4L &&
      sum(manifest$fit_kind == "truth_neighborhood") == 6L &&
      all(table(manifest$study_order[manifest$fit_kind ==
                                      "truth_neighborhood"]) == 3L) &&
      sum(manifest$oracle_mode == "exact") == 2L &&
      sum(manifest$oracle_mode == "near") == 4L &&
      identical(manifest$fit_seed[manifest$fit_kind == "study_swap_g12"],
                84611101:84611104) &&
      all(manifest$n_cpus == 1L) &&
      all(manifest$planned_ordinary_t1_sweeps == 400L) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_selection) &&
      !any(manifest$cross_dimension_selection) &&
      !any(manifest$formal_paper_mc_result),
    "Stage 6I manifest validation failed."
  )
  invisible(TRUE)
}

s6i_source_paths <- function() list(
  observation = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2",
    "observation_bundle.rds"
  ),
  truth = file.path(
    stage6g_output_root, "data", "s6g_minimal_L111_M2", "sealed_truth",
    "truth_bundle.rds"
  ),
  stage6g_complete = file.path(stage6g_output_root,
                               "STAGE6G_TRUTH_FREE_COMPLETE.txt"),
  stage6h_complete = file.path(stage6h_output_root,
                               "STAGE6H_WORKFLOW_COMPLETE.txt"),
  stage6h_evaluation = file.path(
    stage6h_output_root, "truth_evaluation_20260913_v1",
    "EVALUATION_COMPLETE.txt"
  ),
  stage6h_manifest = file.path(stage6h_output_root,
                               "ANALYSIS_ENDPOINT_MANIFEST.csv")
)

s6i_verify_package_binding <- function() {
  pp_assert(file.exists(package_binding_path),
            "Stage 6G package source binding is missing.")
  source_root <- pp_env("STAGE6I_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
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

s6i_verify_environment <- function() {
  s6i_assert_authorized()
  source_paths <- s6i_source_paths()
  pp_assert(all(file.exists(unlist(source_paths))),
            "A required Stage 6G/6H artifact is missing.")
  pp_load_package()
  pp_source_runtime(FALSE)
  s6i_verify_package_binding()
  stage6h_manifest <- pp_read_csv(source_paths$stage6h_manifest)
  imported <- stage6h_manifest[
    stage6h_manifest$condition_id == "full_data_6_9" &
      stage6h_manifest$seed_index %in% 1:4, , drop = FALSE
  ]
  pp_assert(nrow(imported) == 4L &&
              all(file.exists(file.path(
                stage6h_output_root, "fits", imported$fit_id, "fit.rds"
              ))), "Stage 6H paired read-only endpoints are incomplete.")
  invisible(TRUE)
}

s6i_template_fit <- function(observation, seed) {
  captured <- c2_capture_conditions(do.call(
    getExportedValue("multiFSYNC", "bayesSYNC_multi"), list(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 5L, anneal = NULL,
      list_hyper = NULL, n_g = 51L, time_g = NULL,
      tol_abs = 0, tol_rel = 0, maxit = 1L, n_cpus = 1L,
      verbose = FALSE, seed = as.integer(seed), bool_scale = FALSE,
      bool_var_spec_prob = FALSE, d_0 = 500L,
      convergence_rule = "parameters", lambda_orth = 0,
      initialization = "random", continuation_state = NULL
    )
  ), phase = "stage6i_state_template")
  pp_assert(is.null(captured$error) && is.list(captured$value),
            paste0("Stage 6I template failed: ",
                   captured$error$message %||% "unknown"))
  captured
}

s6i_oracle_result <- function(observation, state, row) {
  control <- multiFSYNC::g12_stopping_control()
  captured <- c2_capture_conditions(do.call(
    getExportedValue("multiFSYNC", "bayesSYNC_multi"), list(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 5L, anneal = NULL,
      list_hyper = NULL, n_g = 51L, time_g = NULL,
      tol_abs = 1e-3, tol_rel = 1e-5, maxit = 400L,
      n_cpus = 1L, verbose = FALSE,
      seed = as.integer(row$fit_seed[[1L]]), bool_scale = FALSE,
      bool_var_spec_prob = FALSE, d_0 = 500L,
      convergence_rule = "practical", practical_control = control,
      lambda_orth = 0, initialization = "random",
      continuation_state = state
    )
  ), phase = "stage6i_truth_neighborhood_fixed400")
  fit <- captured$value
  objective <- list(eligible = FALSE, checks = logical(),
                    invalid_reasons = "terminal_error")
  if (is.null(captured$error) && is.list(fit)) {
    objective <- r_route_v3_candidate2_objective_status(fit)
    checks <- c(
      objective_eligible = objective$eligible,
      fixed_400 = identical(as.integer(fit$t1_sweeps), 400L),
      no_annealing = identical(as.integer(fit$annealing_sweeps), 0L),
      one_cpu = identical(as.integer(fit$n_cpus_used), 1L),
      ordinary_elbo = identical(as.numeric(fit$lambda_orth), 0)
    )
    if (!all(checks)) {
      captured$error <- list(
        class = "stage6i_validation_error",
        message = paste(names(checks)[!checks], collapse = ";"),
        call = "s6i_oracle_result", trace_summary = "route validation"
      )
      fit <- NULL
      objective$eligible <- FALSE
      objective$invalid_reasons <- names(checks)[!checks]
    }
  }
  list(fit = fit, captured = captured, objective = objective)
}

s6i_prepare <- function(manifest) {
  s6i_assert_authorized()
  pp_assert(nzchar(output_root) && !dir.exists(output_root),
            "Stage 6I output must be a new path.")
  s6i_verify_environment()
  source_paths <- s6i_source_paths()
  observation <- readRDS(source_paths$observation)
  truth <- readRDS(source_paths$truth)
  pp_assert(identical(truth$data_id, "s6g_minimal_L111_M2") &&
              identical(as.integer(truth$data$true_params$L_f), 1L) &&
              all(as.integer(truth$data$true_params$M_f) == 2L),
            "Stage 6I truth identity failed.")
  swapped_observation <- s6i_swap_observation(observation)
  swapped_truth <- s6i_swap_truth_bundle(truth)
  double_swap <- s6i_swap_truth_bundle(swapped_truth)
  pp_assert(
    identical(double_swap$data$true_params$a_true,
              truth$data$true_params$a_true) &&
      identical(double_swap$data$true_params$b_true,
                truth$data$true_params$b_true) &&
      identical(double_swap$data$Y, truth$data$Y),
    "Study swap is not involutive for scientific inputs."
  )
  dir.create(file.path(output_root, "inputs"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  pp_save_rds(swapped_observation,
              file.path(output_root, "inputs", "swapped_observation.rds"))
  pp_save_rds(swapped_truth,
              file.path(output_root, "inputs", "swapped_truth.rds"))
  pp_write_csv(manifest, file.path(output_root, "FIT_MANIFEST.csv"))

  original_template <- s6i_template_fit(observation, 84612901L)
  swapped_template <- s6i_template_fit(swapped_observation, 84612902L)
  maker <- getFromNamespace(".make_reduced_continuation_state", "multiFSYNC")
  templates <- list(
    original = list(fit = original_template$value,
                    state = maker(original_template$value), truth = truth),
    swapped = list(fit = swapped_template$value,
                   state = maker(swapped_template$value), truth = swapped_truth)
  )
  state_root <- file.path(output_root, "oracle_states")
  dir.create(state_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  state_rows <- list()
  oracle_rows <- manifest[manifest$fit_kind == "truth_neighborhood", ,
                          drop = FALSE]
  for (index in seq_len(nrow(oracle_rows))) {
    row <- oracle_rows[index, , drop = FALSE]
    template <- templates[[row$study_order[[1L]]]]
    state <- s6i_map_truth_state(
      template$state, template$truth, template$fit,
      mode = row$oracle_mode[[1L]], seed = row$state_seed[[1L]]
    )
    s6i_validate_state(state, template$state)
    state_path <- file.path(state_root, paste0(row$fit_id[[1L]], ".rds"))
    pp_save_rds(state, state_path)
    diagnostic <- state$diagnostics
    state_rows[[index]] <- data.frame(
      fit_id = row$fit_id[[1L]], study_order = row$study_order[[1L]],
      oracle_mode = row$oracle_mode[[1L]],
      loading_target_cosine = diagnostic$loading_target_cosine,
      projected_mean_rmse_study1 =
        diagnostic$projected_mean_relative_rmse[[1L]],
      projected_mean_rmse_study2 =
        diagnostic$projected_mean_relative_rmse[[2L]],
      projected_shared_function_rmse =
        diagnostic$projected_shared_function_relative_rmse,
      projected_specific_function_rmse_study1 =
        diagnostic$projected_specific_function_relative_rmse[[1L]],
      projected_specific_function_rmse_study2 =
        diagnostic$projected_specific_function_relative_rmse[[2L]],
      state_sha256 = pp_sha256(state_path), stringsAsFactors = FALSE
    )
  }
  pp_write_csv(do.call(rbind, state_rows),
               file.path(output_root, "ORACLE_STATE_DIAGNOSTICS.csv"))
  pp_write_csv(rbind(
    data.frame(study_order = "original", template_seed = 84612901L,
               warnings = nrow(original_template$warnings)),
    data.frame(study_order = "swapped", template_seed = 84612902L,
               warnings = nrow(swapped_template$warnings))
  ), file.path(output_root, "STATE_TEMPLATE_DIAGNOSTICS.csv"))
  bound <- unlist(source_paths, use.names = FALSE)
  pp_write_csv(data.frame(
    role = names(source_paths), path = bound,
    bytes = as.numeric(file.info(bound)$size),
    sha256 = vapply(bound, pp_sha256, character(1L)),
    mutable_by_stage6i = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "SOURCE_ARTIFACT_BINDING.csv"))
  pp_write_lines(c(
    "status=PREPARED", paste0("prepared_utc=", pp_iso_time()),
    "new_scientific_fits=10", "study_swap_g12_fits=4",
    "truth_neighborhood_fits=6", "state_templates=2",
    "read_only_stage6h_paired_endpoints=4",
    "truth_used_for_oracle_state_construction=TRUE",
    "truth_used_for_study_swap_g12_fit=FALSE",
    "model_prior_cavi_elbo_changed=FALSE",
    "new_dimension_criterion_used=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6I_PREPARED.txt"))
  cat("STAGE6I_PREPARED new_fits=10 templates=2 imported=4\n")
  invisible(TRUE)
}

s6i_load_manifest <- function() {
  path <- file.path(output_root, "FIT_MANIFEST.csv")
  pp_assert(file.exists(path), "Stage 6I is not prepared.")
  value <- pp_read_csv(path)
  s6i_validate_manifest(value)
  value
}

s6i_smoke <- function() {
  s6i_verify_environment()
  manifest <- s6i_load_manifest()
  swapped_observation <- readRDS(file.path(
    output_root, "inputs", "swapped_observation.rds"
  ))
  g12_row <- manifest[manifest$fit_kind == "study_swap_g12", ,
                      drop = FALSE][1L, , drop = FALSE]
  g12 <- pp_fit_g12(swapped_observation, g12_row, smoke = TRUE)
  pp_assert(!is.null(g12$fit) && isTRUE(g12$objective$eligible),
            "Stage 6I swapped G12 micro-smoke failed.")

  oracle_row <- manifest[
    manifest$fit_kind == "truth_neighborhood" &
      manifest$study_order == "original" & manifest$oracle_mode == "exact",
    , drop = FALSE
  ][1L, , drop = FALSE]
  state <- readRDS(file.path(
    output_root, "oracle_states", paste0(oracle_row$fit_id[[1L]], ".rds")
  ))
  original_observation <- readRDS(s6i_source_paths()$observation)
  captured <- c2_capture_conditions(do.call(
    getExportedValue("multiFSYNC", "bayesSYNC_multi"), list(
      Y = original_observation$Y, Z = original_observation$Z,
      time_obs = original_observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 5L, anneal = NULL,
      list_hyper = NULL, n_g = 51L, time_g = NULL,
      tol_abs = 0, tol_rel = 0, maxit = 3L, n_cpus = 1L,
      verbose = FALSE, seed = 84612801L, bool_scale = FALSE,
      bool_var_spec_prob = FALSE, d_0 = 500L,
      convergence_rule = "parameters", lambda_orth = 0,
      initialization = "random", continuation_state = state
    )
  ), phase = "stage6i_oracle_micro_smoke")
  pp_assert(is.null(captured$error) && is.list(captured$value) &&
              length(captured$value$ELBO) == 3L &&
              all(is.finite(captured$value$ELBO)),
            "Stage 6I oracle-state micro-smoke failed.")
  cat("STAGE6I_SMOKE_PASS swapped_g12=TRUE oracle_state=TRUE\n")
  invisible(TRUE)
}

s6i_run_fit <- function(fit_id) {
  s6i_verify_environment()
  manifest <- s6i_load_manifest()
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Unknown Stage 6I fit_id.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    pp_assert(xor(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")),
                  file.exists(file.path(fit_dir, "FIT_ERROR.txt"))),
              "Existing Stage 6I fit directory is incomplete.")
    cat("FIT_ALREADY_TERMINAL fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  original_observation <- readRDS(s6i_source_paths()$observation)
  observation <- if (row$study_order[[1L]] == "original") {
    original_observation
  } else readRDS(file.path(output_root, "inputs", "swapped_observation.rds"))
  result <- if (row$fit_kind[[1L]] == "study_swap_g12") {
    pp_fit_g12(observation, row, smoke = FALSE)
  } else {
    state_path <- file.path(output_root, "oracle_states", paste0(fit_id, ".rds"))
    pp_assert(file.exists(state_path), "Oracle state is missing.")
    s6i_oracle_result(observation, readRDS(state_path), row)
  }
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6I_CAUSAL_LOCALIZATION_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$fit_kind <- row$fit_kind[[1L]]
  record$study_order <- row$study_order[[1L]]
  record$oracle_mode <- row$oracle_mode[[1L]]
  record$truth_used_for_fit_stopping_or_selection <-
    isTRUE(row$truth_available_to_fit[[1L]])
  record$truth_used_for_fit <- isTRUE(row$truth_available_to_fit[[1L]])
  record$truth_used_for_stopping <- FALSE
  record$truth_used_for_selection <- FALSE
  record$continuation_used <- isTRUE(row$continuation_state_used[[1L]])
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
  pp_write_csv(data.frame(
    t1_sweep = seq_along(result$fit$ELBO),
    ordinary_elbo = as.numeric(result$fit$ELBO)
  ), file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  if (is.data.frame(result$fit$driver_trace)) {
    pp_write_csv(result$fit$driver_trace,
                 file.path(fit_dir, "DRIVER_TRACE.csv"))
  }
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("fit_kind=", row$fit_kind[[1L]]),
    paste0("study_order=", row$study_order[[1L]]),
    paste0("truth_used_for_fit=", row$truth_available_to_fit[[1L]]),
    "truth_used_for_stopping_or_selection=FALSE",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    "ordinary_T1_sweeps=400", "n_cpus=1",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

s6i_write_status <- function(state, detail = "") {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_new_fits = 10L, complete_new_fits = complete,
    error_new_fits = errors, outer_workers = workers,
    per_fit_n_cpus = 1L, formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

s6i_run_supervisor <- function() {
  s6i_verify_environment()
  manifest <- s6i_load_manifest()
  pp_assert(.Platform$OS.type == "unix", "Supervisor requires Unix.")
  log_root <- file.path(output_root, "supervisor_logs")
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  s6i_write_status("RUNNING", "10 new diagnostic fits")
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
    pp_write_csv(data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      exit_status = as.integer(status), stringsAsFactors = FALSE
    ), file.path(log_root, paste0(key, ".exit.csv")))
    list(fit_id = fit_id, status = as.integer(status))
  }
  results <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = workers,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    s6i_write_status("INCOMPLETE_PROCESS_FAILURE", paste(failed, collapse = ";"))
    pp_abort("Stage 6I child failures retained: ", paste(failed, collapse = ","))
  }
  records <- lapply(manifest$fit_id, function(fit_id) {
    readRDS(file.path(output_root, "fits", fit_id, "terminal_record.rds"))
  })
  terminals <- pp_rbind_fill(lapply(records, pp_terminal_row))
  terminals$fit_kind <- manifest$fit_kind
  terminals$study_order <- manifest$study_order
  terminals$oracle_mode <- manifest$oracle_mode
  terminals$truth_used_for_fit <- manifest$truth_available_to_fit
  pp_write_csv(terminals, file.path(output_root, "ALL_NEW_TERMINALS.csv"))
  pp_assert(nrow(terminals) == 10L &&
              all(terminals$terminal_status == "fixed_400_complete") &&
              all(terminals$objective_eligible),
            "Stage 6I endpoints are not all complete and eligible.")
  pp_write_lines(c(
    "status=FITS_COMPLETE", paste0("completed_utc=", pp_iso_time()),
    "new_fits=10", "objective_eligible=10",
    "scientific_evaluation_not_yet_started=TRUE",
    "new_dimension_criterion_used=FALSE", "formal_paper_mc=FALSE"
  ), file.path(output_root, "STAGE6I_FITS_COMPLETE.txt"))
  s6i_write_status("FITS_COMPLETE", "10/10 complete and eligible")
  cat("STAGE6I_FITS_COMPLETE new=10 eligible=10\n")
  invisible(TRUE)
}

s6i_status <- function() {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  evaluated <- file.exists(file.path(output_root,
                                     "STAGE6I_EVALUATION_COMPLETE.txt"))
  cat("STAGE6I_STATUS complete=", complete, "/10 errors=", errors,
      " evaluated=", evaluated, "\n", sep = "")
}

manifest <- s6i_build_manifest()
s6i_validate_manifest(manifest)

if (action == "check") {
  s6i_assert_authorized()
  cat("STAGE6I_CHECK_PASS new=10 swap_g12=4 oracle=6 imported=4\n")
} else if (action == "preflight") {
  s6i_verify_environment()
  cat("STAGE6I_PREFLIGHT_PASS package_binding=29 sources=6\n")
} else if (action == "prepare") {
  s6i_prepare(manifest)
} else if (action == "smoke") {
  s6i_smoke()
} else if (action == "run-fit") {
  fit_id <- cli$`fit-id` %||% pp_abort("--fit-id is required.")
  s6i_run_fit(fit_id)
} else if (action == "run") {
  s6i_run_supervisor()
} else if (action == "status") {
  s6i_status()
} else {
  pp_abort("Unknown Stage 6I action: ", action)
}
