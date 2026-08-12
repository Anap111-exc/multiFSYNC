# Complete non-formal candidate2 micro chain. This script never reads or
# modifies the formal manifest start fields and never uses formal data seeds.

candidate2_smoke_script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run this smoke with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/",
                mustWork = TRUE)
}

script_path <- candidate2_smoke_script_path()
base_dir <- dirname(script_path)
project_root <- normalizePath(file.path(base_dir, "..", ".."),
                              winslash = "/", mustWork = TRUE)
.libPaths(c(normalizePath(file.path(project_root, "_r_test_lib"),
                          winslash = "/", mustWork = TRUE), .libPaths()))
source(file.path(base_dir, "v0lv_candidate2_runtime.R"), local = FALSE)
source(file.path(base_dir, "r_route_v3_candidate2_runner.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_controller.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_evaluation.R"), local = FALSE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[grepl("^--output-dir=", arguments)]
output_dir <- if (length(output_argument)) {
  sub("^--output-dir=", "", output_argument[[1L]])
} else file.path(base_dir, "candidate2_e2e_smoke_20260812_v1")
if (dir.exists(output_dir)) stop("Smoke output directory already exists.")
dir.create(output_dir, recursive = TRUE)

status_rows <- list(); status_index <- 0L
record_status <- function(stage, status, detail = "") {
  status_index <<- status_index + 1L
  status_rows[[status_index]] <<- data.frame(
    stage_index = status_index, stage = stage, status = status,
    detail = detail, formal_experiment = FALSE,
    observed_at_utc = c2_iso_time(), stringsAsFactors = FALSE
  )
}

smoke_error <- NULL
tryCatch({
  data_dir <- file.path(output_dir, "data", "candidate2_micro")
  generated <- candidate2_generate_and_seal(
    data_id = "candidate2_micro", data_seed = 82649101L,
    output_dir = data_dir, mode = "smoke", project_root = project_root
  )
  observation <- generated$observation
  observation_path <- file.path(data_dir, "observation_bundle.rds")
  hashes <- list(
    input = c2_sha256(observation_path),
    config = c2_function_signature_sha256(candidate2_generator_config),
    source = c2_sha256(file.path(base_dir, "r_route_v3_candidate2_runner.R"))
  )
  record_status("generate_and_seal", "PASS",
                "micro data sealed; observation-only bundle emitted")

  truth_guard <- tryCatch({
    candidate2_unseal_truth(
      generated$sealed_path, file.path(output_dir, "selection", "all"),
      file.path(output_dir, "unseal_authorization")
    ); FALSE
  }, error = function(error) TRUE)
  if (!truth_guard) stop("Truth was readable before selection authorization.")
  record_status("preselection_truth_lock", "PASS")

  fit_manifest_rows <- list(); fit_manifest_index <- 0L
  add_manifest <- function(fit_id, path_id, seed, seed_index,
                           initialization_id) {
    fit_manifest_index <<- fit_manifest_index + 1L
    fit_manifest_rows[[fit_manifest_index]] <<- data.frame(
      manifest_schema_id = "SMOKE_NOT_FORMAL",
      total_protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
      r_only_protocol_id = if (path_id == "R") {
        R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID
      } else "not_applicable",
      fit_id = fit_id, data_id = observation$data_id, path_id = path_id,
      route_id = if (path_id == "R") {
        r_route_v3_candidate2_spec("smoke")$route_id
      } else path_id,
      fit_seed = seed, seed_index = seed_index,
      initialization_independence_id = initialization_id,
      terminal_record_required = TRUE, actual_racing = FALSE,
      truth_available_to_fit = FALSE, formal_fit_started = FALSE,
      stringsAsFactors = FALSE
    )
  }

  fit_objects <- list(); terminal_records <- list()
  for (method in c("B0", "A")) {
    fit_id <- paste0("candidate2_micro_", method, "_01")
    seed <- if (method == "B0") 82649111L else 82649112L
    independence <- paste(observation$data_id, method, seed, sep = "__")
    result <- candidate2_run_multi_main(
      observation, method = method, fit_id = fit_id, fit_seed = seed,
      seed_index = 1L, initialization_independence_id = independence,
      mode = "smoke", project_root = project_root, hashes = hashes
    )
    c2_write_terminal_bundle(result$terminal_record, result$fit,
                             file.path(output_dir, "fits", fit_id))
    fit_objects[[fit_id]] <- result$fit
    terminal_records[[fit_id]] <- result$terminal_record
    add_manifest(fit_id, method, seed, 1L, independence)
  }
  record_status("B0_A_fit_and_terminal", "PASS")

  r_ids <- character()
  for (seed_index in seq_len(2L)) {
    fit_id <- sprintf("candidate2_micro_R_%02d", seed_index)
    seed <- 82649120L + seed_index
    independence <- paste(observation$data_id, "R", seed, sep = "__")
    result <- fit_R_route_v3_candidate2(
      Y = observation$Y, Z = observation$Z, time_obs = observation$time_obs,
      L_f = observation$dimensions$L_f, L_s = observation$dimensions$L_s,
      M_f = observation$dimensions$M_f, M_s = observation$dimensions$M_s,
      K = observation$dimensions$K, fit_id = fit_id,
      data_id = observation$data_id, fit_seed = seed, seed_index = seed_index,
      initialization_independence_id = independence, mode = "smoke",
      project_root = project_root, hashes = hashes
    )
    c2_write_terminal_bundle(result$terminal_record, result$fit,
                             file.path(output_dir, "fits", fit_id))
    fit_objects[[fit_id]] <- result$fit
    terminal_records[[fit_id]] <- result$terminal_record
    add_manifest(fit_id, "R", seed, seed_index, independence)
    r_ids <- c(r_ids, fit_id)
  }
  record_status("R_public_pre_score_multistart_fit_and_terminal", "PASS")

  pooled_id <- "candidate2_micro_pooled_bayesSYNC_01"
  pooled_seed <- 82649131L
  pooled_independence <- paste(
    observation$data_id, "pooled_bayesSYNC", pooled_seed, sep = "__"
  )
  pooled <- candidate2_run_pooled_main(
    observation, fit_id = pooled_id, fit_seed = pooled_seed,
    seed_index = 1L, initialization_independence_id = pooled_independence,
    mode = "smoke", project_root = project_root, hashes = hashes
  )
  c2_write_terminal_bundle(pooled$terminal_record, pooled$fit,
                           file.path(output_dir, "fits", pooled_id))
  fit_objects[[pooled_id]] <- pooled$fit
  terminal_records[[pooled_id]] <- pooled$terminal_record
  add_manifest(pooled_id, "pooled_bayesSYNC", pooled_seed, 1L,
               pooled_independence)
  record_status("pooled_fit_once_without_study_labels", "PASS")

  warning_capture <- c2_capture_conditions({
    warning("candidate2 e2e warning sentinel", call. = FALSE); 1L
  }, phase = "e2e_warning_terminal")
  warning_record <- list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
    runner_version = "candidate2_e2e_smoke", fit_id = "diagnostic_warning",
    data_id = observation$data_id, method = "diagnostic", route = "warning",
    fit_seed = 82649141L, seed_index = 1L,
    initialization_independence_id = "diagnostic_warning_iid",
    started_at_utc = warning_capture$started_at_utc,
    ended_at_utc = warning_capture$ended_at_utc,
    elapsed_seconds = warning_capture$elapsed_seconds,
    peak_memory_bytes = warning_capture$peak_memory_bytes,
    terminal_status = "max_budget_reached", terminal_reason = "warning_smoke",
    error = list(class = "", message = "", call = "", trace_summary = ""),
    warnings = warning_capture$warnings, actual_annealing_sweeps = 0L,
    actual_T1_sweeps = 1L, first_strict_convergence_sweep = NA_integer_,
    strict_practical_converged = FALSE, convergence_diagnostics = list(),
    objective_eligible = TRUE, objective_checks = logical(),
    objective_invalid_reasons = character(), final_elbo = -1,
    hashes = hashes, formal_experiment = FALSE,
    truth_used_for_fit_or_selection = FALSE
  )
  c2_write_terminal_bundle(warning_record, NULL,
                           file.path(output_dir, "diagnostics", "warning"))
  error_result <- fit_R_route_v3_candidate2(
    Y = observation$Y, Z = observation$Z, time_obs = observation$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L), K = 3L,
    fit_id = "diagnostic_error", data_id = observation$data_id,
    fit_seed = 82649142L, seed_index = 1L,
    initialization_independence_id = "diagnostic_error_iid", mode = "smoke",
    fit_function = function(...) stop("candidate2 e2e error sentinel"),
    hashes = hashes
  )
  c2_write_terminal_bundle(error_result$terminal_record, NULL,
                           file.path(output_dir, "diagnostics", "error"))
  if (nrow(warning_record$warnings) != 1L ||
      error_result$terminal_record$terminal_status != "error") {
    stop("Warning/error terminal smoke failed.")
  }
  record_status("warning_and_error_terminal_records", "PASS")

  fit_manifest <- do.call(rbind, fit_manifest_rows)
  c2_atomic_write_csv(fit_manifest,
                      file.path(output_dir, "SMOKE_FIT_MANIFEST.csv"))
  selection_root <- file.path(output_dir, "selection")
  component_dirs <- list(B0 = file.path(selection_root, "B0"),
                         A = file.path(selection_root, "A"),
                         R = file.path(selection_root, "R"))
  for (method in c("B0", "A")) {
    rows <- fit_manifest[fit_manifest$path_id == method, , drop = FALSE]
    candidate2_freeze_generic_truth_free_selection(
      fit_manifest, terminal_records[rows$fit_id], method,
      component_dirs[[method]], expected_starts = 1L,
      formal_experiment = FALSE
    )
  }
  r_manifest <- fit_manifest[fit_manifest$path_id == "R", , drop = FALSE]
  r_selection <- candidate2_freeze_truth_free_selection(
    r_manifest, terminal_records[r_manifest$fit_id], component_dirs$R,
    expected_starts = 2L, formal_experiment = FALSE
  )
  candidate2_freeze_combined_truth_free_selection(
    component_dirs, pooled$terminal_record, file.path(selection_root, "all"),
    formal_experiment = FALSE
  )
  record_status("truth_free_all_method_selection_and_hash_freeze", "PASS")

  r_summaries <- do.call(rbind, lapply(
    terminal_records[r_manifest$fit_id], c2_terminal_summary
  ))
  r_traces <- lapply(r_manifest$fit_id, function(id) {
    fit <- fit_objects[[id]]
    data.frame(t1_sweep = seq_along(fit$ELBO), elbo = as.numeric(fit$ELBO))
  }); names(r_traces) <- r_manifest$fit_id
  replay <- r_route_v3_candidate2_offline_racing_replay(
    r_summaries, r_traces, r_selection$winner,
    screen_t1 = 100L, minimum_keep = 2L, elbo_margin = 500
  )
  racing_dir <- file.path(selection_root, "offline_racing")
  dir.create(racing_dir, recursive = TRUE)
  c2_atomic_write_csv(replay$per_fit, file.path(racing_dir, "PER_FIT.csv"))
  c2_atomic_write_csv(replay$per_data, file.path(racing_dir, "PER_DATA.csv"))
  if (!all(replay$per_fit$completed_before_screen) ||
      !all(replay$per_fit$retained)) stop("Early completed candidates not retained.")
  record_status("offline_racing_real_trajectory_replay", "PASS")

  r_winner_id <- r_selection$winner$fit_id[[1L]]
  r_winner_fit <- fit_objects[[r_winner_id]]
  r_winner_record <- terminal_records[[r_winner_id]]
  audit_target <- as.integer(r_winner_record$actual_T1_sweeps + 2L)
  audit_id <- paste0(r_winner_id, "_micro_to_", audit_target)
  audit_dir <- file.path(output_dir, "audit", audit_id)
  audit <- candidate2_run_R_audit_horizon(
    observation, r_winner_fit, r_winner_record,
    target_cumulative_T1 = audit_target, audit_id = audit_id,
    output_dir = audit_dir, mode = "smoke", project_root = project_root,
    hashes = hashes
  )
  if (audit$record$terminal_status != "audit_horizon_reached") {
    stop("Micro fixed-horizon continuation did not reach target.")
  }
  audit_fit <- readRDS(file.path(audit_dir, "fit.rds"))
  audit_manifest <- data.frame(
    audit_id = audit_id, data_id = observation$data_id,
    seed_index = r_winner_record$seed_index,
    original_fit_seed = r_winner_record$fit_seed,
    initialization_independence_id =
      r_winner_record$initialization_independence_id,
    target_cumulative_T1 = audit_target, stringsAsFactors = FALSE
  )
  candidate2_freeze_audit_horizon_truth_free(
    audit_manifest, list(audit$record), setNames(list(audit_fit), audit_id),
    audit_target, file.path(output_dir, "audit_selection"),
    expected_starts = 1L, formal_experiment = FALSE
  )
  record_status("micro_fixed_horizon_continuation_and_truth_free_freeze", "PASS")

  candidate2_authorize_unseal(
    file.path(selection_root, "all"), generated$sealed_path,
    file.path(output_dir, "unseal_authorization"),
    formal_experiment = FALSE, operator_authorized = FALSE
  )
  truth <- candidate2_unseal_truth(
    generated$sealed_path, file.path(selection_root, "all"),
    file.path(output_dir, "unseal_authorization")
  )
  record_status("authorized_truth_unseal", "PASS")

  selected_ids <- c(
    utils::read.csv(file.path(component_dirs$B0, "TRUTH_FREE_SELECTION.csv"),
                    stringsAsFactors = FALSE)$fit_id[[1L]],
    utils::read.csv(file.path(component_dirs$A, "TRUTH_FREE_SELECTION.csv"),
                    stringsAsFactors = FALSE)$fit_id[[1L]],
    r_winner_id
  )
  names(selected_ids) <- c("B0", "A", "R")
  multi_evaluations <- lapply(selected_ids, function(id) {
    candidate2_evaluate_multi(
      fit_objects[[id]], truth, project_root, mode = "smoke"
    )
  })
  pooled_evaluation <- candidate2_evaluate_pooled(
    pooled$fit, truth, project_root, mode = "smoke"
  )
  evaluation_dir <- file.path(output_dir, "evaluation")
  dir.create(evaluation_dir, recursive = TRUE)
  c2_atomic_save_rds(multi_evaluations,
                     file.path(evaluation_dir, "MULTIFSYNC_EVALUATIONS.rds"))
  c2_atomic_save_rds(pooled_evaluation,
                     file.path(evaluation_dir, "POOLED_EVALUATION.rds"))
  c2_atomic_write_csv(pooled_evaluation$direction_matches,
                      file.path(evaluation_dir, "POOLED_DIRECTION_MATCHES.csv"))
  c2_atomic_write_csv(pooled_evaluation$direction_metrics,
                      file.path(evaluation_dir, "POOLED_DIRECTION_METRICS.csv"))
  if (nrow(pooled_evaluation$direction_matches) != 3L ||
      !all(pooled_evaluation$direction_matches$evaluated_in_primary)) {
    stop("Pooled evaluation did not retain all three oracle candidates.")
  }
  record_status("multiFSYNC_and_pooled_evaluation", "PASS")

  metric_field <- "dense_signal_nrmse_ppi_selected"
  paired <- data.frame(
    data_id = observation$data_id, comparison = "R_minus_B0",
    metric = metric_field,
    difference = as.numeric(
      multi_evaluations$R$summary[[metric_field]] -
        multi_evaluations$B0$summary[[metric_field]]
    ), stringsAsFactors = FALSE
  )
  summary <- candidate2_summarize_paired_differences(paired)
  c2_atomic_write_csv(paired, file.path(evaluation_dir, "PAIRED_METRICS.csv"))
  c2_atomic_write_csv(summary,
                      file.path(evaluation_dir, "SUCCESS_GATE_SUMMARY.csv"))
  record_status("cluster_bootstrap_and_success_gate_summary", "PASS")
}, error = function(error) {
  smoke_error <<- error
  record_status("terminal_smoke_status", "FAIL", conditionMessage(error))
})

status <- do.call(rbind, status_rows)
status_write <- c2_atomic_write_csv(status, file.path(output_dir, "SMOKE_STATUS.csv"))
c2_atomic_write_lines(c(
  "V0LV_CANDIDATE2_END_TO_END_SMOKE_COMPLETE",
  paste0("status=", if (is.null(smoke_error)) "PASS" else "FAIL"),
  "formal_experiment=FALSE", "formal_data_generated=FALSE",
  "formal_fits_started=FALSE", "formal_continuations_started=FALSE",
  "frozen_bundle_created=FALSE",
  paste0("status_sha256=", status_write$sha256),
  capture.output(sessionInfo())
), file.path(output_dir, "SMOKE_COMPLETE.txt"))
print(status, row.names = FALSE)
if (!is.null(smoke_error)) stop(conditionMessage(smoke_error))
