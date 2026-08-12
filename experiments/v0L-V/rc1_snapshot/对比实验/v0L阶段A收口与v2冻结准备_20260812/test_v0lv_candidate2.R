# Candidate2 unit and regression tests. These tests use only generated micro
# data and synthetic records; formal_experiment is always FALSE.

candidate2_test_script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run candidate2 tests with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/",
                mustWork = TRUE)
}

script_path <- candidate2_test_script_path()
base_dir <- dirname(script_path)
project_root <- normalizePath(file.path(base_dir, "..", ".."),
                              winslash = "/", mustWork = TRUE)
local_library <- normalizePath(file.path(project_root, "_r_test_lib"),
                               winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))
source(file.path(base_dir, "v0lv_candidate2_runtime.R"), local = FALSE)
source(file.path(base_dir, "r_route_v3_candidate2_runner.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_controller.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_evaluation.R"), local = FALSE)

if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
  stop("The current multiFSYNC source must be installed in _r_test_lib.")
}

assert_true <- function(value, message = "assertion failed") {
  if (!isTRUE(value)) stop(message)
  invisible(TRUE)
}

expect_error_message <- function(expression, pattern = NULL) {
  condition <- tryCatch({ force(expression); NULL }, error = identity)
  if (is.null(condition)) stop("Expected an error but none was raised.")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(condition),
                                  fixed = TRUE)) {
    stop("Unexpected error: ", conditionMessage(condition))
  }
  invisible(condition)
}

candidate2_test_record <- function(
    fit_id = "test_fit", data_id = "test_data", method = "R",
    route = "test", fit_seed = 1L, seed_index = 1L,
    initialization_independence_id = "test_independence",
    terminal_status = "max_budget_reached", warnings = c2_empty_warnings(),
    error = list(class = "", message = "", call = "", trace_summary = ""),
    objective_eligible = TRUE, final_elbo = -1,
    formal_experiment = FALSE) {
  now <- c2_iso_time()
  list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID,
    runner_version = "candidate2_test", fit_id = fit_id, data_id = data_id,
    method = method, route = route, fit_seed = fit_seed,
    seed_index = seed_index,
    initialization_independence_id = initialization_independence_id,
    started_at_utc = now, ended_at_utc = now, elapsed_seconds = 0,
    peak_memory_bytes = 0, terminal_status = terminal_status,
    terminal_reason = "test", error = error, warnings = warnings,
    actual_annealing_sweeps = 0L, actual_T1_sweeps = 1L,
    first_strict_convergence_sweep = NA_integer_,
    strict_practical_converged = FALSE, convergence_diagnostics = list(),
    objective_eligible = objective_eligible, objective_checks = logical(),
    objective_invalid_reasons = if (objective_eligible) character() else "test",
    final_elbo = final_elbo,
    hashes = list(input = "test", config = "test", source = "test"),
    formal_experiment = formal_experiment,
    truth_used_for_fit_or_selection = FALSE
  )
}

test_results <- list(); test_index <- 0L
run_test <- function(name, expression) {
  test_index <<- test_index + 1L
  started <- Sys.time()
  condition <- NULL
  tryCatch(force(expression), error = function(error) condition <<- error)
  test_results[[test_index]] <<- data.frame(
    test_id = test_index, test_name = name,
    status = if (is.null(condition)) "PASS" else "FAIL",
    message = if (is.null(condition)) "" else conditionMessage(condition),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    formal_experiment = FALSE, stringsAsFactors = FALSE
  )
  invisible(is.null(condition))
}

run_test("cross_platform_sha256_matches_known_value", {
  path <- tempfile("candidate2_sha256_")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, open = "wb")
  writeBin(charToRaw("abc"), connection)
  close(connection)
  assert_true(identical(
    c2_sha256(path),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
  ))
})

generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
micro_config <- candidate2_generator_config("smoke")
generator_call <- micro_config[setdiff(names(micro_config), c(
  "loading_calibration", "specific_direction_abs_cosine_max",
  "formal_experiment"
))]
generator_call$seed <- 82649001L
micro_data <- do.call(generator, generator_call)
micro_observation <- list(
  bundle_class = "v0lv_observation_only_candidate2", data_id = "micro_test",
  data_seed = 82649001L, formal_experiment = FALSE,
  Y = micro_data$Y, Z = micro_data$Z, time_obs = micro_data$time_obs,
  dimensions = micro_config[c("S", "n_s", "p", "d", "L_f", "L_s",
                              "M_f", "M_s", "K")]
)
micro_truth <- list(
  bundle_class = "v0lv_sealed_truth_candidate2", data_id = "micro_test",
  data_seed = 82649001L, formal_experiment = FALSE, data = micro_data,
  sealed_generation_provenance = list(mean_amp = 0.6)
)

run_test("truth_assertions_match_generator_contract_and_allow_B_overlap", {
  checks <- candidate2_truth_assertions(micro_data, micro_config)
  assert_true(all(checks$passed), paste(
    checks$assertion_id[!checks$passed], collapse = ";"
  ))
  support_a <- which(micro_data$true_params$a_true[, 1L] != 0)
  support_b1 <- which(micro_data$true_params$b_true[[1L]][, 1L] != 0)
  support_b2 <- which(micro_data$true_params$b_true[[2L]][, 1L] != 0)
  assert_true(!length(intersect(support_a, support_b1)))
  assert_true(!length(intersect(support_a, support_b2)))
  assert_true(is.finite(candidate2_abs_cosine(
    micro_data$true_params$b_true[[1L]][, 1L],
    micro_data$true_params$b_true[[2L]][, 1L]
  )))
})

common_fit_args <- list(
  Y = micro_data$Y, Z = micro_data$Z, time_obs = micro_data$time_obs,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L), K = 3L,
  anneal = c(1, 1.2, 3), list_hyper = NULL, n_g = 31L, time_g = NULL,
  tol_abs = 0, tol_rel = 0, maxit = 6L, n_cpus = 1L, verbose = FALSE,
  seed = 82649011L, bool_scale = FALSE, bool_var_spec_prob = FALSE,
  d_0 = 20L, convergence_rule = "parameters", lambda_orth = 0,
  practical_control = NULL, initialization = "random",
  initialization_control = list(grid_size = 81L, perturb_sd = 0.05,
                                rank_tol = 1e-8),
  continuation_state = NULL
)

run_test("public_pre_score_equals_historical_private_driver", {
  public <- do.call(
    getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score"),
    c(common_fit_args, list(pre_score_sweeps = 1L, trace_sweeps = 1L))
  )
  private <- do.call(
    getFromNamespace(".bayesSYNC_multi_driver_diagnostic", "multiFSYNC"),
    c(common_fit_args, list(control = list(
      trace_sweeps = 1L, dense_gate_sweeps = 0L,
      random_scale_calibration = "none", pre_score_sweeps = 1L
    )))
  )
  for (field in c("ELBO", "mu_q_a_hat", "mu_q_b_specific_hat",
                  "list_Zeta_hat", "factor_ppi_shared",
                  "factor_ppi_specific")) {
    assert_true(isTRUE(all.equal(public[[field]], private[[field]],
                                 tolerance = 1e-12)), field)
  }
  assert_true(public$pre_score_interface$pre_score_sweeps == 1L)
})

run_test("warning_is_preserved_in_atomic_terminal_record", {
  captured <- c2_capture_conditions({
    warning("candidate2 warning sentinel", call. = FALSE); 1L
  }, phase = "warning_test")
  assert_true(nrow(captured$warnings) == 1L)
  record <- candidate2_test_record(warnings = captured$warnings)
  directory <- tempfile("candidate2_warning_terminal_")
  written <- c2_write_terminal_bundle(record, NULL, directory)
  stored <- readRDS(written$terminal_path)
  assert_true(stored$warnings$message[[1L]] == "candidate2 warning sentinel")
  assert_true(stored$warnings$phase[[1L]] == "warning_test")
})

error_result <- NULL
run_test("error_path_produces_complete_terminal_row", {
  error_result <<- fit_R_route_v3_candidate2(
    Y = micro_data$Y, Z = micro_data$Z, time_obs = micro_data$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L), K = 3L,
    fit_id = "micro_error_string_id", data_id = "micro_test",
    fit_seed = 82649021L, seed_index = 1L,
    initialization_independence_id = "micro_test__R__82649021",
    mode = "smoke", fit_function = function(...) stop("sentinel fit error")
  )
  assert_true(error_result$terminal_record$terminal_status == "error")
  assert_true(grepl("sentinel fit error",
                    error_result$terminal_record$error$message, fixed = TRUE))
  directory <- tempfile("candidate2_error_terminal_")
  c2_write_terminal_bundle(error_result$terminal_record, NULL, directory)
  assert_true(file.exists(file.path(directory, "TERMINAL_COMPLETE.txt")))
})

run_test("string_fit_id_exact_manifest_join_and_seed_mismatch_rejection", {
  manifest <- data.frame(
    fit_id = c("fit_alpha_001", "fit_beta_X"), data_id = "d",
    r_only_protocol_id = R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID,
    route_id = R_ROUTE_V3_CANDIDATE2_ID, fit_seed = c(11L, 12L),
    seed_index = c(1L, 2L),
    initialization_independence_id = c("iid_alpha", "iid_beta"),
    terminal_record_required = TRUE, actual_racing = FALSE,
    truth_available_to_fit = FALSE, formal_fit_started = FALSE,
    stringsAsFactors = FALSE
  )
  terminals <- manifest[c("fit_id", "data_id", "fit_seed", "seed_index",
                          "initialization_independence_id")]
  joined <- r_route_v3_candidate2_exact_join(manifest, terminals[2:1, ])
  assert_true(identical(joined$fit_id, manifest$fit_id))
  wrong <- terminals; wrong$fit_seed[[2L]] <- 99L
  expect_error_message(
    r_route_v3_candidate2_exact_join(manifest, wrong),
    "do not exactly join"
  )
})

run_test("formal_selection_rejects_seed_identity_mismatch", {
  manifest <- data.frame(
    fit_id = c("formal_fit_A", "formal_fit_B"), data_id = "formal_d",
    r_only_protocol_id = R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID,
    route_id = R_ROUTE_V3_CANDIDATE2_ID, fit_seed = c(101L, 102L),
    seed_index = c(1L, 2L),
    initialization_independence_id = c("formal_iid_A", "formal_iid_B"),
    terminal_record_required = TRUE, actual_racing = FALSE,
    truth_available_to_fit = FALSE, formal_fit_started = FALSE,
    stringsAsFactors = FALSE
  )
  records <- list(
    candidate2_test_record(
      fit_id = "formal_fit_A", data_id = "formal_d", route = R_ROUTE_V3_CANDIDATE2_ID,
      fit_seed = 101L, seed_index = 1L,
      initialization_independence_id = "formal_iid_A",
      final_elbo = -2, formal_experiment = TRUE
    ),
    candidate2_test_record(
      fit_id = "formal_fit_B", data_id = "formal_d", route = R_ROUTE_V3_CANDIDATE2_ID,
      fit_seed = 999L, seed_index = 2L,
      initialization_independence_id = "formal_iid_B",
      final_elbo = -1, formal_experiment = TRUE
    )
  )
  summaries <- do.call(rbind, lapply(records, c2_terminal_summary))
  expect_error_message(select_R_route_v3_candidate2_winner(
    manifest, summaries, expected_starts = 2L, formal_experiment = TRUE
  ), "do not exactly join")
})

run_test("formal_rejects_fake_fit_injection", {
  expect_error_message(fit_R_route_v3_candidate2(
    Y = micro_data$Y, Z = micro_data$Z, time_obs = micro_data$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L), K = 3L,
    fit_id = "formal_fake", data_id = "formal_data", fit_seed = 1L,
    seed_index = 1L, initialization_independence_id = "formal_iid",
    mode = "formal", fit_function = function(...) NULL
  ), "forbids fit_function injection")
})

run_test("formal_source_hash_mismatch_is_rejected", {
  directory <- tempfile("candidate2_bad_binding_"); dir.create(directory)
  utils::write.csv(data.frame(
    status = "frozen", formal_execution_authorized = TRUE
  ), file.path(directory, "BINDING_STATUS.csv"), row.names = FALSE)
  relative <- c2_relative_path(file.path(base_dir, "v0lv_candidate2_runtime.R"),
                               project_root)
  utils::write.csv(data.frame(relative_path = relative,
                              sha256 = paste(rep("0", 64L), collapse = "")),
                   file.path(directory, "FULL_RUNTIME_SOURCE_BINDINGS.csv"),
                   row.names = FALSE)
  expect_error_message(c2_verify_formal_binding(directory, project_root),
                       "Formal source hash mismatch")
})

run_test("endpoint_alignment_compares_all_four_outputs", {
  control <- r_route_v3_candidate2_practical_control()
  winner <- list(practical_control = control, stability_snapshot = list(
    fitted = c(1, 2), rss = c(3, 4), ppi = c(0.1, 0.2),
    factor_ppi = c(0.8, 0.7)
  ))
  same <- winner
  aligned <- r_route_v3_candidate2_endpoint_output_alignment(same, winner)
  assert_true(aligned$endpoint_output_aligned[[1L]])
  changed <- same; changed$stability_snapshot$factor_ppi[[1L]] <- 0.75
  not_aligned <- r_route_v3_candidate2_endpoint_output_alignment(changed, winner)
  assert_true(not_aligned$factor_ppi_max_abs_to_winner[[1L]] > 0.01)
  assert_true(!not_aligned$endpoint_output_aligned[[1L]])
  assert_true(all(c("fitted_nrmse_to_winner", "rss_relative_l2_to_winner",
                    "variable_ppi_max_abs_to_winner",
                    "factor_ppi_max_abs_to_winner") %in% names(not_aligned)))
})

run_test("ever_stable_but_endpoint_unstable_counterexample", {
  diagnostic <- data.frame(
    t1_sweep = 1:7, objective_pass = c(rep(TRUE, 5L), FALSE, FALSE),
    finite = TRUE, monotone = TRUE, long_monotone = TRUE,
    fitted_nrmse = 0, rss_rel = 0, ppi_quantile_abs = 0,
    factor_ppi_max_abs = 0, long_fitted_nrmse = 0, long_rss_rel = 0,
    long_ppi_quantile_abs = 0, long_factor_ppi_max_abs = 0
  )
  stability <- r_route_v3_candidate2_output_stability(list(
    practical_diagnostics = diagnostic,
    practical_control = r_route_v3_candidate2_practical_control()
  ))
  assert_true(stability$ever_reached_5_consecutive)
  assert_true(stability$first_reached_sweep == 5L)
  assert_true(stability$maximum_streak == 5L)
  assert_true(!stability$endpoint_pass)
  assert_true(stability$endpoint_streak == 0L)
  assert_true(!stability$endpoint_stable_5_consecutive)
})

run_test("truth_cannot_be_unsealed_before_selection_authorization", {
  root <- tempfile("candidate2_truth_lock_"); dir.create(root)
  truth_path <- file.path(root, "truth.rds"); saveRDS(micro_truth, truth_path)
  selection <- file.path(root, "selection"); dir.create(selection)
  authorization <- file.path(root, "authorization")
  expect_error_message(candidate2_unseal_truth(
    truth_path, selection, authorization
  ), "sealed until selection freeze")
})

run_test("terminal_resume_is_idempotent_and_non_overwriting", {
  directory <- tempfile("candidate2_resume_")
  record <- candidate2_test_record(fit_id = "resume_string_fit")
  first <- c2_write_terminal_bundle(record, list(value = 1L), directory)
  terminal_hash <- c2_sha256(first$terminal_path)
  second <- c2_write_terminal_bundle(record, list(value = 999L), directory)
  assert_true(second$reused)
  assert_true(c2_sha256(first$terminal_path) == terminal_hash)
  assert_true(readRDS(file.path(directory, "fit.rds"))$value == 1L)
})

run_test("audit_main_error_is_terminal_and_never_rerun", {
  directory <- tempfile("candidate2_audit_error_")
  audit <- candidate2_run_R_audit_horizon(
    observation = micro_observation, source_fit = NULL,
    source_record = error_result$terminal_record,
    target_cumulative_T1 = 6L, audit_id = "micro_audit_error",
    output_dir = directory, mode = "smoke",
    fit_function = function(...) stop("must not be called"),
    state_function = function(...) stop("must not be called")
  )
  assert_true(audit$record$terminal_status ==
                "audit_unavailable_due_to_main_error")
  assert_true(file.exists(file.path(directory, "TERMINAL_COMPLETE.txt")))
})

run_test("audit_horizon_freeze_retains_unavailable_error_row", {
  eligible_id <- "audit_eligible"
  unavailable_id <- "audit_unavailable"
  eligible <- candidate2_test_record(
    fit_id = eligible_id, data_id = "audit_data", method = "R_audit",
    route = "fixed_horizon_continuation", fit_seed = 201L, seed_index = 1L,
    initialization_independence_id = "audit_iid_1",
    terminal_status = "audit_horizon_reached", final_elbo = -10
  )
  unavailable <- candidate2_test_record(
    fit_id = unavailable_id, data_id = "audit_data", method = "R_audit",
    route = "fixed_horizon_continuation", fit_seed = 202L, seed_index = 2L,
    initialization_independence_id = "audit_iid_2",
    terminal_status = "audit_unavailable_due_to_main_error",
    objective_eligible = FALSE, final_elbo = NA_real_,
    error = list(class = "main_fit_error", message = "unavailable",
                 call = "", trace_summary = "")
  )
  manifest <- data.frame(
    audit_id = c(eligible_id, unavailable_id), data_id = "audit_data",
    seed_index = c(1L, 2L), original_fit_seed = c(201L, 202L),
    initialization_independence_id = c("audit_iid_1", "audit_iid_2"),
    target_cumulative_T1 = 200L, stringsAsFactors = FALSE
  )
  fit <- list(
    practical_control = r_route_v3_candidate2_practical_control(),
    stability_snapshot = list(
      fitted = c(1, 2), rss = c(3, 4), ppi = c(0.1, 0.2),
      factor_ppi = c(0.8, 0.7)
    )
  )
  fits <- setNames(list(fit, NULL), c(eligible_id, unavailable_id))
  output <- tempfile("candidate2_audit_freeze_")
  frozen <- candidate2_freeze_audit_horizon_truth_free(
    manifest, list(eligible, unavailable), fits, 200L, output,
    expected_starts = 2L, formal_experiment = FALSE
  )
  assert_true(nrow(frozen$endpoints) == 2L)
  assert_true(frozen$winner$fit_id[[1L]] == eligible_id)
  assert_true(any(frozen$endpoints$terminal_status ==
                    "audit_unavailable_due_to_main_error"))
})

run_test("pooled_matches_all_three_oracles_and_uses_ppi_auxiliarily", {
  truth <- micro_data$true_params
  pooled_fit <- list(
    B_hat = cbind(truth$b_true[[1L]][, 1L], truth$a_true[, 1L],
                  truth$b_true[[2L]][, 1L]),
    factor_ppi = c(0.1, 0.4, 0.9), list_h_hat = list()
  )
  global <- 0L
  for (study in seq_len(2L)) for (subject in seq_len(truth$n_s[[study]])) {
    global <- global + 1L
    zeros <- rep(0, length(micro_data$time_obs[[study]][[subject]]))
    pooled_fit$list_h_hat[[global]] <- list(
      if (study == 1L) candidate2_truth_path(
        micro_data, study, "B10", subject
      ) else zeros,
      candidate2_truth_path(micro_data, study, "A0", subject),
      if (study == 2L) candidate2_truth_path(
        micro_data, study, "B20", subject
      ) else zeros
    )
  }
  matches <- candidate2_match_pooled_directions(pooled_fit, micro_data)
  assert_true(nrow(matches) == 3L && all(matches$evaluated_in_primary))
  assert_true(identical(matches$pooled_factor, c(2L, 1L, 3L)))
  assert_true(!matches$auxiliary_factor_retained[matches$oracle_id == "A0"])
  metric <- candidate2_pooled_direction_metrics(pooled_fit, micro_data, matches)
  assert_true(all(is.na(metric$misplaced_status)))
  assert_true(all(is.na(metric$block_purity)))
  assert_true(is.na(metric$inactive_study_leakage_rms[metric$oracle_id == "A0"]))
  assert_true(all(metric$complete_contribution_nrmse < 1e-12))
})

run_test("statistics_use_median_gate_and_descriptive_worsening_count", {
  paired <- data.frame(
    data_id = paste0("d", 1:4), comparison = "R_minus_A",
    metric = "observed_nrmse", difference = c(-0.2, -0.1, 0.1, 0.2)
  )
  summary <- candidate2_summarize_paired_differences(paired)
  assert_true(summary$success_gate_name[[1L]] == "中位无恶化/优势门槛")
  assert_true(summary$median_no_worsening_or_advantage_gate[[1L]])
  assert_true(summary$worsening_count_is_descriptive_only[[1L]])
  assert_true(!summary$invented_noninferiority_margin_used[[1L]])
  assert_true(summary$bootstrap_replicates[[1L]] == 9999L)
})

run_test("candidate2_manifest_counts_flags_and_seed_collisions", {
  manifest_dir <- file.path(
    base_dir, "protocol_manifests_20260812_v3_candidate2"
  )
  fit <- utils::read.csv(file.path(
    manifest_dir, "V0LV_FIT_MANIFEST_CANDIDATE2.csv"
  ), stringsAsFactors = FALSE)
  audit <- utils::read.csv(file.path(
    manifest_dir, "V0LV_AUDIT_MANIFEST_CANDIDATE2.csv"
  ), stringsAsFactors = FALSE)
  status <- utils::read.csv(file.path(manifest_dir, "BINDING_STATUS.csv"),
                            stringsAsFactors = FALSE)
  collision <- utils::read.csv(file.path(
    manifest_dir, "SEED_COLLISION_AUDIT.csv"
  ), stringsAsFactors = FALSE)
  assert_true(nrow(fit) == 210L && nrow(audit) == 72L)
  counts <- table(fit$path_id)
  assert_true(counts[["B0"]] == 40L && counts[["A"]] == 40L &&
                counts[["R"]] == 120L &&
                counts[["pooled_bayesSYNC"]] == 10L)
  assert_true(all(!fit$formal_fit_started) &&
                all(!audit$formal_continuation_started))
  assert_true(status$status[[1L]] == "candidate" &&
                !status$formal_execution_authorized[[1L]])
  assert_true(all(collision$intentional_or_unique))
})

results <- do.call(rbind, test_results)
requested_output <- commandArgs(trailingOnly = TRUE)
requested_output <- requested_output[grepl("^--output-dir=", requested_output)]
output_dir <- if (length(requested_output)) {
  sub("^--output-dir=", "", requested_output[[1L]])
} else file.path(base_dir, "candidate2_test_results_20260812_v1")
if (dir.exists(output_dir)) stop("Test output directory already exists: ", output_dir)
dir.create(output_dir, recursive = TRUE)
results_write <- c2_atomic_write_csv(results, file.path(output_dir, "TEST_RESULTS.csv"))
c2_atomic_write_lines(c(
  "V0LV_CANDIDATE2_TEST_RUN_COMPLETE",
  paste0("tests=", nrow(results)), paste0("passed=", sum(results$status == "PASS")),
  paste0("failed=", sum(results$status == "FAIL")),
  "formal_experiment=FALSE", "formal_data_generated=FALSE",
  "formal_fits_started=FALSE", "formal_continuations_started=FALSE",
  paste0("test_results_sha256=", results_write$sha256),
  capture.output(sessionInfo())
), file.path(output_dir, "TEST_RUN_COMPLETE.txt"))
print(results[, c("test_id", "test_name", "status", "message")], row.names = FALSE)
if (any(results$status == "FAIL")) stop("Candidate2 regression tests failed.")
