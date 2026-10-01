#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(
  script_root, "..", "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)

PP_ROOT <- script_root
PP_EXPERIMENT_ID <- "JAOUA_SPLINETS_PAIRED_CONFIRMATION_V1_20260909"
PP_MODEL_PARENT_COMMIT <- "5ec35267def4ccbc1a71c6811a9b06acccb19e9c"
PP_DATA_IDS <- c(
  "jspc_01__legacy", "jspc_01__jaoua",
  "jspc_02__legacy", "jspc_02__jaoua",
  "jspc_03__legacy", "jspc_03__jaoua"
)
PP_FIT_CONFIG_IDS <- "correct_L1_M2"
PP_EXPECTED_FITS <- 72L
PP_EXPECTED_WINNERS <- 6L
PP_WORKERS <- as.integer(pp_env("JSPC_OUTER_WORKERS", "7"))
pp_assert(length(PP_WORKERS) == 1L && is.finite(PP_WORKERS) &&
            PP_WORKERS >= 1L && PP_WORKERS <= 8L,
          "JSPC_OUTER_WORKERS must be an integer from 1 through 8.")

parent_make_terminal_record <- pp_make_terminal_record
parent_terminal_row <- pp_terminal_row

pp_validate_data_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "data_id", "logical_data_id", "logical_data_index",
    "basis_arm", "truth_basis_method", "data_seed", "truth_basis_seed",
    "scenario_id", "S", "n_s", "p", "d", "truth_L_f", "truth_L_s",
    "truth_M_f", "truth_M_s", "K", "n_obs_min", "n_obs_max",
    "sigma_eps", "mean_amp", "paired_nonbasis_rng_required",
    "new_unseen_at_registration", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "truth_unseal_authorized", "development_only", "formal_v0lv_result",
    "formal_paper_mc_result"
  )
  pp_assert(is.data.frame(manifest) && nrow(manifest) == 6L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              identical(as.character(manifest$data_id), PP_DATA_IDS) &&
              !anyDuplicated(manifest$data_id) &&
              all(table(manifest$logical_data_id) == 2L) &&
              all(table(manifest$data_seed) == 2L) &&
              all(table(manifest$truth_basis_seed) == 2L) &&
              setequal(manifest$truth_basis_method,
                       c("legacy_dense_gram", "jaoua_splinets")) &&
              all(manifest$paired_nonbasis_rng_required) &&
              all(manifest$new_unseen_at_registration) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              !any(manifest$truth_unseal_authorized) &&
              all(manifest$development_only) &&
              !any(manifest$formal_v0lv_result) &&
              !any(manifest$formal_paper_mc_result),
            "Paired data manifest contract failed.")
  pp_assert(all(manifest$S == 2L) && all(manifest$n_s == "30;30") &&
              all(manifest$p == 500L) && all(manifest$d == 0L) &&
              all(manifest$truth_L_f == 1L) &&
              all(manifest$truth_L_s == "1;1") &&
              all(manifest$truth_M_f == "2") &&
              all(manifest$truth_M_s == "2;2") &&
              all(manifest$K == 5L) && all(manifest$n_obs_min == 6L) &&
              all(manifest$n_obs_max == 9L) &&
              all(manifest$sigma_eps == 0.3) &&
              all(manifest$mean_amp == 0.6),
            "Paired data dimensions or generation parameters changed.")
  invisible(TRUE)
}

pp_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "logical_data_id",
    "logical_data_index", "basis_arm", "truth_basis_method", "scenario_id",
    "method_id", "fit_config_id", "selection_stratum_id",
    "paired_start_id", "seed_index", "fit_seed", "fit_L_f",
    "fit_L_s_1", "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2",
    "route_id", "initialization", "function_initialization",
    "pre_score_sweeps", "anneal", "planned_annealing_sweeps",
    "planned_ordinary_t1_sweeps", "n_cpus", "technical_gate",
    "objective_eligibility_required", "truth_free_winner_candidate",
    "cross_arm_elbo_comparison_forbidden", "actual_racing", "continuation",
    "automatic_800", "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_selection", "development_only", "formal_v0lv_result",
    "formal_paper_mc_result"
  )
  pp_assert(is.data.frame(manifest) && nrow(manifest) == 72L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              !anyDuplicated(manifest$fit_id) &&
              setequal(manifest$data_id, PP_DATA_IDS) &&
              all(table(manifest$selection_stratum_id) == 12L) &&
              length(unique(manifest$paired_start_id)) == 36L &&
              all(table(manifest$paired_start_id) == 2L),
            "Paired fit-manifest identity failed.")
  pairing_ok <- vapply(split(manifest, manifest$paired_start_id), function(x) {
    nrow(x) == 2L && setequal(x$basis_arm, c("legacy", "jaoua")) &&
      length(unique(x$fit_seed)) == 1L &&
      length(unique(x$logical_data_id)) == 1L
  }, logical(1L))
  pp_assert(all(pairing_ok) && sum(manifest$technical_gate) == 8L &&
              all(manifest$method_id == "G12") &&
              all(manifest$fit_config_id == "correct_L1_M2") &&
              all(manifest$fit_L_f == 1L) &&
              all(manifest$fit_L_s_1 == 1L) &&
              all(manifest$fit_L_s_2 == 1L) &&
              all(manifest$fit_M_f == "2") &&
              all(manifest$fit_M_s_1 == "2") &&
              all(manifest$fit_M_s_2 == "2") &&
              all(manifest$route_id ==
                    "random_g12_prescore1_jaoua_fixed400") &&
              all(manifest$initialization == "random") &&
              all(manifest$function_initialization == "gram_unit_energy") &&
              all(manifest$pre_score_sweeps == 1L) &&
              all(manifest$anneal == "1;1.9;100") &&
              all(manifest$planned_annealing_sweeps == 99L) &&
              all(manifest$planned_ordinary_t1_sweeps == 400L) &&
              all(manifest$n_cpus == 1L) &&
              all(manifest$objective_eligibility_required) &&
              all(manifest$truth_free_winner_candidate) &&
              all(manifest$cross_arm_elbo_comparison_forbidden) &&
              !any(manifest$actual_racing) && !any(manifest$continuation) &&
              !any(manifest$automatic_800) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              all(manifest$development_only) &&
              !any(manifest$formal_v0lv_result) &&
              !any(manifest$formal_paper_mc_result),
            "Paired fit route, budget, pairing or truth isolation changed.")
  invisible(TRUE)
}

pp_make_terminal_record <- function(result, fit_row, smoke = FALSE) {
  value <- parent_make_terminal_record(result, fit_row, smoke)
  value$terminal_schema <- "JAOUA_SPLINETS_PAIRED_CONFIRMATION_TERMINAL_V1"
  value$logical_data_id <- fit_row$logical_data_id[[1L]]
  value$basis_arm <- fit_row$basis_arm[[1L]]
  value$truth_basis_method <- fit_row$truth_basis_method[[1L]]
  value$paired_start_id <- fit_row$paired_start_id[[1L]]
  value
}

pp_terminal_row <- function(record) {
  value <- parent_terminal_row(record)
  value$logical_data_id <- record$logical_data_id
  value$basis_arm <- record$basis_arm
  value$truth_basis_method <- record$truth_basis_method
  value$paired_start_id <- record$paired_start_id
  value
}

pc_verify_environment <- function() {
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  pp_assert(requireNamespace("Splinets", quietly = TRUE),
            "Splinets is required for the Jaoua generator arm.")
  binding <- pp_verify_source_binding()
  pp_load_package()
  pp_source_runtime(FALSE)
  manifests <- pp_load_manifests()
  qc <- pp_read_csv(file.path(PP_ROOT, "MANIFEST_QC.csv"))
  pp_assert(nrow(qc) == 14L && all(qc$passed),
            "Manifest QC is not 14/14 PASS.")
  invisible(list(manifests = manifests, binding = binding))
}

pc_pair_checks <- function(legacy, jaoua) {
  top_fields <- c("time_obs", "C", "Z")
  truth_fields <- c(
    "phase_mu_mat", "phase_beta", "zeta_true", "xi_true", "eta_true",
    "chi_true", "zeta_fpca_true", "xi_fpca_true", "a_true", "b_true",
    "gamma_a_true", "gamma_b_true", "omega_a_true", "omega_b_true",
    "lambda_phi_true", "lambda_psi_true", "mu_true_values",
    "beta_true_values", "noise_true_values", "n_obs_by_subject",
    "int_knots", "sigma2_eps_true"
  )
  checks <- c(
    setNames(vapply(top_fields, function(field) {
      identical(legacy[[field]], jaoua[[field]])
    }, logical(1L)), paste0("identical_", top_fields)),
    setNames(vapply(truth_fields, function(field) {
      identical(legacy$true_params[[field]], jaoua$true_params[[field]])
    }, logical(1L)), paste0("identical_", truth_fields)),
    basis_methods_distinct =
      !identical(legacy$true_params$truth_basis_method,
                 jaoua$true_params$truth_basis_method),
    shared_basis_distinct = max(abs(
      legacy$true_params$phi_dense_list[[1L]] -
        jaoua$true_params$phi_dense_list[[1L]]
    )) > 1e-6,
    observed_responses_distinct = max(abs(
      unlist(legacy$Y, recursive = TRUE, use.names = FALSE) -
        unlist(jaoua$Y, recursive = TRUE, use.names = FALSE)
    )) > 1e-6
  )
  data.frame(
    check_id = names(checks), passed = unname(checks),
    detail = ifelse(checks, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

pc_save_arm <- function(data, data_row, config, captured, output_root) {
  data_id <- data_row$data_id[[1L]]
  output_dir <- file.path(output_root, "data", data_id)
  pp_assert(!dir.exists(output_dir),
            paste0("Refusing to overwrite data directory: ", output_dir))
  observation <- list(
    bundle_class = "v0lv_observation_only_candidate2",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_id,
    logical_data_id = data_row$logical_data_id[[1L]],
    data_seed = as.integer(data_row$data_seed[[1L]]),
    basis_arm = data_row$basis_arm[[1L]],
    scenario_id = data_row$scenario_id[[1L]],
    scenario_label_available_to_fit = FALSE,
    engineering_pilot = TRUE, development_only = TRUE,
    formal_experiment = FALSE,
    Y = data$Y, time_obs = data$time_obs, Z = data$Z,
    dimensions = config[c("S", "n_s", "p", "d", "L_f", "L_s",
                          "M_f", "M_s", "K")]
  )
  pp_assert(!length(candidate2_forbidden_observation_names(observation)),
            paste0("Observation bundle leaks truth: ", data_id))
  assertions <- pp_truth_assertions(data, data_row, config)
  pp_assert(all(assertions$passed), paste0(
    "Truth assertion failed for ", data_id, ": ",
    paste(assertions$assertion_id[!assertions$passed], collapse = ";")
  ))
  sealed <- list(
    bundle_class = "v0lv_sealed_truth_candidate2",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_id,
    logical_data_id = data_row$logical_data_id[[1L]],
    data_seed = as.integer(data_row$data_seed[[1L]]),
    truth_basis_seed = as.integer(data_row$truth_basis_seed[[1L]]),
    basis_arm = data_row$basis_arm[[1L]],
    truth_basis_method = data_row$truth_basis_method[[1L]],
    scenario_id = data_row$scenario_id[[1L]],
    engineering_pilot = TRUE, development_only = TRUE,
    formal_experiment = FALSE, data = data,
    sealed_generation_provenance = list(
      generator_config = config,
      truth_basis_method = data_row$truth_basis_method[[1L]],
      truth_basis_seed = as.integer(data_row$truth_basis_seed[[1L]]),
      generation_warnings = captured$warnings,
      truth_assertions = assertions,
      generated_at_utc = captured$ended_at_utc
    )
  )
  dir.create(file.path(output_dir, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  observation_path <- file.path(output_dir, "observation_bundle.rds")
  truth_path <- file.path(output_dir, "sealed_truth", "truth_bundle.rds")
  pp_save_rds(observation, observation_path)
  pp_save_rds(sealed, truth_path)
  pp_write_csv(assertions, file.path(output_dir, "sealed_truth",
                                     "truth_assertions.csv"))
  pp_write_lines(c(
    "status=DATA_SEAL_COMPLETE", paste0("data_id=", data_id),
    paste0("observation_sha256=", pp_sha256(observation_path)),
    paste0("sealed_truth_sha256=", pp_sha256(truth_path)),
    "truth_unseal_authorized=FALSE", "development_only=TRUE",
    "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"
  ), file.path(output_dir, "DATA_SEAL_COMPLETE.txt"))
  invisible(TRUE)
}

pc_generate_pair <- function(pair_rows, output_root, smoke = FALSE) {
  pp_assert(nrow(pair_rows) == 2L &&
              setequal(pair_rows$basis_arm, c("legacy", "jaoua")) &&
              length(unique(pair_rows$data_seed)) == 1L &&
              length(unique(pair_rows$truth_basis_seed)) == 1L,
            "Malformed logical data pair.")
  arm_dirs <- file.path(output_root, "data", pair_rows$data_id)
  if (any(dir.exists(arm_dirs))) {
    complete <- file.exists(file.path(arm_dirs, "DATA_SEAL_COMPLETE.txt")) &
      file.exists(file.path(arm_dirs, "observation_bundle.rds")) &
      file.exists(file.path(arm_dirs, "sealed_truth", "truth_bundle.rds"))
    pp_assert(all(complete), "Incomplete pre-existing paired data directory.")
    return(invisible(FALSE))
  }
  config <- pp_base_config(pair_rows[1L, , drop = FALSE], smoke)
  generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
  calibration <- new.env(parent = globalenv())
  sys.source(pp_find_one(PP_SNAPSHOT, "^v0f_data_calibration[.]R$"),
             envir = calibration)
  generated <- vector("list", 2L)
  captures <- vector("list", 2L)
  names(generated) <- names(captures) <- pair_rows$basis_arm
  for (index in seq_len(2L)) {
    row <- pair_rows[index, , drop = FALSE]
    arguments <- config
    arguments$seed <- as.integer(row$data_seed[[1L]])
    arguments$truth_basis_seed <- as.integer(row$truth_basis_seed[[1L]])
    arguments$truth_basis_method <- row$truth_basis_method[[1L]]
    captured <- c2_capture_conditions(
      do.call(generator, arguments), phase = paste0(
        "paired_generation_", row$basis_arm[[1L]]
      )
    )
    pp_assert(is.null(captured$error), paste0(
      "Generation failed for ", row$data_id[[1L]], ": ",
      captured$error$message %||% "unknown"
    ))
    generated[[row$basis_arm[[1L]]]] <-
      calibration$calibrate_loading_information_v0f(
        captured$value, regime = "constant_per_active_loading",
        anchor_p = 100L, active_fraction = 0.1, anchor_norm = 3
      )
    captures[[row$basis_arm[[1L]]]] <- captured
  }
  pairing <- pc_pair_checks(generated$legacy, generated$jaoua)
  pp_assert(all(pairing$passed), paste0(
    "Non-basis pairing failed: ",
    paste(pairing$check_id[!pairing$passed], collapse = ";")
  ))
  for (index in seq_len(2L)) {
    row <- pair_rows[index, , drop = FALSE]
    pc_save_arm(generated[[row$basis_arm[[1L]]]], row, config,
                captures[[row$basis_arm[[1L]]]], output_root)
  }
  qc_dir <- file.path(output_root, "pairing_qc")
  pp_write_csv(pairing, file.path(
    qc_dir, paste0(pair_rows$logical_data_id[[1L]], "_PAIRING_QC.csv")
  ))
  invisible(TRUE)
}

pc_generate_all <- function(output_root, smoke = FALSE) {
  verification <- pc_verify_environment()
  manifests <- verification$manifests
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  pp_assert(!file.exists(file.path(output_root,
                                  "FITS_AND_SELECTION_COMPLETE.txt")),
            "Refusing generation after fit/selection completion.")
  groups <- split(manifests$data, manifests$data$logical_data_id)
  created <- vapply(groups, pc_generate_pair, logical(1L),
                    output_root = output_root, smoke = smoke)
  status <- data.frame(
    logical_data_id = names(groups), newly_created = unname(created),
    paired_arms_complete = TRUE, truth_unseal_authorized = FALSE,
    smoke = smoke, checked_at_utc = pp_iso_time(),
    stringsAsFactors = FALSE
  )
  pp_write_csv(status, file.path(output_root, "DATA_GENERATION_STATUS.csv"))
  pp_write_lines(c(
    "status=PAIRED_DATA_GENERATION_AND_SEAL_COMPLETE",
    paste0("completed_utc=", pp_iso_time()),
    "logical_data_sets=3", "observation_bundles=6",
    "sealed_truth_bundles=6", "pairing_qc_passed=TRUE",
    "truth_unseal_authorized=FALSE",
    paste0("smoke=", toupper(as.character(smoke))),
    "development_only=TRUE", "formal_v0lv_result=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("PAIRED_DATA_GENERATION_COMPLETE logical=3 arms=6 smoke=",
      smoke, "\n", sep = "")
  invisible(status)
}

pc_run_fit <- function(fit_id, output_root, smoke = FALSE) {
  manifests <- pc_verify_environment()$manifests
  pp_assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
            "Paired data generation/sealing is incomplete.")
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
  result <- pp_fit_g12(observation, row, smoke)
  record <- pp_make_terminal_record(result, row, smoke)
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
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
  }
  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(), call = character(),
               phase = character(), stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (record$terminal_status == "error") {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", record$error$message %||% "unknown"),
      paste0("ended_utc=", pp_iso_time())
    ), file.path(fit_dir, "FIT_ERROR.txt"))
    pp_abort("Fit failed and was retained: ", fit_id, ": ",
             record$error$message %||% "unknown")
  }
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("basis_arm=", row$basis_arm[[1L]]),
    paste0("terminal_status=", record$terminal_status),
    paste0("objective_eligible=", record$objective_eligible),
    "truth_used_for_fit_stopping_or_selection=FALSE", "n_cpus=1",
    "continuation_used=FALSE", "automatic_800_used=FALSE",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(0L)
}

pc_collect_terminals <- function(manifest, output_root) {
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  pp_assert(all(file.exists(paths)), "One or more terminal records are missing.")
  records <- lapply(paths, readRDS)
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  pp_assert(nrow(rows) == nrow(manifest) && !anyDuplicated(rows$fit_id) &&
              identical(rows$fit_id, manifest$fit_id),
            "Terminal records do not exactly join to the manifest.")
  rows
}

pc_freeze_selection <- function(manifest, output_root, smoke = FALSE) {
  rows <- pc_collect_terminals(manifest, output_root)
  pp_assert(all(rows$terminal_status != "error") &&
              all(rows$objective_eligible) &&
              !any(rows$truth_used_for_fit_stopping_or_selection) &&
              !any(rows$continuation_used),
            "Selection pool contains an invalid or truth-bearing endpoint.")
  winners <- pp_select_g12_winners(rows)
  expected <- length(unique(manifest$selection_stratum_id))
  pp_assert(nrow(winners) == expected &&
              all(table(winners$selection_stratum_id) == 1L),
            "Truth-free paired winner table is invalid.")
  selection_dir <- file.path(output_root, "truth_free_selection")
  pp_assert(!dir.exists(selection_dir),
            "Truth-free selection already exists; refusing overwrite.")
  dir.create(selection_dir, recursive = TRUE, mode = "0700")
  endpoints_path <- file.path(selection_dir, "ALL_TRUTH_FREE_ENDPOINTS.csv")
  winners_path <- file.path(selection_dir, "G12_TRUTH_FREE_WINNERS.csv")
  pp_write_csv(rows, endpoints_path)
  pp_write_csv(winners, winners_path)
  selection_hash <- pp_sha256(winners_path)
  pp_write_lines(c(
    "status=TRUTH_FREE_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("endpoints=", nrow(rows)), paste0("winners=", nrow(winners)),
    paste0("selection_sha256=", selection_hash),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_arm",
    "tie_break_rule=fit_id_lexicographic",
    "cross_arm_ELBO_comparison=FALSE", "selection_used_truth=FALSE",
    "truth_unseal_authorized=FALSE", "continuation_started=FALSE",
    paste0("smoke=", toupper(as.character(smoke)))
  ), file.path(selection_dir, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  if (!smoke) {
    pp_write_lines(c(
      "status=FITS_COMPLETE_AND_TRUTH_FREE_SELECTION_FROZEN",
      paste0("completed_utc=", pp_iso_time()), "fits=72", "winners=6",
      paste0("outer_workers=", PP_WORKERS), "per_fit_n_cpus=1",
      "truth_unseal_authorized=FALSE", "evaluation_started=FALSE",
      "continuation_started=FALSE", "development_only=TRUE",
      "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"
    ), file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
    pp_write_lines(c(
      "status=WAITING_FOR_EXPLICIT_TRUTH_UNSEAL_AUTHORIZATION",
      paste0("selection_sha256=", selection_hash),
      "truth_unseal_authorized=FALSE", "evaluation_started=FALSE",
      "new_fits_authorized=FALSE", "continuation_authorized=FALSE"
    ), file.path(output_root,
                 "WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt"))
  }
  cat("TRUTH_FREE_SELECTION_FROZEN endpoints=", nrow(rows),
      " winners=", nrow(winners), "\n", sep = "")
  invisible(list(rows = rows, winners = winners, hash = selection_hash))
}

pc_count_markers <- function(output_root, marker) {
  root <- file.path(output_root, "fits")
  if (!dir.exists(root)) return(0L)
  length(list.files(root, pattern = paste0("^", marker, "$"),
                    recursive = TRUE))
}

pc_write_status <- function(output_root, state, detail = "") {
  status <- data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    total_fits = 72L,
    complete_fits = pc_count_markers(output_root, "FIT_COMPLETE[.]txt"),
    error_fits = pc_count_markers(output_root, "FIT_ERROR[.]txt"),
    outer_workers = PP_WORKERS, per_fit_n_cpus = 1L,
    truth_read_by_fit_or_selection = FALSE,
    development_only = TRUE, formal_v0lv_result = FALSE,
    stringsAsFactors = FALSE
  )
  pp_write_csv(status, file.path(output_root, "CURRENT_STATUS.csv"))
  invisible(status)
}

pc_run_batch <- function(manifest, output_root, label) {
  pp_assert(.Platform$OS.type == "unix",
            "Concurrent batches require a Unix/Linux host.")
  log_dir <- file.path(output_root, "supervisor_logs", label)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    started <- Sys.time()
    status <- system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=fit", paste0("--fit-id=", fit_id),
                    paste0("--output-root=", output_root)),
      stdout = stdout, stderr = stderr
    )
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
    manifest$fit_id, run_one, mc.cores = min(PP_WORKERS, nrow(manifest)),
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    pc_write_status(output_root, "ERROR",
                    paste0("retained=", paste(failed, collapse = ";")))
    pp_abort("Fit errors retained without selective rerun: ",
             paste(failed, collapse = ", "))
  }
  rows <- pc_collect_terminals(manifest, output_root)
  pp_assert(all(rows$objective_eligible) &&
              all(rows$terminal_status != "error") &&
              all(rows$actual_T1_sweeps == 400L) &&
              !any(rows$truth_used_for_fit_stopping_or_selection),
            paste0("Batch endpoint validation failed: ", label))
  invisible(rows)
}

pc_technical_gate <- function(output_root) {
  manifests <- pc_verify_environment()$manifests
  pp_assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
            "Generate paired data before the technical gate.")
  gate <- manifests$fits[manifests$fits$technical_gate, , drop = FALSE]
  pp_assert(nrow(gate) == 8L, "Technical gate must contain eight fits.")
  pc_write_status(output_root, "TECHNICAL_GATE_RUNNING", "8 registered fits")
  rows <- pc_run_batch(gate, output_root, "technical_gate")
  pp_write_lines(c(
    "status=TECHNICAL_GATE_PASS", paste0("completed_utc=", pp_iso_time()),
    "registered_fits=8", "complete_fits=8", "objective_eligible=8",
    "fits_are_part_of_final_72=TRUE", "truth_read=FALSE",
    "selective_rerun=FALSE", "development_only=TRUE"
  ), file.path(output_root, "TECHNICAL_GATE_COMPLETE.txt"))
  pc_write_status(output_root, "TECHNICAL_GATE_PASS", "8/8 eligible")
  cat("TECHNICAL_GATE_PASS fits=8\n")
  invisible(rows)
}

pc_run_remaining <- function(output_root) {
  manifests <- pc_verify_environment()$manifests
  pp_assert(file.exists(file.path(output_root, "TECHNICAL_GATE_COMPLETE.txt")),
            "Technical gate has not passed.")
  remaining <- manifests$fits[!manifests$fits$technical_gate, , drop = FALSE]
  pp_assert(nrow(remaining) == 64L, "Remaining batch must contain 64 fits.")
  pc_write_status(output_root, "REMAINING_FITS_RUNNING", "64 registered fits")
  pc_run_batch(remaining, output_root, "remaining_64")
  pc_freeze_selection(manifests$fits, output_root, FALSE)
  pc_write_status(output_root, "WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION",
                  "72/72 complete; six truth-free winners frozen")
  invisible(TRUE)
}

pc_micro_smoke <- function(output_root) {
  pp_assert(!dir.exists(output_root),
            "Micro-smoke output root must not already exist.")
  pc_generate_all(output_root, TRUE)
  manifests <- pp_load_manifests()
  smoke_rows <- manifests$fits[
    manifests$fits$logical_data_index == 1L &
      manifests$fits$seed_index == 1L, , drop = FALSE]
  pp_assert(nrow(smoke_rows) == 2L &&
              setequal(smoke_rows$basis_arm, c("legacy", "jaoua")),
            "Micro-smoke pair is not exactly two arms.")
  for (fit_id in smoke_rows$fit_id) {
    pc_run_fit(fit_id, output_root, TRUE)
  }
  selection <- pc_freeze_selection(smoke_rows, output_root, TRUE)
  pp_assert(nrow(selection$rows) == 2L &&
              all(selection$rows$objective_eligible) &&
              nrow(selection$winners) == 2L,
            "Micro-smoke endpoint contract failed.")
  pp_write_lines(c(
    "status=MICRO_SMOKE_PASS", paste0("completed_utc=", pp_iso_time()),
    "logical_data_pairs_generated=3", "observation_bundles=6",
    "fits=2", "objective_eligible=2", "truth_free_winners=2",
    "pairing_qc_passed=TRUE", "truth_unseal_authorized=FALSE",
    "formal_data_generated=FALSE", "formal_fits_started=FALSE",
    "development_only=TRUE"
  ), file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
  cat("MICRO_SMOKE_PASS data_arms=6 fits=2 winners=2\n")
  invisible(TRUE)
}

pc_status <- function(output_root) {
  complete <- pc_count_markers(output_root, "FIT_COMPLETE[.]txt")
  errors <- pc_count_markers(output_root, "FIT_ERROR[.]txt")
  gate <- file.exists(file.path(output_root, "TECHNICAL_GATE_COMPLETE.txt"))
  selection <- file.exists(file.path(
    output_root, "truth_free_selection", "TRUTH_FREE_SELECTION_FROZEN.txt"
  ))
  cat("complete_fits=", complete, "\n", sep = "")
  cat("error_fits=", errors, "\n", sep = "")
  cat("technical_gate_complete=", toupper(as.character(gate)), "\n", sep = "")
  cat("truth_free_selection_frozen=",
      toupper(as.character(selection)), "\n", sep = "")
  cat("truth_unseal_authorized=FALSE\n")
  invisible(TRUE)
}

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  verification <- pc_verify_environment()
  cat("CHECK_PASS source_rows=", nrow(verification$binding),
      " data=6 fits=72 gate=8 workers=", PP_WORKERS, "\n", sep = "")
} else if (action == "micro-smoke") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_micro_smoke(output_root)
} else if (action == "generate") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_generate_all(output_root, FALSE)
} else if (action == "fit") {
  pp_assert(nzchar(output_root) && nzchar(cli$`fit-id` %||% ""),
            "--output-root and --fit-id are required.")
  pc_run_fit(cli$`fit-id`, output_root, pp_bool(cli$smoke, FALSE))
} else if (action == "technical-gate") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_technical_gate(output_root)
} else if (action == "remaining") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_run_remaining(output_root)
} else if (action == "full") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_generate_all(output_root, FALSE)
  if (!file.exists(file.path(output_root, "TECHNICAL_GATE_COMPLETE.txt"))) {
    pc_technical_gate(output_root)
  }
  pc_run_remaining(output_root)
} else if (action == "status") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  pc_status(output_root)
} else {
  pp_abort("Unknown action: ", action)
}
