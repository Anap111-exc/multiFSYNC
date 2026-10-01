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
PP_EXPERIMENT_ID <- "G12_LOADING_DGP_ROBUSTNESS_STAGE6L_V1_20260915"
PP_DATA_IDS <- c("s6l_01__controlled", "s6l_01__jaoua_beta10")
PP_FIT_CONFIG_IDS <- "correct_L1_M2"
PP_EXPECTED_FITS <- 8L
PP_EXPECTED_WINNERS <- 2L
PP_WORKERS <- as.integer(pp_env("STAGE6L_OUTER_WORKERS", "8"))
pp_assert(length(PP_WORKERS) == 1L && is.finite(PP_WORKERS) &&
            PP_WORKERS >= 1L && PP_WORKERS <= 8L,
          "STAGE6L_OUTER_WORKERS must be an integer from 1 through 8.")

parent_make_terminal_record <- pp_make_terminal_record
parent_terminal_row <- pp_terminal_row

pp_validate_data_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "data_id", "logical_data_id", "loading_arm",
    "scenario_id", "data_seed", "truth_basis_seed", "truth_loading_seed",
    "truth_basis_method", "loading_calibration", "S", "n_s", "p", "d",
    "truth_L_f", "truth_L_s", "truth_M_f", "truth_M_s", "K",
    "n_obs_min", "n_obs_max", "sigma_eps", "mean_amp",
    "new_unseen_at_registration", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "truth_unseal_authorized", "development_only", "formal_v0lv_result",
    "formal_paper_mc_result")
  pp_assert(is.data.frame(manifest) && nrow(manifest) == 2L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              identical(as.character(manifest$data_id), PP_DATA_IDS) &&
              !anyDuplicated(manifest$data_id) &&
              all(manifest$logical_data_id == "s6l_01") &&
              setequal(manifest$loading_arm, c("controlled", "jaoua_beta10")) &&
              length(unique(manifest$data_seed)) == 1L &&
              identical(unique(manifest$data_seed), 96151001L) &&
              length(unique(manifest$truth_basis_seed)) == 1L &&
              identical(unique(manifest$truth_basis_seed), 97151001L) &&
              length(unique(manifest$truth_loading_seed)) == 1L &&
              identical(unique(manifest$truth_loading_seed), 98151001L) &&
              all(manifest$truth_basis_method == "jaoua_splinets") &&
              all(manifest$loading_calibration ==
                    "v0f_constant_per_active_loading") &&
              all(manifest$new_unseen_at_registration) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              !any(manifest$truth_unseal_authorized) &&
              all(manifest$development_only) &&
              !any(manifest$formal_v0lv_result) &&
              !any(manifest$formal_paper_mc_result),
            "Stage 6L data-manifest identity or truth-isolation contract failed.")
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
            "Stage 6L data dimensions or generation parameters changed.")
  invisible(TRUE)
}

pp_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "logical_data_id",
    "loading_arm", "scenario_id", "method_id", "fit_config_id",
    "selection_stratum_id", "paired_start_id", "seed_index", "fit_seed",
    "fit_L_f", "fit_L_s_1", "fit_L_s_2", "fit_M_f", "fit_M_s_1",
    "fit_M_s_2", "route_id", "initialization", "function_initialization",
    "pre_score_sweeps", "anneal", "planned_annealing_sweeps",
    "planned_ordinary_t1_sweeps", "n_cpus",
    "objective_eligibility_required", "truth_free_winner_candidate",
    "cross_arm_elbo_comparison_forbidden", "actual_racing", "continuation",
    "automatic_800", "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_selection", "development_only", "formal_v0lv_result",
    "formal_paper_mc_result")
  pp_assert(is.data.frame(manifest) && nrow(manifest) == PP_EXPECTED_FITS &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              !anyDuplicated(manifest$fit_id) &&
              setequal(manifest$data_id, PP_DATA_IDS) &&
              all(table(manifest$data_id) == 4L) &&
              all(table(manifest$paired_start_id) == 2L),
            "Stage 6L fit-manifest identity or pairing failed.")
  pairing_ok <- vapply(split(manifest, manifest$paired_start_id), function(x) {
    nrow(x) == 2L && setequal(x$loading_arm, c("controlled", "jaoua_beta10")) &&
      length(unique(x$fit_seed)) == 1L
  }, logical(1L))
  pp_assert(all(pairing_ok) && all(manifest$method_id == "G12") &&
              all(manifest$fit_config_id == "correct_L1_M2") &&
              all(manifest$fit_L_f == 1L) &&
              all(manifest$fit_L_s_1 == 1L) &&
              all(manifest$fit_L_s_2 == 1L) &&
              all(manifest$fit_M_f == "2") &&
              all(manifest$fit_M_s_1 == "2") &&
              all(manifest$fit_M_s_2 == "2") &&
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
            "Stage 6L route, dimension, isolation or execution contract failed.")
  invisible(TRUE)
}

pp_base_config <- function(data_row, smoke = FALSE) {
  config <- if (smoke) list(
    S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
    K = 3L, n_obs = 5L, common_grid = FALSE,
    n_obs_range = c(4L, 6L), sigma_eps = 0.3,
    bool_sparse_loadings = TRUE, prop_sparse = 0.8,
    score_var_decay = TRUE, n_dense = 51L, mean_amp = 0.6,
    beta_amp = 0, bs_degree = c(2L, 3L), omega_beta = c(1, 10),
    truth_basis_method = "jaoua_splinets"
  ) else list(
    S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
    K = 5L, n_obs = 8L, common_grid = FALSE,
    n_obs_range = c(6L, 9L), sigma_eps = 0.3,
    bool_sparse_loadings = TRUE, prop_sparse = 0.9,
    score_var_decay = TRUE, n_dense = 201L, mean_amp = 0.6,
    beta_amp = 0, bs_degree = c(2L, 3L), omega_beta = c(1, 10),
    truth_basis_method = "jaoua_splinets"
  )
  if (data_row$loading_arm[[1L]] == "controlled") {
    config$identified_loadings <- TRUE
    config$sparsity_mode <- "fixed"
  } else {
    config$identified_loadings <- FALSE
    config$sparsity_mode <- "jaoua_beta_iid_normal"
  }
  config
}

s6l_calibrate_constant_per_active <- function(data, anchor_active = 10L,
                                              anchor_norm = 3) {
  truth <- data$true_params
  truth$generator_omega_a_true <- truth$omega_a_true
  truth$generator_omega_b_true <- truth$omega_b_true
  calibrate <- function(value) {
    active <- sum(value != 0)
    norm <- sqrt(sum(value^2))
    pp_assert(active >= 1L && is.finite(norm) && norm > 1e-14,
              "Cannot calibrate an empty/non-finite loading column.")
    anchor_norm * sqrt(active / anchor_active) * value / norm
  }
  truth$a_true[, 1L] <- calibrate(truth$a_true[, 1L])
  for (study in seq_len(truth$S)) {
    truth$b_true[[study]][, 1L] <- calibrate(truth$b_true[[study]][, 1L])
  }
  truth$gamma_a_true <- 1L * (truth$a_true != 0)
  truth$gamma_b_true <- lapply(truth$b_true, function(x) 1L * (x != 0))
  truth$omega_a_true <- colMeans(truth$gamma_a_true)
  truth$omega_b_true <- lapply(truth$gamma_b_true, colMeans)
  truth$loading_calibration <- "v0f_constant_per_active_loading"
  truth$v0f_anchor_active <- as.integer(anchor_active)
  truth$v0f_anchor_loading_norm <- anchor_norm
  for (study in seq_len(truth$S)) for (subject in seq_len(truth$n_s[[study]])) {
    signal <- truth$mu_true_values[[study]][[subject]] +
      truth$beta_true_values[[study]][[subject]] +
      truth$f_true_values[[study]][[subject]] %*% t(truth$a_true) +
      truth$g_true_values[[study]][[subject]] %*% t(truth$b_true[[study]])
    truth$signal_true_values[[study]][[subject]] <- signal
    noise <- truth$noise_true_values[[study]][[subject]]
    data$Y[[study]][[subject]] <- lapply(
      seq_len(truth$p),
      function(variable) as.numeric(signal[, variable] + noise[, variable]))
  }
  data$true_params <- truth
  data
}

s6l_abs_cosine <- function(left, right) {
  denom <- sqrt(sum(left^2) * sum(right^2))
  if (!is.finite(denom) || denom <= 1e-14) return(NA_real_)
  abs(sum(left * right) / denom)
}

pp_truth_assertions <- function(data, data_row, config) {
  truth <- data$true_params
  counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  loading <- c(list(as.numeric(truth$a_true[, 1L])),
               lapply(truth$b_true, function(x) as.numeric(x[, 1L])))
  active <- vapply(loading, function(x) sum(x != 0), integer(1L))
  norms <- vapply(loading, function(x) sqrt(sum(x^2)), numeric(1L))
  target <- 3 * sqrt(active / 10)
  cosines <- vapply(1:2, function(study) {
    s6l_abs_cosine(loading[[1L]], loading[[study + 1L]])
  }, numeric(1L))
  min_eigen <- min(vapply(1:2, function(study) {
    min(eigen(crossprod(cbind(loading[[1L]], loading[[study + 1L]])),
              symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1L)))
  maximum_error <- 0
  for (study in seq_len(truth$S)) for (subject in seq_len(truth$n_s[[study]])) {
    signal <- truth$mu_true_values[[study]][[subject]] +
      truth$beta_true_values[[study]][[subject]] +
      truth$f_true_values[[study]][[subject]] %*% t(truth$a_true) +
      truth$g_true_values[[study]][[subject]] %*% t(truth$b_true[[study]])
    observed <- do.call(cbind, data$Y[[study]][[subject]])
    maximum_error <- max(
      maximum_error,
      abs(signal - truth$signal_true_values[[study]][[subject]]),
      abs(observed - signal - truth$noise_true_values[[study]][[subject]]))
  }
  controlled <- data_row$loading_arm[[1L]] == "controlled"
  checks <- c(
    dimensions = truth$S == 2L && all(truth$n_s == config$n_s) &&
      truth$p == config$p && truth$d_use == 0L,
    factor_fpca_counts = truth$L_f == 1L &&
      all(truth$L_s_by_study == 1L) && all(truth$M_f == 2L) &&
      all(vapply(truth$M_s, function(x) all(x == 2L), logical(1L))),
    observation_counts = all(counts >= config$n_obs_range[[1L]] &
                               counts <= config$n_obs_range[[2L]]),
    nonempty_loading_blocks = all(active >= 1L),
    calibration_exact = max(abs(norms - target)) < 1e-12,
    controlled_geometry = if (controlled) {
      all(active == round(config$p * (1 - config$prop_sparse))) &&
        max(cosines) < 1e-12
    } else TRUE,
    beta_geometry_full_rank = if (!controlled) min_eigen > 1e-8 else TRUE,
    signal_noise_reconstruction = maximum_error < 1e-12)
  data.frame(
    assertion_id = names(checks), passed = unname(checks),
    detail = c(
      paste0("S=", truth$S, ";n=", paste(truth$n_s, collapse = ";"),
             ";p=", truth$p),
      "L=1;1;1;M=2;2;2", paste(range(counts), collapse = ";"),
      paste(active, collapse = ";"),
      format(max(abs(norms - target)), scientific = TRUE),
      paste(format(cosines, digits = 8), collapse = ";"),
      format(min_eigen, scientific = TRUE),
      format(maximum_error, scientific = TRUE)),
    stringsAsFactors = FALSE)
}

pp_make_terminal_record <- function(result, fit_row, smoke = FALSE) {
  value <- parent_make_terminal_record(result, fit_row, smoke)
  value$terminal_schema <- "G12_LOADING_DGP_ROBUSTNESS_STAGE6L_TERMINAL_V1"
  value$logical_data_id <- fit_row$logical_data_id[[1L]]
  value$loading_arm <- fit_row$loading_arm[[1L]]
  value$paired_start_id <- fit_row$paired_start_id[[1L]]
  value
}

pp_terminal_row <- function(record) {
  value <- parent_terminal_row(record)
  value$logical_data_id <- record$logical_data_id
  value$loading_arm <- record$loading_arm
  value$paired_start_id <- record$paired_start_id
  value
}

s6l_verify_environment <- function() {
  pp_assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
  pp_assert(requireNamespace("Splinets", quietly = TRUE), "Splinets is required.")
  binding <- pp_verify_source_binding()
  pp_load_package()
  pp_source_runtime(FALSE)
  manifests <- pp_load_manifests()
  invisible(list(binding = binding, manifests = manifests))
}

s6l_pair_checks <- function(controlled, jaoua) {
  top_fields <- c("time_obs", "C", "Z")
  truth_fields <- c(
    "phase_mu_mat", "phase_beta", "zeta_true", "xi_true", "eta_true",
    "chi_true", "zeta_fpca_true", "xi_fpca_true", "phi_dense_list",
    "psi_dense_list", "theta_dense_list", "kappa_dense_list",
    "lambda_phi_true", "lambda_psi_true", "mu_true_values",
    "f_true_values", "g_true_values", "beta_true_values",
    "noise_true_values", "n_obs_by_subject", "int_knots", "sigma2_eps_true")
  checks <- c(
    setNames(vapply(top_fields, function(field) {
      identical(controlled[[field]], jaoua[[field]])
    }, logical(1L)), paste0("identical_", top_fields)),
    setNames(vapply(truth_fields, function(field) {
      identical(controlled$true_params[[field]], jaoua$true_params[[field]])
    }, logical(1L)), paste0("identical_", truth_fields)),
    loading_differs = !identical(controlled$true_params$a_true,
                                 jaoua$true_params$a_true),
    response_differs = !identical(controlled$Y, jaoua$Y))
  data.frame(check_id = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}

s6l_save_arm <- function(data, row, config, captured, output_root) {
  data_id <- row$data_id[[1L]]
  output_dir <- file.path(output_root, "data", data_id)
  pp_assert(!dir.exists(output_dir), paste0("Refusing overwrite: ", output_dir))
  assertions <- pp_truth_assertions(data, row, config)
  pp_assert(all(assertions$passed), paste0(
    "Truth assertion failure for ", data_id, ": ",
    paste(assertions$assertion_id[!assertions$passed], collapse = ";")))
  observation <- list(
    bundle_class = "v0lv_observation_only_candidate2",
    experiment_id = PP_EXPERIMENT_ID, data_id = data_id,
    logical_data_id = row$logical_data_id[[1L]],
    data_seed = as.integer(row$data_seed[[1L]]),
    scenario_label_available_to_fit = FALSE,
    engineering_pilot = TRUE, development_only = TRUE,
    formal_experiment = FALSE, Y = data$Y, time_obs = data$time_obs,
    Z = data$Z,
    dimensions = config[c("S", "n_s", "p", "d", "L_f", "L_s",
                          "M_f", "M_s", "K")])
  pp_assert(!length(candidate2_forbidden_observation_names(observation)),
            paste0("Observation bundle leaks truth: ", data_id))
  sealed <- list(
    bundle_class = "v0lv_sealed_truth_candidate2",
    experiment_id = PP_EXPERIMENT_ID, data_id = data_id,
    logical_data_id = row$logical_data_id[[1L]],
    loading_arm = row$loading_arm[[1L]], scenario_id = row$scenario_id[[1L]],
    data_seed = as.integer(row$data_seed[[1L]]),
    truth_basis_seed = as.integer(row$truth_basis_seed[[1L]]),
    truth_loading_seed = as.integer(row$truth_loading_seed[[1L]]),
    development_only = TRUE, formal_experiment = FALSE, data = data,
    sealed_generation_provenance = list(
      generator_config = config,
      loading_calibration = "v0f_constant_per_active_loading",
      generation_warnings = captured$warnings,
      truth_assertions = assertions,
      generated_at_utc = captured$ended_at_utc))
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
    "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE"),
    file.path(output_dir, "DATA_SEAL_COMPLETE.txt"))
}

s6l_generate <- function(output_root, smoke = FALSE) {
  verification <- s6l_verify_environment()
  manifest <- verification$manifests$data
  pp_assert(!dir.exists(output_root), "Output root already exists.")
  dir.create(output_root, recursive = TRUE, mode = "0700")
  generated <- captures <- vector("list", nrow(manifest))
  names(generated) <- names(captures) <- manifest$loading_arm
  generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    config <- pp_base_config(row, smoke)
    arguments <- c(config, list(
      seed = as.integer(row$data_seed[[1L]]),
      truth_basis_seed = as.integer(row$truth_basis_seed[[1L]]),
      truth_loading_seed = as.integer(row$truth_loading_seed[[1L]])))
    captured <- c2_capture_conditions(
      do.call(generator, arguments),
      phase = paste0("stage6l_generation_", row$loading_arm[[1L]]))
    pp_assert(is.null(captured$error), paste0(
      "Generation failed for ", row$data_id[[1L]], ": ",
      captured$error$message %||% "unknown"))
    generated[[row$loading_arm[[1L]]]] <-
      s6l_calibrate_constant_per_active(captured$value)
    captures[[row$loading_arm[[1L]]]] <- captured
  }
  pairing <- s6l_pair_checks(generated$controlled, generated$jaoua_beta10)
  pp_assert(all(pairing$passed), paste0(
    "Non-loading pairing failed: ",
    paste(pairing$check_id[!pairing$passed], collapse = ";")))
  for (index in seq_len(nrow(manifest))) {
    row <- manifest[index, , drop = FALSE]
    s6l_save_arm(generated[[row$loading_arm[[1L]]]], row,
                 pp_base_config(row, smoke), captures[[row$loading_arm[[1L]]]],
                 output_root)
  }
  pp_write_csv(pairing, file.path(output_root, "PAIRING_QC.csv"))
  pp_write_lines(c(
    "status=DATA_GENERATION_AND_SEAL_COMPLETE", paste0("smoke=", smoke),
    "logical_data_sets=1", "observation_bundles=2",
    "sealed_truth_bundles=2", "pairing_qc_passed=TRUE",
    "truth_unseal_authorized=FALSE", "formal_v0lv_result=FALSE",
    "formal_paper_mc_result=FALSE", paste0("completed_utc=", pp_iso_time())),
    file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  invisible(TRUE)
}

s6l_run_fit <- function(fit_id, output_root, smoke = FALSE) {
  manifest <- s6l_verify_environment()$manifests$fits
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  pp_assert(nrow(row) == 1L, "fit_id must join to exactly one manifest row.")
  observation <- pp_load_observation(row$data_id[[1L]], output_root)
  fit_dir <- file.path(output_root, "fits", fit_id)
  pp_assert(!dir.exists(fit_dir), paste0("Fit directory exists: ", fit_id))
  dir.create(fit_dir, recursive = TRUE, mode = "0700")
  result <- pp_fit_g12(observation, row, smoke)
  record <- pp_make_terminal_record(result, row, smoke)
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  warnings <- if (is.data.frame(record$warnings)) record$warnings else
    data.frame(class = character(), message = character(), call = character(),
               phase = character(), stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", record$error$message %||% "unknown")),
      file.path(fit_dir, "FIT_ERROR.txt"))
    pp_abort("Fit failed and was retained: ", fit_id)
  }
  pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  pp_write_csv(data.frame(
    t1_sweep = seq_along(result$fit$ELBO),
    ordinary_elbo = as.numeric(result$fit$ELBO)),
    file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  pp_write_csv(result$fit$practical_diagnostics,
               file.path(fit_dir, "PRACTICAL_DIAGNOSTICS.csv"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("loading_arm=", row$loading_arm[[1L]]),
    paste0("terminal_status=", record$terminal_status),
    paste0("objective_eligible=", record$objective_eligible),
    "truth_used_for_fit_stopping_or_selection=FALSE", "n_cpus=1",
    "continuation_used=FALSE", "automatic_800_used=FALSE",
    paste0("ended_utc=", pp_iso_time())),
    file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      record$elapsed_seconds, "\n", sep = "")
  invisible(TRUE)
}

s6l_collect <- function(manifest, output_root) {
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  pp_assert(all(file.exists(paths)), "One or more terminal records are missing.")
  records <- lapply(paths, readRDS)
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  pp_assert(identical(rows$fit_id, manifest$fit_id) &&
              all(rows$objective_eligible) &&
              all(rows$terminal_status != "error") &&
              !any(rows$truth_used_for_fit_stopping_or_selection) &&
              !any(rows$continuation_used),
            "Terminal join or eligibility validation failed.")
  rows
}

s6l_freeze_selection <- function(manifest, output_root, smoke = FALSE) {
  rows <- s6l_collect(manifest, output_root)
  winners <- pp_select_g12_winners(rows)
  pp_assert(nrow(winners) == 2L && all(table(winners$loading_arm) == 1L),
            "Stage 6L winner table is invalid.")
  selection_root <- file.path(output_root, "truth_free_selection")
  pp_assert(!dir.exists(selection_root), "Selection already exists.")
  dir.create(selection_root, recursive = TRUE, mode = "0700")
  endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
  winners_path <- file.path(selection_root, "G12_TRUTH_FREE_WINNERS.csv")
  pp_write_csv(rows, endpoints_path)
  pp_write_csv(winners, winners_path)
  pp_write_lines(c(
    "status=TRUTH_FREE_SELECTION_FROZEN",
    paste0("endpoints_sha256=", pp_sha256(endpoints_path)),
    paste0("winners_sha256=", pp_sha256(winners_path)),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_loading_arm",
    "tie_break_rule=fit_id_lexicographic",
    "cross_arm_elbo_comparison=FALSE", "selection_used_truth=FALSE",
    "truth_unseal_authorized=FALSE", paste0("smoke=", smoke),
    paste0("completed_utc=", pp_iso_time())),
    file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  invisible(list(rows = rows, winners = winners))
}

s6l_run_batch <- function(manifest, output_root) {
  pp_assert(.Platform$OS.type == "unix",
            "The full concurrent bridge requires Linux.")
  log_root <- file.path(output_root, "supervisor_logs")
  dir.create(log_root, recursive = TRUE, mode = "0700")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    status <- system2(
      "Rscript",
      c(PP_SCRIPT_FILE, "--action=fit", paste0("--fit-id=", fit_id),
        paste0("--output-root=", output_root)),
      stdout = file.path(log_root, paste0(key, ".stdout.log")),
      stderr = file.path(log_root, paste0(key, ".stderr.log")))
    list(fit_id = fit_id, status = as.integer(status))
  }
  results <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = min(PP_WORKERS, nrow(manifest)),
    mc.preschedule = FALSE, mc.set.seed = FALSE)
  status <- vapply(results, `[[`, integer(1L), "status")
  if (any(status != 0L)) {
    failed <- vapply(results[status != 0L], `[[`, character(1L), "fit_id")
    pp_abort("Fit errors retained without selective rerun: ",
             paste(failed, collapse = ";"))
  }
  rows <- s6l_collect(manifest, output_root)
  pp_assert(all(rows$actual_T1_sweeps == 400L),
            "A full endpoint did not complete 400 ordinary sweeps.")
  rows
}

s6l_full <- function(output_root) {
  manifests <- s6l_verify_environment()$manifests
  s6l_generate(output_root, FALSE)
  s6l_run_batch(manifests$fits, output_root)
  selection <- s6l_freeze_selection(manifests$fits, output_root, FALSE)
  pp_write_lines(c(
    "status=FITS_COMPLETE_AND_TRUTH_FREE_SELECTION_FROZEN",
    "fits=8", "winners=2", paste0("outer_workers=", PP_WORKERS),
    "per_fit_n_cpus=1", "truth_unseal_authorized=FALSE",
    "evaluation_started=FALSE", "continuation_started=FALSE",
    "development_only=TRUE", "formal_v0lv_result=FALSE",
    "formal_paper_mc_result=FALSE", paste0("completed_utc=", pp_iso_time())),
    file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  invisible(selection)
}

s6l_smoke <- function(output_root) {
  manifests <- s6l_verify_environment()$manifests
  s6l_generate(output_root, TRUE)
  rows <- manifests$fits[manifests$fits$seed_index == 1L, , drop = FALSE]
  for (fit_id in rows$fit_id) s6l_run_fit(fit_id, output_root, TRUE)
  selection <- s6l_freeze_selection(rows, output_root, TRUE)
  pp_assert(nrow(selection$rows) == 2L && nrow(selection$winners) == 2L,
            "Stage 6L smoke result is incomplete.")
  pp_write_lines(c(
    "status=MICRO_SMOKE_PASS", "generated_arms=2", "fits=2", "winners=2",
    "truth_unseal_authorized=FALSE", "formal_v0lv_result=FALSE",
    paste0("completed_utc=", pp_iso_time())),
    file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
}

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  verification <- s6l_verify_environment()
  cat("CHECK_PASS source_rows=", nrow(verification$binding),
      " data=2 fits=8 workers=", PP_WORKERS, "\n", sep = "")
} else if (action == "smoke") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  s6l_smoke(output_root)
} else if (action == "generate") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  s6l_generate(output_root, FALSE)
} else if (action == "fit") {
  pp_assert(nzchar(output_root) && nzchar(cli$`fit-id` %||% ""),
            "--output-root and --fit-id are required.")
  s6l_run_fit(cli$`fit-id`, output_root, pp_bool(cli$smoke, FALSE))
} else if (action == "full") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  s6l_full(output_root)
} else if (action == "status") {
  pp_assert(nzchar(output_root), "--output-root is required.")
  complete <- length(list.files(file.path(output_root, "fits"),
                                "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  errors <- length(list.files(file.path(output_root, "fits"),
                              "^FIT_ERROR[.]txt$", recursive = TRUE))
  cat("complete_fits=", complete, "\nerror_fits=", errors,
      "\ntruth_free_selection_frozen=",
      file.exists(file.path(output_root, "truth_free_selection",
                            "TRUTH_FREE_SELECTION_FROZEN.txt")),
      "\n", sep = "")
} else {
  pp_abort("Unknown action: ", action)
}
