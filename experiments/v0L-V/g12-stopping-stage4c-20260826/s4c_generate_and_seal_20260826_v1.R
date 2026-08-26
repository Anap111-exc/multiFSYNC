#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

s4c_verify_environment(load_runtime = TRUE, package_role = "frozen_generation")
manifest <- uc_read_csv(file.path(UC_ROOT, "DATA_MANIFEST.csv"))
s4c_validate_data_manifest(manifest)

calibration_env <- new.env(parent = globalenv())
sys.source(uc_find_one(UC_SNAPSHOT, "^v0f_data_calibration[.]R$"),
           envir = calibration_env)
v0f_rebuild_observations <- calibration_env$v0f_rebuild_observations
calibrate_loading_information_v0f <-
  calibration_env$calibrate_loading_information_v0f
generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")

data_root <- file.path(UC_ROOT, "data")
dir.create(data_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")
status_rows <- vector("list", nrow(manifest))

for (index in seq_len(nrow(manifest))) {
  row <- manifest[index, , drop = FALSE]
  output_dir <- file.path(data_root, row$data_id[[1L]])
  if (dir.exists(output_dir)) {
    complete <- file.path(output_dir, "DATA_SEAL_COMPLETE.txt")
    observation_path <- file.path(output_dir, "observation_bundle.rds")
    truth_path <- file.path(output_dir, "sealed_truth", "truth_bundle.rds")
    uc_assert(file.exists(complete) && file.exists(observation_path) &&
                file.exists(truth_path),
              paste0("Incomplete pre-existing data directory: ",
                     row$data_id[[1L]]))
    observation <- readRDS(observation_path)
    uc_assert(identical(observation$data_id, row$data_id[[1L]]) &&
                identical(as.integer(observation$data_seed),
                          as.integer(row$data_seed[[1L]])) &&
                identical(observation$scenario_id, row$scenario_id[[1L]]) &&
                !length(candidate2_forbidden_observation_names(observation)),
              "Pre-existing observation bundle failed isolation.")
    reused <- TRUE
  } else {
    config <- s4c_base_config(row, smoke = FALSE)
    call <- config
    call$seed <- as.integer(row$data_seed[[1L]])
    captured <- c2_capture_conditions(
      do.call(generator, call), phase = "data_generation"
    )
    if (!is.null(captured$error)) {
      stop("Data generation failed: ", captured$error$message)
    }
    data <- calibrate_loading_information_v0f(
      captured$value,
      regime = "constant_per_active_loading",
      anchor_p = 100L, active_fraction = 0.1, anchor_norm = 3
    )
    if (row$scenario_id[[1L]] == "weak_loading_separation") {
      data <- s4c_apply_weak_loading_transform(data, rho = 0.6)
    }
    assertions <- s4c_truth_assertions(data, row, config)
    if (!all(assertions$passed)) {
      stop("Truth assertions failed for ", row$data_id[[1L]], ": ",
           paste(assertions$assertion_id[!assertions$passed], collapse = ", "))
    }
    observation <- list(
      bundle_class = "v0lv_observation_only_candidate2",
      data_id = row$data_id[[1L]],
      data_seed = as.integer(row$data_seed[[1L]]),
      scenario_id = row$scenario_id[[1L]],
      stopping_semantics_development = TRUE,
      formal_experiment = FALSE,
      Y = data$Y, time_obs = data$time_obs, Z = data$Z,
      dimensions = config[c(
        "S", "n_s", "p", "d", "L_f", "L_s", "M_f", "M_s", "K"
      )]
    )
    if (length(candidate2_forbidden_observation_names(observation))) {
      stop("Observation bundle leaks truth fields.")
    }
    sealed <- list(
      bundle_class = "v0lv_sealed_truth_candidate2",
      data_id = row$data_id[[1L]],
      data_seed = as.integer(row$data_seed[[1L]]),
      scenario_id = row$scenario_id[[1L]],
      stopping_semantics_development = TRUE,
      formal_experiment = FALSE,
      data = data,
      sealed_generation_provenance = list(
        generator_config = config,
        weak_loading_target_cosine =
          row$target_shared_specific_abs_cosine[[1L]],
        generation_warnings = captured$warnings,
        truth_assertions = assertions,
        generated_at_utc = captured$ended_at_utc
      )
    )
    dir.create(output_dir, recursive = TRUE, mode = "0700")
    observation_write <- c2_atomic_save_rds(
      observation, file.path(output_dir, "observation_bundle.rds")
    )
    sealed_dir <- file.path(output_dir, "sealed_truth")
    dir.create(sealed_dir, mode = "0700")
    truth_write <- c2_atomic_save_rds(
      sealed, file.path(sealed_dir, "truth_bundle.rds")
    )
    c2_atomic_write_csv(
      assertions, file.path(sealed_dir, "truth_assertions.csv")
    )
    c2_atomic_write_lines(c(
      "G12SS4C_DATA_SEAL_COMPLETE",
      paste0("data_id=", row$data_id[[1L]]),
      paste0("scenario_id=", row$scenario_id[[1L]]),
      "formal_v0lv_result=FALSE",
      paste0("observation_sha256=", observation_write$sha256),
      paste0("sealed_truth_sha256=", truth_write$sha256),
      "truth_unseal_authorized=FALSE"
    ), file.path(output_dir, "DATA_SEAL_COMPLETE.txt"))
    reused <- FALSE
  }
  status_rows[[index]] <- data.frame(
    data_id = row$data_id[[1L]], scenario_id = row$scenario_id[[1L]],
    data_seed = row$data_seed[[1L]], observation_created = TRUE,
    sealed_truth_created = TRUE, truth_unseal_authorized = FALSE,
    truth_available_to_fit_stopping_or_analysis = FALSE,
    reused_complete_bundle = reused,
    generated_or_verified_utc = uc_iso_time(),
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )
  utils::write.csv(
    do.call(rbind, status_rows[seq_len(index)]),
    file.path(UC_ROOT, "DATA_GENERATION_STATUS.csv"),
    row.names = FALSE, quote = TRUE, na = ""
  )
}

status <- do.call(rbind, status_rows)
uc_assert(nrow(status) == 6L && all(status$observation_created) &&
            all(status$sealed_truth_created) &&
            !any(status$truth_unseal_authorized) &&
            !any(status$truth_available_to_fit_stopping_or_analysis),
          "Data generation/sealing completion contract failed.")
writeLines(c(
  "status=DATA_GENERATION_AND_SEAL_COMPLETE",
  paste0("completed_utc=", uc_iso_time()),
  "data_sets=6", "scenarios=3", "replicates_per_scenario=2",
  "observation_bundles=6", "sealed_truth_bundles=6",
  "truth_unseal_authorized=FALSE",
  "truth_available_to_fit_stopping_or_analysis=FALSE",
  "formal_v0lv_result=FALSE"
), file.path(UC_ROOT, "DATA_GENERATION_COMPLETE.txt"), useBytes = TRUE)
cat("DATA_GENERATION_AND_SEAL_COMPLETE data_sets=6 scenarios=3 truth=0\n")
