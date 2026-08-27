argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
s4e_script_path <- normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
)

source(file.path(dirname(s4e_script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

UC_EXPERIMENT_ID <- "G12_FOCUSED_REACHABILITY_STAGE4E_V1_20260827"
UC_DATA_IDS <- c(
  "g12ss4c_base_01", "g12ss4c_sparse_01", "g12ss4c_sparse_02"
)
UC_SCENARIO_IDS <- c("baseline_strong", "sparse_observation")
UC_EXPECTED_FITS <- 18L
UC_OUTER_WORKERS <- 4L
S4E_NEW_SEED_INDICES <- 3:8
S4E_PARENT_STAGE4C_ROOT <- uc_env(
  "G12SS4E_PARENT_STAGE4C_ROOT",
  "/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1"
)

s4e_validate_data_binding <- function(binding) {
  required <- c(
    "experiment_id", "data_id", "source_experiment_id",
    "source_data_index", "data_seed", "scenario_id",
    "observation_sha256", "adaptive_target_reason",
    "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_selection", "truth_previously_unsealed_stage4d",
    "development_only", "formal_v0lv_result"
  )
  expected_hash <- c(
    g12ss4c_base_01 =
      "f1a183dd568ad0838f8c037c5de526304550f9c9feef72d0f35f918670ca8257",
    g12ss4c_sparse_01 =
      "bcfe2f856d9e8c41c6b6810dd52abbc76bb1e3664d9701178138beac588f131b",
    g12ss4c_sparse_02 =
      "895afc709321aad5f8916316c0569c69d86970ca0013f31009d791c41e858c81"
  )
  uc_assert(
    is.data.frame(binding) && nrow(binding) == 3L &&
      all(required %in% names(binding)) &&
      all(binding$experiment_id == UC_EXPERIMENT_ID) &&
      identical(as.character(binding$data_id), UC_DATA_IDS) &&
      identical(as.integer(binding$source_data_index), c(1L, 5L, 6L)) &&
      identical(as.integer(binding$data_seed), c(82726001L, 82726005L,
                                                  82726006L)) &&
      identical(as.character(binding$scenario_id),
                c("baseline_strong", rep("sparse_observation", 2L))) &&
      identical(as.character(binding$observation_sha256),
                unname(expected_hash[UC_DATA_IDS])) &&
      !any(binding$truth_available_to_fit) &&
      !any(binding$truth_available_to_stopping) &&
      !any(binding$truth_available_to_selection) &&
      all(binding$truth_previously_unsealed_stage4d) &&
      all(binding$development_only) && !any(binding$formal_v0lv_result),
    "Stage-4E target-data binding failed."
  )
  invisible(TRUE)
}

s4e_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "source_data_index",
    "data_seed", "scenario_id", "method_id", "function_initialization",
    "calibration_mode", "seed_index", "fit_seed", "n_cpus",
    "initialization", "pre_score_sweeps", "anneal",
    "planned_annealing_sweeps", "fixed_ordinary_T1_sweeps", "total_maxit",
    "stopping_profile", "stopping_min_t1", "stopping_max_t1",
    "stopping_gate", "stopping_consecutive", "checkpoint_sweeps",
    "fixed_400_endpoint", "objective_eligibility_required",
    "outer_parallel_fit_limit", "actual_racing", "continuation",
    "automatic_800", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "adaptive_post_truth_design", "development_only", "formal_v0lv_result"
  )
  expected_ids <- unlist(lapply(UC_DATA_IDS, function(data_id) {
    sprintf("%s__G12__%02d", data_id, S4E_NEW_SEED_INDICES)
  }), use.names = FALSE)
  source_index <- c(g12ss4c_base_01 = 1L, g12ss4c_sparse_01 = 5L,
                    g12ss4c_sparse_02 = 6L)
  expected_seed <- 87266000L +
    10L * unname(source_index[manifest$data_id]) +
    as.integer(manifest$seed_index)
  uc_assert(
    is.data.frame(manifest) && nrow(manifest) == UC_EXPECTED_FITS &&
      all(required %in% names(manifest)) &&
      all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
      identical(as.character(manifest$fit_id), expected_ids) &&
      !anyDuplicated(manifest$fit_id) && !anyDuplicated(manifest$fit_seed) &&
      all(table(manifest$data_id) == 6L) &&
      all(manifest$seed_index %in% S4E_NEW_SEED_INDICES) &&
      identical(as.integer(manifest$fit_seed), expected_seed) &&
      all(manifest$method_id == "G12") &&
      all(manifest$function_initialization == "gram_unit_energy") &&
      all(manifest$calibration_mode == "function") &&
      all(manifest$n_cpus == 1L) &&
      all(manifest$initialization == "random") &&
      all(manifest$pre_score_sweeps == 1L) &&
      all(manifest$anneal == "c(1,1.9,100)") &&
      all(manifest$planned_annealing_sweeps == 99L) &&
      all(manifest$fixed_ordinary_T1_sweeps == S4C_FIXED_T1) &&
      all(manifest$total_maxit == S4C_MAXIT) &&
      all(manifest$stopping_profile == "g12_stopping_control") &&
      all(manifest$stopping_min_t1 == S4C_STOPPING_MIN_T1) &&
      all(manifest$stopping_max_t1 == S4C_FIXED_T1) &&
      all(manifest$stopping_gate == "quantile_factor") &&
      all(manifest$stopping_consecutive == S4C_STOPPING_CONSECUTIVE) &&
      all(manifest$checkpoint_sweeps ==
            paste(S4C_CHECKPOINTS, collapse = ";")) &&
      all(manifest$fixed_400_endpoint) &&
      all(manifest$objective_eligibility_required) &&
      all(manifest$outer_parallel_fit_limit == UC_OUTER_WORKERS) &&
      !any(manifest$actual_racing) && !any(manifest$continuation) &&
      !any(manifest$automatic_800) &&
      !any(manifest$truth_available_to_fit) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_selection) &&
      all(manifest$adaptive_post_truth_design) &&
      all(manifest$development_only) && !any(manifest$formal_v0lv_result),
    "Stage-4E fit route, seed, budget, or truth-isolation contract failed."
  )
  invisible(TRUE)
}

s4e_validate_old_endpoint_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "seed_index", "fit_seed",
    "source_experiment_id", "source_terminal_relative_path",
    "source_fit_relative_path", "objective_eligibility_required",
    "truth_available_to_selection", "development_only"
  )
  expected_ids <- unlist(lapply(UC_DATA_IDS, function(data_id) {
    sprintf("%s__G12__%02d", data_id, 1:2)
  }), use.names = FALSE)
  uc_assert(
    is.data.frame(manifest) && nrow(manifest) == 6L &&
      all(required %in% names(manifest)) &&
      all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
      identical(as.character(manifest$fit_id), expected_ids) &&
      !anyDuplicated(manifest$fit_id) && all(table(manifest$data_id) == 2L) &&
      identical(as.integer(manifest$seed_index), rep(1:2, 3L)) &&
      all(manifest$objective_eligibility_required) &&
      !any(manifest$truth_available_to_selection) &&
      all(manifest$development_only),
    "Stage-4E inherited endpoint manifest failed."
  )
  invisible(TRUE)
}

s4e_verify_inputs <- function(load_runtime = TRUE) {
  s4c_verify_environment(load_runtime = load_runtime,
                         package_role = "development")
  data_binding <- uc_read_csv(file.path(UC_ROOT, "TARGET_DATA_BINDING.csv"))
  fit_manifest <- uc_read_csv(file.path(UC_ROOT, "FIT_MANIFEST.csv"))
  old_manifest <- uc_read_csv(file.path(UC_ROOT, "OLD_ENDPOINT_MANIFEST.csv"))
  s4e_validate_data_binding(data_binding)
  s4e_validate_fit_manifest(fit_manifest)
  s4e_validate_old_endpoint_manifest(old_manifest)
  uc_assert(dir.exists(S4E_PARENT_STAGE4C_ROOT),
            "The bound Stage-4C parent root is missing.")

  local_paths <- file.path(
    UC_ROOT, "data", data_binding$data_id, "observation_bundle.rds"
  )
  parent_paths <- file.path(
    S4E_PARENT_STAGE4C_ROOT, "data", data_binding$data_id,
    "observation_bundle.rds"
  )
  uc_assert(all(file.exists(local_paths)) && all(file.exists(parent_paths)),
            "A target observation-only bundle is missing.")
  local_hash <- vapply(local_paths, uc_sha256, character(1L))
  parent_hash <- vapply(parent_paths, uc_sha256, character(1L))
  uc_assert(identical(unname(local_hash),
                      as.character(data_binding$observation_sha256)) &&
              identical(unname(parent_hash),
                        as.character(data_binding$observation_sha256)),
            "Target observation-only bundle hash mismatch.")
  leaked <- list.files(file.path(UC_ROOT, "data"), recursive = TRUE,
                       full.names = FALSE)
  uc_assert(!any(grepl("sealed_truth|truth_bundle", leaked)),
            "Stage-4E fit root contains a truth bundle.")

  old_terminal <- file.path(
    S4E_PARENT_STAGE4C_ROOT, old_manifest$source_terminal_relative_path
  )
  old_fit <- file.path(
    S4E_PARENT_STAGE4C_ROOT, old_manifest$source_fit_relative_path
  )
  uc_assert(all(file.exists(old_terminal)) && all(file.exists(old_fit)),
            "An inherited Stage-4C endpoint is missing.")
  old_records <- lapply(old_terminal, readRDS)
  uc_assert(all(vapply(old_records, function(record) {
    identical(record$terminal_status, "fixed_400_complete") &&
      isTRUE(record$objective_eligible) &&
      identical(as.integer(record$actual_T1_sweeps), 400L) &&
      !isTRUE(record$truth_used) && !isTRUE(record$continuation_used) &&
      !isTRUE(record$automatic_800_used)
  }, logical(1L))), "An inherited Stage-4C endpoint is not reusable.")
  invisible(list(
    data_binding = data_binding, fit_manifest = fit_manifest,
    old_manifest = old_manifest, old_terminal_paths = old_terminal,
    old_fit_paths = old_fit
  ))
}

s4e_terminal_summary <- function(record, source_stage) {
  value <- s4c_terminal_summary(record)
  value$source_stage <- source_stage
  value
}
