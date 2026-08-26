#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale(
    "LC_CTYPE", "Chinese (Simplified)_China.utf8"
  )))
}

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
))
experiment_id <- "G12_STOPPING_SEMANTICS_STAGE4C_V1_20260826"
registered_utc <- "2026-08-26T08:22:45Z"

scenarios <- data.frame(
  scenario_id = c(
    "baseline_strong", "weak_loading_separation", "sparse_observation"
  ),
  scenario_label = c(
    "strong-identification independent replication",
    "full shared-specific support overlap with loading cosine 0.6",
    "strong-identification with 3-5 irregular observations"
  ),
  n_obs_min = c(6L, 6L, 3L),
  n_obs_max = c(9L, 9L, 5L),
  target_shared_specific_abs_cosine = c(0, 0.6, 0),
  shared_specific_support_overlap = c(0, 1, 0),
  otherwise_frozen_generator_configuration = TRUE,
  stringsAsFactors = FALSE
)
utils::write.csv(scenarios, file.path(root, "SCENARIOS.csv"),
                 row.names = FALSE, quote = TRUE, na = "")

data_ids <- c(
  "g12ss4c_base_01", "g12ss4c_base_02",
  "g12ss4c_weak_01", "g12ss4c_weak_02",
  "g12ss4c_sparse_01", "g12ss4c_sparse_02"
)
scenario_ids <- rep(scenarios$scenario_id, each = 2L)
scenario_match <- match(scenario_ids, scenarios$scenario_id)
data_manifest <- data.frame(
  experiment_id = experiment_id,
  data_id = data_ids,
  data_index = seq_along(data_ids),
  scenario_id = scenario_ids,
  scenario_replicate = rep(1:2, times = 3L),
  data_seed = 82726001:82726006,
  S = 2L, n_s = "30;30", p = 500L, d = 0L,
  L_f = 1L, L_s = "1;1", M_f = 2L, M_s = "2;2", K = 5L,
  n_obs_min = scenarios$n_obs_min[scenario_match],
  n_obs_max = scenarios$n_obs_max[scenario_match],
  sigma_eps = 0.3, mean_amp = 0.6,
  target_shared_specific_abs_cosine =
    scenarios$target_shared_specific_abs_cosine[scenario_match],
  shared_specific_support_overlap =
    scenarios$shared_specific_support_overlap[scenario_match],
  new_unseen_at_registration = TRUE,
  truth_available_to_fit = FALSE,
  truth_available_to_stopping = FALSE,
  truth_available_to_analysis = FALSE,
  truth_unseal_authorized = FALSE,
  development_only = TRUE,
  formal_v0lv_result = FALSE,
  registered_utc = registered_utc,
  stringsAsFactors = FALSE
)
utils::write.csv(data_manifest, file.path(root, "DATA_MANIFEST.csv"),
                 row.names = FALSE, quote = TRUE, na = "")

fit_rows <- list()
index <- 0L
for (data_index in seq_len(nrow(data_manifest))) {
  for (seed_index in 1:2) {
    index <- index + 1L
    fit_seed <- 87266000L + 10L * data_index + seed_index
    fit_rows[[index]] <- data.frame(
      experiment_id = experiment_id,
      fit_id = sprintf("%s__G12__%02d",
                       data_manifest$data_id[[data_index]], seed_index),
      data_id = data_manifest$data_id[[data_index]],
      data_index = data_index,
      data_seed = data_manifest$data_seed[[data_index]],
      scenario_id = data_manifest$scenario_id[[data_index]],
      method_id = "G12",
      function_initialization = "gram_unit_energy",
      calibration_mode = "function",
      seed_index = seed_index,
      fit_seed = fit_seed,
      n_cpus = 1L,
      initialization = "random",
      pre_score_sweeps = 1L,
      anneal = "c(1,1.9,100)",
      planned_annealing_sweeps = 99L,
      fixed_ordinary_T1_sweeps = 400L,
      total_maxit = 499L,
      stopping_profile = "g12_stopping_control",
      stopping_min_t1 = 396L,
      stopping_max_t1 = 400L,
      stopping_gate = "quantile_factor",
      stopping_consecutive = 5L,
      checkpoint_sweeps = paste(seq(20L, 400L, 20L), collapse = ";"),
      fixed_400_endpoint = TRUE,
      objective_eligibility_required = TRUE,
      outer_parallel_fit_limit = 4L,
      actual_racing = FALSE,
      continuation = FALSE,
      automatic_800 = FALSE,
      truth_available_to_fit = FALSE,
      truth_available_to_stopping = FALSE,
      truth_available_to_analysis = FALSE,
      development_only = TRUE,
      formal_v0lv_result = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
fit_manifest <- do.call(rbind, fit_rows)
utils::write.csv(fit_manifest, file.path(root, "FIT_MANIFEST.csv"),
                 row.names = FALSE, quote = TRUE, na = "")

source_binding <- data.frame(
  experiment_id = experiment_id,
  dev_branch = "dev/v0lv-init-geometry-20260818",
  dev_commit = "72d9a53f0ef5e9d3e5f49d9cf837958207b3c980",
  package_code_commit = "72d9a53f0ef5e9d3e5f49d9cf837958207b3c980",
  runner_commit_recorded_at_deployment = TRUE,
  parent_frozen_tag = "v0.3.0-v0L-V-frozen1",
  parent_frozen_commit = "aee98b79a80a6535a0ad7a7c5187b29db2f66176",
  package_version = "0.3.0.9000",
  package_source_sha256 =
    "2e48b3d936196b9e2ddc3146afba6f8cc347d8797bcc5ef1ca340a3a936e9af5",
  public_interface_sha256 =
    "b8e7d74d318aabfc7e936089b75f64ed2e4c5c6182ecd0df06af984c1359662f",
  generator_source_sha256 =
    "4e459c5a1bac12491afb6f4f63a74ca0c4551027fae1f5f037080ebe2a7d6aa4",
  calibration_source_sha256 =
    "7eec6166fac3678126944bfc6bae14824ce26a4a9829fcd0ed49aa1dd9d1ea93",
  scale_trace_sha256 =
    "2e6a6e02c23300f2f8151c68a6521cb29ec30a55f563c2b990d7b9867f432da3",
  multi_core_sha256 =
    "f06a14c853dd16c3695ef43447ef2fe341534070a54cd3d627b0608465816b52",
  practical_stopping_sha256 =
    "8cf9a5a119fbaa23dfeb3abb81a90ea31e6565977629327821e3538ed93bf5fc",
  initialization_sha256 =
    "b69f4f1537222aecfca02103a03dc1add3c156529b9f4bbce88d78bad0face1b",
  g12_stopping_control_sha256 =
    "90b3a2597547fdc926caa0045dcfa353396cdd48f6cc0cad81216dd795c20f34",
  description_sha256 =
    "41438e12379ca645d585c2b755eaeafa5d26677a4b1248d29575bd5915e213b8",
  namespace_sha256 =
    "fd33c1d9e9ae138e3ee1e5ee80dc67c3a9e261d576af09cd8ae19df5d3033ee3",
  upstream_common_sha256 =
    "5dff355a66f3900326e585a959800cb46a2c274705b1baa96cc211fb888c2e26",
  model_prior_cavi_elbo_changed = FALSE,
  public_stopping_default_changed = FALSE,
  formal_v0lv_result = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(source_binding, file.path(root, "SOURCE_BINDING.csv"),
                 row.names = FALSE, quote = TRUE, na = "")

qc <- data.frame(
  check_id = c(
    "six_new_data", "three_scenarios_two_replicates",
    "twelve_fixed_400_fits", "two_fit_seeds_per_data",
    "unique_ids_and_seeds", "per_fit_single_cpu", "G12_route_fixed",
    "public_fixed_400_profile", "checkpoint_capture_preregistered",
    "truth_isolated", "no_racing_continuation_or_800", "not_formal_v0lv"
  ),
  passed = c(
    nrow(data_manifest) == 6L && all(data_manifest$new_unseen_at_registration),
    all(table(data_manifest$scenario_id) == 2L),
    nrow(fit_manifest) == 12L && all(fit_manifest$fixed_ordinary_T1_sweeps == 400L),
    all(table(fit_manifest$data_id) == 2L),
    !anyDuplicated(data_manifest$data_id) &&
      !anyDuplicated(data_manifest$data_seed) &&
      !anyDuplicated(fit_manifest$fit_id) &&
      !anyDuplicated(fit_manifest$fit_seed),
    all(fit_manifest$n_cpus == 1L),
    all(fit_manifest$method_id == "G12") &&
      all(fit_manifest$function_initialization == "gram_unit_energy") &&
      all(fit_manifest$pre_score_sweeps == 1L) &&
      all(fit_manifest$planned_annealing_sweeps == 99L),
    all(fit_manifest$stopping_profile == "g12_stopping_control") &&
      all(fit_manifest$stopping_min_t1 == 396L) &&
      all(fit_manifest$stopping_max_t1 == 400L) &&
      all(fit_manifest$stopping_consecutive == 5L) &&
      all(fit_manifest$fixed_400_endpoint),
    all(fit_manifest$checkpoint_sweeps ==
          paste(seq(20L, 400L, 20L), collapse = ";")),
    !any(fit_manifest$truth_available_to_fit) &&
      !any(fit_manifest$truth_available_to_stopping) &&
      !any(fit_manifest$truth_available_to_analysis),
    !any(fit_manifest$actual_racing) && !any(fit_manifest$continuation) &&
      !any(fit_manifest$automatic_800),
    !any(data_manifest$formal_v0lv_result) &&
      !any(fit_manifest$formal_v0lv_result)
  ),
  detail = c(
    "6", "2 per baseline/weak/sparse", "12", "2", "unique", "1",
    "random + G12 + one pre-score + Jaoua anneal",
    "public min396/max400 quantile-factor five consecutive",
    "20-sweep checkpoints through fixed endpoint 400",
    "fit/stopping/analysis flags FALSE", "FALSE/FALSE/FALSE",
    "post-v0L-V development only"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(qc, file.path(root, "MANIFEST_QC.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
if (!all(qc$passed)) stop("Manifest QC failed.", call. = FALSE)
cat("MANIFEST_BUILD_PASS data=6 fits=12 fixed_t1=400 truth=0\n")
