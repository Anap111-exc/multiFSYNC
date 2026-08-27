#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
if (!grepl("^([A-Za-z]:|/)", script_file)) {
  script_file <- file.path(getwd(), script_file)
}
root <- dirname(script_file)

experiment_id <- "G12_FOCUSED_REACHABILITY_STAGE4E_V1_20260827"
data_id <- c("g12ss4c_base_01", "g12ss4c_sparse_01", "g12ss4c_sparse_02")
source_index <- c(1L, 5L, 6L)
data_seed <- c(82726001L, 82726005L, 82726006L)
scenario_id <- c("baseline_strong", "sparse_observation",
                 "sparse_observation")
observation_sha256 <- c(
  "f1a183dd568ad0838f8c037c5de526304550f9c9feef72d0f35f918670ca8257",
  "bcfe2f856d9e8c41c6b6810dd52abbc76bb1e3664d9701178138beac588f131b",
  "895afc709321aad5f8916316c0569c69d86970ca0013f31009d791c41e858c81"
)
target_reason <- c(
  "stage4d_two_start_endpoints_had_correct_counts_but_wrong_shared_specific_identity",
  rep("stage4d_two_start_endpoints_both_missed_study1_specific_factor", 2L)
)

data_binding <- data.frame(
  experiment_id = experiment_id, data_id = data_id,
  source_experiment_id = "G12_STOPPING_SEMANTICS_STAGE4C_V1_20260826",
  source_data_index = source_index, data_seed = data_seed,
  scenario_id = scenario_id, observation_sha256 = observation_sha256,
  adaptive_target_reason = target_reason,
  truth_available_to_fit = FALSE, truth_available_to_stopping = FALSE,
  truth_available_to_selection = FALSE,
  truth_previously_unsealed_stage4d = TRUE,
  development_only = TRUE, formal_v0lv_result = FALSE,
  stringsAsFactors = FALSE
)

fit_rows <- do.call(rbind, lapply(seq_along(data_id), function(index) {
  seed_index <- 3:8
  data.frame(
    experiment_id = experiment_id,
    fit_id = sprintf("%s__G12__%02d", data_id[[index]], seed_index),
    data_id = data_id[[index]], source_data_index = source_index[[index]],
    data_seed = data_seed[[index]], scenario_id = scenario_id[[index]],
    method_id = "G12", function_initialization = "gram_unit_energy",
    calibration_mode = "function", seed_index = seed_index,
    fit_seed = 87266000L + 10L * source_index[[index]] + seed_index,
    n_cpus = 1L, initialization = "random", pre_score_sweeps = 1L,
    anneal = "c(1,1.9,100)", planned_annealing_sweeps = 99L,
    fixed_ordinary_T1_sweeps = 400L, total_maxit = 499L,
    stopping_profile = "g12_stopping_control", stopping_min_t1 = 396L,
    stopping_max_t1 = 400L, stopping_gate = "quantile_factor",
    stopping_consecutive = 5L,
    checkpoint_sweeps = paste(seq.int(20L, 400L, by = 20L), collapse = ";"),
    fixed_400_endpoint = TRUE, objective_eligibility_required = TRUE,
    outer_parallel_fit_limit = 4L, actual_racing = FALSE,
    continuation = FALSE, automatic_800 = FALSE,
    truth_available_to_fit = FALSE, truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    adaptive_post_truth_design = TRUE, development_only = TRUE,
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )
}))
rownames(fit_rows) <- NULL

old_rows <- do.call(rbind, lapply(seq_along(data_id), function(index) {
  seed_index <- 1:2
  fit_id <- sprintf("%s__G12__%02d", data_id[[index]], seed_index)
  data.frame(
    experiment_id = experiment_id, fit_id = fit_id,
    data_id = data_id[[index]], seed_index = seed_index,
    fit_seed = 87266000L + 10L * source_index[[index]] + seed_index,
    source_experiment_id = "G12_STOPPING_SEMANTICS_STAGE4C_V1_20260826",
    source_terminal_relative_path = file.path(
      "fits", fit_id, "terminal_record.rds"
    ),
    source_fit_relative_path = file.path("fits", fit_id, "fit.rds"),
    objective_eligibility_required = TRUE,
    truth_available_to_selection = FALSE, development_only = TRUE,
    stringsAsFactors = FALSE
  )
}))
rownames(old_rows) <- NULL

parent_binding_path <- file.path(
  root, "..", "g12-stopping-stage4c-20260826", "SOURCE_BINDING.csv"
)
source_binding <- utils::read.csv(
  parent_binding_path, stringsAsFactors = FALSE, check.names = FALSE
)
source_binding$experiment_id <- experiment_id
source_binding$runner_commit_recorded_at_deployment <- TRUE
source_binding$formal_v0lv_result <- FALSE

write_csv <- function(value, name) {
  utils::write.csv(value, file.path(root, name), row.names = FALSE,
                   quote = TRUE, na = "")
}
write_csv(data_binding, "TARGET_DATA_BINDING.csv")
write_csv(fit_rows, "FIT_MANIFEST.csv")
write_csv(old_rows, "OLD_ENDPOINT_MANIFEST.csv")
write_csv(source_binding, "SOURCE_BINDING.csv")

checks <- c(
  three_target_data = nrow(data_binding) == 3L,
  eighteen_new_fits = nrow(fit_rows) == 18L,
  six_new_starts_per_data = all(table(fit_rows$data_id) == 6L),
  six_inherited_endpoints = nrow(old_rows) == 6L,
  eight_total_starts_per_data = all(
    table(c(fit_rows$data_id, old_rows$data_id)) == 8L
  ),
  seeds_unique = !anyDuplicated(c(fit_rows$fit_seed, old_rows$fit_seed)),
  per_fit_single_cpu = all(fit_rows$n_cpus == 1L),
  fixed400_without_continuation = all(fit_rows$fixed_400_endpoint) &&
    !any(fit_rows$continuation) && !any(fit_rows$automatic_800),
  truth_free_fit_stop_selection =
    !any(fit_rows$truth_available_to_fit) &&
    !any(fit_rows$truth_available_to_stopping) &&
    !any(fit_rows$truth_available_to_selection),
  adaptive_development_not_formal = all(fit_rows$adaptive_post_truth_design) &&
    all(fit_rows$development_only) && !any(fit_rows$formal_v0lv_result)
)
write_csv(data.frame(
  check_id = names(checks), passed = unname(checks),
  stringsAsFactors = FALSE
), "MANIFEST_QC.csv")
if (!all(checks)) stop("Stage-4E manifest QC failed.", call. = FALSE)
cat("STAGE4E_MANIFEST_PASS data=3 new_fits=18 inherited=6 total=24\n")
