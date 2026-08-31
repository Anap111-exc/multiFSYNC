options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_arg[[1L]])
root <- dirname(script_file)
if (!dir.exists(root)) {
  root <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
}

write_csv_lf <- function(value, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  text <- capture.output(utils::write.csv(value, row.names = FALSE,
                                          quote = TRUE, na = ""))
  writeBin(charToRaw(paste0(paste(text, collapse = "\n"), "\n")), con)
}

experiment_id <- "G12_DIMENSION_OVERSPEC_STAGE6A_V1_20260831"
data_ids <- "g12s6a_01"

data_manifest <- data.frame(
  experiment_id = experiment_id,
  data_id = data_ids,
  data_index = 1L,
  data_seed = 83131001L,
  scenario_id = "baseline_strong",
  scenario_replicate = 1L,
  S = 2L, n_s = "30;30", p = 500L, d = 0L,
  truth_L_f = 1L, truth_L_s = "1;1",
  truth_M_f = "2", truth_M_s = "2;2", K = 5L,
  n_obs_min = 6L, n_obs_max = 9L,
  sigma_eps = 0.3, mean_amp = 0.6,
  target_shared_specific_abs_cosine = 0,
  new_unseen_at_registration = TRUE,
  truth_available_to_fit = FALSE,
  truth_available_to_stopping = FALSE,
  truth_available_to_selection = FALSE,
  truth_unseal_authorized = FALSE,
  development_only = TRUE,
  engineering_pilot = TRUE,
  formal_v0lv_result = FALSE,
  formal_paper_mc_result = FALSE
)

configs <- data.frame(
  fit_config_index = 1:3,
  fit_config_id = c("truth_L1_M2", "factor_L3_M2", "fpca_L1_M4"),
  fit_L_f = c(1L, 3L, 1L),
  fit_L_s_1 = c(1L, 3L, 1L),
  fit_L_s_2 = c(1L, 3L, 1L),
  fit_M_f = c("2", "2;2;2", "4"),
  fit_M_s_1 = c("2", "2;2;2", "4"),
  fit_M_s_2 = c("2", "2;2;2", "4"),
  dimension_change = c("none", "factor_count_only", "fpca_rank_cap_only")
)

fit_rows <- list()
row_index <- 0L
for (data_index in seq_len(nrow(data_manifest))) {
  data <- data_manifest[data_index, , drop = FALSE]
  for (config_index in seq_len(nrow(configs))) {
    config <- configs[config_index, , drop = FALSE]
    for (seed_index in 1:12) {
      row_index <- row_index + 1L
      fit_id <- paste(data$data_id, config$fit_config_id,
                      sprintf("%02d", seed_index), sep = "__")
      fit_rows[[row_index]] <- data.frame(
        experiment_id = experiment_id,
        fit_id = fit_id,
        data_id = data$data_id,
        data_index = data$data_index,
        data_seed = data$data_seed,
        scenario_id = data$scenario_id,
        method_id = "G12",
        fit_config_id = config$fit_config_id,
        fit_config_index = config$fit_config_index,
        selection_stratum_id = paste(data$data_id, config$fit_config_id,
                                     sep = "__"),
        route_id = paste0(
          "random__gram_unit_energy__pre1__jaoua_default__multistart12__",
          "fixed400__", config$fit_config_id
        ),
        seed_index = seed_index,
        fit_seed = 89300000L + 1000L * data$data_index +
          100L * config$fit_config_index + seed_index,
        fit_L_f = config$fit_L_f,
        fit_L_s_1 = config$fit_L_s_1,
        fit_L_s_2 = config$fit_L_s_2,
        fit_M_f = config$fit_M_f,
        fit_M_s_1 = config$fit_M_s_1,
        fit_M_s_2 = config$fit_M_s_2,
        dimension_change = config$dimension_change,
        initialization = "random",
        function_initialization = "gram_unit_energy",
        pre_score_sweeps = 1L,
        anneal = "c(1,1.9,100)",
        planned_annealing_sweeps = 99L,
        planned_ordinary_t1_sweeps = 400L,
        total_maxit = 499L,
        n_cpus = 1L,
        objective_eligibility_required = TRUE,
        truth_free_winner_candidate = TRUE,
        cross_model_elbo_comparison_forbidden = TRUE,
        cross_dimension_elbo_comparison_forbidden = TRUE,
        actual_racing = FALSE,
        continuation = FALSE,
        automatic_800 = FALSE,
        capacity_probe = data$data_index == 1L && seed_index == 1L,
        truth_available_to_fit = FALSE,
        truth_available_to_stopping = FALSE,
        truth_available_to_selection = FALSE,
        development_only = TRUE,
        engineering_pilot = TRUE,
        formal_v0lv_result = FALSE,
        formal_paper_mc_result = FALSE
      )
    }
  }
}
fit_manifest <- do.call(rbind, fit_rows)

metric_manifest <- data.frame(
  priority = 1:7,
  metric_group = c(
    "loading_first_factor_identity",
    "extra_factor_activity",
    "fpca_effective_rank",
    "reconstruction",
    "scientific_quality",
    "optimization_support",
    "runtime_resource"
  ),
  required_outputs = c(
    "matched;missing;misplaced;duplicate;extra;R;P;L",
    "factor_ppi;sum_loading_ppi;loading_energy;process_energy;contribution_energy",
    "effective_M;component_eigenvalue_share;true_subspace_recall;kernel_error",
    "observed_signal_nrmse;dense_signal_nrmse;observed_dense_gap",
    "feature_ise;loading_error;score_correlation;process_error;contribution_error;covariance_error",
    "eligible_ordinary_elbo;within_stratum_winner;cross_start_identity_rate",
    "elapsed;peak_memory;warning;practical_status;380_to_400"
  ),
  selection_input = FALSE,
  new_success_threshold = FALSE
)

strata <- interaction(fit_manifest$data_id, fit_manifest$fit_config_id,
                      drop = TRUE)
checks <- c(
  one_new_data_row = nrow(data_manifest) == 1L,
  thirty_six_fit_rows = nrow(fit_manifest) == 36L,
  three_registered_fit_configs =
    setequal(fit_manifest$fit_config_id, configs$fit_config_id),
  twelve_starts_per_stratum = all(table(strata) == 12L),
  three_truth_free_winners_planned = length(levels(strata)) == 3L,
  unique_fit_ids = !anyDuplicated(fit_manifest$fit_id),
  unique_fit_seeds = !anyDuplicated(fit_manifest$fit_seed),
  all_single_core = all(fit_manifest$n_cpus == 1L),
  fixed_g12_route = all(fit_manifest$pre_score_sweeps == 1L) &&
    all(fit_manifest$planned_annealing_sweeps == 99L) &&
    all(fit_manifest$planned_ordinary_t1_sweeps == 400L),
  no_cross_dimension_elbo_selection =
    all(fit_manifest$cross_dimension_elbo_comparison_forbidden),
  no_racing_or_continuation = !any(fit_manifest$actual_racing) &&
    !any(fit_manifest$continuation) && !any(fit_manifest$automatic_800),
  all_truth_flags_false = !any(fit_manifest$truth_available_to_fit) &&
    !any(fit_manifest$truth_available_to_stopping) &&
    !any(fit_manifest$truth_available_to_selection),
  three_capacity_probe_rows = sum(fit_manifest$capacity_probe) == 3L,
  no_formal_result_claim = !any(fit_manifest$formal_paper_mc_result)
)
qc <- data.frame(
  check_id = names(checks),
  passed = unname(checks),
  detail = c(
    paste(data_manifest$data_seed, collapse = ";"),
    nrow(fit_manifest),
    paste(configs$fit_config_id, collapse = ";"),
    paste(as.integer(table(strata)), collapse = ";"),
    length(levels(strata)),
    "fit_id unique", "fit_seed unique", "n_cpus=1",
    "G12;pre1;Jaoua99;fixed400",
    "winner selected within data_id and fit_config_id only",
    "racing=FALSE;continuation=FALSE;automatic_800=FALSE",
    "fit/stopping/selection truth flags FALSE",
    paste(fit_manifest$fit_id[fit_manifest$capacity_probe], collapse = ";"),
    "development pilot only"
  )
)
if (!all(qc$passed)) stop("Stage 6A manifest QC failed.", call. = FALSE)

write_csv_lf(data_manifest, file.path(root, "DATA_MANIFEST.csv"))
write_csv_lf(configs, file.path(root, "FIT_CONFIGS.csv"))
write_csv_lf(fit_manifest, file.path(root, "FIT_MANIFEST.csv"))
write_csv_lf(metric_manifest, file.path(root, "METRIC_MANIFEST.csv"))
write_csv_lf(qc, file.path(root, "MANIFEST_QC.csv"))
cat("STAGE6A_MANIFESTS_PASS data=1 configs=3 fits=36 strata=3 qc=14/14\n")
