#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
# Keep an existing relative ASCII path on Windows. Some R installations start
# with a broken C.UTF-8 locale and corrupt non-ASCII absolute workspace paths
# during normalizePath(), while the same relative path remains fully usable.
root <- dirname(script_file)
if (!dir.exists(root)) {
  root <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
}
experiment_id <- "G12_PAPER_PIPELINE_PILOT_V1_20260829"

write_csv_lf <- function(value, path) {
  lines <- capture.output(utils::write.csv(
    value, row.names = FALSE, quote = TRUE, na = ""
  ))
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, sep = "\n", useBytes = TRUE)
  invisible(path)
}

data_manifest <- data.frame(
  experiment_id = experiment_id,
  data_id = c("g12pp_base_01", "g12pp_weak_01"),
  data_index = 1:2,
  data_seed = 82929001:82929002,
  scenario_id = c("baseline_strong", "weak_loading_separation"),
  scenario_replicate = 1L,
  n_obs_min = 6L,
  n_obs_max = 9L,
  target_shared_specific_abs_cosine = c(0, 0.6),
  new_unseen_at_registration = TRUE,
  truth_available_to_fit = FALSE,
  truth_available_to_stopping = FALSE,
  truth_available_to_selection = FALSE,
  truth_unseal_authorized = FALSE,
  engineering_pilot = TRUE,
  formal_v0lv_result = FALSE,
  formal_paper_mc_result = FALSE,
  stringsAsFactors = FALSE
)

g12 <- do.call(rbind, lapply(seq_len(nrow(data_manifest)), function(i) {
  data <- data_manifest[i, , drop = FALSE]
  seed_index <- 1:12
  data.frame(
    experiment_id = experiment_id,
    fit_id = sprintf("%s__G12__%02d", data$data_id, seed_index),
    data_id = data$data_id,
    data_index = data$data_index,
    data_seed = data$data_seed,
    scenario_id = data$scenario_id,
    method_id = "G12",
    route_id =
      "random__gram_unit_energy__pre1__jaoua_default__multistart12__fixed400",
    seed_index = seed_index,
    fit_seed = 88290000L + 100L * data$data_index + seed_index,
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
    actual_racing = FALSE,
    continuation = FALSE,
    automatic_800 = FALSE,
    truth_available_to_fit = FALSE,
    truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    engineering_pilot = TRUE,
    formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
}))

pooled <- do.call(rbind, lapply(seq_len(nrow(data_manifest)), function(i) {
  data <- data_manifest[i, , drop = FALSE]
  data.frame(
    experiment_id = experiment_id,
    fit_id = paste0(data$data_id, "__pooled_bayesSYNC__01"),
    data_id = data$data_id,
    data_index = data$data_index,
    data_seed = data$data_seed,
    scenario_id = data$scenario_id,
    method_id = "pooled_bayesSYNC",
    route_id = "pooled_two_studies_fit_once_Q3_L2",
    seed_index = 1L,
    fit_seed = 88390000L + data$data_index,
    initialization = "reference_random",
    function_initialization = "reference_default",
    pre_score_sweeps = 0L,
    anneal = "c(1,1.9,100)",
    planned_annealing_sweeps = 99L,
    planned_ordinary_t1_sweeps = 200L,
    total_maxit = 299L,
    n_cpus = 1L,
    objective_eligibility_required = TRUE,
    truth_free_winner_candidate = FALSE,
    cross_model_elbo_comparison_forbidden = TRUE,
    actual_racing = FALSE,
    continuation = FALSE,
    automatic_800 = FALSE,
    truth_available_to_fit = FALSE,
    truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    engineering_pilot = TRUE,
    formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
}))

fit_manifest <- rbind(g12, pooled)
fit_manifest <- fit_manifest[order(fit_manifest$data_index,
                                   fit_manifest$method_id,
                                   fit_manifest$seed_index), , drop = FALSE]
rownames(fit_manifest) <- NULL

metric_manifest <- data.frame(
  priority = 1:8,
  metric_group = c(
    "reconstruction", "loading_first_identity", "features",
    "loadings", "scores_and_processes", "contribution_and_covariance",
    "runtime_objective_stability", "auxiliary_dimensions"
  ),
  required_outputs = c(
    "observed_signal_nrmse;dense_signal_nrmse;observed_dense_gap",
    "matched;missing;misplaced;duplicate;extra;R;P;L",
    "total_ise;projection_floor_ise;excess_ise;matched_total_coverage",
    "direction;support;raw_scale;canonical_scale",
    "fpca_score_correlation;factor_process_error",
    "complete_contribution_error;factor_kernel;covariance_operator",
    "ordinary_elbo;elapsed;peak_memory;warning;objective;practical;380_to_400",
    "selected_factor_count;retained_M"
  ),
  selection_input = FALSE,
  new_success_threshold = FALSE,
  stringsAsFactors = FALSE
)

qc <- data.frame(
  check_id = c(
    "two_new_data_rows", "twenty_six_fit_rows", "method_counts",
    "unique_ids_and_seeds", "g12_twelve_per_data", "pooled_once_per_data",
    "all_single_core", "all_truth_flags_false", "no_racing_or_continuation",
    "no_formal_result_claim"
  ),
  passed = c(
    nrow(data_manifest) == 2L && all(data_manifest$new_unseen_at_registration),
    nrow(fit_manifest) == 26L,
    identical(as.integer(table(fit_manifest$method_id)[c("G12", "pooled_bayesSYNC")]),
              c(24L, 2L)),
    !anyDuplicated(fit_manifest$fit_id) && !anyDuplicated(fit_manifest$fit_seed),
    all(table(g12$data_id) == 12L),
    all(table(pooled$data_id) == 1L),
    all(fit_manifest$n_cpus == 1L),
    !any(data_manifest$truth_available_to_fit) &&
      !any(data_manifest$truth_available_to_stopping) &&
      !any(data_manifest$truth_available_to_selection) &&
      !any(fit_manifest$truth_available_to_fit) &&
      !any(fit_manifest$truth_available_to_stopping) &&
      !any(fit_manifest$truth_available_to_selection),
    !any(fit_manifest$actual_racing) && !any(fit_manifest$continuation) &&
      !any(fit_manifest$automatic_800),
    !any(fit_manifest$formal_v0lv_result) &&
      !any(fit_manifest$formal_paper_mc_result)
  ),
  detail = c(
    paste(data_manifest$data_seed, collapse = ";"),
    as.character(nrow(fit_manifest)),
    paste(names(table(fit_manifest$method_id)),
          as.integer(table(fit_manifest$method_id)), collapse = ";"),
    "fit_id and fit_seed are unique",
    "12 G12 starts for each data set",
    "both studies pooled once for each data set",
    "per-fit n_cpus=1",
    "truth unavailable to fit/stopping/selection",
    "racing=FALSE;continuation=FALSE;automatic_800=FALSE",
    "engineering pilot only"
  ),
  stringsAsFactors = FALSE
)

if (!all(qc$passed)) {
  stop("Manifest QC failed: ", paste(qc$check_id[!qc$passed], collapse = ", "))
}

write_csv_lf(data_manifest, file.path(root, "DATA_MANIFEST.csv"))
write_csv_lf(fit_manifest, file.path(root, "FIT_MANIFEST.csv"))
write_csv_lf(metric_manifest, file.path(root, "METRIC_MANIFEST.csv"))
write_csv_lf(qc, file.path(root, "MANIFEST_QC.csv"))
cat("MANIFESTS_COMPLETE data=2 fits=26 G12=24 pooled=2 qc=10/10\n")
