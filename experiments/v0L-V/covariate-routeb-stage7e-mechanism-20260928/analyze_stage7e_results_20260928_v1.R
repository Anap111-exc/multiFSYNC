#!/usr/bin/env Rscript

# Read-only synthesis of completed Stage 7E diagnostics and frozen Stage 7C.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript SCRIPT STAGE7E_ROOT STAGE7C_ALL_START_SCIENCE_CSV")
}
root <- args[[1L]]
science_path <- args[[2L]]
if (!dir.exists(root) || !file.exists(science_path)) {
  stop("Input experiment or Stage 7C science table is missing.")
}
result_root <- file.path(root, "results", "derived_v1")
if (file.exists(result_root) || dir.exists(result_root)) {
  stop("Refusing to overwrite derived_v1.")
}
read <- function(...) utils::read.csv(file.path(...), stringsAsFactors = FALSE)
oracle_root <- file.path(root, "results", "oracle_summary_v1")
evaluation_root <- file.path(oracle_root, "evaluation_v2")
if (!file.exists(file.path(evaluation_root, "EVALUATION_COMPLETE.txt")) ||
    !file.exists(file.path(oracle_root, "FITS_COMPLETE.txt"))) {
  stop("Six-fit evaluation is incomplete.")
}
endpoints <- read(evaluation_root, "ORACLE_ENDPOINT_COMPARISON.csv")
factors <- read(evaluation_root, "ORACLE_FACTOR_DIAGNOSTICS.csv")
trace <- read(evaluation_root, "ORACLE_DRIVER_TRACE.csv")
scales <- read(root, "results", "ORACLE_B2_ENDPOINT_SCALES.csv")
stage7c <- utils::read.csv(science_path, stringsAsFactors = FALSE)
stage7c <- stage7c[stage7c$beta_arm == "C1" & stage7c$truth_free_winner, ]
expected_ids <- unlist(lapply(1:3, function(index) {
  sprintf("s7e_%02d_C1_%s", index, c("exact", "near"))
}), use.names = FALSE)
stopifnot(nrow(endpoints) == 6L,
          identical(endpoints$fit_id, expected_ids),
          nrow(factors) == 18L, nrow(scales) == 6L,
          nrow(stage7c) == 3L,
          all(endpoints$objective_eligible),
          all(endpoints$ordinary_sweeps == 400L),
          all(endpoints$annealing_sweeps == 0L),
          !any(endpoints$truth_used_for_stopping_or_selection),
          !any(endpoints$formal_candidate_pool))
comparison <- do.call(rbind, lapply(seq_len(nrow(endpoints)), function(i) {
  row <- endpoints[i, ]
  reference <- stage7c[stage7c$data_id == row$data_id, ]
  factor <- factors[factors$fit_id == row$fit_id, ]
  scale <- scales[scales$fit_id == row$fit_id, ]
  stopifnot(nrow(reference) == 1L, nrow(factor) == 3L,
            nrow(scale) == 1L,
            identical(as.character(row$stage7c_winner_fit_id),
                      as.character(reference$fit_id)),
            isTRUE(all.equal(row$stage7c_winner_ordinary_elbo,
                             reference$final_ordinary_elbo,
                             tolerance = 1e-9)),
            identical(row$stage7c_winner_core_truth_recovered,
                      reference$core_truth_recovered),
            identical(row$stage7c_winner_exact_structure,
                      reference$exact_thresholded_factor_structure))
  factor <- factor[match(c("shared", "specific_1", "specific_2"),
                         factor$block), ]
  data.frame(
    data_id = row$data_id, oracle_fit_id = row$fit_id, mode = row$mode,
    stage7c_winner_fit_id = reference$fit_id,
    stage7c_winner_elbo = reference$final_ordinary_elbo,
    oracle_elbo = row$final_ordinary_elbo,
    oracle_minus_winner_elbo = row$elbo_minus_stage7c_winner,
    stage7c_winner_core = reference$core_truth_recovered,
    oracle_core = row$frozen_core_truth_recovered,
    stage7c_winner_exact = reference$exact_thresholded_factor_structure,
    oracle_exact = row$frozen_exact_thresholded_structure,
    stage7c_winner_B2_selected =
      reference$selected_specific_study_2 == 1L,
    oracle_B2_selected = factor$factor_selected[[3L]],
    oracle_B2_true_positive_support =
      factor$support_true_positive[[3L]],
    oracle_B2_false_positive_support =
      factor$support_false_positive[[3L]],
    oracle_B2_loading_abs_cosine =
      factor$canonical_loading_abs_cosine[[3L]],
    oracle_B2_raw_loading_norm_ratio = scale$raw_loading_norm_ratio,
    oracle_B2_time_l2_ratio = scale$time_function_l2_ratio,
    oracle_B2_score_rms_ratio = scale$score_rms_ratio,
    stage7c_winner_strict = reference$strict_practical_converged,
    oracle_strict = row$strict_practical_converged,
    stringsAsFactors = FALSE)
}))
pair <- do.call(rbind, lapply(unique(comparison$data_id), function(data_id) {
  group <- comparison[comparison$data_id == data_id, ]
  exact <- group[group$mode == "exact", ]
  near <- group[group$mode == "near", ]
  stopifnot(nrow(exact) == 1L, nrow(near) == 1L)
  data.frame(
    data_id = data_id,
    exact_minus_near_elbo = exact$oracle_elbo - near$oracle_elbo,
    both_core = exact$oracle_core && near$oracle_core,
    both_exact = exact$oracle_exact && near$oracle_exact,
    B2_loading_cosine_absolute_difference =
      abs(exact$oracle_B2_loading_abs_cosine -
            near$oracle_B2_loading_abs_cosine),
    B2_support_count_difference =
      exact$oracle_B2_true_positive_support -
      near$oracle_B2_true_positive_support,
    stringsAsFactors = FALSE)
}))
elbo <- do.call(rbind, lapply(endpoints$fit_id, function(fit_id) {
  path <- file.path(oracle_root, "fits", fit_id, "ORDINARY_ELBO.csv")
  series <- read.csv(path)$ordinary_elbo
  stopifnot(length(series) == 400L, all(is.finite(series)))
  data.frame(
    fit_id = fit_id, elbo_sweep_1 = series[[1L]],
    elbo_sweep_20 = series[[20L]],
    elbo_sweep_100 = series[[100L]],
    elbo_sweep_380 = series[[380L]],
    elbo_sweep_400 = series[[400L]],
    elbo_gain_380_to_400 = series[[400L]] - series[[380L]],
    negative_ordinary_elbo_steps = sum(diff(series) < 0),
    minimum_ordinary_elbo_step = min(diff(series)),
    stringsAsFactors = FALSE)
}))
selected_sweeps <- c(1L, 2L, 5L, 10L, 20L, 50L, 100L, 200L, 400L)
b2_trace <- do.call(rbind, lapply(endpoints$fit_id, function(fit_id) {
  group <- trace[trace$fit_id == fit_id &
    trace$iteration %in% selected_sweeps, ]
  do.call(rbind, lapply(selected_sweeps, function(sweep) {
    at <- group[group$iteration == sweep & !is.na(group$study) &
                  group$study == 2L, ]
    loading <- at[at$block == "specific_loading", ]
    score <- at[at$block == "specific_score", ]
    function_rows <- at[at$block == "specific_time_function", ]
    stopifnot(nrow(loading) == 1L, nrow(score) == 1L,
              nrow(function_rows) == 2L)
    data.frame(
      fit_id = fit_id, iteration = sweep,
      B2_mean_loading_ppi = loading$ppi_mean,
      B2_loading_posterior_mean_norm = loading$posterior_mean_norm,
      B2_loading_slab_mean_norm = loading$slab_mean_norm,
      B2_score_posterior_mean_norm = score$posterior_mean_norm,
      B2_time_coefficient_norm =
        sqrt(sum(function_rows$posterior_mean_norm^2)),
      stringsAsFactors = FALSE)
  }))
}))
dir.create(result_root)
write.csv(comparison,
          file.path(result_root, "STAGE7C_ORACLE_CONTRAST.csv"),
          row.names = FALSE)
write.csv(pair, file.path(result_root, "ORACLE_PAIR_STABILITY.csv"),
          row.names = FALSE)
write.csv(elbo, file.path(result_root, "ORACLE_ELBO_PATH_SUMMARY.csv"),
          row.names = FALSE)
write.csv(b2_trace, file.path(result_root, "ORACLE_B2_TRACE_SNAPSHOT.csv"),
          row.names = FALSE)
writeLines(c("STAGE7E_DERIVED_COMPLETE", "endpoints=6", "pairs=3",
             "truth_used_for_scientific_diagnosis=TRUE",
             "truth_used_for_fit_selection=FALSE", "new_fit_started=FALSE"),
           file.path(result_root, "DERIVED_COMPLETE.txt"))
cat("STAGE7E_DERIVED_COMPLETE endpoints=6 pairs=3\n")
