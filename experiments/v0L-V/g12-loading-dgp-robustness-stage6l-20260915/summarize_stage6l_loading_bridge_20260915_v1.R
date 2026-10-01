#!/usr/bin/env Rscript
# Post-unseal descriptive summarizer. It reads existing outputs only; it has
# no generator, fitting, continuation or winner-selection call.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L, dir.exists(arguments[[1L]]))
output_root <- arguments[[1L]]
evaluation_root <- file.path(output_root, "scientific_evaluation_20260915_v1")
stopifnot(file.exists(file.path(evaluation_root,
                               "STAGE6L_EVALUATION_COMPLETE.txt")))
summary_root <- file.path(evaluation_root, "descriptive_summary_20260915_v1")
if (file.exists(summary_root)) stop("Refusing to overwrite descriptive summary.")
dir.create(summary_root, recursive = TRUE)
read_csv <- function(path) utils::read.csv(path, check.names = FALSE,
                                          stringsAsFactors = FALSE)
write_csv <- function(value, name) utils::write.csv(
  value, file.path(summary_root, name), row.names = FALSE, na = "")
science <- read_csv(file.path(evaluation_root,
                              "ALL_STAGE6L_SCIENTIFIC_RESULTS.csv"))
terminals <- read_csv(file.path(output_root, "truth_free_selection",
                                "ALL_TRUTH_FREE_ENDPOINTS.csv"))
stopifnot(nrow(science) == 8L, !anyDuplicated(science$fit_id),
          setequal(science$fit_id, terminals$fit_id),
          !any(science$truth_used_for_fit_stopping_or_selection))

# These context formulas are the same as the registered generation-only audit.
rms <- function(value) sqrt(mean(as.numeric(value)^2))
flatten <- function(value) unlist(value, recursive = TRUE, use.names = FALSE)
jaccard <- function(left, right) sum(left & right) / sum(left | right)
cosine <- function(left, right) abs(sum(left * right)) /
  sqrt(sum(left^2) * sum(right^2))
blocks <- geometry <- datasets <- list()
for (data_id in unique(science$data_id)) {
  truth <- readRDS(file.path(output_root, "data", data_id, "sealed_truth",
                            "truth_bundle.rds"))
  data <- truth$data
  tp <- data$true_params
  arm <- truth$loading_arm
  values <- c(list(shared_1 = as.numeric(tp$a_true[, 1L])),
              list(specific_s1_1 = as.numeric(tp$b_true[[1L]][, 1L])),
              list(specific_s2_1 = as.numeric(tp$b_true[[2L]][, 1L])))
  omega <- c(tp$generator_omega_a_true[[1L]],
             tp$generator_omega_b_true[[1L]][[1L]],
             tp$generator_omega_b_true[[2L]][[1L]])
  blocks[[arm]] <- do.call(rbind, lapply(seq_along(values), function(index) {
    value <- values[[index]]
    active <- sum(value != 0)
    data.frame(data_id = data_id, loading_arm = arm,
               block_id = names(values)[[index]],
               sampled_generator_omega = omega[[index]],
               active_count = active, active_fraction = active / tp$p,
               loading_l2 = sqrt(sum(value^2)),
               loading_l2_per_sqrt_active = sqrt(sum(value^2) / active))
  }))
  local_geometry <- do.call(rbind, lapply(1:2, function(study) {
    left <- values[[1L]]; right <- values[[study + 1L]]
    eigenvalues <- eigen(crossprod(cbind(left, right)), symmetric = TRUE,
                         only.values = TRUE)$values
    data.frame(data_id = data_id, loading_arm = arm, study = study,
               shared_specific_support_jaccard = jaccard(left != 0, right != 0),
               shared_specific_abs_cosine = cosine(left, right),
               loading_gram_min_eigenvalue = min(eigenvalues),
               loading_gram_condition = max(eigenvalues) / min(eigenvalues))
  }))
  geometry[[arm]] <- local_geometry
  signal <- flatten(tp$signal_true_values)
  mean_part <- flatten(tp$mu_true_values)
  beta_part <- flatten(tp$beta_true_values)
  noise <- flatten(tp$noise_true_values)
  factor_part <- signal - mean_part - beta_part
  datasets[[arm]] <- data.frame(
    data_id = data_id, loading_arm = arm, observation_scalar_count = length(signal),
    signal_rms = rms(signal), factor_signal_rms = rms(factor_part),
    mean_rms = rms(mean_part), noise_rms = rms(noise),
    factor_to_noise_rms_ratio = rms(factor_part) / rms(noise),
    shared_specific_jaccard_mean =
      mean(local_geometry$shared_specific_support_jaccard),
    shared_specific_abs_cosine_mean =
      mean(local_geometry$shared_specific_abs_cosine),
    loading_gram_min_eigenvalue =
      min(local_geometry$loading_gram_min_eigenvalue),
    loading_gram_condition_max = max(local_geometry$loading_gram_condition),
    response_identity_max_abs_error = max(abs(flatten(data$Y) - signal - noise)))
}
write_csv(do.call(rbind, blocks), "LOADING_BLOCK_CONTEXT.csv")
write_csv(do.call(rbind, geometry), "STUDY_LOADING_GEOMETRY.csv")
write_csv(do.call(rbind, datasets), "BRIDGE_DATASET_CONTEXT.csv")

columns <- intersect(c(
  "loading_arm", "seed_index", "fit_seed", "truth_free_winner", "final_elbo",
  "selected_counts", "core_truth_recovered", "exact_thresholded_factor_structure",
  "global_correct", "global_misplaced", "global_missing", "global_extra",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "feature_components_matched", "feature_components_total", "shared_effective_M",
  "specific_effective_M_study_1", "specific_effective_M_study_2",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_total_ise_mean", "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "loading_support_ppi_auc_mean", "canonical_loading_abs_cosine_mean",
  "correctly_placed_function_recall_mean", "correctly_placed_score_abs_cor_mean",
  "correctly_placed_trajectory_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero", "R_A", "R_B1", "R_B2",
  "P_A", "P_B1", "P_B2", "max_cross_block_L", "practical_converged",
  "fit_warning_count", "fit_elapsed_seconds"
), names(science))
write_csv(science[, columns, drop = FALSE], "ALL_ENDPOINT_PRESENTATION.csv")
write_csv(science[science$truth_free_winner, columns, drop = FALSE],
          "WINNER_PRESENTATION.csv")
runtime <- do.call(rbind, lapply(split(terminals, terminals$loading_arm), function(x) {
  data.frame(loading_arm = x$loading_arm[[1L]], fits = nrow(x),
             elapsed_minutes_mean = mean(x$elapsed_seconds) / 60,
             elapsed_minutes_min = min(x$elapsed_seconds) / 60,
             elapsed_minutes_max = max(x$elapsed_seconds) / 60,
             objective_eligible = sum(x$objective_eligible),
             practical_converged = sum(x$practical_converged),
             fit_warnings = sum(x$warning_count))
}))
write_csv(runtime, "FIT_RUNTIME_SUMMARY.csv")
invisible(lapply(c("multiFSYNC", "Splinets", "digest"), function(package) {
  loadNamespace(package)
}))
writeLines(capture.output(sessionInfo()),
           file.path(summary_root, "SESSION_INFO_SUMMARY_RUNTIME.txt"))
writeLines(c("status=DESCRIPTIVE_SUMMARY_COMPLETE", "existing_fits_readonly=TRUE",
             "winner_selection_changed=FALSE", "new_fit_or_continuation=FALSE",
             "data_replicates=one_paired_logical_dataset",
             "fit_seeds_not_treated_as_independent_dataset_replicates=TRUE"),
           file.path(summary_root, "DESCRIPTIVE_SUMMARY_COMPLETE.txt"))
cat("DESCRIPTIVE_SUMMARY_COMPLETE\n")
