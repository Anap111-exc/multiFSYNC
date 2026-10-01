#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
source(file.path(
  script_root, "..", "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)
source(file.path(
  script_root, "..", "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R"
), local = FALSE)

PP_EXPERIMENT_ID <- "G12_LOADING_DGP_ROBUSTNESS_STAGE6L_V1_20260915"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
smoke <- pp_bool(cli$smoke, FALSE)
workers <- as.integer(cli$workers %||% "1")
pp_assert(nzchar(output_root) && dir.exists(output_root) &&
            workers %in% 1:8,
          "A prepared --output-root and 1--8 workers are required.")
pp_assert(pp_bool(cli$`truth-unseal-authorized`, FALSE),
          "Scientific evaluation requires explicit truth-unseal authorization.")
pp_verify_source_binding()

selection_root <- file.path(output_root, "truth_free_selection")
frozen_path <- file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt")
pp_assert(file.exists(frozen_path), "Truth-free selection has not been frozen.")
frozen <- readLines(frozen_path, warn = FALSE)
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "G12_TRUTH_FREE_WINNERS.csv")
pp_assert(all(c(
  "status=TRUTH_FREE_SELECTION_FROZEN", "selection_used_truth=FALSE",
  "cross_arm_elbo_comparison=FALSE",
  paste0("endpoints_sha256=", pp_sha256(endpoints_path)),
  paste0("winners_sha256=", pp_sha256(winners_path))
) %in% frozen), "The frozen selection content or hashes changed.")

manifest <- pp_read_csv(file.path(script_root, "FIT_MANIFEST.csv"))
if (smoke) manifest <- manifest[manifest$seed_index == 1L, , drop = FALSE]
terminals <- pp_read_csv(endpoints_path)
winners <- pp_read_csv(winners_path)
pp_assert(nrow(manifest) == if (smoke) 2L else 8L,
          "The registered evaluation fit count changed.")
pp_assert(nrow(terminals) == nrow(manifest) && nrow(winners) == 2L &&
            setequal(manifest$fit_id, terminals$fit_id) &&
            all(terminals$objective_eligible) &&
            !any(terminals$truth_used_for_fit_stopping_or_selection) &&
            !any(terminals$continuation_used) &&
            all(table(winners$loading_arm) == 1L),
          "Frozen endpoint accounting or truth isolation failed.")
if (!smoke) pp_assert(all(terminals$actual_T1_sweeps == 400L),
                     "A full endpoint did not complete 400 ordinary sweeps.")

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
evaluation_root <- file.path(output_root, "scientific_evaluation_20260915_v1")
pp_assert(!dir.exists(evaluation_root), "Evaluation output already exists.")
dir.create(evaluation_root, recursive = TRUE, mode = "0700")

# The first truth read is below, after both selection tables have been verified.
data_ids <- unique(manifest$data_id)
truths <- setNames(lapply(data_ids, function(data_id) {
  path <- file.path(output_root, "data", data_id, "sealed_truth",
                    "truth_bundle.rds")
  truth <- readRDS(path)
  pp_assert(identical(truth$experiment_id, PP_EXPERIMENT_ID) &&
              identical(truth$data_id, data_id), "Sealed truth identity failed.")
  truth
}), data_ids)

add_recovery_flags <- function(value) {
  value$ppi_factor_counts_exact <- with(value,
    selected_shared == truth_L_f &
      selected_specific_study_1 == truth_L_s_1 &
      selected_specific_study_2 == truth_L_s_2)
  value$core_truth_recovered <- with(value,
    global_correct == truth_global_direction_count &
      global_misplaced == 0L & global_missing == 0L &
      study_correct == truth_study_role_count &
      study_misplaced == 0L & study_missing == 0L)
  value$exact_thresholded_factor_structure <- with(value,
    core_truth_recovered & ppi_factor_counts_exact &
      global_extra == 0L & global_duplicate == 0L &
      study_candidate_extra == 0L & study_candidate_duplicate == 0L)
  value$full_true_feature_component_coverage <- with(value,
    feature_components_matched == feature_components_total)
  value
}

evaluate_one <- function(index) {
  row <- manifest[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, mode = "0700")
  captured <- c2_capture_conditions({
    fit <- readRDS(file.path(output_root, "fits", fit_id, "fit.rds"))
    result <- s6h_evaluate_fit(
      fit, truths[[row$data_id[[1L]]]], evaluator_environment)
    result$evaluation_contract$oracle_L_M <- TRUE
    result$evaluation_contract$stage6l_loading_arm <- row$loading_arm[[1L]]
    result
  }, phase = "stage6l_post_freeze_evaluation")
  pp_write_csv(captured$warnings,
               file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  status <- data.frame(
    fit_id = fit_id, ok = is.null(captured$error),
    elapsed_seconds = captured$elapsed_seconds,
    error = captured$error$message %||% "", stringsAsFactors = FALSE)
  if (!status$ok) {
    pp_write_lines(c("status=ERROR", paste0("message=", status$error)),
                   file.path(result_dir, "EVALUATION_ERROR.txt"))
    return(list(status = status, compact = NULL))
  }
  result <- captured$value
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(result_dir, paste0(name, ".csv")))
  }
  terminal <- terminals[terminals$fit_id == fit_id, , drop = FALSE]
  base <- data.frame(
    fit_id = fit_id, data_id = row$data_id[[1L]],
    loading_arm = row$loading_arm[[1L]],
    paired_start_id = row$paired_start_id[[1L]],
    seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
    final_elbo = terminal$final_elbo[[1L]],
    practical_converged = terminal$practical_converged[[1L]],
    convergence_status = terminal$convergence_status[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    fit_warning_count = terminal$warning_count[[1L]],
    truth_free_winner = fit_id %in% winners$fit_id,
    truth_used_for_fit_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE)
  activity <- result$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  compact <- cbind(base, result$summary,
                   result$process_contribution_summary,
                   result$loading_scale_summary, activity)
  pp_write_lines(c("status=COMPLETE", paste0("fit_id=", fit_id),
                   "truth_used_for_fit_stopping_or_selection=FALSE",
                   "new_fit_or_continuation_started=FALSE"),
                 file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(status = status, compact = compact)
}

results <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(manifest)), evaluate_one,
                    mc.cores = workers, mc.preschedule = FALSE,
                    mc.set.seed = FALSE)
} else lapply(seq_len(nrow(manifest)), evaluate_one)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok), "Evaluation errors retained; no selective rerun.")
science <- add_recovery_flags(pp_rbind_fill(lapply(results, `[[`, "compact")))
science <- science[match(manifest$fit_id, science$fit_id), , drop = FALSE]
pp_write_csv(science, file.path(evaluation_root,
                                "ALL_STAGE6L_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(science[science$truth_free_winner, , drop = FALSE],
               file.path(evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))

pair_columns <- intersect(c(
  "paired_start_id", "fit_id", "core_truth_recovered",
  "exact_thresholded_factor_structure", "full_true_feature_component_coverage",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "covariance_operator_relative_error_mean", "loading_relative_l2_error_mean",
  "correctly_placed_score_abs_cor_mean",
  "factor_process_nrmse_mean_missing_as_zero",
  "complete_contribution_nrmse_mean_missing_as_zero", "fit_elapsed_seconds"
), names(science))
paired <- merge(
  science[science$loading_arm == "controlled", pair_columns, drop = FALSE],
  science[science$loading_arm == "jaoua_beta10", pair_columns, drop = FALSE],
  by = "paired_start_id", suffixes = c("_controlled", "_beta"), sort = TRUE)
pp_write_csv(paired, file.path(evaluation_root,
                               "PAIRED_CONTROLLED_VS_BETA_SCIENCE.csv"))
summary <- do.call(rbind, lapply(split(science, science$loading_arm), function(x) {
  data.frame(loading_arm = x$loading_arm[[1L]], starts = nrow(x),
             core_hits = sum(x$core_truth_recovered),
             exact_structure_hits = sum(x$exact_thresholded_factor_structure),
             full_feature_coverage_hits =
               sum(x$full_true_feature_component_coverage),
             winner_fit_id = x$fit_id[x$truth_free_winner],
             winner_core = x$core_truth_recovered[x$truth_free_winner],
             paired_logical_datasets = 1L, broad_generalization_claim = FALSE,
             stringsAsFactors = FALSE)
}))
pp_write_csv(summary, file.path(evaluation_root, "LOADING_ARM_SUMMARY.csv"))
pp_write_lines(c("status=STAGE6L_SCIENTIFIC_EVALUATION_COMPLETE",
                 paste0("fits_evaluated=", nrow(science)),
                 "truth_read_after_both_within_arm_selection_freezes=TRUE",
                 "truth_used_for_fit_stopping_or_selection=FALSE",
                 "cross_arm_elbo_comparison=FALSE",
                 "new_fit_or_continuation_started=FALSE",
                 "broad_generalization_claim=FALSE",
                 paste0("smoke=", smoke),
                 paste0("completed_utc=", pp_iso_time())),
               file.path(evaluation_root, "STAGE6L_EVALUATION_COMPLETE.txt"))
cat("STAGE6L_EVALUATION_COMPLETE fits=", nrow(science), "\n", sep = "")
