#!/usr/bin/env Rscript

# Post-selection scientific evaluation for the Route-B K=8/9 independent confirmation.
# This script reads sealed truth only after the 48 endpoints and six within-
# data-by-K winners have been verified and their hashes recorded.  It starts
# no fit and no continuation, and reuses the frozen loading-first evaluator.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
support_repo <- Sys.getenv("ROUTE_B_EVALUATION_SUPPORT_REPO", unset = "")
support_v0lv <- if (nzchar(support_repo)) {
  file.path(support_repo, "experiments", "v0L-V")
} else file.path(script_root, "..")
stage6a_support_path <- file.path(
  support_v0lv, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
)
stage6h_support_path <- file.path(
  support_v0lv, "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R"
)
source(stage6a_support_path, local = FALSE)
source(stage6h_support_path, local = FALSE)
if (nzchar(support_repo)) {
  PP_SNAPSHOT <- file.path(support_repo, "experiments", "v0L-V",
                           "rc1_snapshot")
}

PP_EXPERIMENT_ID <- "ROUTE_B_K8_K9_INDEPENDENT_CONFIRMATION_V1_20260925"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% "1")
pp_assert(nzchar(output_root) && dir.exists(output_root) &&
            workers %in% 1:8,
          "A prepared --output-root and 1--8 workers are required.")
pp_assert(pp_bool(cli$`truth-unseal-authorized`, FALSE),
          "Scientific evaluation requires explicit truth-unseal authorization.")
pp_assert(
  identical(tolower(pp_sha256(stage6a_support_path)),
            "bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6") &&
    identical(tolower(pp_sha256(stage6h_support_path)),
              "5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221") &&
    dir.exists(PP_SNAPSHOT),
  "The bound Stage 6A/6H evaluation support or frozen snapshot changed."
)

evaluation_root <- file.path(
  output_root, "scientific_evaluation_20260926_v1")
pp_assert(!dir.exists(evaluation_root),
          "Evaluation output already exists; refusing overwrite.")

data_manifest_path <- file.path(output_root, "DATA_MANIFEST.csv")
fit_manifest_path <- file.path(output_root, "FIT_MANIFEST.csv")
selection_root <- file.path(output_root, "truth_free_selection")
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "TRUTH_FREE_WINNERS.csv")
selection_status_path <- file.path(selection_root, "SELECTION_STATUS.csv")
frozen_path <- file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt")

for (path in c(
  data_manifest_path, fit_manifest_path, endpoints_path, winners_path,
  selection_status_path, frozen_path,
  file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"),
  file.path(output_root, "SUPERVISOR_EXIT_STATUS.txt")
)) pp_assert(file.exists(path), paste0("Missing required artifact: ", path))

data_manifest <- pp_read_csv(data_manifest_path)
manifest <- pp_read_csv(fit_manifest_path)
endpoints <- pp_read_csv(endpoints_path)
winners <- pp_read_csv(winners_path)
selection_status <- pp_read_csv(selection_status_path)
frozen <- readLines(frozen_path, warn = FALSE)

pp_assert(
  nrow(data_manifest) == 3L && nrow(manifest) == 48L &&
    nrow(endpoints) == 48L && nrow(winners) == 6L &&
    nrow(selection_status) == 6L &&
    all(data_manifest$experiment_id == PP_EXPERIMENT_ID) &&
    all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
    setequal(data_manifest$data_id, c("routebc_01", "routebc_02", "routebc_03")) &&
    setequal(manifest$fit_K, c(8L, 9L)) &&
    all(table(manifest$data_id, manifest$fit_K) == 8L) &&
    !anyDuplicated(manifest$fit_id) &&
    setequal(manifest$fit_id, endpoints$fit_id) &&
    all(endpoints$status == "COMPLETE") &&
    all(endpoints$objective_eligible) &&
    all(is.finite(endpoints$final_ordinary_elbo)) &&
    all(endpoints$annealing_sweeps == 99L) &&
    all(endpoints$ordinary_t1_sweeps == 400L) &&
    !any(endpoints$truth_used_for_fit_or_selection) &&
    !any(manifest$truth_available_to_fit) &&
    !any(manifest$truth_available_to_stopping) &&
    !any(manifest$truth_available_to_selection) &&
    !any(manifest$continuation) && !any(manifest$actual_racing) &&
    !any(manifest$automatic_800) &&
    all(manifest$n_cpus == 1L) &&
    all(selection_status$n_registered == 8L) &&
    all(selection_status$n_eligible == 8L) &&
    all(selection_status$winner_status == "SELECTED"),
  "Manifest, endpoint, budget, or truth-isolation contract failed."
)

pp_assert(all(c(
  "status=TRUTH_FREE_SELECTION_FROZEN",
  "rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_K",
  "cross_K_ELBO_comparison=FALSE",
  "truth_used_for_selection=FALSE",
  "truth_unseal_authorized=FALSE",
  "registered_fits=48", "selected_winners=6", "smoke=FALSE"
) %in% frozen), "Frozen selection marker is incomplete or changed.")

for (stratum in unique(manifest$selection_stratum_id)) {
  group <- endpoints[
    endpoints$selection_stratum_id == stratum &
      endpoints$objective_eligible &
      is.finite(endpoints$final_ordinary_elbo), , drop = FALSE]
  chosen <- winners[winners$selection_stratum_id == stratum, , drop = FALSE]
  pp_assert(nrow(group) == 8L && nrow(chosen) == 1L &&
              identical(chosen$fit_id[[1L]],
                        group$fit_id[[which.max(group$final_ordinary_elbo)]]),
            paste0("Winner is not maximum eligible ELBO in ", stratum, "."))
}

fit_complete <- file.path(output_root, "fits", manifest$fit_id,
                          "FIT_COMPLETE.txt")
fit_error <- file.path(output_root, "fits", manifest$fit_id,
                       "FIT_ERROR.txt")
pp_assert(all(file.exists(fit_complete)) && !any(file.exists(fit_error)),
          "Fit completion accounting failed.")

# Verify the three source files bound before generation/fitting.  The new
# evaluator itself is recorded in the authorization file below.
provenance <- pp_read_csv(file.path(output_root, "SOURCE_PROVENANCE.csv"))
pp_assert(nrow(provenance) == 3L &&
            all(c("relative_path", "sha256") %in% names(provenance)) &&
            !anyDuplicated(provenance$relative_path),
          "SOURCE_PROVENANCE.csv is malformed.")
source_matches <- vapply(seq_len(nrow(provenance)), function(index) {
  path <- file.path(PP_REPO_ROOT, provenance$relative_path[[index]])
  file.exists(path) && identical(
    tolower(pp_sha256(path)), tolower(provenance$sha256[[index]]))
}, logical(1L))
pp_assert(all(source_matches), paste0(
  "Bound source mismatch: ",
  paste(provenance$relative_path[!source_matches], collapse = ";")))

dir.create(evaluation_root, recursive = TRUE, mode = "0700")
pp_write_lines(c(
  "status=POST_SELECTION_TRUTH_UNSEAL_AUTHORIZED_BY_USER",
  "authorization_date=2026-09-26",
  paste0("frozen_endpoints_sha256=", pp_sha256(endpoints_path)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  paste0("evaluator_sha256=", pp_sha256(script_file)),
  paste0("stage6a_support_sha256=", pp_sha256(stage6a_support_path)),
  paste0("stage6h_support_sha256=", pp_sha256(stage6h_support_path)),
  "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_K",
  "cross_K_ELBO_comparison=FALSE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE"
), file.path(evaluation_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))

# This is the first sealed-truth read.  All endpoint, winner, provenance, and
# authorization checks above have already completed.
data_ids <- unique(manifest$data_id)
truths <- setNames(lapply(data_ids, function(data_id) {
  truth <- readRDS(file.path(
    output_root, "data", data_id, "sealed_truth", "truth_bundle.rds"))
  pp_assert(identical(truth$bundle_class, "route_b_sealed_truth") &&
              identical(truth$experiment_id, PP_EXPERIMENT_ID) &&
              identical(truth$data_id, data_id) &&
              isFALSE(truth$truth_unseal_authorized) &&
              identical(truth$data$true_params$truth_basis_method,
                        "iid_common_splinet_stiefel") &&
              identical(as.integer(truth$data$true_params$truth_basis_dim), 8L),
            paste0("Sealed truth identity mismatch: ", data_id))
  truth
}), data_ids)

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

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
    result$evaluation_contract$route_b_fit_K <- row$fit_K[[1L]]
    result$evaluation_contract$route_b_truth_basis_method <-
      "iid_common_splinet_stiefel"
    result$evaluation_contract$oracle_L_M <- TRUE
    result
  }, phase = "route_b_k8_k9_confirmation_post_selection_scientific_evaluation")
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
    pp_write_csv(result[[name]],
                 file.path(result_dir, paste0(name, ".csv")))
  }
  terminal <- endpoints[endpoints$fit_id == fit_id, , drop = FALSE]
  activity <- result$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  base <- data.frame(
    fit_id = fit_id, data_id = row$data_id[[1L]],
    data_index = row$data_index[[1L]], fit_K = row$fit_K[[1L]],
    fit_config_id = row$fit_config_id[[1L]],
    selection_stratum_id = row$selection_stratum_id[[1L]],
    paired_start_id = row$paired_start_id[[1L]],
    start_index = row$start_index[[1L]], fit_seed = row$fit_seed[[1L]],
    final_ordinary_elbo = terminal$final_ordinary_elbo[[1L]],
    objective_eligible = terminal$objective_eligible[[1L]],
    strict_practical_converged =
      terminal$strict_practical_converged[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    truth_free_winner = fit_id %in% winners$fit_id,
    truth_used_for_fit_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE)
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
pp_write_csv(science, file.path(
  evaluation_root, "ALL_ROUTE_B_SCIENTIFIC_RESULTS.csv"))
winner_science <- science[science$truth_free_winner, , drop = FALSE]
pp_write_csv(winner_science, file.path(
  evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))

summary <- do.call(rbind, lapply(split(science, science$fit_K), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(
    fit_K = x$fit_K[[1L]], starts = nrow(x), datasets = length(unique(x$data_id)),
    core_hits = sum(x$core_truth_recovered),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    full_feature_coverage_hits = sum(x$full_true_feature_component_coverage),
    winner_core_hits = sum(winner$core_truth_recovered),
    winner_exact_structure_hits =
      sum(winner$exact_thresholded_factor_structure),
    winner_full_feature_coverage_hits =
      sum(winner$full_true_feature_component_coverage),
    mean_fit_elapsed_seconds = mean(x$fit_elapsed_seconds),
    broad_generalization_claim = FALSE, stringsAsFactors = FALSE)
}))
pp_write_csv(summary, file.path(evaluation_root, "K_ARM_SUMMARY.csv"))
pp_write_csv(science, file.path(evaluation_root,
                                "PAIRED_ENDPOINT_SCIENCE_LONG.csv"))

pp_write_lines(c(
  "status=ROUTE_B_K8_K9_CONFIRMATION_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("fits_evaluated=", nrow(science)),
  paste0("winners_evaluated=", nrow(winner_science)),
  "truth_read_after_all_within_data_by_K_selection_freezes=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "cross_K_ELBO_comparison=FALSE",
  "new_fit_or_continuation_started=FALSE",
  "broad_generalization_claim=FALSE",
  paste0("completed_utc=", pp_iso_time())
), file.path(evaluation_root, "ROUTE_B_CONFIRMATION_EVALUATION_COMPLETE.txt"))

cat("ROUTE_B_CONFIRMATION_EVALUATION_COMPLETE fits=", nrow(science),
    " winners=", nrow(winner_science), "\n", sep = "")


