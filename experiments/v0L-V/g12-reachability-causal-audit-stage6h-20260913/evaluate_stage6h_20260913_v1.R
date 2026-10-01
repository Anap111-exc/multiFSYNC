#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
experiments_root <- dirname(script_root)
source(file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)
source(file.path(
  script_root, "stage6h_evaluation_helpers_20260913_v1.R"
), local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6H_REACHABILITY_CAUSAL_AUDIT_V1_20260913"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6H_OUTPUT_ROOT")
workers <- as.integer(cli$workers %||% pp_env("STAGE6H_EVAL_WORKERS", "7"))
stage6g_output_root <- pp_env(
  "STAGE6H_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6g_evaluation_root <- pp_env(
  "STAGE6H_STAGE6G_EVALUATION_ROOT",
  file.path(stage6g_output_root, "truth_evaluation_20260912_v1")
)
stage6h_library <- pp_env("STAGE6H_R_LIB")
if (nzchar(stage6h_library)) Sys.setenv(STAGE6A_R_LIB = stage6h_library)
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify Stage 6H output.")
pp_assert(.Platform$OS.type == "unix" && length(workers) == 1L &&
            is.finite(workers) && workers %in% 1:7,
          "Stage 6H evaluation requires Unix and 1--7 workers.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

authorization_path <- file.path(
  script_root, "STAGE6H_EXECUTION_AUTHORIZATION_20260913.txt"
)
authorization <- readLines(authorization_path, warn = FALSE)
pp_assert(all(c(
  "status=AUTHORIZED_POST_UNSEAL_REACHABILITY_CAUSAL_AUDIT",
  "scientific_evaluation_after_endpoint_freeze_authorized=TRUE",
  "truth_used_only_after_all_endpoints_frozen=TRUE",
  "continuation_authorized=FALSE", "cross_dimension_selection_authorized=FALSE"
) %in% authorization), "Stage 6H evaluation authorization is incomplete.")

freeze_marker <- file.path(
  output_root, "STAGE6H_ENDPOINTS_AND_ELBO_SELECTION_FROZEN.txt"
)
fit_complete_marker <- file.path(output_root, "STAGE6H_FITS_COMPLETE.txt")
pp_assert(all(file.exists(c(freeze_marker, fit_complete_marker))),
          "Stage 6H endpoints must be frozen before truth evaluation.")
all_manifest <- pp_read_csv(file.path(output_root,
                                      "ANALYSIS_ENDPOINT_MANIFEST.csv"))
new_manifest <- pp_read_csv(file.path(output_root, "NEW_FIT_MANIFEST.csv"))
terminals <- pp_read_csv(file.path(output_root, "ALL_ENDPOINT_TERMINALS.csv"))
winners <- pp_read_csv(file.path(output_root, "CONDITION_ELBO_WINNERS.csv"))
pp_assert(
  nrow(all_manifest) == 36L && nrow(new_manifest) == 28L &&
    nrow(terminals) == 36L && nrow(winners) == 3L &&
    !anyDuplicated(all_manifest$fit_id) &&
    setequal(all_manifest$fit_id, terminals$fit_id) &&
    all(table(all_manifest$condition_id) == 12L) &&
    !any(terminals$truth_used_for_fit_stopping_or_selection) &&
    !any(terminals$continuation_used),
  "Stage 6H frozen endpoint accounting failed."
)

source_binding <- pp_read_csv(file.path(output_root,
                                        "SOURCE_ARTIFACT_BINDING.csv"))
pp_assert(nrow(source_binding) == 8L &&
            !anyDuplicated(source_binding$artifact_role),
          "Stage 6H source-artifact binding is malformed.")
source_info <- file.info(source_binding$path)
source_matches <- file.exists(source_binding$path) &
  is.finite(source_info$size) &
  as.numeric(source_info$size) == as.numeric(source_binding$bytes) &
  vapply(source_binding$path, pp_sha256, character(1L)) ==
    source_binding$sha256
pp_assert(all(source_matches), paste0(
  "A read-only Stage 6G source artifact changed: ",
  paste(source_binding$artifact_role[!source_matches], collapse = ";")
))

evaluation_root <- file.path(output_root, "truth_evaluation_20260913_v1")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6H truth evaluation already exists; refusing overwrite.")
dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")
pp_write_lines(c(
  "status=EVALUATION_STARTED", paste0("started_utc=", pp_iso_time()),
  "all_endpoints_frozen_before_truth_read=TRUE", "analysis_endpoints=36",
  "new_fits_to_evaluate=28", "stage6g_evaluations_reused_read_only=8",
  paste0("workers=", workers), "cross_dimension_selection=FALSE",
  "continuation=FALSE", "formal_paper_mc=FALSE"
), file.path(evaluation_root, "EVALUATION_STARTED.txt"))

# Truth is read only after endpoint and within-condition ELBO selection freeze.
truth_path <- source_binding$path[
  source_binding$artifact_role == "truth"
]
pp_assert(length(truth_path) == 1L, "Stage 6H truth binding is ambiguous.")
truth <- readRDS(truth_path)
pp_assert(
  identical(truth$bundle_class, "stage6g_sealed_truth_v1") &&
    identical(truth$data_id, "s6g_minimal_L111_M2") &&
    identical(as.integer(truth$data$true_params$L_f), 1L) &&
    identical(as.integer(unlist(truth$data$true_params$L_s_by_study,
                                use.names = FALSE)), c(1L, 1L)) &&
    all(as.integer(truth$data$true_params$M_f) == 2L) &&
    all(as.integer(unlist(truth$data$true_params$M_s,
                          use.names = FALSE)) == 2L),
  "Stage 6H truth identity or dimension check failed."
)
pp_write_csv(data.frame(
  data_id = truth$data_id, truth_sha256 = pp_sha256(truth_path),
  truth_read_after_endpoint_freeze = TRUE, stringsAsFactors = FALSE
), file.path(evaluation_root, "UNSEALED_TRUTH_BINDING.csv"))

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

s6h_new_record <- function(fit_id) {
  path <- file.path(output_root, "fits", fit_id, "terminal_record.rds")
  pp_assert(file.exists(path), paste0("Missing Stage 6H record: ", fit_id))
  record <- readRDS(path)
  pp_assert(
    identical(record$terminal_schema, "G12_STAGE6H_TERMINAL_V1") &&
      identical(record$experiment_id, PP_EXPERIMENT_ID) &&
      identical(record$fit_id, fit_id) &&
      identical(record$terminal_status, "fixed_400_complete") &&
      isTRUE(record$objective_eligible) &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$continuation_used) &&
      identical(as.integer(record$n_cpus), 1L),
    paste0("Ineligible Stage 6H endpoint: ", fit_id)
  )
  record
}

records <- lapply(new_manifest$fit_id, s6h_new_record)
names(records) <- new_manifest$fit_id
expected_tables <- c(
  "feature_component_errors", "factor_kernel_errors",
  "dense_signal_metrics", "loading_recovery", "threshold_sensitivity",
  "unthresholded_assignment", "covariance_operator_error",
  "ppi_diagnostics", "shared_candidate_metrics", "summary",
  "evaluation_contract", "process_contribution_by_role",
  "process_contribution_summary", "loading_scale_by_block",
  "loading_scale_summary", "candidate_activity",
  "candidate_activity_summary", "fpca_rank_activity",
  "fpca_component_activity"
)

evaluate_one <- function(index) {
  row <- new_manifest[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  record <- records[[fit_id]]
  fit_dir <- file.path(output_root, "fits", fit_id)
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    s6h_evaluate_fit(readRDS(file.path(fit_dir, "fit.rds")), truth,
                     evaluator_environment)
  }, warning = function(condition) {
    warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
      class = paste(class(condition), collapse = ";"),
      message = conditionMessage(condition),
      call = paste(deparse(conditionCall(condition)), collapse = " "),
      stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  }), error = function(condition) condition)
  ended <- Sys.time()
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(class = character(), message = character(), call = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  elapsed <- as.numeric(difftime(ended, started, units = "secs"))
  if (inherits(result, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(result))
    ), file.path(result_dir, "EVALUATION_ERROR.txt"))
    return(list(
      status = data.frame(
        fit_id = fit_id, ok = FALSE, elapsed_seconds = elapsed,
        warning_count = nrow(warnings), error = conditionMessage(result),
        stringsAsFactors = FALSE
      ), compact = NULL, metric = NULL
    ))
  }
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(result_dir, paste0(name, ".csv")))
  }
  base <- data.frame(
    fit_id = fit_id, data_id = record$data_id,
    scenario_id = record$scenario_id,
    condition_id = row$condition_id[[1L]], split_id = row$split_id[[1L]],
    input_kind = row$input_kind[[1L]], fit_source = "stage6h_new",
    fit_config_id = record$fit_config_id,
    generating_dimension_configuration = TRUE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, seed_index = record$seed_index,
    fit_seed = record$fit_seed, fit_L_f = record$fit_L_f,
    fit_L_s_1 = record$fit_L_s_1, fit_L_s_2 = record$fit_L_s_2,
    fit_M_f = record$fit_M_f, fit_M_s_1 = record$fit_M_s_1,
    fit_M_s_2 = record$fit_M_s_2, final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    convergence_status = record$convergence_status,
    fit_elapsed_seconds = record$elapsed_seconds,
    fit_warning_count = nrow(record$warnings),
    stringsAsFactors = FALSE
  )
  activity_summary <- result$candidate_activity_summary
  names(activity_summary) <- paste0("activity_", names(activity_summary))
  compact <- cbind(
    base, result$summary, result$process_contribution_summary,
    result$loading_scale_summary, activity_summary
  )
  tables_present <- all(expected_tables %in% names(result)) &&
    all(vapply(result[expected_tables], function(value) {
      is.data.frame(value) && nrow(value) > 0L
    }, logical(1L)))
  fit_object <- readRDS(file.path(fit_dir, "fit.rds"))
  checkpoint_complete <- identical(
    names(fit_object$practical_checkpoints), c("380", "400")
  )
  metric <- data.frame(
    fit_id = fit_id, expected_table_count = length(expected_tables),
    expected_tables_present = tables_present,
    checkpoint_380_400_complete = checkpoint_complete,
    metric_contract_passed = tables_present && checkpoint_complete,
    stringsAsFactors = FALSE
  )
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "all_endpoints_frozen_before_truth_read=TRUE",
    "fit_stopping_and_elbo_selection_used_truth=FALSE",
    "continuation=FALSE", "formal_paper_mc=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(
    status = data.frame(
      fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed,
      warning_count = nrow(warnings), error = "", stringsAsFactors = FALSE
    ), compact = compact, metric = metric
  )
}

results <- parallel::mclapply(
  seq_len(nrow(new_manifest)), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
if (!all(status$ok)) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("failed_fit_ids=", paste(status$fit_id[!status$ok], collapse = ";"))
  ), file.path(evaluation_root, "EVALUATION_ERROR.txt"))
  pp_abort("Stage 6H evaluation errors retained; no selective rerun.")
}
metric_qc <- do.call(rbind, lapply(results, `[[`, "metric"))
pp_assert(nrow(metric_qc) == 28L && all(metric_qc$metric_contract_passed),
          "A Stage 6H evaluation metric contract failed.")
new_science <- pp_rbind_fill(lapply(results, `[[`, "compact"))

stage6g_science_path <- file.path(
  stage6g_evaluation_root, "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv"
)
stage6g_science <- pp_read_csv(stage6g_science_path)
import_manifest <- all_manifest[
  all_manifest$fit_source == "stage6g_read_only", , drop = FALSE
]
imported_science <- stage6g_science[
  stage6g_science$fit_id %in% import_manifest$source_fit_id, , drop = FALSE
]
imported_science <- imported_science[
  match(import_manifest$source_fit_id, imported_science$fit_id), , drop = FALSE
]
pp_assert(nrow(imported_science) == 8L &&
            identical(as.character(imported_science$fit_id),
                      as.character(import_manifest$source_fit_id)),
          "Stage 6G scientific endpoint import failed.")
imported_science$condition_id <- import_manifest$condition_id
imported_science$input_kind <- import_manifest$input_kind
imported_science$fit_source <- "stage6g_read_only"

all_science <- pp_rbind_fill(list(imported_science, new_science))
all_science <- all_science[
  match(all_manifest$fit_id, all_science$fit_id), , drop = FALSE
]
pp_assert(nrow(all_science) == 36L && !anyDuplicated(all_science$fit_id) &&
            identical(as.character(all_science$fit_id),
                      as.character(all_manifest$fit_id)),
          "Stage 6H scientific endpoint merge failed.")
all_science$within_condition_elbo_winner <-
  all_science$fit_id %in% winners$fit_id
all_science$original_stage6g_four_start <-
  all_science$fit_source == "stage6g_read_only"
all_science$stage6h_new_fit <- all_science$fit_source == "stage6h_new"
all_science$ppi_factor_counts_exact <- with(
  all_science,
  selected_shared == truth_L_f &
    selected_specific_study_1 == truth_L_s_1 &
    selected_specific_study_2 == truth_L_s_2
)
all_science$core_truth_recovered <- with(
  all_science,
  global_correct == truth_global_direction_count &
    global_misplaced == 0L & global_missing == 0L &
    study_correct == truth_study_role_count &
    study_misplaced == 0L & study_missing == 0L
)
all_science$exact_thresholded_factor_structure <- with(
  all_science,
  core_truth_recovered & ppi_factor_counts_exact &
    global_extra == 0L & global_duplicate == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L
)
all_science$full_true_feature_component_coverage <- with(
  all_science, feature_components_matched == feature_components_total
)
all_science$core_with_full_feature_coverage <- with(
  all_science, core_truth_recovered & full_true_feature_component_coverage
)

condition_summary <- do.call(rbind, lapply(
  split(all_science, all_science$condition_id), function(rows) {
    winner <- rows[rows$within_condition_elbo_winner, , drop = FALSE]
    first_four <- rows[rows$seed_index <= 4L, , drop = FALSE]
    data.frame(
      condition_id = rows$condition_id[[1L]], input_kind = rows$input_kind[[1L]],
      endpoints = nrow(rows), stage6g_read_only = sum(rows$original_stage6g_four_start),
      stage6h_new = sum(rows$stage6h_new_fit),
      practical_converged = sum(rows$practical_converged),
      ppi_factor_count_exact_hits = sum(rows$ppi_factor_counts_exact),
      core_hits = sum(rows$core_truth_recovered),
      exact_structure_hits = sum(rows$exact_thresholded_factor_structure),
      full_feature_coverage_hits =
        sum(rows$full_true_feature_component_coverage),
      core_full_feature_hits = sum(rows$core_with_full_feature_coverage),
      first_four_core_hits = sum(first_four$core_truth_recovered),
      added_starts_5_12_core_hits =
        sum(rows$core_truth_recovered & rows$seed_index >= 5L),
      winner_fit_id = winner$fit_id[[1L]],
      winner_seed_index = winner$seed_index[[1L]],
      winner_core_hit = winner$core_truth_recovered[[1L]],
      winner_exact_hit = winner$exact_thresholded_factor_structure[[1L]],
      elbo_missed_available_core =
        any(rows$core_truth_recovered) && !winner$core_truth_recovered[[1L]],
      elbo_missed_available_exact =
        any(rows$exact_thresholded_factor_structure) &&
          !winner$exact_thresholded_factor_structure[[1L]],
      observed_signal_nrmse_mean =
        s6h_finite_mean(rows$signal_nrmse_ppi_selected),
      dense_signal_nrmse_mean =
        s6h_finite_mean(rows$dense_signal_nrmse_ppi_selected),
      feature_ise_mean = s6h_finite_mean(rows$feature_total_ise_mean),
      stringsAsFactors = FALSE
    )
  }
))
rownames(condition_summary) <- NULL

holdout <- all_science[all_science$input_kind == "holdout_train", , drop = FALSE]
full <- all_science[all_science$condition_id == "full_data_6_9", , drop = FALSE]
baseline_four <- holdout[holdout$seed_index <= 4L, , drop = FALSE]
baseline_four_any_core <- any(baseline_four$core_truth_recovered)
holdout_twelve_any_core <- any(holdout$core_truth_recovered)
full_twelve_any_core <- any(full$core_truth_recovered)
decision_code <- if (holdout_twelve_any_core && full_twelve_any_core) {
  "BOTH_HOLDOUT_AND_FULL_REACHABLE"
} else if (holdout_twelve_any_core && !full_twelve_any_core) {
  "HOLDOUT_REACHABLE_BUT_FULL_NOT_REACHABLE_UNEXPECTED"
} else if (!holdout_twelve_any_core && full_twelve_any_core) {
  "ONLY_FULL_DATA_REACHABLE_HOLDOUT_DATA_LOSS_SUPPORTED"
} else {
  "NEITHER_REACHABLE_G12_OR_SIGNAL_GEOMETRY_REMAINS_BLOCKER"
}
causal_decision <- data.frame(
  original_holdout_four_start_any_core = baseline_four_any_core,
  holdout_twelve_start_any_core = holdout_twelve_any_core,
  full_data_twelve_start_any_core = full_twelve_any_core,
  start_budget_effect_supported =
    !baseline_four_any_core && holdout_twelve_any_core,
  holdout_data_loss_effect_supported =
    !holdout_twelve_any_core && full_twelve_any_core,
  decision_code = decision_code,
  dimension_criterion_evaluable_next =
    holdout_twelve_any_core && full_twelve_any_core,
  stringsAsFactors = FALSE
)

paired_rows <- lapply(c("holdout_split1_middle", "holdout_split2_upper"),
                      function(condition_id) {
  held <- all_science[all_science$condition_id == condition_id, , drop = FALSE]
  held <- held[order(held$seed_index), , drop = FALSE]
  full_ordered <- full[match(held$seed_index, full$seed_index), , drop = FALSE]
  pp_assert(nrow(held) == 12L && nrow(full_ordered) == 12L &&
              identical(as.integer(held$fit_seed),
                        as.integer(full_ordered$fit_seed)),
            paste0("Stage 6H paired-start join failed: ", condition_id))
  data.frame(
    holdout_condition_id = condition_id,
    seed_index = held$seed_index, fit_seed = held$fit_seed,
    holdout_core = held$core_truth_recovered,
    full_core = full_ordered$core_truth_recovered,
    holdout_exact = held$exact_thresholded_factor_structure,
    full_exact = full_ordered$exact_thresholded_factor_structure,
    holdout_ppi_counts_exact = held$ppi_factor_counts_exact,
    full_ppi_counts_exact = full_ordered$ppi_factor_counts_exact,
    holdout_observed_signal_nrmse = held$signal_nrmse_ppi_selected,
    full_observed_signal_nrmse = full_ordered$signal_nrmse_ppi_selected,
    full_minus_holdout_observed_nrmse =
      full_ordered$signal_nrmse_ppi_selected - held$signal_nrmse_ppi_selected,
    holdout_dense_signal_nrmse = held$dense_signal_nrmse_ppi_selected,
    full_dense_signal_nrmse = full_ordered$dense_signal_nrmse_ppi_selected,
    full_minus_holdout_dense_nrmse =
      full_ordered$dense_signal_nrmse_ppi_selected -
        held$dense_signal_nrmse_ppi_selected,
    stringsAsFactors = FALSE
  )
})
paired <- do.call(rbind, paired_rows)
rownames(paired) <- NULL

runtime <- data.frame(
  new_fits = nrow(new_science),
  single_core_hours = sum(new_science$fit_elapsed_seconds) / 3600,
  median_fit_minutes = stats::median(new_science$fit_elapsed_seconds) / 60,
  max_fit_minutes = max(new_science$fit_elapsed_seconds) / 60,
  peak_memory_max_bytes = max(terminals$peak_memory_bytes[terminals$fit_source ==
                                                           "stage6h_new"]),
  evaluation_elapsed_seconds = sum(status$elapsed_seconds),
  evaluation_warning_count = sum(status$warning_count),
  stringsAsFactors = FALSE
)

pp_write_csv(metric_qc,
             file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
pp_write_csv(new_science,
             file.path(evaluation_root, "NEW_STAGE6H_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(imported_science,
             file.path(evaluation_root, "IMPORTED_STAGE6G_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(all_science,
             file.path(evaluation_root, "ALL_STAGE6H_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(condition_summary,
             file.path(evaluation_root, "CONDITION_REACHABILITY_SUMMARY.csv"))
pp_write_csv(causal_decision,
             file.path(evaluation_root, "CAUSAL_DECISION.csv"))
pp_write_csv(paired,
             file.path(evaluation_root, "PAIRED_START_CONTRASTS.csv"))
pp_write_csv(runtime,
             file.path(evaluation_root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "analysis_endpoints=36", "new_endpoints_evaluated=28",
  "stage6g_endpoint_evaluations_reused=8", "evaluation_errors=0",
  paste0("core_hits=", sum(all_science$core_truth_recovered)),
  paste0("exact_structure_hits=",
         sum(all_science$exact_thresholded_factor_structure)),
  paste0("decision_code=", decision_code),
  "all_endpoints_frozen_before_truth_read=TRUE",
  "fit_stopping_and_elbo_selection_used_truth=FALSE",
  "cross_dimension_selection=FALSE", "continuation=FALSE",
  "formal_paper_mc=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE6H_EVALUATION_COMPLETE endpoints=36 new=28 decision=",
    decision_code, "\n", sep = "")
