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
source(file.path(experiments_root,
                 "g12-dimension-overspec-stage6a-20260831",
                 "stage6a_common_20260831_v1.R"), local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6G_DIMENSION_CRITERION_FALSIFICATION_V1_20260911"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6G_OUTPUT_ROOT")
workers <- as.integer(cli$workers %||% pp_env("STAGE6G_EVAL_WORKERS", "7"))
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "--output-root must identify the completed Stage 6G output.")
pp_assert(.Platform$OS.type == "unix" && length(workers) == 1L &&
            is.finite(workers) && workers %in% 1:7,
          "Stage 6G evaluation requires Unix and 1--7 workers.")
stage6g_library <- pp_env("STAGE6G_R_LIB")
if (nzchar(stage6g_library)) Sys.setenv(STAGE6A_R_LIB = stage6g_library)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

selection_root <- file.path(output_root, "truth_free_selection")
selection_path <- file.path(
  selection_root, "FINAL_TRUTH_FREE_DIMENSION_SELECTIONS.csv"
)
winner_path <- file.path(selection_root, "CONFIG_SPLIT_WINNERS.csv")
artifact_hash_path <- file.path(
  selection_root, "TRUTH_FREE_ARTIFACT_HASHES.csv"
)
selection_marker <- file.path(selection_root,
                              "DIMENSION_SELECTION_FROZEN.txt")
start_marker <- file.path(selection_root, "START_SELECTION_FROZEN.txt")
experiment_complete <- file.path(output_root,
                                 "STAGE6G_TRUTH_FREE_COMPLETE.txt")
authorization_path <- file.path(output_root,
                                "TRUTH_UNSEAL_AUTHORIZATION.txt")
evaluation_root <- file.path(output_root,
                             "truth_evaluation_20260912_v1")

pp_assert(all(file.exists(c(
  selection_path, winner_path, artifact_hash_path, selection_marker,
  start_marker, experiment_complete
))), "The frozen truth-free Stage 6G evidence is incomplete.")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6G evaluation already exists; refusing overwrite.")

artifact_manifest <- pp_read_csv(artifact_hash_path)
pp_assert(all(c("relative_to_output_root", "bytes", "sha256") %in%
                names(artifact_manifest)) && nrow(artifact_manifest) > 0L &&
            !anyDuplicated(artifact_manifest$relative_to_output_root),
          "Stage 6G truth-free artifact hash manifest is malformed.")
artifact_paths <- file.path(output_root,
                            artifact_manifest$relative_to_output_root)
artifact_info <- file.info(artifact_paths)
artifact_matches <- file.exists(artifact_paths) &
  is.finite(artifact_info$size) &
  as.numeric(artifact_info$size) == as.numeric(artifact_manifest$bytes) &
  vapply(artifact_paths, pp_sha256, character(1L)) ==
    artifact_manifest$sha256
pp_assert(all(artifact_matches), paste0(
  "Truth-free artifact changed after freezing: ",
  paste(artifact_manifest$relative_to_output_root[!artifact_matches],
        collapse = ";")
))

selection_hash <- pp_sha256(selection_path)
winner_hash <- pp_sha256(winner_path)
artifact_manifest_hash <- pp_sha256(artifact_hash_path)
pp_assert(file.exists(authorization_path),
          "Explicit Stage 6G truth-unseal authorization is missing.")
authorization <- readLines(authorization_path, warn = FALSE)
required_authorization <- c(
  "status=TRUTH_UNSEAL_AND_SCIENTIFIC_EVALUATION_AUTHORIZED",
  paste0("final_dimension_selection_sha256=", selection_hash),
  paste0("configuration_split_winner_sha256=", winner_hash),
  paste0("truth_free_artifact_hash_manifest_sha256=",
         artifact_manifest_hash),
  "truth_unseal_authorized=TRUE",
  "scientific_evaluation_authorized=TRUE",
  "new_fits_authorized_after_unseal=FALSE",
  "continuation_authorized=FALSE",
  "selection_recomputation_after_unseal=FALSE",
  "authorization_scope=stage6g_existing_120_endpoints_only"
)
pp_assert(all(required_authorization %in% authorization),
          "Stage 6G authorization is incomplete or hash-mismatched.")

manifest <- pp_read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
pp_assert(nrow(manifest) == 120L && !anyDuplicated(manifest$fit_id) &&
            length(unique(manifest$selection_stratum_id)) == 30L &&
            all(table(manifest$selection_stratum_id) == 4L),
          "Stage 6G runtime fit manifest is malformed.")
fit_root <- file.path(output_root, "fits")
fit_dirs <- file.path(fit_root, manifest$fit_id)
pp_assert(all(dir.exists(fit_dirs)) &&
            all(file.exists(file.path(fit_dirs, "FIT_COMPLETE.txt"))) &&
            all(file.exists(file.path(fit_dirs, "fit.rds"))) &&
            all(file.exists(file.path(fit_dirs, "terminal_record.rds"))),
          "Stage 6G evaluation requires all 120 complete fit artifacts.")
pp_assert(length(list.files(fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE,
                            full.names = TRUE)) == 0L,
          "A Stage 6G fit error is present.")

records <- lapply(seq_len(nrow(manifest)), function(index) {
  record <- readRDS(file.path(fit_dirs[[index]], "terminal_record.rds"))
  pp_assert(
    identical(record$terminal_schema, "G12_STAGE6G_TERMINAL_V1") &&
      identical(record$experiment_id, PP_EXPERIMENT_ID) &&
      identical(record$fit_id, manifest$fit_id[[index]]) &&
      identical(record$terminal_status, "fixed_400_complete") &&
      isTRUE(record$objective_eligible) &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$truth_used_for_holdout_scoring) &&
      !isTRUE(record$continuation_used) &&
      identical(as.integer(record$n_cpus), 1L),
    paste0("Ineligible Stage 6G endpoint: ", manifest$fit_id[[index]])
  )
  record
})
names(records) <- manifest$fit_id

winners <- pp_read_csv(winner_path)
selections <- pp_read_csv(selection_path)
pp_assert(nrow(winners) == 30L && !anyDuplicated(winners$fit_id) &&
            all(table(winners$data_id, winners$fit_config_id) == 2L) &&
            nrow(selections) == 3L && !anyDuplicated(selections$data_id) &&
            !any(selections$selection_used_truth),
          "Malformed frozen Stage 6G selection tables.")
winner_ids <- as.character(winners$fit_id)

map_selected <- merge(
  selections[c("data_id", "selected_config_id")],
  winners[c("data_id", "split_id", "fit_config_id", "fit_id")],
  by.x = c("data_id", "selected_config_id"),
  by.y = c("data_id", "fit_config_id"), sort = FALSE
)
map_minimum <- merge(
  selections[c("data_id", "minimum_loss_config_id")],
  winners[c("data_id", "split_id", "fit_config_id", "fit_id")],
  by.x = c("data_id", "minimum_loss_config_id"),
  by.y = c("data_id", "fit_config_id"), sort = FALSE
)
pp_assert(nrow(map_selected) == 6L && nrow(map_minimum) == 6L,
          "Stage 6G dimension selections do not map to split winners.")
one_se_ids <- as.character(map_selected$fit_id)
minimum_ids <- as.character(map_minimum$fit_id)

winner_scores <- lapply(winner_ids, function(fit_id) {
  path <- file.path(output_root, "truth_free_scores", fit_id,
                    "HOLDOUT_SUMMARY.csv")
  value <- pp_read_csv(path)
  pp_assert(nrow(value) == 1L && value$factor_mode[[1L]] == "all" &&
              !isTRUE(value$truth_used[[1L]]) &&
              identical(as.character(value$fit_id[[1L]]), fit_id),
            paste0("Malformed holdout score for ", fit_id))
  value
})
winner_scores <- pp_rbind_fill(winner_scores)
rownames(winner_scores) <- winner_scores$fit_id

# Truth-bearing registration is deliberately read only after the authorization
# and all truth-free artifacts have been verified.
scenario_truth <- pp_read_csv(file.path(
  output_root, "sealed_registration", "SCENARIOS_WITH_GENERATING_TRUTH.csv"
))
pp_assert(nrow(scenario_truth) == 3L &&
            !anyDuplicated(scenario_truth$scenario_id) &&
            setequal(scenario_truth$scenario_id, selections$data_id),
          "Stage 6G sealed scenario registration is malformed.")

pp_load_package()
pp_source_runtime(TRUE)
truth_paths <- file.path(output_root, "data", selections$data_id,
                         "sealed_truth", "truth_bundle.rds")
pp_assert(all(file.exists(truth_paths)), "A Stage 6G truth bundle is missing.")
truths <- lapply(truth_paths, readRDS)
names(truths) <- selections$data_id
for (data_id in names(truths)) {
  truth <- truths[[data_id]]
  row <- scenario_truth[scenario_truth$scenario_id == data_id,
                        , drop = FALSE]
  true_l_s <- as.integer(unlist(truth$data$true_params$L_s_by_study,
                                use.names = FALSE))
  pp_assert(
    identical(truth$bundle_class, "stage6g_sealed_truth_v1") &&
      identical(truth$experiment_id, PP_EXPERIMENT_ID) &&
      identical(truth$data_id, data_id) && nrow(row) == 1L &&
      identical(as.integer(truth$data$true_params$L_f),
                as.integer(row$truth_L_f[[1L]])) &&
      identical(true_l_s,
                as.integer(c(row$truth_L_s_1[[1L]],
                             row$truth_L_s_2[[1L]]))) &&
      all(as.integer(truth$data$true_params$M_f) ==
            as.integer(row$truth_M_common[[1L]])) &&
      all(as.integer(unlist(truth$data$true_params$M_s,
                            use.names = FALSE)) ==
            as.integer(row$truth_M_common[[1L]])),
    paste0("Stage 6G truth identity/dimension check failed: ", data_id)
  )
}
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

s6g_finite_mean <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}

s6g_multi_additional_metrics <- function(
    fit, truth, evaluation, evaluator_environment) {
  data <- truth$data
  reported_paths <- evaluator_environment$v0g_reported_paths_by_study(
    fit, data
  )
  true_paths <- evaluator_environment$v0h_truth_paths_by_study(data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected",
                            , drop = FALSE]
  expected_roles <- data$true_params$S * data$true_params$L_f +
    sum(as.integer(unlist(data$true_params$L_s_by_study,
                          use.names = FALSE)))
  pp_assert(nrow(functional) == expected_roles,
            "Unexpected Stage 6G study-role truth-row count.")
  process_rows <- lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_role <- as.character(item$true_role[[1L]])
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_role == "shared") {
      data$true_params$a_true[, true_factor]
    } else {
      data$true_params$b_true[[study]][, true_factor]
    }
    true_path <- true_paths[[true_role]][[study]][, true_factor]
    true_contribution <- as.numeric(outer(true_path, true_loading))
    matched <- isTRUE(item$matched_by_loading[[1L]])
    if (matched) {
      estimated_role <- as.character(item$estimated_block[[1L]])
      estimated_factor <- as.integer(item$estimated_factor[[1L]])
      estimated_loading <- if (estimated_role == "shared") {
        fit$mu_q_a_hat[, estimated_factor]
      } else {
        fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      }
      estimated_path <- reported_paths[[estimated_role]][[study]][,
                                                                  estimated_factor]
      aligned_path <- as.numeric(item$loading_sign[[1L]]) * estimated_path
      process_error <- pp_nrmse(aligned_path, true_path)
      contribution_error <- pp_nrmse(
        as.numeric(outer(estimated_path, estimated_loading)),
        true_contribution
      )
      loading_norm_ratio <- pp_vector_l2(estimated_loading) /
        pp_vector_l2(true_loading)
      loading_cosine <- pp_abs_cosine(estimated_loading, true_loading)
      estimated_path_rms <- pp_vector_rms(estimated_path)
    } else {
      estimated_role <- NA_character_
      estimated_factor <- NA_integer_
      process_error <- contribution_error <- loading_norm_ratio <-
        loading_cosine <- estimated_path_rms <- NA_real_
    }
    data.frame(
      study = study, true_role = true_role, true_factor = true_factor,
      true_id = item$true_id[[1L]], estimated_block = estimated_role,
      estimated_factor = estimated_factor,
      loading_structure_status = item$loading_structure_status[[1L]],
      matched_by_loading = matched,
      loading_norm_ratio = loading_norm_ratio,
      loading_abs_cosine = loading_cosine,
      true_path_rms = pp_vector_rms(true_path),
      estimated_path_rms = estimated_path_rms,
      factor_process_nrmse_matched = process_error,
      factor_process_nrmse_missing_as_zero = if (matched) process_error else 1,
      complete_contribution_nrmse_matched = contribution_error,
      complete_contribution_nrmse_missing_as_zero =
        if (matched) contribution_error else 1,
      stringsAsFactors = FALSE
    )
  })
  process <- do.call(rbind, process_rows)
  process_summary <- data.frame(
    study_role_matched = sum(process$matched_by_loading),
    study_role_truth_total = nrow(process),
    factor_process_nrmse_mean_matched =
      s6g_finite_mean(process$factor_process_nrmse_matched),
    factor_process_nrmse_mean_missing_as_zero =
      mean(process$factor_process_nrmse_missing_as_zero),
    complete_contribution_nrmse_mean_matched =
      s6g_finite_mean(process$complete_contribution_nrmse_matched),
    complete_contribution_nrmse_mean_missing_as_zero =
      mean(process$complete_contribution_nrmse_missing_as_zero),
    stringsAsFactors = FALSE
  )

  loading_scale_rows <- lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_role <- as.character(item$true_role[[1L]])
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_role == "shared") {
      data$true_params$a_true[, true_factor]
    } else {
      data$true_params$b_true[[study]][, true_factor]
    }
    matched <- isTRUE(item$matched_by_loading[[1L]])
    estimated_block <- if (matched) {
      as.character(item$estimated_block[[1L]])
    } else NA_character_
    estimated_factor <- if (matched) {
      as.integer(item$estimated_factor[[1L]])
    } else NA_integer_
    if (!matched) {
      raw <- canonical <- rep(NA_real_, length(true_loading))
      factor_scale <- factor_ppi <- NA_real_
    } else if (estimated_block == "shared") {
      original_factor <- as.integer(
        fit$factor_order_shared[[estimated_factor]]
      )
      raw <- fit$mu_q_a_original[, original_factor]
      canonical <- fit$mu_q_a_hat[, estimated_factor]
      factor_scale <- fit$factor_scale_shared[[estimated_factor]]
      factor_ppi <- fit$factor_ppi_shared[[estimated_factor]]
    } else {
      original_factor <- as.integer(
        fit$factor_order_specific[[study]][[estimated_factor]]
      )
      raw <- fit$mu_q_b_specific_original[[study]][, original_factor]
      canonical <- fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      factor_scale <- fit$factor_scale_specific[[study]][[estimated_factor]]
      factor_ppi <- fit$factor_ppi_specific[[study]][[estimated_factor]]
    }
    data.frame(
      study = study, true_role = true_role, true_factor = true_factor,
      estimated_block = estimated_block,
      estimated_factor = estimated_factor,
      matched_by_loading = matched, factor_ppi = factor_ppi,
      factor_scale = factor_scale,
      raw_loading_norm_ratio = pp_vector_l2(raw) /
        pp_vector_l2(true_loading),
      canonical_loading_norm_ratio = pp_vector_l2(canonical) /
        pp_vector_l2(true_loading),
      raw_loading_abs_cosine = pp_abs_cosine(raw, true_loading),
      canonical_loading_abs_cosine = pp_abs_cosine(canonical, true_loading),
      stringsAsFactors = FALSE
    )
  })
  loading_scale <- do.call(rbind, loading_scale_rows)
  loading_scale_summary <- data.frame(
    raw_loading_norm_ratio_mean =
      s6g_finite_mean(loading_scale$raw_loading_norm_ratio),
    canonical_loading_norm_ratio_mean =
      s6g_finite_mean(loading_scale$canonical_loading_norm_ratio),
    raw_loading_abs_cosine_mean =
      s6g_finite_mean(loading_scale$raw_loading_abs_cosine),
    canonical_loading_abs_cosine_mean =
      s6g_finite_mean(loading_scale$canonical_loading_abs_cosine),
    factor_scale_mean = s6g_finite_mean(loading_scale$factor_scale),
    stringsAsFactors = FALSE
  )
  list(
    process_contribution_by_role = process,
    process_contribution_summary = process_summary,
    loading_scale_by_block = loading_scale,
    loading_scale_summary = loading_scale_summary
  )
}

s6g_evaluate_fit <- function(fit, record, truth) {
  # candidate2's frozen evaluator checks its historical bundle-class name.
  # This in-memory adapter changes no truth payload or file bytes.
  adapter <- truth
  adapter$bundle_class <- "v0lv_sealed_truth_candidate2"
  result <- candidate2_evaluate_multi(
    fit, adapter, PP_SNAPSHOT,
    evaluator_environment = evaluator_environment, mode = "smoke"
  )
  result$evaluation_contract$stage6g_original_truth_bundle_class <-
    truth$bundle_class
  result$evaluation_contract$truth_dimensions_used_for_evaluation_only <- TRUE
  result$evaluation_contract$truth_free_selection_frozen_before_unseal <- TRUE
  result$evaluation_contract$cross_dimension_elbo_selection <- FALSE
  result$evaluation_contract$new_fit_after_unseal <- FALSE
  result$evaluation_contract$continuation_after_unseal <- FALSE
  additional <- s6g_multi_additional_metrics(
    fit, truth, result, evaluator_environment
  )
  dimension_activity <- pp_dimension_activity_metrics(
    fit, adapter, result, evaluator_environment
  )
  c(result, additional, dimension_activity)
}

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

dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")
truth_hashes <- vapply(truth_paths, pp_sha256, character(1L))
pp_write_csv(data.frame(
  data_id = names(truths), sealed_truth_sha256 = truth_hashes,
  stringsAsFactors = FALSE
), file.path(evaluation_root, "UNSEALED_TRUTH_BINDING.csv"))
pp_write_lines(c(
  "status=EVALUATION_STARTED", paste0("started_utc=", pp_iso_time()),
  paste0("workers=", workers), "registered_endpoints=120",
  paste0("final_dimension_selection_sha256=", selection_hash),
  paste0("configuration_split_winner_sha256=", winner_hash),
  paste0("truth_free_artifact_hash_manifest_sha256=",
         artifact_manifest_hash),
  "truth_free_artifacts_verified_before_truth_read=TRUE",
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_STARTED.txt"))

evaluate_one <- function(index) {
  record <- records[[index]]
  fit_id <- record$fit_id
  fit_dir <- file.path(fit_root, fit_id)
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  started <- Sys.time()
  warning_rows <- list()
  result <- tryCatch(withCallingHandlers({
    s6g_evaluate_fit(readRDS(file.path(fit_dir, "fit.rds")), record,
                     truths[[record$data_id]])
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
    pp_write_csv(result[[name]],
                 file.path(result_dir, paste0(name, ".csv")))
  }
  scenario <- scenario_truth[
    scenario_truth$scenario_id == record$data_id, , drop = FALSE
  ]
  score <- if (fit_id %in% winner_ids) {
    winner_scores[fit_id, , drop = FALSE]
  } else NULL
  base <- data.frame(
    fit_id = fit_id, data_id = record$data_id,
    scenario_id = record$scenario_id, split_id = record$split_id,
    fit_config_id = record$fit_config_id,
    expected_truth_fit_config_id = scenario$expected_truth_fit_config_id,
    generating_dimension_configuration =
      record$fit_config_id == scenario$expected_truth_fit_config_id,
    truth_L_f = as.integer(scenario$truth_L_f),
    truth_L_s_1 = as.integer(scenario$truth_L_s_1),
    truth_L_s_2 = as.integer(scenario$truth_L_s_2),
    truth_M_common = as.integer(scenario$truth_M_common),
    seed_index = record$seed_index, fit_seed = record$fit_seed,
    fit_L_f = record$fit_L_f, fit_L_s_1 = record$fit_L_s_1,
    fit_L_s_2 = record$fit_L_s_2, fit_M_f = record$fit_M_f,
    fit_M_s_1 = record$fit_M_s_1, fit_M_s_2 = record$fit_M_s_2,
    final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    convergence_status = record$convergence_status,
    fit_elapsed_seconds = record$elapsed_seconds,
    fit_warning_count = nrow(record$warnings),
    heldout_rmse = if (is.null(score)) NA_real_ else score$heldout_rmse,
    heldout_nrmse = if (is.null(score)) NA_real_ else score$heldout_nrmse,
    heldout_mae = if (is.null(score)) NA_real_ else score$heldout_mae,
    within_configuration_split_winner = fit_id %in% winner_ids,
    selected_by_one_se = fit_id %in% one_se_ids,
    selected_by_minimum_loss = fit_id %in% minimum_ids,
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
    paste0("final_dimension_selection_sha256=", selection_hash),
    "selection_frozen_before_truth_read=TRUE",
    "fit_stopping_and_selection_used_truth=FALSE",
    "new_fits_after_unseal=0", "continuation_started=FALSE",
    "formal_paper_mc_result=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(
    status = data.frame(
      fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed,
      warning_count = nrow(warnings), error = "", stringsAsFactors = FALSE
    ), compact = compact, metric = metric
  )
}

results <- parallel::mclapply(
  seq_along(records), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
if (!all(status$ok)) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
    paste0("failed_fit_ids=", paste(status$fit_id[!status$ok],
                                    collapse = ";"))
  ), file.path(evaluation_root, "EVALUATION_ERROR.txt"))
  pp_abort("Stage 6G evaluation errors retained; no selective rerun.")
}

metric_qc <- do.call(rbind, lapply(results, `[[`, "metric"))
pp_assert(nrow(metric_qc) == 120L && all(metric_qc$metric_contract_passed),
          "One or more Stage 6G evaluation metric contracts failed.")
all_science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
config_split_winners <- all_science[
  all_science$within_configuration_split_winner, , drop = FALSE
]
one_se_science <- all_science[all_science$selected_by_one_se,
                              , drop = FALSE]
minimum_science <- all_science[all_science$selected_by_minimum_loss,
                               , drop = FALSE]
pp_assert(nrow(all_science) == 120L &&
            nrow(config_split_winners) == 30L &&
            nrow(one_se_science) == 6L && nrow(minimum_science) == 6L,
          "Stage 6G evaluated selection mapping failed.")

pp_write_csv(metric_qc,
             file.path(evaluation_root, "METRIC_COMPLETENESS_QC.csv"))
pp_write_csv(all_science,
             file.path(evaluation_root,
                       "ALL_STAGE6G_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(config_split_winners,
             file.path(evaluation_root,
                       "CONFIG_SPLIT_WINNER_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(one_se_science,
             file.path(evaluation_root,
                       "ONE_SE_SELECTED_ENDPOINT_SCIENTIFIC_RESULTS.csv"))
pp_write_csv(minimum_science,
             file.path(evaluation_root,
                       "MINIMUM_LOSS_ENDPOINT_SCIENTIFIC_RESULTS.csv"))
pp_write_lines(c(
  "status=EVALUATION_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "fits_evaluated=120", "evaluation_errors=0",
  "configuration_split_winners_evaluated=30",
  "one_se_selected_endpoints_evaluated=6",
  "minimum_loss_endpoints_evaluated=6",
  paste0("final_dimension_selection_sha256=", selection_hash),
  paste0("configuration_split_winner_sha256=", winner_hash),
  paste0("truth_free_artifact_hash_manifest_sha256=",
         artifact_manifest_hash),
  "truth_free_artifacts_verified_before_truth_read=TRUE",
  "selection_frozen_before_truth_read=TRUE", "new_fits_after_unseal=0",
  "continuation_started=FALSE", "formal_paper_mc_result=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE6G_EVALUATION_COMPLETE fits=120 errors=0 one_se_endpoints=6\n")
