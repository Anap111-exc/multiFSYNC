#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
))
output <- file.path(root, "evaluation_20260828_v1")
stage4e_root <- Sys.getenv(
  "STAGE4G_STAGE4E_ROOT",
  unset = "/root/v0lv-g12-focused-reachability-stage4e-20260827-v1"
)
stage4c_root <- Sys.getenv(
  "STAGE4G_STAGE4C_ROOT",
  unset = "/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1"
)
snapshot_root <- Sys.getenv(
  "STAGE4G_SNAPSHOT_ROOT",
  unset = "/root/multiFSYNC-git-rc1/experiments/v0L-V/rc1_snapshot"
)
development_library <- Sys.getenv(
  "STAGE4G_DEVELOPMENT_LIBRARY",
  unset = "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5-lib"
)

abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(value, message) if (!isTRUE(value)) abort(message)
read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing input: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(
    value, file.path(output, name), row.names = FALSE,
    quote = TRUE, na = ""
  )
}
find_one <- function(path, pattern) {
  result <- list.files(
    path, pattern = pattern, recursive = TRUE, full.names = TRUE
  )
  assert(length(result) == 1L,
         paste0("Expected one ", pattern, "; found ", length(result), "."))
  normalizePath(result[[1L]], winslash = "/", mustWork = TRUE)
}
safe_mean <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
vector_rms <- function(value) {
  value <- as.numeric(value)
  if (length(value) && all(is.finite(value))) sqrt(mean(value^2)) else NA_real_
}
vector_l2 <- function(value) {
  value <- as.numeric(value)
  if (length(value) && all(is.finite(value))) sqrt(sum(value^2)) else NA_real_
}
abs_cosine <- function(left, right) {
  left <- as.numeric(left); right <- as.numeric(right)
  denominator <- vector_l2(left) * vector_l2(right)
  if (length(left) != length(right) || !length(left) ||
      any(!is.finite(c(left, right)))) return(NA_real_)
  if (!is.finite(denominator) || denominator <= 0) return(0)
  abs(sum(left * right) / denominator)
}
nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate); truth <- as.numeric(truth)
  assert(length(estimate) == length(truth) && length(truth) > 0L,
         "NRMSE vectors are not aligned.")
  denominator <- sqrt(mean(truth^2))
  assert(is.finite(denominator) && denominator > 0,
         "Truth NRMSE denominator is invalid.")
  sqrt(mean((estimate - truth)^2)) / denominator
}
iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

assert(!dir.exists(output),
       "Refusing to overwrite evaluation_20260828_v1.")
assert(file.exists(file.path(root, "CONTINUATIONS_COMPLETE.txt")),
       "Stage-4G continuations are incomplete.")
manifest <- read_csv(file.path(root, "CONTINUATION_MANIFEST.csv"))
terminals <- read_csv(file.path(root, "ALL_2_TERMINALS.csv"))
assert(nrow(manifest) == 2L && nrow(terminals) == 2L &&
         identical(manifest$continuation_id, terminals$continuation_id) &&
         all(terminals$terminal_status == "fixed_800_complete") &&
         all(terminals$objective_eligible) &&
         !any(terminals$truth_used) &&
         !any(terminals$winner_selection_reopened) &&
         !any(terminals$automatic_further_extension),
       "Stage-4G terminal set is invalid.")

.libPaths(unique(c(development_library, .libPaths())))
wrapper <- find_one(snapshot_root, "^v0lv_candidate2_evaluation[.]R$")
sys.source(wrapper, envir = globalenv(), keep.source = TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(snapshot_root)
dir.create(output, recursive = TRUE, mode = "0700")

structure_exact <- function(summary) {
  with(
    summary,
    selected_counts == "1,1,1" &
      global_correct == truth_global_direction_count &
      global_misplaced == 0L & global_missing == 0L &
      global_extra == 0L & global_duplicate == 0L &
      study_correct == truth_study_role_count &
      study_misplaced == 0L & study_missing == 0L &
      study_candidate_extra == 0L & study_candidate_duplicate == 0L
  )
}
decompose_gate <- function(diagnostics, control) {
  result <- diagnostics
  result$objective_rate_pass <-
    result$elbo_abs_rate <= control$elbo_abs_rate |
    result$elbo_per_response_rate <= control$elbo_per_response_rate
  result$objective_monotone_pass <- result$monotone & result$long_monotone
  result$short_fitted_pass <-
    result$fitted_nrmse <= control$fitted_nrmse
  result$short_rss_pass <- result$rss_rel <= control$rss_rel
  result$short_ppi_quantile_pass <-
    result$ppi_quantile_abs <= control$ppi_quantile_max_abs
  result$short_factor_ppi_pass <-
    result$factor_ppi_max_abs <= control$factor_ppi_max_abs
  result$long_fitted_pass <-
    result$long_fitted_nrmse <= control$long_fitted_nrmse
  result$long_rss_pass <- result$long_rss_rel <= control$long_rss_rel
  result$long_ppi_quantile_pass <-
    result$long_ppi_quantile_abs <= control$long_ppi_quantile_max_abs
  result$long_factor_ppi_pass <-
    result$long_factor_ppi_max_abs <= control$long_factor_ppi_max_abs
  result
}
direction_payload <- function(snapshot, role) {
  if (role == "A") {
    loading <- snapshot$shared_loading[, 1L]
    feature <- snapshot$shared_feature[[1L]]
    score <- lapply(snapshot$shared_score_mean, `[[`, 1L)
  } else {
    study <- if (role == "B1") 1L else 2L
    loading <- snapshot$specific_loading[[study]][, 1L]
    feature <- snapshot$specific_feature[[study]][[1L]]
    score <- list(snapshot$specific_score_mean[[study]][[1L]])
  }
  process <- unlist(lapply(score, function(value) {
    as.numeric(as.matrix(value) %*% t(as.matrix(feature)))
  }), use.names = FALSE)
  list(
    loading = as.numeric(loading), feature = as.matrix(feature),
    score = score, process = process
  )
}
process_recovery <- function(fit, evaluation, truth, metadata) {
  data <- truth$data
  reported_paths <-
    evaluator_environment$v0g_reported_paths_by_study(fit, data)
  truth_paths <- evaluator_environment$v0h_truth_paths_by_study(data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected", , drop = FALSE]
  assert(nrow(functional) == 4L, "Expected four truth study-role rows.")
  do.call(rbind, lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_block <- item$true_role[[1L]]
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_block == "shared") {
      data$true_params$a_true[, true_factor]
    } else {
      data$true_params$b_true[[study]][, true_factor]
    }
    true_path <- truth_paths[[true_block]][[study]][, true_factor]
    matched <- isTRUE(item$matched_by_loading[[1L]])
    if (matched) {
      estimated_block <- item$estimated_block[[1L]]
      estimated_factor <- as.integer(item$estimated_factor[[1L]])
      estimated_loading <- if (estimated_block == "shared") {
        fit$mu_q_a_hat[, estimated_factor]
      } else {
        fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      }
      estimated_path <-
        reported_paths[[estimated_block]][[study]][, estimated_factor]
      process_error <- nrmse(
        item$loading_sign[[1L]] * estimated_path, true_path
      )
      contribution_error <- nrmse(
        as.numeric(outer(estimated_path, estimated_loading)),
        as.numeric(outer(true_path, true_loading))
      )
      loading_norm_ratio <-
        vector_l2(estimated_loading) / vector_l2(true_loading)
      loading_direction <- abs_cosine(estimated_loading, true_loading)
    } else {
      process_error <- contribution_error <- loading_norm_ratio <-
        loading_direction <- NA_real_
    }
    data.frame(
      metadata, study = study, true_id = item$true_id[[1L]],
      true_role = true_block,
      estimated_block = item$estimated_block[[1L]],
      loading_structure_status = item$loading_structure_status[[1L]],
      matched_by_loading = matched,
      loading_norm_ratio = loading_norm_ratio,
      loading_abs_cosine = loading_direction,
      factor_process_nrmse_matched = process_error,
      factor_process_nrmse_missing_as_zero =
        if (matched) process_error else 1,
      complete_contribution_nrmse_matched = contribution_error,
      complete_contribution_nrmse_missing_as_zero =
        if (matched) contribution_error else 1,
      stringsAsFactors = FALSE
    )
  }))
}

endpoint_rows <- list()
change_rows <- list()
gate_rows <- list()
timeline_rows <- list()
process_rows <- list()
evaluation_action_rows <- list()
roles <- c("A", "B1", "B2")
metric_names <- c(
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "dense_mean_only_signal_nrmse", "loading_relative_l2_error_mean",
  "feature_components_matched", "feature_components_total",
  "matched_function_recall_mean", "matched_trajectory_abs_cor_mean",
  "matched_score_abs_cor_mean", "feature_total_ise_mean",
  "feature_projection_lower_bound_ise_mean",
  "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean", "P_A", "P_B1", "P_B2"
)

for (index in seq_len(nrow(manifest))) {
  row <- manifest[index, , drop = FALSE]
  continuation_id <- row$continuation_id[[1L]]
  source_fit_path <- file.path(
    stage4e_root, row$source_fit_relative_path[[1L]]
  )
  continuation_dir <- file.path(root, "continuations", continuation_id)
  continuation_fit_path <- file.path(continuation_dir, "fit.rds")
  source_evaluation_path <- file.path(
    stage4e_root, "truth_evaluation_20260827_v1", "endpoints",
    row$source_fit_id[[1L]], "evaluation.rds"
  )
  truth_path <- file.path(
    stage4c_root, "data", row$data_id[[1L]],
    "sealed_truth", "truth_bundle.rds"
  )
  source_direction_path <- file.path(
    dirname(source_fit_path), "DIRECTION_SNAPSHOTS.rds"
  )
  continuation_direction_path <- file.path(
    continuation_dir, "DIRECTION_SNAPSHOTS.rds"
  )
  paths <- c(
    source_fit_path, continuation_fit_path, source_evaluation_path,
    truth_path, source_direction_path, continuation_direction_path
  )
  assert(all(file.exists(paths)),
         paste0("A Stage-4G evaluation input is missing: ", continuation_id))
  source_fit <- readRDS(source_fit_path)
  continuation_fit <- readRDS(continuation_fit_path)
  source_evaluation <- readRDS(source_evaluation_path)
  truth <- readRDS(truth_path)
  source_directions <- readRDS(source_direction_path)
  continuation_directions <- readRDS(continuation_direction_path)
  assert(identical(truth$data_id, row$data_id[[1L]]) &&
           identical(names(continuation_directions),
                     as.character(c(400L, seq.int(420L, 800L, 20L)))),
         "Truth or continuation checkpoint identity mismatch.")

  started <- Sys.time()
  warnings <- character()
  continuation_evaluation <- withCallingHandlers(
    candidate2_evaluate_multi(
      continuation_fit, truth, snapshot_root,
      evaluator_environment = evaluator_environment,
      mode = "smoke", binding_dir = NULL
    ), warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  ended <- Sys.time()
  endpoint_dir <- file.path(output, "endpoints", continuation_id)
  dir.create(endpoint_dir, recursive = TRUE, mode = "0700")
  saveRDS(
    continuation_evaluation,
    file.path(endpoint_dir, "evaluation.rds"), version = 3
  )
  for (name in names(continuation_evaluation)) {
    if (is.data.frame(continuation_evaluation[[name]])) {
      utils::write.csv(
        continuation_evaluation[[name]],
        file.path(endpoint_dir, paste0(name, ".csv")),
        row.names = FALSE, quote = TRUE, na = ""
      )
    }
  }
  evaluation_action_rows[[index]] <- data.frame(
    continuation_id = continuation_id,
    source_fit_id = row$source_fit_id[[1L]],
    data_id = row$data_id[[1L]], status = "COMPLETE",
    warning_count = length(warnings),
    warnings = paste(unique(warnings), collapse = " | "),
    started_utc = iso_time(started), ended_utc = iso_time(ended),
    elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
    truth_used_for_continuation = FALSE,
    truth_used_for_evaluation = TRUE,
    winner_selection_reopened = FALSE,
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )

  source_summary <- source_evaluation$summary
  continuation_summary <- continuation_evaluation$summary
  assert(nrow(source_summary) == 1L && nrow(continuation_summary) == 1L,
         "An evaluator summary is not one row.")
  source_exact <- structure_exact(source_summary)
  continuation_exact <- structure_exact(continuation_summary)
  available_metrics <- intersect(
    metric_names, intersect(names(source_summary), names(continuation_summary))
  )
  metadata <- data.frame(
    continuation_id = continuation_id,
    source_fit_id = row$source_fit_id[[1L]],
    data_id = row$data_id[[1L]],
    scenario_id = row$scenario_id[[1L]],
    selection_role = row$selection_role[[1L]],
    stringsAsFactors = FALSE
  )
  endpoint_rows[[length(endpoint_rows) + 1L]] <- rbind(
    data.frame(
      metadata, cumulative_t1_sweeps = 400L,
      selected_counts = source_summary$selected_counts,
      structure_exact = source_exact,
      final_elbo = as.numeric(row$source_final_elbo[[1L]]),
      source_summary[, available_metrics, drop = FALSE],
      stringsAsFactors = FALSE
    ),
    data.frame(
      metadata, cumulative_t1_sweeps = 800L,
      selected_counts = continuation_summary$selected_counts,
      structure_exact = continuation_exact,
      final_elbo = as.numeric(tail(continuation_fit$ELBO, 1L)),
      continuation_summary[, available_metrics, drop = FALSE],
      stringsAsFactors = FALSE
    )
  )
  numeric_metrics <- available_metrics[
    vapply(source_summary[, available_metrics, drop = FALSE],
           is.numeric, logical(1L))
  ]
  for (metric in c("final_elbo", numeric_metrics)) {
    source_value <- if (metric == "final_elbo") {
      as.numeric(row$source_final_elbo[[1L]])
    } else as.numeric(source_summary[[metric]][[1L]])
    continuation_value <- if (metric == "final_elbo") {
      as.numeric(tail(continuation_fit$ELBO, 1L))
    } else as.numeric(continuation_summary[[metric]][[1L]])
    change_rows[[length(change_rows) + 1L]] <- data.frame(
      metadata, metric = metric, source_400 = source_value,
      continuation_800 = continuation_value,
      delta_800_minus_400 = continuation_value - source_value,
      relative_change = (continuation_value - source_value) /
        (abs(source_value) + .Machine$double.eps),
      stringsAsFactors = FALSE
    )
  }
  source_factor_scale <- c(
    A = source_fit$factor_scale_shared[[1L]],
    B1 = source_fit$factor_scale_specific[[1L]][[1L]],
    B2 = source_fit$factor_scale_specific[[2L]][[1L]]
  )
  continuation_factor_scale <- c(
    A = continuation_fit$factor_scale_shared[[1L]],
    B1 = continuation_fit$factor_scale_specific[[1L]][[1L]],
    B2 = continuation_fit$factor_scale_specific[[2L]][[1L]]
  )
  for (role in names(source_factor_scale)) {
    source_value <- as.numeric(source_factor_scale[[role]])
    continuation_value <- as.numeric(continuation_factor_scale[[role]])
    change_rows[[length(change_rows) + 1L]] <- data.frame(
      metadata, metric = paste0("canonical_factor_scale_", role),
      source_400 = source_value,
      continuation_800 = continuation_value,
      delta_800_minus_400 = continuation_value - source_value,
      relative_change = (continuation_value - source_value) /
        (abs(source_value) + .Machine$double.eps),
      stringsAsFactors = FALSE
    )
  }

  source_process <- process_recovery(
    source_fit, source_evaluation, truth,
    data.frame(metadata, cumulative_t1_sweeps = 400L)
  )
  continuation_process <- process_recovery(
    continuation_fit, continuation_evaluation, truth,
    data.frame(metadata, cumulative_t1_sweeps = 800L)
  )
  process_rows[[index]] <- rbind(source_process, continuation_process)

  source_gate <- decompose_gate(
    source_fit$practical_diagnostics, source_fit$practical_control
  )
  continuation_gate <- decompose_gate(
    continuation_fit$practical_diagnostics,
    continuation_fit$practical_control
  )
  for (gate_item in list(
    list(stage = 400L, value = source_gate,
         actual_converged = source_fit$practical_converged),
    list(stage = 800L, value = continuation_gate,
         actual_converged = FALSE)
  )) {
    value <- gate_item$value
    final <- tail(value, 1L)
    last_five <- tail(value, 5L)
    gate_rows[[length(gate_rows) + 1L]] <- data.frame(
      metadata, cumulative_t1_sweeps = gate_item$stage,
      actual_run_practical_converged = gate_item$actual_converged,
      final_gate_pass = final$pass[[1L]],
      last_five_all_pass = all(last_five$pass),
      counterfactual_public_g12_converged = all(last_five$pass),
      final_objective_pass = final$objective_pass[[1L]],
      final_output_pass = final$output_pass[[1L]],
      fitted_nrmse = final$fitted_nrmse[[1L]],
      rss_rel = final$rss_rel[[1L]],
      ppi_quantile_abs = final$ppi_quantile_abs[[1L]],
      factor_ppi_max_abs = final$factor_ppi_max_abs[[1L]],
      long_fitted_nrmse = final$long_fitted_nrmse[[1L]],
      long_rss_rel = final$long_rss_rel[[1L]],
      long_ppi_quantile_abs = final$long_ppi_quantile_abs[[1L]],
      long_factor_ppi_max_abs = final$long_factor_ppi_max_abs[[1L]],
      stringsAsFactors = FALSE
    )
  }

  truth_data <- truth$data
  source_ppi <- source_fit$practical_checkpoints[["400"]]$snapshot$factor_ppi
  continuation_ppi <- lapply(
    continuation_fit$practical_checkpoints,
    function(value) value$snapshot$factor_ppi
  )
  all_ppi <- c(
    list(`400` = as.numeric(source_ppi)),
    stats::setNames(
      lapply(continuation_ppi, as.numeric),
      as.character(400L + as.integer(names(continuation_ppi)))
    )
  )
  for (checkpoint in names(continuation_directions)) {
    snapshot <- continuation_directions[[checkpoint]]
    sweep <- as.integer(checkpoint)
    factor_ppi <- all_ppi[[checkpoint]]
    assert(length(factor_ppi) == 3L,
           "Continuation checkpoint factor-PPI shape changed.")
    for (role_index in seq_along(roles)) {
      role <- roles[[role_index]]
      payload <- direction_payload(snapshot, role)
      true_loading <- if (role == "A") {
        truth_data$true_params$a_true[, 1L]
      } else {
        study <- if (role == "B1") 1L else 2L
        truth_data$true_params$b_true[[study]][, 1L]
      }
      loading_rms <- vector_rms(payload$loading)
      process_rms <- vector_rms(payload$process)
      timeline_rows[[length(timeline_rows) + 1L]] <- data.frame(
        metadata, cumulative_t1_sweeps = sweep, block_role = role,
        factor_ppi = factor_ppi[[role_index]],
        loading_norm_ratio =
          vector_l2(payload$loading) / vector_l2(true_loading),
        loading_abs_cosine = abs_cosine(payload$loading, true_loading),
        feature_grid_rms = vector_rms(payload$feature),
        score_mean_rms = vector_rms(unlist(payload$score, use.names = FALSE)),
        factor_process_grid_rms = process_rms,
        complete_contribution_grid_rms_proxy = process_rms * loading_rms,
        stringsAsFactors = FALSE
      )
    }
  }
  rm(source_fit, continuation_fit, source_evaluation,
     continuation_evaluation, truth, source_directions,
     continuation_directions)
  invisible(gc(verbose = FALSE))
}

endpoint <- do.call(rbind, endpoint_rows)
changes <- do.call(rbind, change_rows)
gates <- do.call(rbind, gate_rows)
timeline <- do.call(rbind, timeline_rows)
process <- do.call(rbind, process_rows)
actions <- do.call(rbind, evaluation_action_rows)
write_csv(endpoint, "SOURCE_AND_800_ENDPOINT_METRICS.csv")
write_csv(changes, "SOURCE_TO_800_METRIC_CHANGES.csv")
write_csv(gates, "SOURCE_AND_800_GATE_AUDIT.csv")
write_csv(timeline, "CONTINUATION_MECHANISM_TIMELINE.csv")
write_csv(process, "SOURCE_AND_800_PROCESS_CONTRIBUTION.csv")
write_csv(actions, "EVALUATION_ACTIONS.csv")

process_summary <- do.call(rbind, lapply(
  split(process, interaction(
    process$continuation_id, process$cumulative_t1_sweeps,
    drop = TRUE
  )), function(group) {
    data.frame(
      continuation_id = group$continuation_id[[1L]],
      source_fit_id = group$source_fit_id[[1L]],
      data_id = group$data_id[[1L]],
      cumulative_t1_sweeps = group$cumulative_t1_sweeps[[1L]],
      study_role_matched = sum(group$matched_by_loading),
      study_role_truth_total = nrow(group),
      factor_process_nrmse_mean_missing_as_zero =
        mean(group$factor_process_nrmse_missing_as_zero),
      complete_contribution_nrmse_mean_missing_as_zero =
        mean(group$complete_contribution_nrmse_missing_as_zero),
      mean_loading_norm_ratio_matched =
        safe_mean(group$loading_norm_ratio),
      stringsAsFactors = FALSE
    )
  }
))
rownames(process_summary) <- NULL
write_csv(process_summary, "SOURCE_AND_800_PROCESS_CONTRIBUTION_SUMMARY.csv")

qc <- data.frame(
  check_id = c(
    "two_continuations_evaluated", "source_and_800_rows",
    "four_truth_roles_each_endpoint", "checkpoint_timeline_complete",
    "gate_comparison_complete", "truth_not_used_for_continuation",
    "selection_not_reopened", "no_automatic_further_extension",
    "development_not_formal"
  ),
  passed = c(
    nrow(actions) == 2L && all(actions$status == "COMPLETE"),
    nrow(endpoint) == 4L && all(table(endpoint$continuation_id) == 2L),
    nrow(process) == 16L &&
      all(table(interaction(
        process$continuation_id, process$cumulative_t1_sweeps
      )) == 4L),
    nrow(timeline) == 2L * 21L * 3L,
    nrow(gates) == 4L,
    !any(actions$truth_used_for_continuation) && !any(terminals$truth_used),
    !any(actions$winner_selection_reopened) &&
      !any(terminals$winner_selection_reopened),
    !any(terminals$automatic_further_extension),
    !any(actions$formal_v0lv_result) && !any(terminals$formal_v0lv_result)
  ), stringsAsFactors = FALSE
)
write_csv(qc, "EVALUATION_QC.csv")
assert(all(qc$passed), "Stage-4G evaluation QC failed.")
writeLines(c(
  "status=PASS", paste0("completed_utc=", iso_time()),
  "continuations_evaluated=2", "source_and_800_endpoints=4",
  "truth_study_role_rows=16", "winner_selection_reopened=FALSE",
  "automatic_further_extension=FALSE", "formal_v0lv_result=FALSE"
), file.path(output, "EVALUATION_COMPLETE.txt"), useBytes = TRUE)
utils::capture.output(
  sessionInfo(), file = file.path(output, "SESSION_INFO.txt")
)
cat("STAGE4G_EVALUATION_PASS continuations=2 endpoints=4\n")
