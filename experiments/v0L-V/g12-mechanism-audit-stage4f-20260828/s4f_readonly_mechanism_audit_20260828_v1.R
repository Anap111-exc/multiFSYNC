#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]),
  winslash = "/", mustWork = TRUE
)
root <- dirname(script_path)
output <- file.path(root, "analysis_20260828_v1")

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
stage4e_root <- env_or(
  "STAGE4F_STAGE4E_ROOT",
  "/root/v0lv-g12-focused-reachability-stage4e-20260827-v1"
)
stage4c_root <- env_or(
  "STAGE4F_STAGE4C_ROOT",
  "/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1"
)
snapshot_root <- env_or(
  "STAGE4F_SNAPSHOT_ROOT",
  "/root/multiFSYNC-git-rc1/experiments/v0L-V/rc1_snapshot"
)
development_library <- env_or(
  "STAGE4F_DEVELOPMENT_LIBRARY",
  "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5-lib"
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
safe_max <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) max(value) else NA_real_
}
first_sweep <- function(sweep, condition) {
  keep <- is.finite(sweep) & !is.na(condition) & condition
  if (any(keep)) min(sweep[keep]) else NA_integer_
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
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (length(left) != length(right) || !length(left) ||
      any(!is.finite(c(left, right)))) return(NA_real_)
  denominator <- vector_l2(left) * vector_l2(right)
  if (!is.finite(denominator) || denominator <= 0) return(0)
  abs(sum(left * right) / denominator)
}
nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
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
       "Refusing to overwrite analysis_20260828_v1.")
required_roots <- c(stage4e_root, stage4c_root, snapshot_root)
assert(all(dir.exists(required_roots)), "A bound Stage-4C/4E root is missing.")
required_markers <- c(
  file.path(stage4e_root, "STAGE4E_COMPLETE.txt"),
  file.path(stage4e_root, "TRUTH_FREE_SELECTION_COMPLETE.txt"),
  file.path(stage4e_root, "truth_evaluation_20260827_v1",
            "EVALUATION_COMPLETE.txt"),
  file.path(stage4e_root, "analysis_20260827_v1",
            "ANALYSIS_COMPLETE.txt")
)
assert(all(file.exists(required_markers)), "Stage 4E is not complete.")

inventory <- read_csv(file.path(
  stage4e_root, "ALL_24_TRUTH_FREE_ENDPOINTS.csv"
))
winners <- read_csv(file.path(
  stage4e_root, "TRUTH_FREE_ELBO_WINNERS_8_STARTS.csv"
))
reachability <- read_csv(file.path(
  stage4e_root, "analysis_20260827_v1", "PER_ENDPOINT_REACHABILITY.csv"
))
actions <- read_csv(file.path(
  stage4e_root, "truth_evaluation_20260827_v1", "EVALUATION_ACTIONS.csv"
))
assert(
  nrow(inventory) == 24L && all(table(inventory$data_id) == 8L) &&
    all(inventory$objective_eligible) && !any(inventory$truth_used) &&
    !any(inventory$continuation_used) &&
    !any(inventory$automatic_800_used),
  "The Stage-4E endpoint inventory is not the frozen 24-endpoint input."
)
assert(
  nrow(winners) == 3L && !any(winners$selection_used_truth) &&
    !any(winners$selection_used_structure) &&
    !any(winners$selection_used_NRMSE_ISE_RPL),
  "The Stage-4E truth-free winner table is invalid."
)
assert(nrow(reachability) == 24L && nrow(actions) == 24L &&
         all(actions$status == "COMPLETE"),
       "The Stage-4E truth evaluation is incomplete.")
assert(setequal(inventory$fit_id, reachability$fit_id) &&
         setequal(inventory$fit_id, actions$fit_id),
       "Stage-4E fit identities do not join exactly.")

if (dir.exists(development_library)) {
  .libPaths(unique(c(development_library, .libPaths())))
}
wrapper <- find_one(snapshot_root, "^v0lv_candidate2_evaluation[.]R$")
sys.source(wrapper, envir = globalenv(), keep.source = TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(snapshot_root)

dir.create(output, recursive = TRUE, mode = "0700")

fit_path_for <- function(row) {
  if (identical(row$source_stage[[1L]], "stage4c_inherited")) {
    file.path(stage4c_root, "fits", row$fit_id[[1L]], "fit.rds")
  } else {
    file.path(stage4e_root, "fits", row$fit_id[[1L]], "fit.rds")
  }
}
direction_path_for <- function(row) {
  if (identical(row$source_stage[[1L]], "stage4c_inherited")) {
    file.path(stage4c_root, "fits", row$fit_id[[1L]],
              "DIRECTION_SNAPSHOTS.rds")
  } else {
    file.path(stage4e_root, "fits", row$fit_id[[1L]],
              "DIRECTION_SNAPSHOTS.rds")
  }
}
evaluation_path_for <- function(fit_id) {
  file.path(stage4e_root, "truth_evaluation_20260827_v1",
            "endpoints", fit_id, "evaluation.rds")
}
truth_path_for <- function(data_id) {
  file.path(stage4c_root, "data", data_id,
            "sealed_truth", "truth_bundle.rds")
}
endpoint_role_for <- function(fit_id) {
  if (fit_id %in% winners$selected_fit_id) "winner" else "not_selected"
}

checkpoint_roles <- c("A", "B1", "B2")
ppi_rows <- list()
direction_rows <- list()
global_scale_rows <- list()
gate_rows <- list()
final_five_rows <- list()
process_rows <- list()
process_endpoint_rows <- list()

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
    kernel <- lapply(snapshot$shared_score_kernel, `[[`, 1L)
  } else {
    study <- if (role == "B1") 1L else 2L
    loading <- snapshot$specific_loading[[study]][, 1L]
    feature <- snapshot$specific_feature[[study]][[1L]]
    score <- list(snapshot$specific_score_mean[[study]][[1L]])
    kernel <- list(snapshot$specific_score_kernel[[study]][[1L]])
  }
  process <- unlist(lapply(score, function(score_matrix) {
    as.numeric(as.matrix(score_matrix) %*% t(as.matrix(feature)))
  }), use.names = FALSE)
  list(
    loading = as.numeric(loading), feature = as.matrix(feature),
    score = score, kernel = kernel, process = process
  )
}

for (index in seq_len(nrow(inventory))) {
  row <- inventory[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  data_id <- row$data_id[[1L]]
  fit_path <- fit_path_for(row)
  direction_path <- direction_path_for(row)
  evaluation_path <- evaluation_path_for(fit_id)
  truth_path <- truth_path_for(data_id)
  assert(all(file.exists(c(fit_path, direction_path,
                           evaluation_path, truth_path))),
         paste0("A Stage-4E audit input is missing for ", fit_id, "."))

  fit <- readRDS(fit_path)
  directions <- readRDS(direction_path)
  evaluation <- readRDS(evaluation_path)
  truth <- readRDS(truth_path)
  assert(is.list(fit) && is.list(directions) && is.list(evaluation) &&
           identical(truth$bundle_class, "v0lv_sealed_truth_candidate2") &&
           identical(truth$data_id, data_id),
         paste0("Invalid fit/evaluation/truth identity for ", fit_id, "."))
  assert(identical(as.integer(fit$t1_sweeps), 400L) &&
           length(fit$ELBO) == 400L &&
           identical(names(directions),
                     as.character(c(0L, seq.int(20L, 400L, 20L)))) &&
           identical(names(fit$practical_checkpoints),
                     as.character(seq.int(20L, 400L, 20L))),
         paste0("Incomplete fixed-400 checkpoint path for ", fit_id, "."))

  structure_row <- reachability[reachability$fit_id == fit_id, , drop = FALSE]
  assert(nrow(structure_row) == 1L, "Reachability join is not one-to-one.")
  endpoint_role <- endpoint_role_for(fit_id)
  truth_data <- truth$data
  p <- nrow(truth_data$true_params$a_true)

  for (checkpoint_name in names(fit$practical_checkpoints)) {
    checkpoint <- fit$practical_checkpoints[[checkpoint_name]]
    snapshot <- checkpoint$snapshot
    factor_ppi <- as.numeric(snapshot$factor_ppi)
    variable_ppi <- as.numeric(snapshot$ppi)
    assert(length(factor_ppi) == 3L && length(variable_ppi) == 3L * p,
           paste0("Unexpected PPI checkpoint shape for ", fit_id, "."))
    sweep <- as.integer(checkpoint_name)
    for (role_index in seq_along(checkpoint_roles)) {
      role <- checkpoint_roles[[role_index]]
      segment <- seq.int((role_index - 1L) * p + 1L, role_index * p)
      values <- variable_ppi[segment]
      ppi_rows[[length(ppi_rows) + 1L]] <- data.frame(
        fit_id = fit_id, data_id = data_id,
        scenario_id = row$scenario_id[[1L]],
        seed_index = row$seed_index[[1L]],
        source_stage = row$source_stage[[1L]],
        endpoint_role = endpoint_role,
        structure_exact = structure_row$structure_exact[[1L]],
        structure_class = structure_row$structure_class[[1L]],
        t1_sweep = sweep, block_role = role,
        factor_ppi = factor_ppi[[role_index]],
        variable_ppi_min = min(values),
        variable_ppi_mean = mean(values),
        variable_ppi_median = stats::median(values),
        variable_ppi_max = max(values),
        variable_ppi_gt_half = sum(values > 0.5),
        variable_count = length(values), stringsAsFactors = FALSE
      )
    }
  }

  for (checkpoint_name in names(directions)) {
    snapshot <- directions[[checkpoint_name]]
    sweep <- as.integer(checkpoint_name)
    for (role in checkpoint_roles) {
      payload <- direction_payload(snapshot, role)
      true_loading <- if (role == "A") {
        truth_data$true_params$a_true[, 1L]
      } else {
        study <- if (role == "B1") 1L else 2L
        truth_data$true_params$b_true[[study]][, 1L]
      }
      loading_rms <- vector_rms(payload$loading)
      process_rms <- vector_rms(payload$process)
      score_values <- unlist(payload$score, use.names = FALSE)
      kernel_values <- unlist(payload$kernel, use.names = FALSE)
      true_norm <- vector_l2(true_loading)
      estimated_norm <- vector_l2(payload$loading)
      direction_rows[[length(direction_rows) + 1L]] <- data.frame(
        fit_id = fit_id, data_id = data_id,
        scenario_id = row$scenario_id[[1L]],
        seed_index = row$seed_index[[1L]],
        source_stage = row$source_stage[[1L]],
        endpoint_role = endpoint_role,
        structure_exact = structure_row$structure_exact[[1L]],
        structure_class = structure_row$structure_class[[1L]],
        t1_sweep = sweep, block_role = role,
        loading_l2 = estimated_norm,
        loading_rms = loading_rms,
        true_loading_l2 = true_norm,
        same_role_loading_norm_ratio = estimated_norm / true_norm,
        same_role_loading_abs_cosine =
          abs_cosine(payload$loading, true_loading),
        feature_grid_rms = vector_rms(payload$feature),
        score_mean_rms = vector_rms(score_values),
        score_kernel_rms = vector_rms(kernel_values),
        factor_process_grid_rms = process_rms,
        complete_contribution_grid_rms_proxy = process_rms * loading_rms,
        stringsAsFactors = FALSE
      )
    }
  }

  scale <- fit$scale_trace
  required_scale <- c(
    "stage", "t1_sweep", "elbo", "observed_rms", "fitted_rms",
    "mean_contribution_rms", "shared_contribution_rms",
    "specific_contribution_rms", "expected_rss_per_obs",
    "mean_residual_sse", "noise_precision_median", "phi_rms",
    "psi_rms", "zeta_rms", "xi_rms", "a_rms", "b_rms",
    "shared_effective_loading_q_second_moment",
    "specific_effective_loading_q_second_moment",
    "shared_ppi_mean", "specific_ppi_mean"
  )
  assert(is.data.frame(scale) && all(required_scale %in% names(scale)) &&
           nrow(scale) == 21L,
         paste0("Unexpected scale trace for ", fit_id, "."))
  scale_out <- scale[, required_scale, drop = FALSE]
  scale_out <- cbind(data.frame(
    fit_id = fit_id, data_id = data_id,
    scenario_id = row$scenario_id[[1L]],
    seed_index = row$seed_index[[1L]],
    source_stage = row$source_stage[[1L]],
    endpoint_role = endpoint_role,
    structure_exact = structure_row$structure_exact[[1L]],
    structure_class = structure_row$structure_class[[1L]],
    stringsAsFactors = FALSE
  ), scale_out)
  global_scale_rows[[length(global_scale_rows) + 1L]] <- scale_out

  diagnostics <- fit$practical_diagnostics
  assert(is.data.frame(diagnostics) && nrow(diagnostics) == 400L &&
           tail(diagnostics$t1_sweep, 1L) == 400L,
         paste0("Unexpected practical diagnostics for ", fit_id, "."))
  decomposed <- decompose_gate(diagnostics, fit$practical_control)
  final <- tail(decomposed, 1L)
  component_names <- c(
    "objective_rate_pass", "objective_monotone_pass",
    "short_fitted_pass", "short_rss_pass",
    "short_ppi_quantile_pass", "short_factor_ppi_pass",
    "long_fitted_pass", "long_rss_pass",
    "long_ppi_quantile_pass", "long_factor_ppi_pass"
  )
  failed <- component_names[!unlist(final[1L, component_names])]
  gate_rows[[length(gate_rows) + 1L]] <- data.frame(
    fit_id = fit_id, data_id = data_id,
    scenario_id = row$scenario_id[[1L]],
    seed_index = row$seed_index[[1L]],
    source_stage = row$source_stage[[1L]],
    endpoint_role = endpoint_role,
    structure_exact = structure_row$structure_exact[[1L]],
    structure_class = structure_row$structure_class[[1L]],
    practical_converged = row$practical_converged[[1L]],
    slow_case = row$slow_case[[1L]],
    final_streak = final$streak[[1L]],
    final_pass = final$pass[[1L]],
    final_objective_pass = final$objective_pass[[1L]],
    final_output_pass = final$output_pass[[1L]],
    primary_failed_components = paste(failed, collapse = ";"),
    final[, c(
      "elbo_abs_rate", "elbo_per_response_rate", "elbo_rel_rate",
      "fitted_nrmse", "rss_rel", "ppi_quantile_abs",
      "factor_ppi_max_abs", "long_fitted_nrmse", "long_rss_rel",
      "long_ppi_quantile_abs", "long_factor_ppi_max_abs",
      component_names
    ), drop = FALSE], stringsAsFactors = FALSE
  )
  last_five <- decomposed[decomposed$t1_sweep >= 396L, , drop = FALSE]
  assert(nrow(last_five) == 5L, "Final-five convergence window is missing.")
  final_five_rows[[length(final_five_rows) + 1L]] <- cbind(
    data.frame(
      fit_id = fit_id, data_id = data_id,
      scenario_id = row$scenario_id[[1L]],
      seed_index = row$seed_index[[1L]],
      endpoint_role = endpoint_role,
      structure_exact = structure_row$structure_exact[[1L]],
      stringsAsFactors = FALSE
    ),
    last_five[, c(
      "t1_sweep", "objective_pass", "short_output_pass",
      "long_output_pass", "output_pass", "pass", "streak",
      "elbo_abs_rate", "elbo_per_response_rate", "fitted_nrmse",
      "rss_rel", "ppi_quantile_abs", "factor_ppi_max_abs",
      "long_fitted_nrmse", "long_rss_rel",
      "long_ppi_quantile_abs", "long_factor_ppi_max_abs",
      component_names
    ), drop = FALSE]
  )

  reported_paths <-
    evaluator_environment$v0g_reported_paths_by_study(fit, truth_data)
  truth_paths <- evaluator_environment$v0h_truth_paths_by_study(truth_data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected", , drop = FALSE]
  assert(nrow(functional) == 4L,
         paste0("Expected four truth study-role rows for ", fit_id, "."))
  endpoint_process <- lapply(seq_len(nrow(functional)), function(row_index) {
    item <- functional[row_index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_block <- item$true_role[[1L]]
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_block == "shared") {
      truth_data$true_params$a_true[, true_factor]
    } else {
      truth_data$true_params$b_true[[study]][, true_factor]
    }
    true_path <- truth_paths[[true_block]][[study]][, true_factor]
    true_contribution <- as.numeric(outer(true_path, true_loading))
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
      aligned_path <- item$loading_sign[[1L]] * estimated_path
      process_error <- nrmse(aligned_path, true_path)
      contribution_error <- nrmse(
        as.numeric(outer(estimated_path, estimated_loading)),
        true_contribution
      )
      estimated_loading_norm <- vector_l2(estimated_loading)
      loading_norm_ratio <- estimated_loading_norm / vector_l2(true_loading)
      loading_direction_cosine <-
        abs_cosine(estimated_loading, true_loading)
    } else {
      process_error <- contribution_error <- NA_real_
      estimated_loading_norm <- loading_norm_ratio <-
        loading_direction_cosine <- NA_real_
    }
    data.frame(
      fit_id = fit_id, data_id = data_id,
      scenario_id = row$scenario_id[[1L]],
      seed_index = row$seed_index[[1L]],
      endpoint_role = endpoint_role,
      structure_exact = structure_row$structure_exact[[1L]],
      structure_class = structure_row$structure_class[[1L]],
      study = study, true_id = item$true_id[[1L]],
      true_role = true_block, true_factor = true_factor,
      candidate_id = item$candidate_id[[1L]],
      estimated_block = item$estimated_block[[1L]],
      estimated_factor = item$estimated_factor[[1L]],
      loading_structure_status = item$loading_structure_status[[1L]],
      matched_by_loading = matched,
      true_loading_l2 = vector_l2(true_loading),
      estimated_loading_l2 = estimated_loading_norm,
      loading_norm_ratio = loading_norm_ratio,
      loading_abs_cosine = loading_direction_cosine,
      factor_process_nrmse_matched = process_error,
      factor_process_nrmse_missing_as_zero =
        if (matched) process_error else 1,
      complete_contribution_nrmse_matched = contribution_error,
      complete_contribution_nrmse_missing_as_zero =
        if (matched) contribution_error else 1,
      complete_contribution_grid =
        "subject_observation_times_by_variable",
      missing_candidate_estimate = "zero_contribution",
      stringsAsFactors = FALSE
    )
  })
  endpoint_process <- do.call(rbind, endpoint_process)
  process_rows[[length(process_rows) + 1L]] <- endpoint_process
  correctly_placed <- endpoint_process$matched_by_loading &
    endpoint_process$loading_structure_status == "correct"
  process_endpoint_rows[[length(process_endpoint_rows) + 1L]] <- data.frame(
    fit_id = fit_id, data_id = data_id,
    scenario_id = row$scenario_id[[1L]],
    seed_index = row$seed_index[[1L]],
    endpoint_role = endpoint_role,
    structure_exact = structure_row$structure_exact[[1L]],
    structure_class = structure_row$structure_class[[1L]],
    study_role_truth_total = nrow(endpoint_process),
    study_role_matched = sum(endpoint_process$matched_by_loading),
    study_role_correctly_placed = sum(correctly_placed),
    factor_process_nrmse_mean_matched =
      safe_mean(endpoint_process$factor_process_nrmse_matched),
    factor_process_nrmse_mean_missing_as_zero =
      mean(endpoint_process$factor_process_nrmse_missing_as_zero),
    complete_contribution_nrmse_mean_matched = safe_mean(
      endpoint_process$complete_contribution_nrmse_matched
    ),
    complete_contribution_nrmse_mean_missing_as_zero = mean(
      endpoint_process$complete_contribution_nrmse_missing_as_zero
    ),
    complete_contribution_nrmse_mean_correctly_placed = safe_mean(
      endpoint_process$complete_contribution_nrmse_matched[correctly_placed]
    ), stringsAsFactors = FALSE
  )

  rm(fit, directions, evaluation, truth, truth_data)
  invisible(gc(verbose = FALSE))
}

ppi <- do.call(rbind, ppi_rows)
direction <- do.call(rbind, direction_rows)
global_scale <- do.call(rbind, global_scale_rows)
gate <- do.call(rbind, gate_rows)
final_five <- do.call(rbind, final_five_rows)
process <- do.call(rbind, process_rows)
process_endpoint <- do.call(rbind, process_endpoint_rows)

write_csv(ppi, "CHECKPOINT_BLOCK_PPI.csv")
write_csv(direction, "CHECKPOINT_BLOCK_SCALE.csv")
write_csv(global_scale, "CHECKPOINT_GLOBAL_SCALE.csv")
write_csv(gate, "PER_ENDPOINT_FINAL_GATE_AUDIT.csv")
write_csv(final_five, "FINAL_FIVE_SWEEP_GATE_AUDIT.csv")
write_csv(process, "ALL_PROCESS_AND_CONTRIBUTION_RECOVERY.csv")
write_csv(process_endpoint,
          "PER_ENDPOINT_PROCESS_CONTRIBUTION_SUMMARY.csv")

transition_keys <- unique(ppi[, c(
  "fit_id", "data_id", "scenario_id", "seed_index", "endpoint_role",
  "structure_exact", "structure_class", "block_role"
)])
transition_rows <- lapply(seq_len(nrow(transition_keys)), function(index) {
  key <- transition_keys[index, , drop = FALSE]
  ppi_group <- ppi[ppi$fit_id == key$fit_id &
                     ppi$block_role == key$block_role, , drop = FALSE]
  scale_group <- direction[direction$fit_id == key$fit_id &
                             direction$block_role == key$block_role,
                           , drop = FALSE]
  ppi_group <- ppi_group[order(ppi_group$t1_sweep), , drop = FALSE]
  scale_group <- scale_group[order(scale_group$t1_sweep), , drop = FALSE]
  data.frame(
    key,
    ppi_at_20 = ppi_group$factor_ppi[ppi_group$t1_sweep == 20L],
    ppi_at_200 = ppi_group$factor_ppi[ppi_group$t1_sweep == 200L],
    ppi_at_400 = ppi_group$factor_ppi[ppi_group$t1_sweep == 400L],
    first_checkpoint_ppi_below_0_5 = first_sweep(
      ppi_group$t1_sweep, ppi_group$factor_ppi < 0.5
    ),
    first_checkpoint_ppi_below_0_1 = first_sweep(
      ppi_group$t1_sweep, ppi_group$factor_ppi < 0.1
    ),
    first_checkpoint_ppi_below_0_01 = first_sweep(
      ppi_group$t1_sweep, ppi_group$factor_ppi < 0.01
    ),
    loading_norm_ratio_at_0 = scale_group$same_role_loading_norm_ratio[
      scale_group$t1_sweep == 0L
    ],
    loading_norm_ratio_at_200 = scale_group$same_role_loading_norm_ratio[
      scale_group$t1_sweep == 200L
    ],
    loading_norm_ratio_at_400 = scale_group$same_role_loading_norm_ratio[
      scale_group$t1_sweep == 400L
    ],
    max_loading_norm_ratio =
      safe_max(scale_group$same_role_loading_norm_ratio),
    first_checkpoint_loading_ratio_above_2 = first_sweep(
      scale_group$t1_sweep,
      scale_group$same_role_loading_norm_ratio > 2
    ),
    first_checkpoint_loading_ratio_above_5 = first_sweep(
      scale_group$t1_sweep,
      scale_group$same_role_loading_norm_ratio > 5
    ),
    same_role_loading_abs_cosine_at_400 =
      scale_group$same_role_loading_abs_cosine[
        scale_group$t1_sweep == 400L
      ],
    factor_process_grid_rms_at_400 = scale_group$factor_process_grid_rms[
      scale_group$t1_sweep == 400L
    ],
    complete_contribution_grid_rms_proxy_at_400 =
      scale_group$complete_contribution_grid_rms_proxy[
        scale_group$t1_sweep == 400L
      ], stringsAsFactors = FALSE
  )
})
transitions <- do.call(rbind, transition_rows)
write_csv(transitions, "BLOCK_TRANSITION_SUMMARY.csv")

winner_ids <- winners$selected_fit_id
winner_timeline <- merge(
  direction[direction$fit_id %in% winner_ids, ],
  ppi[ppi$fit_id %in% winner_ids,
      c("fit_id", "t1_sweep", "block_role", "factor_ppi",
        "variable_ppi_mean", "variable_ppi_max")],
  by = c("fit_id", "t1_sweep", "block_role"), sort = FALSE
)
winner_timeline <- winner_timeline[order(
  winner_timeline$data_id, winner_timeline$t1_sweep,
  winner_timeline$block_role
), , drop = FALSE]
write_csv(winner_timeline, "WINNER_MECHANISM_TIMELINE.csv")
write_csv(
  process_endpoint[process_endpoint$fit_id %in% winner_ids, , drop = FALSE],
  "WINNER_PROCESS_CONTRIBUTION_SUMMARY.csv"
)

blocker_names <- c(
  "objective_rate_pass", "objective_monotone_pass",
  "short_fitted_pass", "short_rss_pass",
  "short_ppi_quantile_pass", "short_factor_ppi_pass",
  "long_fitted_pass", "long_rss_pass",
  "long_ppi_quantile_pass", "long_factor_ppi_pass"
)
blocker_summary <- do.call(rbind, lapply(blocker_names, function(name) {
  data.frame(
    gate_component = name,
    failed_all_endpoints = sum(!gate[[name]]),
    failed_slow_endpoints = sum(!gate[[name]] & gate$slow_case),
    failed_converged_endpoints = sum(
      !gate[[name]] & gate$practical_converged
    ), stringsAsFactors = FALSE
  )
}))
write_csv(blocker_summary, "FINAL_GATE_BLOCKER_SUMMARY.csv")

group_names <- c("all", "winner", "exact", "slow")
group_indices <- list(
  rep(TRUE, nrow(process_endpoint)),
  process_endpoint$endpoint_role == "winner",
  process_endpoint$structure_exact,
  gate$slow_case[match(process_endpoint$fit_id, gate$fit_id)]
)
process_group <- do.call(rbind, Map(function(name, keep) {
  group <- process_endpoint[keep, , drop = FALSE]
  data.frame(
    group = name, endpoints = nrow(group),
    study_role_matched = sum(group$study_role_matched),
    study_role_truth_total = sum(group$study_role_truth_total),
    study_role_coverage =
      sum(group$study_role_matched) / sum(group$study_role_truth_total),
    mean_factor_process_nrmse_missing_as_zero = safe_mean(
      group$factor_process_nrmse_mean_missing_as_zero
    ),
    median_factor_process_nrmse_missing_as_zero = stats::median(
      group$factor_process_nrmse_mean_missing_as_zero
    ),
    mean_complete_contribution_nrmse_missing_as_zero = safe_mean(
      group$complete_contribution_nrmse_mean_missing_as_zero
    ),
    median_complete_contribution_nrmse_missing_as_zero = stats::median(
      group$complete_contribution_nrmse_mean_missing_as_zero
    ), stringsAsFactors = FALSE
  )
}, group_names, group_indices))
write_csv(process_group, "PROCESS_CONTRIBUTION_GROUP_SUMMARY.csv")

qc <- data.frame(
  check_id = c(
    "twenty_four_bound_endpoints", "three_frozen_truth_free_winners",
    "checkpoint_ppi_complete", "checkpoint_scale_complete",
    "global_scale_complete", "final_gate_complete",
    "final_five_complete", "four_truth_roles_per_endpoint",
    "process_contribution_complete", "winner_outputs_complete",
    "no_fit_or_continuation", "selection_not_reopened",
    "development_not_formal"
  ),
  passed = c(
    nrow(inventory) == 24L && all(table(inventory$data_id) == 8L),
    nrow(winners) == 3L &&
      sum(process_endpoint$endpoint_role == "winner") == 3L,
    nrow(ppi) == 24L * 20L * 3L,
    nrow(direction) == 24L * 21L * 3L,
    nrow(global_scale) == 24L * 21L,
    nrow(gate) == 24L,
    nrow(final_five) == 24L * 5L,
    all(process_endpoint$study_role_truth_total == 4L),
    nrow(process) == 24L * 4L,
    nrow(winner_timeline) == 3L * 20L * 3L &&
      nrow(process_endpoint[process_endpoint$endpoint_role == "winner", ]) == 3L,
    !any(inventory$continuation_used) &&
      !any(inventory$automatic_800_used),
    !any(winners$selection_used_truth) &&
      setequal(winner_ids, process_endpoint$fit_id[
        process_endpoint$endpoint_role == "winner"
      ]),
    !any(inventory$formal_v0lv_result)
  ), stringsAsFactors = FALSE
)
write_csv(qc, "AUDIT_QC.csv")
assert(all(qc$passed), "Stage-4F read-only mechanism audit QC failed.")

writeLines(c(
  "status=PASS", paste0("completed_utc=", iso_time()),
  "input_endpoints=24", "frozen_truth_free_winners=3",
  "checkpoint_ppi_rows=1440", "checkpoint_scale_rows=1512",
  "process_contribution_rows=96", "new_fit_run=FALSE",
  "continuation_run=FALSE", "selection_reopened=FALSE",
  "formal_v0lv_result=FALSE"
), file.path(output, "AUDIT_COMPLETE.txt"), useBytes = TRUE)
utils::capture.output(sessionInfo(), file = file.path(output, "SESSION_INFO.txt"))
cat("STAGE4F_READONLY_AUDIT_PASS endpoints=24 winners=3\n")
