#!/usr/bin/env Rscript

# Post-selection evaluation for the paired population-mean bridge.  All fit
# endpoints and within-data-by-mean-arm ELBO winners must be frozen before this
# script reads any sealed truth.  No fit or continuation is started here.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
support_repo <- Sys.getenv("MEAN_BRIDGE_EVALUATION_SUPPORT_REPO", unset = "")
support_v0lv <- if (nzchar(support_repo)) {
  file.path(support_repo, "experiments", "v0L-V")
} else file.path(script_root, "..")
stage6a_support_path <- file.path(
  support_v0lv, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R")
stage6h_support_path <- file.path(
  support_v0lv, "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R")
source(stage6a_support_path, local = FALSE)
source(stage6h_support_path, local = FALSE)
if (nzchar(support_repo)) {
  PP_SNAPSHOT <- file.path(support_repo, "experiments", "v0L-V",
                           "rc1_snapshot")
}

PP_EXPERIMENT_ID <- "BSPLINE_MEAN_PAIRED_BRIDGE_V1_20260926"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% "1")
pp_assert(nzchar(output_root) && dir.exists(output_root) && workers %in% 1:8,
          "A prepared --output-root and 1--8 workers are required.")
pp_assert(pp_bool(cli$`truth-unseal-authorized`, FALSE),
          "Scientific evaluation requires explicit truth-unseal authorization.")
pp_assert(
  identical(tolower(pp_sha256(stage6a_support_path)),
            "bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6") &&
    identical(tolower(pp_sha256(stage6h_support_path)),
              "5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221") &&
    dir.exists(PP_SNAPSHOT),
  "The bound Stage 6A/6H evaluation support or frozen snapshot changed.")

evaluation_root <- file.path(
  output_root, "scientific_evaluation_20260927_v1")
pp_assert(!dir.exists(evaluation_root),
          "Evaluation output already exists; refusing overwrite.")

data_manifest_path <- file.path(output_root, "DATA_MANIFEST.csv")
fit_manifest_path <- file.path(output_root, "FIT_MANIFEST.csv")
selection_root <- file.path(output_root, "truth_free_selection")
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "TRUTH_FREE_WINNERS.csv")
selection_status_path <- file.path(selection_root, "SELECTION_STATUS.csv")
frozen_path <- file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt")
required <- c(
  data_manifest_path, fit_manifest_path, endpoints_path, winners_path,
  selection_status_path, frozen_path,
  file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"),
  file.path(output_root, "SUPERVISOR_EXIT_STATUS.txt"))
for (path in required) {
  pp_assert(file.exists(path), paste0("Missing required artifact: ", path))
}

data_manifest <- pp_read_csv(data_manifest_path)
manifest <- pp_read_csv(fit_manifest_path)
endpoints <- pp_read_csv(endpoints_path)
winners <- pp_read_csv(winners_path)
selection_status <- pp_read_csv(selection_status_path)
frozen <- readLines(frozen_path, warn = FALSE)
expected_data_ids <- as.vector(outer(
  c("meanb_01", "meanb_02"),
  c("periodic_phase", "independent_splinet"), paste, sep = "__"))

pp_assert(
  nrow(data_manifest) == 4L && nrow(manifest) == 48L &&
    nrow(endpoints) == 48L && nrow(winners) == 4L &&
    nrow(selection_status) == 4L &&
    all(data_manifest$experiment_id == PP_EXPERIMENT_ID) &&
    all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
    setequal(data_manifest$data_id, expected_data_ids) &&
    setequal(manifest$mean_arm,
             c("periodic_phase", "independent_splinet")) &&
    all(table(manifest$logical_data_id, manifest$mean_arm) == 12L) &&
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
    !any(manifest$automatic_800) && all(manifest$n_cpus == 1L) &&
    all(selection_status$n_registered == 12L) &&
    all(selection_status$n_eligible == 12L) &&
    all(selection_status$winner_status == "SELECTED"),
  "Manifest, endpoint, budget, or truth-isolation contract failed.")

pp_assert(all(c(
  "status=TRUTH_FREE_SELECTION_FROZEN",
  "rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_mean_arm",
  "cross_mean_arm_ELBO_comparison=FALSE",
  "truth_used_for_selection=FALSE", "truth_unseal_authorized=FALSE",
  "registered_fits=48", "selected_winners=4", "smoke=FALSE"
) %in% frozen), "Frozen selection marker is incomplete or changed.")

for (stratum in unique(manifest$selection_stratum_id)) {
  group <- endpoints[
    endpoints$selection_stratum_id == stratum &
      endpoints$objective_eligible &
      is.finite(endpoints$final_ordinary_elbo), , drop = FALSE]
  chosen <- winners[winners$selection_stratum_id == stratum, , drop = FALSE]
  expected <- group[order(-group$final_ordinary_elbo, group$fit_id), ,
                    drop = FALSE][1L, , drop = FALSE]
  pp_assert(nrow(group) == 12L && nrow(chosen) == 1L &&
              identical(chosen$fit_id[[1L]], expected$fit_id[[1L]]),
            paste0("Winner is not maximum eligible ELBO in ", stratum, "."))
}

fit_complete <- file.path(output_root, "fits", manifest$fit_id,
                          "FIT_COMPLETE.txt")
fit_error <- file.path(output_root, "fits", manifest$fit_id,
                       "FIT_ERROR.txt")
pp_assert(all(file.exists(fit_complete)) && !any(file.exists(fit_error)),
          "Fit completion accounting failed.")

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
  "authorization_date=2026-09-27",
  paste0("frozen_endpoints_sha256=", pp_sha256(endpoints_path)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  paste0("evaluator_sha256=", pp_sha256(script_file)),
  paste0("stage6a_support_sha256=", pp_sha256(stage6a_support_path)),
  paste0("stage6h_support_sha256=", pp_sha256(stage6h_support_path)),
  "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_mean_arm",
  "cross_mean_arm_ELBO_comparison=FALSE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE"
), file.path(evaluation_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))

# First sealed-truth read: all endpoint, selection and provenance checks above
# have already completed.
data_ids <- unique(manifest$data_id)
truths <- setNames(lapply(data_ids, function(data_id) {
  truth <- readRDS(file.path(
    output_root, "data", data_id, "sealed_truth", "truth_bundle.rds"))
  row <- data_manifest[data_manifest$data_id == data_id, , drop = FALSE]
  is_b <- identical(row$mean_arm[[1L]], "independent_splinet")
  mean_identity <- if (is_b) {
    identical(truth$data$true_params$mean_truth_method,
              "independent_splinet") &&
      identical(as.integer(truth$data$true_params$mean_basis_dim), 8L)
  } else !"mean_truth_method" %in% names(truth$data$true_params)
  pp_assert(identical(truth$bundle_class,
                      "bspline_mean_bridge_sealed_truth") &&
              identical(truth$experiment_id, PP_EXPERIMENT_ID) &&
              identical(truth$data_id, data_id) &&
              identical(truth$mean_arm, row$mean_arm[[1L]]) &&
              isFALSE(truth$truth_unseal_authorized) && mean_identity &&
              identical(truth$data$true_params$truth_basis_method,
                        "iid_common_splinet_stiefel") &&
              identical(as.integer(
                truth$data$true_params$truth_basis_dim), 8L),
            paste0("Sealed truth identity mismatch: ", data_id))
  truth
}), data_ids)

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

mean_bridge_weights <- function(grid) {
  delta <- diff(as.numeric(grid))
  c(delta[[1L]] / 2,
    (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}

mean_bridge_design <- function(grid, truth) {
  cbind(1, as.numeric(grid), multiFSYNC:::ZOSull(
    as.numeric(grid), range.x = c(0, 1), intKnots = truth$int_knots))
}

mean_bridge_true_dense <- function(data) {
  truth <- data$true_params
  if (identical(truth$mean_truth_method, "independent_splinet")) {
    return(truth$mean_dense_list)
  }
  lapply(seq_len(truth$S), function(study) {
    vapply(seq_len(truth$p), function(variable) {
      0.6 * (sin(2 * pi * truth$t_grid_dense +
                   truth$phase_mu_mat[study, variable]) +
               0.5 * sin(4 * pi * truth$t_grid_dense +
                           2 * truth$phase_mu_mat[study, variable]))
    }, numeric(length(truth$t_grid_dense)))
  })
}

mean_bridge_estimated_dense <- function(fit, data, study) {
  coefficient <- do.call(cbind, fit$mu_q_nu_mu_original[[study]])
  design <- mean_bridge_design(data$true_params$t_grid_dense,
                               data$true_params)
  pp_assert(nrow(coefficient) == ncol(design),
            "Mean coefficient and fitted-basis dimensions differ.")
  design %*% coefficient
}

mean_bridge_mean_recovery <- function(fit, data) {
  truth <- data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  weights <- mean_bridge_weights(grid)
  design <- mean_bridge_design(grid, truth)
  weighted_design <- design * sqrt(weights)
  true_dense <- mean_bridge_true_dense(data)
  rows <- lapply(seq_len(truth$S), function(study) {
    estimate <- mean_bridge_estimated_dense(fit, data, study)
    target <- true_dense[[study]]
    projection_coefficient <- qr.solve(
      weighted_design, target * sqrt(weights))
    projection <- design %*% projection_coefficient
    ise <- colSums(weights * (estimate - target)^2)
    energy <- colSums(weights * target^2)
    floor <- colSums(weights * (projection - target)^2)
    estimate_projection <- colSums(weights * (estimate - projection)^2)
    observed_sse <- observed_energy <- observed_n <- 0
    coefficient <- do.call(cbind, fit$mu_q_nu_mu_original[[study]])
    for (subject in seq_len(truth$n_s[[study]])) {
      time <- data$time_obs[[study]][[subject]]
      observed_estimate <- mean_bridge_design(time, truth) %*% coefficient
      observed_truth <- truth$mu_true_values[[study]][[subject]]
      observed_sse <- observed_sse +
        sum((observed_estimate - observed_truth)^2)
      observed_energy <- observed_energy + sum(observed_truth^2)
      observed_n <- observed_n + length(observed_truth)
    }
    data.frame(
      study = study, variables = truth$p,
      observed_scalar_count = observed_n,
      observed_sse = observed_sse, observed_truth_energy = observed_energy,
      observed_rmse = sqrt(observed_sse / observed_n),
      observed_nrmse = sqrt(observed_sse / max(observed_energy, 1e-12)),
      dense_ise_sum = sum(ise), dense_truth_energy = sum(energy),
      dense_rmse = sqrt(sum(ise) / (truth$p * sum(weights))),
      dense_nrmse = sqrt(sum(ise) / max(sum(energy), 1e-12)),
      dense_ise_mean = mean(ise),
      dense_relative_ise_mean = mean(ise / pmax(energy, 1e-12)),
      projection_floor_ise_sum = sum(floor),
      projection_floor_ise_mean = mean(floor),
      projection_relative_ise_mean =
        mean(floor / pmax(energy, 1e-12)),
      excess_ise_sum = sum(ise - floor),
      excess_ise_mean = mean(ise - floor),
      estimate_to_projection_ise_sum = sum(estimate_projection),
      estimate_to_projection_ise_mean = mean(estimate_projection),
      projection_decomposition_max_abs_error =
        max(abs(ise - floor - estimate_projection)),
      stringsAsFactors = FALSE)
  })
  by_study <- do.call(rbind, rows)
  total_variables <- sum(by_study$variables)
  summary <- data.frame(
    mean_observed_rmse = sqrt(sum(by_study$observed_sse) /
                                sum(by_study$observed_scalar_count)),
    mean_observed_nrmse = sqrt(sum(by_study$observed_sse) /
                                 max(sum(by_study$observed_truth_energy),
                                     1e-12)),
    mean_dense_rmse = sqrt(sum(by_study$dense_ise_sum) /
                             (total_variables * sum(weights))),
    mean_dense_nrmse = sqrt(sum(by_study$dense_ise_sum) /
                              max(sum(by_study$dense_truth_energy), 1e-12)),
    mean_dense_ise_mean = sum(by_study$dense_ise_sum) / total_variables,
    mean_dense_relative_ise_mean =
      mean(by_study$dense_relative_ise_mean),
    mean_projection_floor_ise_mean =
      sum(by_study$projection_floor_ise_sum) / total_variables,
    mean_projection_relative_ise_mean =
      mean(by_study$projection_relative_ise_mean),
    mean_excess_ise_mean = sum(by_study$excess_ise_sum) / total_variables,
    mean_estimate_to_projection_ise_mean =
      sum(by_study$estimate_to_projection_ise_sum) / total_variables,
    mean_projection_decomposition_max_abs_error =
      max(by_study$projection_decomposition_max_abs_error),
    stringsAsFactors = FALSE)
  list(by_study = by_study, summary = summary)
}

# The frozen v0J dense evaluator assumed the historical analytic mean.  Install
# a bridge-local replacement that changes only true/estimated mean assembly;
# all factor, loading, score and identity calculations remain frozen.
mean_bridge_install_dense_override <- function(environment) {
  environment$mean_bridge_true_dense <- mean_bridge_true_dense
  environment$mean_bridge_estimated_dense <- mean_bridge_estimated_dense
  dense_function <- function(fit, data) {
    truth <- data$true_params
    grid <- as.numeric(truth$t_grid_dense)
    true_L_s <- if (!is.null(truth$L_s_by_study)) {
      as.integer(unlist(truth$L_s_by_study, use.names = FALSE))
    } else if (length(truth$L_s) == truth$S) {
      as.integer(unlist(truth$L_s, use.names = FALSE))
    } else rep.int(as.integer(truth$L_s), truth$S)
    if (length(true_L_s) != truth$S || any(!is.finite(true_L_s))) {
      stop("Could not recover study-specific true factor counts.")
    }
    true_means <- mean_bridge_true_dense(data)
    shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
    specific_mask <- lapply(
      fit$factor_ppi_specific, function(value) as.numeric(value >= 0.5))
    accumulators <- list(
      all = c(sse = 0, truth = 0, n = 0),
      selected = c(sse = 0, truth = 0, n = 0),
      mean_only = c(sse = 0, truth = 0, n = 0))
    subject_offset <- 0L
    for (study in seq_len(truth$S)) {
      estimated_mean <- mean_bridge_estimated_dense(fit, data, study)
      true_mean <- true_means[[study]]
      estimated_shared_functions <- lapply(seq_len(fit$L_f), function(factor) {
        v0h_interpolate_matrix(fit$time_g, fit$list_Phi_hat[[factor]], grid)
      })
      estimated_specific_functions <- lapply(
        seq_along(fit$list_Phi_hat_spec[[study]]), function(factor) {
          v0h_interpolate_matrix(
            fit$time_g, fit$list_Phi_hat_spec[[study]][[factor]], grid)
        })
      for (subject in seq_len(truth$n_s[[study]])) {
        global_subject <- subject_offset + subject
        estimated_shared_path <- matrix(0, length(grid), fit$L_f)
        for (factor in seq_len(fit$L_f)) {
          estimated_shared_path[, factor] <-
            estimated_shared_functions[[factor]] %*%
            fit$list_Zeta_hat[[factor]][global_subject, ]
        }
        estimated_specific_path <- matrix(
          0, length(grid), length(estimated_specific_functions))
        if (ncol(estimated_specific_path)) {
          for (factor in seq_len(ncol(estimated_specific_path))) {
            estimated_specific_path[, factor] <-
              estimated_specific_functions[[factor]] %*%
              fit$list_Zeta_hat_spec[[study]][[factor]][subject, ]
          }
        }
        true_shared_path <- matrix(0, length(grid), truth$L_f)
        for (factor in seq_len(truth$L_f)) {
          true_shared_path[, factor] <-
            truth$theta_dense_list[[factor]] %*%
            truth$eta_true[[study]][[factor]][subject, ]
        }
        true_specific_path <- matrix(0, length(grid), true_L_s[[study]])
        if (ncol(true_specific_path)) {
          for (factor in seq_len(ncol(true_specific_path))) {
            true_specific_path[, factor] <-
              truth$kappa_dense_list[[study]][[factor]] %*%
              truth$chi_true[[study]][[factor]][subject, ]
          }
        }
        true_signal <- true_mean +
          true_shared_path %*% t(truth$a_true) +
          true_specific_path %*% t(truth$b_true[[study]])
        estimate_all <- estimated_mean +
          estimated_shared_path %*% t(fit$mu_q_a_hat) +
          estimated_specific_path %*% t(fit$mu_q_b_specific_hat[[study]])
        estimate_selected <- estimated_mean +
          estimated_shared_path %*%
            t(sweep(fit$mu_q_a_hat, 2L, shared_mask, "*")) +
          estimated_specific_path %*%
            t(sweep(fit$mu_q_b_specific_hat[[study]], 2L,
                    specific_mask[[study]], "*"))
        estimates <- list(all = estimate_all, selected = estimate_selected,
                          mean_only = estimated_mean)
        for (name in names(estimates)) {
          error <- estimates[[name]] - true_signal
          accumulators[[name]][["sse"]] <-
            accumulators[[name]][["sse"]] + sum(error^2)
          accumulators[[name]][["truth"]] <-
            accumulators[[name]][["truth"]] + sum(true_signal^2)
          accumulators[[name]][["n"]] <-
            accumulators[[name]][["n"]] + length(true_signal)
        }
      }
      subject_offset <- subject_offset + truth$n_s[[study]]
    }
    rows <- lapply(names(accumulators), function(name) {
      value <- accumulators[[name]]
      data.frame(
        scope = name, grid_size = length(grid),
        squared_error = value[["sse"]],
        truth_squared_energy = value[["truth"]],
        scalar_count = value[["n"]],
        rmse = sqrt(value[["sse"]] / value[["n"]]),
        truth_rms = sqrt(value[["truth"]] / value[["n"]]),
        nrmse = sqrt(value[["sse"]] / max(value[["truth"]], 1e-12)),
        relative_mise = value[["sse"]] / max(value[["truth"]], 1e-12),
        stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }
  environment(dense_function) <- environment
  environment$v0j_dense_signal_metrics <- dense_function
  invisible(environment)
}
mean_bridge_install_dense_override(evaluator_environment)

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
    truth <- truths[[row$data_id[[1L]]]]
    result <- s6h_evaluate_fit(fit, truth, evaluator_environment)
    mean_recovery <- mean_bridge_mean_recovery(fit, truth$data)
    result$mean_recovery_by_study <- mean_recovery$by_study
    result$mean_recovery_summary <- mean_recovery$summary
    result$evaluation_contract$mean_arm <- row$mean_arm[[1L]]
    result$evaluation_contract$mean_bridge_dense_override <- TRUE
    result$evaluation_contract$oracle_L_M <- TRUE
    result
  }, phase = "bspline_mean_bridge_post_selection_scientific_evaluation")
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
    logical_data_id = row$logical_data_id[[1L]],
    data_index = row$data_index[[1L]], mean_arm = row$mean_arm[[1L]],
    fit_K = row$fit_K[[1L]], fit_config_id = row$fit_config_id[[1L]],
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
                   result$loading_scale_summary,
                   result$mean_recovery_summary, activity)
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
  evaluation_root, "ALL_MEAN_BRIDGE_SCIENTIFIC_RESULTS.csv"))
winner_science <- science[science$truth_free_winner, , drop = FALSE]
pp_write_csv(winner_science, file.path(
  evaluation_root, "TRUTH_FREE_WINNER_SCIENCE.csv"))

arm_summary <- do.call(rbind, lapply(split(science, science$mean_arm), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(
    mean_arm = x$mean_arm[[1L]], starts = nrow(x),
    logical_datasets = length(unique(x$logical_data_id)),
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
pp_write_csv(arm_summary, file.path(evaluation_root, "MEAN_ARM_SUMMARY.csv"))

pp_write_lines(c(
  "status=BSPLINE_MEAN_BRIDGE_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("fits_evaluated=", nrow(science)),
  paste0("winners_evaluated=", nrow(winner_science)),
  "truth_read_after_all_within_data_by_mean_arm_selection_freezes=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "cross_mean_arm_ELBO_comparison=FALSE",
  "new_fit_or_continuation_started=FALSE",
  "broad_generalization_claim=FALSE",
  paste0("completed_utc=", pp_iso_time())
), file.path(evaluation_root, "MEAN_BRIDGE_EVALUATION_COMPLETE.txt"))

cat("MEAN_BRIDGE_EVALUATION_COMPLETE fits=", nrow(science),
    " winners=", nrow(winner_science), "\n", sep = "")
