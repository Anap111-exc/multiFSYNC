#!/usr/bin/env Rscript

# Post-selection Stage 7A evaluation. Never used for fitting or winner choice.
# The frozen v0J evaluator supplies the historical residual-factor diagnostics;
# covariate and full-signal metrics below explicitly add the beta contribution.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_root <- dirname(sub("^--file=", "", script_arg[[1L]]))
source(file.path(script_root, "..", "g12-dimension-overspec-stage6a-20260831",
                 "stage6a_common_20260831_v1.R"), local = FALSE)
source(file.path(script_root, "..", "g12-reachability-causal-audit-stage6h-20260913",
                 "stage6h_evaluation_helpers_20260913_v1.R"), local = FALSE)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% "1")
limit <- as.integer(cli$limit %||% "144")
evaluation_name <- cli$`evaluation-name` %||%
  "scientific_evaluation_20260920_v1"
pp_assert(nzchar(output_root) && dir.exists(output_root),
          "A prepared Stage 7A --output-root is required.")
pp_assert(pp_bool(cli$`truth-unseal-authorized`, FALSE),
          "Explicit truth-unseal authorization is required.")
pp_assert(workers %in% 1:8 && limit %in% 1:144 &&
            grepl("^[A-Za-z0-9_-]+$", evaluation_name),
          "Invalid evaluation workers, limit, or output name.")

read_csv <- function(path) {
  pp_assert(file.exists(path), paste0("Missing file: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv_new <- function(value, path) {
  pp_assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
write_lines_new <- function(value, path) {
  pp_assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
}
save_rds_new <- function(value, path) {
  pp_assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3)
}

selection_root <- file.path(output_root, "truth_free_selection")
frozen <- readLines(file.path(selection_root,
                             "TRUTH_FREE_SELECTION_FROZEN.txt"), warn = FALSE)
pp_assert(all(c("status=TRUTH_FREE_SELECTION_FROZEN",
                "truth_used_for_selection=FALSE",
                "selected_winners=12", "registered_fits=144",
                "smoke=FALSE") %in% frozen),
          "Truth-free selection is not completely frozen.")
pp_assert(file.exists(file.path(output_root,
                                "FITS_AND_SELECTION_COMPLETE.txt")),
          "The full main fit batch has not completed.")
manifest <- read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
data_manifest <- read_csv(file.path(output_root, "DATA_MANIFEST.csv"))
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "TRUTH_FREE_WINNERS.csv")
endpoints <- read_csv(endpoints_path)
winners <- read_csv(winners_path)
pp_assert(nrow(manifest) == 144L && nrow(data_manifest) == 12L &&
            nrow(endpoints) == 144L && nrow(winners) == 12L &&
            !anyDuplicated(manifest$fit_id) &&
            !anyDuplicated(winners$data_id) &&
            setequal(manifest$fit_id, endpoints$fit_id) &&
            setequal(data_manifest$data_id, winners$data_id) &&
            all(endpoints$status == "COMPLETE") &&
            all(endpoints$objective_eligible) &&
            all(endpoints$ordinary_t1_sweeps == 400L) &&
            all(endpoints$annealing_sweeps == 99L) &&
            !any(endpoints$truth_used_for_fit_or_selection),
          "Frozen manifest, endpoint, or truth-isolation contract failed.")
for (data_id in data_manifest$data_id) {
  group <- endpoints[endpoints$data_id == data_id, , drop = FALSE]
  chosen <- winners[winners$data_id == data_id, , drop = FALSE]
  pp_assert(nrow(group) == 12L && nrow(chosen) == 1L &&
              identical(chosen$fit_id[[1L]],
                        group$fit_id[[which.max(group$final_ordinary_elbo)]]),
            paste0("Frozen winner is not maximum eligible ELBO for ", data_id))
}

evaluation_root <- file.path(output_root, evaluation_name)
pp_assert(!dir.exists(evaluation_root),
          "Evaluation directory exists; refusing overwrite.")
dir.create(evaluation_root, recursive = TRUE, mode = "0700")
write_lines_new(c(
  "status=POST_SELECTION_TRUTH_UNSEAL_AUTHORIZED_BY_USER",
  "authorization_date=2026-09-20",
  paste0("frozen_endpoints_sha256=", pp_sha256(endpoints_path)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  "truth_used_for_fit_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE"
), file.path(evaluation_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))

# All manifest/selection checks above precede the first sealed-truth read.
manifest <- manifest[seq_len(limit), , drop = FALSE]
data_ids <- unique(manifest$data_id)
truths <- setNames(lapply(data_ids, function(data_id) {
  x <- readRDS(file.path(output_root, "data", data_id,
                         "sealed_truth", "truth_bundle.rds"))
  pp_assert(identical(x$bundle_class, "stage7a_sealed_truth") &&
              identical(x$data_id, data_id) &&
              identical(x$experiment_id,
                        "MULTIFSYNC_COVARIATE_EFFECT_STAGE7A_V1_20260919"),
            paste0("Sealed truth identity mismatch: ", data_id))
  x
}), data_ids)

pp_load_package()
pp_source_runtime(TRUE)
evaluator <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

interp <- function(grid, value, target) {
  stats::approx(grid, as.numeric(value), xout = target,
                rule = 2, ties = mean)$y
}
trapezoid_weights <- function(grid) {
  delta <- diff(grid)
  c(delta[[1L]] / 2, (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}
accumulate <- function(acc, estimate, truth) {
  pp_assert(all(dim(estimate) == dim(truth)) &&
              all(is.finite(estimate)) && all(is.finite(truth)),
            "Signal evaluation dimensions or finite-value checks failed.")
  acc + c(sse = sum((estimate - truth)^2),
          truth_sq = sum(truth^2), n = length(truth))
}
finish_nrmse <- function(acc) sqrt(acc[["sse"]] / acc[["truth_sq"]])
finish_rmse <- function(acc) sqrt(acc[["sse"]] / acc[["n"]])

covariate_metrics <- function(fit, sealed, evaluator_environment) {
  data <- sealed$data
  truth <- data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  p <- truth$p
  pp_assert(length(grid) == 201L && length(fit$list_beta_hat) == p &&
              isTRUE(all.equal(as.numeric(truth$d), 1)) &&
              ncol(truth$beta_true_grid) == p,
            "Expected Stage 7A d=1 beta and 201-point truth grid.")
  weights <- trapezoid_weights(grid)
  beta_hat <- vapply(seq_len(p), function(j) {
    pp_assert(length(fit$list_beta_hat[[j]]) == 1L,
              "Expected one shared covariate coefficient function.")
    interp(fit$time_g, fit$list_beta_hat[[j]][[1L]], grid)
  }, numeric(length(grid)))
  beta_truth <- as.matrix(truth$beta_true_grid)
  beta_ise <- colSums(sweep((beta_hat - beta_truth)^2, 1L, weights, `*`))
  beta_energy <- colSums(sweep(beta_truth^2, 1L, weights, `*`))
  beta_errors <- data.frame(
    variable = seq_len(p), beta_ise = beta_ise,
    truth_beta_energy = beta_energy,
    beta_relative_ise = beta_ise / pmax(beta_energy, 1e-12))
  pp_assert(all(is.finite(beta_ise)) && all(beta_energy > 1e-12),
            "Beta coefficient ISE is non-finite or denominator vanished.")

  paths_obs <- evaluator_environment$v0g_reported_paths_by_study(fit, data)
  names_acc <- c("obs_full_all", "obs_full_selected",
                 "obs_residual_selected", "obs_beta_contribution",
                 "dense_full_all", "dense_full_selected",
                 "dense_residual_selected")
  acc <- setNames(replicate(length(names_acc),
                            c(sse = 0, truth_sq = 0, n = 0),
                            simplify = FALSE), names_acc)
  shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
  specific_mask <- lapply(fit$factor_ppi_specific,
                          function(x) as.numeric(x >= 0.5))
  global_subject_offset <- 0L
  for (s in seq_len(truth$S)) {
    fit_ls <- length(fit$list_Phi_hat_spec[[s]])
    true_ls <- length(truth$kappa_dense_list[[s]])
    pp_assert(fit_ls == 1L && true_ls == 1L,
              "Expected one study-specific factor in each study.")
    beta_obs_offset <- 0L
    estimated_mean_grid <- vapply(seq_len(p), function(j) {
      interp(fit$time_g, fit$list_mu_hat[[s]][[j]], grid)
    }, numeric(length(grid)))
    phase <- truth$phase_mu_mat[s, ]
    true_mean_grid <- vapply(seq_len(p), function(j) {
      0.6 * (sin(2 * pi * grid + phase[[j]]) +
               0.5 * sin(4 * pi * grid + 2 * phase[[j]]))
    }, numeric(length(grid)))
    estimated_shared_functions <- lapply(seq_len(fit$L_f), function(l) {
      evaluator_environment$v0h_interpolate_matrix(
        fit$time_g, fit$list_Phi_hat[[l]], grid)
    })
    estimated_specific_functions <- lapply(seq_len(fit_ls), function(l) {
      evaluator_environment$v0h_interpolate_matrix(
        fit$time_g, fit$list_Phi_hat_spec[[s]][[l]], grid)
    })
    for (i in seq_len(truth$n_s[[s]])) {
      global_i <- global_subject_offset + i
      t_obs <- data$time_obs[[s]][[i]]
      n_obs <- length(t_obs)
      rows <- beta_obs_offset + seq_len(n_obs)
      z <- data$Z[[s]][i, 1L]
      estimated_mean_obs <- vapply(seq_len(p), function(j) {
        interp(fit$time_g, fit$list_mu_hat[[s]][[j]], t_obs)
      }, numeric(n_obs))
      estimated_beta_obs <- vapply(seq_len(p), function(j) {
        interp(fit$time_g, fit$list_beta_hat[[j]][[1L]], t_obs)
      }, numeric(n_obs)) * z
      estimated_shared_obs <- paths_obs$shared[[s]][rows, , drop = FALSE] %*%
        t(fit$mu_q_a_hat)
      estimated_specific_obs <- paths_obs$specific[[s]][rows, , drop = FALSE] %*%
        t(fit$mu_q_b_specific_hat[[s]])
      estimated_residual_obs <- estimated_mean_obs +
        paths_obs$shared[[s]][rows, , drop = FALSE] %*%
          t(sweep(fit$mu_q_a_hat, 2L, shared_mask, `*`)) +
        paths_obs$specific[[s]][rows, , drop = FALSE] %*%
          t(sweep(fit$mu_q_b_specific_hat[[s]], 2L,
                  specific_mask[[s]], `*`))
      true_residual_obs <- truth$signal_true_values[[s]][[i]] -
        truth$beta_true_values[[s]][[i]]
      true_signal_obs <- truth$signal_true_values[[s]][[i]]
      pp_assert(max(abs(true_mean_grid[1L, ] - true_mean_grid[1L, ])) == 0 &&
                  max(abs(vapply(seq_len(p), function(j) {
                    0.6 * (sin(2 * pi * t_obs + phase[[j]]) +
                             0.5 * sin(4 * pi * t_obs + 2 * phase[[j]]))
                  }, numeric(n_obs)) - truth$mu_true_values[[s]][[i]])) < 1e-10,
                "Dense-grid true mean formula differs from observed truth.")
      acc$obs_beta_contribution <- accumulate(
        acc$obs_beta_contribution, estimated_beta_obs,
        truth$beta_true_values[[s]][[i]])
      acc$obs_residual_selected <- accumulate(
        acc$obs_residual_selected, estimated_residual_obs,
        true_residual_obs)
      acc$obs_full_all <- accumulate(
        acc$obs_full_all,
        estimated_mean_obs + estimated_beta_obs +
          estimated_shared_obs + estimated_specific_obs, true_signal_obs)
      acc$obs_full_selected <- accumulate(
        acc$obs_full_selected,
        estimated_residual_obs + estimated_beta_obs, true_signal_obs)

      estimated_shared_dense <- matrix(0, length(grid), fit$L_f)
      true_shared_dense <- matrix(0, length(grid), truth$L_f)
      for (l in seq_len(fit$L_f)) {
        estimated_shared_dense[, l] <- estimated_shared_functions[[l]] %*%
          fit$list_Zeta_hat[[l]][global_i, ]
        true_shared_dense[, l] <- truth$theta_dense_list[[l]] %*%
          truth$eta_true[[s]][[l]][i, ]
      }
      estimated_specific_dense <- matrix(0, length(grid), fit_ls)
      true_specific_dense <- matrix(0, length(grid), true_ls)
      for (l in seq_len(fit_ls)) {
        estimated_specific_dense[, l] <-
          estimated_specific_functions[[l]] %*%
          fit$list_Zeta_hat_spec[[s]][[l]][i, ]
        true_specific_dense[, l] <- truth$kappa_dense_list[[s]][[l]] %*%
          truth$chi_true[[s]][[l]][i, ]
      }
      true_residual_dense <- true_mean_grid +
        true_shared_dense %*% t(truth$a_true) +
        true_specific_dense %*% t(truth$b_true[[s]])
      true_full_dense <- true_residual_dense + z * beta_truth
      estimate_all_dense <- estimated_mean_grid + z * beta_hat +
        estimated_shared_dense %*% t(fit$mu_q_a_hat) +
        estimated_specific_dense %*% t(fit$mu_q_b_specific_hat[[s]])
      estimate_residual_selected_dense <- estimated_mean_grid +
        estimated_shared_dense %*%
          t(sweep(fit$mu_q_a_hat, 2L, shared_mask, `*`)) +
        estimated_specific_dense %*%
          t(sweep(fit$mu_q_b_specific_hat[[s]], 2L,
                  specific_mask[[s]], `*`))
      acc$dense_full_all <- accumulate(
        acc$dense_full_all, estimate_all_dense, true_full_dense)
      acc$dense_full_selected <- accumulate(
        acc$dense_full_selected,
        estimate_residual_selected_dense + z * beta_hat, true_full_dense)
      acc$dense_residual_selected <- accumulate(
        acc$dense_residual_selected,
        estimate_residual_selected_dense, true_residual_dense)
      beta_obs_offset <- beta_obs_offset + n_obs
    }
    global_subject_offset <- global_subject_offset + truth$n_s[[s]]
  }
  summary <- data.frame(
    beta_ise_mean = mean(beta_ise),
    beta_ise_median = median(beta_ise),
    beta_ise_p90 = as.numeric(quantile(beta_ise, 0.9)),
    beta_relative_ise = sum(beta_ise) / sum(beta_energy),
    beta_contribution_obs_rmse = finish_rmse(acc$obs_beta_contribution),
    beta_contribution_obs_nrmse = finish_nrmse(acc$obs_beta_contribution),
    obs_full_signal_nrmse_all = finish_nrmse(acc$obs_full_all),
    obs_full_signal_nrmse_selected = finish_nrmse(acc$obs_full_selected),
    dense_full_signal_nrmse_all = finish_nrmse(acc$dense_full_all),
    dense_full_signal_nrmse_selected = finish_nrmse(acc$dense_full_selected),
    obs_residual_signal_nrmse_selected =
      finish_nrmse(acc$obs_residual_selected),
    dense_residual_signal_nrmse_selected =
      finish_nrmse(acc$dense_residual_selected),
    dense_minus_observed_signal_nrmse =
      finish_nrmse(acc$dense_full_selected) -
        finish_nrmse(acc$obs_full_selected),
    n_observed_values = acc$obs_full_selected[["n"]],
    n_dense_values = acc$dense_full_selected[["n"]])
  list(summary = summary, beta_errors = beta_errors)
}

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
  current_phase <- "load_fit"
  captured <- c2_capture_conditions({
    fit <- readRDS(file.path(output_root, "fits", fit_id, "fit.rds"))
    sealed <- truths[[row$data_id[[1L]]]]
    current_phase <- "historical_evaluator"
    historical <- s6h_evaluate_fit(fit, sealed, evaluator)
    current_phase <- "covariate_evaluator"
    covariate <- covariate_metrics(fit, sealed, evaluator)
    list(historical = historical, covariate = covariate)
  }, phase = "stage7a_post_selection_scientific_evaluation")
  write_csv_new(captured$warnings,
                file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  status <- data.frame(fit_id = fit_id, ok = is.null(captured$error),
                       elapsed_seconds = captured$elapsed_seconds,
                       error = captured$error$message %||% "")
  if (!status$ok) {
    write_lines_new(c("status=ERROR", paste0("phase=", current_phase),
                      paste0("message=", status$error),
                      paste0("call=", captured$error$call %||% ""),
                      paste0("trace=", captured$error$trace_summary %||% "")),
                    file.path(result_dir, "EVALUATION_ERROR.txt"))
    return(list(status = status, compact = NULL))
  }
  result <- captured$value
  save_rds_new(result, file.path(result_dir, "evaluation.rds"))
  write_csv_new(result$covariate$beta_errors,
                file.path(result_dir, "BETA_COEFFICIENT_ERRORS.csv"))
  for (name in names(result$historical)) {
    if (is.data.frame(result$historical[[name]])) {
      write_csv_new(result$historical[[name]],
                    file.path(result_dir, paste0(name, ".csv")))
    }
  }
  terminal <- endpoints[endpoints$fit_id == fit_id, , drop = FALSE]
  legacy <- result$historical$summary
  # The old observed-signal metric omits beta. It is not a valid d>0 metric.
  legacy <- legacy[, !names(legacy) %in% c(
    "signal_nrmse_all_candidates", "signal_nrmse_ppi_selected",
    "signal_nrmse_mean_only"), drop = FALSE]
  names(legacy)[grepl("^dense_signal_", names(legacy))] <- paste0(
    "legacy_residual_", names(legacy)[grepl("^dense_signal_", names(legacy))])
  activity <- result$historical$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  base <- data.frame(
    fit_id = fit_id, data_id = row$data_id[[1L]],
    logical_data_id = row$logical_data_id[[1L]],
    beta_arm = row$beta_arm[[1L]], start_index = row$start_index[[1L]],
    fit_seed = row$fit_seed[[1L]],
    final_ordinary_elbo = terminal$final_ordinary_elbo[[1L]],
    objective_eligible = terminal$objective_eligible[[1L]],
    strict_practical_converged = terminal$strict_practical_converged[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    truth_free_winner = fit_id %in% winners$fit_id,
    truth_used_for_fit_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L)
  compact <- cbind(base, result$covariate$summary, legacy,
                   result$historical$process_contribution_summary,
                   result$historical$loading_scale_summary, activity)
  write_lines_new(c("status=COMPLETE", paste0("fit_id=", fit_id),
                    "truth_used_for_fit_or_selection=FALSE",
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
write_csv_new(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok),
          "Evaluation errors retained; inspect them before new evaluation.")
science <- add_recovery_flags(pp_rbind_fill(lapply(results, `[[`, "compact")))
science <- science[match(manifest$fit_id, science$fit_id), , drop = FALSE]
write_csv_new(science, file.path(evaluation_root, "ALL_START_SCIENCE.csv"))
winner_science <- science[science$truth_free_winner, , drop = FALSE]
write_csv_new(winner_science, file.path(evaluation_root, "WINNER_SCIENCE.csv"))

if (limit == 144L) {
  paired <- merge(
    winner_science[winner_science$beta_arm == "sine", ],
    winner_science[winner_science$beta_arm == "bspline", ],
    by = "logical_data_id", suffixes = c("_sine", "_bspline"),
    sort = TRUE)
  pp_assert(nrow(paired) == 6L, "Expected six paired winner comparisons.")
  write_csv_new(paired, file.path(evaluation_root, "PAIRED_WINNERS.csv"))
}
write_lines_new(c(
  "status=STAGE7A_POST_SELECTION_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("fits_evaluated=", nrow(science)),
  paste0("truth_free_winners_evaluated=", nrow(winner_science)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  "truth_read_after_truth_free_selection_freeze=TRUE",
  "truth_used_for_fit_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE",
  "legacy_observed_signal_metric_excluded_because_it_omits_beta=TRUE",
  "cross_arm_elbo_comparison=FALSE",
  paste0("completed_utc=", pp_iso_time())
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE7A_EVALUATION_COMPLETE fits=", nrow(science),
    " winners=", nrow(winner_science), "\n", sep = "")
