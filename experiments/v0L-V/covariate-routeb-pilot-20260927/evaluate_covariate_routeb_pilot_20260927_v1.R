#!/usr/bin/env Rscript

# Post-selection scientific evaluation for the 36-fit covariate Route-B pilot.
# This script never fits or continues a model and never changes a frozen winner.
# It reads sealed truth only after validating all registered endpoints and the
# truth-free maximum-eligible-ordinary-ELBO selections.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
support_repo <- Sys.getenv("COVARIATE_ROUTEB_EVALUATION_SUPPORT_REPO",
                           unset = "")
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

experiment_id <- "COVARIATE_ROUTEB_PILOT_V1_20260927"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% ""
workers <- as.integer(cli$workers %||% "1")
smoke <- pp_bool(cli$smoke, FALSE)
evaluation_name <- cli$`evaluation-name` %||%
  if (smoke) "scientific_evaluation_smoke_20260927_v1" else
    "scientific_evaluation_20260927_v1"
pp_assert(nzchar(output_root) && dir.exists(output_root) && workers %in% 1:8,
          "A prepared --output-root and 1--8 workers are required.")
pp_assert(pp_bool(cli$`truth-unseal-authorized`, FALSE),
          "Scientific evaluation requires explicit truth-unseal authorization.")
pp_assert(grepl("^[A-Za-z0-9_-]+$", evaluation_name),
          "The evaluation name contains unsupported characters.")
pp_assert(
  identical(tolower(pp_sha256(stage6a_support_path)),
            "bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6") &&
    identical(tolower(pp_sha256(stage6h_support_path)),
              "5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221") &&
    dir.exists(PP_SNAPSHOT),
  "The bound Stage 6A/6H support or frozen evaluation snapshot changed.")

read_csv <- function(path) {
  pp_assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv_new <- function(value, path) {
  pp_assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
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

evaluation_root <- file.path(output_root, evaluation_name)
pp_assert(!dir.exists(evaluation_root),
          "Evaluation output already exists; refusing overwrite.")

data_manifest_path <- file.path(output_root, "DATA_MANIFEST.csv")
fit_manifest_path <- file.path(output_root, "FIT_MANIFEST.csv")
selection_root <- file.path(output_root, "truth_free_selection")
endpoints_path <- file.path(selection_root, "ALL_TRUTH_FREE_ENDPOINTS.csv")
winners_path <- file.path(selection_root, "TRUTH_FREE_WINNERS.csv")
selection_status_path <- file.path(selection_root, "SELECTION_STATUS.csv")
frozen_path <- file.path(selection_root, "TRUTH_FREE_SELECTION_FROZEN.txt")
completion_path <- file.path(
  output_root, if (smoke) "MICRO_SMOKE_COMPLETE.txt" else
    "FITS_AND_SELECTION_COMPLETE.txt")
required <- c(data_manifest_path, fit_manifest_path, endpoints_path,
              winners_path, selection_status_path, frozen_path,
              completion_path)
if (!smoke) required <- c(required,
                          file.path(output_root, "SUPERVISOR_EXIT_STATUS.txt"))
for (path in required) {
  pp_assert(file.exists(path), paste0("Missing required artifact: ", path))
}

data_manifest <- read_csv(data_manifest_path)
manifest <- read_csv(fit_manifest_path)
endpoints <- read_csv(endpoints_path)
winners <- read_csv(winners_path)
selection_status <- read_csv(selection_status_path)
frozen <- readLines(frozen_path, warn = FALSE)
expected_fits <- if (smoke) 3L else 36L
expected_per_arm <- if (smoke) 1L else 12L
expected_data_ids <- paste0("covrbp_01__", c("C0", "C1", "C2"))

pp_assert(
  nrow(data_manifest) == 3L && nrow(manifest) == expected_fits &&
    nrow(endpoints) == expected_fits && nrow(winners) == 3L &&
    nrow(selection_status) == 3L &&
    identical(as.character(data_manifest$data_id), expected_data_ids) &&
    identical(as.character(data_manifest$beta_arm), c("C0", "C1", "C2")) &&
    all(table(manifest$data_id) == expected_per_arm) &&
    !anyDuplicated(manifest$fit_id) && setequal(manifest$fit_id, endpoints$fit_id) &&
    all(endpoints$status == "COMPLETE") && all(endpoints$objective_eligible) &&
    all(is.finite(endpoints$final_ordinary_elbo)) &&
    !any(endpoints$truth_used_for_fit_or_selection) &&
    !any(data_manifest$truth_available_to_fit) &&
    !any(data_manifest$truth_available_to_selection) &&
    !any(data_manifest$truth_unseal_authorized) &&
    !any(manifest$truth_available_to_fit) &&
    !any(manifest$truth_available_to_selection) &&
    all(manifest$n_cpus == 1L) &&
    all(selection_status$n_registered == expected_per_arm) &&
    all(selection_status$n_eligible == expected_per_arm) &&
    all(selection_status$winner_status == "SELECTED"),
  "Manifest, endpoint, or truth-isolation contract failed.")
if (!smoke) {
  pp_assert(all(endpoints$annealing_sweeps == 99L) &&
              all(endpoints$ordinary_t1_sweeps == 400L),
            "The formal endpoints are not the registered 99+400-sweep fits.")
}
pp_assert(all(c(
  "status=TRUTH_FREE_SELECTION_FROZEN",
  "rule=maximum_eligible_ordinary_T1_ELBO_within_data_arm",
  "tie_break=fit_id_lexicographic", "cross_arm_ELBO_comparison=FALSE",
  "truth_used_for_selection=FALSE", "truth_unseal_authorized=FALSE",
  paste0("registered_fits=", expected_fits), "selected_winners=3",
  paste0("smoke=", smoke)
) %in% frozen), "Frozen selection marker is incomplete or changed.")

for (data_id in expected_data_ids) {
  group <- endpoints[endpoints$data_id == data_id &
                       endpoints$objective_eligible &
                       is.finite(endpoints$final_ordinary_elbo), , drop = FALSE]
  chosen <- winners[winners$data_id == data_id, , drop = FALSE]
  expected <- group[order(-group$final_ordinary_elbo, group$fit_id), ,
                    drop = FALSE][1L, , drop = FALSE]
  pp_assert(nrow(group) == expected_per_arm && nrow(chosen) == 1L &&
              identical(chosen$fit_id[[1L]], expected$fit_id[[1L]]),
            paste0("Winner is not maximum eligible ELBO for ", data_id, "."))
}

fit_complete <- file.path(output_root, "fits", manifest$fit_id,
                          "FIT_COMPLETE.txt")
fit_error <- file.path(output_root, "fits", manifest$fit_id,
                       "FIT_ERROR.txt")
pp_assert(all(file.exists(fit_complete)) && !any(file.exists(fit_error)),
          "Fit completion accounting failed.")

obs_c1 <- readRDS(file.path(output_root, "data", "covrbp_01__C1",
                            "observation_bundle.rds"))
obs_c2 <- readRDS(file.path(output_root, "data", "covrbp_01__C2",
                            "observation_bundle.rds"))
pp_assert(identical(obs_c1$Y, obs_c2$Y) &&
            identical(obs_c1$time_obs, obs_c2$time_obs) &&
            is.list(obs_c1$Z) && is.null(obs_c2$Z) &&
            identical(as.integer(obs_c1$dimensions$d), 2L) &&
            identical(as.integer(obs_c2$dimensions$d), 0L),
          "C1/C2 paired observation or omitted-covariate contract failed.")

if (!smoke) {
  binding_path <- file.path(script_root, "EXECUTION_SOURCE_BINDING_V1.csv")
  binding <- read_csv(binding_path)
  pp_assert(nrow(binding) == 11L &&
              all(c("relative_path_or_artifact", "sha256", "role") %in%
                    names(binding)),
            "Execution source binding is malformed.")
  source_rows <- binding[binding$role !=
                           "uploaded_source_bundle_before_launcher_dependency_path_update",
                         , drop = FALSE]
  source_matches <- vapply(seq_len(nrow(source_rows)), function(index) {
    path <- file.path(PP_REPO_ROOT,
                      source_rows$relative_path_or_artifact[[index]])
    file.exists(path) && identical(tolower(pp_sha256(path)),
                                   tolower(source_rows$sha256[[index]]))
  }, logical(1L))
  pp_assert(nrow(source_rows) == 10L && all(source_matches), paste0(
    "Bound fitting source mismatch: ", paste(
      source_rows$relative_path_or_artifact[!source_matches], collapse = ";")))
}

dir.create(evaluation_root, recursive = TRUE, mode = "0700")
write_lines_new(c(
  "status=POST_SELECTION_TRUTH_UNSEAL_AUTHORIZED_BY_USER",
  "authorization_date=2026-09-27",
  paste0("frozen_endpoints_sha256=", pp_sha256(endpoints_path)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  paste0("evaluator_sha256=", pp_sha256(script_file)),
  paste0("stage6a_support_sha256=", pp_sha256(stage6a_support_path)),
  paste0("stage6h_support_sha256=", pp_sha256(stage6h_support_path)),
  "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_arm",
  "cross_arm_ELBO_comparison=FALSE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE",
  paste0("smoke=", smoke)
), file.path(evaluation_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))

# First sealed-truth read. Every endpoint, selection and source check above has
# already completed.
truths <- setNames(lapply(expected_data_ids, function(data_id) {
  truth <- readRDS(file.path(output_root, "data", data_id,
                             "sealed_truth", "truth_bundle.rds"))
  arm <- sub("^.*__", "", data_id)
  pp_assert(identical(truth$bundle_class,
                      "covariate_routeb_pilot_sealed_truth") &&
              identical(truth$experiment_id, experiment_id) &&
              identical(truth$data_id, data_id) &&
              isFALSE(truth$truth_unseal_authorized) &&
              identical(truth$data$true_params$truth_basis_method,
                        "iid_common_splinet_stiefel") &&
              identical(truth$data$true_params$mean_truth_method,
                        "independent_splinet") &&
              identical(truth$data$true_params$beta_truth_method,
                        if (arm == "C0") "zero_effect_d2" else
                          "independent_cubic_bspline_d2_matched_mean_centered_energy"),
            paste0("Sealed truth identity mismatch: ", data_id))
  truth
}), expected_data_ids)
pp_assert(max(abs(truths[["covrbp_01__C0"]]$data$true_params$beta_true_grid)) == 0 &&
            max(abs(truths[["covrbp_01__C1"]]$data$true_params$beta_true_grid)) > 0 &&
            identical(truths[["covrbp_01__C1"]]$data,
                      truths[["covrbp_01__C2"]]$data),
          "C0/C1/C2 sealed-truth pairing failed after authorized unsealing.")

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)

trapezoid_weights <- function(grid) {
  delta <- diff(as.numeric(grid))
  c(delta[[1L]] / 2,
    (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}
interp <- function(grid, value, target) {
  stats::approx(grid, as.numeric(value), xout = target,
                rule = 2, ties = mean)$y
}
mean_design <- function(grid, truth) {
  cbind(1, as.numeric(grid), multiFSYNC:::ZOSull(
    as.numeric(grid), range.x = c(0, 1), intKnots = truth$int_knots))
}
true_mean_dense <- function(data) {
  pp_assert(identical(data$true_params$mean_truth_method,
                      "independent_splinet"),
            "This evaluator requires the registered independent Splinets mean.")
  data$true_params$mean_dense_list
}
estimated_mean_dense <- function(fit, data, study) {
  coefficient <- do.call(cbind, fit$mu_q_nu_mu_original[[study]])
  design <- mean_design(data$true_params$t_grid_dense, data$true_params)
  pp_assert(nrow(coefficient) == ncol(design),
            "Mean coefficient and fitted-basis dimensions differ.")
  design %*% coefficient
}
mean_recovery <- function(fit, data) {
  truth <- data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  weights <- trapezoid_weights(grid)
  design <- mean_design(grid, truth)
  weighted_design <- design * sqrt(weights)
  target_list <- true_mean_dense(data)
  rows <- lapply(seq_len(truth$S), function(study) {
    estimate <- estimated_mean_dense(fit, data, study)
    target <- target_list[[study]]
    projection_coefficient <- qr.solve(weighted_design,
                                       target * sqrt(weights))
    projection <- design %*% projection_coefficient
    ise <- colSums(weights * (estimate - target)^2)
    energy <- colSums(weights * target^2)
    floor <- colSums(weights * (projection - target)^2)
    estimate_projection <- colSums(weights * (estimate - projection)^2)
    observed_sse <- observed_energy <- observed_n <- 0
    coefficient <- do.call(cbind, fit$mu_q_nu_mu_original[[study]])
    for (subject in seq_len(truth$n_s[[study]])) {
      time <- data$time_obs[[study]][[subject]]
      observed_estimate <- mean_design(time, truth) %*% coefficient
      observed_truth <- truth$mu_true_values[[study]][[subject]]
      observed_sse <- observed_sse + sum((observed_estimate - observed_truth)^2)
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
      projection_relative_ise_mean = mean(floor / pmax(energy, 1e-12)),
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
                                 max(sum(by_study$observed_truth_energy), 1e-12)),
    mean_dense_rmse = sqrt(sum(by_study$dense_ise_sum) /
                             (total_variables * sum(weights))),
    mean_dense_nrmse = sqrt(sum(by_study$dense_ise_sum) /
                              max(sum(by_study$dense_truth_energy), 1e-12)),
    mean_dense_ise_mean = sum(by_study$dense_ise_sum) / total_variables,
    mean_dense_relative_ise_mean = mean(by_study$dense_relative_ise_mean),
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

# Replace only the historical analytic-mean assembly. The frozen identity,
# loading, score, process and feature matching logic is unchanged.
install_dense_residual_override <- function(environment) {
  environment$true_mean_dense <- true_mean_dense
  environment$estimated_mean_dense <- estimated_mean_dense
  dense_function <- function(fit, data) {
    truth <- data$true_params
    grid <- as.numeric(truth$t_grid_dense)
    true_L_s <- if (!is.null(truth$L_s_by_study)) {
      as.integer(unlist(truth$L_s_by_study, use.names = FALSE))
    } else if (length(truth$L_s) == truth$S) {
      as.integer(unlist(truth$L_s, use.names = FALSE))
    } else rep.int(as.integer(truth$L_s), truth$S)
    true_means <- true_mean_dense(data)
    shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
    specific_mask <- lapply(fit$factor_ppi_specific,
                            function(value) as.numeric(value >= 0.5))
    accumulators <- list(
      all = c(sse = 0, truth = 0, n = 0),
      selected = c(sse = 0, truth = 0, n = 0),
      mean_only = c(sse = 0, truth = 0, n = 0))
    subject_offset <- 0L
    for (study in seq_len(truth$S)) {
      estimated_mean <- estimated_mean_dense(fit, data, study)
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
          true_shared_path[, factor] <- truth$theta_dense_list[[factor]] %*%
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
        true_residual <- true_mean +
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
          error <- estimates[[name]] - true_residual
          accumulators[[name]][["sse"]] <-
            accumulators[[name]][["sse"]] + sum(error^2)
          accumulators[[name]][["truth"]] <-
            accumulators[[name]][["truth"]] + sum(true_residual^2)
          accumulators[[name]][["n"]] <-
            accumulators[[name]][["n"]] + length(true_residual)
        }
      }
      subject_offset <- subject_offset + truth$n_s[[study]]
    }
    do.call(rbind, lapply(names(accumulators), function(name) {
      value <- accumulators[[name]]
      data.frame(
        scope = name, grid_size = length(grid),
        squared_error = value[["sse"]],
        truth_squared_energy = value[["truth"]], scalar_count = value[["n"]],
        rmse = sqrt(value[["sse"]] / value[["n"]]),
        truth_rms = sqrt(value[["truth"]] / value[["n"]]),
        nrmse = sqrt(value[["sse"]] / max(value[["truth"]], 1e-12)),
        relative_mise = value[["sse"]] / max(value[["truth"]], 1e-12),
        stringsAsFactors = FALSE)
    }))
  }
  environment(dense_function) <- environment
  environment$v0j_dense_signal_metrics <- dense_function
  invisible(environment)
}
install_dense_residual_override(evaluator_environment)

new_error_accumulator <- function() c(sse = 0, truth_sq = 0, n = 0)
accumulate_error <- function(accumulator, estimate, truth) {
  pp_assert(identical(dim(estimate), dim(truth)) &&
              all(is.finite(estimate)) && all(is.finite(truth)),
            "Signal dimensions or finite-value checks failed.")
  accumulator + c(sse = sum((estimate - truth)^2),
                  truth_sq = sum(truth^2), n = length(truth))
}
finish_rmse <- function(accumulator) {
  sqrt(accumulator[["sse"]] / accumulator[["n"]])
}
finish_nrmse <- function(accumulator) {
  if (accumulator[["truth_sq"]] <= 1e-12) NA_real_ else
    sqrt(accumulator[["sse"]] / accumulator[["truth_sq"]])
}
new_projection_accumulator <- function() c(cross = 0, dev_sq = 0, beta_sq = 0)
accumulate_projection <- function(accumulator, deviation, beta) {
  accumulator + c(cross = sum(deviation * beta),
                  dev_sq = sum(deviation^2), beta_sq = sum(beta^2))
}
finish_projection <- function(accumulator) {
  c(
    coefficient = if (accumulator[["beta_sq"]] <= 1e-12) NA_real_ else
      accumulator[["cross"]] / accumulator[["beta_sq"]],
    cosine = if (accumulator[["beta_sq"]] <= 1e-12 ||
                   accumulator[["dev_sq"]] <= 1e-12) NA_real_ else
      accumulator[["cross"]] /
        sqrt(accumulator[["beta_sq"]] * accumulator[["dev_sq"]]))
}

covariate_and_signal_metrics <- function(fit, sealed, evaluator_environment,
                                         beta_arm) {
  data <- sealed$data
  truth <- data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  p <- truth$p
  d_truth <- 2L
  weights <- trapezoid_weights(grid)
  beta_truth <- truth$beta_true_grid
  pp_assert(length(dim(beta_truth)) == 3L &&
              identical(as.integer(dim(beta_truth)),
                        c(length(grid), p, d_truth)),
            "Expected the registered G by p by 2 beta truth array.")
  beta_available <- !identical(beta_arm, "C2")
  beta_hat <- array(0, dim = dim(beta_truth))
  if (beta_available) {
    pp_assert(length(fit$list_beta_hat) == p &&
                all(vapply(fit$list_beta_hat, length, integer(1L)) == d_truth),
              "Expected two fitted shared covariate coefficient functions.")
    for (variable in seq_len(p)) for (covariate in seq_len(d_truth)) {
      beta_hat[, variable, covariate] <- interp(
        fit$time_g, fit$list_beta_hat[[variable]][[covariate]], grid)
    }
  }
  coefficient_rows <- list()
  for (covariate in seq_len(d_truth)) for (variable in seq_len(p)) {
    target <- beta_truth[, variable, covariate]
    estimate <- beta_hat[, variable, covariate]
    truth_energy <- sum(weights * target^2)
    estimated_energy <- sum(weights * estimate^2)
    zero_imputed_ise <- sum(weights * (estimate - target)^2)
    inner <- sum(weights * estimate * target)
    coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
      variable = variable, covariate = covariate,
      beta_estimate_available = beta_available,
      beta_ise = if (beta_available) zero_imputed_ise else NA_real_,
      beta_zero_imputed_ise = zero_imputed_ise,
      truth_beta_energy = truth_energy,
      estimated_beta_energy = if (beta_available) estimated_energy else NA_real_,
      beta_l2_cosine = if (!beta_available || truth_energy <= 1e-12 ||
                              estimated_energy <= 1e-12) NA_real_ else
        inner / sqrt(truth_energy * estimated_energy),
      stringsAsFactors = FALSE)
  }
  beta_errors <- do.call(rbind, coefficient_rows)
  true_beta_total <- sum(beta_errors$truth_beta_energy)
  actual_beta_ise <- if (beta_available) sum(beta_errors$beta_ise) else NA_real_

  error_names <- c(
    "obs_full_all", "obs_full_selected", "obs_residual_selected",
    "obs_factor_selected", "obs_mean", "obs_beta_contribution",
    "dense_full_all", "dense_full_selected", "dense_residual_selected",
    "dense_factor_selected", "dense_mean", "dense_beta_contribution")
  acc <- setNames(replicate(length(error_names), new_error_accumulator(),
                            simplify = FALSE), error_names)
  projection_names <- c("obs_factor", "obs_residual",
                        "dense_factor", "dense_residual")
  projection <- setNames(replicate(length(projection_names),
                                   new_projection_accumulator(),
                                   simplify = FALSE), projection_names)
  paths_obs <- evaluator_environment$v0g_reported_paths_by_study(fit, data)
  shared_mask <- as.numeric(fit$factor_ppi_shared >= 0.5)
  specific_mask <- lapply(fit$factor_ppi_specific,
                          function(value) as.numeric(value >= 0.5))
  true_means <- true_mean_dense(data)
  global_subject_offset <- 0L
  for (study in seq_len(truth$S)) {
    fit_ls <- length(fit$list_Phi_hat_spec[[study]])
    true_ls <- length(truth$kappa_dense_list[[study]])
    pp_assert(fit_ls == 1L && true_ls == 1L,
              "Expected one study-specific factor in each study.")
    beta_obs_offset <- 0L
    estimated_mean_grid <- estimated_mean_dense(fit, data, study)
    true_mean_grid <- true_means[[study]]
    estimated_shared_functions <- lapply(seq_len(fit$L_f), function(factor) {
      evaluator_environment$v0h_interpolate_matrix(
        fit$time_g, fit$list_Phi_hat[[factor]], grid)
    })
    estimated_specific_functions <- lapply(seq_len(fit_ls), function(factor) {
      evaluator_environment$v0h_interpolate_matrix(
        fit$time_g, fit$list_Phi_hat_spec[[study]][[factor]], grid)
    })
    for (subject in seq_len(truth$n_s[[study]])) {
      global_subject <- global_subject_offset + subject
      time_obs <- data$time_obs[[study]][[subject]]
      n_obs <- length(time_obs)
      rows <- beta_obs_offset + seq_len(n_obs)
      z <- as.numeric(data$Z[[study]][subject, seq_len(d_truth)])
      estimated_mean_obs <- vapply(seq_len(p), function(variable) {
        interp(fit$time_g, fit$list_mu_hat[[study]][[variable]], time_obs)
      }, numeric(n_obs))
      estimated_beta_obs <- matrix(0, n_obs, p)
      if (beta_available) {
        for (covariate in seq_len(d_truth)) {
          estimated_beta_obs <- estimated_beta_obs + z[[covariate]] *
            vapply(seq_len(p), function(variable) {
              interp(fit$time_g,
                     fit$list_beta_hat[[variable]][[covariate]], time_obs)
            }, numeric(n_obs))
        }
      }
      estimated_shared_obs <- paths_obs$shared[[study]][rows, , drop = FALSE] %*%
        t(fit$mu_q_a_hat)
      estimated_specific_obs <-
        paths_obs$specific[[study]][rows, , drop = FALSE] %*%
        t(fit$mu_q_b_specific_hat[[study]])
      estimated_factor_all_obs <- estimated_shared_obs + estimated_specific_obs
      estimated_factor_selected_obs <-
        paths_obs$shared[[study]][rows, , drop = FALSE] %*%
          t(sweep(fit$mu_q_a_hat, 2L, shared_mask, "*")) +
        paths_obs$specific[[study]][rows, , drop = FALSE] %*%
          t(sweep(fit$mu_q_b_specific_hat[[study]], 2L,
                  specific_mask[[study]], "*"))
      true_mean_obs <- truth$mu_true_values[[study]][[subject]]
      true_beta_obs <- truth$beta_true_values[[study]][[subject]]
      true_full_obs <- truth$signal_true_values[[study]][[subject]]
      true_residual_obs <- true_full_obs - true_beta_obs
      true_factor_obs <- true_residual_obs - true_mean_obs
      estimated_residual_selected_obs <-
        estimated_mean_obs + estimated_factor_selected_obs
      acc$obs_mean <- accumulate_error(acc$obs_mean,
                                       estimated_mean_obs, true_mean_obs)
      acc$obs_beta_contribution <- accumulate_error(
        acc$obs_beta_contribution, estimated_beta_obs, true_beta_obs)
      acc$obs_factor_selected <- accumulate_error(
        acc$obs_factor_selected, estimated_factor_selected_obs, true_factor_obs)
      acc$obs_residual_selected <- accumulate_error(
        acc$obs_residual_selected, estimated_residual_selected_obs,
        true_residual_obs)
      acc$obs_full_all <- accumulate_error(
        acc$obs_full_all,
        estimated_mean_obs + estimated_beta_obs + estimated_factor_all_obs,
        true_full_obs)
      acc$obs_full_selected <- accumulate_error(
        acc$obs_full_selected,
        estimated_residual_selected_obs + estimated_beta_obs, true_full_obs)
      projection$obs_factor <- accumulate_projection(
        projection$obs_factor,
        estimated_factor_selected_obs - true_factor_obs, true_beta_obs)
      projection$obs_residual <- accumulate_projection(
        projection$obs_residual,
        estimated_residual_selected_obs - true_residual_obs, true_beta_obs)

      estimated_shared_dense <- matrix(0, length(grid), fit$L_f)
      true_shared_dense <- matrix(0, length(grid), truth$L_f)
      for (factor in seq_len(fit$L_f)) {
        estimated_shared_dense[, factor] <-
          estimated_shared_functions[[factor]] %*%
          fit$list_Zeta_hat[[factor]][global_subject, ]
        true_shared_dense[, factor] <- truth$theta_dense_list[[factor]] %*%
          truth$eta_true[[study]][[factor]][subject, ]
      }
      estimated_specific_dense <- matrix(0, length(grid), fit_ls)
      true_specific_dense <- matrix(0, length(grid), true_ls)
      for (factor in seq_len(fit_ls)) {
        estimated_specific_dense[, factor] <-
          estimated_specific_functions[[factor]] %*%
          fit$list_Zeta_hat_spec[[study]][[factor]][subject, ]
        true_specific_dense[, factor] <-
          truth$kappa_dense_list[[study]][[factor]] %*%
          truth$chi_true[[study]][[factor]][subject, ]
      }
      estimated_factor_all_dense <-
        estimated_shared_dense %*% t(fit$mu_q_a_hat) +
        estimated_specific_dense %*% t(fit$mu_q_b_specific_hat[[study]])
      estimated_factor_selected_dense <-
        estimated_shared_dense %*%
          t(sweep(fit$mu_q_a_hat, 2L, shared_mask, "*")) +
        estimated_specific_dense %*%
          t(sweep(fit$mu_q_b_specific_hat[[study]], 2L,
                  specific_mask[[study]], "*"))
      true_factor_dense <- true_shared_dense %*% t(truth$a_true) +
        true_specific_dense %*% t(truth$b_true[[study]])
      true_residual_dense <- true_mean_grid + true_factor_dense
      true_beta_dense <- matrix(0, length(grid), p)
      estimated_beta_dense <- matrix(0, length(grid), p)
      for (covariate in seq_len(d_truth)) {
        true_beta_dense <- true_beta_dense +
          z[[covariate]] * beta_truth[, , covariate]
        estimated_beta_dense <- estimated_beta_dense +
          z[[covariate]] * beta_hat[, , covariate]
      }
      true_full_dense <- true_residual_dense + true_beta_dense
      estimated_residual_selected_dense <-
        estimated_mean_grid + estimated_factor_selected_dense
      acc$dense_mean <- accumulate_error(acc$dense_mean,
                                         estimated_mean_grid, true_mean_grid)
      acc$dense_beta_contribution <- accumulate_error(
        acc$dense_beta_contribution, estimated_beta_dense, true_beta_dense)
      acc$dense_factor_selected <- accumulate_error(
        acc$dense_factor_selected, estimated_factor_selected_dense,
        true_factor_dense)
      acc$dense_residual_selected <- accumulate_error(
        acc$dense_residual_selected, estimated_residual_selected_dense,
        true_residual_dense)
      acc$dense_full_all <- accumulate_error(
        acc$dense_full_all,
        estimated_mean_grid + estimated_beta_dense + estimated_factor_all_dense,
        true_full_dense)
      acc$dense_full_selected <- accumulate_error(
        acc$dense_full_selected,
        estimated_residual_selected_dense + estimated_beta_dense,
        true_full_dense)
      projection$dense_factor <- accumulate_projection(
        projection$dense_factor,
        estimated_factor_selected_dense - true_factor_dense, true_beta_dense)
      projection$dense_residual <- accumulate_projection(
        projection$dense_residual,
        estimated_residual_selected_dense - true_residual_dense,
        true_beta_dense)
      beta_obs_offset <- beta_obs_offset + n_obs
    }
    global_subject_offset <- global_subject_offset + truth$n_s[[study]]
  }
  projection_result <- lapply(projection, finish_projection)
  finite_cosines <- beta_errors$beta_l2_cosine[
    is.finite(beta_errors$beta_l2_cosine)]
  summary <- data.frame(
    beta_estimate_available = beta_available,
    beta_true_energy_total = true_beta_total,
    beta_estimated_energy_total = if (beta_available)
      sum(beta_errors$estimated_beta_energy) else NA_real_,
    beta_absolute_ise_total = actual_beta_ise,
    beta_absolute_ise_mean = if (beta_available)
      mean(beta_errors$beta_ise) else NA_real_,
    beta_relative_ise_when_identifiable = if (beta_available &&
                                                 true_beta_total > 1e-12)
      actual_beta_ise / true_beta_total else NA_real_,
    beta_zero_imputed_relative_ise = if (true_beta_total > 1e-12)
      sum(beta_errors$beta_zero_imputed_ise) / true_beta_total else NA_real_,
    beta_function_l2_cosine_mean = if (length(finite_cosines))
      mean(finite_cosines) else NA_real_,
    beta_contribution_obs_rmse = finish_rmse(acc$obs_beta_contribution),
    beta_contribution_obs_nrmse = finish_nrmse(acc$obs_beta_contribution),
    beta_contribution_dense_rmse = finish_rmse(acc$dense_beta_contribution),
    beta_contribution_dense_nrmse = finish_nrmse(acc$dense_beta_contribution),
    obs_full_signal_nrmse_all = finish_nrmse(acc$obs_full_all),
    obs_full_signal_nrmse_selected = finish_nrmse(acc$obs_full_selected),
    dense_full_signal_nrmse_all = finish_nrmse(acc$dense_full_all),
    dense_full_signal_nrmse_selected = finish_nrmse(acc$dense_full_selected),
    obs_residual_signal_nrmse_selected =
      finish_nrmse(acc$obs_residual_selected),
    dense_residual_signal_nrmse_selected =
      finish_nrmse(acc$dense_residual_selected),
    obs_factor_contribution_nrmse_selected =
      finish_nrmse(acc$obs_factor_selected),
    dense_factor_contribution_nrmse_selected =
      finish_nrmse(acc$dense_factor_selected),
    obs_mean_nrmse = finish_nrmse(acc$obs_mean),
    dense_mean_nrmse = finish_nrmse(acc$dense_mean),
    dense_minus_observed_signal_nrmse =
      finish_nrmse(acc$dense_full_selected) -
      finish_nrmse(acc$obs_full_selected),
    factor_beta_absorption_projection_obs =
      projection_result$obs_factor[["coefficient"]],
    factor_beta_absorption_cosine_obs =
      projection_result$obs_factor[["cosine"]],
    residual_beta_absorption_projection_obs =
      projection_result$obs_residual[["coefficient"]],
    residual_beta_absorption_cosine_obs =
      projection_result$obs_residual[["cosine"]],
    factor_beta_absorption_projection_dense =
      projection_result$dense_factor[["coefficient"]],
    factor_beta_absorption_cosine_dense =
      projection_result$dense_factor[["cosine"]],
    residual_beta_absorption_projection_dense =
      projection_result$dense_residual[["coefficient"]],
    residual_beta_absorption_cosine_dense =
      projection_result$dense_residual[["cosine"]],
    n_observed_values = acc$obs_full_selected[["n"]],
    n_dense_values = acc$dense_full_selected[["n"]],
    stringsAsFactors = FALSE)
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
    current_phase <- "frozen_factor_evaluator"
    historical <- s6h_evaluate_fit(fit, sealed, evaluator_environment)
    current_phase <- "mean_evaluator"
    mean_result <- mean_recovery(fit, sealed$data)
    current_phase <- "covariate_and_full_signal_evaluator"
    covariate <- covariate_and_signal_metrics(
      fit, sealed, evaluator_environment, row$beta_arm[[1L]])
    list(historical = historical, mean = mean_result,
         covariate = covariate)
  }, phase = "covariate_routeb_post_selection_scientific_evaluation")
  write_csv_new(captured$warnings,
                file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  status <- data.frame(
    fit_id = fit_id, ok = is.null(captured$error),
    elapsed_seconds = captured$elapsed_seconds,
    error = captured$error$message %||% "", stringsAsFactors = FALSE)
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
  write_csv_new(result$mean$by_study,
                file.path(result_dir, "mean_recovery_by_study.csv"))
  write_csv_new(result$mean$summary,
                file.path(result_dir, "mean_recovery_summary.csv"))
  for (name in names(result$historical)) {
    if (is.data.frame(result$historical[[name]])) {
      write_csv_new(result$historical[[name]],
                    file.path(result_dir, paste0(name, ".csv")))
    }
  }
  terminal <- endpoints[endpoints$fit_id == fit_id, , drop = FALSE]
  legacy <- result$historical$summary
  legacy <- legacy[, !names(legacy) %in% c(
    "signal_nrmse_all_candidates", "signal_nrmse_ppi_selected",
    "signal_nrmse_mean_only"), drop = FALSE]
  names(legacy)[grepl("^dense_signal_", names(legacy))] <- paste0(
    "legacy_residual_",
    names(legacy)[grepl("^dense_signal_", names(legacy))])
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
    truth_used_for_fit_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE)
  compact <- cbind(base, result$covariate$summary, legacy,
                   result$historical$process_contribution_summary,
                   result$historical$loading_scale_summary,
                   result$mean$summary, activity)
  write_lines_new(c("status=COMPLETE", paste0("fit_id=", fit_id),
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
write_csv_new(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok),
          "Evaluation errors retained; no selective evaluation rerun.")
science <- add_recovery_flags(pp_rbind_fill(lapply(results, `[[`, "compact")))
science <- science[match(manifest$fit_id, science$fit_id), , drop = FALSE]
write_csv_new(science, file.path(evaluation_root, "ALL_START_SCIENCE.csv"))
winner_science <- science[science$truth_free_winner, , drop = FALSE]
write_csv_new(winner_science, file.path(evaluation_root, "WINNER_SCIENCE.csv"))

arm_summary <- do.call(rbind, lapply(split(science, science$beta_arm), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  data.frame(
    beta_arm = x$beta_arm[[1L]], starts = nrow(x),
    core_hits = sum(x$core_truth_recovered),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    full_feature_coverage_hits = sum(x$full_true_feature_component_coverage),
    strict_convergence_hits = sum(x$strict_practical_converged),
    winner_core_hit = winner$core_truth_recovered[[1L]],
    winner_exact_structure_hit = winner$exact_thresholded_factor_structure[[1L]],
    winner_full_feature_coverage_hit =
      winner$full_true_feature_component_coverage[[1L]],
    winner_strict_convergence = winner$strict_practical_converged[[1L]],
    winner_obs_full_signal_nrmse = winner$obs_full_signal_nrmse_selected[[1L]],
    winner_dense_full_signal_nrmse = winner$dense_full_signal_nrmse_selected[[1L]],
    winner_beta_relative_ise =
      winner$beta_relative_ise_when_identifiable[[1L]],
    winner_mean_dense_nrmse = winner$mean_dense_nrmse[[1L]],
    broad_generalization_claim = FALSE, stringsAsFactors = FALSE)
}))
write_csv_new(arm_summary, file.path(evaluation_root, "ARM_SUMMARY.csv"))

paired_c1_c2 <- merge(
  science[science$beta_arm == "C1", ],
  science[science$beta_arm == "C2", ],
  by = c("logical_data_id", "start_index", "fit_seed"),
  suffixes = c("_C1", "_C2"), sort = TRUE)
pp_assert(nrow(paired_c1_c2) == expected_per_arm,
          "C1/C2 paired-start scientific join failed.")
write_csv_new(paired_c1_c2,
              file.path(evaluation_root, "PAIRED_START_C1_C2.csv"))

write_lines_new(c(
  "status=COVARIATE_ROUTEB_SCIENTIFIC_EVALUATION_COMPLETE",
  paste0("fits_evaluated=", nrow(science)),
  paste0("truth_free_winners_evaluated=", nrow(winner_science)),
  paste0("frozen_winners_sha256=", pp_sha256(winners_path)),
  "truth_read_after_truth_free_selection_freeze=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE",
  "legacy_observed_signal_metric_excluded_because_it_omits_beta=TRUE",
  "cross_arm_elbo_comparison=FALSE",
  "broad_generalization_claim=FALSE",
  paste0("smoke=", smoke), paste0("completed_utc=", pp_iso_time())
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("COVARIATE_ROUTEB_EVALUATION_COMPLETE fits=", nrow(science),
    " winners=", nrow(winner_science), " smoke=", smoke, "\n", sep = "")
