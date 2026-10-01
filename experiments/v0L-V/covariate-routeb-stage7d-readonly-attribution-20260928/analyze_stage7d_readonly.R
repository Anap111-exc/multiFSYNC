#!/usr/bin/env Rscript
# Stage 7D: post-selection, read-only attribution of Stage 7C structure failures.
# Reads frozen fits and already unsealed truth. Never calls a fitting routine.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analyze_stage7d_readonly.R STAGE7C_ROOT NEW_OUTPUT_ROOT")
}
input_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_root <- args[[2L]]
if (file.exists(output_root)) stop("Output path already exists: ", output_root)

suppressPackageStartupMessages(library(multiFSYNC))
stopifnot(requireNamespace("Splinets", quietly = TRUE))

read_csv <- function(...) utils::read.csv(file.path(...), stringsAsFactors = FALSE)
manifest <- read_csv(input_root, "FIT_MANIFEST.csv")
science <- read_csv(input_root, "scientific_evaluation_20260928_v1",
                    "ALL_START_SCIENCE.csv")
winners <- read_csv(input_root, "truth_free_selection", "TRUTH_FREE_WINNERS.csv")
stopifnot(nrow(manifest) == 108L, nrow(science) == 108L,
          nrow(winners) == 9L, !anyDuplicated(manifest$fit_id),
          identical(sort(manifest$fit_id), sort(science$fit_id)),
          all(science$objective_eligible),
          all(!science$truth_used_for_fit_stopping_or_selection),
          file.exists(file.path(input_root, "FITS_AND_SELECTION_COMPLETE.txt")),
          file.exists(file.path(input_root, "truth_free_selection",
                                "TRUTH_FREE_SELECTION_FROZEN.txt")),
          file.exists(file.path(input_root, "scientific_evaluation_20260928_v1",
                                "EVALUATION_COMPLETE.txt")))

trap_weights <- function(t) {
  dt <- diff(t)
  c(dt[[1L]] / 2, (dt[-1L] + dt[-length(dt)]) / 2,
    dt[[length(dt)]] / 2)
}
cosine <- function(x, y) {
  nx <- sqrt(sum(x * x)); ny <- sqrt(sum(y * y))
  if (nx <= 0 || ny <= 0) return(NA_real_)
  sum(x * y) / (nx * ny)
}
eval_truth_features <- function(meta, times) {
  system <- Splinets::splinet(meta$knots, degree = meta$degree,
                             periodic = meta$periodic, norm = TRUE)
  splines <- Splinets::lincomb(system$os,
                              t(meta$stiefel_coefficients), reduced = TRUE)
  as.matrix(Splinets::evspline(splines, x = times)[, -1L, drop = FALSE])
}
subspace_correlation <- function(A, B, weights) {
  QA <- qr.Q(qr(sqrt(weights) * A))[, seq_len(ncol(A)), drop = FALSE]
  QB <- qr.Q(qr(sqrt(weights) * B))[, seq_len(ncol(B)), drop = FALSE]
  min(1, max(svd(crossprod(QA, QB), nu = 0L, nv = 0L)$d))
}
project_floor <- function(features, design, weights) {
  qr_design <- qr(sqrt(weights) * design)
  stopifnot(qr_design$rank == ncol(design))
  coefficients <- qr.coef(qr_design, sqrt(weights) * features)
  error <- features - design %*% coefficients
  sum(weights * error^2) / sum(weights * features^2)
}
key_row <- function(df, block, study, iteration, phase = "main") {
  z <- df[df$block == block & !is.na(df$study) &
            df$study == study & df$iteration == iteration &
            df$phase == phase, , drop = FALSE]
  z
}
get_scalar <- function(df, block, study, iteration, field,
                       phase = "main") {
  z <- key_row(df, block, study, iteration, phase)
  if (nrow(z) != 1L) return(NA_real_)
  as.numeric(z[[field]][[1L]])
}
get_time_norm <- function(df, study, iteration) {
  z <- key_row(df, "specific_time_function", study, iteration)
  if (nrow(z) != 2L) return(NA_real_)
  sqrt(sum(z$posterior_mean_norm^2))
}

data_ids <- paste0("covrb7c_0", 1:3)
truth_rows <- list(); design_rows <- list(); geometry_rows <- list()
loading_rows <- list(); contribution_rows <- list(); binding_rows <- list()
for (data_id in data_ids) {
  c1_id <- paste0(data_id, "__C1")
  c0_id <- paste0(data_id, "__C0")
  c1 <- readRDS(file.path(input_root, "data", c1_id,
                          "sealed_truth", "truth_bundle.rds"))
  c0 <- readRDS(file.path(input_root, "data", c0_id,
                          "sealed_truth", "truth_bundle.rds"))
  stopifnot(identical(c1$data_id, c1_id), identical(c0$data_id, c0_id),
            identical(c1$data$true_params$a_true,
                      c0$data$true_params$a_true),
            identical(c1$data$true_params$b_true,
                      c0$data$true_params$b_true),
            identical(c1$data$true_params$phi_dense_list,
                      c0$data$true_params$phi_dense_list),
            identical(c1$data$true_params$psi_dense_list,
                      c0$data$true_params$psi_dense_list))
  data <- c1$data
  truth <- data$true_params
  stopifnot(truth$S == 2L, truth$p == 500L, truth$L_f == 1L)
  grid <- as.numeric(truth$t_grid_dense)
  w <- trap_weights(grid)
  design_dense <- cbind(1, grid,
                        multiFSYNC:::ZOSull(grid, range.x = c(0, 1),
                                             intKnots = truth$int_knots))
  roles <- list(
    A1 = list(study = 1L, loading = truth$a_true[, 1L],
              score = truth$eta_true[[1L]][[1L]],
              feature = truth$phi_dense_list[[1L]],
              process_feature = truth$theta_dense_list[[1L]],
              meta = truth$phi_basis_meta[[1L]], type = "shared"),
    A2 = list(study = 2L, loading = truth$a_true[, 1L],
              score = truth$eta_true[[2L]][[1L]],
              feature = truth$phi_dense_list[[1L]],
              process_feature = truth$theta_dense_list[[1L]],
              meta = truth$phi_basis_meta[[1L]], type = "shared"),
    B1 = list(study = 1L, loading = truth$b_true[[1L]][, 1L],
              score = truth$chi_true[[1L]][[1L]],
              feature = truth$psi_dense_list[[1L]][[1L]],
              process_feature = truth$kappa_dense_list[[1L]][[1L]],
              meta = truth$psi_basis_meta[[1L]][[1L]], type = "specific"),
    B2 = list(study = 2L, loading = truth$b_true[[2L]][, 1L],
              score = truth$chi_true[[2L]][[1L]],
              feature = truth$psi_dense_list[[2L]][[1L]],
              process_feature = truth$kappa_dense_list[[2L]][[1L]],
              meta = truth$psi_basis_meta[[2L]][[1L]], type = "specific")
  )
  for (role in names(roles)) {
    item <- roles[[role]]
    s <- item$study
    time <- unlist(data$time_obs[[s]], use.names = FALSE)
    obs_features <- eval_truth_features(item$meta, time)
    dense_check <- eval_truth_features(item$meta, grid)
    stopifnot(max(abs(dense_check - item$feature)) < 1e-8)
    observed_gram <- crossprod(obs_features) / length(time)
    continuous_gram <- crossprod(item$feature, w * item$feature)
    obs_process <- unlist(lapply(seq_along(data$time_obs[[s]]), function(i) {
      if (identical(item$type, "shared"))
        truth$f_true_values[[s]][[i]][, 1L]
      else truth$g_true_values[[s]][[i]][, 1L]
    }), use.names = FALSE)
    score_scaled <- item$score
    dense_process <- item$process_feature %*% t(score_scaled)
    loading_sq <- sum(item$loading^2)
    observed_rms <- sqrt(loading_sq * mean(obs_process^2) / truth$p)
    dense_rms <- sqrt(loading_sq * sum(w * dense_process^2) /
                        (nrow(item$score) * truth$p))
    noise_expected_rms <- sqrt(mean(truth$sigma2_eps_true[s, ]))
    noise_actual_sq <- sum(vapply(truth$noise_true_values[[s]],
                                  function(x) sum(x^2), numeric(1L)))
    noise_actual_rms <- sqrt(noise_actual_sq / (length(time) * truth$p))
    boundary <- grid <= 0.1 | grid >= 0.9
    boundary_fraction <- sum(w[boundary] *
                               rowSums(item$feature[boundary, , drop = FALSE]^2)) /
      sum(w * rowSums(item$feature^2))
    Q <- item$meta$stiefel_coefficients
    roughness <- sum(diag(crossprod(Q, item$meta$second_derivative_gram %*% Q)))
    design_obs <- cbind(1, time,
                        multiFSYNC:::ZOSull(time, range.x = c(0, 1),
                                             intKnots = truth$int_knots))
    truth_rows[[length(truth_rows) + 1L]] <- data.frame(
      logical_data_id = data_id, role = role, study = s,
      observed_values = length(time) * truth$p,
      observed_times = length(time),
      loading_support = sum(abs(item$loading) > 1e-12),
      loading_l2 = sqrt(loading_sq),
      score_variance_1 = stats::var(item$score[, 1L]),
      score_variance_2 = stats::var(item$score[, 2L]),
      observed_gram_min_eigen = min(eigen(observed_gram, symmetric = TRUE)$values),
      observed_gram_max_eigen = max(eigen(observed_gram, symmetric = TRUE)$values),
      continuous_gram_min_eigen = min(eigen(continuous_gram, symmetric = TRUE)$values),
      observed_contribution_rms = observed_rms,
      dense_contribution_rms = dense_rms,
      observed_to_dense_energy_ratio = (observed_rms / dense_rms)^2,
      expected_noise_rms = noise_expected_rms,
      realized_noise_rms = noise_actual_rms,
      observed_snr_rms_expected = observed_rms / noise_expected_rms,
      observed_snr_rms_realized = observed_rms / noise_actual_rms,
      K9_feature_projection_floor_relative =
        project_floor(item$feature, design_dense, w),
      K9_process_projection_floor_relative =
        project_floor(item$process_feature, design_dense, w),
      K9_observed_feature_projection_floor_relative =
        project_floor(obs_features, design_obs, rep(1 / length(time), length(time))),
      boundary_10pct_feature_energy_fraction = boundary_fraction,
      feature_roughness_sum = roughness)
  }
  for (s in 1:2) {
    sizes <- lengths(data$time_obs[[s]])
    design_rows[[length(design_rows) + 1L]] <- data.frame(
      logical_data_id = data_id, study = s, subjects = length(sizes),
      total_times = sum(sizes), min_times = min(sizes),
      median_times = stats::median(sizes), max_times = max(sizes),
      time_min = min(unlist(data$time_obs[[s]])),
      time_max = max(unlist(data$time_obs[[s]])))
  }
  pair_specs <- list(c("A1", "B1"), c("A2", "B2"), c("B1", "B2"))
  for (pair in pair_specs) {
    lhs <- roles[[pair[[1L]]]]; rhs <- roles[[pair[[2L]]]]
    lc <- cosine(lhs$loading, rhs$loading)
    support_l <- abs(lhs$loading) > 1e-12
    support_r <- abs(rhs$loading) > 1e-12
    same_study <- lhs$study == rhs$study
    obs_corr <- if (same_study) {
      tt <- unlist(data$time_obs[[lhs$study]], use.names = FALSE)
      subspace_correlation(eval_truth_features(lhs$meta, tt),
                           eval_truth_features(rhs$meta, tt),
                           rep(1 / length(tt), length(tt)))
    } else NA_real_
    geometry_rows[[length(geometry_rows) + 1L]] <- data.frame(
      logical_data_id = data_id, pair = paste(pair, collapse = "-"),
      loading_abs_cosine = abs(lc),
      loading_angle_degrees = acos(min(1, abs(lc))) * 180 / pi,
      support_overlap = sum(support_l & support_r),
      support_jaccard = sum(support_l & support_r) /
        sum(support_l | support_r),
      continuous_function_canonical_correlation =
        subspace_correlation(lhs$feature, rhs$feature, w),
      observed_function_canonical_correlation = obs_corr)
  }
  # Full observed contribution directions in study 2, including beta.
  contributions <- lapply(seq_along(data$time_obs[[2L]]), function(i) {
    a <- tcrossprod(truth$f_true_values[[2L]][[i]][, 1L],
                    truth$a_true[, 1L])
    b <- tcrossprod(truth$g_true_values[[2L]][[i]][, 1L],
                    truth$b_true[[2L]][, 1L])
    cbind(as.vector(a), as.vector(b),
          as.vector(truth$mu_true_values[[2L]][[i]]),
          as.vector(truth$beta_true_values[[2L]][[i]]))
  })
  contribution_matrix <- do.call(rbind, contributions)
  colnames(contribution_matrix) <- c("A2", "B2", "mean", "beta")
  b2 <- contribution_matrix[, "B2"]
  q <- qr(contribution_matrix[, c("A2", "mean", "beta"), drop = FALSE])
  b2_residual <- qr.resid(q, b2)
  contribution_rows[[length(contribution_rows) + 1L]] <- data.frame(
    logical_data_id = data_id,
    B2_A2_observed_contribution_cosine =
      cosine(b2, contribution_matrix[, "A2"]),
    B2_mean_observed_contribution_cosine =
      cosine(b2, contribution_matrix[, "mean"]),
    B2_beta_observed_contribution_cosine =
      cosine(b2, contribution_matrix[, "beta"]),
    B2_residual_energy_fraction_after_A2_mean_beta =
      sum(b2_residual^2) / sum(b2^2))
  # Verify that the truth knots reproduce the frozen fitted K9 design.
  fit_path <- file.path(input_root, "fits",
                          paste0(data_id, "__C1__start06"), "fit.rds")
  reference_fit <- readRDS(fit_path)
  reconstructed <- cbind(1, reference_fit$time_g,
                         multiFSYNC:::ZOSull(reference_fit$time_g,
                                              range.x = c(0, 1),
                                              intKnots = truth$int_knots))
  binding_rows[[length(binding_rows) + 1L]] <- data.frame(
    logical_data_id = data_id,
    max_abs_K9_design_difference = max(abs(reconstructed - reference_fit$C_g)),
    C0_C1_nonbeta_truth_identical = TRUE)
  stopifnot(max(abs(reconstructed - reference_fit$C_g)) < 1e-10)
  rm(c1, c0, data, truth, roles, reference_fit, contribution_matrix)
  gc(verbose = FALSE)
}

trace_rows <- list(); fit_rows <- list()
fit_manifest <- manifest[manifest$beta_arm %in% c("C0", "C1"), , drop = FALSE]
for (k in seq_len(nrow(fit_manifest))) {
  entry <- fit_manifest[k, , drop = FALSE]
  fit_id <- entry$fit_id[[1L]]
  fit <- readRDS(file.path(input_root, "fits", fit_id, "fit.rds"))
  science_row <- science[science$fit_id == fit_id, , drop = FALSE]
  stopifnot(nrow(science_row) == 1L, fit$K == 9L,
            fit$fit_seed == entry$fit_seed[[1L]],
            is.data.frame(fit$driver_trace),
            all(c(1L, 380L, 400L) %in% fit$driver_trace$iteration))
  tr <- fit$driver_trace
  keep <- tr$block %in% c("shared_loading", "specific_loading",
                         "specific_score", "specific_time_function")
  selected <- tr[keep, , drop = FALSE]
  selected$fit_id <- fit_id
  selected$logical_data_id <- entry$logical_data_id[[1L]]
  selected$beta_arm <- entry$beta_arm[[1L]]
  selected$start_index <- entry$start_index[[1L]]
  selected$core_truth_recovered <- science_row$core_truth_recovered[[1L]]
  selected$truth_free_winner <- science_row$truth_free_winner[[1L]]
  trace_rows[[length(trace_rows) + 1L]] <- selected
  pp <- fit$mu_q_gamma_b[[2L]][, 1L]
  load <- fit$mu_q_b_specific[[2L]][, 1L]
  coefficient <- fit$mu_q_nu_psi[[2L]][[1L]]
  score <- fit$mu_q_xi[[2L]][[1L]]
  time_values <- fit$C_g %*% coefficient
  weights <- trap_weights(fit$time_g)
  calibration <- fit$random_scale_calibration_diagnostics
  calibration <- calibration[calibration$block == "specific_function_l2" &
                               calibration$study == 2L, , drop = FALSE]
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(
    fit_id = fit_id, logical_data_id = entry$logical_data_id[[1L]],
    beta_arm = entry$beta_arm[[1L]], start_index = entry$start_index[[1L]],
    truth_free_winner = science_row$truth_free_winner[[1L]],
    core_truth_recovered = science_row$core_truth_recovered[[1L]],
    selected_specific_study_2 = science_row$selected_specific_study_2[[1L]],
    initial_B2_function_l2_before_component_1 = calibration$before[calibration$component == 1L],
    initial_B2_function_l2_before_component_2 = calibration$before[calibration$component == 2L],
    pre_score_B2_data_coherence =
      get_scalar(tr, "specific_score", 2L, 1L, "data_coherence", "pre_score"),
    B2_PPI_mean_sweep_1 =
      get_scalar(tr, "specific_loading", 2L, 1L, "ppi_mean"),
    B2_PPI_mean_sweep_380 =
      get_scalar(tr, "specific_loading", 2L, 380L, "ppi_mean"),
    B2_PPI_mean_sweep_400 =
      get_scalar(tr, "specific_loading", 2L, 400L, "ppi_mean"),
    B2_loading_coefficient_norm_sweep_1 =
      get_scalar(tr, "specific_loading", 2L, 1L, "posterior_mean_norm"),
    B2_loading_coefficient_norm_sweep_380 =
      get_scalar(tr, "specific_loading", 2L, 380L, "posterior_mean_norm"),
    B2_loading_coefficient_norm_sweep_400 =
      get_scalar(tr, "specific_loading", 2L, 400L, "posterior_mean_norm"),
    B2_time_coefficient_norm_sweep_1 = get_time_norm(tr, 2L, 1L),
    B2_time_coefficient_norm_sweep_380 = get_time_norm(tr, 2L, 380L),
    B2_time_coefficient_norm_sweep_400 = get_time_norm(tr, 2L, 400L),
    B2_final_sum_loading_PPI = sum(pp),
    B2_final_loading_count_PPI_ge_half = sum(pp >= 0.5),
    B2_final_raw_loading_norm = sqrt(sum(load^2)),
    B2_final_time_function_L2 = sqrt(sum(weights * time_values^2)),
    B2_final_score_RMS = sqrt(mean(score^2)),
    ordinary_ELBO_first_T1 = fit$ELBO[[1L]],
    ordinary_ELBO_final = tail(fit$ELBO, 1L))
  if (k %% 12L == 0L) cat("Read", k, "of", nrow(fit_manifest), "fits\n")
  rm(fit, tr)
  gc(verbose = FALSE)
}

truth_table <- do.call(rbind, truth_rows)
design_table <- do.call(rbind, design_rows)
geometry_table <- do.call(rbind, geometry_rows)
contribution_table <- do.call(rbind, contribution_rows)
binding_table <- do.call(rbind, binding_rows)
trace_table <- do.call(rbind, trace_rows)
fit_table <- do.call(rbind, fit_rows)

summaries <- list()
for (logical_data_id in data_ids) for (arm in c("C0", "C1")) {
  z <- fit_table[fit_table$logical_data_id == logical_data_id &
                   fit_table$beta_arm == arm, , drop = FALSE]
  stopifnot(nrow(z) == 12L)
  summaries[[length(summaries) + 1L]] <- data.frame(
    logical_data_id = logical_data_id, beta_arm = arm,
    core_hits = sum(z$core_truth_recovered),
    B2_selected_count = sum(z$selected_specific_study_2 == 1L),
    B2_first_sweep_PPI_median = median(z$B2_PPI_mean_sweep_1),
    B2_first_sweep_PPI_below_0_1 = sum(z$B2_PPI_mean_sweep_1 < 0.1),
    B2_first_sweep_loading_norm_median =
      median(z$B2_loading_coefficient_norm_sweep_1),
    B2_sweep_380_PPI_median = median(z$B2_PPI_mean_sweep_380),
    B2_sweep_400_PPI_median = median(z$B2_PPI_mean_sweep_400),
    B2_final_loading_norm_median = median(z$B2_final_raw_loading_norm))
}
summary_table <- do.call(rbind, summaries)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(output_root))
utils::write.csv(truth_table, file.path(output_root, "TRUTH_BLOCK_VISIBILITY.csv"),
                 row.names = FALSE)
utils::write.csv(design_table, file.path(output_root, "OBSERVATION_DESIGN.csv"),
                 row.names = FALSE)
utils::write.csv(geometry_table, file.path(output_root, "BLOCK_PAIR_GEOMETRY.csv"),
                 row.names = FALSE)
utils::write.csv(contribution_table,
                 file.path(output_root, "B2_CONTRIBUTION_COMPETITION.csv"),
                 row.names = FALSE)
utils::write.csv(binding_table, file.path(output_root, "DESIGN_BINDING_CHECK.csv"),
                 row.names = FALSE)
utils::write.csv(trace_table, file.path(output_root, "DRIVER_TRACE_SELECTED.csv"),
                 row.names = FALSE)
utils::write.csv(fit_table, file.path(output_root, "FIT_B2_TRAJECTORY.csv"),
                 row.names = FALSE)
utils::write.csv(summary_table, file.path(output_root, "POOL_B2_TRAJECTORY.csv"),
                 row.names = FALSE)
writeLines(c(
  "status=STAGE7D_READONLY_ATTRIBUTION_COMPLETE",
  "new_fit_or_continuation_started=FALSE",
  "existing_stage7c_files_modified=FALSE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  paste0("truth_block_rows=", nrow(truth_table)),
  paste0("fit_rows=", nrow(fit_table)),
  paste0("driver_trace_rows=", nrow(trace_table)),
  "driver_trace_sweeps=1,380,400",
  "trace_gap=No saved intermediate annealing or early ordinary sweeps"
), file.path(output_root, "ANALYSIS_COMPLETE.txt"))
cat("STAGE7D_READONLY_ATTRIBUTION_COMPLETE\n")
