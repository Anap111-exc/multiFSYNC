#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- sub("^--file=", "", script_arg[[1L]])
script_root <- dirname(script_path)
experiment_id <- "COVARIATE_ROUTEB_PILOT_V1_20260927"
`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
save_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3)
}
csv_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
lines_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
}
read_manifest <- function(path) {
  assert(file.exists(path), paste0("Missing manifest: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
parse_cli <- function(x) {
  ans <- list()
  for (item in x) {
    assert(grepl("^--[^=]+=", item), "Use --key=value arguments.")
    key <- sub("^--([^=]+)=.*$", "\\1", item)
    assert(is.null(ans[[key]]), paste0("Duplicate argument: ", key))
    ans[[key]] <- sub("^[^=]+=", "", item)
  }
  ans
}

data_manifest <- function(smoke = FALSE) {
  data.frame(
    logical_data_id = rep("covrbp_01", 3L),
    data_id = paste0("covrbp_01__", c("C0", "C1", "C2")),
    data_index = rep(1L, 3L),
    beta_arm = c("C0", "C1", "C2"),
    observation_source_data_id = c(
      "covrbp_01__C0", "covrbp_01__C1", "covrbp_01__C1"),
    data_seed = rep(99471001L, 3L),
    truth_basis_seed = rep(99472001L, 3L),
    truth_loading_seed = rep(99473001L, 3L),
    mean_truth_seed = rep(99474001L, 3L),
    beta_bspline_seed = rep(99475001L, 3L),
    d_generated = rep(2L, 3L),
    d_fit = c(2L, 2L, 0L),
    truth_available_to_fit = FALSE,
    truth_available_to_selection = FALSE,
    truth_unseal_authorized = FALSE,
    stringsAsFactors = FALSE
  )
}

fit_manifest <- function(data_rows, smoke = FALSE) {
  starts <- if (smoke) 1L else seq_len(12L)
  out <- do.call(rbind, lapply(seq_len(nrow(data_rows)), function(i) {
    row <- data_rows[i, , drop = FALSE]
    data.frame(
      fit_id = sprintf("%s__start%02d", row$data_id[[1L]], starts),
      data_id = row$data_id[[1L]],
      logical_data_id = row$logical_data_id[[1L]],
      beta_arm = row$beta_arm[[1L]],
      start_index = starts,
      fit_seed = 99476000L + starts,
      n_cpus = 1L, truth_available_to_fit = FALSE,
      truth_available_to_selection = FALSE,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

generator_config <- function(smoke = FALSE) {
  list(
    S = 2L, n_s = if (smoke) c(4L, 4L) else c(30L, 30L),
    p = if (smoke) 20L else 500L, d = 2L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
    K = if (smoke) 3L else 9L,
    n_obs = if (smoke) 5L else 8L, common_grid = FALSE,
    n_obs_range = if (smoke) c(4L, 6L) else c(6L, 9L),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE,
    prop_sparse = if (smoke) 0.8 else 0.9,
    score_var_decay = TRUE, n_dense = if (smoke) 51L else 201L,
    mean_amp = 0.6, beta_amp = 0.3, bs_degree = c(2L, 3L),
    identified_loadings = TRUE, sparsity_mode = "fixed",
    truth_basis_method = "iid_common_splinet_stiefel",
    truth_basis_dim = 8L,
    mean_truth_method = "independent_splinet",
    mean_basis_dim = 8L,
    mean_smoothness_alpha = 1.519826404906808
  )
}

calibrate_loadings <- function(data) {
  truth <- data$true_params
  calibrate <- function(x) {
    active <- sum(x != 0)
    norm <- sqrt(sum(x^2))
    assert(active >= 1L && is.finite(norm) && norm > 1e-14,
           "Empty or non-finite loading block.")
    3 * sqrt(active / 10) * x / norm
  }
  truth$a_true[, 1L] <- calibrate(truth$a_true[, 1L])
  for (s in seq_len(truth$S)) {
    truth$b_true[[s]][, 1L] <- calibrate(truth$b_true[[s]][, 1L])
  }
  truth$loading_calibration <- "v0f_constant_per_active_loading"
  for (s in seq_len(truth$S)) for (i in seq_len(truth$n_s[[s]])) {
    signal <- truth$mu_true_values[[s]][[i]] +
      truth$beta_true_values[[s]][[i]] +
      truth$f_true_values[[s]][[i]] %*% t(truth$a_true) +
      truth$g_true_values[[s]][[i]] %*% t(truth$b_true[[s]])
    truth$signal_true_values[[s]][[i]] <- signal
    observed <- signal + truth$noise_true_values[[s]][[i]]
    data$Y[[s]][[i]] <- lapply(seq_len(truth$p), function(j) observed[, j])
  }
  data$true_params <- truth
  data
}

bspline_design <- function(t) {
  splines::bs(t, degree = 3L, knots = c(0.25, 0.5, 0.75),
              Boundary.knots = c(0, 1), intercept = TRUE)
}

paired_beta_data <- function(row, smoke = FALSE) {
  config <- generator_config(smoke)
  base <- do.call(multiFSYNC::simulate_multi_study_structured,
                  c(config, list(
                    seed = row$data_seed[[1L]],
                    truth_basis_seed = row$truth_basis_seed[[1L]],
                    truth_loading_seed = row$truth_loading_seed[[1L]],
                    mean_truth_seed = row$mean_truth_seed[[1L]])))
  base <- calibrate_loadings(base)
  grid <- base$true_params$t_grid_dense
  phase <- base$true_params$phase_beta
  p <- base$true_params$p
  d <- 2L
  reference_grid <- array(NA_real_, dim = c(length(grid), p, d))
  for (r in seq_len(d)) for (j in seq_len(p)) {
    phase_index <- (j - 1L) * d + r
    reference_grid[, j, r] <-
      0.3 * sin(pi * r * grid + phase[[phase_index]])
  }

  B <- bspline_design(grid)
  assert(ncol(B) == 7L, "Unexpected B-spline basis dimension.")
  set.seed(row$beta_bspline_seed[[1L]])
  coefficient <- array(rnorm(7L * p * d), dim = c(7L, p, d))
  weights <- rep(1 / (length(grid) - 1L), length(grid))
  weights[c(1L, length(grid))] <- weights[[1L]] / 2
  weighted_mean <- function(x) colSums(sweep(x, 1L, weights, `*`))
  target_mean <- raw_mean <- scale <- matrix(NA_real_, p, d)
  bspline_grid <- array(NA_real_, dim = c(length(grid), p, d))
  for (r in seq_len(d)) {
    raw <- B %*% coefficient[, , r]
    target_mean[, r] <- weighted_mean(reference_grid[, , r])
    raw_mean[, r] <- weighted_mean(raw)
    target_sd <- sqrt(weighted_mean(sweep(
      reference_grid[, , r], 2L, target_mean[, r], `-`)^2))
    raw_sd <- sqrt(weighted_mean(sweep(raw, 2L, raw_mean[, r], `-`)^2))
    assert(all(is.finite(raw_sd)) && all(raw_sd > 1e-10),
           "Degenerate B-spline beta draw.")
    scale[, r] <- target_sd / raw_sd
    value <- sweep(sweep(raw, 2L, raw_mean[, r], `-`),
                   2L, scale[, r], `*`)
    bspline_grid[, , r] <- sweep(value, 2L, target_mean[, r], `+`)
  }

  c1 <- base
  truth <- c1$true_params
  for (s in seq_len(truth$S)) for (i in seq_len(truth$n_s[[s]])) {
    B_obs <- bspline_design(c1$time_obs[[s]][[i]])
    new_contribution <- matrix(0, nrow(B_obs), p)
    for (r in seq_len(d)) {
      raw_obs <- B_obs %*% coefficient[, , r]
      beta_obs <- sweep(sweep(raw_obs, 2L, raw_mean[, r], `-`),
                        2L, scale[, r], `*`)
      beta_obs <- sweep(beta_obs, 2L, target_mean[, r], `+`)
      new_contribution <- new_contribution +
        c1$Z[[s]][i, r] * beta_obs
    }
    signal <- truth$signal_true_values[[s]][[i]] -
      truth$beta_true_values[[s]][[i]] + new_contribution
    truth$beta_true_values[[s]][[i]] <- new_contribution
    truth$signal_true_values[[s]][[i]] <- signal
    observed <- signal + truth$noise_true_values[[s]][[i]]
    c1$Y[[s]][[i]] <- lapply(seq_len(p), function(j) observed[, j])
  }
  truth$beta_truth_method <-
    "independent_cubic_bspline_d2_matched_mean_centered_energy"
  truth$beta_true_grid <- bspline_grid
  truth$beta_reference_grid <- reference_grid
  truth$beta_bspline_coefficients <- coefficient
  truth$beta_bspline_seed <- row$beta_bspline_seed[[1L]]
  truth$beta_bspline_raw_mean <- raw_mean
  truth$beta_bspline_scale <- scale
  truth$beta_bspline_target_mean <- target_mean
  c1$true_params <- truth

  c0 <- c1
  truth0 <- c0$true_params
  for (s in seq_len(truth0$S)) for (i in seq_len(truth0$n_s[[s]])) {
    zero <- truth0$beta_true_values[[s]][[i]]
    zero[,] <- 0
    signal <- truth0$signal_true_values[[s]][[i]] -
      truth0$beta_true_values[[s]][[i]]
    truth0$beta_true_values[[s]][[i]] <- zero
    truth0$signal_true_values[[s]][[i]] <- signal
    observed <- signal + truth0$noise_true_values[[s]][[i]]
    c0$Y[[s]][[i]] <- lapply(seq_len(p), function(j) observed[, j])
  }
  truth0$beta_truth_method <- "zero_effect_d2"
  truth0$beta_true_grid <- array(0, dim = dim(bspline_grid))
  c0$true_params <- truth0
  list(C0 = c0, C1 = c1)
}

load_package <- function() {
  if (identical(Sys.getenv("COVARIATE_ROUTEB_USE_LOAD_ALL"), "TRUE")) {
    assert(requireNamespace("pkgload", quietly = TRUE),
           "pkgload is unavailable for the local smoke.")
    root <- file.path(script_root, "..", "..", "..")
    pkgload::load_all(root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
  }
  assert(requireNamespace("multiFSYNC", quietly = TRUE),
         "multiFSYNC is not installed or loaded.")
  assert(requireNamespace("Splinets", quietly = TRUE),
         "Splinets is required for the registered truth generator.")
  invisible(TRUE)
}

pair_checks <- function(pair, index) {
  c0 <- pair$C0
  c1 <- pair$C1
  t0 <- c0$true_params
  t1 <- c1$true_params
  grid <- t1$t_grid_dense
  w <- rep(1 / (length(grid) - 1L), length(grid))
  w[c(1L, length(grid))] <- w[[1L]] / 2
  weighted_mean <- function(x) colSums(sweep(x, 1L, w, `*`))
  mean_error <- rms_error <- numeric(2L)
  for (r in seq_len(2L)) {
    reference <- t1$beta_reference_grid[, , r]
    value <- t1$beta_true_grid[, , r]
    m_reference <- weighted_mean(reference)
    m_value <- weighted_mean(value)
    sd_reference <- sqrt(weighted_mean(
      sweep(reference, 2L, m_reference, `-`)^2))
    sd_value <- sqrt(weighted_mean(
      sweep(value, 2L, m_value, `-`)^2))
    mean_error[[r]] <- max(abs(m_reference - m_value))
    rms_error[[r]] <- max(abs(sd_reference - sd_value))
  }
  errors <- numeric()
  for (s in seq_len(t1$S)) for (i in seq_len(t1$n_s[[s]])) {
    obs0 <- do.call(cbind, c0$Y[[s]][[i]])
    obs1 <- do.call(cbind, c1$Y[[s]][[i]])
    errors <- c(errors,
      max(abs(obs0 - t0$signal_true_values[[s]][[i]] -
              t0$noise_true_values[[s]][[i]])),
      max(abs(obs1 - t1$signal_true_values[[s]][[i]] -
              t1$noise_true_values[[s]][[i]])),
      max(abs((obs1 - obs0) - t1$beta_true_values[[s]][[i]])))
  }
  checks <- c(
    identical_observation_times = identical(c0$time_obs, c1$time_obs),
    identical_covariates = identical(c0$Z, c1$Z),
    identical_design = identical(c0$C, c1$C),
    identical_loadings = identical(t0$a_true, t1$a_true) &&
      identical(t0$b_true, t1$b_true),
    identical_nonbeta_truth = all(vapply(c(
      "mu_true_values", "f_true_values", "g_true_values",
      "noise_true_values", "theta_dense_list", "kappa_dense_list",
      "zeta_true", "xi_true"), function(name) {
        identical(t0[[name]], t1[[name]])
      }, logical(1L))),
    route_b_truth = identical(
      t1$truth_basis_method, "iid_common_splinet_stiefel"),
    periodic_bspline_mean = identical(
      t1$mean_truth_method, "independent_splinet"),
    controlled_identifiable_loadings = isTRUE(t1$identified_loadings),
    c0_beta_exactly_zero = max(abs(t0$beta_true_grid)) == 0,
    c1_beta_finite = all(is.finite(t1$beta_true_grid)),
    matched_beta_continuous_mean = max(mean_error) < 1e-12,
    matched_beta_centered_energy = max(rms_error) < 1e-12,
    response_decomposition = max(errors) < 1e-10,
    beta_curves_differ_from_reference = max(abs(
      t1$beta_true_grid - t1$beta_reference_grid)) > 0.05
  )
  data.frame(data_index = index, check_id = names(checks),
             passed = unname(checks), stringsAsFactors = FALSE)
}

save_arm <- function(data, row, output_root, smoke) {
  root <- file.path(output_root, "data", row$data_id[[1L]])
  assert(!dir.exists(root), paste0("Data directory exists: ", root))
  dir.create(file.path(root, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  dimensions <- generator_config(smoke)[c(
    "S", "n_s", "p", "d", "L_f", "L_s", "M_f", "M_s", "K")]
  dimensions$d <- as.integer(row$d_fit[[1L]])
  fit_z <- if (dimensions$d == 0L) NULL else data$Z
  observation <- list(
    bundle_class = "covariate_routeb_pilot_observation_only",
    experiment_id = experiment_id,
    data_id = row$data_id[[1L]],
    beta_arm = row$beta_arm[[1L]],
    observation_source_data_id = row$observation_source_data_id[[1L]],
    Y = data$Y, Z = fit_z, time_obs = data$time_obs,
    dimensions = dimensions,
    generated_d = as.integer(row$d_generated[[1L]]),
    truth_available_to_fit = FALSE)
  sealed <- list(
    bundle_class = "covariate_routeb_pilot_sealed_truth",
    experiment_id = experiment_id,
    data_id = row$data_id[[1L]],
    data_seed = row$data_seed[[1L]],
    truth_basis_seed = row$truth_basis_seed[[1L]],
    truth_loading_seed = row$truth_loading_seed[[1L]],
    mean_truth_seed = row$mean_truth_seed[[1L]],
    beta_bspline_seed = row$beta_bspline_seed[[1L]],
    truth_unseal_authorized = FALSE,
    data = data)
  save_new(observation, file.path(root, "observation_bundle.rds"))
  save_new(sealed, file.path(root, "sealed_truth", "truth_bundle.rds"))
  lines_new(c("status=DATA_SEAL_COMPLETE", "truth_unseal_authorized=FALSE",
              "scientific_evaluation_started=FALSE"),
            file.path(root, "DATA_SEAL_COMPLETE.txt"))
}

generate <- function(output_root, smoke = FALSE) {
  load_package()
  assert(!dir.exists(output_root),
         paste0("Refusing to overwrite output root: ", output_root))
  dir.create(output_root, recursive = TRUE, mode = "0700")
  data_rows <- data_manifest(smoke)
  fit_rows <- fit_manifest(data_rows, smoke)
  assert(!anyDuplicated(data_rows$data_id) &&
           !anyDuplicated(fit_rows$fit_id) &&
           nrow(fit_rows) == if (smoke) 3L else 36L,
         "Manifest dimensions or identifiers are invalid.")
  csv_new(data_rows, file.path(output_root, "DATA_MANIFEST.csv"))
  csv_new(fit_rows, file.path(output_root, "FIT_MANIFEST.csv"))
  checks <- list()
  for (index in unique(data_rows$data_index)) {
    row <- data_rows[data_rows$data_index == index, , drop = FALSE]
    assert(nrow(row) == 3L &&
             identical(as.character(row$beta_arm), c("C0", "C1", "C2")),
           "Paired data manifest is malformed.")
    pair <- paired_beta_data(row[1L, , drop = FALSE], smoke)
    qc <- pair_checks(pair, index)
    assert(all(qc$passed), paste0("Pair QC failed for dataset ", index,
                                  ": ", paste(qc$check_id[!qc$passed],
                                               collapse = ";")))
    checks[[length(checks) + 1L]] <- qc
    save_arm(pair$C0, row[row$beta_arm == "C0", , drop = FALSE],
             output_root, smoke)
    save_arm(pair$C1, row[row$beta_arm == "C1", , drop = FALSE],
             output_root, smoke)
    save_arm(pair$C1, row[row$beta_arm == "C2", , drop = FALSE],
             output_root, smoke)
    obs1 <- readRDS(file.path(
      output_root, "data", "covrbp_01__C1", "observation_bundle.rds"))
    obs2 <- readRDS(file.path(
      output_root, "data", "covrbp_01__C2", "observation_bundle.rds"))
    c2_checks <- c(
      c1_c2_identical_Y = identical(obs1$Y, obs2$Y),
      c1_c2_identical_time = identical(obs1$time_obs, obs2$time_obs),
      c1_has_d2_Z = is.list(obs1$Z) && obs1$dimensions$d == 2L,
      c2_omits_Z = is.null(obs2$Z) && obs2$dimensions$d == 0L)
    checks[[length(checks) + 1L]] <- data.frame(
      data_index = index, check_id = names(c2_checks),
      passed = unname(c2_checks), stringsAsFactors = FALSE)
  }
  csv_new(do.call(rbind, checks), file.path(output_root, "PAIRING_QC.csv"))
  lines_new(c("status=DATA_GENERATION_AND_SEAL_COMPLETE",
              paste0("logical_datasets=", length(unique(data_rows$data_index))),
              paste0("observation_bundles=", nrow(data_rows)),
              paste0("planned_fits=", nrow(fit_rows)),
              "truth_unseal_authorized=FALSE",
              "scientific_evaluation_started=FALSE"),
            file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("GENERATE_PASS data=", nrow(data_rows), " fits=", nrow(fit_rows),
      " smoke=", smoke, "\n", sep = "")
  invisible(TRUE)
}

objective_status <- function(fit) {
  finite_elbo <- length(fit$ELBO) > 0L && all(is.finite(fit$ELBO))
  finite_components <- length(fit$ELBO_finite) == length(fit$ELBO) &&
    all(vapply(fit$ELBO_finite, function(value) {
      is.list(value) && isTRUE(value$all)
    }, logical(1L)))
  checks <- c(
    annealing_completed = isTRUE(fit$annealing_completed),
    final_temperature_one = isTRUE(all.equal(
      as.numeric(fit$final_temperature), 1)),
    ordinary_elbo_present_and_finite = finite_elbo,
    ordinary_elbo_components_finite = finite_components,
    ordinary_objective_valid = isTRUE(fit$elbo_objective_valid),
    no_T1_elbo_decrease = identical(as.integer(fit$ELBO_decrease_count), 0L),
    no_T1_solver_jitter = identical(as.integer(fit$elbo_t1_jitter_count), 0L))
  list(eligible = all(checks), checks = checks,
       invalid_reasons = names(checks)[!checks])
}

capture_fit <- function(observation, row, mode = "main") {
  smoke <- identical(mode, "smoke")
  benchmark <- startsWith(mode, "benchmark")
  if (smoke) {
    anneal <- c(1, 1.2, 3); maxit <- 7L; n_g <- 31L
    convergence_rule <- "parameters"; practical_control <- NULL
    trace_sweeps <- 1L
  } else if (benchmark) {
    anneal <- NULL
    maxit <- as.integer(sub("^benchmark", "", mode))
    n_g <- 51L; convergence_rule <- "parameters"
    practical_control <- NULL; trace_sweeps <- 1L
  } else {
    anneal <- c(1, 1.9, 100); maxit <- 499L; n_g <- 51L
    convergence_rule <- "practical"
    practical_control <- multiFSYNC::g12_stopping_control()
    practical_control$checkpoints <- c(380L, 400L)
    trace_sweeps <- c(1L, 380L, 400L)
  }
  p <- observation$dimensions$p
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(
    multiFSYNC::bayesSYNC_multi_pre_score(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
      K = observation$dimensions$K, n_g = n_g,
      anneal = anneal, maxit = maxit,
      n_cpus = 1L, verbose = FALSE, seed = row$fit_seed[[1L]],
      tol_abs = if (smoke || benchmark) 0 else 1e-3,
      tol_rel = if (smoke || benchmark) 0 else 1e-5,
      bool_scale = FALSE, bool_var_spec_prob = FALSE,
      d_0 = as.integer(p), convergence_rule = convergence_rule,
      lambda_orth = 0, practical_control = practical_control,
      initialization = "random",
      initialization_control = list(
        grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8),
      continuation_state = NULL, pre_score_sweeps = 1L,
      trace_sweeps = trace_sweeps,
      function_initialization = "gram_unit_energy"
    ), warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  error <- if (inherits(value, "error")) conditionMessage(value) else NULL
  fit <- if (is.null(error)) value else NULL
  if (is.null(error) && !smoke && !benchmark) {
    valid <- tryCatch({
      assert(fit$annealing_sweeps == 99L && fit$t1_sweeps == 400L &&
               fit$n_g == 51L,
             "Incomplete fixed-400 fit endpoint.")
      assert(fit$n_cpus_used == 1L &&
               fit$pre_score_interface$pre_score_sweeps == 1L &&
               identical(fit$pre_score_interface$function_initialization,
                         "gram_unit_energy"),
             "Route provenance differs from protocol.")
      TRUE
    }, error = function(e) e)
    if (inherits(valid, "error")) {
      error <- conditionMessage(valid)
      fit <- NULL
    }
  }
  list(fit = fit, error = error, warnings = warnings,
       elapsed_seconds = elapsed)
}

fit_one <- function(output_root, fit_id, smoke = FALSE) {
  load_package()
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  assert(nrow(row) == 1L, "fit_id must match exactly one row.")
  data_root <- file.path(output_root, "data", row$data_id[[1L]])
  assert(file.exists(file.path(data_root, "DATA_SEAL_COMPLETE.txt")),
         "Sealed dataset is incomplete.")
  observation <- readRDS(file.path(data_root, "observation_bundle.rds"))
  assert(identical(
           observation$bundle_class,
           "covariate_routeb_pilot_observation_only") &&
           isFALSE(observation$truth_available_to_fit),
         "Fit input is not observation-only.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  assert(!dir.exists(fit_dir), paste0("Refusing fit rerun: ", fit_id))
  dir.create(fit_dir, recursive = TRUE, mode = "0700")
  outcome <- capture_fit(observation, row, if (smoke) "smoke" else "main")
  fit <- outcome$fit
  objective <- if (is.null(fit)) {
    list(eligible = FALSE, invalid_reasons = "fit_error")
  } else if (smoke) {
    list(eligible = all(is.finite(fit$ELBO)), invalid_reasons = character())
  } else objective_status(fit)
  final_elbo <- if (is.null(fit) || !length(fit$ELBO)) NA_real_ else
    tail(fit$ELBO, 1L)
  terminal <- list(
    experiment_id = experiment_id,
    fit_id = fit_id, data_id = row$data_id[[1L]],
    logical_data_id = row$logical_data_id[[1L]],
    beta_arm = row$beta_arm[[1L]], start_index = row$start_index[[1L]],
    fit_seed = row$fit_seed[[1L]], n_cpus = 1L,
    status = if (is.null(fit)) "ERROR" else "COMPLETE",
    error = outcome$error,
    warnings = outcome$warnings,
    elapsed_seconds = outcome$elapsed_seconds,
    objective_eligible = objective$eligible,
    invalid_reasons = paste(objective$invalid_reasons, collapse = ";"),
    final_ordinary_elbo = final_elbo,
    annealing_sweeps = if (is.null(fit)) NA_integer_ else
      fit$annealing_sweeps,
    ordinary_t1_sweeps = if (is.null(fit)) NA_integer_ else fit$t1_sweeps,
    strict_practical_converged = if (is.null(fit)) FALSE else
      isTRUE(fit$practical_converged) && isTRUE(fit$converged) &&
      identical(as.character(fit$convergence_status),
                "converged_practical"),
    truth_used_for_fit_or_selection = FALSE)
  save_new(terminal, file.path(fit_dir, "terminal_record.rds"))
  csv_new(data.frame(message = outcome$warnings),
          file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(fit)) {
    lines_new(c("status=ERROR_RETAINED_NO_SELECTIVE_RERUN",
                paste0("fit_id=", fit_id),
                paste0("message=", outcome$error)),
              file.path(fit_dir, "FIT_ERROR.txt"))
    stop("Fit failed: ", fit_id, ": ", outcome$error, call. = FALSE)
  }
  save_new(fit, file.path(fit_dir, "fit.rds"))
  csv_new(data.frame(t1_sweep = seq_along(fit$ELBO),
                     ordinary_elbo = as.numeric(fit$ELBO)),
          file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  lines_new(c("status=COMPLETE", paste0("fit_id=", fit_id),
              paste0("objective_eligible=", objective$eligible),
              paste0("elapsed_seconds=", outcome$elapsed_seconds),
              "n_cpus=1", "truth_used_for_fit_or_selection=FALSE"),
            file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id,
      " elapsed_seconds=", outcome$elapsed_seconds,
      " eligible=", objective$eligible, "\n", sep = "")
  invisible(terminal)
}

benchmark <- function(output_root) {
  load_package()
  marker <- file.path(output_root, "DATA_GENERATION_COMPLETE.txt")
  assert(file.exists(marker), "Generate and seal before benchmarking.")
  benchmark_root <- file.path(output_root, "preflight")
  assert(!dir.exists(benchmark_root), "Refusing benchmark rerun.")
  dir.create(benchmark_root, recursive = TRUE, mode = "0700")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  row <- manifest[manifest$data_id == "covrbp_01__C1" &
                    manifest$start_index == 1L, , drop = FALSE]
  assert(nrow(row) == 1L, "Benchmark row missing.")
  observation <- readRDS(file.path(
    output_root, "data", "covrbp_01__C1", "observation_bundle.rds"))
  assert(observation$dimensions$p == 500L &&
           observation$dimensions$d == 2L &&
           observation$dimensions$K == 9L,
         "Benchmark must use full-size covariate data.")
  short <- capture_fit(observation, row, "benchmark3")
  long <- capture_fit(observation, row, "benchmark8")
  assert(is.null(short$error) && is.null(long$error), paste0(
    "Full-size benchmark fit failed: ", short$error %||% long$error))
  assert(short$fit$t1_sweeps == 3L && long$fit$t1_sweeps == 8L &&
           isTRUE(all.equal(short$fit$ELBO, long$fit$ELBO[1:3],
                            tolerance = 1e-10)),
         "Benchmark ordinary T=1 trajectory is not reproducible.")
  marginal <- (long$elapsed_seconds - short$elapsed_seconds) / 5
  assert(is.finite(marginal) && marginal > 0,
         "Benchmark cannot estimate a positive marginal sweep time.")
  result <- data.frame(
    scenario = "covrbp_01__C1", p = 500L, d = 2L, K = 9L,
    short_t1_sweeps = 3L, long_t1_sweeps = 8L,
    short_elapsed_seconds = short$elapsed_seconds,
    long_elapsed_seconds = long$elapsed_seconds,
    estimated_marginal_seconds_per_T1_sweep = marginal,
    projected_seconds_per_499_sweep_fit = 499 * marginal,
    projected_hours_36_fits_8_workers = 36 * 499 * marginal / 8 / 3600,
    stringsAsFactors = FALSE)
  csv_new(result, file.path(benchmark_root, "FULL_SIZE_SWEEP_BENCHMARK.csv"))
  lines_new(c("status=FULL_SIZE_SWEEP_BENCHMARK_COMPLETE",
              "main_fits_started=FALSE", "truth_unseal_authorized=FALSE",
              paste0("marginal_T1_seconds=", marginal)),
            file.path(benchmark_root, "BENCHMARK_COMPLETE.txt"))
  print(result)
  invisible(result)
}

select_winners <- function(output_root, smoke = FALSE) {
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  assert(all(file.exists(paths)),
         "All registered fit terminal records are required for selection.")
  records <- lapply(paths, readRDS)
  endpoints <- do.call(rbind, lapply(records, function(x) {
    data.frame(
      fit_id = x$fit_id, data_id = x$data_id,
      logical_data_id = x$logical_data_id, beta_arm = x$beta_arm,
      start_index = x$start_index, fit_seed = x$fit_seed,
      status = x$status, objective_eligible = x$objective_eligible,
      invalid_reasons = x$invalid_reasons,
      final_ordinary_elbo = x$final_ordinary_elbo,
      annealing_sweeps = x$annealing_sweeps,
      ordinary_t1_sweeps = x$ordinary_t1_sweeps,
      strict_practical_converged = x$strict_practical_converged,
      elapsed_seconds = x$elapsed_seconds,
      truth_used_for_fit_or_selection = x$truth_used_for_fit_or_selection,
      stringsAsFactors = FALSE)
  }))
  assert(identical(as.character(endpoints$fit_id),
                   as.character(manifest$fit_id)) &&
           !any(endpoints$truth_used_for_fit_or_selection),
         "Endpoint join or truth isolation failed.")
  winners <- list()
  status <- list()
  for (data_id in unique(manifest$data_id)) {
    group <- endpoints[endpoints$data_id == data_id, , drop = FALSE]
    eligible <- group$status == "COMPLETE" & group$objective_eligible &
      is.finite(group$final_ordinary_elbo)
    status[[length(status) + 1L]] <- data.frame(
      data_id = data_id, n_registered = nrow(group),
      n_eligible = sum(eligible),
      winner_status = if (any(eligible)) "SELECTED" else "NO_ELIGIBLE_START",
      stringsAsFactors = FALSE)
    if (any(eligible)) {
      candidate <- group[eligible, , drop = FALSE]
      candidate <- candidate[order(-candidate$final_ordinary_elbo,
                                   candidate$fit_id), , drop = FALSE]
      winners[[length(winners) + 1L]] <- candidate[1L, , drop = FALSE]
    }
  }
  selected <- if (length(winners)) do.call(rbind, winners) else
    endpoints[FALSE, , drop = FALSE]
  root <- file.path(output_root, "truth_free_selection")
  assert(!dir.exists(root), "Refusing to overwrite frozen selection.")
  dir.create(root, recursive = TRUE, mode = "0700")
  csv_new(endpoints, file.path(root, "ALL_TRUTH_FREE_ENDPOINTS.csv"))
  csv_new(selected, file.path(root, "TRUTH_FREE_WINNERS.csv"))
  csv_new(do.call(rbind, status),
          file.path(root, "SELECTION_STATUS.csv"))
  lines_new(c("status=TRUTH_FREE_SELECTION_FROZEN",
              "rule=maximum_eligible_ordinary_T1_ELBO_within_data_arm",
              "tie_break=fit_id_lexicographic",
              "cross_arm_ELBO_comparison=FALSE",
              "truth_used_for_selection=FALSE",
              "truth_unseal_authorized=FALSE",
              paste0("registered_fits=", nrow(endpoints)),
              paste0("selected_winners=", nrow(selected)),
              paste0("smoke=", smoke)),
            file.path(root, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  cat("SELECTION_FROZEN endpoints=", nrow(endpoints),
      " winners=", nrow(selected), "\n", sep = "")
  invisible(selected)
}

batch <- function(output_root, workers) {
  assert(.Platform$OS.type == "unix", "Batch concurrency requires Linux.")
  assert(file.exists(file.path(output_root, "preflight",
                               "BENCHMARK_COMPLETE.txt")),
         "Full-size sweep benchmark must finish before batch fitting.")
  assert(is.finite(workers) && workers >= 1L && workers <= 8L,
         "Outer worker count must be from 1 through 8.")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  assert(nrow(manifest) == 36L && !anyDuplicated(manifest$fit_id),
         "Main fit manifest must contain exactly 36 unique fits.")
  logs <- file.path(output_root, "supervisor_logs")
  dir.create(logs, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  run_one <- function(fit_id) {
    status <- system2(
      "Rscript",
      c(script_path, "--action=fit", paste0("--fit-id=", fit_id),
        paste0("--output-root=", output_root)),
      stdout = file.path(logs, paste0(fit_id, ".stdout.log")),
      stderr = file.path(logs, paste0(fit_id, ".stderr.log")))
    data.frame(fit_id = fit_id, exit_status = as.integer(status))
  }
  statuses <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = workers,
    mc.preschedule = FALSE, mc.set.seed = FALSE)
  csv_new(do.call(rbind, statuses),
          file.path(output_root, "BATCH_PROCESS_STATUS.csv"))
  complete <- length(list.files(file.path(output_root, "fits"),
                                "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  errors <- length(list.files(file.path(output_root, "fits"),
                              "^FIT_ERROR[.]txt$", recursive = TRUE))
  if (complete + errors == 36L) {
    select_winners(output_root, FALSE)
    lines_new(c("status=MAIN_FITS_TERMINAL_AND_SELECTION_FROZEN",
                paste0("complete=", complete),
                paste0("errors=", errors),
                paste0("outer_workers=", workers),
                "per_fit_n_cpus=1", "truth_unseal_authorized=FALSE",
                "scientific_evaluation_started=FALSE"),
              file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  } else {
    lines_new(c("status=BATCH_INCOMPLETE_NO_SELECTIVE_RERUN",
                paste0("complete=", complete), paste0("errors=", errors)),
              file.path(output_root, "BATCH_INCOMPLETE.txt"))
    stop("Batch ended without all 36 terminal records.", call. = FALSE)
  }
  invisible(TRUE)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  load_package()
  d <- data_manifest(FALSE)
  f <- fit_manifest(d, FALSE)
  registered <- read_manifest(file.path(
    script_root, "FIT_MANIFEST_CANDIDATE_V1.csv"))
  assert(nrow(d) == 3L && nrow(f) == 36L &&
           all(table(f$data_id) == 12L) &&
           all(table(f$fit_seed) == 3L) &&
           identical(as.character(d$beta_arm), c("C0", "C1", "C2")) &&
           all(d$truth_available_to_fit == FALSE) &&
           all(d$truth_available_to_selection == FALSE) &&
           all(d$truth_unseal_authorized == FALSE),
         "Covariate Route-B pilot manifest or seed pairing is invalid.")
  assert(nrow(registered) == 36L &&
           identical(as.character(f$fit_id),
                     as.character(registered$fit_id)) &&
           identical(as.character(f$logical_data_id),
                     as.character(registered$logical_data_id)) &&
           identical(as.character(f$beta_arm),
                     as.character(registered$arm)) &&
           identical(as.integer(f$start_index),
                     as.integer(registered$start_index)) &&
           identical(as.integer(f$fit_seed),
                     as.integer(registered$fit_seed)) &&
           identical(as.integer(f$n_cpus),
                     as.integer(registered$n_cpus)),
         "Runtime fit IDs or seeds do not exactly match the candidate manifest.")
  cat("CHECK_PASS data=3 fits=36 package=",
      as.character(utils::packageVersion("multiFSYNC")), "\n", sep = "")
} else if (action == "generate") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, FALSE)
} else if (action == "benchmark") {
  assert(nzchar(output_root), "--output-root is required.")
  benchmark(output_root)
} else if (action == "fit") {
  assert(nzchar(output_root) && nzchar(cli$`fit-id` %||% ""),
         "--output-root and --fit-id are required.")
  fit_one(output_root, cli$`fit-id`, identical(cli$smoke, "true"))
} else if (action == "select") {
  assert(nzchar(output_root), "--output-root is required.")
  select_winners(output_root, identical(cli$smoke, "true"))
} else if (action == "batch") {
  assert(nzchar(output_root), "--output-root is required.")
  batch(output_root,
        as.integer(Sys.getenv("COVARIATE_ROUTEB_OUTER_WORKERS",
                              unset = "8")))
} else if (action == "smoke") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, TRUE)
  rows <- fit_manifest(data_manifest(TRUE), TRUE)
  for (fit_id in rows$fit_id) fit_one(output_root, fit_id, TRUE)
  select_winners(output_root, TRUE)
  lines_new(c("status=MICRO_SMOKE_PASS", "generated_arms=3", "fits=3",
              "truth_unseal_authorized=FALSE"),
            file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
  cat("MICRO_SMOKE_PASS\n")
} else if (action == "status") {
  assert(nzchar(output_root), "--output-root is required.")
  cat("complete=", length(list.files(file.path(output_root, "fits"),
                                       "^FIT_COMPLETE[.]txt$", recursive = TRUE)),
      " errors=", length(list.files(file.path(output_root, "fits"),
                                     "^FIT_ERROR[.]txt$", recursive = TRUE)),
      " selected=", file.exists(file.path(output_root,
                                           "truth_free_selection",
                                           "TRUTH_FREE_SELECTION_FROZEN.txt")),
      "\n", sep = "")
} else stop("Unknown action: ", action, call. = FALSE)

