#!/usr/bin/env Rscript

# Small paired bridge: historical analytic periodic mean versus independent
# periodic B-spline mean truth. Truth is sealed through winner selection.
options(warn = 1, stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")

`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
assert(length(script_arg) == 1L, "Run this file through Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
script_root <- dirname(script_path)
repo_root <- normalizePath(file.path(script_root, "..", "..", ".."),
                           winslash = "/", mustWork = TRUE)
experiment_id <- "BSPLINE_MEAN_PAIRED_BRIDGE_V1_20260926"
expected_generator_sha256 <-
  "226131949adf52a28fc5e8e91b765a8cdb98a4e2c3a2ffbfa7ee40fc7bc932dd"
expected_stage7a_sha256 <-
  "7b99b9030433d7033eec74cd3fa12736a4952b22a315d7b358ce881c9d5e3732"

assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
sha256 <- function(path) tolower(digest::digest(
  file = path, algo = "sha256", serialize = FALSE))
generator_path <- file.path(repo_root, "R", "generate_data_structured.R")
stage7a_path <- file.path(
  repo_root, "experiments", "v0L-V", "covariate-effect-stage7a-20260919",
  "run_covariate_stage7a_20260919_v1.R")
assert(file.exists(generator_path) && file.exists(stage7a_path),
       "Required package or Stage-7A source is missing.")
assert(identical(sha256(generator_path), expected_generator_sha256),
       "Generator SHA256 differs from the registered bridge source.")
assert(identical(sha256(stage7a_path), expected_stage7a_sha256),
       "Stage-7A runner SHA256 differs from the registered source.")

# Reuse tested I/O, objective and calibrated-loading declarations without
# executing the Stage-7A CLI.
stage7a <- new.env(parent = globalenv())
expressions <- parse(stage7a_path)
cli_index <- which(vapply(expressions, function(x) {
  is.call(x) && length(x) >= 3L && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name("cli"))
}, logical(1L)))
assert(length(cli_index) == 1L && cli_index > 1L,
       "Cannot locate the Stage-7A dispatch boundary.")
for (i in seq_len(cli_index - 1L)) eval(expressions[[i]], envir = stage7a)
stage7a$script_path <- normalizePath(stage7a_path, winslash = "/")
stage7a$script_root <- dirname(stage7a$script_path)

save_new <- stage7a$save_new
csv_new <- stage7a$csv_new
lines_new <- stage7a$lines_new
read_manifest <- stage7a$read_manifest
parse_cli <- stage7a$parse_cli

load_package <- function() {
  stage7a$load_package()
  assert(as.character(utils::packageVersion("multiFSYNC")) == "0.3.0.9000",
         "The installed multiFSYNC version must be 0.3.0.9000.")
  invisible(TRUE)
}

data_manifest <- function(smoke = FALSE) {
  index <- if (smoke) 1L else 1:2
  arm <- c("periodic_phase", "independent_splinet")
  rows <- do.call(rbind, lapply(index, function(k) data.frame(
    experiment_id = experiment_id,
    logical_data_id = sprintf("meanb_%02d", k),
    data_id = paste0(sprintf("meanb_%02d", k), "__", arm),
    data_index = k, mean_arm = arm,
    data_seed = 99311000L + k,
    truth_basis_seed = 99312000L + k,
    truth_loading_seed = 99313000L + k,
    mean_truth_seed = 99314000L + k,
    truth_basis_method = "iid_common_splinet_stiefel",
    truth_basis_dim = 8L, truth_basis_smoothness_alpha = 3,
    mean_basis_dim = ifelse(arm == "independent_splinet", 8L, NA_integer_),
    mean_smoothness_alpha = ifelse(
      arm == "independent_splinet", 1.519826404906808, NA_real_),
    loading_scenario = "controlled_identifiable_fixed90_calibrated",
    truth_available_to_fit = FALSE,
    truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    truth_unseal_authorized = FALSE,
    development_only = TRUE, formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE)))
  rownames(rows) <- NULL
  rows
}

fit_manifest <- function(data_rows, smoke = FALSE) {
  starts <- if (smoke) 1L else 1:12
  fit_K <- if (smoke) 3L else 9L
  rows <- do.call(rbind, lapply(seq_len(nrow(data_rows)), function(i) {
    row <- data_rows[i, , drop = FALSE]
    data.frame(
      experiment_id = experiment_id,
      fit_id = sprintf("%s__start%02d", row$data_id[[1L]], starts),
      data_id = row$data_id[[1L]],
      logical_data_id = row$logical_data_id[[1L]],
      data_index = row$data_index[[1L]], mean_arm = row$mean_arm[[1L]],
      fit_K = fit_K, fit_config_id = sprintf("K%02d", fit_K),
      selection_stratum_id = row$data_id[[1L]],
      paired_start_id = sprintf("%s__start%02d",
                                row$logical_data_id[[1L]], starts),
      start_index = starts, fit_seed = 99315000L + starts,
      method_id = "G12", initialization = "random",
      function_initialization = "gram_unit_energy", pre_score_sweeps = 1L,
      anneal = "1;1.9;100", planned_annealing_sweeps = 99L,
      planned_ordinary_t1_sweeps = 400L, n_cpus = 1L,
      objective_eligibility_required = TRUE,
      truth_free_winner_candidate = TRUE,
      cross_mean_arm_elbo_comparison_forbidden = TRUE,
      actual_racing = FALSE, continuation = FALSE, automatic_800 = FALSE,
      truth_available_to_fit = FALSE,
      truth_available_to_stopping = FALSE,
      truth_available_to_selection = FALSE,
      development_only = TRUE, formal_v0lv_result = FALSE,
      formal_paper_mc_result = FALSE, stringsAsFactors = FALSE)
  }))
  rownames(rows) <- NULL
  rows
}

validate_manifests <- function(data_rows, fit_rows, smoke = FALSE) {
  expected_logical <- if (smoke) 1L else 2L
  expected_starts <- if (smoke) 1L else 12L
  expected_data <- 2L * expected_logical
  expected_fits <- expected_data * expected_starts
  strata <- interaction(fit_rows$data_id, drop = TRUE)
  pairs <- interaction(
    fit_rows$logical_data_id, fit_rows$start_index, drop = TRUE)
  assert(nrow(data_rows) == expected_data && nrow(fit_rows) == expected_fits &&
           !anyDuplicated(data_rows$data_id) && !anyDuplicated(fit_rows$fit_id),
         "Bridge manifest dimensions or identifiers are invalid.")
  assert(all(table(data_rows$logical_data_id) == 2L) &&
           all(table(strata) == expected_starts) && all(table(pairs) == 2L),
         "Arm/start pairing in the manifests is invalid.")
  paired_seed_ok <- vapply(split(fit_rows, pairs), function(x) {
    nrow(x) == 2L && length(unique(x$fit_seed)) == 1L &&
      setequal(x$mean_arm, c("periodic_phase", "independent_splinet"))
  }, logical(1L))
  assert(all(paired_seed_ok) && all(fit_rows$n_cpus == 1L) &&
           all(fit_rows$method_id == "G12") &&
           all(fit_rows$planned_ordinary_t1_sweeps == 400L) &&
           !any(fit_rows$actual_racing) && !any(fit_rows$continuation) &&
           !any(fit_rows$automatic_800) &&
           !any(fit_rows$truth_available_to_fit) &&
           !any(fit_rows$truth_available_to_selection),
         "Route, pairing or truth-isolation contract failed.")
  invisible(TRUE)
}

generator_config <- function(smoke = FALSE) list(
  S = 2L, n_s = if (smoke) c(4L, 4L) else c(30L, 30L),
  p = if (smoke) 20L else 500L, d = 0L,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
  K = if (smoke) 3L else 9L,
  n_obs = if (smoke) 5L else 8L, common_grid = FALSE,
  n_obs_range = if (smoke) c(4L, 6L) else c(6L, 9L),
  sigma_eps = 0.3, bool_sparse_loadings = TRUE,
  prop_sparse = if (smoke) 0.8 else 0.9, score_var_decay = TRUE,
  n_dense = if (smoke) 101L else 201L, mean_amp = 0.6, beta_amp = 0,
  bs_degree = c(2L, 3L), identified_loadings = TRUE,
  sparsity_mode = "fixed", truth_basis_method = "iid_common_splinet_stiefel",
  truth_basis_dim = 8L)

generate_one <- function(row, config) {
  extra <- list(
    seed = row$data_seed[[1L]],
    truth_basis_seed = row$truth_basis_seed[[1L]],
    truth_loading_seed = row$truth_loading_seed[[1L]],
    mean_truth_method = row$mean_arm[[1L]])
  if (row$mean_arm[[1L]] == "independent_splinet") {
    extra$mean_truth_seed <- row$mean_truth_seed[[1L]]
    extra$mean_basis_dim <- 8L
    extra$mean_smoothness_alpha <- 1.519826404906808
  }
  data <- do.call(multiFSYNC::simulate_multi_study_structured,
                  c(config, extra))
  stage7a$calibrate_loadings(data)
}

mean_dense_truth <- function(data, arm) {
  truth <- data$true_params
  if (arm == "independent_splinet") return(truth$mean_dense_list)
  lapply(seq_len(truth$S), function(s) vapply(seq_len(truth$p), function(j) {
    multiFSYNC:::mean_function_periodic(
      truth$t_grid_dense, truth$phase_mu_mat[s, j], 0.6)
  }, numeric(length(truth$t_grid_dense))))
}

projection_summary <- function(data, arm) {
  truth <- data$true_params
  grid <- truth$t_grid_dense
  dt <- diff(grid)
  w <- c(dt[[1L]] / 2,
         (dt[-1L] + dt[-length(dt)]) / 2,
         dt[[length(dt)]] / 2)
  design <- cbind(
    1, grid,
    multiFSYNC:::ZOSull(
      grid, range.x = c(0, 1), intKnots = truth$int_knots))
  weighted_design <- design * sqrt(w)
  means <- mean_dense_truth(data, arm)
  rows <- lapply(seq_len(truth$S), function(s) {
    weighted_truth <- means[[s]] * sqrt(w)
    coefficient <- qr.solve(weighted_design, weighted_truth)
    residual <- means[[s]] - design %*% coefficient
    ise <- colSums(w * residual^2)
    energy <- colSums(w * means[[s]]^2)
    relative <- ise / pmax(energy, .Machine$double.eps)
    data.frame(
      study = s, mean_arm = arm, n_variables = truth$p,
      median_projection_relative_ise = stats::median(relative),
      q90_projection_relative_ise = as.numeric(stats::quantile(
        relative, 0.9, names = FALSE)),
      max_projection_relative_ise = max(relative),
      median_grid_energy = stats::median(energy),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

truth_qc <- function(data, row, config) {
  truth <- data$true_params
  meta <- c(truth$phi_basis_meta,
            unlist(truth$psi_basis_meta, recursive = FALSE))
  gram_error <- max(vapply(meta, function(x) {
    max(abs(x$gram - diag(nrow(x$gram))))
  }, numeric(1L)))
  loadings <- c(list(truth$a_true[, 1L]),
                lapply(truth$b_true, function(x) x[, 1L]))
  loading_cosine <- vapply(loadings[-1L], function(x) {
    abs(sum(loadings[[1L]] * x)) /
      sqrt(sum(loadings[[1L]]^2) * sum(x^2))
  }, numeric(1L))
  counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  reconstruction_error <- 0
  for (s in seq_len(truth$S)) for (i in seq_len(truth$n_s[[s]])) {
    observed <- do.call(cbind, data$Y[[s]][[i]])
    reconstruction_error <- max(reconstruction_error, abs(
      observed - truth$signal_true_values[[s]][[i]] -
        truth$noise_true_values[[s]][[i]]))
  }
  is_splinet <- row$mean_arm[[1L]] == "independent_splinet"
  mean_ok <- if (!is_splinet) {
    !"mean_truth_method" %in% names(truth)
  } else {
    energy <- unlist(lapply(
      truth$mean_basis_meta, `[[`, "continuous_energy"))
    integral <- unlist(lapply(
      truth$mean_basis_meta, `[[`, "continuous_integral"))
    identical(truth$mean_truth_method, "independent_splinet") &&
      identical(truth$mean_truth_seed,
                as.integer(row$mean_truth_seed[[1L]])) &&
      max(abs(energy - 0.225)) < 1e-10 &&
      max(abs(integral)) < 1e-10
  }
  checks <- c(
    dimensions = truth$S == 2L && all(truth$n_s == config$n_s) &&
      truth$p == config$p && truth$d_use == 0L,
    truth_counts = truth$L_f == 1L &&
      identical(as.integer(truth$L_s_by_study), c(1L, 1L)) &&
      identical(as.integer(truth$M_f), 2L) &&
      all(vapply(truth$M_s, function(x) identical(as.integer(x), 2L),
                 logical(1L))),
    observation_counts = all(counts >= config$n_obs_range[[1L]] &
                               counts <= config$n_obs_range[[2L]]),
    route_b_factor_truth = identical(
      truth$truth_basis_method, "iid_common_splinet_stiefel") &&
      identical(truth$truth_basis_dim, 8L),
    factor_continuous_orthonormality = gram_error < 1e-8,
    strong_loading_identification = max(loading_cosine) < 1e-10,
    registered_mean_truth = mean_ok,
    response_reconstruction = reconstruction_error < 1e-10)
  data.frame(
    check_id = names(checks), passed = unname(checks),
    detail = c(
      sprintf("S=%d;n=%s;p=%d;d=%d", truth$S,
              paste(truth$n_s, collapse = ";"), truth$p, truth$d_use),
      "L=1;1;1;M=2;2;2", paste(range(counts), collapse = ";"),
      "method=iid_common_splinet_stiefel;J=8;alpha=3",
      format(gram_error, scientific = TRUE),
      paste(format(loading_cosine, digits = 8), collapse = ";"),
      row$mean_arm[[1L]], format(reconstruction_error, scientific = TRUE)),
    stringsAsFactors = FALSE)
}

pair_qc <- function(periodic, splinet, logical_data_id) {
  left <- periodic$true_params
  right <- splinet$true_params
  top_fields <- c("time_obs", "C", "Z")
  truth_fields <- c(
    "phase_mu_mat", "phase_beta", "zeta_true", "xi_true", "eta_true",
    "chi_true", "zeta_fpca_true", "xi_fpca_true", "phi_dense_list",
    "psi_dense_list", "theta_dense_list", "kappa_dense_list",
    "lambda_phi_true", "lambda_psi_true", "a_true", "b_true",
    "gamma_a_true", "gamma_b_true", "omega_a_true", "omega_b_true",
    "f_true_values", "g_true_values", "beta_true_values",
    "noise_true_values", "n_obs_by_subject", "int_knots",
    "sigma2_eps_true", "loading_calibration")
  exact_top <- vapply(top_fields, function(x) {
    identical(periodic[[x]], splinet[[x]])
  }, logical(1L))
  exact_truth <- vapply(truth_fields, function(x) {
    identical(left[[x]], right[[x]])
  }, logical(1L))
  response_identity <- 0
  for (s in seq_len(left$S)) for (i in seq_len(left$n_s[[s]])) {
    y_left <- do.call(cbind, periodic$Y[[s]][[i]])
    y_right <- do.call(cbind, splinet$Y[[s]][[i]])
    response_identity <- max(response_identity, abs(
      (y_right - y_left) -
        (right$mu_true_values[[s]][[i]] - left$mu_true_values[[s]][[i]])))
  }
  checks <- c(
    exact_nonmean_top_fields = all(exact_top),
    exact_nonmean_truth_fields = all(exact_truth),
    response_difference_is_mean_difference = response_identity < 1e-10,
    mean_truth_differs = max(abs(
      right$mu_true_values[[1L]][[1L]] -
        left$mu_true_values[[1L]][[1L]])) > 1e-6)
  data.frame(
    logical_data_id = logical_data_id,
    check_id = names(checks), passed = unname(checks),
    detail = c(
      paste(names(exact_top)[!exact_top], collapse = ";"),
      paste(names(exact_truth)[!exact_truth], collapse = ";"),
      format(response_identity, scientific = TRUE), "mu differs"),
    stringsAsFactors = FALSE)
}

save_dataset <- function(data, row, config, qc, projection, output_root) {
  root <- file.path(output_root, "data", row$data_id[[1L]])
  assert(!dir.exists(root), paste0("Data directory exists: ", root))
  dir.create(file.path(root, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  observation <- list(
    bundle_class = "bspline_mean_bridge_observation_only",
    experiment_id = experiment_id, data_id = row$data_id[[1L]],
    logical_data_id = row$logical_data_id[[1L]],
    mean_arm = row$mean_arm[[1L]],
    Y = data$Y, Z = data$Z, time_obs = data$time_obs,
    dimensions = config[c("S", "n_s", "p", "d", "L_f", "L_s",
                          "M_f", "M_s")],
    generator_cache_K = config$K, truth_available_to_fit = FALSE)
  sealed <- list(
    bundle_class = "bspline_mean_bridge_sealed_truth",
    experiment_id = experiment_id, data_id = row$data_id[[1L]],
    logical_data_id = row$logical_data_id[[1L]],
    mean_arm = row$mean_arm[[1L]], data_seed = row$data_seed[[1L]],
    truth_basis_seed = row$truth_basis_seed[[1L]],
    truth_loading_seed = row$truth_loading_seed[[1L]],
    mean_truth_seed = row$mean_truth_seed[[1L]],
    truth_unseal_authorized = FALSE, data = data)
  save_new(observation, file.path(root, "observation_bundle.rds"))
  save_new(sealed, file.path(root, "sealed_truth", "truth_bundle.rds"))
  csv_new(qc, file.path(root, "sealed_truth", "TRUTH_GENERATION_QC.csv"))
  csv_new(projection, file.path(
    root, "sealed_truth", "MEAN_PROJECTION_FLOOR_SUMMARY.csv"))
  lines_new(c("status=DATA_SEAL_COMPLETE", "truth_unseal_authorized=FALSE",
              "scientific_evaluation_started=FALSE"),
            file.path(root, "DATA_SEAL_COMPLETE.txt"))
}

generate <- function(output_root, smoke = FALSE) {
  load_package()
  assert(!dir.exists(output_root),
         paste0("Refusing to overwrite output root: ", output_root))
  data_rows <- data_manifest(smoke)
  fit_rows <- fit_manifest(data_rows, smoke)
  validate_manifests(data_rows, fit_rows, smoke)
  dir.create(output_root, recursive = TRUE, mode = "0700")
  csv_new(data_rows, file.path(output_root, "DATA_MANIFEST.csv"))
  csv_new(fit_rows, file.path(output_root, "FIT_MANIFEST.csv"))
  csv_new(data.frame(
    relative_path = c(
      "R/generate_data_structured.R",
      file.path("experiments", "v0L-V", "covariate-effect-stage7a-20260919",
                basename(stage7a_path)),
      file.path("experiments", "v0L-V",
                "bspline-mean-paired-bridge-20260926", basename(script_path))),
    sha256 = c(sha256(generator_path), sha256(stage7a_path),
               sha256(script_path)), stringsAsFactors = FALSE),
    file.path(output_root, "SOURCE_PROVENANCE.csv"))

  config <- generator_config(smoke)
  pair_checks <- list()
  for (logical_id in unique(data_rows$logical_data_id)) {
    pair_rows <- data_rows[data_rows$logical_data_id == logical_id, , drop = FALSE]
    generated <- list()
    for (index in seq_len(nrow(pair_rows))) {
      row <- pair_rows[index, , drop = FALSE]
      data <- generate_one(row, config)
      qc <- truth_qc(data, row, config)
      assert(all(qc$passed), paste0(
        "Truth generation QC failed for ", row$data_id[[1L]], ": ",
        paste(qc$check_id[!qc$passed], collapse = ";")))
      projection <- projection_summary(data, row$mean_arm[[1L]])
      save_dataset(data, row, config, qc, projection, output_root)
      generated[[row$mean_arm[[1L]]]] <- data
    }
    one_pair <- pair_qc(
      generated$periodic_phase, generated$independent_splinet, logical_id)
    assert(all(one_pair$passed), paste0(
      "Paired generation QC failed for ", logical_id, ": ",
      paste(one_pair$check_id[!one_pair$passed], collapse = ";")))
    pair_checks[[length(pair_checks) + 1L]] <- one_pair
  }
  csv_new(do.call(rbind, pair_checks),
          file.path(output_root, "PAIRED_GENERATION_QC.csv"))
  lines_new(c("status=DATA_GENERATION_AND_SEAL_COMPLETE",
              paste0("datasets=", nrow(data_rows)),
              paste0("planned_fits=", nrow(fit_rows)),
              "truth_unseal_authorized=FALSE",
              "scientific_evaluation_started=FALSE"),
            file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("MEAN_BRIDGE_GENERATION_COMPLETE datasets=", nrow(data_rows),
      " fits=", nrow(fit_rows), " smoke=", smoke, "\n", sep = "")
}

capture_fit <- function(observation, row, mode = "main") {
  smoke <- identical(mode, "smoke")
  benchmark <- startsWith(mode, "benchmark")
  if (smoke) {
    anneal <- c(1, 1.2, 3); maxit <- 7L; n_g <- 31L
    convergence_rule <- "parameters"; practical_control <- NULL
    trace_sweeps <- 1L
  } else if (benchmark) {
    anneal <- NULL; maxit <- as.integer(sub("^benchmark", "", mode))
    n_g <- 51L; convergence_rule <- "parameters"; practical_control <- NULL
    trace_sweeps <- 1L
  } else {
    anneal <- c(1, 1.9, 100); maxit <- 499L; n_g <- 51L
    convergence_rule <- "practical"
    practical_control <- multiFSYNC::g12_stopping_control()
    practical_control$checkpoints <- c(380L, 400L)
    trace_sweeps <- c(1L, 380L, 400L)
  }
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(
    multiFSYNC::bayesSYNC_multi_pre_score(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
      K = as.integer(row$fit_K[[1L]]), n_g = n_g,
      anneal = anneal, maxit = maxit, n_cpus = 1L, verbose = FALSE,
      seed = as.integer(row$fit_seed[[1L]]),
      tol_abs = if (smoke || benchmark) 0 else 1e-3,
      tol_rel = if (smoke || benchmark) 0 else 1e-5,
      bool_scale = FALSE, bool_var_spec_prob = FALSE,
      d_0 = as.integer(observation$dimensions$p),
      convergence_rule = convergence_rule, lambda_orth = 0,
      practical_control = practical_control, initialization = "random",
      initialization_control = list(
        grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8),
      continuation_state = NULL, pre_score_sweeps = 1L,
      trace_sweeps = trace_sweeps,
      function_initialization = "gram_unit_energy"),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  error <- if (inherits(value, "error")) conditionMessage(value) else NULL
  fit <- if (is.null(error)) value else NULL
  if (is.null(error) && !smoke && !benchmark) {
    valid <- tryCatch({
      assert(fit$annealing_sweeps == 99L && fit$t1_sweeps == 400L &&
               fit$n_g == 51L, "Incomplete fixed-400 fit endpoint.")
      assert(fit$n_cpus_used == 1L &&
               fit$pre_score_interface$pre_score_sweeps == 1L &&
               identical(fit$pre_score_interface$function_initialization,
                         "gram_unit_energy") &&
               identical(as.integer(fit$K), as.integer(row$fit_K[[1L]])),
             "Returned route or K differs from the manifest.")
      TRUE
    }, error = function(e) e)
    if (inherits(valid, "error")) {
      error <- conditionMessage(valid); fit <- NULL
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
  assert(identical(observation$bundle_class,
                   "bspline_mean_bridge_observation_only") &&
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
  } else stage7a$objective_status(fit)
  final_elbo <- if (is.null(fit) || !length(fit$ELBO)) NA_real_ else
    as.numeric(tail(fit$ELBO, 1L))
  terminal <- list(
    experiment_id = experiment_id, fit_id = fit_id,
    data_id = row$data_id[[1L]], logical_data_id = row$logical_data_id[[1L]],
    data_index = row$data_index[[1L]], mean_arm = row$mean_arm[[1L]],
    fit_K = row$fit_K[[1L]], fit_config_id = row$fit_config_id[[1L]],
    selection_stratum_id = row$selection_stratum_id[[1L]],
    paired_start_id = row$paired_start_id[[1L]],
    start_index = row$start_index[[1L]], fit_seed = row$fit_seed[[1L]],
    n_cpus = 1L, status = if (is.null(fit)) "ERROR" else "COMPLETE",
    error = outcome$error, warnings = outcome$warnings,
    elapsed_seconds = outcome$elapsed_seconds,
    objective_eligible = objective$eligible,
    invalid_reasons = paste(objective$invalid_reasons, collapse = ";"),
    final_ordinary_elbo = final_elbo,
    annealing_sweeps = if (is.null(fit)) NA_integer_ else fit$annealing_sweeps,
    ordinary_t1_sweeps = if (is.null(fit)) NA_integer_ else fit$t1_sweeps,
    strict_practical_converged = if (is.null(fit)) FALSE else
      isTRUE(fit$practical_converged) && isTRUE(fit$converged) &&
      identical(as.character(fit$convergence_status), "converged_practical"),
    truth_used_for_fit_or_selection = FALSE,
    continuation_used = FALSE, formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE)
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
  if (is.data.frame(fit$practical_diagnostics)) {
    csv_new(fit$practical_diagnostics,
            file.path(fit_dir, "PRACTICAL_DIAGNOSTICS.csv"))
  }
  lines_new(c("status=COMPLETE", paste0("fit_id=", fit_id),
              paste0("mean_arm=", row$mean_arm[[1L]]),
              paste0("objective_eligible=", objective$eligible),
              paste0("elapsed_seconds=", outcome$elapsed_seconds),
              "n_cpus=1", "truth_used_for_fit_or_selection=FALSE"),
            file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id,
      " elapsed_seconds=", outcome$elapsed_seconds,
      " eligible=", objective$eligible, "\n", sep = "")
  invisible(terminal)
}

select_winners <- function(output_root, smoke = FALSE) {
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  paths <- file.path(output_root, "fits", manifest$fit_id,
                     "terminal_record.rds")
  assert(all(file.exists(paths)),
         "All registered fit terminal records are required for selection.")
  records <- lapply(paths, readRDS)
  endpoints <- do.call(rbind, lapply(records, function(x) data.frame(
    fit_id = x$fit_id, data_id = x$data_id,
    logical_data_id = x$logical_data_id, data_index = x$data_index,
    mean_arm = x$mean_arm, fit_K = x$fit_K,
    selection_stratum_id = x$selection_stratum_id,
    paired_start_id = x$paired_start_id, start_index = x$start_index,
    fit_seed = x$fit_seed, status = x$status,
    objective_eligible = x$objective_eligible,
    invalid_reasons = x$invalid_reasons,
    final_ordinary_elbo = x$final_ordinary_elbo,
    annealing_sweeps = x$annealing_sweeps,
    ordinary_t1_sweeps = x$ordinary_t1_sweeps,
    strict_practical_converged = x$strict_practical_converged,
    elapsed_seconds = x$elapsed_seconds,
    truth_used_for_fit_or_selection = x$truth_used_for_fit_or_selection,
    stringsAsFactors = FALSE)))
  assert(identical(as.character(endpoints$fit_id),
                   as.character(manifest$fit_id)) &&
           !any(endpoints$truth_used_for_fit_or_selection),
         "Endpoint join or truth isolation failed.")
  winners <- list(); selection_status <- list()
  for (stratum in unique(manifest$selection_stratum_id)) {
    group <- endpoints[endpoints$selection_stratum_id == stratum, , drop = FALSE]
    eligible <- group$status == "COMPLETE" & group$objective_eligible &
      is.finite(group$final_ordinary_elbo)
    selection_status[[length(selection_status) + 1L]] <- data.frame(
      selection_stratum_id = stratum, data_id = group$data_id[[1L]],
      logical_data_id = group$logical_data_id[[1L]],
      mean_arm = group$mean_arm[[1L]], n_registered = nrow(group),
      n_eligible = sum(eligible),
      winner_status = if (any(eligible)) "SELECTED" else "NO_ELIGIBLE_START")
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
  csv_new(do.call(rbind, selection_status),
          file.path(root, "SELECTION_STATUS.csv"))
  lines_new(c("status=TRUTH_FREE_SELECTION_FROZEN",
              "rule=maximum_eligible_ordinary_T1_ELBO_within_data_and_mean_arm",
              "tie_break=fit_id_lexicographic",
              "cross_mean_arm_ELBO_comparison=FALSE",
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

benchmark <- function(output_root) {
  load_package()
  assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
         "Generate and seal before benchmarking.")
  root <- file.path(output_root, "preflight")
  assert(!dir.exists(root), "Refusing benchmark rerun.")
  dir.create(root, recursive = TRUE, mode = "0700")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  row <- manifest[manifest$data_index == 1L &
                    manifest$mean_arm == "independent_splinet" &
                    manifest$start_index == 1L, , drop = FALSE]
  assert(nrow(row) == 1L, "B-spline benchmark row is missing.")
  observation <- readRDS(file.path(output_root, "data", row$data_id[[1L]],
                                   "observation_bundle.rds"))
  short <- capture_fit(observation, row, "benchmark3")
  long <- capture_fit(observation, row, "benchmark8")
  assert(is.null(short$error) && is.null(long$error), paste0(
    "Full-size benchmark failed: ", short$error %||% long$error))
  assert(short$fit$t1_sweeps == 3L && long$fit$t1_sweeps == 8L &&
           isTRUE(all.equal(short$fit$ELBO, long$fit$ELBO[1:3],
                            tolerance = 1e-10)),
         "Benchmark trajectory is not reproducible.")
  marginal <- (long$elapsed_seconds - short$elapsed_seconds) / 5
  result <- data.frame(
    data_id = row$data_id[[1L]], fit_K = row$fit_K[[1L]],
    short_t1_sweeps = 3L, long_t1_sweeps = 8L,
    short_elapsed_seconds = short$elapsed_seconds,
    long_elapsed_seconds = long$elapsed_seconds,
    estimated_marginal_seconds_per_T1_sweep = marginal,
    projected_hours_48_fits_8_workers = 48 * 499 * marginal / 8 / 3600)
  csv_new(result, file.path(root, "FULL_SIZE_K9_SWEEP_BENCHMARK.csv"))
  lines_new(c("status=FULL_SIZE_K9_SWEEP_BENCHMARK_COMPLETE",
              "main_fits_started=FALSE", "truth_unseal_authorized=FALSE",
              paste0("marginal_T1_seconds=", marginal)),
            file.path(root, "BENCHMARK_COMPLETE.txt"))
  print(result)
}

batch <- function(output_root, workers) {
  assert(.Platform$OS.type == "unix", "Batch concurrency requires Linux.")
  assert(file.exists(file.path(output_root, "preflight",
                               "BENCHMARK_COMPLETE.txt")),
         "Full-size benchmark must finish before the batch.")
  assert(is.finite(workers) && workers >= 1L && workers <= 8L,
         "Outer worker count must be from 1 through 8.")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  validate_manifests(read_manifest(file.path(output_root, "DATA_MANIFEST.csv")),
                     manifest, FALSE)
  logs <- file.path(output_root, "supervisor_logs")
  assert(!dir.exists(logs), "Refusing to rerun the batch.")
  dir.create(logs, recursive = TRUE, mode = "0700")
  run_one <- function(fit_id) {
    status <- system2(
      "Rscript", c(script_path, "--action=fit", paste0("--fit-id=", fit_id),
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
  if (complete + errors == nrow(manifest)) {
    select_winners(output_root, FALSE)
    lines_new(c("status=FITS_TERMINAL_AND_TRUTH_FREE_SELECTION_FROZEN",
                paste0("complete=", complete), paste0("errors=", errors),
                paste0("outer_workers=", workers), "per_fit_n_cpus=1",
                "truth_unseal_authorized=FALSE",
                "scientific_evaluation_started=FALSE"),
              file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  } else {
    lines_new(c("status=BATCH_INCOMPLETE_NO_SELECTIVE_RERUN",
                paste0("complete=", complete), paste0("errors=", errors)),
              file.path(output_root, "BATCH_INCOMPLETE.txt"))
    stop("Batch ended without all registered terminal records.", call. = FALSE)
  }
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  load_package()
  data_rows <- data_manifest(FALSE); fit_rows <- fit_manifest(data_rows, FALSE)
  validate_manifests(data_rows, fit_rows, FALSE)
  cat("MEAN_BRIDGE_CHECK_PASS logical_data=2 arms=2 fits=48 K=9 starts=12 package=",
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
  batch(output_root, as.integer(cli$workers %||% "8"))
} else if (action == "smoke") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, TRUE)
  rows <- fit_manifest(data_manifest(TRUE), TRUE)
  for (fit_id in rows$fit_id) fit_one(output_root, fit_id, TRUE)
  select_winners(output_root, TRUE)
  lines_new(c("status=MICRO_SMOKE_PASS", "logical_data=1", "arms=2",
              "fits=2", "truth_unseal_authorized=FALSE"),
            file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
  cat("MEAN_BRIDGE_MICRO_SMOKE_PASS\n")
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
