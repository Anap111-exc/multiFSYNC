options(warn = 1, stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

pp_abort <- function(...) stop(paste0(...), call. = FALSE)
pp_assert <- function(value, message) if (!isTRUE(value)) pp_abort(message)
`%||%` <- function(left, right) if (is.null(left)) right else left
pp_iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}
pp_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
pp_parse_cli <- function(arguments) {
  result <- list()
  for (argument in arguments) {
    pp_assert(startsWith(argument, "--") && grepl("=", argument, fixed = TRUE),
              "Arguments must use --key=value syntax.")
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    pp_assert(nzchar(key) && !key %in% names(result),
              "Duplicate or empty CLI key.")
    result[[key]] <- paste(pieces[-1L], collapse = "=")
  }
  result
}
pp_bool <- function(value, default = FALSE) {
  if (is.null(value)) return(default)
  value <- tolower(as.character(value))
  pp_assert(value %in% c("true", "false"), "Boolean values must be true/false.")
  identical(value, "true")
}

pp_script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
pp_assert(length(pp_script_argument) == 1L, "Run with Rscript.")
PP_SCRIPT_FILE <- sub("^--file=", "", pp_script_argument[[1L]])
PP_ROOT <- dirname(PP_SCRIPT_FILE)
if (!dir.exists(PP_ROOT)) {
  PP_ROOT <- dirname(normalizePath(PP_SCRIPT_FILE, winslash = "/",
                                  mustWork = TRUE))
}
PP_REPO_ROOT <- file.path(PP_ROOT, "..", "..", "..")
if (!dir.exists(file.path(PP_REPO_ROOT, "R"))) {
  PP_REPO_ROOT <- normalizePath(PP_REPO_ROOT, winslash = "/", mustWork = TRUE)
}
PP_SNAPSHOT <- file.path(PP_REPO_ROOT, "experiments", "v0L-V", "rc1_snapshot")
PP_EXPERIMENT_ID <- "G12_PAPER_PIPELINE_PILOT_V1_20260829"
PP_MODEL_PARENT_COMMIT <- "6b10c796c0ad2d53435c9a134d9476313bf94586"
PP_DATA_IDS <- c("g12pp_base_01", "g12pp_weak_01")
PP_EXPECTED_FITS <- 26L
PP_WORKERS <- 3L

pp_read_csv <- function(path) {
  pp_assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
pp_sha256 <- function(path) {
  pp_assert(file.exists(path), paste0("Cannot hash missing file: ", path))
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
pp_find_one <- function(root, pattern) {
  hits <- list.files(root, pattern = pattern, recursive = TRUE,
                     full.names = TRUE)
  pp_assert(length(hits) == 1L,
            paste0("Expected one file matching ", pattern,
                   "; found ", length(hits), "."))
  hits[[1L]]
}
pp_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}
pp_write_lines <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
  invisible(path)
}
pp_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3)
  invisible(path)
}

pp_load_package <- function() {
  library_path <- pp_env("PAPER_PILOT_R_LIB")
  if (nzchar(library_path)) {
    pp_assert(dir.exists(library_path),
              paste0("PAPER_PILOT_R_LIB is missing: ", library_path))
    .libPaths(unique(c(library_path, .libPaths())))
  }
  if (pp_bool(pp_env("PAPER_PILOT_USE_LOAD_ALL", "false"))) {
    pp_assert(requireNamespace("pkgload", quietly = TRUE),
              "pkgload is required for local load_all smoke.")
    pkgload::load_all(PP_REPO_ROOT, quiet = TRUE, export_all = FALSE,
                      helpers = FALSE)
  }
  pp_assert(requireNamespace("multiFSYNC", quietly = TRUE),
            "The current multiFSYNC package is unavailable.")
  pp_assert(as.character(utils::packageVersion("multiFSYNC")) == "0.3.0.9000",
            "The installed multiFSYNC version must be 0.3.0.9000.")
  for (name in c("simulate_multi_study_structured",
                 "bayesSYNC_multi_pre_score", "g12_stopping_control")) {
    pp_assert(is.function(getExportedValue("multiFSYNC", name)),
              paste0("Missing public multiFSYNC function: ", name))
  }
  invisible(TRUE)
}

pp_source_runtime <- function(include_evaluation = FALSE) {
  pp_assert(dir.exists(PP_SNAPSHOT), "The rc1 snapshot is missing.")
  source(pp_find_one(PP_SNAPSHOT, "^v0lv_candidate2_runtime[.]R$"),
         local = FALSE)
  source(pp_find_one(PP_SNAPSHOT, "^r_route_v3_candidate2_runner[.]R$"),
         local = FALSE)
  source(pp_find_one(PP_SNAPSHOT, "^v0lv_candidate2_controller[.]R$"),
         local = FALSE)
  if (include_evaluation) {
    source(pp_find_one(PP_SNAPSHOT, "^v0lv_candidate2_evaluation[.]R$"),
           local = FALSE)
  }
  invisible(TRUE)
}

pp_verify_source_binding <- function() {
  path <- file.path(PP_ROOT, "SOURCE_BINDING.csv")
  pp_assert(file.exists(path), "SOURCE_BINDING.csv is missing.")
  binding <- pp_read_csv(path)
  pp_assert(all(c("relative_path", "sha256") %in% names(binding)) &&
              nrow(binding) >= 10L && !anyDuplicated(binding$relative_path),
            "SOURCE_BINDING.csv is malformed.")
  matches <- vapply(seq_len(nrow(binding)), function(index) {
    file <- file.path(PP_REPO_ROOT, binding$relative_path[[index]])
    file.exists(file) && identical(pp_sha256(file), binding$sha256[[index]])
  }, logical(1L))
  pp_assert(all(matches), paste0(
    "Source binding mismatch: ",
    paste(binding$relative_path[!matches], collapse = ";")
  ))
  invisible(binding)
}

pp_validate_data_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "data_id", "data_index", "data_seed", "scenario_id",
    "n_obs_min", "n_obs_max", "target_shared_specific_abs_cosine",
    "new_unseen_at_registration", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "truth_unseal_authorized", "engineering_pilot", "formal_v0lv_result",
    "formal_paper_mc_result"
  )
  pp_assert(is.data.frame(manifest) && nrow(manifest) == 2L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              identical(as.character(manifest$data_id), PP_DATA_IDS) &&
              identical(as.integer(manifest$data_seed), 82929001:82929002) &&
              !anyDuplicated(manifest$data_id) &&
              !anyDuplicated(manifest$data_seed) &&
              all(manifest$new_unseen_at_registration) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              !any(manifest$truth_unseal_authorized) &&
              all(manifest$engineering_pilot) &&
              !any(manifest$formal_v0lv_result) &&
              !any(manifest$formal_paper_mc_result),
            "Pilot data manifest contract failed.")
  pp_assert(manifest$scenario_id[[1L]] == "baseline_strong" &&
              manifest$scenario_id[[2L]] == "weak_loading_separation" &&
              all(manifest$n_obs_min == 6L) && all(manifest$n_obs_max == 9L) &&
              identical(as.numeric(manifest$target_shared_specific_abs_cosine),
                        c(0, 0.6)),
            "Pilot scenario registration failed.")
  invisible(TRUE)
}

pp_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "data_index", "data_seed",
    "scenario_id", "method_id", "route_id", "seed_index", "fit_seed",
    "n_cpus", "planned_annealing_sweeps", "planned_ordinary_t1_sweeps",
    "objective_eligibility_required", "truth_free_winner_candidate",
    "cross_model_elbo_comparison_forbidden", "actual_racing",
    "continuation", "automatic_800", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "engineering_pilot", "formal_v0lv_result", "formal_paper_mc_result"
  )
  pp_assert(is.data.frame(manifest) && nrow(manifest) == PP_EXPECTED_FITS &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == PP_EXPERIMENT_ID) &&
              setequal(manifest$data_id, PP_DATA_IDS) &&
              !anyDuplicated(manifest$fit_id) && !anyDuplicated(manifest$fit_seed),
            "Pilot fit manifest identity failed.")
  g12 <- manifest$method_id == "G12"
  pooled <- manifest$method_id == "pooled_bayesSYNC"
  pp_assert(sum(g12) == 24L && sum(pooled) == 2L &&
              all(table(manifest$data_id[g12]) == 12L) &&
              all(table(manifest$data_id[pooled]) == 1L) &&
              all(manifest$n_cpus == 1L) &&
              all(manifest$planned_annealing_sweeps == 99L) &&
              all(manifest$planned_ordinary_t1_sweeps[g12] == 400L) &&
              all(manifest$planned_ordinary_t1_sweeps[pooled] == 200L) &&
              all(manifest$objective_eligibility_required) &&
              all(manifest$truth_free_winner_candidate == g12) &&
              all(manifest$cross_model_elbo_comparison_forbidden) &&
              !any(manifest$actual_racing) && !any(manifest$continuation) &&
              !any(manifest$automatic_800) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              all(manifest$engineering_pilot) &&
              !any(manifest$formal_v0lv_result) &&
              !any(manifest$formal_paper_mc_result),
            "Pilot fit route, budget, or isolation contract failed.")
  invisible(TRUE)
}

pp_load_manifests <- function() {
  data <- pp_read_csv(file.path(PP_ROOT, "DATA_MANIFEST.csv"))
  fits <- pp_read_csv(file.path(PP_ROOT, "FIT_MANIFEST.csv"))
  pp_validate_data_manifest(data)
  pp_validate_fit_manifest(fits)
  list(data = data, fits = fits)
}

pp_base_config <- function(data_row, smoke = FALSE) {
  if (smoke) {
    return(list(
      S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 3L, n_obs = 5L,
      common_grid = FALSE, n_obs_range = c(4L, 6L),
      sigma_eps = 0.3, bool_sparse_loadings = TRUE,
      prop_sparse = 0.8, score_var_decay = TRUE,
      n_dense = 51L, mean_amp = 0.6, beta_amp = 0,
      bs_degree = c(2L, 3L), identified_loadings = TRUE,
      sparsity_mode = "fixed"
    ))
  }
  list(
    S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = 5L, n_obs = 8L,
    common_grid = FALSE, n_obs_range = c(6L, 9L),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE,
    prop_sparse = 0.9, score_var_decay = TRUE,
    n_dense = 201L, mean_amp = 0.6, beta_amp = 0,
    bs_degree = c(2L, 3L), identified_loadings = TRUE,
    sparsity_mode = "fixed"
  )
}

pp_abs_cosine <- function(left, right) {
  denominator <- sqrt(sum(left^2) * sum(right^2))
  if (!is.finite(denominator) || denominator <= 1e-14) return(NA_real_)
  abs(sum(left * right) / denominator)
}

pp_apply_weak_loading_transform <- function(data, rebuild, rho = 0.6) {
  truth <- data$true_params
  pp_assert(truth$L_f == 1L && all(truth$L_s_by_study == 1L) &&
              isTRUE(truth$identified_loadings),
            "Weak-loading transform requires identified 1/1/1 truth.")
  a <- as.numeric(truth$a_true[, 1L])
  support <- which(a != 0)
  a_unit <- a / sqrt(sum(a^2))
  observed <- numeric(truth$S)
  for (study in seq_len(truth$S)) {
    old <- as.numeric(truth$b_true[[study]][, 1L])
    values <- old[old != 0]
    pp_assert(length(values) == length(support),
              "Weak transform requires equal active counts.")
    raw <- numeric(truth$p)
    raw[support] <- values
    orthogonal <- raw - a_unit * sum(a_unit * raw)
    norm <- sqrt(sum(orthogonal^2))
    if (!is.finite(norm) || norm <= 1e-12) {
      raw[support] <- c(tail(values, -1L), values[[1L]])
      orthogonal <- raw - a_unit * sum(a_unit * raw)
      norm <- sqrt(sum(orthogonal^2))
    }
    pp_assert(is.finite(norm) && norm > 1e-12,
              "Could not construct weak loading direction.")
    new <- sqrt(sum(old^2)) *
      (rho * a_unit + sqrt(1 - rho^2) * orthogonal / norm)
    new[-support] <- 0
    truth$b_true[[study]][, 1L] <- new
    observed[[study]] <- pp_abs_cosine(a, new)
  }
  truth$gamma_b_true <- lapply(truth$b_true, function(value) 1L * (value != 0))
  truth$omega_b_true <- lapply(truth$gamma_b_true, colMeans)
  truth$identified_loadings <- FALSE
  truth$paper_pilot_loading_scenario <- "weak_loading_separation"
  truth$paper_pilot_target_abs_cosine <- rho
  truth$paper_pilot_observed_abs_cosine <- observed
  data$true_params <- truth
  rebuild(data)
}

pp_truth_assertions <- function(data, data_row, config) {
  truth <- data$true_params
  counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  a <- as.numeric(truth$a_true[, 1L])
  b <- lapply(truth$b_true, function(value) as.numeric(value[, 1L]))
  cosines <- vapply(b, function(value) pp_abs_cosine(a, value), numeric(1L))
  expected_active <- as.integer(round((1 - config$prop_sparse) * config$p))
  checks <- data.frame(
    assertion_id = c(
      "dimensions", "factor_fpca_counts", "observation_counts",
      "active_loading_counts", "loading_norms", "scenario_geometry",
      "signal_noise_reconstruction"
    ),
    passed = FALSE,
    detail = "",
    stringsAsFactors = FALSE
  )
  checks$passed[1L] <- truth$S == 2L && all(truth$n_s == config$n_s) &&
    truth$p == config$p && truth$d_use == 0L
  checks$detail[1L] <- "S=2;n=(30,30) or smoke;p and d match"
  checks$passed[2L] <- truth$L_f == 1L &&
    all(truth$L_s_by_study == 1L) && identical(as.integer(truth$M_f), 2L) &&
    all(vapply(truth$M_s, function(value) as.integer(value) == 2L, logical(1L)))
  checks$detail[2L] <- "L=1/1/1;M=2"
  checks$passed[3L] <- all(counts >= config$n_obs_range[[1L]] &
                            counts <= config$n_obs_range[[2L]])
  checks$detail[3L] <- paste(range(counts), collapse = "-")
  active <- c(sum(a != 0), vapply(b, function(value) sum(value != 0), integer(1L)))
  checks$passed[4L] <- all(active == expected_active)
  checks$detail[4L] <- paste(active, collapse = ";")
  norms <- c(sqrt(sum(a^2)), vapply(b, function(x) sqrt(sum(x^2)), numeric(1L)))
  checks$passed[5L] <- max(abs(norms - norms[[1L]])) <= 1e-10
  checks$detail[5L] <- paste(format(norms, digits = 8), collapse = ";")
  weak <- data_row$scenario_id[[1L]] == "weak_loading_separation"
  checks$passed[6L] <- if (weak) max(abs(cosines - 0.6)) <= 1e-10 else
    max(cosines) <= 1e-12
  checks$detail[6L] <- paste(format(cosines, digits = 8), collapse = ";")
  maximum_error <- 0
  for (study in seq_len(truth$S)) for (subject in seq_len(truth$n_s[[study]])) {
    signal <- truth$mu_true_values[[study]][[subject]] +
      truth$beta_true_values[[study]][[subject]] +
      truth$f_true_values[[study]][[subject]] %*% t(truth$a_true) +
      truth$g_true_values[[study]][[subject]] %*% t(truth$b_true[[study]])
    observed <- do.call(cbind, data$Y[[study]][[subject]])
    maximum_error <- max(maximum_error,
      abs(signal - truth$signal_true_values[[study]][[subject]]),
      abs(observed - signal - truth$noise_true_values[[study]][[subject]]))
  }
  checks$passed[7L] <- maximum_error <= 1e-12
  checks$detail[7L] <- format(maximum_error, scientific = TRUE)
  checks
}

pp_generate_one <- function(data_row, output_root, smoke = FALSE) {
  pp_load_package()
  pp_source_runtime(FALSE)
  output_dir <- file.path(output_root, "data", data_row$data_id[[1L]])
  if (dir.exists(output_dir)) {
    pp_assert(file.exists(file.path(output_dir, "DATA_SEAL_COMPLETE.txt")) &&
                file.exists(file.path(output_dir, "observation_bundle.rds")) &&
                file.exists(file.path(output_dir, "sealed_truth", "truth_bundle.rds")),
              paste0("Incomplete pre-existing data directory: ", output_dir))
    return(invisible(FALSE))
  }
  calibration <- new.env(parent = globalenv())
  sys.source(pp_find_one(PP_SNAPSHOT, "^v0f_data_calibration[.]R$"),
             envir = calibration)
  config <- pp_base_config(data_row, smoke)
  arguments <- config
  arguments$seed <- as.integer(data_row$data_seed[[1L]])
  generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
  captured <- c2_capture_conditions(do.call(generator, arguments),
                                    phase = "pilot_data_generation")
  pp_assert(is.null(captured$error), paste0(
    "Data generation failed: ", captured$error$message %||% "unknown"
  ))
  data <- calibration$calibrate_loading_information_v0f(
    captured$value, regime = "constant_per_active_loading",
    anchor_p = 100L, active_fraction = 0.1, anchor_norm = 3
  )
  if (data_row$scenario_id[[1L]] == "weak_loading_separation") {
    data <- pp_apply_weak_loading_transform(
      data, calibration$v0f_rebuild_observations, rho = 0.6
    )
  }
  assertions <- pp_truth_assertions(data, data_row, config)
  pp_assert(all(assertions$passed), paste0(
    "Truth assertions failed: ",
    paste(assertions$assertion_id[!assertions$passed], collapse = ";")
  ))
  observation <- list(
    bundle_class = "v0lv_observation_only_candidate2",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_row$data_id[[1L]],
    data_seed = as.integer(data_row$data_seed[[1L]]),
    scenario_id = data_row$scenario_id[[1L]],
    scenario_label_available_to_fit = FALSE,
    engineering_pilot = TRUE,
    formal_experiment = FALSE,
    Y = data$Y, time_obs = data$time_obs, Z = data$Z,
    dimensions = config[c("S", "n_s", "p", "d", "L_f", "L_s",
                          "M_f", "M_s", "K")]
  )
  pp_assert(!length(candidate2_forbidden_observation_names(observation)),
            "Observation bundle leaks a forbidden truth field.")
  sealed <- list(
    bundle_class = "v0lv_sealed_truth_candidate2",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_row$data_id[[1L]],
    data_seed = as.integer(data_row$data_seed[[1L]]),
    scenario_id = data_row$scenario_id[[1L]],
    engineering_pilot = TRUE,
    formal_experiment = FALSE,
    data = data,
    sealed_generation_provenance = list(
      generator_config = config,
      generation_warnings = captured$warnings,
      truth_assertions = assertions,
      generated_at_utc = captured$ended_at_utc
    )
  )
  dir.create(file.path(output_dir, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  observation_path <- file.path(output_dir, "observation_bundle.rds")
  truth_path <- file.path(output_dir, "sealed_truth", "truth_bundle.rds")
  pp_save_rds(observation, observation_path)
  pp_save_rds(sealed, truth_path)
  pp_write_csv(assertions, file.path(output_dir, "sealed_truth",
                                     "truth_assertions.csv"))
  pp_write_lines(c(
    "status=DATA_SEAL_COMPLETE",
    paste0("data_id=", data_row$data_id[[1L]]),
    paste0("observation_sha256=", pp_sha256(observation_path)),
    paste0("sealed_truth_sha256=", pp_sha256(truth_path)),
    "truth_unseal_authorized=FALSE", "engineering_pilot=TRUE",
    "formal_paper_mc_result=FALSE"
  ), file.path(output_dir, "DATA_SEAL_COMPLETE.txt"))
  invisible(TRUE)
}

pp_load_observation <- function(data_id, output_root) {
  path <- file.path(output_root, "data", data_id, "observation_bundle.rds")
  pp_assert(file.exists(path), paste0("Missing observation bundle: ", data_id))
  value <- readRDS(path)
  pp_assert(identical(value$bundle_class, "v0lv_observation_only_candidate2") &&
              identical(value$data_id, data_id) &&
              !length(candidate2_forbidden_observation_names(value)),
            paste0("Observation identity/isolation failed: ", data_id))
  value
}

pp_fit_g12 <- function(observation, fit_row, smoke = FALSE) {
  p <- length(observation$Y[[1L]][[1L]])
  if (smoke) {
    anneal <- c(1, 1.2, 3); maxit <- 6L; n_g <- 31L
    tol_abs <- 0; tol_rel <- 0; convergence_rule <- "parameters"
    practical_control <- NULL; trace_sweeps <- 1L
  } else {
    anneal <- c(1, 1.9, 100); maxit <- 499L; n_g <- 51L
    tol_abs <- 1e-3; tol_rel <- 1e-5; convergence_rule <- "practical"
    practical_control <- multiFSYNC::g12_stopping_control()
    practical_control$checkpoints <- c(380L, 400L)
    trace_sweeps <- c(1L, 380L, 400L)
  }
  arguments <- list(
    Y = observation$Y, Z = observation$Z, time_obs = observation$time_obs,
    L_f = observation$dimensions$L_f, L_s = observation$dimensions$L_s,
    M_f = observation$dimensions$M_f, M_s = observation$dimensions$M_s,
    K = if (smoke) 3L else observation$dimensions$K,
    anneal = anneal, list_hyper = NULL, n_g = n_g, time_g = NULL,
    tol_abs = tol_abs, tol_rel = tol_rel, maxit = maxit,
    n_cpus = 1L, verbose = FALSE, seed = as.integer(fit_row$fit_seed[[1L]]),
    bool_scale = FALSE, bool_var_spec_prob = FALSE, d_0 = as.integer(p),
    convergence_rule = convergence_rule, lambda_orth = 0,
    practical_control = practical_control, initialization = "random",
    initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ),
    continuation_state = NULL, pre_score_sweeps = 1L,
    trace_sweeps = trace_sweeps,
    function_initialization = "gram_unit_energy"
  )
  captured <- c2_capture_conditions(
    do.call(getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score"),
            arguments),
    phase = if (smoke) "pilot_smoke_g12" else "pilot_g12_fixed400"
  )
  fit <- captured$value
  validation_error <- NULL
  objective <- list(eligible = FALSE, checks = logical(),
                    invalid_reasons = "terminal_error")
  if (is.null(captured$error)) {
    validation_error <- tryCatch({
      pp_assert(is.list(fit), "G12 did not return a fit list.")
      pp_assert(identical(fit$initialization, "random") &&
                  fit$n_cpus_used == 1L &&
                  identical(as.numeric(fit$lambda_orth), 0) &&
                  identical(fit$pre_score_interface$function_initialization,
                            "gram_unit_energy") &&
                  fit$pre_score_interface$pre_score_sweeps == 1L,
                "G12 returned route provenance differs from manifest.")
      calibration <- fit$random_scale_calibration_diagnostics
      pp_assert(is.data.frame(calibration) && nrow(calibration) == 6L &&
                  all(calibration$success) &&
                  max(abs(calibration$after - 1)) < 1e-10,
                "G12 unit Gram-energy calibration failed.")
      if (!smoke) {
        pp_assert(fit$annealing_sweeps == 99L && fit$t1_sweeps == 400L &&
                    fit$n_g == 51L &&
                    identical(fit$practical_control, practical_control) &&
                    identical(names(fit$practical_checkpoints), c("380", "400")),
                  "G12 fixed-400 endpoint/checkpoints are incomplete.")
        objective <- r_route_v3_candidate2_objective_status(fit)
      } else {
        checks <- c(
          finite_elbo = length(fit$ELBO) > 0L && all(is.finite(fit$ELBO)),
          finite_state = all(is.finite(fit$mu_q_a))
        )
        objective <- list(eligible = all(checks), checks = checks,
                          invalid_reasons = names(checks)[!checks])
      }
      pp_assert(objective$eligible, paste0(
        "G12 endpoint ineligible: ",
        paste(objective$invalid_reasons, collapse = ";")
      ))
      NULL
    }, error = function(condition) condition)
  }
  if (!is.null(validation_error)) {
    captured$error <- list(
      class = paste(class(validation_error), collapse = ";"),
      message = conditionMessage(validation_error),
      call = "pp_fit_g12", trace_summary = "route validation"
    )
    fit <- NULL
  }
  list(fit = fit, captured = captured, objective = objective)
}

pp_fit_pooled <- function(observation, fit_row, smoke = FALSE) {
  reference <- candidate2_load_bayes_reference(PP_SNAPSHOT)
  pooled_Y <- unname(unlist(observation$Y, recursive = FALSE))
  pooled_time <- unname(unlist(observation$time_obs, recursive = FALSE))
  pp_assert(length(pooled_Y) == sum(observation$dimensions$n_s) &&
              length(pooled_time) == length(pooled_Y),
            "Pooled subject concatenation failed.")
  spec <- if (smoke) list(
    anneal = c(1, 1.2, 3), annealing_sweeps = 2L,
    maxit = 6L, n_g = 31L, K = 3L, tol_abs = 0, tol_rel = 0
  ) else list(
    anneal = c(1, 1.9, 100), annealing_sweeps = 99L,
    maxit = 299L, n_g = 51L, K = 5L, tol_abs = 1e-3, tol_rel = 1e-5
  )
  captured <- c2_capture_conditions(reference$bayesSYNC(
    time_obs = pooled_time, Y = pooled_Y, Q = 3L, L = 2L, K = spec$K,
    anneal = spec$anneal, list_hyper = NULL, n_g = spec$n_g, time_g = NULL,
    tol_abs = spec$tol_abs, tol_rel = spec$tol_rel, maxit = spec$maxit,
    n_cpus = 1L, verbose = FALSE,
    seed = as.integer(fit_row$fit_seed[[1L]]), bool_scale = FALSE,
    bool_var_spec_prob = FALSE, show_factor_ppi_progress = FALSE
  ), phase = if (smoke) "pilot_smoke_pooled" else "pilot_pooled_bayesSYNC")
  fit <- captured$value
  checks <- c(
    no_error = is.null(captured$error),
    valid_dimensions = is.list(fit) && identical(as.integer(fit$Q), 3L) &&
      identical(as.integer(fit$L), 2L),
    finite_elbo = is.list(fit) && length(fit$ELBO_iter) > 0L &&
      all(is.finite(fit$ELBO_iter)),
    finite_loadings = is.list(fit) && all(is.finite(fit$B_hat)),
    finite_factor_ppi = is.list(fit) && all(is.finite(fit$factor_ppi))
  )
  objective <- list(eligible = all(checks), checks = checks,
                    invalid_reasons = names(checks)[!checks])
  if (!objective$eligible) fit <- NULL
  list(fit = fit, captured = captured, objective = objective,
       spec = spec)
}

pp_make_terminal_record <- function(result, fit_row, smoke = FALSE) {
  fit <- result$fit
  pooled <- fit_row$method_id[[1L]] == "pooled_bayesSYNC"
  captured <- result$captured
  actual_t1 <- if (is.null(fit)) NA_integer_ else if (pooled) {
    max(0L, as.integer(fit$i_iter) - result$spec$annealing_sweeps)
  } else as.integer(fit$t1_sweeps)
  practical <- if (is.null(fit) || pooled || smoke) FALSE else
    isTRUE(fit$practical_converged)
  list(
    terminal_schema = "G12_PAPER_PILOT_TERMINAL_V1",
    experiment_id = PP_EXPERIMENT_ID,
    fit_id = fit_row$fit_id[[1L]], data_id = fit_row$data_id[[1L]],
    scenario_id = fit_row$scenario_id[[1L]],
    method_id = fit_row$method_id[[1L]], route_id = fit_row$route_id[[1L]],
    seed_index = as.integer(fit_row$seed_index[[1L]]),
    fit_seed = as.integer(fit_row$fit_seed[[1L]]),
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = if (is.null(fit)) "error" else if (pooled) {
      "max_budget_reached"
    } else if (smoke) "smoke_complete" else "fixed_400_complete",
    error = captured$error %||% list(class = "", message = "", call = "",
                                     trace_summary = ""),
    warnings = captured$warnings,
    actual_annealing_sweeps = if (is.null(fit)) NA_integer_ else if (pooled) {
      result$spec$annealing_sweeps
    } else as.integer(fit$annealing_sweeps),
    actual_T1_sweeps = actual_t1,
    objective_eligible = isTRUE(result$objective$eligible),
    objective_checks = result$objective$checks,
    objective_invalid_reasons = result$objective$invalid_reasons,
    final_elbo = if (is.null(fit)) NA_real_ else if (pooled) {
      as.numeric(tail(fit$ELBO_iter, 1L))
    } else as.numeric(tail(fit$ELBO, 1L)),
    practical_converged = practical,
    convergence_status = if (is.null(fit)) "error" else if (pooled) {
      "reference_endpoint"
    } else as.character(fit$convergence_status),
    studies_combined_once = pooled,
    study_labels_passed_to_model = FALSE,
    cross_model_elbo_comparison_forbidden = TRUE,
    truth_used_for_fit_stopping_or_selection = FALSE,
    continuation_used = FALSE, automatic_800_used = FALSE,
    n_cpus = 1L, engineering_pilot = TRUE,
    smoke = smoke, formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE
  )
}

pp_terminal_row <- function(record) {
  data.frame(
    fit_id = record$fit_id, data_id = record$data_id,
    scenario_id = record$scenario_id, method_id = record$method_id,
    route_id = record$route_id, seed_index = record$seed_index,
    fit_seed = record$fit_seed, terminal_status = record$terminal_status,
    actual_annealing_sweeps = record$actual_annealing_sweeps,
    actual_T1_sweeps = record$actual_T1_sweeps,
    objective_eligible = record$objective_eligible,
    final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    convergence_status = record$convergence_status,
    warning_count = if (is.data.frame(record$warnings)) nrow(record$warnings)
      else 0L,
    elapsed_seconds = record$elapsed_seconds,
    peak_memory_bytes = record$peak_memory_bytes,
    studies_combined_once = record$studies_combined_once,
    study_labels_passed_to_model = record$study_labels_passed_to_model,
    truth_used_for_fit_stopping_or_selection =
      record$truth_used_for_fit_stopping_or_selection,
    continuation_used = record$continuation_used,
    engineering_pilot = record$engineering_pilot,
    formal_paper_mc_result = record$formal_paper_mc_result,
    stringsAsFactors = FALSE
  )
}

pp_select_g12_winners <- function(terminals) {
  g12 <- terminals[terminals$method_id == "G12", , drop = FALSE]
  pp_assert(nrow(g12) > 0L &&
              !any(g12$truth_used_for_fit_stopping_or_selection),
            "G12 selection input is empty or truth-bearing.")
  groups <- split(g12, g12$data_id)
  winners <- lapply(groups, function(group) {
    eligible <- group$objective_eligible & is.finite(group$final_elbo) &
      group$terminal_status %in% c("fixed_400_complete", "smoke_complete")
    pp_assert(any(eligible), paste0(
      "No eligible G12 endpoint for ", group$data_id[[1L]]
    ))
    candidates <- group[eligible, , drop = FALSE]
    candidates <- candidates[order(-candidates$final_elbo,
                                   candidates$fit_id), , drop = FALSE]
    winner <- candidates[1L, , drop = FALSE]
    winner$selection_rule <- "maximum_eligible_ordinary_T1_ELBO_within_data"
    winner$tie_break_rule <- "fit_id_lexicographic"
    winner
  })
  result <- do.call(rbind, winners)
  rownames(result) <- NULL
  result[order(result$data_id), , drop = FALSE]
}

pp_rbind_fill <- function(values) {
  values <- values[vapply(values, is.data.frame, logical(1L))]
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(value) {
    missing <- setdiff(columns, names(value))
    for (name in missing) value[[name]] <- NA
    value[columns]
  })
  do.call(rbind, values)
}

pp_nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate); truth <- as.numeric(truth)
  if (length(estimate) != length(truth) || !length(truth) ||
      any(!is.finite(c(estimate, truth)))) return(NA_real_)
  denominator <- sqrt(mean(truth^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sqrt(mean((estimate - truth)^2)) / denominator
}

pp_vector_l2 <- function(value) sqrt(sum(as.numeric(value)^2))
pp_vector_rms <- function(value) sqrt(mean(as.numeric(value)^2))

pp_multi_additional_metrics <- function(
    fit, truth, evaluation, evaluator_environment) {
  data <- truth$data
  reported_paths <- evaluator_environment$v0g_reported_paths_by_study(fit, data)
  true_paths <- evaluator_environment$v0h_truth_paths_by_study(data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected", , drop = FALSE]
  pp_assert(nrow(functional) == 4L,
            "Expected four study-role truth rows for multiFSYNC.")
  roles <- lapply(seq_len(nrow(functional)), function(index) {
    item <- functional[index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_role <- item$true_role[[1L]]
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_role == "shared") {
      data$true_params$a_true[, true_factor]
    } else data$true_params$b_true[[study]][, true_factor]
    true_path <- true_paths[[true_role]][[study]][, true_factor]
    true_contribution <- as.numeric(outer(true_path, true_loading))
    matched <- isTRUE(item$matched_by_loading[[1L]])
    if (matched) {
      estimated_role <- item$estimated_block[[1L]]
      estimated_factor <- as.integer(item$estimated_factor[[1L]])
      estimated_loading <- if (estimated_role == "shared") {
        fit$mu_q_a_hat[, estimated_factor]
      } else fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      estimated_path <-
        reported_paths[[estimated_role]][[study]][, estimated_factor]
      aligned_path <- item$loading_sign[[1L]] * estimated_path
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
      estimated_role <- NA_character_; estimated_factor <- NA_integer_
      process_error <- contribution_error <- loading_norm_ratio <-
        loading_cosine <- estimated_path_rms <- NA_real_
    }
    data.frame(
      study = study, true_role = true_role,
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
  process <- do.call(rbind, roles)
  finite_mean <- function(value) {
    value <- as.numeric(value); value <- value[is.finite(value)]
    if (length(value)) mean(value) else NA_real_
  }
  process_summary <- data.frame(
    study_role_matched = sum(process$matched_by_loading),
    study_role_truth_total = nrow(process),
    factor_process_nrmse_mean_matched =
      finite_mean(process$factor_process_nrmse_matched),
    factor_process_nrmse_mean_missing_as_zero =
      mean(process$factor_process_nrmse_missing_as_zero),
    complete_contribution_nrmse_mean_matched =
      finite_mean(process$complete_contribution_nrmse_matched),
    complete_contribution_nrmse_mean_missing_as_zero =
      mean(process$complete_contribution_nrmse_missing_as_zero),
    stringsAsFactors = FALSE
  )
  loading_scale <- do.call(rbind, lapply(c("A", "B1", "B2"), function(role) {
    if (role == "A") {
      raw <- fit$mu_q_a_original[, 1L]
      canonical <- fit$mu_q_a_hat[, 1L]
      true <- data$true_params$a_true[, 1L]
      factor_scale <- fit$factor_scale_shared[[1L]]
      factor_ppi <- fit$factor_ppi_shared[[1L]]
    } else {
      study <- if (role == "B1") 1L else 2L
      raw <- fit$mu_q_b_specific_original[[study]][, 1L]
      canonical <- fit$mu_q_b_specific_hat[[study]][, 1L]
      true <- data$true_params$b_true[[study]][, 1L]
      factor_scale <- fit$factor_scale_specific[[study]][[1L]]
      factor_ppi <- fit$factor_ppi_specific[[study]][[1L]]
    }
    data.frame(
      block_role = role, factor_ppi = factor_ppi,
      factor_scale = factor_scale,
      raw_loading_norm_ratio = pp_vector_l2(raw) / pp_vector_l2(true),
      canonical_loading_norm_ratio =
        pp_vector_l2(canonical) / pp_vector_l2(true),
      raw_loading_abs_cosine = pp_abs_cosine(raw, true),
      canonical_loading_abs_cosine = pp_abs_cosine(canonical, true),
      stringsAsFactors = FALSE
    )
  }))
  loading_scale_summary <- data.frame(
    raw_loading_norm_ratio_mean = mean(loading_scale$raw_loading_norm_ratio),
    canonical_loading_norm_ratio_mean =
      mean(loading_scale$canonical_loading_norm_ratio),
    raw_loading_abs_cosine_mean = mean(loading_scale$raw_loading_abs_cosine),
    canonical_loading_abs_cosine_mean =
      mean(loading_scale$canonical_loading_abs_cosine),
    factor_scale_mean = mean(loading_scale$factor_scale),
    stringsAsFactors = FALSE
  )
  list(
    process_contribution_by_role = process,
    process_contribution_summary = process_summary,
    loading_scale_by_block = loading_scale,
    loading_scale_summary = loading_scale_summary
  )
}

pp_evaluate_fit <- function(
    fit, record, truth, evaluator_environment, bayes_environment) {
  if (record$method_id == "G12") {
    result <- candidate2_evaluate_multi(
      fit, truth, PP_SNAPSHOT, evaluator_environment = evaluator_environment,
      mode = "smoke"
    )
    additional <- pp_multi_additional_metrics(
      fit, truth, result, evaluator_environment
    )
    c(result, additional)
  } else {
    candidate2_evaluate_pooled(
      fit, truth, PP_SNAPSHOT, bayes_environment = bayes_environment,
      mode = "smoke"
    )
  }
}

pp_pooled_compact_summary <- function(evaluation) {
  direction <- evaluation$direction_metrics
  feature <- evaluation$feature_and_score_metrics
  reconstruction <- evaluation$reconstruction_metrics
  data.frame(
    observed_signal_nrmse = reconstruction$observed_signal_nrmse[[1L]],
    dense_grid_signal_nrmse = reconstruction$dense_grid_signal_nrmse[[1L]],
    loading_abs_cosine_mean = mean(direction$loading_abs_cosine, na.rm = TRUE),
    loading_relative_l2_error_mean =
      mean(direction$loading_relative_l2_error, na.rm = TRUE),
    factor_process_abs_cor_mean =
      mean(direction$factor_process_abs_cor, na.rm = TRUE),
    factor_process_nrmse_mean =
      mean(direction$factor_process_nrmse, na.rm = TRUE),
    complete_contribution_nrmse_mean =
      mean(direction$complete_contribution_nrmse, na.rm = TRUE),
    inactive_study_leakage_energy_fraction_mean =
      mean(direction$inactive_study_leakage_energy_fraction, na.rm = TRUE),
    feature_total_ise_mean = mean(feature$feature_total_ise, na.rm = TRUE),
    feature_projection_floor_ise_mean =
      mean(feature$feature_projection_floor_ise, na.rm = TRUE),
    feature_excess_ise_mean = mean(feature$feature_excess_ise, na.rm = TRUE),
    fpca_score_abs_cor_mean = mean(feature$fpca_score_abs_cor, na.rm = TRUE),
    loading_directions_matched = nrow(direction),
    feature_components_matched = nrow(feature),
    stringsAsFactors = FALSE
  )
}
