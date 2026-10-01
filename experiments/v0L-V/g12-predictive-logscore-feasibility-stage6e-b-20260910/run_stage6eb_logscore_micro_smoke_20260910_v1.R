#!/usr/bin/env Rscript

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- command_line[startsWith(command_line, "--file=")]
script_file <- if (length(file_argument)) {
  sub("--file=", "", file_argument[[1L]], fixed = TRUE)
} else {
  "run_stage6eb_logscore_micro_smoke_20260910_v1.R"
}
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
repo_root <- file.path(script_root, "..", "..", "..")
if (!dir.exists(repo_root)) {
  stop("Cannot resolve the multiFSYNC source-tree root.", call. = FALSE)
}
stage6b_root <- file.path(
  repo_root,
  "experiments", "v0L-V",
  "g12-truth-free-dimension-selection-stage6b-20260901"
)
source(file.path(
  script_root, "stage6eb_logscore_tools_20260910_v1.R"
), local = TRUE)
source(file.path(
  stage6b_root, "stage6b_truth_free_tools_20260901_v1.R"
), local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
if (!length(output_argument)) {
  stop("--output-root is required.", call. = FALSE)
}
output_root <- sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Smoke output root already exists; refusing to overwrite it.",
       call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

library_argument <- arguments[startsWith(arguments, "--library-root=")]
if (length(library_argument)) {
  library_root <- sub(
    "--library-root=", "", library_argument[[1L]], fixed = TRUE
  )
  s6eb_assert(dir.exists(library_root),
               "--library-root does not identify an installed R library.")
  .libPaths(c(library_root, .libPaths()))
  s6eb_assert(requireNamespace("multiFSYNC", quietly = TRUE),
               "multiFSYNC is not installed in --library-root.")
  package_loading_mode <- "preinstalled_temporary_library"
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required when --library-root is not supplied.",
         call. = FALSE)
  }
  pkgload::load_all(
    repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE
  )
  package_loading_mode <- "pkgload_source_tree"
}

started <- proc.time()[["elapsed"]]
generated <- multiFSYNC::simulate_multi_study_structured(
  S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
  K = 3L, n_obs = 6L, common_grid = FALSE, n_obs_range = c(6L, 7L),
  sigma_eps = 0.3, bool_sparse_loadings = TRUE, prop_sparse = 0.8,
  score_var_decay = TRUE, n_dense = 51L, mean_amp = 0.6, beta_amp = 0,
  bs_degree = c(2L, 3L), identified_loadings = TRUE,
  sparsity_mode = "fixed", seed = 91061001L
)
observation <- list(
  bundle_class = "stage6eb_micro_observation_only_v1",
  Y = generated$Y,
  time_obs = generated$time_obs,
  Z = generated$Z,
  dimensions = list(
    S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = 3L
  )
)
rm(generated)
invisible(gc(FALSE))
s6b_assert(!length(s6b_forbidden_observation_names(observation)),
           "Micro observation bundle leaked a truth field.")

split <- s6b_make_middle_time_holdout(observation, min_train_times = 5L)
rm(observation)
invisible(gc(FALSE))

fit_warnings <- character()
fit <- withCallingHandlers(
  multiFSYNC::bayesSYNC_multi_pre_score(
    Y = split$train_observation$Y,
    Z = split$train_observation$Z,
    time_obs = split$train_observation$time_obs,
    L_f = 1L,
    L_s = c(1L, 1L),
    M_f = 2L,
    M_s = list(2L, 2L),
    K = 3L,
    anneal = c(1, 1.2, 3),
    list_hyper = NULL,
    n_g = 31L,
    time_g = NULL,
    tol_abs = 0,
    tol_rel = 0,
    maxit = 6L,
    n_cpus = 1L,
    verbose = FALSE,
    seed = 91061101L,
    bool_scale = FALSE,
    bool_var_spec_prob = FALSE,
    d_0 = 20L,
    convergence_rule = "parameters",
    lambda_orth = 0,
    practical_control = NULL,
    initialization = "random",
    initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ),
    continuation_state = NULL,
    pre_score_sweeps = 1L,
    trace_sweeps = 1L,
    function_initialization = "gram_unit_energy"
  ),
  warning = function(condition) {
    fit_warnings <<- c(fit_warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)

s6eb_assert(is.list(fit) && fit$n_cpus_used == 1L,
              "Micro fit did not complete on one CPU.")
s6eb_assert(identical(fit$pre_score_interface$function_initialization,
                       "gram_unit_energy"),
              "Micro fit did not retain the G12 initialization provenance.")

# The preregistered Stage 6B point-risk implementation is called without any
# modification.  Stage 6E-B adds a parallel experimental score only.
registered_rmse <- s6b_score_holdout(fit, split$holdout_plan, "all")
experimental_logscore <- s6eb_score_holdout(fit, split$holdout_plan)

first_record <- split$holdout_plan$records[[1L]]
first_moments <- s6eb_predictive_moments(
  fit, first_record$study, first_record$subject, first_record$time
)
efficient <- s6eb_lowrank_gaussian_nlpd(
  first_record$y,
  first_moments$mean,
  first_moments$diagonal_variance,
  first_moments$low_rank
)
dense <- s6eb_dense_gaussian_nlpd(
  first_record$y,
  first_moments$mean,
  first_moments$diagonal_variance,
  first_moments$low_rank
)
woodbury_dense_difference <- abs(
  efficient$total_negative_log_density -
    dense$total_negative_log_density
)

registered_subject <- registered_rmse$subject_time
experimental_subject <- experimental_logscore$subject_time
mean_loss_difference <- max(abs(
  registered_subject$mse - experimental_subject$mse
))
rmse_difference <- abs(
  registered_rmse$aggregate$heldout_rmse -
    experimental_logscore$aggregate$heldout_rmse
)

aggregate <- experimental_logscore$aggregate
checks <- data.frame(
  check_id = c(
    "fit_completed",
    "one_cpu",
    "g12_initialization",
    "eight_subject_time_vectors",
    "twenty_variables_per_vector",
    "registered_rmse_unchanged",
    "posterior_mean_preserved",
    "woodbury_matches_dense",
    "finite_log_scores",
    "positive_predictive_variances",
    "no_variance_floor_or_clamp",
    "coverage_computation_bounded",
    "moment_gaussian_label_honest",
    "not_promoted_to_primary",
    "truth_free"
  ),
  passed = c(
    is.list(fit),
    fit$n_cpus_used == 1L,
    identical(fit$pre_score_interface$function_initialization,
              "gram_unit_energy"),
    nrow(experimental_subject) == 8L,
    all(experimental_subject$variables == 20L),
    identical(
      registered_rmse$aggregate$predictive_variance_scope,
      "residual_noise_only_not_calibrated_full_posterior"
    ),
    mean_loss_difference < 1e-10 && rmse_difference < 1e-10,
    woodbury_dense_difference < 1e-10,
    all(is.finite(experimental_subject$joint_nlpd_per_variable)) &&
      all(is.finite(experimental_subject$marginal_nlpd)),
    all(experimental_subject$minimum_predictive_variance > 0),
    aggregate$total_diagonal_floor_count == 0L &&
      aggregate$total_numerical_clamp_count == 0L,
    is.finite(aggregate$marginal_interval_coverage) &&
      aggregate$marginal_interval_coverage >= 0 &&
      aggregate$marginal_interval_coverage <= 1,
    identical(aggregate$exact_posterior_predictive_density, FALSE),
    identical(aggregate$primary_dimension_score_ready, FALSE) &&
      identical(aggregate$replaces_preregistered_rmse, FALSE),
    identical(aggregate$truth_used, FALSE)
  ),
  detail = c(
    "micro fit object returned",
    "n_cpus_used=1",
    "gram_unit_energy + one pre-score + annealing",
    as.character(nrow(experimental_subject)),
    paste(range(experimental_subject$variables), collapse = "-"),
    "historical Stage 6B function sourced and called as-is",
    paste0("max_mse_difference=", format(mean_loss_difference, digits = 6),
           "; rmse_difference=", format(rmse_difference, digits = 6)),
    paste0("absolute_nld_difference=",
           format(woodbury_dense_difference, digits = 6)),
    "joint and marginal negative log predictive densities finite",
    paste0("minimum=",
           format(min(experimental_subject$minimum_predictive_variance),
                  digits = 6)),
    paste0("floors=", aggregate$total_diagonal_floor_count,
           "; clamps=", aggregate$total_numerical_clamp_count),
    paste0("descriptive_coverage=",
           format(aggregate$marginal_interval_coverage, digits = 6)),
    "full mean-field moments with Gaussian approximation; not exact density",
    "experimental_not_preregistered_primary",
    "holdout_plan truth_used=FALSE"
  ),
  stringsAsFactors = FALSE
)

s6b_write_csv(checks, file.path(output_root, "SMOKE_STATUS.csv"))
s6b_write_csv(registered_rmse$aggregate,
              file.path(output_root, "REGISTERED_RMSE_SUMMARY.csv"))
s6b_write_csv(aggregate,
              file.path(output_root, "EXPERIMENTAL_LOGSCORE_SUMMARY.csv"))
s6b_write_csv(experimental_subject,
              file.path(output_root, "EXPERIMENTAL_SUBJECT_TIME.csv"))
s6b_write_csv(data.frame(
  warning_index = seq_along(fit_warnings),
  warning = fit_warnings,
  stringsAsFactors = FALSE
), file.path(output_root, "SMOKE_WARNINGS.csv"))

elapsed <- proc.time()[["elapsed"]] - started
status <- if (all(checks$passed)) {
  "STAGE6E_B_LOGSCORE_MICRO_SMOKE_PASS"
} else {
  "STAGE6E_B_LOGSCORE_MICRO_SMOKE_FAIL"
}
s6b_write_lines(c(
  paste0("status=", status),
  "formal_experiment=FALSE",
  "new_full_size_fit_started=FALSE",
  "truth_accessed=FALSE",
  "n_cpus=1",
  "subject_time_vectors=8",
  "variables_per_vector=20",
  paste0("warnings_captured=", length(fit_warnings)),
  paste0("elapsed_seconds=", format(elapsed, digits = 8)),
  paste0("package_loading_mode=", package_loading_mode),
  "preregistered_rmse_interface_modified=FALSE",
  "logscore_replaces_rmse=FALSE",
  "primary_dimension_score_ready=FALSE"
), file.path(output_root, "SMOKE_COMPLETE.txt"))

if (!all(checks$passed)) {
  failed <- paste(checks$check_id[!checks$passed], collapse = ", ")
  stop(paste0("Stage 6E-B micro smoke failed: ", failed), call. = FALSE)
}
cat(status, " output_root=", output_root, "\n", sep = "")
