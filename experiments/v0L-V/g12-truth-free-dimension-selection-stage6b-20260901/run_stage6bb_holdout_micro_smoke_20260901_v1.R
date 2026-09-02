#!/usr/bin/env Rscript

script_argument <- commandArgs(trailingOnly = FALSE)
file_argument <- script_argument[startsWith(script_argument, "--file=")]
script_file <- if (length(file_argument)) {
  sub("--file=", "", file_argument[[1L]], fixed = TRUE)
} else {
  "run_stage6bb_holdout_micro_smoke_20260901_v1.R"
}
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
repo_root <- file.path(script_root, "..", "..", "..")
if (!dir.exists(repo_root)) {
  stop("Cannot resolve the multiFSYNC source-tree root.", call. = FALSE)
}
source(file.path(script_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = TRUE)

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

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required for the source-tree micro smoke.", call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)

generated <- multiFSYNC::simulate_multi_study_structured(
  S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
  K = 3L, n_obs = 6L, common_grid = FALSE, n_obs_range = c(6L, 7L),
  sigma_eps = 0.3, bool_sparse_loadings = TRUE, prop_sparse = 0.8,
  score_var_decay = TRUE, n_dense = 51L, mean_amp = 0.6, beta_amp = 0,
  bs_degree = c(2L, 3L), identified_loadings = TRUE,
  sparsity_mode = "fixed", seed = 83136001L
)
observation <- list(
  bundle_class = "stage6bb_micro_observation_only_v1",
  Y = generated$Y,
  time_obs = generated$time_obs,
  Z = generated$Z,
  dimensions = list(S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
                    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
                    M_s = list(2L, 2L), K = 3L)
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
    seed = 83136101L,
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

s6b_assert(is.list(fit) && fit$n_cpus_used == 1L,
           "Micro fit did not complete on one CPU.")
s6b_assert(identical(fit$pre_score_interface$function_initialization,
                     "gram_unit_energy"),
           "Micro fit did not retain the G12 initialization provenance.")

score_all <- s6b_score_holdout(fit, split$holdout_plan, "all")
score_ppi <- s6b_score_holdout(fit, split$holdout_plan, "ppi_0.5")
aggregate <- rbind(score_all$aggregate, score_ppi$aggregate)
s6b_assert(nrow(score_all$subject_time) == 8L &&
             all(score_all$subject_time$variables == 20L) &&
             all(is.finite(aggregate$heldout_rmse)) &&
             all(is.finite(aggregate$heldout_nrmse)) &&
             all(is.finite(aggregate$heldout_mae)) &&
             all(is.finite(aggregate$residual_only_nlpd)) &&
             all(!aggregate$primary_dimension_score_ready) &&
             all(!aggregate$truth_used),
           "Micro held-out prediction contract failed.")

s6b_write_csv(aggregate, file.path(output_root, "SMOKE_HOLDOUT_SUMMARY.csv"))
s6b_write_csv(score_all$subject_time,
              file.path(output_root, "SMOKE_SUBJECT_TIME_ALL.csv"))
s6b_write_csv(data.frame(
  warning_index = seq_along(fit_warnings),
  warning = fit_warnings,
  stringsAsFactors = FALSE
), file.path(output_root, "SMOKE_WARNINGS.csv"))
s6b_write_csv(data.frame(
  check_id = c(
    "fit_completed", "one_cpu", "g12_initialization",
    "eight_subject_time_vectors", "twenty_variables_per_vector",
    "finite_heldout_point_losses", "truth_free",
    "nlpd_not_primary_or_calibrated"
  ),
  passed = rep(TRUE, 8L),
  detail = c(
    "micro fit object returned", "n_cpus_used=1",
    "gram_unit_energy + one pre-score + annealing",
    "8", "20", "RMSE/NRMSE/MAE finite", "truth_used=FALSE",
    "residual-noise-only NLPD; primary_dimension_score_ready=FALSE"
  ),
  stringsAsFactors = FALSE
), file.path(output_root, "SMOKE_STATUS.csv"))
s6b_write_lines(c(
  "status=STAGE6B_B_HOLDOUT_MICRO_SMOKE_PASS",
  "formal_experiment=FALSE",
  "new_full_size_fit_started=FALSE",
  "truth_accessed=FALSE",
  "n_cpus=1",
  "subject_time_vectors=8",
  "variables_per_vector=20",
  paste0("warnings_captured=", length(fit_warnings)),
  "primary_dimension_score_ready=FALSE"
), file.path(output_root, "SMOKE_COMPLETE.txt"))

cat("STAGE6B_B_HOLDOUT_MICRO_SMOKE_PASS output_root=", output_root,
    "\n", sep = "")
