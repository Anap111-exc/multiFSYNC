#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

s6g_abort <- function(...) stop(paste0(...), call. = FALSE)
s6g_assert <- function(value, message) {
  if (!isTRUE(value)) s6g_abort(message)
  invisible(TRUE)
}
s6g_write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, na = "")
  invisible(path)
}

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
s6g_assert(length(script_argument) == 1L, "Run this smoke with Rscript.")
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(
    script_file, winslash = "/", mustWork = TRUE
  ))
}
arguments <- commandArgs(trailingOnly = TRUE)
package_argument <- arguments[startsWith(arguments, "--package-root=")]
s6g_assert(length(package_argument) <= 1L,
           "Supply at most one --package-root argument.")
repo_root <- if (length(package_argument)) {
  sub("--package-root=", "", package_argument[[1L]], fixed = TRUE)
} else {
  file.path(script_root, "..", "..", "..")
}
s6g_assert(file.exists(file.path(repo_root, "DESCRIPTION")) &&
             dir.exists(file.path(repo_root, "R")),
           "Cannot resolve the multiFSYNC source-tree root.")

output_argument <- arguments[startsWith(arguments, "--output-root=")]
s6g_assert(length(output_argument) == 1L, "--output-root is required.")
output_root <- sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
s6g_assert(!dir.exists(output_root) && !file.exists(output_root),
           "Smoke output root already exists; refusing to overwrite it.")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

stage6b_root <- file.path(
  repo_root, "experiments", "v0L-V",
  "g12-truth-free-dimension-selection-stage6b-20260901"
)
stage6f_root <- file.path(
  repo_root, "experiments", "v0L-V",
  "g12-joint-dimension-selection-stage6f-20260911"
)
source(file.path(stage6b_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = TRUE)
source(file.path(stage6f_root,
                 "stage6f_joint_selection_tools_20260911_v1.R"),
       local = TRUE)

library_argument <- arguments[startsWith(arguments, "--library-root=")]
s6g_assert(length(library_argument) <= 1L,
           "Supply at most one --library-root argument.")
if (length(library_argument)) {
  library_root <- sub(
    "--library-root=", "", library_argument[[1L]], fixed = TRUE
  )
  s6g_assert(dir.exists(library_root),
             "--library-root does not identify an installed R library.")
  .libPaths(c(library_root, .libPaths()))
  s6g_assert(requireNamespace("multiFSYNC", quietly = TRUE),
             "multiFSYNC is not installed in --library-root.")
  package_loading_mode <- "preinstalled_temporary_library"
} else {
  s6g_assert(requireNamespace("pkgload", quietly = TRUE),
             "Package 'pkgload' is required for a source-tree smoke.")
  pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE)
  package_loading_mode <- "pkgload_source_tree"
}
s6g_assert(requireNamespace("Splinets", quietly = TRUE),
           "Package 'Splinets' is required for Jaoua-style truth generation.")

scenarios <- utils::read.csv(
  file.path(script_root, "SCENARIOS_STAGE6G_V1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
configs <- utils::read.csv(
  file.path(script_root, "FIT_CONFIGS_STAGE6G_V1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
seeds <- utils::read.csv(
  file.path(script_root, "START_SEEDS_STAGE6G_V1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
s6g_assert(nrow(scenarios) == 3L && nrow(configs) == 5L && nrow(seeds) == 4L,
           "Stage 6G registration tables have unexpected row counts.")

parse_integer_vector <- function(value, expected_length, label) {
  pieces <- strsplit(as.character(value), ";", fixed = TRUE)[[1L]]
  result <- suppressWarnings(as.integer(pieces))
  s6g_assert(length(result) == expected_length && all(is.finite(result)) &&
               all(result >= 1L), paste0("Invalid ", label, "."))
  result
}

started <- proc.time()[["elapsed"]]
generator_rows <- vector("list", nrow(scenarios))
fit_observation <- NULL

for (index in seq_len(nrow(scenarios))) {
  row <- scenarios[index, , drop = FALSE]
  L_f <- as.integer(row$truth_L_f[[1L]])
  L_s <- c(as.integer(row$truth_L_s_1[[1L]]),
           as.integer(row$truth_L_s_2[[1L]]))
  M <- as.integer(row$truth_M_common[[1L]])
  generated <- multiFSYNC::simulate_multi_study_structured(
    S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
    L_f = L_f, L_s = L_s,
    M_f = rep(M, L_f),
    M_s = list(rep(M, L_s[[1L]]), rep(M, L_s[[2L]])),
    K = 3L, n_obs = 6L, common_grid = FALSE,
    n_obs_range = c(6L, 7L), sigma_eps = 0.3,
    bool_sparse_loadings = TRUE, prop_sparse = 0.8,
    score_var_decay = TRUE, n_dense = 51L,
    mean_amp = 0.6, beta_amp = 0,
    bs_degree = c(2L, 3L), identified_loadings = TRUE,
    sparsity_mode = "fixed", seed = as.integer(row$data_seed[[1L]]),
    truth_basis_method = "jaoua_splinets",
    truth_basis_seed = as.integer(row$truth_basis_seed[[1L]])
  )
  truth <- generated$true_params
  metadata <- c(
    truth$phi_basis_meta,
    unlist(truth$psi_basis_meta, recursive = FALSE)
  )
  continuous_error <- max(vapply(metadata, function(value) {
    max(abs(value$gram - diag(nrow(value$gram))))
  }, numeric(1L)))
  response_error <- 0
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      observed <- do.call(cbind, generated$Y[[study]][[subject]])
      response_error <- max(
        response_error,
        abs(observed - truth$signal_true_values[[study]][[subject]] -
              truth$noise_true_values[[study]][[subject]])
      )
    }
  }
  n_time <- unlist(truth$n_obs_by_subject, use.names = FALSE)
  expected_blocks <- L_f + sum(L_s)
  passed <- identical(truth$truth_basis_method, "jaoua_splinets") &&
    length(metadata) == expected_blocks &&
    continuous_error < 1e-12 && response_error < 1e-12 &&
    all(n_time >= 6L & n_time <= 7L)
  generator_rows[[index]] <- data.frame(
    scenario_id = as.character(row$scenario_id[[1L]]),
    truth_L_f = L_f, truth_L_s = paste(L_s, collapse = ";"),
    truth_M = M, total_factor_blocks = length(metadata),
    continuous_gram_max_abs_error = continuous_error,
    response_identity_max_abs_error = response_error,
    min_times = min(n_time), max_times = max(n_time),
    passed = passed, stringsAsFactors = FALSE
  )
  s6g_assert(passed, paste0(
    "Generator smoke failed for ", row$scenario_id[[1L]], "."
  ))

  observation <- list(
    bundle_class = "stage6g_micro_observation_only_v1",
    Y = generated$Y, time_obs = generated$time_obs, Z = generated$Z,
    dimensions = list(
      S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
      L_f = L_f, L_s = L_s, M_f = rep(M, L_f),
      M_s = list(rep(M, L_s[[1L]]), rep(M, L_s[[2L]])), K = 3L
    )
  )
  s6g_assert(!length(s6b_forbidden_observation_names(observation)),
             "Micro observation bundle leaked a truth field.")
  if (index == 1L) fit_observation <- observation
  rm(generated, truth, metadata, observation)
  invisible(gc(FALSE))
}

middle <- s6f_make_time_holdout(fit_observation, "middle")
upper <- s6f_make_time_holdout(fit_observation, "upper")
middle_index <- vapply(middle$holdout_plan$records, `[[`, integer(1L),
                       "original_index")
upper_index <- vapply(upper$holdout_plan$records, `[[`, integer(1L),
                      "original_index")
s6g_assert(
  length(middle_index) == 8L && length(upper_index) == 8L &&
    all(upper_index == middle_index + 1L) &&
    all(vapply(middle$holdout_plan$records, `[[`, integer(1L),
               "train_time_count") >= 5L) &&
    all(vapply(upper$holdout_plan$records, `[[`, integer(1L),
               "train_time_count") >= 5L) &&
    identical(middle$holdout_plan$truth_used, FALSE) &&
    identical(upper$holdout_plan$truth_used, FALSE),
  "The two registered Stage 6G holdout splits failed their contract."
)
rm(fit_observation)
invisible(gc(FALSE))

registered_splits <- list(middle = middle, upper = upper)
fit_rows <- vector("list", nrow(configs) * length(registered_splits))
loss_rows <- vector("list", length(fit_rows))
warning_rows <- list()
fit_index <- 0L
for (split_name in names(registered_splits)) {
  split <- registered_splits[[split_name]]
  for (index in seq_len(nrow(configs))) {
    fit_index <- fit_index + 1L
    row <- configs[index, , drop = FALSE]
    L_f <- as.integer(row$fit_L_f[[1L]])
    L_s <- c(as.integer(row$fit_L_s_1[[1L]]),
             as.integer(row$fit_L_s_2[[1L]]))
    M_f <- parse_integer_vector(row$fit_M_f[[1L]], L_f, "fit_M_f")
    M_s <- list(
      parse_integer_vector(row$fit_M_s_1[[1L]], L_s[[1L]], "fit_M_s_1"),
      parse_integer_vector(row$fit_M_s_2[[1L]], L_s[[2L]], "fit_M_s_2")
    )
    fit_warnings <- character()
    fit <- withCallingHandlers(
      multiFSYNC::bayesSYNC_multi_pre_score(
        Y = split$train_observation$Y,
        Z = split$train_observation$Z,
        time_obs = split$train_observation$time_obs,
        L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s,
        K = 3L, anneal = c(1, 1.2, 3), list_hyper = NULL,
        n_g = 31L, time_g = NULL, tol_abs = 0, tol_rel = 0,
        maxit = 6L, n_cpus = 1L, verbose = FALSE,
        seed = as.integer(seeds$fit_seed[[1L]]),
        bool_scale = FALSE, bool_var_spec_prob = FALSE, d_0 = 20L,
        convergence_rule = "parameters", lambda_orth = 0,
        practical_control = NULL, initialization = "random",
        initialization_control = list(
          grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
        ),
        continuation_state = NULL, pre_score_sweeps = 1L,
        trace_sweeps = 1L,
        function_initialization = "gram_unit_energy"
      ),
      warning = function(condition) {
        fit_warnings <<- c(fit_warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    returned_L_s <- as.integer(fit$L_s_by_study)
    if (length(returned_L_s) == 1L) returned_L_s <- rep(returned_L_s, 2L)
    calibration <- fit$random_scale_calibration_diagnostics
    expected_calibrations <- sum(M_f) +
      sum(vapply(M_s, sum, integer(1L)))
    score <- s6b_score_holdout(fit, split$holdout_plan, "all")
    passed <- is.list(fit) && fit$n_cpus_used == 1L &&
      identical(fit$pre_score_interface$function_initialization,
                "gram_unit_energy") &&
      fit$pre_score_interface$pre_score_sweeps == 1L &&
      fit$L_f == L_f && identical(returned_L_s, L_s) &&
      identical(as.integer(fit$M_f), M_f) &&
      identical(lapply(fit$M_s, as.integer), M_s) &&
      length(fit$ELBO) > 0L && all(is.finite(fit$ELBO)) &&
      is.data.frame(calibration) &&
      nrow(calibration) == expected_calibrations &&
      all(calibration$success) && max(abs(calibration$after - 1)) < 1e-10 &&
      all(is.finite(score$aggregate$heldout_rmse)) &&
      all(is.finite(score$aggregate$heldout_nrmse)) &&
      identical(score$aggregate$truth_used, FALSE)
    fit_rows[[fit_index]] <- data.frame(
      split_id = as.character(split$holdout_plan$split_id),
      fit_config_id = as.character(row$fit_config_id[[1L]]),
      fit_L_f = L_f, fit_L_s = paste(L_s, collapse = ";"),
      fit_M_f = paste(M_f, collapse = ";"),
      fit_M_s = paste(vapply(M_s, paste, character(1L), collapse = ";"),
                      collapse = "|"),
      warnings_captured = length(fit_warnings),
      final_elbo = tail(fit$ELBO, 1L),
      heldout_rmse = score$aggregate$heldout_rmse,
      heldout_nrmse = score$aggregate$heldout_nrmse,
      passed = passed, stringsAsFactors = FALSE
    )
    score_rows <- score$subject_time
    score_rows$fit_config_id <- as.character(row$fit_config_id[[1L]])
    score_rows$split_id <- as.character(split$holdout_plan$split_id)
    loss_rows[[fit_index]] <- score_rows
    if (length(fit_warnings)) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        split_id = as.character(split$holdout_plan$split_id),
        fit_config_id = as.character(row$fit_config_id[[1L]]),
        warning_index = seq_along(fit_warnings), warning = fit_warnings,
        stringsAsFactors = FALSE
      )
    }
    s6g_assert(passed, paste0(
      "G12 fit/holdout smoke failed for ", row$fit_config_id[[1L]],
      " on ", split_name, "."
    ))
    rm(fit, score, calibration, score_rows)
    invisible(gc(FALSE))
  }
}

generator_checks <- do.call(rbind, generator_rows)
fit_checks <- do.call(rbind, fit_rows)
all_losses <- do.call(rbind, loss_rows)
subject_losses <- s6f_subject_cluster_losses(
  all_losses, "mse", expected_splits = 2L
)
selection_summary <- do.call(rbind, lapply(
  as.character(configs$fit_config_id), function(config_id) {
    rows <- subject_losses[
      subject_losses$fit_config_id == config_id, , drop = FALSE
    ]
    metadata <- configs[configs$fit_config_id == config_id, , drop = FALSE]
    s6g_assert(nrow(rows) == 8L, paste0(
      "Paired subject losses are incomplete for ", config_id, "."
    ))
    data.frame(
      fit_config_id = config_id,
      complexity_proxy = as.integer(metadata$complexity_proxy[[1L]]),
      subjects = nrow(rows), paired_splits_per_subject = 2L,
      mean_mse = mean(rows$mse),
      stratified_subject_se = s6f_stratified_se(rows$mse, rows$study),
      stringsAsFactors = FALSE
    )
  }
))
minimum_order <- order(
  selection_summary$mean_mse, selection_summary$complexity_proxy,
  selection_summary$fit_config_id
)
minimum_index <- minimum_order[[1L]]
one_se_threshold <- selection_summary$mean_mse[[minimum_index]] +
  selection_summary$stratified_subject_se[[minimum_index]]
selection_summary$within_one_se <-
  selection_summary$mean_mse <= one_se_threshold
selected_order <- order(
  selection_summary$complexity_proxy[selection_summary$within_one_se],
  selection_summary$mean_mse[selection_summary$within_one_se],
  selection_summary$fit_config_id[selection_summary$within_one_se]
)
selected_index <- which(selection_summary$within_one_se)[selected_order[[1L]]]
selection_summary$minimum_loss <- seq_len(nrow(selection_summary)) ==
  minimum_index
selection_summary$one_se_selected <- seq_len(nrow(selection_summary)) ==
  selected_index
s6g_assert(
  sum(selection_summary$minimum_loss) == 1L &&
    sum(selection_summary$one_se_selected) == 1L &&
    all(is.finite(selection_summary$mean_mse)) &&
    all(is.finite(selection_summary$stratified_subject_se)),
  "Five-configuration paired one-SE selector failed its micro contract."
)
warnings <- if (length(warning_rows)) {
  do.call(rbind, warning_rows)
} else {
  data.frame(split_id = character(), fit_config_id = character(),
             warning_index = integer(), warning = character(),
             stringsAsFactors = FALSE)
}
elapsed <- proc.time()[["elapsed"]] - started
s6g_write_csv(generator_checks,
              file.path(output_root, "GENERATOR_CHECKS.csv"))
s6g_write_csv(fit_checks, file.path(output_root, "FIT_CHECKS.csv"))
s6g_write_csv(selection_summary,
              file.path(output_root, "ONE_SE_SELECTION_CHECK.csv"))
s6g_write_csv(warnings, file.path(output_root, "WARNINGS.csv"))
writeLines(c(
  "status=STAGE6G_LOCAL_MICRO_SMOKE_PASS",
  "formal_experiment=FALSE",
  "data_scenarios_generated=3",
  "micro_fits=10",
  "registered_holdout_splits_checked=2",
  "five_configuration_one_se_selection=PASS",
  paste0("micro_one_se_selected_config=",
         selection_summary$fit_config_id[[selected_index]]),
  "micro_selection_has_no_scientific_interpretation=TRUE",
  "truth_used_for_fit_or_holdout_scoring=FALSE",
  "n_cpus_per_fit=1",
  paste0("package_loading_mode=", package_loading_mode),
  paste0("warnings_captured=", nrow(warnings)),
  paste0("elapsed_seconds=", format(elapsed, digits = 8))
), file.path(output_root, "SMOKE_COMPLETE.txt"), useBytes = TRUE)

cat("STAGE6G_LOCAL_MICRO_SMOKE_PASS scenarios=3 fits=10 splits=2 ",
    "truth_used_for_fit_or_scoring=FALSE elapsed_seconds=",
    format(elapsed, digits = 8), "\n", sep = "")
