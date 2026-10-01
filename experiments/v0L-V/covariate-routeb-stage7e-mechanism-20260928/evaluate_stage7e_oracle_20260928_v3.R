#!/usr/bin/env Rscript

# Stage 7E scientific diagnosis only. These oracle-initialized endpoints never
# enter Stage 7C selection or a formal performance denominator.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript evaluate_stage7e_oracle_20260928_v3.R STAGE7C_ROOT STAGE7E_ROOT SUPPORT_REPO")
}
stage7c_root <- normalizePath(args[[1L]], mustWork = TRUE)
stage7e_root <- normalizePath(args[[2L]], mustWork = TRUE)
support_repo <- normalizePath(args[[3L]], mustWork = TRUE)
evaluation_root <- file.path(stage7e_root, "evaluation_v2")
if (file.exists(evaluation_root) || dir.exists(evaluation_root)) {
  stop("Refusing to overwrite evaluation_v2.")
}

hash_file <- function(path) {
  if (!file.exists(path)) stop("Missing bound file: ", path)
  result <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(result[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
check_hash <- function(path, expected) {
  if (!identical(tolower(hash_file(path)), tolower(expected))) {
    stop("Bound file hash changed: ", path)
  }
}
write_csv_new <- function(x, name) {
  path <- file.path(evaluation_root, name)
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
write_lines_new <- function(x, name) {
  path <- file.path(evaluation_root, name)
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  writeLines(x, path, useBytes = TRUE)
}
abs_cosine <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  abs(sum(x * y) / denominator)
}

check_hash("/root/multiFSYNC-stage7c-20260927-v1.tar.gz",
           "93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924")
if (as.character(utils::packageVersion("multiFSYNC")) != "0.3.0.9000" ||
    normalizePath(system.file(package = "multiFSYNC")) !=
      "/root/multiFSYNC-stage7c-lib-20260927-v1/multiFSYNC") {
  stop("Stage 7C isolated package binding failed.")
}
support_v0lv <- file.path(support_repo, "experiments", "v0L-V")
stage6a_path <- file.path(support_v0lv,
  "g12-dimension-overspec-stage6a-20260831", "stage6a_common_20260831_v1.R")
stage6h_path <- file.path(support_v0lv,
  "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R")
check_hash(stage6a_path,
           "bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6")
check_hash(stage6h_path,
           "5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221")
source(stage6a_path, local = FALSE)
source(stage6h_path, local = FALSE)
PP_SNAPSHOT <- file.path(support_v0lv, "rc1_snapshot")
if (!dir.exists(PP_SNAPSHOT)) stop("Frozen evaluator snapshot is missing.")
pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
stage7c_evaluator_path <- file.path(
  "/root/multiFSYNC-stage7c-20260927-v1", "experiments", "v0L-V",
  "covariate-routeb-stage7c-independent-confirmation-20260927",
  "evaluate_covariate_routeb_stage7c_20260928_v1.R")
check_hash(stage7c_evaluator_path,
  "1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8")
stage7c_lines <- readLines(stage7c_evaluator_path, warn = FALSE)
override_start <- which(grepl("^mean_design <- function", stage7c_lines))
override_end <- which(grepl(
  "^install_dense_residual_override\\(evaluator_environment\\)",
  stage7c_lines))
if (length(override_start) != 1L || length(override_end) != 1L ||
    override_start >= override_end) stop("Stage 7C override markers changed.")
eval(parse(text = paste(stage7c_lines[override_start:override_end],
                        collapse = "\n")), envir = environment())

manifest <- read.csv(file.path(stage7e_root, "FIT_MANIFEST.csv"))
binding <- read.csv(file.path(stage7e_root, "INPUT_BINDING.csv"))
expected_ids <- unlist(lapply(1:3, function(index) {
  sprintf("s7e_%02d_C1_%s", index, c("exact", "near"))
}), use.names = FALSE)
if (!identical(manifest$fit_id, expected_ids) || nrow(binding) != 3L ||
    !file.exists(file.path(stage7e_root, "FITS_COMPLETE.txt")) ||
    !file.exists(file.path(stage7e_root, "MICRO_SMOKE_PASS.txt"))) {
  stop("Six-fit registration or completion contract failed.")
}
for (index in seq_len(nrow(binding))) {
  row <- binding[index, ]
  observation <- file.path(stage7c_root, "data", row$data_id,
                           "observation_bundle.rds")
  truth <- file.path(stage7c_root, "data", row$data_id,
                     "sealed_truth", "truth_bundle.rds")
  check_hash(observation, row$observation_sha256)
  check_hash(truth, row$truth_sha256)
}
for (fit_id in manifest$fit_id) {
  directory <- file.path(stage7e_root, "fits", fit_id)
  if (!file.exists(file.path(directory, "FIT_COMPLETE.txt")) ||
      file.exists(file.path(directory, "FIT_ERROR.txt")) ||
      !file.exists(file.path(directory, "fit.rds"))) {
    stop("Incomplete diagnostic fit: ", fit_id)
  }
}

winners <- read.csv(file.path(stage7c_root, "truth_free_selection",
                              "TRUTH_FREE_WINNERS.csv"))
stage7c_science <- read.csv(file.path(stage7c_root,
  "scientific_evaluation_20260928_v1", "ALL_START_SCIENCE.csv"))
reference <- stage7c_science[stage7c_science$fit_id %in%
  winners$fit_id & stage7c_science$data_id %in% binding$data_id, ]
if (nrow(reference) != 3L || !setequal(reference$data_id, binding$data_id)) {
  stop("Stage 7C frozen winner/science reference mismatch.")
}

factor_rows <- list()
endpoint_rows <- list()
trace_rows <- list()
science_rows <- list()
for (index in seq_len(nrow(manifest))) {
  row <- manifest[index, ]
  fit_id <- row$fit_id
  fit <- readRDS(file.path(stage7e_root, "fits", fit_id, "fit.rds"))
  if (length(fit$ELBO) != 400L || any(!is.finite(fit$ELBO)) ||
      fit$t1_sweeps != 400L || fit$annealing_sweeps != 0L) {
    stop("Diagnostic objective or sweep contract failed: ", fit_id)
  }
  truth_bundle <- readRDS(file.path(stage7c_root, "data", row$data_id,
                                    "sealed_truth", "truth_bundle.rds"))
  if (!identical(truth_bundle$data_id, row$data_id) ||
      !identical(truth_bundle$bundle_class,
                 "covariate_routeb_pilot_sealed_truth")) {
    stop("Truth bundle identity failed: ", fit_id)
  }
  truth <- truth_bundle$data$true_params
  frozen <- s6h_evaluate_fit(fit, truth_bundle, evaluator_environment)
  summary <- frozen$summary
  if (nrow(summary) != 1L) stop("Frozen evaluator returned multiple rows.")
  expected_columns <- c("selected_shared", "selected_specific_study_1",
    "selected_specific_study_2", "global_correct", "global_misplaced",
    "global_missing", "global_extra", "global_duplicate",
    "study_correct", "study_misplaced", "study_missing",
    "study_candidate_extra", "study_candidate_duplicate",
    "truth_global_direction_count", "truth_study_role_count")
  if (!all(expected_columns %in% names(summary))) {
    stop("Frozen evaluator output lacks Stage 7C recovery fields.")
  }
  exact_counts <- summary$selected_shared == 1L &&
    summary$selected_specific_study_1 == 1L &&
    summary$selected_specific_study_2 == 1L
  core <- summary$global_correct == summary$truth_global_direction_count &&
    summary$global_misplaced == 0L && summary$global_missing == 0L &&
    summary$study_correct == summary$truth_study_role_count &&
    summary$study_misplaced == 0L && summary$study_missing == 0L
  exact <- core && exact_counts && summary$global_extra == 0L &&
    summary$global_duplicate == 0L &&
    summary$study_candidate_extra == 0L &&
    summary$study_candidate_duplicate == 0L

  names_block <- c("shared", "specific_1", "specific_2")
  raw <- list(fit$mu_q_a[, 1L],
              fit$mu_q_b_specific[[1L]][, 1L],
              fit$mu_q_b_specific[[2L]][, 1L])
  canonical <- list(fit$mu_q_a_hat[, 1L],
                    fit$mu_q_b_specific_hat[[1L]][, 1L],
                    fit$mu_q_b_specific_hat[[2L]][, 1L])
  target <- list(truth$a_true[, 1L], truth$b_true[[1L]][, 1L],
                 truth$b_true[[2L]][, 1L])
  gamma <- list(fit$mu_q_gamma_a[, 1L],
                fit$mu_q_gamma_b[[1L]][, 1L],
                fit$mu_q_gamma_b[[2L]][, 1L])
  selected <- c(fit$factor_ppi_shared[[1L]],
                fit$factor_ppi_specific[[1L]][[1L]],
                fit$factor_ppi_specific[[2L]][[1L]]) >= 0.5
  cosines <- numeric(3L)
  for (block in seq_len(3L)) {
    true_support <- abs(target[[block]]) > 0
    estimated_support <- gamma[[block]] >= 0.5
    cosines[[block]] <- abs_cosine(canonical[[block]], target[[block]])
    factor_rows[[length(factor_rows) + 1L]] <- data.frame(
      fit_id = fit_id, data_id = row$data_id, mode = row$mode,
      block = names_block[[block]], factor_selected = selected[[block]],
      factor_ppi = c(fit$factor_ppi_shared[[1L]],
        fit$factor_ppi_specific[[1L]][[1L]],
        fit$factor_ppi_specific[[2L]][[1L]])[[block]],
      raw_loading_abs_cosine = abs_cosine(raw[[block]], target[[block]]),
      canonical_loading_abs_cosine = cosines[[block]],
      true_support = sum(true_support),
      estimated_support = sum(estimated_support),
      support_true_positive = sum(true_support & estimated_support),
      support_false_positive = sum(!true_support & estimated_support),
      loading_ppi_mean = mean(gamma[[block]]),
      stringsAsFactors = FALSE)
  }
  stage7c <- reference[reference$data_id == row$data_id, ]
  if (nrow(stage7c) != 1L) stop("Stage 7C reference row missing.")
  winner <- winners[winners$data_id == row$data_id, ]
  if (nrow(winner) != 1L || winner$fit_id != stage7c$fit_id) {
    stop("Stage 7C truth-free winner identity changed.")
  }
  objective <- tail(fit$ELBO, 1L)
  objective_eligible <- isTRUE(fit$elbo_objective_valid) &&
    identical(fit$elbo_role, "objective") &&
    identical(as.numeric(tail(fit$ELBO_temperature, 1L)), 1)
  endpoint_rows[[index]] <- data.frame(
    fit_id = fit_id, data_id = row$data_id, mode = row$mode,
    final_ordinary_elbo = objective,
    objective_eligible = objective_eligible,
    elbo_invalid_reasons = paste(fit$elbo_invalid_reasons,
                                 collapse = ";"),
    stage7c_winner_fit_id = winner$fit_id,
    stage7c_winner_ordinary_elbo = winner$final_ordinary_elbo,
    elbo_minus_stage7c_winner = if (objective_eligible) {
      objective - winner$final_ordinary_elbo
    } else NA_real_,
    stage7c_winner_core_truth_recovered = stage7c$core_truth_recovered,
    stage7c_winner_exact_structure =
      stage7c$exact_thresholded_factor_structure,
    selected_shared = summary$selected_shared,
    selected_specific_1 = summary$selected_specific_study_1,
    selected_specific_2 = summary$selected_specific_study_2,
    frozen_core_truth_recovered = core,
    frozen_exact_thresholded_structure = exact,
    all_three_factors_selected = all(selected),
    minimum_loading_cosine = min(cosines),
    diagnostic_three_direction_core = all(selected) &&
      all(is.finite(cosines) & cosines >= 0.9),
    strict_practical_converged = isTRUE(fit$practical_converged),
    ordinary_sweeps = fit$t1_sweeps,
    annealing_sweeps = fit$annealing_sweeps,
    truth_used_for_initialization = TRUE,
    truth_used_for_stopping_or_selection = FALSE,
    formal_candidate_pool = FALSE,
    stringsAsFactors = FALSE)
  science_rows[[index]] <- cbind(
    data.frame(fit_id = fit_id, data_id = row$data_id, mode = row$mode),
    summary, frozen$candidate_activity_summary)
  trace <- fit$driver_trace
  trace_rows[[index]] <- cbind(
    data.frame(fit_id = fit_id, data_id = row$data_id, mode = row$mode),
    trace)
  cat("STAGE7E_EVALUATED ", fit_id, "\n", sep = "")
}

dir.create(evaluation_root, recursive = FALSE)
write_csv_new(do.call(rbind, endpoint_rows), "ORACLE_ENDPOINT_COMPARISON.csv")
write_csv_new(do.call(rbind, factor_rows), "ORACLE_FACTOR_DIAGNOSTICS.csv")
write_csv_new(do.call(rbind, trace_rows), "ORACLE_DRIVER_TRACE.csv")
write_csv_new(do.call(rbind, science_rows), "ORACLE_FROZEN_EVALUATOR_SCIENCE.csv")
write_lines_new(c(
  "STAGE7E_EVALUATION_COMPLETE", "endpoints=6",
  "evaluator=frozen_Stage7C_s6h_plus_explicit_loading_diagnostics",
  "truth_used_for_initialization=TRUE",
  "truth_used_for_stopping_or_selection=FALSE",
  "formal_candidate_pool=FALSE",
  paste0("stage6a_support_sha256=", hash_file(stage6a_path)),
  paste0("stage6h_support_sha256=", hash_file(stage6h_path)),
  paste0("stage7c_evaluator_sha256=", hash_file(stage7c_evaluator_path))),
  "EVALUATION_COMPLETE.txt")
cat("STAGE7E_EVALUATION_COMPLETE endpoints=6\n")
