#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(sub("^--file=", "", file_argument[[1L]]))
if (!dir.exists(root)) {
  root <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]),
                                winslash = "/", mustWork = TRUE))
}
stage6b_root <- file.path(dirname(root),
                          "g12-truth-free-dimension-selection-stage6b-20260901")
source(file.path(stage6b_root, "stage6b_truth_free_tools_20260901_v1.R"))
source(file.path(stage6b_root,
                 "stage6b_predictive_pilot_tools_20260902_v1.R"))
source(file.path(root, "stage6d_tools_20260906_v1.R"))

tests <- list()
run_test <- function(id, expression) {
  started <- Sys.time()
  error <- tryCatch({ force(expression); NULL }, error = identity)
  tests[[length(tests) + 1L]] <<- data.frame(
    test_id = id, passed = is.null(error),
    detail = if (is.null(error)) "PASS" else conditionMessage(error),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE
  )
}

data <- s6d_read_csv(file.path(root, "DATA_SEEDS_STAGE6D_V1.csv"))
configs <- s6d_read_csv(file.path(root, "FIT_CONFIGS_STAGE6D_V1.csv"))
seeds <- s6d_read_csv(file.path(root, "START_SEEDS_STAGE6D_V1.csv"))

run_test("registration", s6d_validate_registration(data, configs, seeds))
manifest <- s6d_build_fit_manifest(data, configs, seeds)
run_test("manifest_180", stopifnot(nrow(manifest) == 180L))
run_test("strata_15_by_12", stopifnot(
  length(unique(manifest$selection_stratum_id)) == 15L,
  all(table(manifest$selection_stratum_id) == 12L)
))
run_test("paired_seeds", stopifnot(
  all(vapply(split(manifest$fit_seed, manifest$seed_index),
             function(value) length(unique(value)) == 1L, logical(1L)))
))
run_test("truth_isolation", stopifnot(
  !any(manifest$truth_available_to_fit),
  !any(manifest$truth_available_to_stopping),
  !any(manifest$truth_available_to_start_selection),
  !any(manifest$truth_available_to_dimension_selection)
))
run_test("per_fit_single_cpu", stopifnot(all(manifest$n_cpus == 1L)))
run_test("factor_axis", stopifnot(
  identical(s6d_axis_metadata(configs, "factor")$complexity_value, 1:3)
))
run_test("fpca_axis", stopifnot(
  identical(s6d_axis_metadata(configs, "fpca")$complexity_value, 2:4)
))

make_losses <- function(config_id, mse, data_id = "g12s6d_01") {
  data.frame(
    study = rep(1:2, each = 30L), subject = rep(1:30, 2L),
    mse = rep(mse, 60L), mae = rep(sqrt(mse), 60L),
    observed_second_moment = rep(1, 60L), data_id = data_id,
    fit_config_id = config_id, winner_fit_id = paste0(config_id, "__winner"),
    stringsAsFactors = FALSE
  )
}
losses <- do.call(rbind, list(
  make_losses("factor_L1_M2", 0.10),
  make_losses("factor_L2_M2", 0.11),
  make_losses("factor_L3_M2", 0.12),
  make_losses("fpca_L1_M3", 0.115),
  make_losses("fpca_L1_M4", 0.13)
))
axes <- s6d_select_axes(losses, configs, "g12s6d_01")
run_test("factor_one_se", stopifnot(
  identical(axes$factor$selected_config_id, "factor_L1_M2")
))
run_test("fpca_one_se", stopifnot(
  identical(axes$fpca$selected_config_id, "factor_L1_M2")
))

empty_warnings <- data.frame(
  class = character(), message = character(), call = character(),
  phase = character(), stringsAsFactors = FALSE
)
run_test("zero_warning_rows", stopifnot(
  nrow(s6d_set_warning_phase(empty_warnings, "stage6d")) == 0L
))

results <- do.call(rbind, tests)
output <- commandArgs(trailingOnly = TRUE)
output <- sub("^--output=", "", output[grepl("^--output=", output)])
if (length(output) == 1L) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, output, row.names = FALSE, quote = TRUE)
}
if (!all(results$passed)) {
  print(results[!results$passed, , drop = FALSE])
  stop("Stage 6D targeted tests failed.", call. = FALSE)
}
cat("STAGE6D_TEST_PASS ", sum(results$passed), "/", nrow(results), "\n",
    sep = "")
