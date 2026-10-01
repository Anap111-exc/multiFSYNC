#!/usr/bin/env Rscript

# Read-only binding check against a previously evaluated Stage 7C winner.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript SCRIPT STAGE7C_ROOT SUPPORT_REPO")
stage7c_root <- normalizePath(args[[1L]], mustWork = TRUE)
support_v0lv <- file.path(normalizePath(args[[2L]], mustWork = TRUE),
                           "experiments", "v0L-V")
source(file.path(support_v0lv, "g12-dimension-overspec-stage6a-20260831",
                 "stage6a_common_20260831_v1.R"), local = FALSE)
source(file.path(support_v0lv,
                 "g12-reachability-causal-audit-stage6h-20260913",
                 "stage6h_evaluation_helpers_20260913_v1.R"), local = FALSE)
PP_SNAPSHOT <- file.path(support_v0lv, "rc1_snapshot")
pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
stage7c_evaluator_path <- file.path(
  "/root/multiFSYNC-stage7c-20260927-v1", "experiments", "v0L-V",
  "covariate-routeb-stage7c-independent-confirmation-20260927",
  "evaluate_covariate_routeb_stage7c_20260928_v1.R")
expected_hash <- "1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8"
stopifnot(identical(tolower(pp_sha256(stage7c_evaluator_path)), expected_hash))
stage7c_lines <- readLines(stage7c_evaluator_path, warn = FALSE)
override_start <- which(grepl("^mean_design <- function", stage7c_lines))
override_end <- which(grepl(
  "^install_dense_residual_override\\(evaluator_environment\\)",
  stage7c_lines))
stopifnot(length(override_start) == 1L, length(override_end) == 1L,
          override_start < override_end)
eval(parse(text = paste(stage7c_lines[override_start:override_end],
                        collapse = "\n")), envir = environment())
data_id <- "covrb7c_03__C1"
winners <- read.csv(file.path(stage7c_root, "truth_free_selection",
                              "TRUTH_FREE_WINNERS.csv"))
fit_id <- winners$fit_id[winners$data_id == data_id]
stopifnot(length(fit_id) == 1L)
fit <- readRDS(file.path(stage7c_root, "fits", fit_id, "fit.rds"))
truth <- readRDS(file.path(stage7c_root, "data", data_id,
                           "sealed_truth", "truth_bundle.rds"))
result <- s6h_evaluate_fit(fit, truth, evaluator_environment)
reference <- read.csv(file.path(stage7c_root,
  "scientific_evaluation_20260928_v1", "ALL_START_SCIENCE.csv"))
reference <- reference[reference$fit_id == fit_id, ]
actual <- result$summary
fields <- c("selected_shared", "selected_specific_study_1",
            "selected_specific_study_2", "global_correct",
            "global_misplaced", "global_missing", "global_extra",
            "global_duplicate", "study_correct", "study_misplaced",
            "study_missing", "study_candidate_extra",
            "study_candidate_duplicate", "truth_global_direction_count",
            "truth_study_role_count")
stopifnot(nrow(actual) == 1L, nrow(reference) == 1L,
          all(fields %in% names(actual)),
          all(vapply(fields, function(field) {
            identical(as.numeric(actual[[field]]),
                      as.numeric(reference[[field]]))
          }, logical(1L))))
cat("STAGE7E_FROZEN_EVALUATOR_BINDING_PASS fit_id=", fit_id, "\n", sep = "")
