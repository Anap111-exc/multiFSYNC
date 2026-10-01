#!/usr/bin/env Rscript

# Targeted no-fit test for the Stage 7C evaluator/summarizer adaptation.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_argument) == 1L)
script_path <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_path)
evaluator_path <- file.path(
  script_root, "evaluate_covariate_routeb_stage7c_20260928_v1.R")
summarizer_path <- file.path(
  script_root, "summarize_covariate_routeb_stage7c_20260928_v1.R")
pilot_evaluation <- file.path(
  script_root, "..", "covariate-routeb-pilot-20260927",
  "scientific_evaluation_20260927_v1")

stopifnot(file.exists(evaluator_path), file.exists(summarizer_path),
          file.exists(file.path(pilot_evaluation, "ALL_START_SCIENCE.csv")),
          file.exists(file.path(pilot_evaluation, "WINNER_SCIENCE.csv")))
invisible(parse(file = evaluator_path))
invisible(parse(file = summarizer_path))

evaluator <- paste(readLines(evaluator_path, warn = FALSE), collapse = "\n")
checks <- c(
  stage7c_experiment_id = grepl(
    'experiment_id <- "COVARIATE_ROUTEB_STAGE7C_V1_20260927"',
    evaluator, fixed = TRUE),
  formal_108 = grepl("expected_fits <- if (smoke) 3L else 108L",
                     evaluator, fixed = TRUE),
  winner_9 = grepl('paste0("selected_winners=", length(expected_data_ids))',
                   evaluator, fixed = TRUE),
  three_logical_seeds = grepl('sprintf("covrb7c_%02d", 1:3)', evaluator,
                              fixed = TRUE),
  no_pilot_data_id = !grepl("covrbp_01", evaluator, fixed = TRUE),
  no_fit_or_continuation = grepl("new_fit_or_continuation_started=FALSE",
                                 evaluator, fixed = TRUE))
stopifnot(all(checks))

read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)
science_one <- read_csv(file.path(pilot_evaluation, "ALL_START_SCIENCE.csv"))
winners_one <- read_csv(file.path(pilot_evaluation, "WINNER_SCIENCE.csv"))

replicate_seed <- function(value, index) {
  logical_id <- sprintf("covrb7c_%02d", index)
  value$logical_data_id <- logical_id
  value$data_id <- paste0(logical_id, "__", value$beta_arm)
  value$fit_id <- paste0(logical_id, "__", value$fit_id)
  value
}
science <- do.call(rbind, lapply(1:3, function(index) {
  replicate_seed(science_one, index)
}))
winners <- do.call(rbind, lapply(1:3, function(index) {
  replicate_seed(winners_one, index)
}))
stopifnot(nrow(science) == 108L, nrow(winners) == 9L,
          all(table(science$data_id) == 12L))

test_root <- tempfile("stage7c_scientific_summary_test_")
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)
evaluation_root <- file.path(test_root, "scientific_evaluation_20260928_v1")
dir.create(evaluation_root, recursive = TRUE)
utils::write.csv(science, file.path(evaluation_root, "ALL_START_SCIENCE.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
utils::write.csv(winners, file.path(evaluation_root, "WINNER_SCIENCE.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
writeLines("status=SYNTHETIC_NO_FIT_TEST",
           file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")
status <- system2(
  rscript,
  c(shQuote(normalizePath(summarizer_path, winslash = "/")),
    shQuote(paste0("--output-root=",
                   normalizePath(test_root, winslash = "/")))),
  stdout = TRUE, stderr = TRUE)
stopifnot(is.null(attr(status, "status")) || attr(status, "status") == 0L)

summary_root <- file.path(
  evaluation_root, "descriptive_summary_20260928_v1")
pool <- read_csv(file.path(summary_root, "POOL_BASIN_AND_CONVERGENCE.csv"))
arm <- read_csv(file.path(summary_root, "ARM_BASIN_AND_CONVERGENCE.csv"))
gates <- read_csv(file.path(summary_root,
                            "PREREGISTERED_CONFIRMATION_GATES.csv"))
stopifnot(nrow(pool) == 9L, nrow(arm) == 3L, nrow(gates) == 5L,
          file.exists(file.path(summary_root,
                                "DESCRIPTIVE_SUMMARY_COMPLETE.txt")))

cat("STAGE7C_SCIENTIFIC_EVALUATION_TEST_PASS checks=", length(checks),
    " synthetic_fits=", nrow(science), " pools=", nrow(pool), "\n",
    sep = "")
