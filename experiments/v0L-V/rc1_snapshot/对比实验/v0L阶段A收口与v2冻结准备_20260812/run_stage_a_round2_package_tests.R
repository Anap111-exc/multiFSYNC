#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_stage_a_round2_package_tests.R <project-root> <ascii-lib>")
}
project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
ascii_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
local_test_lib <- file.path(project_root, "_r_test_lib")
.libPaths(c(
  ascii_lib,
  if (dir.exists(local_test_lib)) {
    normalizePath(local_test_lib, winslash = "/")
  } else character(),
  .libPaths()
))
library(multiFSYNC)
test_dir <- file.path(ascii_lib, "multiFSYNC", "tests", "testthat")
testthat::test_dir(
  test_dir,
  package = "multiFSYNC",
  load_package = "installed",
  reporter = "summary",
  stop_on_failure = TRUE
)
cat("STAGE_A_ROUND2_FULL_PACKAGE_TESTS_PASS\n")
