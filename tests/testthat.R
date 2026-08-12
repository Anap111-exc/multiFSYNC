if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(multiFSYNC)
  test_check("multiFSYNC")
} else {
  message("Skipping testthat suite because the suggested package 'testthat' is unavailable.")
}
