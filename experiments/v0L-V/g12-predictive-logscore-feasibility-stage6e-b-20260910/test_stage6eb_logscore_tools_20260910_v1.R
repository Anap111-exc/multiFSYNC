#!/usr/bin/env Rscript

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- command_line[startsWith(command_line, "--file=")]
script_file <- if (length(file_argument)) {
  sub("--file=", "", file_argument[[1L]], fixed = TRUE)
} else {
  "test_stage6eb_logscore_tools_20260910_v1.R"
}
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
source(file.path(script_root, "stage6eb_logscore_tools_20260910_v1.R"),
       local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
if (!length(output_argument)) {
  stop("--output-root is required.", call. = FALSE)
}
output_root <- sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Test output root already exists; refusing to overwrite it.",
       call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
write_lines <- function(value, path) {
  writeLines(enc2utf8(value), path, useBytes = TRUE)
}

test_results <- list()
run_test <- function(test_id, expression) {
  started <- proc.time()[["elapsed"]]
  failure <- NULL
  tryCatch(
    force(expression),
    error = function(condition) failure <<- conditionMessage(condition)
  )
  passed <- is.null(failure)
  test_results[[length(test_results) + 1L]] <<- data.frame(
    test_id = test_id,
    passed = passed,
    detail = if (passed) "PASS" else failure,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  )
  invisible(passed)
}

run_test("independent_gaussian_matches_closed_form", {
  observed <- c(-0.1, 0.5, 1.2, -0.8)
  mean <- c(0.2, 0.1, 0.9, -0.2)
  variance <- c(0.4, 0.7, 1.1, 0.6)
  result <- s6eb_lowrank_gaussian_nlpd(
    observed, mean, variance,
    matrix(numeric(), nrow = length(observed), ncol = 0L)
  )
  expected <- sum(
    0.5 * log(2 * pi * variance) +
      0.5 * (observed - mean)^2 / variance
  )
  stopifnot(abs(result$total_negative_log_density - expected) < 1e-12,
            result$diagonal_floor_count == 0L,
            result$covariance_rank_update == 0L)
})

run_test("woodbury_matches_dense_cholesky", {
  set.seed(910601L)
  variable_count <- 17L
  rank <- 4L
  observed <- stats::rnorm(variable_count)
  mean <- stats::rnorm(variable_count, sd = 0.3)
  diagonal <- stats::runif(variable_count, 0.2, 1.4)
  low_rank <- matrix(
    stats::rnorm(variable_count * rank, sd = 0.25),
    nrow = variable_count,
    ncol = rank
  )
  efficient <- s6eb_lowrank_gaussian_nlpd(
    observed, mean, diagonal, low_rank
  )
  dense <- s6eb_dense_gaussian_nlpd(
    observed, mean, diagonal, low_rank
  )
  stopifnot(
    abs(efficient$total_negative_log_density -
          dense$total_negative_log_density) < 1e-10,
    efficient$diagonal_floor_count == 0L,
    efficient$covariance_rank_update == rank
  )
})

run_test("bilinear_factor_moments_match_monte_carlo", {
  set.seed(910602L)
  design <- c(0.7, -0.2, 0.5)
  coefficient_mean <- matrix(
    c(0.4, -0.1, 0.2,
      -0.3, 0.6, 0.1),
    nrow = 3L,
    ncol = 2L
  )
  coefficient_covariances <- list(
    matrix(c(0.08, 0.01, 0.00,
             0.01, 0.05, 0.01,
             0.00, 0.01, 0.04), 3L, 3L),
    matrix(c(0.06, 0.00, 0.01,
             0.00, 0.04, 0.00,
             0.01, 0.00, 0.05), 3L, 3L)
  )
  score_mean <- c(0.35, -0.20)
  score_covariance <- matrix(c(0.30, 0.05, 0.05, 0.25), 2L, 2L)
  analytic <- s6eb_factor_process_moments(
    score_mean, score_covariance, coefficient_mean,
    coefficient_covariances, design
  )

  draw_count <- 160000L
  score_draw <- matrix(stats::rnorm(draw_count * 2L), ncol = 2L) %*%
    chol(score_covariance)
  score_draw <- sweep(score_draw, 2L, score_mean, "+")
  function_mean <- as.numeric(design %*% coefficient_mean)
  function_variance <- vapply(coefficient_covariances, function(covariance) {
    drop(crossprod(design, covariance %*% design))
  }, numeric(1L))
  function_draw <- matrix(stats::rnorm(draw_count * 2L), ncol = 2L)
  function_draw <- sweep(function_draw, 2L, sqrt(function_variance), "*")
  function_draw <- sweep(function_draw, 2L, function_mean, "+")
  process_draw <- rowSums(score_draw * function_draw)
  stopifnot(
    abs(mean(process_draw) - analytic$mean) < 0.006,
    abs(mean(process_draw^2) - analytic$second_moment) < 0.012,
    abs(stats::var(process_draw) - analytic$variance) < 0.012,
    analytic$numerical_clamp_count == 0L
  )
})

run_test("loading_times_process_covariance_identity", {
  loading_mean <- c(0.5, -0.2, 0.1)
  loading_variance <- c(0.08, 0.03, 0.05)
  process <- list(
    mean = -0.4,
    second_moment = 0.9,
    variance = 0.74,
    numerical_clamp_count = 0L
  )
  contribution <- s6eb_factor_contribution(
    loading_mean,
    loading_variance + loading_mean^2,
    process,
    "unit test factor"
  )
  reconstructed_covariance <- diag(contribution$diagonal) +
    tcrossprod(contribution$low_rank)
  expected_covariance <-
    process$variance * tcrossprod(loading_mean) +
    process$second_moment * diag(loading_variance)
  stopifnot(
    max(abs(reconstructed_covariance - expected_covariance)) < 1e-12,
    max(abs(contribution$mean - loading_mean * process$mean)) < 1e-12,
    contribution$loading_variance_clamp_count == 0L
  )
})

run_test("known_gaussian_marginal_coverage", {
  set.seed(910603L)
  variable_count <- 8L
  rank <- 3L
  draw_count <- 50000L
  diagonal <- seq(0.3, 0.9, length.out = variable_count)
  low_rank <- matrix(
    stats::rnorm(variable_count * rank, sd = 0.22),
    nrow = variable_count,
    ncol = rank
  )
  independent_draw <- matrix(
    stats::rnorm(draw_count * variable_count),
    nrow = draw_count,
    ncol = variable_count
  )
  independent_draw <- sweep(independent_draw, 2L, sqrt(diagonal), "*")
  latent_draw <- matrix(
    stats::rnorm(draw_count * rank),
    nrow = draw_count,
    ncol = rank
  )
  draws <- independent_draw + latent_draw %*% t(low_rank)
  marginal_sd <- sqrt(diagonal + rowSums(low_rank^2))
  standardized <- sweep(draws, 2L, marginal_sd, "/")
  coverage <- mean(abs(standardized) <= stats::qnorm(0.975))
  mean_square <- mean(standardized^2)
  stopifnot(abs(coverage - 0.95) < 0.004,
            abs(mean_square - 1) < 0.015)
})

run_test("variance_guards_are_explicit", {
  negative_error <- tryCatch({
    s6eb_lowrank_gaussian_nlpd(
      observed = c(0, 0), mean = c(0, 0),
      diagonal_variance = c(1, -0.1),
      low_rank = matrix(numeric(), 2L, 0L)
    )
    FALSE
  }, error = function(condition) TRUE)
  floored <- s6eb_lowrank_gaussian_nlpd(
    observed = c(0, 0), mean = c(0, 0),
    diagonal_variance = c(1, 1e-14),
    low_rank = matrix(numeric(), 2L, 0L),
    variance_floor = 1e-10
  )
  stopifnot(negative_error, floored$diagonal_floor_count == 1L)
})

run_test("source_has_no_truth_or_formal_fit_path", {
  source_text <- paste(
    readLines(
      file.path(script_root, "stage6eb_logscore_tools_20260910_v1.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  forbidden_path_tokens <- c(
    "sealed_truth", "truth_bundle[.]rds", "[/\\\\]evaluation[/\\\\]",
    "formal_fits"
  )
  stopifnot(!any(vapply(
    forbidden_path_tokens,
    grepl,
    logical(1L),
    x = source_text,
    ignore.case = TRUE
  )))
})

results <- do.call(rbind, test_results)
write_csv(results, file.path(output_root, "TEST_RESULTS.csv"))
write_lines(c(
  paste0("status=", if (all(results$passed)) "PASS" else "FAIL"),
  paste0("tests_passed=", sum(results$passed)),
  paste0("tests_total=", nrow(results)),
  "truth_accessed=FALSE",
  "fit_started=FALSE",
  "preregistered_rmse_interface_modified=FALSE"
), file.path(output_root, "TEST_COMPLETE.txt"))

if (!all(results$passed)) {
  stop("One or more Stage 6E-B targeted tests failed.", call. = FALSE)
}
cat("STAGE6E_B_LOGSCORE_TESTS_PASS tests=", nrow(results), "\n", sep = "")
