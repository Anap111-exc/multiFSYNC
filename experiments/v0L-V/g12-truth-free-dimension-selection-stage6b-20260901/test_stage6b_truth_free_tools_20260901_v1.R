#!/usr/bin/env Rscript

script_argument <- commandArgs(trailingOnly = FALSE)
file_argument <- script_argument[startsWith(script_argument, "--file=")]
script_file <- if (length(file_argument)) {
  sub("--file=", "", file_argument[[1L]], fixed = TRUE)
} else {
  "test_stage6b_truth_free_tools_20260901_v1.R"
}
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
source(file.path(script_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
output_root <- if (length(output_argument)) {
  sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
} else {
  tempfile("stage6b-targeted-tests-")
}
observation_argument <- arguments[startsWith(arguments,
                                             "--observation-bundle=")]
observation_path <- if (length(observation_argument)) {
  sub("--observation-bundle=", "", observation_argument[[1L]], fixed = TRUE)
} else ""

test_results <- list()
test_index <- 0L
run_test <- function(test_id, expression) {
  test_index <<- test_index + 1L
  started <- proc.time()[["elapsed"]]
  error <- tryCatch({
    force(expression)
    NULL
  }, error = identity)
  test_results[[test_index]] <<- data.frame(
    test_id = test_id,
    passed = is.null(error),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    detail = if (is.null(error)) "PASS" else conditionMessage(error),
    stringsAsFactors = FALSE
  )
  invisible(is.null(error))
}

run_test("permutation_and_one_to_one_matching", {
  stopifnot(nrow(s6b_permutations(1L)) == 1L,
            nrow(s6b_permutations(2L)) == 2L,
            nrow(s6b_permutations(3L)) == 6L)
  similarity <- matrix(c(0.1, 0.9, 0.8, 0.2), nrow = 2L, byrow = TRUE)
  stopifnot(identical(as.integer(s6b_best_assignment(similarity)), c(2L, 1L)))
})

run_test("cosine_zero_and_sign_contract", {
  stopifnot(identical(s6b_abs_cosine(c(0, 0), c(1, 2)), 0),
            abs(s6b_abs_cosine(c(1, 2), c(-1, -2)) - 1) < 1e-12)
})

run_test("loading_and_fpca_subspace_contract", {
  left <- cbind(c(1, 0, 0), c(0, 1, 0))
  rotated <- left %*% matrix(c(0, 1, -1, 0), nrow = 2L)
  loading <- s6b_euclidean_subspace_similarity(left, rotated)
  stopifnot(abs(loading[["mean_squared_canonical_correlation"]] - 1) < 1e-12,
            abs(loading[["minimum_canonical_correlation"]] - 1) < 1e-12)

  time <- c(0, 0.25, 0.5, 0.75, 1)
  functions <- cbind(sin(pi * time), cos(pi * time))
  function_rotation <- functions %*% matrix(c(0, 1, -1, 0), nrow = 2L)
  fpca <- s6b_subspace_similarity(functions, function_rotation, time)
  stopifnot(abs(fpca[["mean_squared_canonical_correlation"]] - 1) < 1e-12,
            abs(fpca[["minimum_canonical_correlation"]] - 1) < 1e-12)
})

fake_subject <- function(offset) {
  list(offset + seq_len(6L), offset + 10 + seq_len(6L))
}
observation <- list(
  bundle_class = "observation_only_test",
  Y = list(list(fake_subject(0)), list(fake_subject(100))),
  time_obs = list(list(seq(0, 1, length.out = 6L)),
                  list(seq(0, 1, length.out = 6L))),
  Z = NULL,
  dimensions = list(S = 2L, n_s = c(1L, 1L), p = 2L, d = 0L)
)

run_test("whole_middle_time_holdout_contract", {
  heldout <- s6b_make_middle_time_holdout(observation, min_train_times = 5L)
  stopifnot(heldout$holdout_plan$subject_time_count == 2L,
            identical(heldout$holdout_plan$truth_used, FALSE),
            all(vapply(heldout$holdout_plan$records, `[[`, integer(1L),
                       "original_index") == 3L),
            all(vapply(heldout$holdout_plan$records, `[[`, integer(1L),
                       "train_time_count") == 5L),
            length(heldout$train_observation$time_obs[[1L]][[1L]]) == 5L,
            all(vapply(heldout$train_observation$Y[[1L]][[1L]], length,
                       integer(1L)) == 5L),
            identical(heldout$holdout_plan$records[[1L]]$y, c(3, 13)))
})

run_test("observation_truth_rejection", {
  contaminated <- observation
  contaminated$true_params <- list(secret = 1)
  error <- tryCatch({
    s6b_make_middle_time_holdout(contaminated, min_train_times = 5L)
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl("forbidden truth", conditionMessage(error),
                  ignore.case = TRUE))
})

if (nzchar(observation_path)) {
  run_test("stage6a_full_observation_split_contract", {
    full_observation <- readRDS(observation_path)
    stopifnot(!length(s6b_forbidden_observation_names(full_observation)))
    full_split <- s6b_make_middle_time_holdout(
      full_observation, min_train_times = 5L
    )
    records <- full_split$holdout_plan$records
    stopifnot(length(records) == sum(full_observation$dimensions$n_s),
              length(records) == 60L,
              all(vapply(records, `[[`, integer(1L), "p") == 500L),
              all(vapply(records, `[[`, integer(1L),
                         "train_time_count") %in% 5:8),
              identical(full_split$holdout_plan$truth_used, FALSE),
              isTRUE(full_split$holdout_plan$complete_variable_vector_held_out))
  })
}

fake_fit <- list(
  S = 2L,
  n_s = c(1L, 1L),
  p = 2L,
  d = 0L,
  time_g = c(0, 0.5, 1),
  list_mu_hat = list(
    list(c(1, 1, 1), c(2, 2, 2)),
    list(c(1, 1, 1), c(2, 2, 2))
  ),
  L_f = 1L,
  factor_ppi_shared = list(0.9),
  list_Phi_hat = list(matrix(c(0, 1, 0), ncol = 1L)),
  list_Zeta_hat = list(matrix(c(2, 3), nrow = 2L, ncol = 1L)),
  mu_q_a_hat = matrix(c(1, -1), nrow = 2L, ncol = 1L),
  L_s_by_study = c(1L, 1L),
  factor_ppi_specific = list(list(0.4), list(0.8)),
  list_Phi_hat_spec = list(
    list(matrix(c(1, 1, 1), ncol = 1L)),
    list(matrix(c(1, 1, 1), ncol = 1L))
  ),
  list_Zeta_hat_spec = list(
    list(matrix(0.5, nrow = 1L, ncol = 1L)),
    list(matrix(0.25, nrow = 1L, ncol = 1L))
  ),
  mu_q_b_specific_hat = list(
    matrix(c(2, 0), nrow = 2L, ncol = 1L),
    matrix(c(0, 4), nrow = 2L, ncol = 1L)
  ),
  sigsq_eps_hat = matrix(1, nrow = 2L, ncol = 2L)
)

run_test("posterior_mean_prediction_all_and_ppi", {
  all_prediction <- s6b_predict_subject_time(fake_fit, 1L, 1L, 0.5, "all")
  ppi_prediction <- s6b_predict_subject_time(fake_fit, 1L, 1L, 0.5,
                                             "ppi_0.5")
  study_two <- s6b_predict_subject_time(fake_fit, 2L, 1L, 0.5, "ppi_0.5")
  stopifnot(max(abs(all_prediction - c(4, 0))) < 1e-12,
            max(abs(ppi_prediction - c(3, 0))) < 1e-12,
            max(abs(study_two - c(4, 0))) < 1e-12)
})

run_test("holdout_score_semantics", {
  plan <- list(
    truth_used = FALSE,
    records = list(
      list(study = 1L, subject = 1L, time = 0.5, original_index = 3L,
           train_time_count = 5L, y = c(4, 0)),
      list(study = 2L, subject = 1L, time = 0.5, original_index = 3L,
           train_time_count = 5L, y = c(4, 0))
    )
  )
  score <- s6b_score_holdout(fake_fit, plan, "all")
  stopifnot(nrow(score$subject_time) == 2L,
            score$aggregate$heldout_rmse == 0,
            score$aggregate$heldout_nrmse == 0,
            identical(score$aggregate$truth_used, FALSE),
            identical(score$aggregate$primary_dimension_score_ready, FALSE),
            identical(
              score$aggregate$predictive_variance_scope,
              "residual_noise_only_not_calibrated_full_posterior"
            ))
})

run_test("audit_source_static_truth_path_guard", {
  audit_source <- c(
    file.path(script_root, "stage6b_truth_free_tools_20260901_v1.R"),
    file.path(script_root, "run_stage6ba_readonly_audit_20260901_v1.R")
  )
  source_text <- paste(unlist(lapply(audit_source, readLines, warn = FALSE)),
                       collapse = "\n")
  forbidden_path_tokens <- c("sealed_truth", "truth_bundle[.]rds",
                             "[/\\\\]evaluation[/\\\\]")
  stopifnot(!any(vapply(forbidden_path_tokens, grepl, logical(1L),
                        x = source_text, ignore.case = TRUE)))
})

results <- do.call(rbind, test_results)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
s6b_write_csv(results, file.path(output_root, "TEST_RESULTS.csv"))
s6b_write_lines(c(
  paste0("status=", if (all(results$passed)) "PASS" else "FAIL"),
  paste0("tests_passed=", sum(results$passed)),
  paste0("tests_total=", nrow(results)),
  "truth_accessed=FALSE",
  "fit_started=FALSE"
), file.path(output_root, "TEST_COMPLETE.txt"))

if (!all(results$passed)) {
  stop(paste0("Stage 6B targeted tests failed: ",
              paste(results$test_id[!results$passed], collapse = ";")),
       call. = FALSE)
}
cat("STAGE6B_TARGETED_TESTS_PASS tests=", nrow(results),
    " output_root=", output_root, "\n", sep = "")
