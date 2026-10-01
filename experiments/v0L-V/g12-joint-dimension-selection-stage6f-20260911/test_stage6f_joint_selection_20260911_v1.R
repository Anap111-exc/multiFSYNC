#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[[1L]])
} else "test_stage6f_joint_selection_20260911_v1.R"
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(script_root,
                 "stage6f_joint_selection_tools_20260911_v1.R"),
       local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
if (!length(output_argument)) stop("--output-root is required.", call. = FALSE)
output_root <- sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Test output root already exists; refusing overwrite.", call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

results <- list()
run_test <- function(test_id, expression) {
  started <- proc.time()[["elapsed"]]
  failure <- NULL
  tryCatch(force(expression),
           error = function(condition) failure <<- conditionMessage(condition))
  passed <- is.null(failure)
  results[[length(results) + 1L]] <<- data.frame(
    test_id = test_id, passed = passed,
    detail = if (passed) "PASS" else failure,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  )
  invisible(passed)
}

configs <- s6f_read_csv(file.path(script_root, "FIT_CONFIGS_STAGE6F_V1.csv"))
seeds <- s6f_read_csv(file.path(script_root, "START_SEEDS_STAGE6F_V1.csv"))
splits <- s6f_read_csv(file.path(script_root, "HOLDOUT_SPLITS_STAGE6F_V1.csv"))

run_test("registration_and_20_fit_accounting", {
  s6f_validate_registration(configs, seeds, splits)
  manifest <- s6f_build_manifest(configs, seeds, splits)
  stopifnot(nrow(manifest) == 32L, sum(manifest$reuse_stage6b) == 12L,
            sum(manifest$execute_new_fit) == 20L,
            all(table(manifest$selection_stratum_id) == 4L),
            all(manifest$n_cpus == 1L))
})

run_test("paired_holdout_indices_are_distinct_and_interior", {
  make_subject <- function(n_time, shift) {
    time <- seq(0, 1, length.out = n_time)
    list(time = time, y = lapply(seq_len(3L), function(variable) {
      sin(time * (variable + 1L)) + shift
    }))
  }
  subjects <- lapply(6:9, make_subject, shift = 0.1)
  observation <- list(
    Y = list(lapply(subjects, `[[`, "y"), lapply(subjects, `[[`, "y")),
    time_obs = list(lapply(subjects, `[[`, "time"),
                    lapply(subjects, `[[`, "time")),
    Z = NULL, dimensions = list(d = 0L)
  )
  middle <- s6f_make_time_holdout(observation, "middle")
  upper <- s6f_make_time_holdout(observation, "upper")
  middle_index <- vapply(middle$holdout_plan$records, `[[`, integer(1L),
                           "original_index")
  upper_index <- vapply(upper$holdout_plan$records, `[[`, integer(1L),
                          "original_index")
  stopifnot(all(upper_index == middle_index + 1L),
            all(vapply(middle$holdout_plan$records, `[[`, integer(1L),
                       "train_time_count") >= 5L),
            all(vapply(upper$holdout_plan$records, `[[`, integer(1L),
                       "train_time_count") >= 5L))
})

run_test("subject_clustered_one_se_selects_simplest_supported", {
  set.seed(91101L)
  score <- expand.grid(
    fit_config_id = configs$fit_config_id,
    split_id = splits$split_id,
    study = 1:2, subject = 1:30,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  base <- c(joint_L1_M2 = 1.00, joint_L1_M4 = 1.01,
            joint_L2_M2 = 1.03, joint_L2_M4 = 1.05)
  score$mse <- unname(base[score$fit_config_id]) +
    rep(seq(-0.04, 0.04, length.out = 60L), times = 8L)
  clustered <- s6f_subject_cluster_losses(score, "mse")
  selection <- s6f_select_one_se(
    clustered, configs, "mse", "paired_holdout_mse"
  )
  stopifnot(selection$minimum_loss_config_id == "joint_L1_M2",
            selection$selected_config_id == "joint_L1_M2",
            nrow(selection$summary) == 4L)
})

run_test("cross_split_ranking_is_measured", {
  score <- expand.grid(
    fit_config_id = configs$fit_config_id,
    split_id = splits$split_id,
    study = 1:2, subject = 1:30,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  loss <- c(joint_L1_M2 = 1, joint_L1_M4 = 2,
            joint_L2_M2 = 3, joint_L2_M4 = 4)
  score$loss <- unname(loss[score$fit_config_id])
  stability <- s6f_split_ranking_stability(score, "loss", "unit_test")
  stopifnot(stability$summary$spearman_rank_correlation == 1,
            stability$summary$minimum_loss_winner_agreement,
            stability$summary$pairwise_ordering_concordance == 1,
            nrow(stability$pairwise) == 6L)
})

run_test("runner_keeps_truth_out_of_fit_and_selection_actions", {
  text <- paste(readLines(file.path(
    script_root, "run_stage6f_joint_dimension_selection_20260911_v1.R"
  ), warn = FALSE), collapse = "\n")
  stopifnot(!grepl("truth_bundle[.]rds", text, ignore.case = TRUE),
            !grepl("candidate2_evaluate_multi", text, fixed = TRUE))
})

result <- do.call(rbind, results)
utils::write.csv(result, file.path(output_root, "TEST_RESULTS.csv"),
                 row.names = FALSE, na = "")
writeLines(c(
  paste0("status=", if (all(result$passed)) "PASS" else "FAIL"),
  paste0("tests_passed=", sum(result$passed)),
  paste0("tests_total=", nrow(result)),
  "fit_started=FALSE", "truth_accessed=FALSE"
), file.path(output_root, "TEST_COMPLETE.txt"), useBytes = TRUE)
if (!all(result$passed)) stop("Stage 6F targeted tests failed.", call. = FALSE)
cat("STAGE6F_TARGETED_TESTS_PASS tests=", nrow(result), "\n", sep = "")
