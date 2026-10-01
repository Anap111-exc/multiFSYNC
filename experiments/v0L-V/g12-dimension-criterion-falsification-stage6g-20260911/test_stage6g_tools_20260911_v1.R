#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[[1L]])
} else "test_stage6g_tools_20260911_v1.R"
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(
    script_file, winslash = "/", mustWork = TRUE
  ))
}
experiments_root <- dirname(script_root)
source(file.path(
  experiments_root, "g12-joint-dimension-selection-stage6f-20260911",
  "stage6f_joint_selection_tools_20260911_v1.R"
), local = TRUE)
source(file.path(script_root, "stage6g_tools_20260911_v1.R"), local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
if (length(output_argument) != 1L) {
  stop("--output-root is required exactly once.", call. = FALSE)
}
output_root <- sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Test output root already exists; refusing overwrite.", call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

results <- list()
run_test <- function(test_id, expression) {
  started <- proc.time()[["elapsed"]]
  failure <- NULL
  tryCatch(force(expression), error = function(condition) {
    failure <<- conditionMessage(condition)
  })
  passed <- is.null(failure)
  results[[length(results) + 1L]] <<- data.frame(
    test_id = test_id, passed = passed,
    detail = if (passed) "PASS" else failure,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  )
  invisible(passed)
}

scenarios <- s6g_read_csv(file.path(script_root, "SCENARIOS_STAGE6G_V1.csv"))
configs <- s6g_read_csv(file.path(script_root, "FIT_CONFIGS_STAGE6G_V1.csv"))
splits <- s6g_read_csv(file.path(script_root, "HOLDOUT_SPLITS_STAGE6G_V1.csv"))
starts <- s6g_read_csv(file.path(script_root, "START_SEEDS_STAGE6G_V1.csv"))

run_test("registration_and_120_endpoint_accounting", {
  s6g_validate_registration(scenarios, configs, splits, starts)
  manifest <- s6g_build_fit_manifest(scenarios, configs, splits, starts)
  truth_names <- c(
    "truth_basis_seed", "expected_truth_fit_config_id", "truth_L_f",
    "truth_L_s_1", "truth_L_s_2", "truth_M_common"
  )
  stopifnot(
    nrow(manifest) == 120L,
    length(unique(manifest$selection_stratum_id)) == 30L,
    all(table(manifest$selection_stratum_id) == 4L),
    all(table(manifest$data_id) == 40L),
    all(manifest$n_cpus == 1L),
    !length(intersect(names(manifest), truth_names))
  )
})

run_test("package_source_binding_is_exact", {
  repo_root <- Sys.getenv(
    "STAGE6G_PACKAGE_SOURCE_ROOT",
    unset = file.path(script_root, "..", "..", "..")
  )
  binding <- s6g_read_csv(file.path(
    script_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
  ))
  paths <- file.path(repo_root, binding$relative_path)
  hashes <- vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L))
  stopifnot(
    nrow(binding) == 29L,
    all(file.exists(paths)),
    identical(as.numeric(file.info(paths)$size), as.numeric(binding$bytes)),
    requireNamespace("digest", quietly = TRUE),
    identical(unname(hashes), as.character(binding$sha256))
  )
})

make_subject_losses <- function(base, spread = 0.01) {
  value <- expand.grid(
    fit_config_id = as.character(configs$fit_config_id),
    study = 1:2, subject = 1:30,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  subject_effect <- seq(-spread, spread, length.out = 30L)
  value$mse <- unname(base[value$fit_config_id]) +
    subject_effect[value$subject]
  value
}

run_test("one_se_prefers_simplest_when_supported", {
  base <- setNames(c(1.000, 0.999, 1.002, 1.001, 0.998),
                   configs$fit_config_id)
  selection <- s6g_select_one_se(
    make_subject_losses(base, 0.04), configs, "mse"
  )
  stopifnot(
    selection$minimum_loss_config_id == "cand_L211_M2",
    selection$selected_config_id == "cand_L111_M2",
    nrow(selection$summary) == 5L
  )
})

run_test("one_se_can_select_nonminimal_factor_truth", {
  base <- setNames(c(1.00, 1.01, 1.02, 0.95, 0.50),
                   configs$fit_config_id)
  selection <- s6g_select_one_se(
    make_subject_losses(base, 0.01), configs, "mse"
  )
  stopifnot(
    selection$minimum_loss_config_id == "cand_L211_M2",
    selection$selected_config_id == "cand_L211_M2"
  )
})

run_test("one_se_can_select_nonminimal_rank_truth", {
  base <- setNames(c(1.00, 0.50, 1.02, 1.01, 1.03),
                   configs$fit_config_id)
  selection <- s6g_select_one_se(
    make_subject_losses(base, 0.01), configs, "mse"
  )
  stopifnot(
    selection$minimum_loss_config_id == "cand_L111_M3",
    selection$selected_config_id == "cand_L111_M3"
  )
})

run_test("split_ranking_and_pairwise_order_are_measured", {
  score <- expand.grid(
    fit_config_id = as.character(configs$fit_config_id),
    split_id = as.character(splits$split_id), study = 1:2, subject = 1:30,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  base <- setNames(seq_len(nrow(configs)), configs$fit_config_id)
  score$mse <- unname(base[score$fit_config_id])
  stability <- s6g_split_ranking_stability(score, "mse")
  stopifnot(
    abs(stability$summary$spearman_rank_correlation - 1) < 1e-12,
    stability$summary$minimum_loss_winner_agreement,
    stability$summary$pairwise_ordering_concordance == 1,
    nrow(stability$pairwise) == 10L
  )
})

run_test("paired_contrasts_and_leave_one_subject_out_are_complete", {
  base <- setNames(c(1.00, 1.01, 1.02, 1.03, 1.04),
                   configs$fit_config_id)
  losses <- make_subject_losses(base, 0.02)
  selection <- s6g_select_one_se(losses, configs, "mse")
  contrasts <- s6g_paired_contrasts(losses, selection)
  leave_one <- s6g_leave_one_subject_out(losses, configs)
  stopifnot(
    nrow(contrasts) == 5L,
    all(is.finite(contrasts$maximum_leave_one_subject_influence)),
    nrow(leave_one$decisions) == 60L,
    nrow(leave_one$counts) == 5L,
    sum(leave_one$counts$leave_one_subject_out_selections) == 60L
  )
})

run_test("authorization_preserves_truth_and_continuation_lock", {
  authorization <- readLines(file.path(
    script_root, "STAGE6G_EXECUTION_AUTHORIZATION_20260911.txt"
  ), warn = FALSE)
  stopifnot(
    "status=AUTHORIZED_TRUTH_FREE_STAGE6G_SCREEN" %in% authorization,
    "truth_unseal_authorized=FALSE" %in% authorization,
    "continuation_authorized=FALSE" %in% authorization,
    "formal_paper_mc=FALSE" %in% authorization
  )
})

run_test("runner_has_no_unseal_or_evaluation_action", {
  runner <- paste(readLines(file.path(
    script_root, "run_stage6g_20260911_v1.R"
  ), warn = FALSE), collapse = "\n")
  stopifnot(
    !grepl("candidate2_evaluate_multi", runner, fixed = TRUE),
    !grepl("action, \"evaluate\"", runner, fixed = TRUE),
    !grepl("action, \"unseal\"", runner, fixed = TRUE),
    grepl("cross_dimension_elbo_comparison=FALSE", runner, fixed = TRUE)
  )
})

result <- do.call(rbind, results)
utils::write.csv(result, file.path(output_root, "TEST_RESULTS.csv"),
                 row.names = FALSE, na = "")
writeLines(c(
  paste0("status=", if (all(result$passed)) "PASS" else "FAIL"),
  paste0("tests_passed=", sum(result$passed)),
  paste0("tests_total=", nrow(result)),
  "fit_started=FALSE", "truth_unsealed=FALSE", "continuation_started=FALSE"
), file.path(output_root, "TEST_COMPLETE.txt"), useBytes = TRUE)
if (!all(result$passed)) {
  stop("Stage 6G targeted tests failed.", call. = FALSE)
}
cat("STAGE6G_TARGETED_TESTS_PASS tests=", nrow(result), "\n", sep = "")
