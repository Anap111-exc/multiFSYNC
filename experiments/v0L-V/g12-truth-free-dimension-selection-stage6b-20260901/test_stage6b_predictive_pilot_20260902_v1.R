#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(file_argument)) {
  sub("^--file=", "", file_argument[[1L]])
} else "test_stage6b_predictive_pilot_20260902_v1.R"
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, mustWork = TRUE))
}
source(file.path(script_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = TRUE)
source(file.path(script_root,
                 "stage6b_predictive_pilot_tools_20260902_v1.R"),
       local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
output_argument <- arguments[startsWith(arguments, "--output-root=")]
output_root <- if (length(output_argument)) {
  sub("--output-root=", "", output_argument[[1L]], fixed = TRUE)
} else tempfile("stage6b-pilot-tests-")

results <- list()
run_test <- function(test_id, expression) {
  started <- proc.time()[["elapsed"]]
  error <- tryCatch({ force(expression); NULL }, error = identity)
  results[[length(results) + 1L]] <<- data.frame(
    test_id = test_id, passed = is.null(error),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    detail = if (is.null(error)) "PASS" else conditionMessage(error),
    stringsAsFactors = FALSE
  )
  invisible(is.null(error))
}

configs <- s6bp_read_csv(file.path(
  script_root, "FIT_CONFIGS_STAGE6B_PILOT_V1.csv"
))
seeds <- s6bp_read_csv(file.path(
  script_root, "START_SEEDS_STAGE6B_PILOT_V1.csv"
))

run_test("registered_grid_and_seed_contract", {
  s6bp_validate_registration(configs, seeds)
  stopifnot(nrow(configs) == 9L, nrow(seeds) == 4L,
            identical(as.integer(seeds$fit_seed), 83206101:83206104))
})

manifest <- s6bp_build_fit_manifest(configs, seeds)
run_test("conditional_manifest_contract", {
  stopifnot(nrow(manifest) == 36L, !anyDuplicated(manifest$fit_id),
            nrow(s6bp_phase_manifest(manifest, "factor_screen")) == 12L,
            nrow(s6bp_phase_manifest(manifest, "fpca_screen", 1L)) == 8L,
            nrow(s6bp_phase_manifest(manifest, "fpca_screen", 2L)) == 8L,
            nrow(s6bp_phase_manifest(manifest, "fpca_screen", 3L)) == 8L,
            all(manifest$n_cpus == 1L),
            !any(manifest$truth_available_to_fit),
            !any(manifest$truth_available_to_dimension_selection),
            !any(manifest$continuation))
})

run_test("zero_warning_table_is_safe_and_fit_is_saved_first", {
  empty_warnings <- data.frame(
    warning_index = integer(), condition_class = character(),
    message = character(), phase = character(), stringsAsFactors = FALSE
  )
  labelled_empty <- s6bp_set_warning_phase(
    empty_warnings, "stage6b_predictive_pilot_g12_fixed400"
  )
  stopifnot(nrow(labelled_empty) == 0L,
            identical(labelled_empty$phase, character()))

  one_warning <- data.frame(
    warning_index = 1L, condition_class = "simpleWarning",
    message = "test", phase = "old", stringsAsFactors = FALSE
  )
  labelled_one <- s6bp_set_warning_phase(
    one_warning, "stage6b_predictive_pilot_g12_fixed400"
  )
  stopifnot(identical(labelled_one$phase,
                      "stage6b_predictive_pilot_g12_fixed400"))

  runner <- readLines(file.path(
    script_root, "run_stage6b_predictive_pilot_20260902_v1.R"
  ), warn = FALSE)
  save_line <- grep("pp_save_rds\\(result\\$fit", runner)
  label_line <- grep("result\\$captured\\$warnings <- s6bp_set_warning_phase",
                     runner)
  stopifnot(length(save_line) == 1L, length(label_line) == 1L,
            save_line < label_line)
})

make_losses <- function(config_id, mse_values) {
  data.frame(
    fit_config_id = config_id,
    study = rep(1:2, each = 4L),
    subject = rep(seq_len(4L), 2L),
    mse = mse_values,
    mae = sqrt(mse_values),
    observed_second_moment = rep(4, 8L),
    stringsAsFactors = FALSE
  )
}

run_test("one_se_prefers_simplest_tied_configuration", {
  metadata <- data.frame(
    fit_config_id = c("L1", "L2", "L3"), complexity_value = 1:3,
    stringsAsFactors = FALSE
  )
  losses <- rbind(
    make_losses("L1", c(1.00, 1.10, 0.90, 1.00, 1.10, 0.90, 1.00, 1.00)),
    make_losses("L2", c(0.995, 1.095, 0.895, 0.995,
                        1.095, 0.895, 0.995, 0.995)),
    make_losses("L3", c(0.99, 1.09, 0.89, 0.99,
                        1.09, 0.89, 0.99, 0.99))
  )
  selected <- s6bp_select_one_se(losses, metadata, "factor_screen")
  stopifnot(identical(selected$minimum_loss_config_id, "L3"),
            identical(selected$selected_config_id, "L1"),
            sum(selected$summary$selected_by_one_se) == 1L)
})

run_test("one_se_retains_material_predictive_gain", {
  metadata <- data.frame(
    fit_config_id = c("M2", "M3", "M4"), complexity_value = 2:4,
    stringsAsFactors = FALSE
  )
  losses <- rbind(
    make_losses("M2", rep(1.00, 8L)),
    make_losses("M3", rep(0.80, 8L)),
    make_losses("M4", rep(0.25, 8L))
  )
  selected <- s6bp_select_one_se(losses, metadata, "fpca_screen")
  stopifnot(identical(selected$minimum_loss_config_id, "M4"),
            identical(selected$selected_config_id, "M4"))
})

run_test("paired_subject_time_key_guard", {
  metadata <- data.frame(
    fit_config_id = c("A", "B", "C"), complexity_value = 1:3,
    stringsAsFactors = FALSE
  )
  losses <- rbind(
    make_losses("A", rep(1, 8L)),
    make_losses("B", rep(1, 8L)),
    make_losses("C", rep(1, 8L))
  )
  losses$subject[losses$fit_config_id == "C" & losses$study == 2L][1L] <- 99L
  error <- tryCatch({
    s6bp_select_one_se(losses, metadata, "factor_screen")
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl("keys differ", conditionMessage(error), fixed = TRUE))
})

run_test("runner_static_truth_source_guard", {
  runner <- readLines(file.path(
    script_root, "run_stage6b_predictive_pilot_20260902_v1.R"
  ), warn = FALSE)
  source_text <- paste(runner, collapse = "\n")
  stopifnot(!grepl("sealed[_-]truth", source_text, ignore.case = TRUE),
            !grepl("truth[_-]bundle[.]rds", source_text, ignore.case = TRUE),
            !grepl("[/\\\\]evaluation[/\\\\]", source_text,
                   ignore.case = TRUE))
})

run_test("registered_primary_metric_contract", {
  protocol <- paste(readLines(file.path(
    script_root, "PROTOCOL_G12_STAGE6B_PREDICTIVE_PILOT_V1_20260902.md"
  ), warn = FALSE), collapse = "\n")
  stopifnot(grepl("primary score is held-out signal NRMSE", protocol,
                 fixed = TRUE),
            grepl("at most 20 fits", protocol, fixed = TRUE),
            grepl("does not unseal truth", protocol, fixed = TRUE))
})

table <- do.call(rbind, results)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
s6b_write_csv(table, file.path(output_root, "TEST_RESULTS.csv"))
s6b_write_lines(c(
  paste0("status=", if (all(table$passed)) "PASS" else "FAIL"),
  paste0("tests_passed=", sum(table$passed)),
  paste0("tests_total=", nrow(table)),
  "fit_started=FALSE", "truth_accessed=FALSE"
), file.path(output_root, "TEST_COMPLETE.txt"))
if (!all(table$passed)) {
  stop(paste0("Stage 6B pilot tests failed: ",
              paste(table$test_id[!table$passed], collapse = ";")),
       call. = FALSE)
}
cat("STAGE6B_PILOT_TESTS_PASS tests=", nrow(table),
    " output_root=", output_root, "\n", sep = "")
