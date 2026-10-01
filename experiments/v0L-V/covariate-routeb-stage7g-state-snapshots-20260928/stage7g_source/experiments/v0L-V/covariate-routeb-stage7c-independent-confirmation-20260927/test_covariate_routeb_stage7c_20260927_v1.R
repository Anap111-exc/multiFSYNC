#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
parse_cli <- function(x) {
  out <- list()
  for (item in x) {
    assert(grepl("^--[^=]+=", item), "Use --key=value arguments.")
    key <- sub("^--([^=]+)=.*$", "\\1", item)
    out[[key]] <- sub("^[^=]+=", "", item)
  }
  out
}
read_key_value <- function(path) {
  value <- readLines(path, warn = FALSE)
  keep <- grepl("=", value, fixed = TRUE)
  stats::setNames(as.list(sub("^[^=]*=", "", value[keep])),
                  sub("=.*$", "", value[keep]))
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
assert(length(script_arg) == 1L, "Cannot locate this test script.")
script_root <- dirname(sub("^--file=", "", script_arg[[1L]]))
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
smoke_root <- cli$`smoke-root` %||% ""
assert(dir.exists(smoke_root), "Smoke root does not exist.")

checks <- logical()
record <- function(name, ok) {
  checks[[name]] <<- isTRUE(ok)
  assert(ok, paste0("FAILED: ", name))
}

data_candidate <- read.csv(file.path(
  script_root, "DATA_MANIFEST_CANDIDATE_V1.csv"), check.names = FALSE)
fit_candidate <- read.csv(file.path(
  script_root, "FIT_MANIFEST_CANDIDATE_V1.csv"), check.names = FALSE)
record("candidate_has_3_logical_data_and_9_arms",
       nrow(data_candidate) == 9L &&
         length(unique(data_candidate$logical_data_id)) == 3L &&
         all(table(data_candidate$logical_data_id) == 3L))
record("candidate_has_108_unique_fits",
       nrow(fit_candidate) == 108L && !anyDuplicated(fit_candidate$fit_id))
record("candidate_has_12_starts_per_seed_arm",
       all(table(fit_candidate$logical_data_id,
                 fit_candidate$beta_arm) == 12L))
record("candidate_fit_seeds_are_paired_only_within_logical_data",
       length(unique(fit_candidate$fit_seed)) == 36L &&
         all(table(fit_candidate$fit_seed) == 3L))
record("candidate_scientific_route_fixed",
       all(fit_candidate$fit_K == 9L) && all(fit_candidate$n_g == 51L) &&
         all(fit_candidate$lambda_orth == 0) &&
         all(fit_candidate$function_initialization == "gram_unit_energy") &&
         all(fit_candidate$pre_score_sweeps == 1L) &&
         all(fit_candidate$planned_t1_sweeps == 400L) &&
         all(fit_candidate$n_cpus == 1L))
record("candidate_truth_and_execution_flags_false",
       !any(fit_candidate$truth_available_to_fit) &&
         !any(fit_candidate$truth_available_to_stopping) &&
         !any(fit_candidate$truth_available_to_selection) &&
         !any(fit_candidate$execution_authorized) &&
         !any(fit_candidate$fit_started))

assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
auth <- read_key_value(file.path(
  script_root, "STAGE7C_PROTOCOL_SEED_AND_EXECUTION_AUTHORIZATION_20260927.txt"))
digest_file <- function(name) digest::digest(
  file.path(script_root, name), algo = "sha256", file = TRUE,
  serialize = FALSE)
record("authorization_binds_protocol_and_manifests",
       identical(auth$protocol_confirmed, "TRUE") &&
         identical(auth$seed_schedule_confirmed, "TRUE") &&
         identical(digest_file(auth$bound_protocol_file),
                   auth$bound_protocol_sha256) &&
         identical(digest_file(auth$bound_data_manifest_file),
                   auth$bound_data_manifest_sha256) &&
         identical(digest_file(auth$bound_fit_manifest_file),
                   auth$bound_fit_manifest_sha256))
record("authorization_preserves_truth_isolation",
       identical(auth$truth_unseal_authorized, "FALSE") &&
         identical(auth$scientific_evaluation_authorized, "FALSE") &&
         identical(auth$continuation_authorized, "FALSE"))

data_manifest <- read.csv(file.path(smoke_root, "DATA_MANIFEST.csv"))
fit_manifest <- read.csv(file.path(smoke_root, "FIT_MANIFEST.csv"))
record("smoke_uses_three_paired_arms_and_dedicated_seeds",
       nrow(data_manifest) == 3L &&
         identical(data_manifest$beta_arm, c("C0", "C1", "C2")) &&
         all(data_manifest$logical_data_id == "covrb7c_smoke") &&
         all(data_manifest$data_seed == 99719001L))
record("smoke_has_three_unique_fits_with_one_paired_seed",
       nrow(fit_manifest) == 3L && !anyDuplicated(fit_manifest$fit_id) &&
         length(unique(fit_manifest$fit_seed)) == 1L &&
         unique(fit_manifest$fit_seed) == 99719601L)

qc <- read.csv(file.path(smoke_root, "PAIRING_QC.csv"))
record("all_generation_and_pairing_checks_pass",
       nrow(qc) >= 18L && all(qc$passed))
read_observation <- function(arm) readRDS(file.path(
  smoke_root, "data", paste0("covrb7c_smoke__", arm),
  "observation_bundle.rds"))
c0 <- read_observation("C0")
c1 <- read_observation("C1")
c2 <- read_observation("C2")
record("observation_bundles_do_not_expose_truth",
       all(vapply(list(c0, c1, c2), function(x) {
         isFALSE(x$truth_available_to_fit) &&
           !any(c("true_params", "data", "sealed") %in% names(x))
       }, logical(1L))))
record("C0_C1_include_two_covariates_and_C2_omits_them",
       is.list(c0$Z) && is.list(c1$Z) && is.null(c2$Z) &&
         c0$dimensions$d == 2L && c1$dimensions$d == 2L &&
         c2$dimensions$d == 0L && c2$generated_d == 2L)
record("C1_C2_reuse_identical_observations",
       identical(c1$Y, c2$Y) && identical(c1$time_obs, c2$time_obs))
record("C0_C1_pairing_is_preserved",
       identical(c0$Z, c1$Z) && identical(c0$time_obs, c1$time_obs) &&
         !identical(c0$Y, c1$Y))
record("sealed_truth_exists_but_test_does_not_open_it",
       all(file.exists(file.path(
         smoke_root, "data",
         paste0("covrb7c_smoke__", c("C0", "C1", "C2")),
         "sealed_truth", "truth_bundle.rds"))))

terminal_paths <- file.path(smoke_root, "fits", fit_manifest$fit_id,
                            "terminal_record.rds")
record("all_terminal_records_exist", all(file.exists(terminal_paths)))
terminals <- lapply(terminal_paths, readRDS)
record("all_smoke_fits_complete_and_eligible",
       all(vapply(terminals, function(x) {
         identical(x$status, "COMPLETE") && isTRUE(x$objective_eligible)
       }, logical(1L))))
record("warning_and_error_fields_are_retained",
       all(vapply(seq_along(terminals), function(i) {
         is.character(terminals[[i]]$warnings) &&
           "error" %in% names(terminals[[i]]) &&
           file.exists(file.path(dirname(terminal_paths[[i]]), "WARNINGS.csv"))
       }, logical(1L))))
record("fit_truth_isolation_and_single_core",
       all(vapply(terminals, function(x) {
         isFALSE(x$truth_used_for_fit_or_selection) && x$n_cpus == 1L
       }, logical(1L))))

winners <- read.csv(file.path(
  smoke_root, "truth_free_selection", "TRUTH_FREE_WINNERS.csv"))
record("truth_free_selection_has_one_winner_per_arm",
       nrow(winners) == 3L && setequal(winners$beta_arm, c("C0", "C1", "C2")) &&
         !any(winners$truth_used_for_fit_or_selection))
record("no_scientific_evaluation_started",
       !dir.exists(file.path(smoke_root, "evaluation")) &&
         file.exists(file.path(smoke_root, "MICRO_SMOKE_COMPLETE.txt")))

cat("STAGE7C_TARGETED_TEST_PASS checks=", length(checks), "\n", sep = "")
