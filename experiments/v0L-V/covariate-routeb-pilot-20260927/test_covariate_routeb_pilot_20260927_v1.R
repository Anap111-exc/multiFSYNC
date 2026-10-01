#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
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
  key <- sub("=.*$", "", value[keep])
  answer <- sub("^[^=]*=", "", value[keep])
  stats::setNames(as.list(answer), key)
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
assert(length(script_arg) == 1L, "Cannot locate this test script.")
script_path <- sub("^--file=", "", script_arg[[1L]])
script_root <- dirname(script_path)
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
smoke_root <- normalizePath(cli$`smoke-root` %||% "", mustWork = TRUE)

checks <- logical()
record <- function(name, ok) {
  checks[[name]] <<- isTRUE(ok)
  assert(ok, paste0("FAILED: ", name))
}

candidate <- utils::read.csv(
  file.path(script_root, "FIT_MANIFEST_CANDIDATE_V1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
record("candidate_has_36_unique_fits",
       nrow(candidate) == 36L && !anyDuplicated(candidate$fit_id))
record("candidate_has_12_per_arm",
       identical(as.integer(table(candidate$arm)[c("C0", "C1", "C2")]),
                 c(12L, 12L, 12L)))
record("candidate_has_paired_fit_seeds",
       length(unique(candidate$fit_seed)) == 12L &&
         all(table(candidate$fit_seed) == 3L))
record("candidate_scientific_route_fixed",
       all(candidate$fit_K == 9L) &&
         all(candidate$function_initialization == "gram_unit_energy") &&
         all(candidate$pre_score_sweeps == 1L) &&
         all(candidate$planned_t1_sweeps == 400L) &&
         all(candidate$n_cpus == 1L))
record("candidate_truth_isolation",
       all(candidate$truth_available_to_fit == FALSE))

authorization_path <- file.path(
  script_root, "PILOT_IMPLEMENTATION_AND_EXECUTION_AUTHORIZATION_20260927.txt")
authorization <- read_key_value(authorization_path)
record("execution_authorized",
       identical(authorization$ecs_36_fit_pilot_authorized_after_local_acceptance,
                 "TRUE"))
record("truth_unseal_not_authorized",
       identical(authorization$truth_unseal_authorized, "FALSE") &&
         identical(authorization$scientific_evaluation_authorized, "FALSE"))
assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
protocol_path <- file.path(script_root, authorization$bound_protocol_file)
record("authorization_protocol_hash_matches",
       identical(digest::digest(protocol_path, algo = "sha256", file = TRUE,
                                serialize = FALSE),
                 authorization$bound_protocol_sha256))

data_manifest <- utils::read.csv(file.path(smoke_root, "DATA_MANIFEST.csv"),
                                 stringsAsFactors = FALSE)
fit_manifest <- utils::read.csv(file.path(smoke_root, "FIT_MANIFEST.csv"),
                                stringsAsFactors = FALSE)
record("smoke_has_three_paired_arms",
       nrow(data_manifest) == 3L &&
         identical(data_manifest$beta_arm, c("C0", "C1", "C2")))
record("smoke_has_three_unique_fits",
       nrow(fit_manifest) == 3L && !anyDuplicated(fit_manifest$fit_id))
record("smoke_fit_seeds_are_paired",
       length(unique(fit_manifest$fit_seed)) == 1L)

qc <- utils::read.csv(file.path(smoke_root, "PAIRING_QC.csv"),
                      stringsAsFactors = FALSE)
record("all_generation_and_pairing_checks_pass", nrow(qc) >= 18L &&
         all(qc$passed))

read_observation <- function(arm) readRDS(file.path(
  smoke_root, "data", paste0("covrbp_01__", arm),
  "observation_bundle.rds"))
c0 <- read_observation("C0")
c1 <- read_observation("C1")
c2 <- read_observation("C2")
record("observation_bundle_class_is_exact",
       all(vapply(list(c0, c1, c2), function(x) identical(
         x$bundle_class, "covariate_routeb_pilot_observation_only"),
         logical(1L))))
record("observation_bundles_do_not_expose_truth",
       all(vapply(list(c0, c1, c2), function(x) {
         isFALSE(x$truth_available_to_fit) &&
           !any(c("true_params", "data", "sealed") %in% names(x))
       }, logical(1L))))
record("C0_and_C1_include_two_covariates",
       is.list(c0$Z) && is.list(c1$Z) &&
         c0$dimensions$d == 2L && c1$dimensions$d == 2L)
record("C2_omits_covariates",
       is.null(c2$Z) && c2$dimensions$d == 0L && c2$generated_d == 2L)
record("C1_C2_reuse_identical_observations",
       identical(c1$Y, c2$Y) && identical(c1$time_obs, c2$time_obs))
record("C0_C1_pairing_is_preserved",
       identical(c0$Z, c1$Z) && identical(c0$time_obs, c1$time_obs) &&
         !identical(c0$Y, c1$Y))
record("sealed_truth_exists_but_is_not_opened_by_test",
       all(file.exists(file.path(
         smoke_root, "data", paste0("covrbp_01__", c("C0", "C1", "C2")),
         "sealed_truth", "truth_bundle.rds"))))

terminal_paths <- file.path(smoke_root, "fits", fit_manifest$fit_id,
                            "terminal_record.rds")
record("all_terminal_records_exist", all(file.exists(terminal_paths)))
terminals <- lapply(terminal_paths, readRDS)
record("all_smoke_fits_complete_and_eligible",
       all(vapply(terminals, function(x) {
         identical(x$status, "COMPLETE") && isTRUE(x$objective_eligible)
       }, logical(1L))))
record("warnings_are_retained",
       all(vapply(terminals, function(x) length(x$warnings) >= 1L,
                  logical(1L))))
record("fit_truth_isolation_and_single_core",
       all(vapply(terminals, function(x) {
         isFALSE(x$truth_used_for_fit_or_selection) && x$n_cpus == 1L
       }, logical(1L))))

winners <- utils::read.csv(file.path(
  smoke_root, "truth_free_selection", "TRUTH_FREE_WINNERS.csv"),
  stringsAsFactors = FALSE)
record("truth_free_selection_has_one_winner_per_arm",
       nrow(winners) == 3L &&
         setequal(winners$beta_arm, c("C0", "C1", "C2")) &&
         all(winners$truth_used_for_fit_or_selection == FALSE))
record("no_scientific_evaluation_started",
       !dir.exists(file.path(smoke_root, "evaluation")) &&
         file.exists(file.path(smoke_root, "MICRO_SMOKE_COMPLETE.txt")))

cat("TARGETED_TEST_PASS checks=", length(checks), "\n", sep = "")
