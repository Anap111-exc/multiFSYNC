#!/usr/bin/env Rscript

# Stage 7C intentionally reuses the already accepted Route-B covariate pilot
# implementation. Only the registered data/fit manifests and batch cardinality
# differ. The inherited file is an explicit execution-source dependency.
wrapper_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(wrapper_arg) != 1L) stop("Cannot locate Stage 7C runner.")
wrapper_path <- sub("^--file=", "", wrapper_arg[[1L]])
wrapper_root <- dirname(wrapper_path)
pilot_runner <- file.path(
  wrapper_root, "..", "covariate-routeb-pilot-20260927",
  "run_covariate_routeb_pilot_20260927_v1.R")
if (!file.exists(pilot_runner)) stop("Missing bound pilot runtime dependency.")
pilot_lines <- readLines(pilot_runner, warn = FALSE)
cli_boundary <- grep("^cli <- parse_cli", pilot_lines)
if (length(cli_boundary) != 1L) stop("Pilot runtime boundary is ambiguous.")
eval(parse(text = pilot_lines[seq_len(cli_boundary - 1L)]),
     envir = .GlobalEnv)

script_path <- wrapper_path
script_root <- wrapper_root
experiment_id <- "COVARIATE_ROUTEB_STAGE7C_V1_20260927"

read_key_value <- function(path) {
  value <- readLines(path, warn = FALSE)
  keep <- grepl("=", value, fixed = TRUE)
  stats::setNames(as.list(sub("^[^=]*=", "", value[keep])),
                  sub("=.*$", "", value[keep]))
}

data_manifest <- function(smoke = FALSE) {
  if (smoke) {
    return(data.frame(
      logical_data_id = rep("covrb7c_smoke", 3L),
      data_id = paste0("covrb7c_smoke__", c("C0", "C1", "C2")),
      data_index = rep(0L, 3L), beta_arm = c("C0", "C1", "C2"),
      observation_source_data_id = c(
        "covrb7c_smoke__C0", "covrb7c_smoke__C1", "covrb7c_smoke__C1"),
      data_seed = rep(99719001L, 3L),
      truth_basis_seed = rep(99719002L, 3L),
      truth_loading_seed = rep(99719003L, 3L),
      mean_truth_seed = rep(99719004L, 3L),
      beta_bspline_seed = rep(99719005L, 3L),
      d_generated = rep(2L, 3L), d_fit = c(2L, 2L, 0L),
      truth_available_to_fit = FALSE,
      truth_available_to_selection = FALSE,
      truth_unseal_authorized = FALSE,
      stringsAsFactors = FALSE))
  }
  registered <- read_manifest(file.path(
    script_root, "DATA_MANIFEST_CANDIDATE_V1.csv"))
  registered$beta_bspline_seed <- registered$beta_truth_seed
  registered$truth_unseal_authorized <- FALSE
  registered
}

fit_manifest <- function(data_rows, smoke = FALSE) {
  if (smoke) {
    return(data.frame(
      fit_id = sprintf("%s__start01", data_rows$data_id),
      data_id = data_rows$data_id,
      logical_data_id = data_rows$logical_data_id,
      beta_arm = data_rows$beta_arm,
      start_index = 1L, fit_seed = 99719601L, n_cpus = 1L,
      truth_available_to_fit = FALSE,
      truth_available_to_selection = FALSE,
      stringsAsFactors = FALSE))
  }
  read_manifest(file.path(script_root, "FIT_MANIFEST_CANDIDATE_V1.csv"))
}

generate <- function(output_root, smoke = FALSE) {
  load_package()
  assert(!dir.exists(output_root),
         paste0("Refusing to overwrite output root: ", output_root))
  dir.create(output_root, recursive = TRUE, mode = "0700")
  data_rows <- data_manifest(smoke)
  fit_rows <- fit_manifest(data_rows, smoke)
  expected_data <- if (smoke) 3L else 9L
  expected_fits <- if (smoke) 3L else 108L
  assert(!anyDuplicated(data_rows$data_id) &&
           !anyDuplicated(fit_rows$fit_id) &&
           nrow(data_rows) == expected_data && nrow(fit_rows) == expected_fits,
         "Stage 7C manifest dimensions or identifiers are invalid.")
  csv_new(data_rows, file.path(output_root, "DATA_MANIFEST.csv"))
  csv_new(fit_rows, file.path(output_root, "FIT_MANIFEST.csv"))
  checks <- list()
  for (logical_id in unique(data_rows$logical_data_id)) {
    row <- data_rows[data_rows$logical_data_id == logical_id, , drop = FALSE]
    assert(nrow(row) == 3L &&
             identical(as.character(row$beta_arm), c("C0", "C1", "C2")),
           "Paired Stage 7C data manifest is malformed.")
    pair <- paired_beta_data(row[1L, , drop = FALSE], smoke)
    qc <- pair_checks(pair, row$data_index[[1L]])
    qc$logical_data_id <- logical_id
    assert(all(qc$passed), paste0("Pair QC failed for ", logical_id, ": ",
      paste(qc$check_id[!qc$passed], collapse = ";")))
    checks[[length(checks) + 1L]] <- qc
    save_arm(pair$C0, row[row$beta_arm == "C0", , drop = FALSE],
             output_root, smoke)
    save_arm(pair$C1, row[row$beta_arm == "C1", , drop = FALSE],
             output_root, smoke)
    save_arm(pair$C1, row[row$beta_arm == "C2", , drop = FALSE],
             output_root, smoke)
    c1_id <- row$data_id[row$beta_arm == "C1"]
    c2_id <- row$data_id[row$beta_arm == "C2"]
    obs1 <- readRDS(file.path(output_root, "data", c1_id,
                              "observation_bundle.rds"))
    obs2 <- readRDS(file.path(output_root, "data", c2_id,
                              "observation_bundle.rds"))
    c2_checks <- c(
      c1_c2_identical_Y = identical(obs1$Y, obs2$Y),
      c1_c2_identical_time = identical(obs1$time_obs, obs2$time_obs),
      c1_has_d2_Z = is.list(obs1$Z) && obs1$dimensions$d == 2L,
      c2_omits_Z = is.null(obs2$Z) && obs2$dimensions$d == 0L)
    checks[[length(checks) + 1L]] <- data.frame(
      data_index = row$data_index[[1L]], check_id = names(c2_checks),
      passed = unname(c2_checks), logical_data_id = logical_id,
      stringsAsFactors = FALSE)
  }
  qc_all <- do.call(rbind, checks)
  assert(all(qc_all$passed), "Stage 7C cross-arm pairing QC failed.")
  csv_new(qc_all, file.path(output_root, "PAIRING_QC.csv"))
  lines_new(c("status=DATA_GENERATION_AND_SEAL_COMPLETE",
              paste0("logical_datasets=",
                     length(unique(data_rows$logical_data_id))),
              paste0("observation_bundles=", nrow(data_rows)),
              paste0("planned_fits=", nrow(fit_rows)),
              "truth_unseal_authorized=FALSE",
              "scientific_evaluation_started=FALSE"),
            file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("GENERATE_PASS data=", nrow(data_rows), " fits=", nrow(fit_rows),
      " smoke=", smoke, "\n", sep = "")
  invisible(TRUE)
}

benchmark <- function(output_root) {
  load_package()
  assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
         "Generate and seal before benchmarking.")
  benchmark_root <- file.path(output_root, "preflight")
  assert(!dir.exists(benchmark_root), "Refusing benchmark rerun.")
  dir.create(benchmark_root, recursive = TRUE, mode = "0700")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  row <- manifest[manifest$logical_data_id == "covrb7c_01" &
                    manifest$beta_arm == "C1" &
                    manifest$start_index == 1L, , drop = FALSE]
  assert(nrow(row) == 1L, "Stage 7C benchmark row missing.")
  observation <- readRDS(file.path(
    output_root, "data", row$data_id[[1L]], "observation_bundle.rds"))
  assert(observation$dimensions$p == 500L &&
           observation$dimensions$d == 2L &&
           observation$dimensions$K == 9L,
         "Benchmark must use full-size covariate data.")
  short <- capture_fit(observation, row, "benchmark3")
  long <- capture_fit(observation, row, "benchmark8")
  assert(is.null(short$error) && is.null(long$error), paste0(
    "Full-size benchmark failed: ", short$error %||% long$error))
  assert(short$fit$t1_sweeps == 3L && long$fit$t1_sweeps == 8L &&
           isTRUE(all.equal(short$fit$ELBO, long$fit$ELBO[1:3],
                            tolerance = 1e-10)),
         "Benchmark ordinary T=1 trajectory is not reproducible.")
  marginal <- (long$elapsed_seconds - short$elapsed_seconds) / 5
  assert(is.finite(marginal) && marginal > 0,
         "Benchmark cannot estimate a positive marginal sweep time.")
  result <- data.frame(
    scenario = row$data_id[[1L]], p = 500L, d = 2L, K = 9L,
    short_t1_sweeps = 3L, long_t1_sweeps = 8L,
    short_elapsed_seconds = short$elapsed_seconds,
    long_elapsed_seconds = long$elapsed_seconds,
    estimated_marginal_seconds_per_T1_sweep = marginal,
    projected_seconds_per_499_sweep_fit = 499 * marginal,
    projected_hours_108_fits_8_workers = 108 * 499 * marginal / 8 / 3600,
    stringsAsFactors = FALSE)
  csv_new(result, file.path(benchmark_root, "FULL_SIZE_SWEEP_BENCHMARK.csv"))
  lines_new(c("status=FULL_SIZE_SWEEP_BENCHMARK_COMPLETE",
              "main_fits_started=FALSE", "truth_unseal_authorized=FALSE",
              paste0("marginal_T1_seconds=", marginal)),
            file.path(benchmark_root, "BENCHMARK_COMPLETE.txt"))
  print(result)
  invisible(result)
}

batch <- function(output_root, workers) {
  assert(.Platform$OS.type == "unix", "Batch concurrency requires Linux.")
  assert(file.exists(file.path(output_root, "preflight",
                               "BENCHMARK_COMPLETE.txt")),
         "Full-size sweep benchmark must finish before batch fitting.")
  assert(is.finite(workers) && workers >= 1L && workers <= 8L,
         "Outer worker count must be from 1 through 8.")
  manifest <- read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  assert(nrow(manifest) == 108L && !anyDuplicated(manifest$fit_id) &&
           all(table(manifest$logical_data_id, manifest$beta_arm) == 12L),
         "Main Stage 7C manifest must contain 108 registered fits.")
  logs <- file.path(output_root, "supervisor_logs")
  dir.create(logs, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  run_one <- function(fit_id) {
    status <- system2(
      "Rscript",
      c(script_path, "--action=fit", paste0("--fit-id=", fit_id),
        paste0("--output-root=", output_root)),
      stdout = file.path(logs, paste0(fit_id, ".stdout.log")),
      stderr = file.path(logs, paste0(fit_id, ".stderr.log")))
    data.frame(fit_id = fit_id, exit_status = as.integer(status))
  }
  statuses <- parallel::mclapply(
    manifest$fit_id, run_one, mc.cores = workers,
    mc.preschedule = FALSE, mc.set.seed = FALSE)
  csv_new(do.call(rbind, statuses),
          file.path(output_root, "BATCH_PROCESS_STATUS.csv"))
  complete <- length(list.files(file.path(output_root, "fits"),
                                "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  errors <- length(list.files(file.path(output_root, "fits"),
                              "^FIT_ERROR[.]txt$", recursive = TRUE))
  if (complete + errors == 108L) {
    select_winners(output_root, FALSE)
    lines_new(c("status=MAIN_FITS_TERMINAL_AND_SELECTION_FROZEN",
                paste0("complete=", complete), paste0("errors=", errors),
                paste0("outer_workers=", workers), "per_fit_n_cpus=1",
                "truth_unseal_authorized=FALSE",
                "scientific_evaluation_started=FALSE"),
              file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  } else {
    lines_new(c("status=BATCH_INCOMPLETE_NO_SELECTIVE_RERUN",
                paste0("complete=", complete), paste0("errors=", errors)),
              file.path(output_root, "BATCH_INCOMPLETE.txt"))
    stop("Batch ended without all 108 terminal records.", call. = FALSE)
  }
  invisible(TRUE)
}

check_registration <- function() {
  load_package()
  assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
  d <- data_manifest(FALSE)
  f <- fit_manifest(d, FALSE)
  assert(nrow(d) == 9L && nrow(f) == 108L &&
           !anyDuplicated(d$data_id) && !anyDuplicated(f$fit_id) &&
           all(table(f$logical_data_id, f$beta_arm) == 12L) &&
           all(table(f$fit_seed) == 3L) &&
           all(f$n_cpus == 1L) && all(f$fit_K == 9L) &&
           all(f$planned_t1_sweeps == 400L) &&
           all(f$function_initialization == "gram_unit_energy") &&
           !any(f$truth_available_to_fit) &&
           !any(f$truth_available_to_stopping) &&
           !any(f$truth_available_to_selection),
         "Stage 7C candidate registration is invalid.")
  auth <- read_key_value(file.path(
    script_root, "STAGE7C_PROTOCOL_SEED_AND_EXECUTION_AUTHORIZATION_20260927.txt"))
  digest_file <- function(name) digest::digest(
    file.path(script_root, name), algo = "sha256", file = TRUE,
    serialize = FALSE)
  assert(identical(auth$protocol_confirmed, "TRUE") &&
           identical(auth$seed_schedule_confirmed, "TRUE") &&
           identical(auth$stage7c_108_fit_execution_after_acceptance_authorized,
                     "TRUE") &&
           identical(auth$truth_unseal_authorized, "FALSE") &&
           identical(auth$scientific_evaluation_authorized, "FALSE") &&
           identical(digest_file(auth$bound_protocol_file),
                     auth$bound_protocol_sha256) &&
           identical(digest_file(auth$bound_data_manifest_file),
                     auth$bound_data_manifest_sha256) &&
           identical(digest_file(auth$bound_fit_manifest_file),
                     auth$bound_fit_manifest_sha256),
         "Stage 7C authorization or bound hashes are invalid.")
  cat("CHECK_PASS data=9 fits=108 package=",
      as.character(utils::packageVersion("multiFSYNC")), "\n", sep = "")
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  check_registration()
} else if (action == "generate") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, FALSE)
} else if (action == "benchmark") {
  assert(nzchar(output_root), "--output-root is required.")
  benchmark(output_root)
} else if (action == "fit") {
  assert(nzchar(output_root) && nzchar(cli$`fit-id` %||% ""),
         "--output-root and --fit-id are required.")
  fit_one(output_root, cli$`fit-id`, identical(cli$smoke, "true"))
} else if (action == "select") {
  assert(nzchar(output_root), "--output-root is required.")
  select_winners(output_root, identical(cli$smoke, "true"))
} else if (action == "batch") {
  assert(nzchar(output_root), "--output-root is required.")
  batch(output_root, as.integer(Sys.getenv(
    "COVARIATE_ROUTEB_STAGE7C_OUTER_WORKERS", unset = "8")))
} else if (action == "smoke") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, TRUE)
  rows <- fit_manifest(data_manifest(TRUE), TRUE)
  for (fit_id in rows$fit_id) fit_one(output_root, fit_id, TRUE)
  select_winners(output_root, TRUE)
  lines_new(c("status=STAGE7C_MICRO_SMOKE_PASS", "generated_arms=3",
              "fits=3", "formal_candidate_seeds_used=FALSE",
              "truth_unseal_authorized=FALSE"),
            file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
  cat("STAGE7C_MICRO_SMOKE_PASS\n")
} else if (action == "status") {
  assert(nzchar(output_root), "--output-root is required.")
  cat("complete=", length(list.files(file.path(output_root, "fits"),
                                      "^FIT_COMPLETE[.]txt$", recursive = TRUE)),
      " errors=", length(list.files(file.path(output_root, "fits"),
                                     "^FIT_ERROR[.]txt$", recursive = TRUE)),
      " selected=", file.exists(file.path(
        output_root, "truth_free_selection",
        "TRUTH_FREE_SELECTION_FROZEN.txt")), "\n", sep = "")
} else stop("Unknown action: ", action, call. = FALSE)
