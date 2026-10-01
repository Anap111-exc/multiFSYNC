#!/usr/bin/env Rscript

# Stage 7I truth-blind A/S continuation. The U parity gate is mandatory.
# No truth bundle is opened, and no frozen Stage 7C--7H file is written.
# Usage: Rscript this.R {prepare|fit|resume_s|status} U_ROOT OUTPUT_ROOT
# `fit` runs A and then S. `resume_s` is only for an interrupted process after
# an independently complete A; it never reruns either branch or overwrites data.
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 3L || !argv[[1L]] %in%
    c("prepare", "fit", "resume_s", "status")) {
  stop("Usage: Rscript run_stage7i_as_resume_20260929_v1.R ",
       "{prepare|fit|resume_s|status} U_ROOT OUTPUT_ROOT", call. = FALSE)
}
action <- argv[[1L]]
u_root <- normalizePath(argv[[2L]], mustWork = TRUE)
output_root <- argv[[3L]]
file_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_option) != 1L) stop("Cannot identify A/S runner source.")
this_script <- normalizePath(sub("^--file=", "", file_option))
helper_path <- file.path(dirname(this_script),
                         "build_stage7i_uas_states_20260929_v1.R")
protocol_path <- file.path(dirname(this_script),
                           "PROTOCOL_STAGE7I_UAS_20260929_v1.md")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
sha256 <- function(path) {
  assert(file.exists(path) && !dir.exists(path),
         paste0("Missing bound file: ", path))
  lines <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  assert(length(lines) == 1L &&
           grepl("^[0-9a-f]{64}[[:space:]]", lines[[1L]]),
         paste0("sha256sum failed: ", path))
  substr(lines[[1L]], 1L, 64L)
}
new_dir <- function(path) {
  assert(!file.exists(path) && !dir.exists(path),
         paste0("Refusing existing output directory: ", path))
  assert(dir.create(path, recursive = TRUE, showWarnings = FALSE),
         paste0("Cannot create output directory: ", path))
}
new_csv <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
new_rds <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  saveRDS(x, path, version = 3L)
}
new_lines <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  writeLines(x, path, useBytes = TRUE)
}
all_finite <- function(x) {
  if (is.null(x) || !length(x)) return(TRUE)
  if (is.list(x)) return(all(vapply(x, all_finite, logical(1L))))
  if (is.numeric(x)) return(all(is.finite(x)))
  if (is.logical(x)) return(all(!is.na(x)))
  FALSE
}
as_binding <- function(paths) {
  data.frame(role = names(paths),
             path = unname(vapply(paths, normalizePath, character(1L),
                                  mustWork = TRUE)),
             sha256 = unname(vapply(paths, sha256, character(1L))),
             stringsAsFactors = FALSE)
}

u_binding_path <- file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv")
u_pass_path <- file.path(u_root, "u_parity_v1", "PARITY_PASS.txt")
u_parity_files <- file.path(u_root, "u_parity_v1", c(
  "ALL_STATE_PARITY.csv", "ORDINARY_ELBO_PARITY.csv",
  "ENDPOINT_PARITY.csv", "CONFIGURATION_PARITY.csv",
  "DRIVER_TRACE_PARITY.csv"))
u_fit_path <- file.path(u_root, "u_full_v1", "fit.rds")
source_binding_path <- file.path(output_root, "SOURCE_BINDING.csv")
injected_binding_path <- file.path(output_root, "INJECTED_STATE_BINDING.csv")

read_u_binding <- function() {
  assert(file.exists(u_binding_path) && file.exists(u_pass_path),
         "U input binding or parity pass marker is missing.")
  frozen <- utils::read.csv(u_binding_path, stringsAsFactors = FALSE)
  assert(identical(names(frozen), c("role", "path", "sha256")) &&
           !anyDuplicated(frozen$role) && nrow(frozen) >= 14L &&
           all(grepl("^[0-9a-f]{64}$", frozen$sha256)),
         "U input binding schema is invalid.")
  required <- c("runner", "stage7i_archive", "stage7i_description",
                "stage7i_namespace", "stage7i_rdb", "stage7i_rdx",
                "stage7g_runner", "stage7g_archive", "stage7g_manifest",
                "stage7g_input_binding", "stage7h_blind_manifest",
                "observation", "stage7g_fit", "stage7g_snapshots")
  assert(setequal(frozen$role, required),
         "U source binding roles have changed.")
  expected_frozen <- c(
    stage7i_archive =
      "b6c2210089061d99e31825e358d9c03b6d867bcd7fdde21e8432e784e4436178",
    stage7g_runner =
      "668788f6cb159badedb15c0630d0754019e2de18aaaa79cf61bff60eb3acd935",
    stage7g_archive =
      "b3e9050f1de7154baff28b7dd986254df0b634bf85d6f50f464d00dbd4c8f17f",
    stage7g_manifest =
      "20a68304cccb3f4176c046b4c5c4f4deb4ca81ffd6d10c99ccc5dfdd20ad4e7a",
    stage7g_input_binding =
      "1d0e6ae8badf33edbb9b2578e615d395f00f5501a7aad2fce3db0f456ce62a54",
    stage7h_blind_manifest =
      "d67bef6e692fa87296f3d0b0e19a6a954d05a0085bd3f3ec4b08350b3f0f3f98",
    observation =
      "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6",
    stage7g_fit =
      "b2166588d932c4432261e026745f1d1d8185dc0b00f741125d6a779de5d7d761",
    stage7g_snapshots =
      "e85ce3fdd8b2770674232e11cdab5453f61b30271977a41cb08d54b6383621c2")
  assert(all(vapply(names(expected_frozen), function(role)
    identical(frozen$sha256[frozen$role == role],
              unname(expected_frozen[[role]])), logical(1L))),
    "U binding no longer matches the frozen Stage 7G/H input hashes.")
  assert(all(vapply(seq_len(nrow(frozen)), function(i)
    identical(sha256(frozen$path[[i]]), frozen$sha256[[i]]), logical(1L))),
    "A file in U's frozen input/source binding has changed.")
  marker <- readLines(u_pass_path, warn = FALSE)
  assert(length(marker) >= 1L &&
           identical(marker[[1L]], "STAGE7I_U_PARITY_PASS"),
         "U parity marker is not a pass.")
  check_csv <- function(path, expected_n = NULL) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    assert("pass" %in% names(x) && nrow(x) > 0L &&
             all(!is.na(x$pass)) && all(as.logical(x$pass)) &&
             (is.null(expected_n) || nrow(x) == expected_n),
           paste0("U parity table failed or changed: ", path))
    x
  }
  state <- check_csv(u_parity_files[[1L]])
  assert(identical(sort(unique(as.integer(state$iteration))),
                   as.integer(51:499)),
         "U state parity must cover every resumed sweep 51:499.")
  check_csv(u_parity_files[[2L]], 400L)
  check_csv(u_parity_files[[3L]])
  check_csv(u_parity_files[[4L]])
  check_csv(u_parity_files[[5L]], 1L)
  assert(file.exists(u_fit_path), "U full continuation fit is missing.")
  stats::setNames(frozen$path, frozen$role)
}
source_paths <- function() {
  read_u_binding()
  c(as_runner = this_script, initializer = helper_path,
    protocol = protocol_path, u_binding = u_binding_path,
    u_parity_pass = u_pass_path,
    u_state_parity = u_parity_files[[1L]],
    u_elbo_parity = u_parity_files[[2L]],
    u_endpoint_parity = u_parity_files[[3L]],
    u_config_parity = u_parity_files[[4L]],
    u_trace_parity = u_parity_files[[5L]],
    u_fit = u_fit_path)
}
check_prepare <- function() {
  assert(file.exists(file.path(output_root, "INJECTION_COMPLETE.txt")),
         "A/S injection completion marker is missing.")
  previous <- utils::read.csv(source_binding_path,
                              stringsAsFactors = FALSE)
  current <- as_binding(source_paths())
  assert(identical(previous, current),
         "A/S helper, protocol, U result or source binding changed.")
  injected <- utils::read.csv(injected_binding_path,
                              stringsAsFactors = FALSE)
  expected <- as_binding(c(
    U = file.path(output_root, "STATE_U_SWEEP50.rds"),
    A = file.path(output_root, "STATE_A_SWEEP50.rds"),
    S = file.path(output_root, "STATE_S_SWEEP50.rds"),
    diagnostics = file.path(output_root,
                            "BLIND_INITIALIZATION_DIAGNOSTICS.rds"),
    summary = file.path(output_root, "BLIND_INITIALIZATION_SUMMARY.csv"),
    driver_coordinates = file.path(output_root,
                                   "BLIND_FIXED_RESIDUAL_DRIVER.csv"),
    driver_subjects = file.path(output_root,
                                "BLIND_SUBJECT_DRIVER.csv")))
  assert(identical(injected, expected),
         "Injected states or blind diagnostics changed after freeze.")
  invisible(TRUE)
}
load_bound_package <- function(u_paths) {
  package_description <- u_paths[["stage7i_description"]]
  library_root <- dirname(dirname(package_description))
  .libPaths(c(library_root, .libPaths()))
  suppressPackageStartupMessages(library(multiFSYNC))
  assert(identical(normalizePath(system.file(package = "multiFSYNC")),
                   normalizePath(dirname(package_description))) &&
           identical(as.character(utils::packageVersion("multiFSYNC")),
                     "0.3.0.9000"),
         "Loaded package differs from the U-bound Stage 7I build.")
}
read_observation <- function(u_paths) {
  x <- readRDS(u_paths[["observation"]])
  assert(identical(x$bundle_class,
                   "covariate_routeb_pilot_observation_only") &&
           isFALSE(x$truth_available_to_fit) &&
           identical(x$data_id, "covrb7c_01__C1") &&
           x$dimensions$p == 500L && x$dimensions$d == 2L &&
           x$dimensions$K == 9L,
         "Observation-only contract changed.")
  x
}
stage7i_fit_args <- function(observation) {
  control <- multiFSYNC::g12_stopping_control()
  control$checkpoints <- c(380L, 400L)
  list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
    K = observation$dimensions$K, n_g = 51L,
    anneal = c(1, 1.9, 100), maxit = 499L,
    n_cpus = 1L, verbose = FALSE, seed = 99716106L,
    tol_abs = 1e-3, tol_rel = 1e-5,
    bool_scale = FALSE, bool_var_spec_prob = FALSE,
    d_0 = as.integer(observation$dimensions$p),
    convergence_rule = "practical", lambda_orth = 0,
    practical_control = control, initialization = "random",
    initialization_control = list(grid_size = 81L,
                                  perturb_sd = 0.05, rank_tol = 1e-8),
    continuation_state = NULL, pre_score_sweeps = 1L,
    trace_sweeps = as.integer(51:499),
    function_initialization = "gram_unit_energy")
}
scale_trace_control <- list(
  include_initial = FALSE, include_annealing = TRUE,
  t1_sweeps = as.integer(c(1:20, 50, 100, 200, 300, 380, 400)),
  block_annealing_sweeps = as.integer(51:60)
)
new_snapshot_sink <- function() {
  x <- new.env(parent = emptyenv())
  x$schema <- "stage7g_snapshot_v1"
  x$sweeps <- as.integer(51:499)
  x$records <- list()
  x
}
check_branch_complete <- function(branch) {
  root <- file.path(output_root, paste0(branch, "_full_v1"))
  marker_path <- file.path(root, "FIT_COMPLETE.txt")
  binding_path <- file.path(root, "BRANCH_OUTPUT_BINDING.csv")
  assert(file.exists(marker_path) && file.exists(binding_path),
         paste0(branch, " completion marker or output binding is missing."))
  marker <- readLines(marker_path, warn = FALSE)
  assert(length(marker) >= 2L &&
           identical(marker[[1L]], "STAGE7I_AS_FIT_COMPLETE") &&
           paste0("branch=", branch) %in% marker &&
           paste0("binding_sha256=", sha256(binding_path)) %in% marker,
         paste0(branch, " completion marker is invalid."))
  previous <- utils::read.csv(binding_path, stringsAsFactors = FALSE)
  assert(identical(names(previous), c("role", "path", "sha256")) &&
           !anyDuplicated(previous$role) && nrow(previous) >= 6L,
         paste0(branch, " output binding schema is invalid."))
  current <- as_binding(stats::setNames(previous$path, previous$role))
  assert(identical(previous, current),
         paste0(branch, " completed output changed."))
  invisible(TRUE)
}

if (action == "prepare") {
  u_paths <- read_u_binding()
  proposed <- normalizePath(output_root, mustWork = FALSE)
  protected_roots <- c(
    u_root,
    dirname(dirname(dirname(u_paths[["observation"]]))),
    dirname(dirname(dirname(u_paths[["stage7g_fit"]]))),
    dirname(dirname(u_paths[["stage7h_blind_manifest"]])),
    dirname(dirname(u_paths[["stage7i_description"]])),
    dirname(this_script))
  assert(!any(vapply(protected_roots, function(root)
    identical(proposed, root) ||
      startsWith(proposed, paste0(root, .Platform$file.sep)),
    logical(1L))) &&
           !file.exists(output_root) && !dir.exists(output_root),
         "A/S output root overlaps a bound input or already exists.")
  pre_binding <- as_binding(source_paths())
  load_bound_package(u_paths)
  observation <- read_observation(u_paths)
  fit <- readRDS(u_paths[["stage7g_fit"]])
  source_states <- readRDS(u_paths[["stage7g_snapshots"]])
  assert(length(source_states) == 500L &&
           identical(vapply(source_states, `[[`, integer(1L), "iteration"),
                     as.integer(0:499)),
         "Frozen Stage 7G states must cover 0:499.")
  snapshot50 <- source_states[[51L]]
  rm(source_states)
  invisible(gc(verbose = FALSE))
  grid <- multiFSYNC:::get_grid_objects(
    unlist(observation$time_obs, recursive = FALSE), fit$K,
    n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
  n1 <- fit$n_s[[1L]]
  n2 <- fit$n_s[[2L]]
  C2 <- grid$C[(n1 + 1L):(n1 + n2)]
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helper_env)
  starts <- helper_env$stage7i_build_uas_states(
    observation, fit, snapshot50, C2)
  assert(identical(starts$U, snapshot50) &&
           identical(names(starts), c("U", "A", "S", "blind_diagnostics")) &&
           isTRUE(starts$blind_diagnostics$amplitudes_matched),
         "Initializer did not preserve U or match the fixed A/S amplitudes.")
  assert(identical(pre_binding, as_binding(source_paths())),
         "A/S source or U parity changed during injection construction.")

  # All construction finishes before the first output path is created.
  new_dir(output_root)
  new_csv(pre_binding, source_binding_path)
  for (branch in c("U", "A", "S")) {
    new_rds(starts[[branch]], file.path(output_root,
                                       paste0("STATE_", branch,
                                              "_SWEEP50.rds")))
  }
  details <- starts$blind_diagnostics
  new_rds(details, file.path(output_root,
                             "BLIND_INITIALIZATION_DIAGNOSTICS.rds"))
  energies <- details$initial_energy
  summary <- data.frame(
    branch = rownames(energies), energies,
    fixed_residual_driver_norm = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L), "h_norm"),
    normalized_fixed_residual_driver = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L), "normalized_h"),
    v_signed_driver = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L),
      "v_signed_driver"),
    v_abs_cosine_to_driver = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L),
      "v_abs_cosine_to_driver"),
    v_subject_cancellation_ratio = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L),
      "v_subject_cancellation_ratio"),
    subject_driver_cancellation = vapply(
      details$fixed_residual_driver, `[[`, numeric(1L),
      "global_subject_cancellation_ratio"),
    stringsAsFactors = FALSE)
  new_csv(summary, file.path(output_root,
                             "BLIND_INITIALIZATION_SUMMARY.csv"))
  coordinate <- do.call(rbind, lapply(c("U", "A", "S"), function(branch) {
    d <- details$fixed_residual_driver[[branch]]
    data.frame(branch = branch, coordinate = seq_along(d$h),
               h = d$h, signed_envelope = d$signed_envelope,
               coherence = d$per_coordinate_coherence)
  }))
  new_csv(coordinate, file.path(output_root,
                                "BLIND_FIXED_RESIDUAL_DRIVER.csv"))
  subject <- do.call(rbind, lapply(c("U", "A", "S"), function(branch) {
    d <- details$fixed_residual_driver[[branch]]
    data.frame(branch = branch,
               subject = seq_along(d$subject_driver_norms),
               driver_norm = d$subject_driver_norms,
               v_signed_contribution = d$v_subject_signed,
               v_absolute_contribution = abs(d$v_subject_signed),
               v_subject_cancellation_ratio =
                 d$v_subject_cancellation_ratio)
  }))
  new_csv(subject, file.path(output_root, "BLIND_SUBJECT_DRIVER.csv"))
  state_paths <- c(
    U = file.path(output_root, "STATE_U_SWEEP50.rds"),
    A = file.path(output_root, "STATE_A_SWEEP50.rds"),
    S = file.path(output_root, "STATE_S_SWEEP50.rds"),
    diagnostics = file.path(output_root,
                            "BLIND_INITIALIZATION_DIAGNOSTICS.rds"),
    summary = file.path(output_root, "BLIND_INITIALIZATION_SUMMARY.csv"),
    driver_coordinates = file.path(output_root,
                                   "BLIND_FIXED_RESIDUAL_DRIVER.csv"),
    driver_subjects = file.path(output_root,
                                "BLIND_SUBJECT_DRIVER.csv"))
  new_csv(as_binding(state_paths), injected_binding_path)
  new_lines(c(
    "STAGE7I_AS_INJECTION_COMPLETE", "seed01_C1_start06_only=TRUE",
    "source_sweep=50", "truth_opened=FALSE",
    "formal_candidate_pool=FALSE",
    paste0("helper_sha256=", sha256(helper_path)),
    paste0("u_parity_sha256=", sha256(u_pass_path))
  ), file.path(output_root, "INJECTION_COMPLETE.txt"))
  cat("STAGE7I_AS_INJECTION_COMPLETE\n")
}

run_branch <- function(branch) {
  assert(branch %in% c("A", "S"), "Only A and S are accepted.")
  check_prepare()
  if (branch == "S") {
    check_branch_complete("A")
  }
  u_paths <- read_u_binding()
  load_bound_package(u_paths)
  observation <- read_observation(u_paths)
  injected <- readRDS(file.path(output_root,
                                paste0("STATE_", branch,
                                       "_SWEEP50.rds")))
  assert(identical(injected$iteration, 50L) &&
           isTRUE(injected$annealing),
         paste0("Injected ", branch, " snapshot is not sweep 50."))
  result_dir <- file.path(output_root, paste0(branch, "_full_v1"))
  new_dir(result_dir)
  resume <- new.env(parent = emptyenv())
  resume$schema <- "stage7i_sweep50_resume_v1"
  resume$state <- injected
  resume$stop_after_iteration <- 499L
  sink <- new_snapshot_sink()
  old_options <- options(
    multiFSYNC.__stage7i_resume_20260929__ = resume,
    multiFSYNC.__stage7g_snapshot_20260928__ = sink,
    multiFSYNC.__private_scale_trace_8d72c46f__ = scale_trace_control)
  on.exit(options(old_options), add = TRUE)
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  outcome <- tryCatch(withCallingHandlers(
    do.call(multiFSYNC::bayesSYNC_multi_pre_score,
            stage7i_fit_args(observation)),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))
  ), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  new_lines(warnings, file.path(result_dir, "WARNINGS.txt"))
  new_rds(sink$records,
          file.path(result_dir, "RAW_STATE_SNAPSHOTS_51_499.rds"))
  if (length(sink$records)) {
    new_csv(data.frame(
      iteration = vapply(sink$records, `[[`, integer(1L), "iteration"),
      temperature = vapply(sink$records, `[[`, numeric(1L), "temperature"),
      annealing = vapply(sink$records, `[[`, logical(1L), "annealing")
    ), file.path(result_dir, "RAW_STATE_INDEX.csv"))
  }
  if (inherits(outcome, "error")) {
    new_lines(c("STAGE7I_AS_FIT_ERROR", paste0("branch=", branch),
                conditionMessage(outcome)),
              file.path(result_dir, "FIT_ERROR.txt"))
    stop(paste0(branch, " continuation failed; partial trace retained."))
  }
  fit <- outcome
  new_rds(fit, file.path(result_dir, "fit.rds"))
  new_csv(fit$driver_trace, file.path(result_dir, "DRIVER_TRACE.csv"))
  if (!is.null(fit$scale_trace)) {
    new_csv(fit$scale_trace, file.path(result_dir, "SCALE_TRACE.csv"))
  }
  new_csv(data.frame(t1_sweep = seq_along(fit$ELBO),
                     ordinary_elbo = as.numeric(fit$ELBO)),
          file.path(result_dir, "ORDINARY_ELBO.csv"))
  valid <- isTRUE(resume$applied) &&
    length(sink$records) == 449L &&
    all_finite(sink$records) &&
    identical(vapply(sink$records, `[[`, integer(1L), "iteration"),
              as.integer(51:499)) &&
    fit$i_iter == 499L && fit$annealing_sweeps == 99L &&
    fit$t1_sweeps == 400L && length(fit$ELBO) == 400L &&
    all(is.finite(fit$ELBO)) &&
    isTRUE(fit$elbo_objective_valid) &&
    identical(as.integer(fit$fit_seed), 99716106L)
  if (!valid) {
    new_lines(c("STAGE7I_AS_FIT_CONTRACT_FAILED",
                paste0("branch=", branch)),
              file.path(result_dir, "FIT_CONTRACT_FAILED.txt"))
    stop(paste0(branch, " violated the fixed 449-sweep contract."))
  }
  output_paths <- c(
    fit = file.path(result_dir, "fit.rds"),
    snapshots = file.path(result_dir, "RAW_STATE_SNAPSHOTS_51_499.rds"),
    state_index = file.path(result_dir, "RAW_STATE_INDEX.csv"),
    driver_trace = file.path(result_dir, "DRIVER_TRACE.csv"),
    elbo = file.path(result_dir, "ORDINARY_ELBO.csv"),
    warnings = file.path(result_dir, "WARNINGS.txt"))
  if (file.exists(file.path(result_dir, "SCALE_TRACE.csv"))) {
    output_paths <- c(output_paths,
                      scale_trace = file.path(result_dir, "SCALE_TRACE.csv"))
  }
  output_binding_path <- file.path(result_dir, "BRANCH_OUTPUT_BINDING.csv")
  new_csv(as_binding(output_paths), output_binding_path)
  new_lines(c(
    "STAGE7I_AS_FIT_COMPLETE", paste0("branch=", branch),
    "source_sweep=50", "resumed_sweeps=449", "ordinary_sweeps=400",
    "truth_opened=FALSE", "formal_candidate_pool=FALSE",
    paste0("final_ordinary_elbo=",
           format(tail(fit$ELBO, 1L), digits = 16L)),
    paste0("binding_sha256=", sha256(output_binding_path)),
    paste0("elapsed_seconds=", format(elapsed, digits = 10L))
  ), file.path(result_dir, "FIT_COMPLETE.txt"))
  cat("STAGE7I_AS_FIT_COMPLETE branch=", branch, "\n", sep = "")
  invisible(TRUE)
}

if (action == "fit") {
  check_prepare()
  assert(!dir.exists(file.path(output_root, "A_full_v1")) &&
           !dir.exists(file.path(output_root, "S_full_v1")),
         "A/S fit outputs already exist; refusing automatic rerun.")
  run_branch("A")
  run_branch("S")
  new_lines(c("STAGE7I_AS_FITS_COMPLETE", "branches=A,S",
              "truth_opened=FALSE", "formal_candidate_pool=FALSE"),
            file.path(output_root, "AS_FITS_COMPLETE.txt"))
  cat("STAGE7I_AS_FITS_COMPLETE\n")
}

if (action == "resume_s") {
  check_prepare()
  check_branch_complete("A")
  assert(!dir.exists(file.path(output_root, "S_full_v1")) &&
           !file.exists(file.path(output_root, "AS_FITS_COMPLETE.txt")),
         "Manual S continuation requires completed A and absent S.")
  run_branch("S")
  new_lines(c("STAGE7I_AS_FITS_COMPLETE", "branches=A,S",
              "manual_S_resume=TRUE", "truth_opened=FALSE",
              "formal_candidate_pool=FALSE"),
            file.path(output_root, "AS_FITS_COMPLETE.txt"))
  cat("STAGE7I_AS_FITS_COMPLETE manual_S_resume=TRUE\n")
}

if (action == "status") {
  cat("injection=", file.exists(file.path(output_root,
                                          "INJECTION_COMPLETE.txt")),
      " A_complete=", file.exists(file.path(output_root, "A_full_v1",
                                           "FIT_COMPLETE.txt")),
      " S_complete=", file.exists(file.path(output_root, "S_full_v1",
                                           "FIT_COMPLETE.txt")),
      " all_complete=", file.exists(file.path(output_root,
                                              "AS_FITS_COMPLETE.txt")),
      " A_error=", file.exists(file.path(output_root, "A_full_v1",
                                        "FIT_ERROR.txt")),
      " S_error=", file.exists(file.path(output_root, "S_full_v1",
                                        "FIT_ERROR.txt")), "\n", sep = "")
}
