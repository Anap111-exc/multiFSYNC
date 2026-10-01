#!/usr/bin/env Rscript
# One read-only-input, truth-blind continuation of the frozen Stage 7J S state.
# Usage: Rscript SCRIPT {preflight|fit|status} RESULTS_PARENT OUTPUT_ROOT
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 3L || !argv[[1L]] %in% c("preflight", "fit", "status"))
  stop("Expected {preflight|fit|status} RESULTS_PARENT OUTPUT_ROOT.", call. = FALSE)
action <- argv[[1L]]
parent <- normalizePath(argv[[2L]], mustWork = TRUE)
out <- argv[[3L]]
script_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_option) != 1L) stop("Cannot identify runner source.")
this_script <- normalizePath(sub("^--file=", "", script_option))
protocol <- file.path(dirname(this_script), "PROTOCOL_STAGE7K_SEED02_S_20260930_v1.md")
u_root <- file.path(parent, "stage7j_seed02_u_20260930_v1")
j_root <- file.path(parent, "stage7j_seed02_a_20260930_v1")
freeze_root <- file.path(parent, "stage7j_seed02_blind_freeze_20260930_v1")
expected_s_sha <- "72a19fd46ed9f096966eca6c5b6f9bd673a6010e51a356ddcaa7e5156ad6d292"
selected <- as.integer(c(51:70, 75, 80, 90, 100, 120, 150, 200, 250,
                         300, 350, 400, 450, 499))
scale_control <- list(include_initial = FALSE, include_annealing = TRUE,
                      t1_sweeps = as.integer(c(1:20, 50, 100, 200, 300,
                                               380, 400)),
                      block_annealing_sweeps = as.integer(51:60))
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
sha <- function(path) {
  assert(file.exists(path) && !dir.exists(path), paste("Missing", path))
  lines <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  assert(length(lines) == 1L && grepl("^[0-9a-f]{64}[[:space:]]", lines),
         paste("SHA256 failed", path))
  substr(lines, 1L, 64L)
}
check <- function(path, expected, label) {
  assert(identical(sha(path), expected), paste(label, "SHA256 changed"))
  invisible(TRUE)
}
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                           check.names = FALSE)
new_csv <- function(x, path) {
  assert(!file.exists(path), paste("Refusing overwrite", path))
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}
new_rds <- function(x, path) {
  assert(!file.exists(path), paste("Refusing overwrite", path))
  saveRDS(x, path, version = 3L)
}
new_lines <- function(x, path) {
  assert(!file.exists(path), paste("Refusing overwrite", path))
  writeLines(x, path, useBytes = TRUE)
}
all_finite <- function(x) {
  if (is.null(x) || !length(x)) return(TRUE)
  if (is.list(x)) return(all(vapply(x, all_finite, logical(1L))))
  if (is.numeric(x)) return(all(is.finite(x)))
  if (is.logical(x)) return(all(!is.na(x)))
  FALSE
}
bind <- function(paths) data.frame(
  role = names(paths),
  path = unname(vapply(paths, normalizePath, character(1L), mustWork = TRUE)),
  sha256 = unname(vapply(paths, sha, character(1L))),
  stringsAsFactors = FALSE)
check_binding <- function(path, expected_roles = NULL) {
  x <- read_csv(path)
  assert(identical(names(x), c("role", "path", "sha256")) &&
           !anyDuplicated(x$role) && nrow(x) > 0L &&
           (is.null(expected_roles) || all(expected_roles %in% x$role)),
         paste("Invalid binding", path))
  for (i in seq_len(nrow(x))) check(x$path[[i]], x$sha256[[i]], x$role[[i]])
  stats::setNames(x$path, x$role)
}
preflight <- function() {
  assert(file.exists(protocol), "Protocol is missing.")
  proposed <- normalizePath(out, mustWork = FALSE)
  assert(identical(dirname(proposed), parent) &&
           identical(basename(proposed), "stage7k_seed02_s_20260930_v1") &&
           !file.exists(out) && !dir.exists(out),
         "Output root must be a new dedicated Stage 7K directory.")
  j_marker <- readLines(file.path(freeze_root, "BLIND_ENDPOINTS_FROZEN.txt"),
                        warn = FALSE)
  assert(all(c("STAGE7J_BLIND_ENDPOINTS_FROZEN", "branches=U,A",
               "truth_bundle_opened=FALSE", "formal_candidate_pool=FALSE") %in%
               j_marker), "Stage 7J blind freeze changed.")
  j_manifest_path <- file.path(freeze_root, "BLIND_ENDPOINT_MANIFEST.csv")
  expected_manifest <- sub("^manifest_sha256=", "",
                           j_marker[startsWith(j_marker, "manifest_sha256=")])
  assert(length(expected_manifest) == 1L, "Stage 7J manifest hash absent.")
  check(j_manifest_path, expected_manifest, "Stage 7J manifest")
  j_manifest <- read_csv(j_manifest_path)
  assert(identical(j_manifest$branch, c("U", "A")) &&
           all(j_manifest$fit_id == "covrb7c_02__C1__start11"),
         "Stage 7J endpoint set changed.")
  for (i in seq_len(nrow(j_manifest))) {
    check(j_manifest$fit_path[[i]], j_manifest$fit_sha256[[i]],
          paste("Stage 7J", j_manifest$branch[[i]], "fit"))
    check(j_manifest$snapshots_path[[i]],
          j_manifest$snapshots_sha256[[i]],
          paste("Stage 7J", j_manifest$branch[[i]], "snapshots"))
  }
  u_marker <- readLines(file.path(u_root, "u_parity_v1", "PARITY_PASS.txt"),
                        warn = FALSE)
  assert(identical(u_marker[[1L]], "STAGE7J_U_PARITY_PASS"),
         "U parity marker changed.")
  u_paths <- check_binding(file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv"),
                           c("observation", "stage7i_description"))
  injection_marker <- readLines(file.path(j_root, "INJECTION_COMPLETE.txt"),
                                warn = FALSE)
  assert(all(c("STAGE7J_A_INJECTION_COMPLETE", "truth_opened=FALSE") %in%
               injection_marker), "Stage 7J initialization changed.")
  check_binding(file.path(j_root, "SOURCE_BINDING.csv"))
  injected <- check_binding(file.path(j_root, "INJECTED_STATE_BINDING.csv"),
                            c("S"))
  s_path <- injected[["S"]]
  assert(identical(normalizePath(s_path),
                   normalizePath(file.path(j_root, "STATE_S_SWEEP50.rds"))) &&
           identical(sha(s_path), expected_s_sha),
         "Frozen S state path or hash differs.")
  state <- readRDS(s_path)
  assert(identical(state$iteration, 50L) && isTRUE(state$annealing) &&
           all_finite(state), "Frozen S state is invalid.")
  observation <- readRDS(u_paths[["observation"]])
  assert(identical(observation$data_id, "covrb7c_02__C1") &&
           isFALSE(observation$truth_available_to_fit) &&
           observation$dimensions$p == 500L && observation$dimensions$K == 9L,
         "Observation-only input changed.")
  library_root <- dirname(dirname(u_paths[["stage7i_description"]]))
  .libPaths(c(library_root, .libPaths()))
  suppressPackageStartupMessages(library(multiFSYNC))
  assert(identical(normalizePath(system.file(package = "multiFSYNC")),
                   normalizePath(dirname(u_paths[["stage7i_description"]]))),
         "Installed Stage 7I package changed.")
  paths <- c(runner = this_script, protocol = protocol,
             j_freeze = file.path(freeze_root, "BLIND_ENDPOINTS_FROZEN.txt"),
             j_manifest = j_manifest_path,
             j_injection = file.path(j_root, "INJECTED_STATE_BINDING.csv"),
             s_state = s_path,
             u_binding = file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv"),
             u_parity = file.path(u_root, "u_parity_v1", "PARITY_PASS.txt"),
             observation = u_paths[["observation"]],
             package_description = u_paths[["stage7i_description"]])
  list(state = state, observation = observation, paths = paths)
}
fit_args <- function(observation) {
  control <- multiFSYNC::g12_stopping_control()
  control$checkpoints <- c(380L, 400L)
  list(Y = observation$Y, Z = observation$Z,
       time_obs = observation$time_obs,
       L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
       K = observation$dimensions$K, n_g = 51L,
       anneal = c(1, 1.9, 100), maxit = 499L,
       n_cpus = 1L, verbose = FALSE, seed = 99716211L,
       tol_abs = 1e-3, tol_rel = 1e-5,
       bool_scale = FALSE, bool_var_spec_prob = FALSE,
       d_0 = as.integer(observation$dimensions$p),
       convergence_rule = "practical", lambda_orth = 0,
       practical_control = control, initialization = "random",
       initialization_control = list(grid_size = 81L,
                                     perturb_sd = 0.05, rank_tol = 1e-8),
       continuation_state = NULL, pre_score_sweeps = 1L,
       trace_sweeps = selected,
       function_initialization = "gram_unit_energy")
}
if (action == "preflight") {
  x <- preflight()
  cat("STAGE7K_PREFLIGHT_PASS s_sha256=", sha(x$paths[["s_state"]]),
      " fit_started=FALSE truth_opened=FALSE\n", sep = "")
}
if (action == "fit") {
  x <- preflight()
  assert(dir.create(out, mode = "0700"), "Cannot create Stage 7K output.")
  new_csv(bind(x$paths), file.path(out, "INPUT_AND_SOURCE_BINDING.csv"))
  result_dir <- file.path(out, "S_full_v1")
  assert(dir.create(result_dir, mode = "0700"), "Cannot create S output.")
  resume <- new.env(parent = emptyenv())
  resume$schema <- "stage7i_sweep50_resume_v1"
  resume$state <- x$state
  resume$stop_after_iteration <- 499L
  sink <- new.env(parent = emptyenv())
  sink$schema <- "stage7g_snapshot_v1"
  sink$sweeps <- selected
  sink$records <- list()
  old_options <- options(
    multiFSYNC.__stage7i_resume_20260929__ = resume,
    multiFSYNC.__stage7g_snapshot_20260928__ = sink,
    multiFSYNC.__private_scale_trace_8d72c46f__ = scale_control)
  on.exit(options(old_options), add = TRUE)
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  outcome <- tryCatch(withCallingHandlers(
    do.call(multiFSYNC::bayesSYNC_multi_pre_score,
            fit_args(x$observation)),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  new_lines(warnings, file.path(result_dir, "WARNINGS.txt"))
  new_rds(sink$records,
          file.path(result_dir, "RAW_STATE_SNAPSHOTS_SELECTED.rds"))
  if (length(sink$records))
    new_csv(data.frame(
      iteration = vapply(sink$records, `[[`, integer(1L), "iteration"),
      temperature = vapply(sink$records, `[[`, numeric(1L), "temperature"),
      annealing = vapply(sink$records, `[[`, logical(1L), "annealing")),
      file.path(result_dir, "RAW_STATE_INDEX.csv"))
  if (inherits(outcome, "error")) {
    new_lines(c("STAGE7K_S_FIT_ERROR", conditionMessage(outcome)),
              file.path(result_dir, "FIT_ERROR.txt"))
    stop("S continuation failed; partial snapshots retained.")
  }
  fit <- outcome
  new_rds(fit, file.path(result_dir, "fit.rds"))
  new_csv(fit$driver_trace, file.path(result_dir, "DRIVER_TRACE.csv"))
  if (!is.null(fit$scale_trace))
    new_csv(fit$scale_trace, file.path(result_dir, "SCALE_TRACE.csv"))
  new_csv(data.frame(t1_sweep = seq_along(fit$ELBO),
                     ordinary_elbo = as.numeric(fit$ELBO)),
          file.path(result_dir, "ORDINARY_ELBO.csv"))
  valid <- isTRUE(resume$applied) &&
    length(sink$records) == length(selected) &&
    all_finite(sink$records) &&
    identical(vapply(sink$records, `[[`, integer(1L), "iteration"),
              selected) &&
    fit$i_iter == 499L && fit$annealing_sweeps == 99L &&
    fit$t1_sweeps == 400L && length(fit$ELBO) == 400L &&
    all(is.finite(fit$ELBO)) && isTRUE(fit$elbo_objective_valid) &&
    identical(as.integer(fit$fit_seed), 99716211L)
  if (!valid) {
    new_lines("STAGE7K_S_FIT_CONTRACT_FAILED",
              file.path(result_dir, "FIT_CONTRACT_FAILED.txt"))
    stop("S continuation violated the fixed contract.")
  }
  output_paths <- c(fit = file.path(result_dir, "fit.rds"),
    snapshots = file.path(result_dir, "RAW_STATE_SNAPSHOTS_SELECTED.rds"),
    state_index = file.path(result_dir, "RAW_STATE_INDEX.csv"),
    driver_trace = file.path(result_dir, "DRIVER_TRACE.csv"),
    elbo = file.path(result_dir, "ORDINARY_ELBO.csv"),
    warnings = file.path(result_dir, "WARNINGS.txt"))
  if (file.exists(file.path(result_dir, "SCALE_TRACE.csv")))
    output_paths <- c(output_paths,
                      scale_trace = file.path(result_dir, "SCALE_TRACE.csv"))
  binding_path <- file.path(result_dir, "BRANCH_OUTPUT_BINDING.csv")
  new_csv(bind(output_paths), binding_path)
  new_lines(c("STAGE7K_S_FIT_COMPLETE", "branch=S",
              "source_sweep=50", "resumed_sweeps=449",
              "ordinary_sweeps=400", "truth_opened=FALSE",
              "formal_candidate_pool=FALSE",
              paste0("saved_snapshots=", length(selected)),
              paste0("final_ordinary_elbo=",
                     format(tail(fit$ELBO, 1L), digits = 16L)),
              paste0("binding_sha256=", sha(binding_path)),
              paste0("elapsed_seconds=", format(elapsed, digits = 10L))),
            file.path(result_dir, "FIT_COMPLETE.txt"))
  cat("STAGE7K_S_FIT_COMPLETE\n")
}
if (action == "status") {
  cat("S_complete=", file.exists(file.path(out, "S_full_v1",
                                         "FIT_COMPLETE.txt")),
      " S_error=", file.exists(file.path(out, "S_full_v1",
                                      "FIT_ERROR.txt")), "\n", sep = "")
}
