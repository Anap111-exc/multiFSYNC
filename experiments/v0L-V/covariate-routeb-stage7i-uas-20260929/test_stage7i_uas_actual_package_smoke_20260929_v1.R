#!/usr/bin/env Rscript

# Focused Stage 7I integration smoke. Run `original` with the installed frozen
# Stage 7G package, then `resume` with the isolated Stage 7I package. The only
# saved input from the simulator is Y/Z/time_obs; its truth is never passed to
# the U/A/S initializer. No long Stage 7I continuation is performed.
# Usage: Rscript this.R original|resume PACKAGE_LIBRARY OUTPUT_DIRECTORY
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
if (.Platform$OS.type == "windows") {
  locale_now <- Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
  if (is.na(locale_now)) stop("A Windows UTF-8 locale is required.")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("original", "resume")) {
  stop("Usage: Rscript test_stage7i_uas_actual_package_smoke_20260929_v1.R ",
       "original|resume PACKAGE_LIBRARY OUTPUT_DIRECTORY", call. = FALSE)
}
mode <- args[[1L]]
.libPaths(c(normalizePath(args[[2L]], mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
loaded_package <- normalizePath(system.file(package = "multiFSYNC"),
                                mustWork = TRUE)
expected_package <- normalizePath(file.path(args[[2L]], "multiFSYNC"),
                                  mustWork = TRUE)
if (!identical(loaded_package, expected_package)) {
  stop("Loaded package differs from the requested isolated library.")
}
out <- args[[3L]]
script_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_option) != 1L) stop("Cannot identify smoke script path.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_option)))

fit_args <- function(observation, maxit) list(
  Y = observation$Y, Z = observation$Z,
  time_obs = observation$time_obs,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
  K = 9L, n_g = 51L, anneal = c(1, 1.9, 100), maxit = maxit,
  n_cpus = 1L, verbose = FALSE, seed = 99152L,
  tol_abs = 0, tol_rel = 0, bool_scale = FALSE,
  bool_var_spec_prob = FALSE, d_0 = 12L,
  convergence_rule = "parameters", lambda_orth = 0,
  initialization = "random", pre_score_sweeps = 1L,
  trace_sweeps = c(51L, 52L),
  function_initialization = "gram_unit_energy"
)
snapshot_sink <- function(sweeps) {
  env <- new.env(parent = emptyenv())
  env$schema <- "stage7g_snapshot_v1"
  env$sweeps <- as.integer(sweeps)
  env$records <- list()
  env
}
run_with_capture <- function(observation, maxit, sink, resume = NULL) {
  old_options <- options(
    multiFSYNC.__stage7g_snapshot_20260928__ = sink,
    multiFSYNC.__stage7i_resume_20260929__ = resume)
  on.exit(options(old_options), add = TRUE)
  suppressWarnings(do.call(multiFSYNC::bayesSYNC_multi_pre_score,
                           fit_args(observation, maxit)))
}
all_finite <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.list(x)) return(all(vapply(x, all_finite, logical(1L))))
  if (is.numeric(x)) return(all(is.finite(x)))
  is.logical(x) && all(!is.na(x))
}
max_diff <- function(x, y) {
  if (is.null(x) || is.null(y)) return(if (is.null(x) && is.null(y)) 0 else Inf)
  if (is.list(x) || is.list(y)) {
    if (!is.list(x) || !is.list(y) || length(x) != length(y) ||
        !identical(names(x), names(y))) return(Inf)
    if (!length(x)) return(0)
    return(max(unlist(Map(max_diff, x, y), use.names = FALSE)))
  }
  if (!identical(dim(x), dim(y)) || length(x) != length(y)) return(Inf)
  if (is.numeric(x) && is.numeric(y)) {
    if (!identical(is.na(x), is.na(y))) return(Inf)
    keep <- !is.na(x)
    if (any(!is.finite(x[keep])) || any(!is.finite(y[keep]))) return(Inf)
    return(max(abs(as.numeric(x[keep]) - as.numeric(y[keep])), 0))
  }
  if (identical(x, y)) 0 else Inf
}

if (mode == "original") {
  if (file.exists(out) || dir.exists(out)) {
    stop("Refusing existing smoke output directory.")
  }
  data <- multiFSYNC::simulate_multi_study_osullivan(
    S = 2L, n_s = c(8L, 8L), p = 12L, d = 1L,
    L_f = 1L, L_s = 1L, M_f = 2L, M_s = list(2L, 2L),
    K = 9L, n_obs = 8L, common_grid = FALSE,
    sigma_eps = 0.10, seed = 99151L, use_explicit = FALSE)
  observation <- list(Y = data$Y, Z = data$Z,
                      time_obs = data$time_obs)
  rm(data)
  sink <- snapshot_sink(50:52)
  fit <- run_with_capture(observation, 103L, sink)
  if (!identical(vapply(sink$records, `[[`, integer(1L), "iteration"),
                 50:52)) stop("Frozen-package toy snapshots 50:52 are missing.")
  if (!dir.create(out, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create smoke output directory.")
  }
  saveRDS(observation, file.path(out, "OBSERVATION.rds"))
  saveRDS(fit, file.path(out, "STAGE7G_TOY_FIT.rds"))
  saveRDS(sink$records, file.path(out, "STAGE7G_TOY_SNAPSHOTS.rds"))
  writeLines(c("ORIGINAL_COMPLETE", paste0("package=", loaded_package)),
             file.path(out, "ORIGINAL_COMPLETE.txt"))
  cat("ORIGINAL_COMPLETE\n")
} else {
  if (!file.exists(file.path(out, "ORIGINAL_COMPLETE.txt"))) {
    stop("Original toy fit is incomplete.")
  }
  observation <- readRDS(file.path(out, "OBSERVATION.rds"))
  frozen_fit <- readRDS(file.path(out, "STAGE7G_TOY_FIT.rds"))
  original <- readRDS(file.path(out, "STAGE7G_TOY_SNAPSHOTS.rds"))
  if (!identical(vapply(original, `[[`, integer(1L), "iteration"), 50:52)) {
    stop("Original toy snapshots are not bound to sweeps 50:52.")
  }
  grid <- multiFSYNC:::get_grid_objects(
    unlist(observation$time_obs, recursive = FALSE), frozen_fit$K,
    n_g = frozen_fit$n_g, time_g = frozen_fit$time_g,
    format_univ = TRUE)
  n1 <- frozen_fit$n_s[[1L]]
  n2 <- frozen_fit$n_s[[2L]]
  C2 <- grid$C[(n1 + 1L):(n1 + n2)]
  source(file.path(script_dir, "build_stage7i_uas_states_20260929_v1.R"))
  starts <- stage7i_build_uas_states(observation, frozen_fit,
                                     original[[1L]], C2)
  results <- vector("list", 3L)
  names(results) <- c("U", "A", "S")
  for (branch in names(results)) {
    resume <- new.env(parent = emptyenv())
    resume$schema <- "stage7i_sweep50_resume_v1"
    resume$state <- starts[[branch]]
    resume$stop_after_iteration <- 52L
    sink <- snapshot_sink(51:52)
    fit <- run_with_capture(observation, 499L, sink, resume)
    if (!isTRUE(resume$applied) || length(sink$records) != 2L ||
        !identical(vapply(sink$records, `[[`, integer(1L), "iteration"),
                   51:52) || !all_finite(sink$records) ||
        !all_finite(fit$mu_q_normal_b) ||
        !all_finite(fit$mu_q_nu_psi)) {
      stop("Stage 7I branch failed two-sweep state smoke: ", branch)
    }
    state_difference <- if (branch == "U") {
      max(max_diff(original[[2L]], sink$records[[1L]]),
          max_diff(original[[3L]], sink$records[[2L]]))
    } else NA_real_
    if (branch == "U" && (!is.finite(state_difference) ||
                           state_difference > 1e-8)) {
      stop("U two-sweep parity failed: ", state_difference)
    }
    results[[branch]] <- data.frame(
      branch = branch, sweeps = length(sink$records),
      resumed = isTRUE(resume$applied),
      finite = all_finite(sink$records),
      U_reference_max_abs = state_difference,
      B2_effective_loading_norm_52 =
        sqrt(sum(sink$records[[2L]]$mu_q_b_specific[[2L]][, 1L]^2)),
      B2_ppi_mean_52 =
        mean(sink$records[[2L]]$mu_q_gamma_b[[2L]][, 1L]),
      stringsAsFactors = FALSE)
  }
  summary <- do.call(rbind, results)
  if (file.exists(file.path(out, "UAS_SMOKE_SUMMARY.csv")) ||
      file.exists(file.path(out, "UAS_SMOKE_PASS.txt"))) {
    stop("Refusing overwrite of smoke results.")
  }
  utils::write.csv(summary, file.path(out, "UAS_SMOKE_SUMMARY.csv"),
                   row.names = FALSE)
  saveRDS(starts$blind_diagnostics,
          file.path(out, "UAS_BLIND_INITIALIZATION_DIAGNOSTICS.rds"))
  writeLines(c("UAS_SMOKE_PASS", paste0("package=", loaded_package),
               paste0("U_reference_max_abs=",
                      summary$U_reference_max_abs[[1L]])),
             file.path(out, "UAS_SMOKE_PASS.txt"))
  cat("UAS_SMOKE_PASS U_reference_max_abs=",
      summary$U_reference_max_abs[[1L]], "\n", sep = "")
}
