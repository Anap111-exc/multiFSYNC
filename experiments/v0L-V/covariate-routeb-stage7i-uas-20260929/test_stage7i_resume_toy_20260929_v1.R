#!/usr/bin/env Rscript

# Small, independent replay check for the private annealing resume hook.
# Usage: Rscript this.R original|resume PACKAGE_LIBRARY OUTPUT_DIRECTORY
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("original", "resume")) {
  stop("Usage: Rscript test_stage7i_resume_toy_20260929_v1.R ",
       "original|resume PACKAGE_LIBRARY OUTPUT_DIRECTORY", call. = FALSE)
}
mode <- args[[1L]]
.libPaths(c(normalizePath(args[[2L]], mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
out <- args[[3L]]
if (mode == "original" && (dir.exists(out) || file.exists(out))) {
  stop("Refusing existing toy output directory.")
}
if (mode == "resume" && !dir.exists(out)) {
  stop("Original toy output directory is missing.")
}

fit_args <- function(observation) list(
  Y = observation$Y, Z = observation$Z,
  time_obs = observation$time_obs,
  L_f = 1L, L_s = c(1L, 1L), M_f = 1L, M_s = list(1L, 1L),
  K = 3L, n_g = 21L, anneal = c(1, 1.9, 100), maxit = 499L,
  n_cpus = 1L, verbose = FALSE, seed = 99102L,
  tol_abs = 0, tol_rel = 0, bool_scale = FALSE,
  bool_var_spec_prob = FALSE, d_0 = 3L,
  convergence_rule = "parameters", lambda_orth = 0,
  initialization = "random", pre_score_sweeps = 1L,
  trace_sweeps = c(51L, 52L),
  function_initialization = "gram_unit_energy")

snapshot_sink <- function(sweeps) {
  env <- new.env(parent = emptyenv())
  env$schema <- "stage7g_snapshot_v1"
  env$sweeps <- as.integer(sweeps)
  env$records <- list()
  env
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
    if (!identical(is.na(x), is.na(y)) ||
        any(!is.finite(x[!is.na(x)])) ||
        any(!is.finite(y[!is.na(y)]))) return(Inf)
    return(max(abs(as.numeric(x[!is.na(x)]) -
                   as.numeric(y[!is.na(y)])), 0))
  }
  if (identical(x, y)) 0 else Inf
}

if (mode == "original") {
  data <- multiFSYNC::simulate_multi_study_osullivan(
    S = 2L, n_s = c(3L, 3L), p = 3L, d = 1L,
    L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L, 1L),
    K = 3L, n_obs = 7L, common_grid = TRUE,
    sigma_eps = 0.10, seed = 99101L, use_explicit = FALSE)
  observation <- list(Y = data$Y, Z = data$Z,
                      time_obs = data$time_obs)
  sink <- snapshot_sink(50:52)
  old_options <- options(multiFSYNC.__stage7g_snapshot_20260928__ = sink)
  on.exit(options(old_options), add = TRUE)
  fit <- suppressWarnings(do.call(multiFSYNC::bayesSYNC_multi_pre_score,
                                  fit_args(observation)))
  if (!identical(vapply(sink$records, `[[`, integer(1), "iteration"),
                 50:52)) stop("Original toy snapshots 50:52 are missing.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(observation, file.path(out, "TOY_OBSERVATION.rds"))
  saveRDS(sink$records, file.path(out, "TOY_ORIGINAL_SNAPSHOTS.rds"))
  saveRDS(fit$driver_trace[fit$driver_trace$iteration %in% 51:52,
                           , drop = FALSE],
          file.path(out, "TOY_ORIGINAL_DRIVER_TRACE.rds"))
  writeLines("TOY_ORIGINAL_COMPLETE", file.path(out, "ORIGINAL_COMPLETE.txt"))
  cat("TOY_ORIGINAL_COMPLETE\n")
} else {
  observation <- readRDS(file.path(out, "TOY_OBSERVATION.rds"))
  original <- readRDS(file.path(out, "TOY_ORIGINAL_SNAPSHOTS.rds"))
  if (!identical(vapply(original, `[[`, integer(1), "iteration"), 50:52)) {
    stop("Original toy snapshots are not bound to 50:52.")
  }
  resume <- new.env(parent = emptyenv())
  resume$schema <- "stage7i_sweep50_resume_v1"
  resume$state <- original[[1L]]
  resume$stop_after_iteration <- 52L
  sink <- snapshot_sink(51:52)
  old_options <- options(multiFSYNC.__stage7i_resume_20260929__ = resume,
                         multiFSYNC.__stage7g_snapshot_20260928__ = sink)
  on.exit(options(old_options), add = TRUE)
  fit <- suppressWarnings(do.call(multiFSYNC::bayesSYNC_multi_pre_score,
                                  fit_args(observation)))
  if (!isTRUE(resume$applied) || length(sink$records) != 2L ||
      !identical(vapply(sink$records, `[[`, integer(1), "iteration"),
                 51:52)) stop("Private resume hook did not capture 51:52.")
  fields <- names(original[[2L]])
  differences <- vapply(fields, function(name) max(
    max_diff(original[[2L]][[name]], sink$records[[1L]][[name]]),
    max_diff(original[[3L]][[name]], sink$records[[2L]][[name]])),
    numeric(1))
  if (any(!is.finite(differences)) || max(differences) > 1e-8) {
    print(sort(differences, decreasing = TRUE)[1:10])
    stop("Toy continuation snapshots diverged.")
  }
  original_trace <- readRDS(file.path(out,
                                      "TOY_ORIGINAL_DRIVER_TRACE.rds"))
  trace_diff <- max_diff(original_trace,
                         fit$driver_trace[fit$driver_trace$iteration %in% 51:52,
                                          , drop = FALSE])
  if (!is.finite(trace_diff) || trace_diff > 1e-8) {
    stop("Toy continuation driver trace diverged: ", trace_diff)
  }
  write.csv(data.frame(field = fields, max_abs_difference = differences),
            file.path(out, "TOY_RESUME_PARITY.csv"), row.names = FALSE)
  writeLines(c("TOY_RESUME_PASS", paste0("max_state_abs=",
              max(differences)), paste0("max_trace_abs=", trace_diff)),
             file.path(out, "TOY_RESUME_PASS.txt"))
  cat("TOY_RESUME_PASS max_state_abs=", max(differences),
      " max_trace_abs=", trace_diff, "\n", sep = "")
}
