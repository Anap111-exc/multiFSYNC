#!/usr/bin/env Rscript
# Tiny toy validation of the sparse A snapshot/driver schedule. No seed02 data.
# Usage: Rscript SCRIPT PACKAGE_LIBRARY EXISTING_TOY_OUTPUT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
.libPaths(c(normalizePath(args[[1L]], mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
root <- normalizePath(args[[2L]], mustWork = TRUE)
observation <- readRDS(file.path(root, "TOY_OBSERVATION.rds"))
original <- readRDS(file.path(root, "TOY_ORIGINAL_SNAPSHOTS.rds"))
original_driver <- readRDS(file.path(root, "TOY_ORIGINAL_DRIVER_TRACE.rds"))
stopifnot(identical(vapply(original, `[[`, integer(1L), "iteration"),
                    as.integer(50:499)))
selected <- as.integer(c(51:70, 75, 80, 90, 100, 120, 150,
                         200, 250, 300, 350, 400, 450, 499))
resume <- new.env(parent = emptyenv())
resume$schema <- "stage7i_sweep50_resume_v1"
resume$state <- original[[1L]]
resume$stop_after_iteration <- 499L
sink <- new.env(parent = emptyenv())
sink$schema <- "stage7g_snapshot_v1"
sink$sweeps <- selected
sink$records <- list()
scale_trace <- list(include_initial = FALSE, include_annealing = TRUE,
                    t1_sweeps = as.integer(c(1:20, 50, 100, 200, 300,
                                             380, 400)),
                    block_annealing_sweeps = as.integer(51:60))
old_options <- options(
  multiFSYNC.__stage7i_resume_20260929__ = resume,
  multiFSYNC.__stage7g_snapshot_20260928__ = sink,
  multiFSYNC.__private_scale_trace_8d72c46f__ = scale_trace)
on.exit(options(old_options), add = TRUE)
fit <- suppressWarnings(multiFSYNC::bayesSYNC_multi_pre_score(
  Y = observation$Y, Z = observation$Z,
  time_obs = observation$time_obs,
  L_f = 1L, L_s = c(1L, 1L), M_f = 1L, M_s = list(1L, 1L),
  K = 3L, n_g = 21L, anneal = c(1, 1.9, 100), maxit = 499L,
  n_cpus = 1L, verbose = FALSE, seed = 99102L,
  tol_abs = 0, tol_rel = 0, bool_scale = FALSE,
  bool_var_spec_prob = FALSE, d_0 = 3L,
  convergence_rule = "parameters", lambda_orth = 0,
  initialization = "random", pre_score_sweeps = 1L,
  trace_sweeps = selected,
  function_initialization = "gram_unit_energy"))
stopifnot(isTRUE(resume$applied), fit$i_iter == 499L,
          identical(vapply(sink$records, `[[`, integer(1L), "iteration"),
                    selected))
reference <- original[match(selected, 50:499)]
for (i in seq_along(selected)) {
  stopifnot(isTRUE(all.equal(sink$records[[i]], reference[[i]],
                             tolerance = 1e-10,
                             check.attributes = TRUE)))
}
reference_driver <- original_driver[
  original_driver$iteration %in% selected, , drop = FALSE]
actual_driver <- fit$driver_trace[
  fit$driver_trace$iteration %in% selected, , drop = FALSE]
reference_driver$sequence <- NULL
actual_driver$sequence <- NULL
rownames(reference_driver) <- NULL
rownames(actual_driver) <- NULL
attr(reference_driver, "diagnostic_control") <- NULL
attr(actual_driver, "diagnostic_control") <- NULL
driver_comparison <- all.equal(actual_driver, reference_driver,
                               tolerance = 1e-10,
                               check.attributes = TRUE)
if (!isTRUE(driver_comparison)) {
  cat("actual_rows=", nrow(actual_driver),
      " reference_rows=", nrow(reference_driver), "\n", sep = "")
  print(driver_comparison)
  stop("Sparse toy driver trace differs from original.")
}
cat("STAGE7J_SPARSE_RESUME_TOY_PASS snapshots=", length(selected),
    " sweeps=449\n", sep = "")
