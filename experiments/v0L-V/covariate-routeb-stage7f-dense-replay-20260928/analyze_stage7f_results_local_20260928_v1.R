#!/usr/bin/env Rscript

# Read-only local derivative of the verified Stage 7F archive.
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript analyze_stage7f_results_local_20260928_v1.R STAGE7F_RESULT_ROOT")
root <- normalizePath(args[[1L]], mustWork = TRUE)
parity_marker <- file.path(root, "endpoint_parity_v1", "PARITY_PASS.txt")
attribution_marker <- file.path(root, "attribution_v1",
                                "ATTRIBUTION_TABLES_COMPLETE.txt")
if (!file.exists(parity_marker) || !file.exists(attribution_marker)) {
  stop("Fit parity and first-stage attribution must be complete: ",
       parity_marker, "=", file.exists(parity_marker), "; ",
       attribution_marker, "=", file.exists(attribution_marker))
}
out <- file.path(dirname(root), "local_analysis_v1")
if (dir.exists(out)) stop("Refusing to overwrite local analysis.")
dir.create(out)
write_new <- function(x, name) {
  path <- file.path(out, name)
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  write.csv(x, path, row.names = FALSE, na = "")
}
trajectory <- read.csv(file.path(root, "attribution_v1",
                                 "SPECIFIC_BLOCK_TRAJECTORY.csv"),
                       stringsAsFactors = FALSE)
trajectory$seed <- as.integer(sub("^covrb7c_([0-9][0-9]).*$", "\\1",
                                  trajectory$fit_id))
stopifnot(setequal(trajectory$seed, 1:3),
          setequal(trajectory$study, 1:2))

key_sweeps <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                55, 60, 65, 70, 75, 80, 90, 99, 100, 110, 150,
                200, 300, 400, 499)
snapshot <- trajectory[trajectory$total_sweep %in% key_sweeps, , drop = FALSE]
snapshot <- snapshot[order(snapshot$seed, snapshot$study,
                           snapshot$total_sweep), , drop = FALSE]
write_new(snapshot, "KEY_SNAPSHOTS.csv")

# First recorded observation of each threshold, as a descriptive event only.
event_definitions <- list(
  loading_ppi_lt_0.1 = function(x) x$ppi_mean < 0.1,
  loading_ppi_lt_0.05 = function(x) x$ppi_mean < 0.05,
  loading_ppi_lt_0.01 = function(x) x$ppi_mean < 0.01,
  loading_ppi_lt_1e_4 = function(x) x$ppi_mean < 1e-4,
  loading_norm_lt_1 = function(x) x$loading_posterior_mean_norm < 1,
  loading_norm_lt_0.01 = function(x) x$loading_posterior_mean_norm < 0.01,
  score_norm_lt_0.01 = function(x) x$score_posterior_mean_norm < 0.01,
  time_coeff_norm_lt_0.01 = function(x) x$time_coefficient_norm_combined < 0.01,
  loading_driver_ratio_lt_0.1 = function(x)
    x$loading_driver_precision_ratio < 0.1
)
events <- list()
for (seed in 1:3) for (study in 1:2) {
  x <- trajectory[trajectory$seed == seed &
                    trajectory$study == study, , drop = FALSE]
  x <- x[order(x$total_sweep), , drop = FALSE]
  for (event in names(event_definitions)) {
    hits <- which(event_definitions[[event]](x))
    j <- if (length(hits)) hits[[1L]] else NA_integer_
    events[[length(events) + 1L]] <- data.frame(
      seed = seed, study = study, event = event,
      first_recorded_sweep = if (is.na(j)) NA_integer_ else x$total_sweep[[j]],
      previous_recorded_sweep = if (is.na(j) || j == 1L)
        NA_integer_ else x$total_sweep[[j - 1L]],
      phase = if (is.na(j)) NA_character_ else x$phase[[j]],
      temperature = if (is.na(j)) NA_real_ else x$temperature[[j]],
      stringsAsFactors = FALSE
    )
  }
}
write_new(do.call(rbind, events), "FIRST_RECORDED_EVENTS.csv")

paired <- merge(
  trajectory[trajectory$study == 1L, c("fit_id", "seed", "total_sweep",
                                      "ppi_mean", "loading_posterior_mean_norm",
                                      "score_posterior_mean_norm",
                                      "time_coefficient_norm_combined")],
  trajectory[trajectory$study == 2L, c("fit_id", "seed", "total_sweep",
                                      "ppi_mean", "loading_posterior_mean_norm",
                                      "score_posterior_mean_norm",
                                      "time_coefficient_norm_combined")],
  by = c("fit_id", "seed", "total_sweep"), suffixes = c("_B1", "_B2")
)
paired$B2_to_B1_loading_norm <-
  paired$loading_posterior_mean_norm_B2 /
  pmax(paired$loading_posterior_mean_norm_B1, .Machine$double.eps)
write_new(paired, "B2_VS_B1_TRAJECTORY.csv")

parity <- read.csv(file.path(root, "endpoint_parity_v1",
                             "PARAMETER_PARITY.csv"))
checks <- read.csv(file.path(root, "endpoint_parity_v1",
                             "CONFIGURATION_PARITY.csv"))
old_trace <- read.csv(file.path(root, "endpoint_parity_v1",
                                "ORIGINAL_TRACE_PARITY.csv"))
stopifnot(nrow(parity) == 54L, all(parity$pass),
          all(parity$max_abs_difference == 0),
          all(checks$configuration_and_sweeps_equal),
          nrow(old_trace) == 9L, all(old_trace$pass),
          all(old_trace$max_abs_metric_difference == 0))
writeLines(c("STAGE7F_LOCAL_ANALYSIS_COMPLETE",
             "parameter_parity_checks=54_all_exact",
             "original_trace_checks=9_all_exact",
             "truth_read=FALSE"),
           file.path(out, "LOCAL_ANALYSIS_COMPLETE.txt"))
cat("STAGE7F_LOCAL_ANALYSIS_COMPLETE\n")
