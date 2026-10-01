#!/usr/bin/env Rscript

# Analyze diagnostic paths only after the frozen Stage 7C endpoints are matched.
# This script reads no truth and performs no fit or selection.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript analyze_stage7f_trace_20260928_v1.R STAGE7F_ROOT")
root <- normalizePath(args[[1L]], mustWork = TRUE)
parity_root <- file.path(root, "endpoint_parity_v1")
if (!file.exists(file.path(parity_root, "PARITY_PASS.txt")) ||
    file.exists(file.path(parity_root, "PARITY_FAILED.txt"))) {
  stop("Endpoint parity must pass before trajectory attribution.")
}
out <- file.path(root, "attribution_v1")
if (dir.exists(out) || file.exists(out)) stop("Refusing attribution overwrite.")
dir.create(out)
write_new <- function(value, name) {
  path <- file.path(out, name)
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  utils::write.csv(value, path, row.names = FALSE, na = "")
}

trace <- utils::read.csv(file.path(parity_root, "TRACE_COMBINED.csv"),
                         stringsAsFactors = FALSE)
required <- c("fit_id", "data_id", "iteration", "temperature",
              "annealing", "phase", "block", "study", "component",
              "ppi_mean", "posterior_mean_norm", "data_driver_norm",
              "driver_precision_ratio")
if (!all(required %in% names(trace))) stop("Trace schema changed.")
trace <- trace[trace$phase == "main", , drop = FALSE]
trace <- trace[order(trace$fit_id, trace$iteration, trace$sequence), , drop = FALSE]
manifest <- utils::read.csv(file.path(root, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
expected_sweeps <- utils::read.csv(file.path(root, "TRACE_SCHEDULE.csv"))$iteration
if (!identical(sort(unique(trace$fit_id)), sort(manifest$fit_id))) {
  stop("Trace fit IDs do not match registered manifest.")
}
if (!all(vapply(split(trace$iteration, trace$fit_id), function(x)
  setequal(unique(x), expected_sweeps), logical(1L)))) {
  stop("Trace schedule is incomplete.")
}

rows <- list()
index <- 0L
for (fit_id in manifest$fit_id) {
  x <- trace[trace$fit_id == fit_id, , drop = FALSE]
  for (iteration in expected_sweeps) {
    current <- x[x$iteration == iteration, , drop = FALSE]
    for (study in 1:2) {
      loading <- current[current$block == "specific_loading" &
                           current$study == study, , drop = FALSE]
      score <- current[current$block == "specific_score" &
                         current$study == study, , drop = FALSE]
      function_rows <- current[current$block == "specific_time_function" &
                                 current$study == study, , drop = FALSE]
      if (nrow(loading) != 1L || nrow(score) != 1L ||
          nrow(function_rows) != 2L ||
          !setequal(function_rows$component, 1:2)) {
        stop("Missing or duplicated specific block at ", fit_id,
             " iteration ", iteration, " study ", study)
      }
      function_rows <- function_rows[order(function_rows$component), , drop = FALSE]
      index <- index + 1L
      rows[[index]] <- data.frame(
        fit_id = fit_id, data_id = loading$data_id,
        total_sweep = iteration,
        phase = if (isTRUE(loading$annealing)) "annealing" else "ordinary_T1",
        temperature = loading$temperature, study = study,
        ppi_mean = loading$ppi_mean,
        ppi_min = loading$ppi_min, ppi_max = loading$ppi_max,
        loading_posterior_mean_norm = loading$posterior_mean_norm,
        loading_slab_mean_norm = loading$slab_mean_norm,
        loading_data_driver_norm = loading$data_driver_norm,
        loading_driver_precision_ratio = loading$driver_precision_ratio,
        loading_log_prior_odds_mean = loading$log_prior_odds_mean,
        score_posterior_mean_norm = score$posterior_mean_norm,
        score_data_driver_norm = score$data_driver_norm,
        score_driver_precision_ratio = score$driver_precision_ratio,
        time_coefficient_norm_m1 = function_rows$posterior_mean_norm[[1L]],
        time_coefficient_norm_m2 = function_rows$posterior_mean_norm[[2L]],
        time_coefficient_norm_combined = sqrt(sum(
          function_rows$posterior_mean_norm^2)),
        time_driver_precision_ratio_m1 =
          function_rows$driver_precision_ratio[[1L]],
        time_driver_precision_ratio_m2 =
          function_rows$driver_precision_ratio[[2L]],
        stringsAsFactors = FALSE
      )
    }
  }
}
trajectory <- do.call(rbind, rows)
write_new(trajectory, "SPECIFIC_BLOCK_TRAJECTORY.csv")

thresholds <- c(1e-2, 1e-4, 1e-6)
crossings <- list()
changes <- list()
first_close <- list()
for (key in unique(paste(trajectory$fit_id, trajectory$study, sep = "__"))) {
  x <- trajectory[paste(trajectory$fit_id, trajectory$study,
                        sep = "__") == key, , drop = FALSE]
  x <- x[order(x$total_sweep), , drop = FALSE]
  if (anyDuplicated(x$total_sweep)) stop("Duplicate total sweep in trajectory.")
  for (threshold in thresholds) {
    hits <- which(x$ppi_mean < threshold)
    first <- if (length(hits)) hits[[1L]] else NA_integer_
    sustained <- which(vapply(seq_len(nrow(x)), function(i)
      all(x$ppi_mean[i:nrow(x)] < threshold), logical(1L)))
    stable <- if (length(sustained)) sustained[[1L]] else NA_integer_
    crossings[[length(crossings) + 1L]] <- data.frame(
      fit_id = x$fit_id[[1L]], study = x$study[[1L]],
      ppi_threshold = threshold,
      first_below = if (is.na(first)) NA_integer_ else x$total_sweep[[first]],
      first_below_previous_record = if (is.na(first) || first == 1L)
        NA_integer_ else x$total_sweep[[first - 1L]],
      first_stably_below = if (is.na(stable)) NA_integer_ else
        x$total_sweep[[stable]],
      first_stably_below_previous_record = if (is.na(stable) || stable == 1L)
        NA_integer_ else x$total_sweep[[stable - 1L]],
      terminal_ppi = tail(x$ppi_mean, 1L),
      stringsAsFactors = FALSE
    )
  }
  closed <- x$ppi_mean < 1e-4 &
    x$loading_posterior_mean_norm < 1e-3
  stable_close <- which(vapply(seq_len(nrow(x)), function(i)
    all(closed[i:nrow(x)]), logical(1L)))
  j <- if (length(stable_close)) stable_close[[1L]] else NA_integer_
  first_close[[length(first_close) + 1L]] <- data.frame(
    fit_id = x$fit_id[[1L]], study = x$study[[1L]],
    first_stable_joint_closure = if (is.na(j)) NA_integer_ else
      x$total_sweep[[j]],
    previous_record = if (is.na(j) || j == 1L) NA_integer_ else
      x$total_sweep[[j - 1L]],
    closure_during_annealing = if (is.na(j)) NA else
      identical(x$phase[[j]], "annealing"),
    terminal_loading_norm = tail(x$loading_posterior_mean_norm, 1L),
    terminal_ppi = tail(x$ppi_mean, 1L),
    stringsAsFactors = FALSE
  )
  if (nrow(x) > 1L) {
    for (j in 2:nrow(x)) {
      changes[[length(changes) + 1L]] <- data.frame(
        fit_id = x$fit_id[[1L]], study = x$study[[1L]],
        from_sweep = x$total_sweep[[j - 1L]],
        to_sweep = x$total_sweep[[j]],
        phase_at_end = x$phase[[j]],
        temperature_at_end = x$temperature[[j]],
        log10_ppi_change = log10(max(x$ppi_mean[[j]], 1e-16)) -
          log10(max(x$ppi_mean[[j - 1L]], 1e-16)),
        loading_norm_change = x$loading_posterior_mean_norm[[j]] -
          x$loading_posterior_mean_norm[[j - 1L]],
        score_norm_change = x$score_posterior_mean_norm[[j]] -
          x$score_posterior_mean_norm[[j - 1L]],
        time_coefficient_norm_change =
          x$time_coefficient_norm_combined[[j]] -
          x$time_coefficient_norm_combined[[j - 1L]],
        stringsAsFactors = FALSE
      )
    }
  }
}
write_new(do.call(rbind, crossings), "PPI_THRESHOLD_CROSSINGS.csv")
write_new(do.call(rbind, first_close), "FIRST_STABLE_JOINT_CLOSURE.csv")
write_new(do.call(rbind, changes), "ADJACENT_RECORDED_CHANGES.csv")
writeLines(c("STAGE7F_TRUTH_BLIND_TRACE_ANALYSIS_COMPLETE",
             "requires_endpoint_parity_pass=TRUE",
             "truth_read=FALSE",
             "joint_closure_definition=ppi_mean<1e-4_and_loading_norm<1e-3_persisting_to_terminal",
             "time_norm_definition=coefficient_Euclidean_norm_not_function_L2"),
           file.path(out, "ATTRIBUTION_TABLES_COMPLETE.txt"))
cat("STAGE7F_ATTRIBUTION_TABLES_COMPLETE fits=3 sweeps=",
    length(expected_sweeps), "\n", sep = "")
