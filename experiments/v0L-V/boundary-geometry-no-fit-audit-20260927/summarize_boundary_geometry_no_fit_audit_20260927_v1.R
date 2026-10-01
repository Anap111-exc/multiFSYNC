#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)

repo_candidates <- unique(c(
  Sys.getenv("MULTIFSYNC_REPO_ROOT", unset = ""),
  getwd(), file.path(getwd(), "Rcode", "multiFSYNC")))
repo_candidates <- repo_candidates[
  nzchar(repo_candidates) &
    file.exists(file.path(repo_candidates, "DESCRIPTION"))]
if (!length(repo_candidates)) stop("Cannot locate package root.")
repo_root <- normalizePath(repo_candidates[[1L]], mustWork = TRUE)
audit_root <- file.path(
  repo_root, "experiments", "v0L-V",
  "boundary-geometry-no-fit-audit-20260927",
  "boundary_geometry_no_fit_audit_20260927_v1")
metrics_path <- file.path(audit_root, "FUNCTION_LEVEL_METRICS.csv")
if (!file.exists(metrics_path)) stop("Function-level metrics are missing.")

output_paths <- file.path(audit_root, c(
  "ENDPOINT_DERIVATIVE_SUMMARY.csv",
  "FLAG_BREAKDOWN.csv",
  "ASSOCIATION_SUMMARY.csv"))
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite an existing derived summary.")
}

metrics <- utils::read.csv(metrics_path, check.names = FALSE)
q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE, type = 7L))

endpoint_rows <- list()
position <- 0L
derivatives <- list(
  function_value = c("f_0", "f_1"),
  first_derivative = c("d1_0", "d1_1"),
  second_derivative = c("d2_0", "d2_1"))
for (family in unique(metrics$family)) {
  x <- metrics[metrics$family == family, , drop = FALSE]
  for (derivative in names(derivatives)) {
    columns <- derivatives[[derivative]]
    left <- abs(x[[columns[[1L]]]])
    right <- abs(x[[columns[[2L]]]])
    difference <- abs(x[[columns[[1L]]]] - x[[columns[[2L]]]])
    position <- position + 1L
    endpoint_rows[[position]] <- data.frame(
      family = family,
      derivative = derivative,
      n_functions = nrow(x),
      left_abs_median = stats::median(left),
      left_abs_q95 = q(left, 0.95),
      left_abs_max = max(left),
      right_abs_median = stats::median(right),
      right_abs_q95 = q(right, 0.95),
      right_abs_max = max(right),
      endpoint_abs_difference_median = stats::median(difference),
      endpoint_abs_difference_q95 = q(difference, 0.95),
      endpoint_abs_difference_max = max(difference),
      stringsAsFactors = FALSE)
  }
}
endpoint_summary <- do.call(rbind, endpoint_rows)

flags <- c(
  "constraint_fail", "edge_energy_high_flag", "edge_energy_low_flag",
  "edge_peak_flag", "projection_l2_high_flag")
flag_rows <- list()
position <- 0L
for (group_name in c("family", "family_role")) {
  grouping <- if (group_name == "family") {
    metrics$family
  } else {
    paste(metrics$family, metrics$role, sep = "::")
  }
  for (group in unique(grouping)) {
    x <- metrics[grouping == group, , drop = FALSE]
    for (flag in flags) {
      position <- position + 1L
      flag_rows[[position]] <- data.frame(
        grouping = group_name,
        group = group,
        flag = flag,
        n_functions = nrow(x),
        flagged_n = sum(x[[flag]]),
        flagged_rate = mean(x[[flag]]),
        stringsAsFactors = FALSE)
    }
  }
}
flag_breakdown <- do.call(rbind, flag_rows)

association_pairs <- list(
  roughness_vs_projection = c(
    "normalized_roughness", "projection_relative_l2"),
  edge_energy_vs_edge_peak = c(
    "edge_energy_fraction", "edge_peak_ratio"),
  endpoint_magnitude_vs_edge_energy = c(
    "f_0", "edge_energy_fraction"))
association_rows <- list()
position <- 0L
for (family in unique(metrics$family)) {
  x <- metrics[metrics$family == family, , drop = FALSE]
  x$f_0 <- abs(x$f_0)
  for (association in names(association_pairs)) {
    columns <- association_pairs[[association]]
    position <- position + 1L
    association_rows[[position]] <- data.frame(
      family = family,
      association = association,
      variable_1 = columns[[1L]],
      variable_2 = columns[[2L]],
      spearman_rho = suppressWarnings(stats::cor(
        x[[columns[[1L]]]], x[[columns[[2L]]]],
        method = "spearman", use = "complete.obs")),
      stringsAsFactors = FALSE)
  }
}
association_summary <- do.call(rbind, association_rows)

utils::write.csv(endpoint_summary, output_paths[[1L]], row.names = FALSE)
utils::write.csv(flag_breakdown, output_paths[[2L]], row.names = FALSE)
utils::write.csv(association_summary, output_paths[[3L]], row.names = FALSE)
cat("BOUNDARY_GEOMETRY_DERIVED_SUMMARY_PASS\n")
