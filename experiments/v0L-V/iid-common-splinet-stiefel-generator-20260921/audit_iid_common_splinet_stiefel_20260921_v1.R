#!/usr/bin/env Rscript

repo_root_env <- Sys.getenv("MULTIFSYNC_REPO_ROOT", unset = "")
repo_candidates <- c(
  repo_root_env,
  getwd(),
  file.path(getwd(), "Rcode", "multiFSYNC"))
repo_candidates <- unique(repo_candidates[nzchar(repo_candidates)])
repo_hits <- repo_candidates[file.exists(file.path(repo_candidates,
                                                    "DESCRIPTION"))]
if (length(repo_hits) > 0L) {
  repo_root <- normalizePath(repo_hits[[1L]], mustWork = TRUE)
} else {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Cannot locate the multiFSYNC package root.")
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  repo_root <- normalizePath(
    file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
}
script_dir <- file.path(
  repo_root, "experiments", "v0L-V",
  "iid-common-splinet-stiefel-generator-20260921")

library(splines)
source(file.path(repo_root, "R", "utils_multi.R"), local = .GlobalEnv)
source(file.path(repo_root, "R", "generate_data_structured.R"),
       local = .GlobalEnv)

if (!requireNamespace("Splinets", quietly = TRUE)) {
  stop("Splinets is required for this audit.")
}

output_dir <- file.path(
  script_dir, "acceptance", "windows", "no_fit_geometry_audit_v1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_draw <- 500L
M_by_block <- c(2L, 2L, 2L)
basis_dim <- 8L
n_grid <- 201L
audit_seed <- 20260921L
set.seed(audit_seed)

t_grid <- seq(0, 1, length.out = n_grid)
dt <- diff(t_grid)
w <- c(dt[1L] / 2,
       (dt[-1L] + dt[-length(dt)]) / 2,
       dt[length(dt)] / 2)
K_fit <- 5L
int_knots <- seq(0, 1, length.out = K_fit)[-c(1L, K_fit)]
C_fit <- cbind(
  1, t_grid,
  ZOSull(t_grid, range.x = c(0, 1), intKnots = int_knots))
C_weighted <- sqrt(w) * C_fit

relative_projection_floor <- function(phi) {
  coef <- qr.solve(C_weighted, sqrt(w) * phi)
  residual <- phi - C_fit %*% coef
  sqrt(colSums(w * residual^2) / colSums(w * phi^2))
}

roughness_fd <- function(phi) {
  h <- t_grid[[2L]] - t_grid[[1L]]
  second <- apply(phi, 2L, diff, differences = 2L) / h^2
  colSums(second^2) * h
}

boundary_energy <- function(phi) {
  boundary <- t_grid <= 0.1 | t_grid >= 0.9
  colSums(w[boundary] * phi[boundary, , drop = FALSE]^2) /
    colSums(w * phi^2)
}

block_rows <- vector("list", n_draw * length(M_by_block))
pair_rows <- vector("list", n_draw * choose(length(M_by_block), 2L))
block_pos <- 0L
pair_pos <- 0L

for (draw_id in seq_len(n_draw)) {
  blocks <- construct_iid_common_splinet_stiefel_blocks(
    M_by_block, n_g = n_grid, basis_dim = basis_dim)

  for (q in seq_along(blocks)) {
    block_pos <- block_pos + 1L
    block <- blocks[[q]]
    rough <- roughness_fd(block$phi)
    edge <- boundary_energy(block$phi)
    projection <- relative_projection_floor(block$phi)
    block_rows[[block_pos]] <- data.frame(
      draw_id = draw_id,
      block_id = q,
      gram_error_continuous = max(abs(block$gram - diag(M_by_block[[q]]))),
      gram_error_grid = max(abs(block$grid_gram - diag(M_by_block[[q]]))),
      roughness_mean = mean(rough),
      roughness_max = max(rough),
      boundary_energy_mean = mean(edge),
      boundary_energy_max = max(edge),
      projection_floor_mean = mean(projection),
      projection_floor_max = max(projection),
      stringsAsFactors = FALSE
    )
  }

  for (q1 in seq_len(length(blocks) - 1L)) {
    for (q2 in (q1 + 1L):length(blocks)) {
      pair_pos <- pair_pos + 1L
      cross_gram <- crossprod(
        blocks[[q1]]$stiefel_coefficients,
        blocks[[q2]]$stiefel_coefficients)
      canonical <- svd(cross_gram, nu = 0, nv = 0)$d
      pair_rows[[pair_pos]] <- data.frame(
        draw_id = draw_id,
        block_1 = q1,
        block_2 = q2,
        max_canonical_correlation = max(canonical),
        min_canonical_correlation = min(canonical),
        stringsAsFactors = FALSE
      )
    }
  }
}

block_metrics <- do.call(rbind, block_rows)
pair_metrics <- do.call(rbind, pair_rows)

metric_names <- c(
  "roughness_mean", "boundary_energy_mean", "projection_floor_mean")
exchangeability <- do.call(rbind, lapply(metric_names, function(metric) {
  block_means <- tapply(block_metrics[[metric]], block_metrics$block_id, mean)
  pooled_sd <- stats::sd(block_metrics[[metric]])
  standardized_range <- if (is.finite(pooled_sd) && pooled_sd > 0) {
    diff(range(block_means)) / pooled_sd
  } else {
    0
  }
  data.frame(
    metric = metric,
    block_1_mean = unname(block_means[[1L]]),
    block_2_mean = unname(block_means[[2L]]),
    block_3_mean = unname(block_means[[3L]]),
    pooled_sd = pooled_sd,
    standardized_mean_range = standardized_range,
    passed = standardized_range <= 0.25,
    stringsAsFactors = FALSE
  )
}))

projection_all <- block_metrics$projection_floor_mean
validity <- data.frame(
  check_id = c(
    "continuous_gram_exact",
    "all_metrics_finite",
    "no_numerically_identical_subspaces",
    "cross_block_not_forced_orthogonal",
    "role_exchangeability",
    "projection_floor_usable"),
  passed = c(
    max(block_metrics$gram_error_continuous) <= 1e-8,
    all(vapply(block_metrics, function(x) all(is.finite(x)), logical(1))) &&
      all(vapply(pair_metrics, function(x) all(is.finite(x)), logical(1))),
    max(pair_metrics$max_canonical_correlation) < 1 - 1e-8,
    stats::median(pair_metrics$max_canonical_correlation) > 0.1,
    all(exchangeability$passed),
    stats::median(projection_all) <= 0.25 &&
      unname(stats::quantile(projection_all, 0.95)) <= 0.60),
  detail = c(
    sprintf("max=%.3e", max(block_metrics$gram_error_continuous)),
    sprintf("block_rows=%d;pair_rows=%d", nrow(block_metrics),
            nrow(pair_metrics)),
    sprintf("max=%.6f", max(pair_metrics$max_canonical_correlation)),
    sprintf("median=%.6f",
            stats::median(pair_metrics$max_canonical_correlation)),
    paste(sprintf("%s=%.3f", exchangeability$metric,
                  exchangeability$standardized_mean_range), collapse = ";"),
    sprintf("median=%.4f;q95=%.4f",
            stats::median(projection_all),
            unname(stats::quantile(projection_all, 0.95)))),
  stringsAsFactors = FALSE
)

overall_summary <- data.frame(
  quantity = c(
    "n_draw", "n_block_rows", "n_pair_rows", "basis_dim",
    "continuous_gram_error_max", "grid_gram_error_max",
    "cross_canonical_median", "cross_canonical_q95",
    "cross_canonical_gt_0_95_rate", "cross_canonical_gt_0_99_rate",
    "roughness_median", "boundary_energy_median",
    "projection_floor_median", "projection_floor_q95"),
  value = c(
    n_draw, nrow(block_metrics), nrow(pair_metrics), basis_dim,
    max(block_metrics$gram_error_continuous),
    max(block_metrics$gram_error_grid),
    stats::median(pair_metrics$max_canonical_correlation),
    unname(stats::quantile(pair_metrics$max_canonical_correlation, 0.95)),
    mean(pair_metrics$max_canonical_correlation > 0.95),
    mean(pair_metrics$max_canonical_correlation > 0.99),
    stats::median(block_metrics$roughness_mean),
    stats::median(block_metrics$boundary_energy_mean),
    stats::median(block_metrics$projection_floor_mean),
    unname(stats::quantile(block_metrics$projection_floor_mean, 0.95))),
  stringsAsFactors = FALSE
)

utils::write.csv(block_metrics, file.path(output_dir, "BLOCK_METRICS.csv"),
                 row.names = FALSE)
utils::write.csv(pair_metrics, file.path(output_dir, "PAIR_METRICS.csv"),
                 row.names = FALSE)
utils::write.csv(exchangeability,
                 file.path(output_dir, "ROLE_EXCHANGEABILITY.csv"),
                 row.names = FALSE)
utils::write.csv(overall_summary,
                 file.path(output_dir, "GEOMETRY_SUMMARY.csv"),
                 row.names = FALSE)
utils::write.csv(validity, file.path(output_dir, "AUDIT_VALIDITY.csv"),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "SESSION_INFO.txt"))

if (!all(validity$passed)) {
  writeLines("FAIL", file.path(output_dir, "AUDIT_STATUS.txt"))
  stop("No-fit geometry audit failed one or more preregistered checks.")
}

writeLines(c(
  "PASS",
  sprintf("draws=%d", n_draw),
  "formal_data_generated=FALSE",
  "model_fits_started=FALSE"),
  file.path(output_dir, "AUDIT_STATUS.txt"))
cat("NO_FIT_GEOMETRY_AUDIT_PASS\n")
