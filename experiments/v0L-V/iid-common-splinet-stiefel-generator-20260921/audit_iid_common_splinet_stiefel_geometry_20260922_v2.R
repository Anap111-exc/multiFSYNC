#!/usr/bin/env Rscript

repo_root_env <- Sys.getenv("MULTIFSYNC_REPO_ROOT", unset = "")
repo_candidates <- unique(c(
  repo_root_env,
  getwd(),
  file.path(getwd(), "Rcode", "multiFSYNC")))
repo_candidates <- repo_candidates[nzchar(repo_candidates)]
repo_hits <- repo_candidates[file.exists(file.path(repo_candidates,
                                                    "DESCRIPTION"))]
if (length(repo_hits) < 1L) {
  stop("Cannot locate the multiFSYNC package root.", call. = FALSE)
}
repo_root <- normalizePath(repo_hits[[1L]], mustWork = TRUE)
script_dir <- file.path(
  repo_root, "experiments", "v0L-V",
  "iid-common-splinet-stiefel-generator-20260921")
output_dir <- file.path(
  script_dir, "acceptance", "windows", "geometry_calibration_v2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Splinets", quietly = TRUE)) {
  stop("Splinets is required for this audit.", call. = FALSE)
}
suppressPackageStartupMessages(library(splines))
source(file.path(repo_root, "R", "utils_multi.R"), local = .GlobalEnv)
source(file.path(repo_root, "R", "generate_data_structured.R"),
       local = .GlobalEnv)

audit_seed <- 20260922L
n_draw <- 500L
n_block <- 3L
M <- 2L
J <- 8L
n_grid <- 201L
K_fit <- 5L

candidates <- data.frame(
  candidate_id = c(
    "NP_U", "NP_S10", "NP_S30", "NP_S100",
    "P_U", "P_S1", "P_S3", "P_S10"),
  periodic = c(FALSE, FALSE, FALSE, FALSE,
               TRUE, TRUE, TRUE, TRUE),
  alpha = c(0, 10, 30, 100, 0, 1, 3, 10),
  stringsAsFactors = FALSE
)

t_grid <- seq(0, 1, length.out = n_grid)
dt <- diff(t_grid)
w_grid <- c(
  dt[[1L]] / 2,
  (dt[-1L] + dt[-length(dt)]) / 2,
  dt[[length(dt)]] / 2)
int_knots <- seq(0, 1, length.out = K_fit)[-c(1L, K_fit)]
C_fit <- cbind(
  1,
  t_grid,
  ZOSull(t_grid, range.x = c(0, 1), intKnots = int_knots))
C_weighted <- sqrt(w_grid) * C_fit

relative_projection_floor <- function(phi) {
  coefficients <- qr.solve(C_weighted, sqrt(w_grid) * phi)
  residual <- phi - C_fit %*% coefficients
  sqrt(colSums(w_grid * residual^2) / colSums(w_grid * phi^2))
}

boundary_energy <- function(phi) {
  boundary <- t_grid <= 0.1 | t_grid >= 0.9
  colSums(w_grid[boundary] * phi[boundary, , drop = FALSE]^2) /
    colSums(w_grid * phi^2)
}

standardize_qr_sign <- function(qr_object, M) {
  Q <- qr.Q(qr_object)[, seq_len(M), drop = FALSE]
  qr_sign <- sign(diag(qr.R(qr_object))[seq_len(M)])
  qr_sign[qr_sign == 0] <- 1
  sweep(Q, 2L, qr_sign, FUN = "*")
}

set.seed(audit_seed)
block_rows <- list()
pair_rows <- list()
block_pos <- 0L
pair_pos <- 0L

for (candidate_row in seq_len(nrow(candidates))) {
  candidate <- candidates[candidate_row, , drop = FALSE]
  knot_count <- if (candidate$periodic[[1L]]) J + 1L else J + 4L
  knots <- seq(0, 1, length.out = knot_count)
  common <- Splinets::splinet(
    knots,
    degree = 3L,
    periodic = candidate$periodic[[1L]],
    norm = TRUE)
  common_gram <- Splinets::gramian(common$os)
  if (!identical(dim(common_gram), c(J, J))) {
    stop("Candidate did not produce the registered J-dimensional system.")
  }

  second_derivative <- Splinets::deriva(Splinets::deriva(common$os))
  P2 <- Splinets::gramian(second_derivative)
  P2 <- (P2 + t(P2)) / 2
  eig <- eigen(P2, symmetric = TRUE)
  d <- pmax(eig$values, 0)
  d_scale <- stats::median(d)
  if (!is.finite(d_scale) || d_scale <= 0) {
    stop("Invalid second-derivative Gram spectrum.")
  }
  spectral_weight <- (1 + candidate$alpha[[1L]] * d / d_scale)^(-0.5)
  coefficient_map <- eig$vectors %*%
    (spectral_weight * t(eig$vectors))

  for (draw_id in seq_len(n_draw)) {
    Q_blocks <- vector("list", n_block)
    for (block_id in seq_len(n_block)) {
      gaussian <- matrix(stats::rnorm(J * M), J, M)
      weighted <- coefficient_map %*% gaussian
      qr_weighted <- qr(weighted)
      if (qr_weighted$rank < M) {
        stop("Unexpected rank-deficient weighted Gaussian draw.")
      }
      Q <- standardize_qr_sign(qr_weighted, M)
      Q_blocks[[block_id]] <- Q
      combined <- Splinets::lincomb(common$os, t(Q), reduced = TRUE)
      continuous_gram <- Splinets::gramian(combined)
      evaluated <- Splinets::evspline(combined, x = t_grid)
      phi <- evaluated[, -1L, drop = FALSE]
      grid_gram <- crossprod(phi, w_grid * phi)
      projection <- relative_projection_floor(phi)
      roughness <- diag(t(Q) %*% P2 %*% Q)
      edge <- boundary_energy(phi)

      block_pos <- block_pos + 1L
      block_rows[[block_pos]] <- data.frame(
        candidate_id = candidate$candidate_id[[1L]],
        periodic = candidate$periodic[[1L]],
        alpha = candidate$alpha[[1L]],
        draw_id = draw_id,
        block_id = block_id,
        gram_error_continuous =
          max(abs(continuous_gram - diag(M))),
        gram_error_grid = max(abs(grid_gram - diag(M))),
        roughness_mean = mean(roughness),
        roughness_max = max(roughness),
        boundary_energy_mean = mean(edge),
        boundary_energy_max = max(edge),
        projection_floor_mean = mean(projection),
        projection_floor_max = max(projection),
        stringsAsFactors = FALSE
      )
    }

    for (block_1 in seq_len(n_block - 1L)) {
      for (block_2 in (block_1 + 1L):n_block) {
        canonical <- svd(
          crossprod(Q_blocks[[block_1]], Q_blocks[[block_2]]),
          nu = 0,
          nv = 0)$d
        pair_pos <- pair_pos + 1L
        pair_rows[[pair_pos]] <- data.frame(
          candidate_id = candidate$candidate_id[[1L]],
          periodic = candidate$periodic[[1L]],
          alpha = candidate$alpha[[1L]],
          draw_id = draw_id,
          block_1 = block_1,
          block_2 = block_2,
          max_canonical_correlation = max(canonical),
          min_canonical_correlation = min(canonical),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  message("completed candidate ", candidate$candidate_id[[1L]])
}

block_metrics <- do.call(rbind, block_rows)
pair_metrics <- do.call(rbind, pair_rows)

role_metrics <- c(
  "roughness_mean", "boundary_energy_mean", "projection_floor_mean")
exchangeability <- do.call(rbind, lapply(
  candidates$candidate_id,
  function(candidate_id) {
    candidate_blocks <- block_metrics[
      block_metrics$candidate_id == candidate_id, , drop = FALSE]
    do.call(rbind, lapply(role_metrics, function(metric) {
      block_means <- tapply(
        candidate_blocks[[metric]], candidate_blocks$block_id, mean)
      pooled_sd <- stats::sd(candidate_blocks[[metric]])
      standardized_range <- if (is.finite(pooled_sd) && pooled_sd > 0) {
        diff(range(block_means)) / pooled_sd
      } else {
        0
      }
      data.frame(
        candidate_id = candidate_id,
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
  }))

candidate_summary <- do.call(rbind, lapply(
  candidates$candidate_id,
  function(candidate_id) {
    b <- block_metrics[block_metrics$candidate_id == candidate_id,
                       , drop = FALSE]
    p <- pair_metrics[pair_metrics$candidate_id == candidate_id,
                      , drop = FALSE]
    data.frame(
      candidate_id = candidate_id,
      periodic = b$periodic[[1L]],
      alpha = b$alpha[[1L]],
      continuous_gram_error_max = max(b$gram_error_continuous),
      grid_gram_error_max = max(b$gram_error_grid),
      roughness_median = stats::median(b$roughness_mean),
      boundary_energy_median = stats::median(b$boundary_energy_mean),
      projection_floor_median = stats::median(b$projection_floor_mean),
      projection_floor_q95 = unname(stats::quantile(
        b$projection_floor_mean, 0.95)),
      cross_canonical_median = stats::median(
        p$max_canonical_correlation),
      cross_canonical_q95 = unname(stats::quantile(
        p$max_canonical_correlation, 0.95)),
      cross_canonical_gt_0_95_rate = mean(
        p$max_canonical_correlation > 0.95),
      cross_canonical_max = max(p$max_canonical_correlation),
      all_metrics_finite =
        all(vapply(
          b[vapply(b, is.numeric, logical(1))],
          function(x) all(is.finite(x)), logical(1))) &&
        all(vapply(
          p[vapply(p, is.numeric, logical(1))],
          function(x) all(is.finite(x)), logical(1))),
      role_exchangeability = all(exchangeability$passed[
        exchangeability$candidate_id == candidate_id]),
      stringsAsFactors = FALSE
    )
  }))

candidate_summary$gram_pass <-
  candidate_summary$continuous_gram_error_max <= 1e-8
candidate_summary$projection_pass <-
  candidate_summary$projection_floor_median <= 0.25 &
  candidate_summary$projection_floor_q95 <= 0.60
candidate_summary$overlap_pass <-
  candidate_summary$cross_canonical_q95 <= 0.98 &
  candidate_summary$cross_canonical_gt_0_95_rate <= 0.10 &
  candidate_summary$cross_canonical_max < 1 - 1e-8 &
  candidate_summary$cross_canonical_median > 0.10
candidate_summary$passed <- with(
  candidate_summary,
  gram_pass & all_metrics_finite & role_exchangeability &
    projection_pass & overlap_pass)

eligible <- candidate_summary[candidate_summary$passed, , drop = FALSE]
if (nrow(eligible) > 0L) {
  eligible <- eligible[order(
    eligible$alpha,
    eligible$periodic,
    eligible$candidate_id), , drop = FALSE]
  selected_candidate <- eligible$candidate_id[[1L]]
  decision <- if (eligible$periodic[[1L]]) {
    "GEOMETRY_PASS_PERIODIC_REQUIRES_SCIENTIFIC_BOUNDARY_DECISION"
  } else {
    "GEOMETRY_PASS_NONPERIODIC_ELIGIBLE_FOR_IMPLEMENTATION"
  }
} else {
  selected_candidate <- NA_character_
  decision <- "NO_CANDIDATE_PASSED_STOP_WITHOUT_IMPLEMENTATION"
}

decision_table <- data.frame(
  audit_id = "IID_COMMON_SPLINET_STIEFEL_GEOMETRY_V2_20260922",
  audit_seed = audit_seed,
  n_draw_per_candidate = n_draw,
  selected_candidate = selected_candidate,
  decision = decision,
  formal_data_generated = FALSE,
  model_fits_started = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(candidates, file.path(output_dir, "CANDIDATES.csv"),
                 row.names = FALSE)
utils::write.csv(block_metrics, file.path(output_dir, "BLOCK_METRICS.csv"),
                 row.names = FALSE)
utils::write.csv(pair_metrics, file.path(output_dir, "PAIR_METRICS.csv"),
                 row.names = FALSE)
utils::write.csv(exchangeability,
                 file.path(output_dir, "ROLE_EXCHANGEABILITY.csv"),
                 row.names = FALSE)
utils::write.csv(candidate_summary,
                 file.path(output_dir, "CANDIDATE_SUMMARY.csv"),
                 row.names = FALSE)
utils::write.csv(decision_table, file.path(output_dir, "DECISION.csv"),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "SESSION_INFO.txt"))
writeLines(c(
  decision,
  paste0("selected_candidate=", ifelse(
    is.na(selected_candidate), "NA", selected_candidate)),
  "formal_data_generated=FALSE",
  "model_fits_started=FALSE"),
  file.path(output_dir, "AUDIT_STATUS.txt"))

cat(decision, "\n", sep = "")
if (!is.na(selected_candidate)) {
  cat("selected_candidate=", selected_candidate, "\n", sep = "")
}
