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
  script_dir, "acceptance", "windows", "route_b_K_audit_v1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Splinets", quietly = TRUE)) {
  stop("Splinets is required for this audit.", call. = FALSE)
}
suppressPackageStartupMessages(library(splines))
source(file.path(repo_root, "R", "utils_multi.R"), local = .GlobalEnv)
source(file.path(repo_root, "R", "generate_data_structured.R"),
       local = .GlobalEnv)

audit_seed <- 20260925L
n_draw <- 1000L
n_block <- 3L
M <- 2L
J <- 8L
alpha <- 3
K_candidates <- 5:10
n_grid <- 401L

t_grid <- seq(0, 1, length.out = n_grid)
dt <- diff(t_grid)
w_grid <- c(
  dt[[1L]] / 2,
  (dt[-1L] + dt[-length(dt)]) / 2,
  dt[[length(dt)]] / 2)

fit_spaces <- lapply(K_candidates, function(K) {
  int_knots <- seq(0, 1, length.out = K)[-c(1L, K)]
  C <- cbind(
    1,
    t_grid,
    ZOSull(t_grid, range.x = c(0, 1), intKnots = int_knots))
  list(
    K = K,
    n_internal_knots = K - 2L,
    total_coefficient_dimension = ncol(C),
    C = C,
    C_weighted = sqrt(w_grid) * C)
})

knots <- seq(0, 1, length.out = J + 4L)
common <- Splinets::splinet(
  knots, degree = 3L, periodic = FALSE, norm = TRUE)
common_gram <- Splinets::gramian(common$os)
if (!identical(dim(common_gram), c(J, J)) ||
    max(abs(common_gram - diag(J))) > 1e-8) {
  stop("Invalid common Splinets system.", call. = FALSE)
}

second_derivative <- Splinets::deriva(Splinets::deriva(common$os))
P2 <- Splinets::gramian(second_derivative)
P2 <- (P2 + t(P2)) / 2
eig <- eigen(P2, symmetric = TRUE)
d <- pmax(eig$values, 0)
d_scale <- stats::median(d)
spectral_weight <- (1 + alpha * d / d_scale)^(-0.5)
coefficient_map <- eig$vectors %*%
  (spectral_weight * t(eig$vectors))

projection_rows <- vector(
  "list", n_draw * n_block * length(K_candidates))
pair_rows <- vector("list", n_draw * choose(n_block, 2L))
projection_pos <- 0L
pair_pos <- 0L
continuous_gram_error_max <- 0

standardize_qr_sign <- function(qr_object) {
  Q <- qr.Q(qr_object)[, seq_len(M), drop = FALSE]
  qr_sign <- sign(diag(qr.R(qr_object))[seq_len(M)])
  qr_sign[qr_sign == 0] <- 1
  sweep(Q, 2L, qr_sign, FUN = "*")
}

set.seed(audit_seed)
for (draw_id in seq_len(n_draw)) {
  Q_blocks <- vector("list", n_block)
  for (block_id in seq_len(n_block)) {
    gaussian <- matrix(stats::rnorm(J * M), J, M)
    qr_weighted <- qr(coefficient_map %*% gaussian)
    if (qr_weighted$rank < M) {
      stop("Unexpected rank-deficient weighted Gaussian draw.")
    }
    Q <- standardize_qr_sign(qr_weighted)
    Q_blocks[[block_id]] <- Q
    combined <- Splinets::lincomb(common$os, t(Q), reduced = TRUE)
    continuous_gram <- Splinets::gramian(combined)
    continuous_gram_error_max <- max(
      continuous_gram_error_max,
      max(abs(continuous_gram - diag(M))))
    phi <- Splinets::evspline(combined, x = t_grid)[, -1L, drop = FALSE]

    for (fit_space in fit_spaces) {
      coefficients <- qr.solve(
        fit_space$C_weighted, sqrt(w_grid) * phi)
      residual <- phi - fit_space$C %*% coefficients
      component_error <- sqrt(
        colSums(w_grid * residual^2) / colSums(w_grid * phi^2))
      projection_pos <- projection_pos + 1L
      projection_rows[[projection_pos]] <- data.frame(
        draw_id = draw_id,
        block_id = block_id,
        K = fit_space$K,
        n_internal_knots = fit_space$n_internal_knots,
        total_coefficient_dimension = fit_space$total_coefficient_dimension,
        projection_floor_mean = mean(component_error),
        projection_floor_max = max(component_error),
        unexplained_energy_mean = mean(component_error^2),
        stringsAsFactors = FALSE
      )
    }
  }

  for (block_1 in seq_len(n_block - 1L)) {
    for (block_2 in (block_1 + 1L):n_block) {
      canonical <- svd(
        crossprod(Q_blocks[[block_1]], Q_blocks[[block_2]]),
        nu = 0,
        nv = 0)$d
      pair_pos <- pair_pos + 1L
      pair_rows[[pair_pos]] <- data.frame(
        draw_id = draw_id,
        block_1 = block_1,
        block_2 = block_2,
        max_canonical_correlation = max(canonical),
        stringsAsFactors = FALSE
      )
    }
  }
}

projection_metrics <- do.call(rbind, projection_rows)
pair_metrics <- do.call(rbind, pair_rows)
overlap_summary <- data.frame(
  cross_canonical_median = stats::median(
    pair_metrics$max_canonical_correlation),
  cross_canonical_q95 = unname(stats::quantile(
    pair_metrics$max_canonical_correlation, 0.95)),
  cross_canonical_gt_0_95_rate = mean(
    pair_metrics$max_canonical_correlation > 0.95),
  cross_canonical_max = max(pair_metrics$max_canonical_correlation),
  stringsAsFactors = FALSE
)
overlap_pass <- with(
  overlap_summary,
  cross_canonical_q95 <= 0.98 &
    cross_canonical_gt_0_95_rate <= 0.10 &
    cross_canonical_max < 1 - 1e-8)

K_summary <- do.call(rbind, lapply(K_candidates, function(K) {
  x <- projection_metrics[projection_metrics$K == K, , drop = FALSE]
  data.frame(
    K = K,
    n_internal_knots = x$n_internal_knots[[1L]],
    total_coefficient_dimension = x$total_coefficient_dimension[[1L]],
    projection_floor_median = stats::median(x$projection_floor_mean),
    projection_floor_q95 = unname(stats::quantile(
      x$projection_floor_mean, 0.95)),
    unexplained_energy_median = stats::median(x$unexplained_energy_mean),
    unexplained_energy_q95 = unname(stats::quantile(
      x$unexplained_energy_mean, 0.95)),
    projection_not_exact_containment =
      max(x$projection_floor_max) > 1e-8,
    gram_pass = continuous_gram_error_max <= 1e-8,
    overlap_pass = overlap_pass,
    stringsAsFactors = FALSE
  )
}))
K_summary$projection_pass <-
  K_summary$projection_floor_median <= 0.25 &
  K_summary$projection_floor_q95 <= 0.60
K_summary$passed <- with(
  K_summary,
  projection_pass & projection_not_exact_containment &
    gram_pass & overlap_pass)

eligible <- K_summary[K_summary$passed, , drop = FALSE]
if (nrow(eligible) > 0L) {
  selected_K <- min(eligible$K)
  decision <- "MINIMUM_K_GEOMETRY_PASS"
} else {
  selected_K <- NA_integer_
  decision <- "NO_K_PASSED_STOP"
}
decision_table <- data.frame(
  audit_id = "ROUTE_B_NONPERIODIC_K_AUDIT_V1_20260922",
  audit_seed = audit_seed,
  n_draw = n_draw,
  J_truth = J,
  alpha = alpha,
  selected_K = selected_K,
  decision = decision,
  formal_data_generated = FALSE,
  model_fits_started = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  projection_metrics,
  file.path(output_dir, "PROJECTION_METRICS.csv"), row.names = FALSE)
utils::write.csv(
  pair_metrics,
  file.path(output_dir, "PAIR_METRICS.csv"), row.names = FALSE)
utils::write.csv(
  overlap_summary,
  file.path(output_dir, "OVERLAP_SUMMARY.csv"), row.names = FALSE)
utils::write.csv(
  K_summary,
  file.path(output_dir, "K_SUMMARY.csv"), row.names = FALSE)
utils::write.csv(
  decision_table,
  file.path(output_dir, "DECISION.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "SESSION_INFO.txt"))
writeLines(c(
  decision,
  paste0("selected_K=", ifelse(is.na(selected_K), "NA", selected_K)),
  "formal_data_generated=FALSE",
  "model_fits_started=FALSE"),
  file.path(output_dir, "AUDIT_STATUS.txt"))

cat(decision, "\n", sep = "")
cat("selected_K=", ifelse(is.na(selected_K), "NA", selected_K),
    "\n", sep = "")

