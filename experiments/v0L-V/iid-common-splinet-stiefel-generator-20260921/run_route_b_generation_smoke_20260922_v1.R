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
  script_dir, "acceptance", "windows", "route_b_generation_smoke_v1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(splines))
source(file.path(repo_root, "R", "utils_multi.R"), local = .GlobalEnv)
source(file.path(repo_root, "R", "generate_data_structured.R"),
       local = .GlobalEnv)

dat <- simulate_multi_study_structured(
  S = 4,
  n_s = rep(4L, 4L),
  p = 8,
  d = 2,
  L_f = 2,
  L_s = rep(1L, 4L),
  M_f = c(2L, 2L),
  M_s = replicate(4L, 2L, simplify = FALSE),
  K = 8,
  n_obs = 7,
  common_grid = FALSE,
  n_obs_range = c(6L, 9L),
  sigma_eps = 0.3,
  bool_sparse_loadings = FALSE,
  score_var_decay = TRUE,
  n_dense = 201,
  mean_amp = 0.6,
  beta_amp = 0.3,
  seed = 20260926L,
  truth_basis_seed = 20261926L,
  truth_basis_method = "iid_common_splinet_stiefel",
  truth_basis_dim = 8L)

tp <- dat$true_params
all_meta <- c(tp$phi_basis_meta,
              unlist(tp$psi_basis_meta, recursive = FALSE))
all_phi <- c(tp$phi_dense_list,
             unlist(tp$psi_dense_list, recursive = FALSE))

n_dense_smoke <- nrow(all_phi[[1L]])
t_grid <- seq(0, 1, length.out = n_dense_smoke)
dt <- diff(t_grid)
w <- c(dt[[1L]] / 2,
       (dt[-1L] + dt[-length(dt)]) / 2,
       dt[[length(dt)]] / 2)
C_fit <- cbind(
  1,
  t_grid,
  ZOSull(t_grid, range.x = c(0, 1), intKnots = tp$int_knots))
C_weighted <- sqrt(w) * C_fit

projection_error <- vapply(all_phi, function(phi) {
  coefficients <- qr.solve(C_weighted, sqrt(w) * phi)
  residual <- phi - C_fit %*% coefficients
  mean(sqrt(colSums(w * residual^2) / colSums(w * phi^2)))
}, numeric(1L))

cross_canonical <- numeric(0)
for (block_1 in seq_len(length(all_meta) - 1L)) {
  for (block_2 in (block_1 + 1L):length(all_meta)) {
    cross_canonical <- c(
      cross_canonical,
      max(svd(crossprod(
        all_meta[[block_1]]$stiefel_coefficients,
        all_meta[[block_2]]$stiefel_coefficients),
        nu = 0, nv = 0)$d))
  }
}

reconstruction_error <- 0
for (s in seq_len(tp$S)) {
  for (i in seq_len(tp$n_s[[s]])) {
    y <- do.call(cbind, dat$Y[[s]][[i]])
    reconstructed <- tp$signal_true_values[[s]][[i]] +
      tp$noise_true_values[[s]][[i]]
    reconstruction_error <- max(
      reconstruction_error, max(abs(y - reconstructed)))
  }
}

checks <- data.frame(
  check_id = c(
    "six_factor_blocks_supported",
    "route_b_metadata",
    "continuous_gram_exact",
    "response_reconstruction_exact",
    "irregular_observation_counts",
    "all_metrics_finite"),
  passed = c(
    length(all_meta) == 6L,
    identical(tp$truth_basis_method, "iid_common_splinet_stiefel") &&
      identical(tp$truth_basis_dim, 8L) &&
      identical(tp$truth_basis_smoothness_alpha, 3) &&
      all(vapply(all_meta, function(x) {
        !x$periodic && identical(x$smoothness_alpha, 3)
      }, logical(1L))),
    max(vapply(all_meta, function(x) {
      max(abs(x$gram - diag(nrow(x$gram))))
    }, numeric(1L))) <= 1e-8,
    reconstruction_error == 0,
    all(unlist(tp$n_obs_by_subject, use.names = FALSE) >= 6L &
        unlist(tp$n_obs_by_subject, use.names = FALSE) <= 9L),
    all(is.finite(c(projection_error, cross_canonical)))),
  detail = c(
    paste0("blocks=", length(all_meta)),
    paste0("J=", tp$truth_basis_dim,
           ";alpha=", tp$truth_basis_smoothness_alpha,
           ";K=", tp$K),
    sprintf("max=%.3e", max(vapply(all_meta, function(x) {
      max(abs(x$gram - diag(nrow(x$gram))))
    }, numeric(1L)))),
    sprintf("max=%.3e", reconstruction_error),
    paste(range(unlist(tp$n_obs_by_subject, use.names = FALSE)),
          collapse = ":"),
    paste0("projection_n=", length(projection_error),
           ";cross_n=", length(cross_canonical))),
  stringsAsFactors = FALSE)

metrics <- data.frame(
  metric = c(
    "projection_error_median",
    "projection_error_max",
    "cross_canonical_median",
    "cross_canonical_max",
    "response_reconstruction_error"),
  value = c(
    stats::median(projection_error),
    max(projection_error),
    stats::median(cross_canonical),
    max(cross_canonical),
    reconstruction_error),
  stringsAsFactors = FALSE)

utils::write.csv(checks, file.path(output_dir, "SMOKE_CHECKS.csv"),
                 row.names = FALSE)
utils::write.csv(metrics, file.path(output_dir, "SMOKE_METRICS.csv"),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "SESSION_INFO.txt"))

if (!all(checks$passed)) {
  writeLines("FAIL", file.path(output_dir, "SMOKE_STATUS.txt"))
  stop("Route-B generation smoke failed.", call. = FALSE)
}
writeLines(c(
  "PASS",
  "formal_data_generated=FALSE",
  "model_fits_started=FALSE"),
  file.path(output_dir, "SMOKE_STATUS.txt"))
cat("ROUTE_B_GENERATION_SMOKE_PASS\n")
