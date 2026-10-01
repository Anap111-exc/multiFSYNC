#!/usr/bin/env Rscript

# Rebuild the study-2 specific loading coordinate at one frozen sweep.
# Read-only: no optimizer, no update function, no truth bundle.
options(warn = 2, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(
  "Usage: Rscript smoke_stage7h_loading_driver_20260929_v1.R ",
  "STAGE7C_ROOT STAGE7G_ROOT PREFLIGHT_ROOT OUTPUT_ROOT FIT_ID",
  call. = FALSE)
stage7c <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g <- normalizePath(args[[2L]], mustWork = TRUE)
preflight <- normalizePath(args[[3L]], mustWork = TRUE)
out <- args[[4L]]
fit_id <- args[[5L]]
if (dir.exists(out)) stop("Refusing existing output directory: ", out)
if (!file.exists(file.path(preflight, "PREFLIGHT_PASS.txt")))
  stop("Stage 7H preflight marker missing.")

binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
row <- binding[binding$fit_id == fit_id, , drop = FALSE]
man <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
if (nrow(row) != 1L || nrow(man) != 1L ||
    !identical(row$data_id, man$data_id)) stop("Unbound fit ID.")
obs_path <- file.path(stage7c, "data", row$data_id[[1L]],
                      "observation_bundle.rds")
fit_path <- file.path(stage7g, "fits", fit_id, "fit.rds")
snapshot_path <- file.path(stage7g, "fits", fit_id,
                           "RAW_STATE_SNAPSHOTS.rds")
trace_path <- file.path(stage7g, "fits", fit_id, "DRIVER_TRACE.csv")
sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(line) != 1L) stop("sha256sum failed: ", path)
  strsplit(line, "[[:space:]]+")[[1L]][[1L]]
}
if (!identical(sha256(obs_path), man$observation_sha256[[1L]]) ||
    !identical(sha256(fit_path), row$fit_sha256[[1L]]) ||
    !identical(sha256(snapshot_path), row$snapshot_sha256[[1L]]))
  stop("Bound input SHA256 mismatch.")

.libPaths(c("/root/multiFSYNC-stage7g-lib-20260928-v1", .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
if (normalizePath(system.file(package = "multiFSYNC")) !=
    "/root/multiFSYNC-stage7g-lib-20260928-v1/multiFSYNC")
  stop("Wrong fitted basis package.")

fit <- readRDS(fit_path)
obs <- readRDS(obs_path)
snapshots <- readRDS(snapshot_path)
if (length(snapshots) != 500L ||
    !identical(vapply(snapshots, `[[`, integer(1L), "iteration"),
               as.integer(0:499))) stop("Snapshot index changed.")
t <- 50L
previous <- snapshots[[t]]        # iteration t-1 = 49
current <- snapshots[[t + 1L]]    # iteration t = 50
rm(snapshots)
invisible(gc())

time_flat <- unlist(obs$time_obs, recursive = FALSE)
basis <- multiFSYNC:::get_grid_objects(
  time_flat, fit$K, n_g = fit$n_g, time_g = fit$time_g,
  format_univ = TRUE)
C <- basis$C[(fit$n_s[[1L]] + 1L):sum(fit$n_s)]
Y <- obs$Y[[2L]]
Z <- obs$Z[[2L]]
p <- fit$p
n <- fit$n_s[[2L]]
K_total <- fit$K_total
M <- fit$M_s[[2L]][[1L]]
c_val <- 1 / current$temperature
if (!is.matrix(Z) || nrow(Z) != n || ncol(Z) != fit$d ||
    length(C) != n || length(Y) != n || M != 2L ||
    !is.finite(c_val) || c_val <= 0) stop("Design contract changed.")

mean_coeff <- do.call(cbind, current$mu_q_nu_mu[[2L]])
beta_coeff <- lapply(seq_len(fit$d), function(r)
  do.call(cbind, lapply(current$mu_q_nu_beta, `[[`, r)))
psi <- current$mu_q_nu_psi[[2L]][[1L]]
psi_cov <- current$Sigma_q_nu_psi[[2L]][[1L]]
xi <- current$mu_q_xi[[2L]][[1L]]
xi_cov <- current$Sigma_q_xi[[2L]][[1L]]
a <- current$mu_q_a[, 1L]
phi <- current$mu_q_nu_phi[[1L]]
zeta <- current$mu_q_zeta[[2L]][[1L]]
tau <- previous$mu_q_recip_sigsq_eps[2L, ]
if (!all(dim(mean_coeff) == c(K_total, p)) ||
    !all(dim(psi) == c(K_total, M)) ||
    !all(dim(xi) == c(n, M)) || length(a) != p ||
    length(tau) != p) stop("State shape mismatch.")

h <- numeric(p)
U <- 0
U_mean <- U_score_cov <- U_function_cov <- U_joint_cov <- 0
data_envelope <- 0
max_design_residual_error <- 0

for (i in seq_len(n)) {
  C_i <- C[[i]]
  cp_C_i <- fit$list_cp_C[[2L]][[i]]
  cpY_i <- fit$list_cp_C_Y[[2L]][[i]]
  Y_i <- do.call(cbind, Y[[i]])
  mean_i <- C_i %*% mean_coeff
  beta_i <- matrix(0, nrow(C_i), p)
  for (r in seq_len(fit$d))
    beta_i <- beta_i + Z[i, r] * (C_i %*% beta_coeff[[r]])
  shared_time_i <- as.vector(C_i %*% phi %*% zeta[i, ])
  shared_i <- tcrossprod(shared_time_i, a)
  residual_i <- Y_i - mean_i - beta_i - shared_i
  ct_res_direct <- crossprod(C_i, residual_i)
  ct_res_suff <- cpY_i - cp_C_i %*% mean_coeff
  for (r in seq_len(fit$d))
    ct_res_suff <- ct_res_suff -
      Z[i, r] * (cp_C_i %*% beta_coeff[[r]])
  ct_res_suff <- ct_res_suff -
    tcrossprod(as.vector(cp_C_i %*% phi %*% zeta[i, ]), a)
  max_design_residual_error <- max(max_design_residual_error,
                                   max(abs(ct_res_direct - ct_res_suff)))

  g_i <- as.vector(C_i %*% psi %*% xi[i, ])
  coordinate_driver_i <- as.vector(crossprod(g_i, residual_i))
  h_i <- tau * coordinate_driver_i
  h <- h + h_i
  data_envelope <- data_envelope + sum(abs(h_i))

  H_mean <- crossprod(psi, cp_C_i %*% psi)
  H_function_cov <- diag(vapply(seq_len(M), function(m)
    sum(diag(cp_C_i %*% psi_cov[[m]])), numeric(1L)), M)
  H <- H_mean + H_function_cov
  xi_mean_outer <- tcrossprod(xi[i, ])
  U_mean <- U_mean + sum(diag(H_mean %*% xi_mean_outer))
  U_score_cov <- U_score_cov + sum(diag(H_mean %*% xi_cov[[i]]))
  U_function_cov <- U_function_cov +
    sum(diag(H_function_cov %*% xi_mean_outer))
  U_joint_cov <- U_joint_cov +
    sum(diag(H_function_cov %*% xi_cov[[i]]))
  U <- U + sum(diag(H %*% (xi_cov[[i]] + xi_mean_outer)))

}

if (max(abs(U - (U_mean + U_score_cov +
                 U_function_cov + U_joint_cov))) > 1e-8)
  stop("Loading precision decomposition mismatch.")
precision <- c_val * (tau * U + 1)
v <- 1 / precision
m <- c_val * v * h
log_omega <- previous$mu_q_log_omega_b[[2L]]
log_1_omega <- previous$mu_q_log_1_omega_b[[2L]]
prior_log_odds <- if (is.matrix(log_omega))
  log_omega[, 1L] - log_1_omega[, 1L] else
    rep(log_omega[[1L]] - log_1_omega[[1L]], p)
logit_prior <- c_val * prior_log_odds
logit_signal <- 0.5 * m^2 / v
logit_variance <- 0.5 * log(v)
logit_measure <- 0.5 * (1 - c_val) * log(2 * pi)
rho <- stats::plogis(logit_prior + logit_signal +
                       logit_variance + logit_measure)
b <- rho * m
second_b <- rho * (v + m^2)

trace <- utils::read.csv(trace_path, stringsAsFactors = FALSE)
reference <- trace[trace$iteration == t & trace$phase == "main" &
                     trace$block == "specific_loading" &
                     trace$study == 2L & trace$factor == 1L, , drop = FALSE]
if (nrow(reference) != 1L) stop("Expected driver trace row missing.")
metrics <- c(
  data_driver_norm = sqrt(sum(h^2)),
  data_envelope = data_envelope,
  data_coherence = sqrt(sum(h^2)) / max(data_envelope,
                                        .Machine$double.eps),
  precision_trace = sum(precision),
  precision_frobenius = sqrt(sum(precision^2)),
  prior_precision_trace = c_val * p,
  driver_precision_ratio = sqrt(sum(h^2)) /
    max(sqrt(sum(precision^2)), .Machine$double.eps),
  posterior_mean_norm = sqrt(sum(b^2)),
  slab_mean_norm = sqrt(sum(m^2)),
  slab_variance_mean = mean(v),
  ppi_min = min(rho), ppi_mean = mean(rho), ppi_max = max(rho),
  log_prior_odds_mean = mean(prior_log_odds))
saved <- list(slab_mean = current$mu_q_normal_b[[2L]][, 1L],
              slab_variance = current$Sigma_q_normal_b[[2L]][, 1L],
              ppi = current$mu_q_gamma_b[[2L]][, 1L],
              effective_loading = current$mu_q_b_specific[[2L]][, 1L],
              loading_second = current$term_b_specific[[2L]][, 1L])
field_error <- c(
  slab_mean = max(abs(m - saved$slab_mean)),
  slab_variance = max(abs(v - saved$slab_variance)),
  ppi = max(abs(rho - saved$ppi)),
  effective_loading = max(abs(b - saved$effective_loading)),
  loading_second = max(abs(second_b - saved$loading_second)))
trace_error <- vapply(names(metrics), function(nm)
  abs(metrics[[nm]] - as.numeric(reference[[nm]][[1L]])), numeric(1L))
tol <- function(ref) 1e-7 + 1e-7 * max(1, abs(ref))
if (any(!is.finite(field_error)) || any(!is.finite(trace_error)) ||
    max_design_residual_error > 1e-7 ||
    any(vapply(names(field_error), function(nm)
      field_error[[nm]] > tol(max(abs(saved[[nm]]))), logical(1L))) ||
    any(vapply(names(metrics), function(nm)
      trace_error[[nm]] > tol(as.numeric(reference[[nm]][[1L]])),
      logical(1L)))) {
  stop("Stage 7H smoke loading driver parity failed: design=",
       signif(max_design_residual_error, 7L),
       " state=", paste(names(field_error), signif(field_error, 7L),
                        collapse = ";"),
       " trace=", paste(names(trace_error), signif(trace_error, 7L),
                        collapse = ";"))
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Cannot create output directory.")
utils::write.csv(data.frame(fit_id = fit_id, iteration = t,
                            max_design_residual_error =
                              max_design_residual_error,
                            U = U, U_mean = U_mean,
                            U_score_cov = U_score_cov,
                            U_function_cov = U_function_cov,
                            U_joint_cov = U_joint_cov,
                            fitted_shared_b2_abs_cosine =
                              abs(sum(a * b)) /
                              sqrt(sum(a^2) * sum(b^2)),
                            stringsAsFactors = FALSE),
                 file.path(out, "SMOKE_SUMMARY.csv"), row.names = FALSE)
utils::write.csv(data.frame(field = names(field_error),
                            max_abs_difference = as.numeric(field_error)),
                 file.path(out, "SMOKE_STATE_PARITY.csv"), row.names = FALSE)
utils::write.csv(data.frame(field = names(trace_error),
                            abs_difference = as.numeric(trace_error)),
                 file.path(out, "SMOKE_TRACE_PARITY.csv"), row.names = FALSE)
writeLines("STAGE7H_LOADING_SMOKE_PASS",
           file.path(out, "SMOKE_PASS.txt"), useBytes = TRUE)
cat("STAGE7H_LOADING_SMOKE_PASS\n")
