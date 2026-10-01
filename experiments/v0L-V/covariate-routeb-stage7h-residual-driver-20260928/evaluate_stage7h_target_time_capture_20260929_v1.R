#!/usr/bin/env Rscript

# Exploratory truth-assisted time-capture check after Stage 7H blind freeze.
# It reads only frozen states and observation bundles; it does not fit.
options(warn = 2, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(
  "Usage: Rscript evaluate_stage7h_target_time_capture_20260929_v1.R ",
  "STAGE7C_ROOT STAGE7G_ROOT BLIND_ROOT OUTPUT_ROOT", call. = FALSE)
stage7c <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g <- normalizePath(args[[2L]], mustWork = TRUE)
blind <- normalizePath(args[[3L]], mustWork = TRUE)
out <- args[[4L]]
if (dir.exists(out)) stop("Refusing existing output directory: ", out)
sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(line) != 1L) stop("sha256sum failed: ", path)
  strsplit(line, "[[:space:]]+")[[1L]][[1L]]
}
hash_lines <- readLines(file.path(blind, "BLIND_TABLES.sha256"),
                        warn = FALSE)
if (length(hash_lines) != 6L ||
    !all(c("mode=all", "point_count=20") %in%
         readLines(file.path(blind, "BLIND_COMPLETE.txt"), warn = FALSE)))
  stop("Blind freeze marker missing.")
for (line in hash_lines) {
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  if (length(fields) != 2L ||
      !identical(sha256(file.path(blind, fields[[2L]])), fields[[1L]]))
    stop("Blind table hash mismatch.")
}
.libPaths(c("/root/multiFSYNC-stage7g-lib-20260928-v1", .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
if (normalizePath(system.file(package = "multiFSYNC")) !=
    "/root/multiFSYNC-stage7g-lib-20260928-v1/multiFSYNC")
  stop("Wrong frozen source package.")
manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
if (nrow(manifest) != 3L || nrow(binding) != 3L ||
    !identical(as.character(manifest$fit_id),
               as.character(binding$fit_id)))
  stop("Input binding changed.")

rows <- list()
index <- 0L
for (k in seq_len(3L)) {
  fit_id <- manifest$fit_id[[k]]
  data_id <- manifest$data_id[[k]]
  obs_path <- file.path(stage7c, "data", data_id,
                        "observation_bundle.rds")
  truth_path <- file.path(stage7c, "data", data_id,
                          "sealed_truth", "truth_bundle.rds")
  fit_path <- file.path(stage7g, "fits", fit_id, "fit.rds")
  snap_path <- file.path(stage7g, "fits", fit_id,
                         "RAW_STATE_SNAPSHOTS.rds")
  if (!identical(sha256(obs_path), manifest$observation_sha256[[k]]) ||
      !identical(sha256(fit_path), binding$fit_sha256[[k]]) ||
      !identical(sha256(snap_path), binding$snapshot_sha256[[k]]) ||
      !identical(sha256(truth_path), binding$truth_sha256[[k]]))
    stop("Bound input SHA256 mismatch: ", fit_id)
  obs <- readRDS(obs_path)
  truth_bundle <- readRDS(truth_path)
  fit <- readRDS(fit_path)
  snapshots <- readRDS(snap_path)
  if (length(snapshots) != 500L ||
      !identical(truth_bundle$data_id, data_id))
    stop("Snapshot/truth contract changed.")
  true_b <- as.numeric(truth_bundle$data$true_params$b_true[[2L]][, 1L])
  time_flat <- unlist(obs$time_obs, recursive = FALSE)
  basis <- multiFSYNC:::get_grid_objects(
    time_flat, fit$K, n_g = fit$n_g, time_g = fit$time_g,
    format_univ = TRUE)
  C <- basis$C[(fit$n_s[[1L]] + 1L):sum(fit$n_s)]
  n <- fit$n_s[[2L]]
  p <- fit$p
  Z <- obs$Z[[2L]]
  times <- if (k < 3L) c(20L, 25L, 40L, 50L, 55L, 60L, 67L)
           else c(20L, 25L, 40L, 50L, 55L, 60L)
  for (t in times) {
    state <- snapshots[[t + 1L]]
    before <- snapshots[[t]]
    mean_coeff <- do.call(cbind, state$mu_q_nu_mu[[2L]])
    beta_coeff <- lapply(seq_len(fit$d), function(r)
      do.call(cbind, lapply(state$mu_q_nu_beta, `[[`, r)))
    phi <- state$mu_q_nu_phi[[1L]]
    zeta <- state$mu_q_zeta[[2L]][[1L]]
    a <- state$mu_q_a[, 1L]
    psi <- state$mu_q_nu_psi[[2L]][[1L]]
    xi <- state$mu_q_xi[[2L]][[1L]]
    tau <- before$mu_q_recip_sigsq_eps[2L, ]
    for (metric in c("raw_target", "tau_weighted_target")) {
      target <- if (metric == "raw_target") true_b else tau * true_b
      target <- target / sqrt(sum(target^2))
      total_resid_energy <- total_span_energy <- total_actual_energy <- 0
      g_norm_sum <- 0
      rank_min <- Inf
      zero_g_count <- 0L
      for (i in seq_len(n)) {
        C_i <- C[[i]]
        Y_i <- do.call(cbind, obs$Y[[2L]][[i]])
        mean_i <- C_i %*% mean_coeff
        beta_i <- matrix(0, nrow(C_i), p)
        for (r in seq_len(fit$d))
          beta_i <- beta_i + Z[i, r] * (C_i %*% beta_coeff[[r]])
        shared_time_i <- as.vector(C_i %*% phi %*% zeta[i, ])
        residual_i <- Y_i - mean_i - beta_i -
          tcrossprod(shared_time_i, a)
        residual_target_i <- as.vector(residual_i %*% target)
        Gi <- C_i %*% psi
        qrGi <- qr(Gi, tol = 1e-8)
        rank_min <- min(rank_min, qrGi$rank)
        total_resid_energy <- total_resid_energy +
          sum(residual_target_i^2)
        if (qrGi$rank > 0L)
          total_span_energy <- total_span_energy +
            sum(qr.fitted(qrGi, residual_target_i)^2)
        g_i <- as.vector(Gi %*% xi[i, ])
        g_sq <- sum(g_i^2)
        g_norm_sum <- g_norm_sum + g_sq
        if (g_sq > 1e-24) {
          total_actual_energy <- total_actual_energy +
            sum(residual_target_i * g_i)^2 / g_sq
        } else {
          zero_g_count <- zero_g_count + 1L
        }
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        fit_id = fit_id, data_id = data_id,
        iteration = t, target_metric = metric,
        true_B2_residual_energy = total_resid_energy,
        current_function_span_capture =
          total_span_energy / max(total_resid_energy,
                                  .Machine$double.eps),
        current_score_function_capture =
          total_actual_energy / max(total_resid_energy,
                                    .Machine$double.eps),
        current_score_function_squared_norm = g_norm_sum,
        function_span_min_rank = rank_min,
        zero_score_function_subjects = zero_g_count,
        stringsAsFactors = FALSE)
    }
  }
  rm(obs, truth_bundle, fit, snapshots, basis, C)
  invisible(gc())
}
result <- do.call(rbind, rows)
if (nrow(result) != 40L) stop("Output count changed: ", nrow(result))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result,
                 file.path(out, "TRUE_B2_TIME_CAPTURE_EXPLORATORY.csv"),
                 row.names = FALSE)
writeLines(c("STAGE7H_TARGET_TIME_CAPTURE_COMPLETE",
             "status=posthoc_exploratory_not_preregistered_threshold",
             paste0("blind_hash_list_sha256=",
                    sha256(file.path(blind, "BLIND_TABLES.sha256")))),
           file.path(out, "TIME_CAPTURE_COMPLETE.txt"),
           useBytes = TRUE)
cat("STAGE7H_TIME_CAPTURE_COMPLETE rows=", nrow(result), "\n",
    sep = "")
