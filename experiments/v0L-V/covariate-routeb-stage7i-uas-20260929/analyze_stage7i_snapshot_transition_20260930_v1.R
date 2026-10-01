#!/usr/bin/env Rscript

# Read-only reconstruction from the frozen sweep-50 injected states and the
# saved end-of-sweep states 51:60. No CAVI update or sealed truth is invoked.
# Usage: Rscript SCRIPT RESULTS_PARENT PACKAGE_LIB OUTPUT_DIR

options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected RESULTS_PARENT PACKAGE_LIB OUTPUT_DIR.", call. = FALSE)
parent <- normalizePath(args[[1L]], mustWork = TRUE)
package_lib <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (file.exists(output) || dir.exists(output))
  stop("Refusing existing output: ", output, call. = FALSE)
.libPaths(unique(c(package_lib, .libPaths())))
stopifnot(requireNamespace("digest", quietly = TRUE),
          requireNamespace("multiFSYNC", quietly = TRUE))
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                           check.names = FALSE)
write_csv <- function(value, filename) {
  path <- file.path(output, filename)
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}

freeze <- file.path(parent, "stage7i_blind_freeze_20260929_v1")
marker_path <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
marker <- readLines(marker_path, warn = FALSE)
assert(identical(marker[[1L]], "STAGE7I_BLIND_ENDPOINTS_FROZEN") &&
         "truth_bundle_opened=FALSE" %in% marker,
       "Blind freeze marker changed.")
marker_field <- function(prefix) {
  hit <- marker[startsWith(marker, prefix)]
  assert(length(hit) == 1L, paste0("Missing marker field: ", prefix))
  substring(hit, nchar(prefix) + 1L)
}
hash_table_path <- file.path(freeze, "BLIND_FILE_HASHES.csv")
assert(identical(sha256(hash_table_path), marker_field("file_hashes_sha256=")),
       "Blind file-hash table changed.")
file_hashes <- read_csv(hash_table_path)
assert(identical(names(file_hashes), c("path", "sha256")),
       "Blind file-hash schema changed.")
check_bound <- function(relative) {
  assert(!grepl("(^/|^[A-Za-z]:|\\.\\.)", relative), "Unsafe relative path.")
  hit <- file_hashes[file_hashes$path == relative, , drop = FALSE]
  assert(nrow(hit) == 1L, paste0("Unbound input: ", relative))
  path <- file.path(parent, relative)
  assert(identical(sha256(path), hit$sha256[[1L]]),
         paste0("Bound input hash mismatch: ", relative))
  path
}
manifest_path <- file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv")
assert(identical(sha256(manifest_path), marker_field("manifest_sha256=")),
       "Blind endpoint manifest changed.")
manifest <- read_csv(manifest_path)
assert(identical(manifest$branch, c("U", "A", "S")), "Branch order changed.")
observation_relative <- file.path("covariate_routeb_stage7c_20260927_v1",
                                  "data", "covrb7c_01__C1", "observation_bundle.rds")
observation_path <- file.path(parent, observation_relative)
assert(identical(sha256(observation_path),
                 "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6"),
       "Frozen observation bundle changed.")
observation <- readRDS(observation_path)
fit_path <- check_bound(manifest$fit_path[[1L]])
fit <- readRDS(fit_path)
grid <- multiFSYNC:::get_grid_objects(
  unlist(observation$time_obs, recursive = FALSE), fit$K,
  n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
n1 <- fit$n_s[[1L]]
n2 <- fit$n_s[[2L]]
C2 <- grid$C[(n1 + 1L):(n1 + n2)]
Nobs <- sum(vapply(C2, nrow, integer(1L)))
p <- fit$p
assert(p == 500L && length(C2) == 30L, "Frozen observation dimensions changed.")
rm(fit, grid)
invisible(gc())

as_root <- file.path(parent, "stage7i_as_20260929_v1")
trajectory_root <- file.path(parent, "stage7i_trajectory_readout_20260929_v1")
trajectory_path <- file.path(trajectory_root, "STAGE7I_B2_TRAJECTORY.csv")
trajectory_marker <- readLines(file.path(trajectory_root,
                                         "TRAJECTORY_COMPLETE.txt"), warn = FALSE)
assert(identical(trajectory_marker[[1L]],
                 "STAGE7I_BLIND_TRAJECTORY_READOUT_COMPLETE"),
       "Blind trajectory completion marker changed.")
trajectory_hash_line <- trajectory_marker[
  startsWith(trajectory_marker, "trajectory_sha256=")]
assert(length(trajectory_hash_line) == 1L &&
         identical(sha256(trajectory_path),
                   substring(trajectory_hash_line, nchar("trajectory_sha256=") + 1L)),
       "Blind trajectory table changed.")
trajectory <- read_csv(trajectory_path)
assert(all(trajectory$iteration %in% 50:499),
       "Blind trajectory index changed.")

g_norm <- function(psi, xi) {
  assert(identical(dim(xi), c(n2, 2L)) && ncol(psi) == 2L,
         "Study-2 B2 mean dimensions changed.")
  sqrt(sum(vapply(seq_len(n2), function(i) {
    z <- as.vector(C2[[i]] %*% psi %*% xi[i, ])
    sum(z^2)
  }, numeric(1L))))
}
stage_row <- function(branch, iteration, block, psi, xi, b) {
  G <- g_norm(psi, xi)
  B <- sqrt(sum(b^2))
  data.frame(branch = branch, iteration = iteration, block = block,
             b2_g_norm = G, b2_loading_norm = B,
             b2_mean_product_rms = B * G / sqrt(Nobs * p),
             psi_coefficient_frobenius = sqrt(sum(psi^2)),
             xi_mean_frobenius = sqrt(sum(xi^2)),
             xi_score_rms_1 = sqrt(mean(xi[, 1L]^2)),
             xi_score_rms_2 = sqrt(mean(xi[, 2L]^2)),
             stringsAsFactors = FALSE)
}
prior_log_odds <- function(state) {
  x <- state$mu_q_log_omega_b[[2L]]
  y <- state$mu_q_log_1_omega_b[[2L]]
  if (is.matrix(x) || is.matrix(y)) {
    assert(is.matrix(x) && is.matrix(y) && nrow(x) == p && nrow(y) == p,
           "Variable-level omega moments changed.")
    return(x[, 1L] - y[, 1L])
  }
  assert(length(x) >= 1L && length(y) >= 1L,
         "Factor-level omega moments changed.")
  rep(x[[1L]] - y[[1L]], p)
}
relative_close <- function(actual, expected, label) {
  difference <- abs(actual - expected)
  assert(all(is.finite(difference)) &&
           all(difference <= 1e-8 + 1e-8 * abs(expected)),
         paste0(label, " parity failed: max_abs=", max(difference)))
  max(difference)
}

branch_roots <- c(
  U = "stage7i_u_20260929_v1/u_full_v1",
  A = "stage7i_as_20260929_v1/A_full_v1",
  S = "stage7i_as_20260929_v1/S_full_v1")
all_stage <- list()
all_score <- list()
all_ppi <- list()
all_binding <- list()
for (branch in names(branch_roots)) {
  row <- manifest[manifest$branch == branch, , drop = FALSE]
  snapshot_path <- check_bound(row$snapshot_path[[1L]])
  initial_relative <- file.path("stage7i_as_20260929_v1",
                                paste0("STATE_", branch, "_SWEEP50.rds"))
  initial_path <- check_bound(initial_relative)
  driver_relative <- file.path(branch_roots[[branch]], "DRIVER_TRACE.csv")
  driver_path <- check_bound(driver_relative)
  driver <- read_csv(driver_path)
  driver <- driver[driver$study == 2L & driver$factor == 1L &
                     driver$iteration %in% 51:60, , drop = FALSE]
  assert(all(table(driver$iteration) == 5L),
         paste0(branch, " B2 driver rows changed."))
  states_full <- readRDS(snapshot_path)
  assert(length(states_full) == 449L &&
           identical(vapply(states_full[1:10], `[[`, integer(1L), "iteration"),
                     as.integer(51:60)),
         paste0(branch, " saved sweep order changed."))
  states <- c(list(readRDS(initial_path)), states_full[1:10])
  rm(states_full)
  invisible(gc())
  assert(identical(states[[1L]]$iteration, 50L),
         paste0(branch, " injected sweep changed."))
  branch_stage <- vector("list", 10L)
  branch_score <- vector("list", 10L)
  branch_ppi <- vector("list", 10L)
  for (step in 1:10) {
    old <- states[[step]]
    new <- states[[step + 1L]]
    iteration <- 50L + step
    psi_old <- old$mu_q_nu_psi[[2L]][[1L]]
    psi_new <- new$mu_q_nu_psi[[2L]][[1L]]
    xi_old <- old$mu_q_xi[[2L]][[1L]]
    xi_new <- new$mu_q_xi[[2L]][[1L]]
    b_old <- as.numeric(old$mu_q_b_specific[[2L]][, 1L])
    b_new <- as.numeric(new$mu_q_b_specific[[2L]][, 1L])
    stages <- rbind(
      stage_row(branch, iteration, "before_updates", psi_old, xi_old, b_old),
      stage_row(branch, iteration, "after_psi", psi_new, xi_old, b_old),
      stage_row(branch, iteration, "after_xi", psi_new, xi_new, b_old),
      stage_row(branch, iteration, "after_specific_loadings",
                psi_new, xi_new, b_new))
    expected <- trajectory[trajectory$branch == branch &
                             trajectory$iteration == iteration,
                           "b2_mean_product_rms"]
    assert(length(expected) == 1L,
           paste0(branch, " missing blind trajectory parity row."))
    stages$final_product_parity_max_abs <- relative_close(
      stages$b2_mean_product_rms[[4L]], expected,
      paste(branch, iteration, "B2 product"))
    branch_stage[[step]] <- stages

    c_val <- 1 / new$temperature
    score_rows <- lapply(seq_len(n2), function(i) {
      Sigma <- new$Sigma_q_xi[[2L]][[1L]][[i]]
      assert(identical(dim(Sigma), c(2L, 2L)), "B2 score covariance changed.")
      Q <- solve(Sigma)
      mu <- as.numeric(xi_new[i, ])
      h <- as.vector(Q %*% mu)
      eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
      assert(all(is.finite(eigenvalues)) && min(eigenvalues) > 0,
             "B2 score precision is not positive definite.")
      hnorm <- sqrt(sum(h^2))
      data.frame(subject = i, score_driver_norm = hnorm / c_val,
                 natural_parameter_norm = hnorm,
                 updated_score_norm = sqrt(sum(mu^2)),
                 Q_trace = sum(diag(Q)), Q_min_eigen = min(eigenvalues),
                 Q_max_eigen = max(eigenvalues),
                 Q_h_direction = if (hnorm > 0)
                   as.numeric(crossprod(h, Q %*% h)) / hnorm^2 else NA_real_,
                 stringsAsFactors = FALSE)
    })
    score_by_subject <- do.call(rbind, score_rows)
    trace_row <- driver[driver$iteration == iteration &
                          driver$block == "specific_score", , drop = FALSE]
    assert(nrow(trace_row) == 1L, "B2 score driver trace is not unique.")
    reconstructed_driver <- sqrt(sum(score_by_subject$score_driver_norm^2))
    reconstructed_mean <- sqrt(sum(score_by_subject$updated_score_norm^2))
    reconstructed_trace <- sum(score_by_subject$Q_trace)
    driver_error <- relative_close(reconstructed_driver,
                                   trace_row$linear_driver_norm,
                                   paste(branch, iteration, "score driver"))
    mean_error <- relative_close(reconstructed_mean,
                                 trace_row$posterior_mean_norm,
                                 paste(branch, iteration, "score mean"))
    precision_error <- relative_close(reconstructed_trace,
                                      trace_row$precision_trace,
                                      paste(branch, iteration, "score precision"))
    branch_score[[step]] <- data.frame(
      branch = branch, iteration = iteration, temperature = new$temperature,
      score_data_driver_norm = trace_row$data_driver_norm,
      score_linear_driver_norm = reconstructed_driver,
      score_natural_parameter_norm = reconstructed_driver * c_val,
      score_updated_mean_norm = reconstructed_mean,
      score_mean_to_natural_ratio = reconstructed_mean /
        max(reconstructed_driver * c_val, .Machine$double.eps),
      Q_trace = reconstructed_trace,
      Q_min_eigen = min(score_by_subject$Q_min_eigen),
      Q_max_eigen = max(score_by_subject$Q_max_eigen),
      Q_median_h_direction = median(score_by_subject$Q_h_direction,
                                    na.rm = TRUE),
      score_driver_parity_max_abs = driver_error,
      score_mean_parity_max_abs = mean_error,
      score_precision_parity_max_abs = precision_error)

    slab_mean <- as.numeric(new$mu_q_normal_b[[2L]][, 1L])
    slab_variance <- as.numeric(new$Sigma_q_normal_b[[2L]][, 1L])
    rho <- as.numeric(new$mu_q_gamma_b[[2L]][, 1L])
    old_rho <- as.numeric(old$mu_q_gamma_b[[2L]][, 1L])
    assert(all(is.finite(slab_variance)) && all(slab_variance > 0),
           "Invalid B2 slab variance.")
    prior <- c_val * prior_log_odds(old)
    signal <- 0.5 * slab_mean^2 / slab_variance
    variance <- 0.5 * log(slab_variance)
    reference <- rep(0.5 * (1 - c_val) * log(2 * pi), p)
    recomputed <- plogis(prior + signal + variance + reference)
    ppi_error <- relative_close(recomputed, rho,
                                paste(branch, iteration, "PPI logit"))
    loading_trace <- driver[driver$iteration == iteration &
                              driver$block == "specific_loading", , drop = FALSE]
    assert(nrow(loading_trace) == 1L &&
             !isTRUE(loading_trace$gate_active[[1L]]),
           "B2 loading trace changed or dense gate was active.")
    branch_ppi[[step]] <- data.frame(
      branch = branch, iteration = iteration,
      old_ppi_mean = mean(old_rho), new_ppi_mean = mean(rho),
      old_ppi_ge_half = sum(old_rho >= 0.5),
      new_ppi_ge_half = sum(rho >= 0.5),
      prior_logit_mean = mean(prior), signal_logit_mean = mean(signal),
      variance_logit_mean = mean(variance),
      reference_logit_mean = mean(reference),
      total_logit_mean = mean(prior + signal + variance + reference),
      signal_logit_selected_mean = if (any(rho >= 0.5))
        mean(signal[rho >= 0.5]) else NA_real_,
      signal_logit_unselected_mean = if (any(rho < 0.5))
        mean(signal[rho < 0.5]) else NA_real_,
      ppi_parity_max_abs = ppi_error,
      stringsAsFactors = FALSE)
  }
  all_stage[[branch]] <- do.call(rbind, branch_stage)
  all_score[[branch]] <- do.call(rbind, branch_score)
  all_ppi[[branch]] <- do.call(rbind, branch_ppi)
  all_binding[[branch]] <- data.frame(
    branch = branch,
    initial_state_sha256 = sha256(initial_path),
    snapshots_sha256 = row$snapshot_sha256[[1L]],
    driver_trace_sha256 = sha256(driver_path), stringsAsFactors = FALSE)
  rm(states, branch_stage, branch_score, branch_ppi, driver)
  invisible(gc())
  cat("STAGE7I_READONLY_SNAPSHOT_BRANCH_COMPLETE branch=", branch, "\n", sep = "")
}

assert(dir.create(output, mode = "0700"), "Cannot create output directory.")
write_csv(do.call(rbind, all_stage), "B2_COORDINATE_PRODUCT_51_60.csv")
write_csv(do.call(rbind, all_score), "B2_SCORE_PRECISION_51_60.csv")
write_csv(do.call(rbind, all_ppi), "B2_PPI_LOGIT_51_60.csv")
write_csv(do.call(rbind, all_binding), "SNAPSHOT_INPUT_SHA256.csv")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                      value = TRUE)[[1L]])
writeLines(c("status=COMPLETE", "branches=U,A,S", "sweeps=51:60",
             "truth_opened=FALSE", "new_fit_started=FALSE",
             "cavi_update_replayed=FALSE",
             paste0("script_sha256=", sha256(script_file)),
             paste0("coordinate_table_sha256=", sha256(file.path(
               output, "B2_COORDINATE_PRODUCT_51_60.csv")))),
           file.path(output, "SNAPSHOT_RECONSTRUCTION_COMPLETE.txt"),
           useBytes = TRUE)
cat("STAGE7I_READONLY_SNAPSHOT_RECONSTRUCTION_COMPLETE\n")
