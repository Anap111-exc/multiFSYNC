#!/usr/bin/env Rscript

# Post-unseal, read-only normalization check for observed/dense signal NRMSE.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, dir.exists(args[[1L]]),
          !file.exists(args[[2L]]))
root <- args[[1L]]
manifest <- utils::read.csv(file.path(root, "DATA_MANIFEST.csv"))
stopifnot(nrow(manifest) == 12L)
rows <- lapply(manifest$data_id, function(data_id) {
  sealed <- readRDS(file.path(root, "data", data_id,
                             "sealed_truth", "truth_bundle.rds"))
  stopifnot(identical(sealed$bundle_class, "stage7a_sealed_truth"),
            identical(sealed$data_id, data_id))
  data <- sealed$data
  truth <- data$true_params
  grid <- truth$t_grid_dense
  beta <- truth$beta_true_grid
  obs <- unlist(truth$signal_true_values, recursive = TRUE,
                use.names = FALSE)
  dense_sq <- 0
  dense_n <- 0L
  for (s in seq_len(truth$S)) {
    mean_matrix <- vapply(seq_len(truth$p), function(j) {
      phase <- truth$phase_mu_mat[s, j]
      0.6 * (sin(2 * pi * grid + phase) +
               0.5 * sin(4 * pi * grid + 2 * phase))
    }, numeric(length(grid)))
    for (i in seq_len(truth$n_s[[s]])) {
      shared <- truth$theta_dense_list[[1L]] %*%
        truth$eta_true[[s]][[1L]][i, ]
      specific <- truth$kappa_dense_list[[s]][[1L]] %*%
        truth$chi_true[[s]][[1L]][i, ]
      signal <- mean_matrix +
        as.numeric(shared) %o% truth$a_true[, 1L] +
        as.numeric(specific) %o% truth$b_true[[s]][, 1L] +
        data$Z[[s]][i, 1L] * beta
      dense_sq <- dense_sq + sum(signal^2)
      dense_n <- dense_n + length(signal)
    }
  }
  data.frame(data_id = data_id,
             observed_truth_rms = sqrt(mean(obs^2)),
             dense_truth_rms = sqrt(dense_sq / dense_n),
             observed_values = length(obs), dense_values = dense_n)
})
result <- do.call(rbind, rows)
result$dense_to_observed_truth_rms_ratio <-
  result$dense_truth_rms / result$observed_truth_rms
utils::write.csv(result, args[[2L]], row.names = FALSE)
cat("TRUTH_ENERGY_DIAGNOSTIC_COMPLETE rows=", nrow(result), "\n", sep = "")
