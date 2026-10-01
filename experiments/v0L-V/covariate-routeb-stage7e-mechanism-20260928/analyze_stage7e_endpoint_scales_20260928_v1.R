#!/usr/bin/env Rscript

# Read-only scale audit for the six completed Stage 7E oracle diagnostics.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript SCRIPT STAGE7C_ROOT STAGE7E_ROOT")
stage7c_root <- normalizePath(args[[1L]], mustWork = TRUE)
stage7e_root <- normalizePath(args[[2L]], mustWork = TRUE)
if (as.character(utils::packageVersion("multiFSYNC")) != "0.3.0.9000" ||
    normalizePath(system.file(package = "multiFSYNC")) !=
      "/root/multiFSYNC-stage7c-lib-20260927-v1/multiFSYNC") {
  stop("Stage 7C isolated package binding failed.")
}
if (!file.exists(file.path(stage7e_root, "evaluation_v2",
                           "EVALUATION_COMPLETE.txt"))) {
  stop("Scientific evaluation must be complete first.")
}
output_root <- file.path(stage7e_root, "endpoint_scales_v1")
if (file.exists(output_root) || dir.exists(output_root)) {
  stop("Refusing to overwrite endpoint_scales_v1.")
}
manifest <- read.csv(file.path(stage7e_root, "FIT_MANIFEST.csv"))
stopifnot(nrow(manifest) == 6L)
trap_weights <- function(grid) {
  delta <- diff(as.numeric(grid))
  c(delta[[1L]] / 2, (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}
function_l2 <- function(value, weights) {
  sqrt(sum(sweep(as.matrix(value)^2, 1L, weights, "*")))
}
rows <- lapply(seq_len(nrow(manifest)), function(index) {
  item <- manifest[index, ]
  fit <- readRDS(file.path(stage7e_root, "fits", item$fit_id, "fit.rds"))
  truth_bundle <- readRDS(file.path(stage7c_root, "data", item$data_id,
                                    "sealed_truth", "truth_bundle.rds"))
  stopifnot(identical(truth_bundle$data_id, item$data_id))
  truth <- truth_bundle$data$true_params
  grid <- as.numeric(truth$t_grid_dense)
  weights <- trap_weights(grid)
  design <- cbind(1, grid, multiFSYNC:::ZOSull(
    grid, range.x = c(0, 1), intKnots = truth$int_knots))
  stopifnot(ncol(design) == nrow(fit$mu_q_nu_psi[[2L]][[1L]]))
  estimated_function <- design %*% fit$mu_q_nu_psi[[2L]][[1L]]
  true_function <- truth$kappa_dense_list[[2L]][[1L]]
  estimated_score <- fit$mu_q_xi[[2L]][[1L]]
  true_score <- truth$chi_true[[2L]][[1L]]
  estimated_loading <- fit$mu_q_b_specific[[2L]][, 1L]
  true_loading <- truth$b_true[[2L]][, 1L]
  true_function_norm <- function_l2(true_function, weights)
  true_score_rms <- sqrt(mean(true_score^2))
  true_loading_norm <- sqrt(sum(true_loading^2))
  data.frame(
    fit_id = item$fit_id, data_id = item$data_id, mode = item$mode,
    raw_loading_norm = sqrt(sum(estimated_loading^2)),
    true_loading_norm = true_loading_norm,
    raw_loading_norm_ratio = sqrt(sum(estimated_loading^2)) /
      true_loading_norm,
    estimated_time_function_l2 = function_l2(estimated_function, weights),
    true_time_function_l2 = true_function_norm,
    time_function_l2_ratio =
      function_l2(estimated_function, weights) / true_function_norm,
    estimated_score_rms = sqrt(mean(estimated_score^2)),
    true_score_rms = true_score_rms,
    score_rms_ratio = sqrt(mean(estimated_score^2)) / true_score_rms,
    scale_product_ratio =
      sqrt(sum(estimated_loading^2)) / true_loading_norm *
      function_l2(estimated_function, weights) / true_function_norm *
      sqrt(mean(estimated_score^2)) / true_score_rms,
    stringsAsFactors = FALSE)
})
dir.create(output_root)
write.csv(do.call(rbind, rows),
          file.path(output_root, "ORACLE_B2_ENDPOINT_SCALES.csv"),
          row.names = FALSE)
writeLines(c("STAGE7E_ENDPOINT_SCALES_COMPLETE", "fits=6",
             "truth_used_for_scientific_diagnosis=TRUE",
             "new_fit_started=FALSE"),
           file.path(output_root, "ENDPOINT_SCALES_COMPLETE.txt"))
cat("STAGE7E_ENDPOINT_SCALES_COMPLETE fits=6\n")
