#!/usr/bin/env Rscript

# Read-only Stage 7C endpoint audit. No truth bundle is opened.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analyze_stage7e_raw_post_20260928_v1.R STAGE7C_ROOT OUTPUT_DIR")
}
stage7c_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
if (dir.exists(output_dir) || file.exists(output_dir)) {
  stop("Refusing to overwrite an existing Stage 7E output directory.")
}
dir.create(output_dir, recursive = TRUE)

fit_ids <- c(
  "covrb7c_01__C0__start06", "covrb7c_01__C1__start06",
  "covrb7c_01__C1__start11", "covrb7c_02__C0__start11",
  "covrb7c_02__C1__start11", "covrb7c_02__C1__start08",
  "covrb7c_03__C0__start07", "covrb7c_03__C1__start06"
)

safe_rms <- function(x) sqrt(mean(x^2))
safe_relative <- function(x, reference) {
  safe_rms(x) / max(safe_rms(reference), .Machine$double.eps)
}
unit_cosine <- function(x, y) {
  norm <- sqrt(sum(x^2) * sum(y^2))
  if (norm <= .Machine$double.eps) return(NA_real_)
  abs(sum(x * y) / norm)
}

factor_rows <- list()
reconstruction_rows <- list()
for (fit_id in fit_ids) {
  fit_path <- file.path(stage7c_root, "fits", fit_id, "fit.rds")
  data_id <- sub("__start[0-9]+$", "", fit_id)
  observation_path <- file.path(
    stage7c_root, "data", data_id, "observation_bundle.rds"
  )
  fit <- readRDS(fit_path)
  observation <- readRDS(observation_path)
  stopifnot(identical(observation$data_id, data_id))
  stopifnot(length(fit$mu_q_zeta) == 2L,
            length(observation$Z) == 2L,
            nrow(fit$C_g) == length(fit$time_g))
  C_g <- fit$C_g
  grid_n <- nrow(C_g)
  p <- length(fit$response_scale)
  study_start <- 1L

  for (study in seq_len(2L)) {
    Z <- observation$Z[[study]]
    n <- nrow(Z)
    rows <- study_start:(study_start + n - 1L)
    study_start <- study_start + n
    stopifnot(ncol(Z) == fit$d, length(fit$list_mu_hat[[study]]) == p)
    raw_factor <- matrix(0, n * grid_n, p)
    post_factor <- matrix(0, n * grid_n, p)

    add_factor <- function(block, factor, raw_score, raw_coef,
                           raw_loading, post_score, post_phi,
                           post_loading, raw_gamma, post_gamma,
                           effective_M) {
      raw_path <- raw_score %*% t(C_g %*% raw_coef)
      post_path <- post_score %*% t(post_phi)
      raw_contribution <- tcrossprod(
        as.vector(t(raw_path)), raw_loading
      )
      post_contribution <- tcrossprod(
        as.vector(t(post_path)), post_loading
      )
      raw_factor <<- raw_factor + raw_contribution
      post_factor <<- post_factor + post_contribution
      factor_rows[[length(factor_rows) + 1L]] <<- data.frame(
        fit_id = fit_id, data_id = data_id, study = study,
        block = block, factor = factor,
        raw_ppi_mean = mean(raw_gamma),
        raw_loading_support_ge_half = sum(raw_gamma >= 0.5),
        post_loading_support_ge_half = sum(post_gamma >= 0.5),
        max_abs_ppi_difference = max(abs(raw_gamma - post_gamma)),
        loading_abs_cosine = unit_cosine(raw_loading, post_loading),
        loading_norm_ratio = sqrt(sum(post_loading^2)) /
          max(sqrt(sum(raw_loading^2)), .Machine$double.eps),
        raw_contribution_rms = safe_rms(raw_contribution),
        post_contribution_rms = safe_rms(post_contribution),
        difference_rms = safe_rms(raw_contribution - post_contribution),
        difference_relative_raw = safe_relative(
          raw_contribution - post_contribution, raw_contribution
        ),
        effective_M = effective_M
      )
    }

    for (factor in seq_along(fit$mu_q_nu_phi)) {
      add_factor(
        "shared", factor, fit$mu_q_zeta[[study]][[factor]],
        fit$mu_q_nu_phi[[factor]],
        fit$mu_q_a[, factor] * fit$response_scale,
        fit$list_Zeta_hat[[factor]][rows, , drop = FALSE],
        fit$list_Phi_hat[[factor]], fit$mu_q_a_hat[, factor],
        fit$mu_q_gamma_a[, factor], fit$mu_q_gamma_a_hat[, factor],
        fit$list_effective_M[[factor]]
      )
    }
    for (factor in seq_along(fit$mu_q_nu_psi[[study]])) {
      add_factor(
        "specific", factor, fit$mu_q_xi[[study]][[factor]],
        fit$mu_q_nu_psi[[study]][[factor]],
        fit$mu_q_b_specific[[study]][, factor] * fit$response_scale,
        fit$list_Zeta_hat_spec[[study]][[factor]],
        fit$list_Phi_hat_spec[[study]][[factor]],
        fit$mu_q_b_specific_hat[[study]][, factor],
        fit$mu_q_gamma_b[[study]][, factor],
        fit$mu_q_gamma_b_hat[[study]][, factor],
        fit$list_effective_M_spec[[study]][[factor]]
      )
    }

    raw_base <- matrix(0, n * grid_n, p)
    post_base <- matrix(0, n * grid_n, p)
    mean_max_diff <- 0
    beta_max_diff <- 0
    for (variable in seq_len(p)) {
      raw_mean <- fit$response_center[[variable]] +
        fit$response_scale[[variable]] *
          as.vector(C_g %*% fit$mu_q_nu_mu[[study]][[variable]])
      post_mean <- fit$list_mu_hat[[study]][[variable]]
      mean_max_diff <- max(mean_max_diff, max(abs(raw_mean - post_mean)))
      raw_base[, variable] <- rep(raw_mean, times = n)
      post_base[, variable] <- rep(post_mean, times = n)
      for (covariate in seq_len(fit$d)) {
        raw_beta <- fit$response_scale[[variable]] *
          as.vector(C_g %*% fit$mu_q_nu_beta[[variable]][[covariate]])
        post_beta <- fit$list_beta_hat[[variable]][[covariate]]
        beta_max_diff <- max(beta_max_diff, max(abs(raw_beta - post_beta)))
        multiplier <- rep(Z[, covariate], each = grid_n)
        raw_base[, variable] <- raw_base[, variable] +
          multiplier * rep(raw_beta, times = n)
        post_base[, variable] <- post_base[, variable] +
          multiplier * rep(post_beta, times = n)
      }
    }
    raw_full <- raw_base + raw_factor
    post_full <- post_base + post_factor
    reconstruction_rows[[length(reconstruction_rows) + 1L]] <- data.frame(
      fit_id = fit_id, data_id = data_id, study = study,
      n_subjects = n, grid_n = grid_n, p = p,
      mean_max_abs_difference = mean_max_diff,
      beta_max_abs_difference = beta_max_diff,
      raw_factor_rms = safe_rms(raw_factor),
      post_factor_rms = safe_rms(post_factor),
      factor_difference_rms = safe_rms(raw_factor - post_factor),
      factor_difference_relative_raw = safe_relative(
        raw_factor - post_factor, raw_factor
      ),
      raw_full_rms = safe_rms(raw_full),
      post_full_rms = safe_rms(post_full),
      full_difference_rms = safe_rms(raw_full - post_full),
      full_difference_relative_raw = safe_relative(
        raw_full - post_full, raw_full
      )
    )
  }
}

factor_table <- do.call(rbind, factor_rows)
reconstruction_table <- do.call(rbind, reconstruction_rows)
stopifnot(nrow(factor_table) == 32L,
          nrow(reconstruction_table) == 16L,
          all(is.finite(reconstruction_table$full_difference_relative_raw)))
write.csv(factor_table, file.path(output_dir, "RAW_POST_FACTOR.csv"),
          row.names = FALSE)
write.csv(reconstruction_table,
          file.path(output_dir, "RAW_POST_RECONSTRUCTION.csv"),
          row.names = FALSE)
writeLines(c(
  "STAGE7E_RAW_POST_READONLY_COMPLETE",
  paste0("fits=", length(fit_ids)),
  paste0("factor_rows=", nrow(factor_table)),
  paste0("reconstruction_rows=", nrow(reconstruction_table)),
  "truth_read=FALSE", "new_fit=FALSE"
), file.path(output_dir, "READONLY_COMPLETE.txt"))
