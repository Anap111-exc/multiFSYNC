#!/usr/bin/env Rscript

# Stage 7H: read-only reconstruction of frozen C1 paths. No truth or fit call.
options(warn = 2, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop(
  "Usage: Rscript run_stage7h_blind_20260929_v1.R ",
  "STAGE7C_ROOT STAGE7G_ROOT PREFLIGHT_ROOT SMOKE_ROOT OUTPUT_ROOT MODE",
  call. = FALSE)
stage7c <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g <- normalizePath(args[[2L]], mustWork = TRUE)
preflight <- normalizePath(args[[3L]], mustWork = TRUE)
smoke <- normalizePath(args[[4L]], mustWork = TRUE)
out <- args[[5L]]
mode <- args[[6L]]
if (!mode %in% c("one", "all")) stop("Mode must be one or all.")
if (dir.exists(out)) stop("Refusing existing output directory: ", out)
if (!file.exists(file.path(preflight, "PREFLIGHT_PASS.txt")) ||
    !file.exists(file.path(smoke, "SMOKE_PASS.txt")))
  stop("Preflight and smoke must pass first.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Cannot create output directory.")

.libPaths(c("/root/multiFSYNC-stage7g-lib-20260928-v1", .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
if (normalizePath(system.file(package = "multiFSYNC")) !=
    "/root/multiFSYNC-stage7g-lib-20260928-v1/multiFSYNC")
  stop("Wrong frozen Stage 7G package.")

sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(line) != 1L) stop("sha256sum failed: ", path)
  strsplit(line, "[[:space:]]+")[[1L]][[1L]]
}
max_abs <- function(a, b) {
  if (length(a) != length(b)) return(Inf)
  max(abs(as.numeric(a) - as.numeric(b)))
}
norm2 <- function(x) sqrt(sum(x^2))
tol <- function(ref) 1e-7 + 1e-7 * max(1, max(abs(ref)))
check <- function(name, found, reference, fit_id, iteration,
                  records) {
  err <- max_abs(found, reference)
  allowed <- tol(reference)
  rows <- data.frame(fit_id = fit_id, iteration = iteration,
                     field = name, max_abs_difference = err,
                     tolerance = allowed,
                     pass = is.finite(err) && err <= allowed,
                     stringsAsFactors = FALSE)
  records[[length(records) + 1L]] <- rows
  if (!rows$pass) stop("Parity failure ", fit_id, " t=", iteration,
                       " ", name, " err=", err, " tol=", allowed)
  records
}
gram_eigen <- function(G, fit_id, iteration, stage, domain) {
  value <- eigen((G + t(G)) / 2, symmetric = TRUE)
  singular <- sqrt(pmax(value$values, 0))
  U <- value$vectors[, 1:5, drop = FALSE]
  gaps <- singular[1:5] / pmax(singular[2:6], .Machine$double.eps)
  names(singular) <- NULL
  row <- data.frame(fit_id = fit_id, iteration = iteration,
                    stage = stage, domain = domain,
                    sigma1 = singular[[1L]], sigma2 = singular[[2L]],
                    sigma3 = singular[[3L]], sigma4 = singular[[4L]],
                    sigma5 = singular[[5L]],
                    gap1 = gaps[[1L]], gap2 = gaps[[2L]],
                    gap3 = gaps[[3L]], gap4 = gaps[[4L]],
                    gap5 = gaps[[5L]], stringsAsFactors = FALSE)
  list(row = row, vectors = U)
}
full_time_second <- function(C, coeff, coeff_cov, score_mean,
                             score_cov) {
  F <- C %*% coeff
  E_score <- score_cov + tcrossprod(score_mean)
  out <- rowSums((F %*% E_score) * F)
  for (m in seq_len(ncol(F))) {
    CF <- C %*% coeff_cov[[m]]
    out <- out + E_score[m, m] * rowSums(CF * C)
  }
  if (any(!is.finite(out)) || min(out) < -1e-6)
    stop("Invalid full time second moment.")
  pmax(out, 0)
}

manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
if (nrow(manifest) != 3L || nrow(binding) != 3L ||
    !identical(as.character(manifest$fit_id),
               as.character(binding$fit_id)))
  stop("Input manifest changed.")
selected <- if (mode == "one") 3L else 1:3
checks <- drivers <- scores <- energy <- geometry <- list()
eigenvectors <- list()
row_c <- row_d <- row_s <- row_e <- row_g <- 0L
fixed_times <- c(20L, 25L, 40L, 50L, 55L, 60L)

for (k in selected) {
  fit_id <- manifest$fit_id[[k]]
  data_id <- manifest$data_id[[k]]
  obs_path <- file.path(stage7c, "data", data_id,
                        "observation_bundle.rds")
  fit_path <- file.path(stage7g, "fits", fit_id, "fit.rds")
  snap_path <- file.path(stage7g, "fits", fit_id,
                         "RAW_STATE_SNAPSHOTS.rds")
  trace_path <- file.path(stage7g, "fits", fit_id, "DRIVER_TRACE.csv")
  if (!identical(sha256(obs_path), manifest$observation_sha256[[k]]) ||
      !identical(sha256(fit_path), binding$fit_sha256[[k]]) ||
      !identical(sha256(snap_path), binding$snapshot_sha256[[k]]))
    stop("Input hash mismatch: ", fit_id)
  obs <- readRDS(obs_path)
  fit <- readRDS(fit_path)
  trace <- utils::read.csv(trace_path, stringsAsFactors = FALSE)
  snapshots <- readRDS(snap_path)
  if (length(snapshots) != 500L ||
      !identical(vapply(snapshots, `[[`, integer(1L), "iteration"),
                 as.integer(0:499))) stop("Snapshot schedule changed.")
  basis <- multiFSYNC:::get_grid_objects(
    unlist(obs$time_obs, recursive = FALSE), fit$K,
    n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
  n1 <- fit$n_s[[1L]]
  n <- fit$n_s[[2L]]
  C <- basis$C[(n1 + 1L):(n1 + n)]
  qrC <- lapply(C, qr, tol = 1e-8)
  design_times <- vapply(C, nrow, integer(1L))
  design_ranks <- vapply(qrC, `[[`, integer(1L), "rank")
  design_conditions <- vapply(C, function(x)
    kappa(x, exact = FALSE), numeric(1L))
  Y <- lapply(obs$Y[[2L]], function(x) do.call(cbind, x))
  Z <- obs$Z[[2L]]
  p <- fit$p
  M <- fit$M_s[[2L]][[1L]]
  if (n != length(C) || n != length(Y) || !is.matrix(Z) ||
      nrow(Z) != n || p != 500L || M != 2L)
    stop("Frozen design shape changed: ", fit_id)
  time_points <- if (k == 3L || mode == "one") fixed_times else
    c(fixed_times, 67L)

  for (iteration in time_points) {
    before <- snapshots[[iteration]]
    state <- snapshots[[iteration + 1L]]
    c_val <- 1 / state$temperature
    tau <- before$mu_q_recip_sigsq_eps[2L, ]
    prior_b <- before$mu_q_b_specific[[2L]][, 1L]
    prior_b2 <- before$term_b_specific[[2L]][, 1L]
    prior_a <- before$mu_q_a[, 1L]
    mean_coeff <- do.call(cbind, state$mu_q_nu_mu[[2L]])
    beta_coeff <- lapply(seq_len(fit$d), function(r)
      do.call(cbind, lapply(state$mu_q_nu_beta, `[[`, r)))
    phi <- state$mu_q_nu_phi[[1L]]
    phi_cov <- state$Sigma_q_nu_phi[[1L]]
    psi <- state$mu_q_nu_psi[[2L]][[1L]]
    psi_cov <- state$Sigma_q_nu_psi[[2L]][[1L]]
    zeta <- state$mu_q_zeta[[2L]][[1L]]
    zeta_cov <- state$Sigma_q_zeta[[2L]][[1L]]
    xi <- state$mu_q_xi[[2L]][[1L]]
    xi_cov <- state$Sigma_q_xi[[2L]][[1L]]
    a <- state$mu_q_a[, 1L]
    b <- state$mu_q_b_specific[[2L]][, 1L]
    if (any(!is.finite(tau)) || any(tau <= 0) ||
        !all(dim(mean_coeff) == c(fit$K_total, p)) ||
        !all(dim(psi) == c(fit$K_total, M)))
      stop("State contract changed: ", fit_id, " t=", iteration)

    # Exact loading-coordinate sufficient statistics and four U terms.
    h_data <- h_mean <- h_beta <- h_shared <- numeric(p)
    U_mean <- U_score_cov <- U_function_cov <- U_joint_cov <- 0
    load_envelope <- 0
    # Exact score-coordinate reconstruction uses the previous loadings.
    score_data_norm_sq <- score_linear_norm_sq <- 0
    score_mean_norm_sq <- score_precision_norm_sq <- 0
    score_envelope <- score_precision_trace <- 0
    score_cov_error <- score_mean_error <- 0
    score_rows_t <- vector("list", n)
    sum_b <- sum(tau * prior_b2)
    sum_b_mean <- sum(tau * prior_b^2)
    sum_b_var <- sum(tau * pmax(prior_b2 - prior_b^2, 0))
    stages <- c("observed", "minus_mean", "minus_beta",
                "minus_shared")
    domains <- c("raw", "basis_projected")
    G <- lapply(stages, function(x) list(raw = matrix(0, p, p),
                                          basis_projected = matrix(0, p, p)))
    names(G) <- stages
    splitG <- lapply(c("odd", "even"), function(x)
      list(raw = matrix(0, p, p),
           basis_projected = matrix(0, p, p)))
    names(splitG) <- c("odd", "even")
    total_raw_resid <- total_basis_resid <- total_func_resid <- 0
    shared_sq <- specific_sq <- block_inner <- 0
    shared_second <- specific_second <- 0
    total_cells <- 0
    max_residual_crossprod_error <- 0
    function_rank_min <- Inf
    function_condition_max <- 0

    for (i in seq_len(n)) {
      C_i <- C[[i]]
      cpC_i <- fit$list_cp_C[[2L]][[i]]
      Yi <- Y[[i]]
      mean_i <- C_i %*% mean_coeff
      beta_i <- matrix(0, nrow(C_i), p)
      for (r in seq_len(fit$d))
        beta_i <- beta_i + Z[i, r] * (C_i %*% beta_coeff[[r]])
      f_shared_i <- as.vector(C_i %*% phi %*% zeta[i, ])
      shared_i <- tcrossprod(f_shared_i, a)
      residuals <- list(
        observed = Yi,
        minus_mean = Yi - mean_i,
        minus_beta = Yi - mean_i - beta_i,
        minus_shared = Yi - mean_i - beta_i - shared_i)
      R <- residuals$minus_shared
      ct_res <- crossprod(C_i, R)
      ct_suff <- fit$list_cp_C_Y[[2L]][[i]] - cpC_i %*% mean_coeff
      for (r in seq_len(fit$d))
        ct_suff <- ct_suff - Z[i, r] * (cpC_i %*% beta_coeff[[r]])
      ct_suff <- ct_suff -
        tcrossprod(as.vector(cpC_i %*% phi %*% zeta[i, ]), a)
      max_residual_crossprod_error <- max(
        max_residual_crossprod_error, max(abs(ct_res - ct_suff)))

      Gi <- C_i %*% psi
      gi <- as.vector(Gi %*% xi[i, ])
      time_products <- lapply(residuals, function(Rs)
        as.vector(crossprod(gi, Rs)))
      h_data <- h_data + tau * time_products$observed
      h_mean <- h_mean - tau *
        (time_products$observed - time_products$minus_mean)
      h_beta <- h_beta - tau *
        (time_products$minus_mean - time_products$minus_beta)
      h_shared <- h_shared - tau *
        (time_products$minus_beta - time_products$minus_shared)
      load_envelope <- load_envelope +
        sum(abs(tau * time_products$minus_shared))

      H_mean <- crossprod(psi, cpC_i %*% psi)
      H_cov <- diag(vapply(seq_len(M), function(m)
        sum(diag(cpC_i %*% psi_cov[[m]])), numeric(1L)), M)
      E_xi <- xi_cov[[i]] + tcrossprod(xi[i, ])
      U_mean <- U_mean + sum(diag(H_mean %*% tcrossprod(xi[i, ])))
      U_score_cov <- U_score_cov + sum(diag(H_mean %*% xi_cov[[i]]))
      U_function_cov <- U_function_cov +
        sum(diag(H_cov %*% tcrossprod(xi[i, ])))
      U_joint_cov <- U_joint_cov + sum(diag(H_cov %*% xi_cov[[i]]))

      # Score update: shared loading is still from t-1 at block 6.
      wb <- tau * prior_b
      wb[abs(wb) < 1e-15] <- 0
      score_data_ct <- crossprod(C_i, residuals$minus_mean)
      score_data_i <- as.vector(score_data_ct %*% wb)
      score_beta_i <- as.vector(crossprod(C_i, beta_i) %*% wb)
      wa <- tau * prior_b * prior_a
      wa[abs(wa) < 1e-15] <- 0
      score_shared_i <- as.vector(crossprod(C_i, f_shared_i)) * sum(wa)
      score_ct <- score_data_i - score_beta_i - score_shared_i
      score_h <- as.vector(crossprod(psi, score_ct))
      score_prec <- c_val * (sum_b * (H_mean + H_cov) + diag(M))
      score_cov_i <- solve(score_prec)
      score_mean_i <- as.vector(c_val * score_cov_i %*% score_h)
      score_cov_error <- max(score_cov_error,
                             max_abs(score_cov_i, xi_cov[[i]]))
      score_mean_error <- max(score_mean_error,
                              max_abs(score_mean_i, xi[i, ]))
      score_data_norm_sq <- score_data_norm_sq + sum(score_data_i^2)
      score_envelope <- score_envelope +
        sum(abs(wb) * sqrt(colSums(score_data_ct^2)))
      score_linear_norm_sq <- score_linear_norm_sq + sum(score_h^2)
      score_mean_norm_sq <- score_mean_norm_sq + sum(score_mean_i^2)
      score_precision_norm_sq <- score_precision_norm_sq +
        sum(score_prec^2)
      score_precision_trace <- score_precision_trace + sum(diag(score_prec))
      score_rows_t[[i]] <- data.frame(
        fit_id = fit_id, iteration = iteration, subject = i,
        driver_data_norm = norm2(as.vector(crossprod(psi, score_data_i))),
        driver_beta_norm = norm2(as.vector(crossprod(psi, score_beta_i))),
        driver_shared_norm =
          norm2(as.vector(crossprod(psi, score_shared_i))),
        driver_net_norm = norm2(score_h),
        driver_cancellation_ratio = norm2(score_h) /
          max(norm2(as.vector(crossprod(psi, score_data_i))) +
                norm2(as.vector(crossprod(psi, score_beta_i))) +
                norm2(as.vector(crossprod(psi, score_shared_i))),
              .Machine$double.eps),
        load_mean_weight = sum_b_mean,
        load_variance_weight = sum_b_var,
        precision_prior_trace = c_val * M,
        precision_mean_mean_trace =
          c_val * sum_b_mean * sum(diag(H_mean)),
        precision_loadvar_mean_trace =
          c_val * sum_b_var * sum(diag(H_mean)),
        precision_mean_funcvar_trace =
          c_val * sum_b_mean * sum(diag(H_cov)),
        precision_jointvar_trace =
          c_val * sum_b_var * sum(diag(H_cov)),
        stringsAsFactors = FALSE)

      # Raw observed-point geometry and geometry after basis projection.
      for (name in stages) {
        Rs <- residuals[[name]]
        Proj <- qr.fitted(qrC[[i]], Rs)
        G[[name]]$raw <- G[[name]]$raw + crossprod(Rs)
        G[[name]]$basis_projected <-
          G[[name]]$basis_projected + crossprod(Proj)
        if (name == "minus_shared") {
          group <- if (i %% 2L == 1L) "odd" else "even"
          splitG[[group]]$raw <- splitG[[group]]$raw + crossprod(Rs)
          splitG[[group]]$basis_projected <-
            splitG[[group]]$basis_projected + crossprod(Proj)
          total_raw_resid <- total_raw_resid + sum(Rs^2)
          total_basis_resid <- total_basis_resid + sum(Proj^2)
        }
      }
      qrGi <- qr(Gi, tol = 1e-8)
      function_rank_min <- min(function_rank_min, qrGi$rank)
      function_condition_max <- max(function_condition_max,
                                     kappa(Gi, exact = FALSE))
      if (qrGi$rank > 0L)
        total_func_resid <- total_func_resid +
          sum(qr.fitted(qrGi, R)^2)

      f_specific_i <- as.vector(Gi %*% xi[i, ])
      F_specific <- tcrossprod(f_specific_i, b)
      shared_sq <- shared_sq + sum(shared_i^2)
      specific_sq <- specific_sq + sum(F_specific^2)
      block_inner <- block_inner + sum(shared_i * F_specific)
      shared_second <- shared_second +
        sum(state$term_a[, 1L]) * sum(full_time_second(
          C_i, phi, phi_cov, zeta[i, ], zeta_cov[[i]]))
      specific_second <- specific_second +
        sum(state$term_b_specific[[2L]][, 1L]) *
        sum(full_time_second(C_i, psi, psi_cov, xi[i, ], xi_cov[[i]]))
      total_cells <- total_cells + length(Yi)
    }

    U <- U_mean + U_score_cov + U_function_cov + U_joint_cov
    h <- h_data + h_mean + h_beta + h_shared
    prec <- c_val * (tau * U + 1)
    v <- 1 / prec
    m <- c_val * v * h
    log_omega <- before$mu_q_log_omega_b[[2L]]
    log_1_omega <- before$mu_q_log_1_omega_b[[2L]]
    prior_log_odds <- if (is.matrix(log_omega))
      log_omega[, 1L] - log_1_omega[, 1L] else
      rep(log_omega[[1L]] - log_1_omega[[1L]], p)
    logit_prior <- c_val * prior_log_odds
    logit_signal <- 0.5 * m^2 / v
    logit_variance <- 0.5 * log(v)
    logit_measure <- 0.5 * (1 - c_val) * log(2 * pi)
    rho <- stats::plogis(logit_prior + logit_signal +
                           logit_variance + logit_measure)
    b_rebuilt <- rho * m
    second_rebuilt <- rho * (v + m^2)
    trace_load <- trace[trace$iteration == iteration &
                          trace$phase == "main" &
                          trace$block == "specific_loading" &
                          trace$study == 2L & trace$factor == 1L,
                        , drop = FALSE]
    trace_score <- trace[trace$iteration == iteration &
                           trace$phase == "main" &
                           trace$block == "specific_score" &
                           trace$study == 2L & trace$factor == 1L,
                         , drop = FALSE]
    trace_scheduled <- nrow(trace_load) == 1L &&
      nrow(trace_score) == 1L
    if (!trace_scheduled &&
        !(iteration == 67L && k %in% 1:2 &&
          nrow(trace_load) == 0L && nrow(trace_score) == 0L))
      stop("Unexpected trace availability: ", fit_id,
           " t=", iteration)
    load_metrics <- c(
      data_driver_norm = norm2(h),
      data_envelope = load_envelope,
      data_coherence = norm2(h) /
        max(load_envelope, .Machine$double.eps),
      linear_driver_norm = norm2(h),
      precision_trace = sum(prec),
      precision_frobenius = norm2(prec),
      prior_precision_trace = c_val * p,
      driver_precision_ratio = norm2(h) /
        max(norm2(prec), .Machine$double.eps),
      posterior_mean_norm = norm2(b_rebuilt),
      slab_mean_norm = norm2(m),
      slab_variance_mean = mean(v),
      unrestricted_ppi_min = min(rho),
      unrestricted_ppi_mean = mean(rho),
      unrestricted_ppi_max = max(rho),
      ppi_min = min(rho), ppi_mean = mean(rho), ppi_max = max(rho),
      log_prior_odds_mean = mean(prior_log_odds))
    score_metrics <- c(
      data_driver_norm = sqrt(score_data_norm_sq),
      data_envelope = score_envelope,
      data_coherence = sqrt(score_data_norm_sq) /
        max(score_envelope, .Machine$double.eps),
      linear_driver_norm = sqrt(score_linear_norm_sq),
      precision_trace = score_precision_trace,
      precision_frobenius = sqrt(score_precision_norm_sq),
      prior_precision_trace = c_val * n * M,
      driver_precision_ratio = sqrt(score_linear_norm_sq) /
        max(sqrt(score_precision_norm_sq), .Machine$double.eps),
      posterior_mean_norm = sqrt(score_mean_norm_sq))
    check_rows <- list()
    check_rows <- check("observation_residual_crossprod",
                        max_residual_crossprod_error, 0,
                        fit_id, iteration, check_rows)
    check_rows <- check("loading_slab_mean", m,
                        state$mu_q_normal_b[[2L]][, 1L],
                        fit_id, iteration, check_rows)
    check_rows <- check("loading_slab_variance", v,
                        state$Sigma_q_normal_b[[2L]][, 1L],
                        fit_id, iteration, check_rows)
    check_rows <- check("loading_ppi", rho,
                        state$mu_q_gamma_b[[2L]][, 1L],
                        fit_id, iteration, check_rows)
    check_rows <- check("loading_effective_mean", b_rebuilt, b,
                        fit_id, iteration, check_rows)
    check_rows <- check("loading_second", second_rebuilt,
                        state$term_b_specific[[2L]][, 1L],
                        fit_id, iteration, check_rows)
    check_rows <- check("score_covariance", score_cov_error, 0,
                        fit_id, iteration, check_rows)
    check_rows <- check("score_mean", score_mean_error, 0,
                        fit_id, iteration, check_rows)
    if (trace_scheduled) {
      for (field in names(load_metrics))
        check_rows <- check(paste0("loading_trace_", field),
                            load_metrics[[field]],
                            as.numeric(trace_load[[field]][[1L]]),
                            fit_id, iteration, check_rows)
      for (field in names(score_metrics))
        check_rows <- check(paste0("score_trace_", field),
                            score_metrics[[field]],
                            as.numeric(trace_score[[field]][[1L]]),
                            fit_id, iteration, check_rows)
    }
    row_c <- row_c + 1L
    checks[[row_c]] <- do.call(rbind, check_rows)

    row_d <- row_d + 1L
    drivers[[row_d]] <- data.frame(
      fit_id = fit_id, iteration = iteration, coordinate = seq_len(p),
      tau = tau, h_data = h_data, h_mean = h_mean,
      h_beta = h_beta, h_shared = h_shared, h_total = h,
      cancellation_ratio = abs(h) /
        pmax(abs(h_data) + abs(h_mean) + abs(h_beta) + abs(h_shared),
             .Machine$double.eps),
      U_mean = U_mean, U_score_cov = U_score_cov,
      U_function_cov = U_function_cov, U_joint_cov = U_joint_cov,
      precision = prec, slab_variance = v, slab_mean = m,
      prior_log_odds = prior_log_odds,
      logit_prior = logit_prior, logit_signal = logit_signal,
      logit_variance = logit_variance,
      logit_measure = logit_measure,
      ppi = rho, effective_loading = b_rebuilt,
      loading_second = second_rebuilt,
      stringsAsFactors = FALSE)
    row_s <- row_s + 1L
    scores[[row_s]] <- do.call(rbind, score_rows_t)
    row_e <- row_e + 1L
    energy[[row_e]] <- data.frame(
      fit_id = fit_id, iteration = iteration,
      shared_mean_rms = sqrt(shared_sq / total_cells),
      B2_mean_rms = sqrt(specific_sq / total_cells),
      shared_B2_mean_cosine = block_inner /
        max(sqrt(shared_sq * specific_sq), .Machine$double.eps),
      shared_B2_mean_inner_per_cell = block_inner / total_cells,
      shared_full_second_rms = sqrt(shared_second / total_cells),
      B2_full_second_rms = sqrt(specific_second / total_cells),
      fitted_shared_B2_loading_abs_cosine =
        abs(sum(a * b)) / max(norm2(a) * norm2(b),
                            .Machine$double.eps),
      U = U, U_mean = U_mean, U_score_cov = U_score_cov,
      U_function_cov = U_function_cov, U_joint_cov = U_joint_cov,
      basis_capture_ratio = total_basis_resid /
        max(total_raw_resid, .Machine$double.eps),
      current_B2_function_capture_ratio = total_func_resid /
        max(total_raw_resid, .Machine$double.eps),
      current_B2_function_min_rank = function_rank_min,
      current_B2_function_max_condition = function_condition_max,
      stringsAsFactors = FALSE)

    split_stability <- list()
    for (domain in domains) {
      eig_odd <- gram_eigen(splitG$odd[[domain]], fit_id,
                            iteration, "minus_shared", domain)
      eig_even <- gram_eigen(splitG$even[[domain]], fit_id,
                             iteration, "minus_shared", domain)
      cosines <- svd(crossprod(eig_odd$vectors,
                               eig_even$vectors))$d
      top2 <- svd(crossprod(eig_odd$vectors[, 1:2, drop = FALSE],
                            eig_even$vectors[, 1:2, drop = FALSE]))$d
      split_stability[[domain]] <- c(
        top1_abs = abs(sum(eig_odd$vectors[, 1L] *
                             eig_even$vectors[, 1L])),
        top2_mean_sq = mean(top2^2),
        min_sq = min(cosines^2),
        mean_sq = mean(cosines^2))
      eigenvectors[[paste(fit_id, iteration, "minus_shared",
                          paste0(domain, "_odd"), sep = "__")]] <-
        eig_odd$vectors
      eigenvectors[[paste(fit_id, iteration, "minus_shared",
                          paste0(domain, "_even"), sep = "__")]] <-
        eig_even$vectors
    }
    for (stage in stages) for (domain in domains) {
      eig <- gram_eigen(G[[stage]][[domain]], fit_id,
                        iteration, stage, domain)
      eig$row$basis_capture_ratio <-
        total_basis_resid / max(total_raw_resid,
                                .Machine$double.eps)
      eig$row$current_B2_function_capture_ratio <-
        total_func_resid / max(total_raw_resid,
                               .Machine$double.eps)
      eig$row$odd_even_min_cosine_sq <-
        if (stage == "minus_shared")
          split_stability[[domain]][["min_sq"]] else NA_real_
      eig$row$odd_even_mean_cosine_sq <-
        if (stage == "minus_shared")
          split_stability[[domain]][["mean_sq"]] else NA_real_
      eig$row$odd_even_top1_abs_cosine <-
        if (stage == "minus_shared")
          split_stability[[domain]][["top1_abs"]] else NA_real_
      eig$row$odd_even_top2_mean_cosine_sq <-
        if (stage == "minus_shared")
          split_stability[[domain]][["top2_mean_sq"]] else NA_real_
      eig$row$design_times_min <- min(design_times)
      eig$row$design_times_median <- stats::median(design_times)
      eig$row$design_times_max <- max(design_times)
      eig$row$design_rank_min <- min(design_ranks)
      eig$row$design_rank_max <- max(design_ranks)
      eig$row$design_full_row_rank_fraction <-
        mean(design_ranks == design_times)
      eig$row$design_max_condition <- max(design_conditions)
      row_g <- row_g + 1L
      geometry[[row_g]] <- eig$row
      key <- paste(fit_id, iteration, stage, domain, sep = "__")
      eigenvectors[[key]] <- eig$vectors
    }
    cat("BLIND_POINT_PASS ", fit_id, " t=", iteration, "\n", sep = "")
  }
  rm(obs, fit, trace, snapshots, basis, C, qrC, Y, Z)
  invisible(gc())
}

files <- c("RECONSTRUCTION_PARITY.csv", "UPDATE_DRIVER_BLIND.csv",
           "SCORE_DRIVER_BLIND.csv", "BLOCK_ENERGY_BLIND.csv",
           "RESIDUAL_GEOMETRY_BLIND.csv",
           "RESIDUAL_EIGENVECTORS_BLIND.rds")
utils::write.csv(do.call(rbind, checks), file.path(out, files[[1L]]),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, drivers), file.path(out, files[[2L]]),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, scores), file.path(out, files[[3L]]),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, energy), file.path(out, files[[4L]]),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, geometry), file.path(out, files[[5L]]),
                 row.names = FALSE)
saveRDS(eigenvectors, file.path(out, files[[6L]]), version = 3L)
hashes <- vapply(file.path(out, files), sha256, character(1L))
writeLines(paste(hashes, files), file.path(out, "BLIND_TABLES.sha256"),
           useBytes = TRUE)
writeLines(c(paste0("mode=", mode),
             paste0("point_count=", row_e),
             paste0("source_sha256=",
                    sha256(sub("^--file=", "",
                               grep("^--file=", commandArgs(FALSE),
                                    value = TRUE))))),
           file.path(out, "BLIND_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7H_BLIND_COMPLETE points=", row_e, "\n", sep = "")
