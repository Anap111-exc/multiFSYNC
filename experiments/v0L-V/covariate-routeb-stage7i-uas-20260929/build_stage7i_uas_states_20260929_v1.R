# Stage 7I private, truth-blind state construction. This file does not fit a
# model, read simulation truth, alter the frozen source package, or write files.
# The caller verifies source/input hashes and supplies the exact frozen C_i.

.stage7i_fail <- function(message) stop(message, call. = FALSE)
.stage7i_check <- function(ok, message) {
  if (!isTRUE(ok)) .stage7i_fail(message)
  invisible(TRUE)
}
.stage7i_sym <- function(x) (x + t(x)) / 2
.stage7i_norm <- function(x) sqrt(sum(x^2))
.stage7i_finite <- function(x) is.numeric(x) && all(is.finite(x))
.stage7i_trap_weights <- function(grid) {
  grid <- as.numeric(grid)
  .stage7i_check(length(grid) >= 2L && all(is.finite(grid)) &&
                   all(diff(grid) > 0), "Invalid frozen time grid.")
  delta <- diff(grid)
  w <- numeric(length(grid))
  w[1L] <- delta[1L] / 2
  w[length(w)] <- delta[length(delta)] / 2
  if (length(w) > 2L) {
    w[2:(length(w) - 1L)] <-
      (delta[1:(length(delta) - 1L)] +
         delta[2:length(delta)]) / 2
  }
  w
}
.stage7i_solve <- function(A, b, label, min_rcond = 1e-12) {
  A <- .stage7i_sym(A)
  rc <- rcond(A)
  .stage7i_check(is.finite(rc) && rc > min_rcond,
                 paste0(label, " is ill-conditioned; rcond=", rc))
  out <- tryCatch(solve(A, b), error = function(e) e)
  .stage7i_check(!inherits(out, "error") && all(is.finite(out)),
                 paste0(label, " could not be solved."))
  out
}
.stage7i_grid_norm <- function(coef, C_g, weights) {
  sqrt(sum(weights * as.vector(C_g %*% coef)^2))
}
.stage7i_sign_grid <- function(coef, score, C_g) {
  curve <- as.vector(C_g %*% coef)
  first_max <- which.max(abs(curve))
  sign <- if (curve[first_max] < 0) -1 else 1
  list(coef = sign * coef, score = sign * score)
}
.stage7i_subject_score <- function(C, r, Psi, B, noise_var) {
  M <- ncol(Psi)
  T <- C %*% Psi
  Q <- B^2 * crossprod(T) + noise_var * diag(M)
  as.vector(.stage7i_solve(Q, B * crossprod(T, r),
                            "S subject score precision"))
}
.stage7i_g_energy <- function(C2, Psi, Xi) {
  sum(vapply(seq_along(C2), function(i) {
    g <- as.vector(C2[[i]] %*% Psi %*% Xi[i, ])
    sum(g^2)
  }, numeric(1L)))
}
.stage7i_fixed_residual_driver <- function(R, C2, Psi, Xi, tau, v) {
  n <- length(C2)
  p <- length(tau)
  contributions <- matrix(0, n, p)
  U_mean <- 0
  for (i in seq_len(n)) {
    gi <- as.vector(C2[[i]] %*% Psi %*% Xi[i, ])
    contributions[i, ] <- tau * as.vector(crossprod(gi, R[[i]]))
    U_mean <- U_mean + sum(gi^2)
  }
  h <- colSums(contributions)
  signed_envelope <- colSums(abs(contributions))
  v_subject_signed <- as.vector(contributions %*% v)
  v_driver <- sum(v_subject_signed)
  per_coordinate_coherence <- ifelse(signed_envelope > 0,
                                     abs(h) / signed_envelope, 0)
  row_norms <- sqrt(rowSums(contributions^2))
  list(
    h = h,
    contributions = contributions,
    U_mean = U_mean,
    h_norm = .stage7i_norm(h),
    normalized_h = .stage7i_norm(h) / sqrt(U_mean),
    v_signed_driver = v_driver,
    v_abs_cosine_to_driver =
      abs(v_driver) / max(.stage7i_norm(h), .Machine$double.eps),
    v_subject_signed = v_subject_signed,
    v_subject_cancellation_ratio =
      abs(v_driver) / max(sum(abs(v_subject_signed)),
                          .Machine$double.eps),
    signed_envelope = signed_envelope,
    per_coordinate_coherence = per_coordinate_coherence,
    global_subject_cancellation_ratio =
      .stage7i_norm(h) / max(sum(row_norms), .Machine$double.eps),
    subject_driver_norms = row_norms
  )
}
.stage7i_state_energy <- function(state, C2, p, N_obs) {
  Psi <- state$mu_q_nu_psi[[2L]][[1L]]
  Xi <- state$mu_q_xi[[2L]][[1L]]
  cov_psi <- state$Sigma_q_nu_psi[[2L]][[1L]]
  cov_xi <- state$Sigma_q_xi[[2L]][[1L]]
  b <- state$mu_q_b_specific[[2L]][, 1L]
  b_second <- state$term_b_specific[[2L]][, 1L]
  M <- ncol(Psi)
  G2 <- U_score_cov <- U_function_cov <- U_joint_cov <- time_second <- 0
  for (i in seq_along(C2)) {
    C_i <- C2[[i]]
    CTC <- crossprod(C_i)
    xi <- Xi[i, ]
    Sigma_xi <- cov_xi[[i]]
    Gi <- C_i %*% Psi
    H <- crossprod(Gi)
    H_cov_diag <- vapply(seq_len(M), function(m)
      sum(diag(CTC %*% cov_psi[[m]])), numeric(1L))
    G2 <- G2 + sum((Gi %*% xi)^2)
    U_score_cov <- U_score_cov + sum(diag(H %*% Sigma_xi))
    U_function_cov <- U_function_cov + sum(xi^2 * H_cov_diag)
    U_joint_cov <- U_joint_cov + sum(diag(Sigma_xi) * H_cov_diag)
    time_second <- time_second +
      sum(diag(H %*% (Sigma_xi + tcrossprod(xi)))) +
      sum(diag(Sigma_xi + tcrossprod(xi)) * H_cov_diag)
  }
  c(
    effective_loading_norm = .stage7i_norm(b),
    slab_mean_norm = .stage7i_norm(state$mu_q_normal_b[[2L]][, 1L]),
    mean_time_product_norm = sqrt(G2),
    mean_product_rms = .stage7i_norm(b) * sqrt(G2 / (p * N_obs)),
    full_posterior_product_second_rms =
      sqrt(sum(b_second) * time_second / (p * N_obs)),
    loading_second_sum = sum(b_second),
    U_mean = G2,
    U_score_cov = U_score_cov,
    U_function_cov = U_function_cov,
    U_joint_cov = U_joint_cov,
    score_rms_1 = sqrt(mean(Xi[, 1L]^2)),
    score_rms_2 = sqrt(mean(Xi[, 2L]^2))
  )
}
.stage7i_install_B2_mean <- function(state, b, Psi, Xi) {
  rho <- state$mu_q_gamma_b[[2L]][, 1L]
  V <- state$Sigma_q_normal_b[[2L]][, 1L]
  .stage7i_check(all(is.finite(rho)) && all(rho > 0) && all(rho <= 1),
                 "Stored B2 PPI must lie in (0, 1].")
  .stage7i_check(all(is.finite(V)) && all(V > 0),
                 "Stored B2 slab variances must be positive.")
  m <- b / rho
  second <- rho * (V + m^2)
  .stage7i_check(all(is.finite(m)) && all(is.finite(second)) &&
                   all(second > 0), "Injected B2 loading moments are invalid.")
  state$mu_q_normal_b[[2L]][, 1L] <- m
  state$mu_q_b_specific[[2L]][, 1L] <- rho * m
  state$term_b_specific[[2L]][, 1L] <- second
  state$mu_q_nu_psi[[2L]][[1L]] <- Psi
  state$mu_q_xi[[2L]][[1L]] <- Xi
  state
}

# Returns state lists with the exact Stage 7G snapshot field schema. The
# caller may persist them only after its manifest and input-hash preflight.
stage7i_build_uas_states <- function(observation, fit, snapshot50, C2) {
  numerical <- list(
    min_rcond = 1e-12,
    min_rank_ratio = 1e-8,
    min_grid_rcond = 1e-10,
    min_mean_norm = sqrt(.Machine$double.eps),
    match_atol = 1e-10,
    match_rtol = 1e-10,
    als_rounds = 3L
  )
  .stage7i_check(is.list(observation) && is.list(fit) &&
                   is.list(snapshot50) && is.list(C2),
                 "Expected observation, fit, snapshot50, and C2 lists.")
  .stage7i_check(identical(as.integer(snapshot50$iteration), 50L) &&
                   isTRUE(snapshot50$annealing),
                 "Expected a Stage 7G sweep-50 annealed snapshot.")
  # The fitted object may compact c(1, 1) to scalar L_s = 1.
  L_s2 <- if (length(fit$L_s) == 1L) fit$L_s[[1L]] else
    fit$L_s[[2L]]
  .stage7i_check(as.integer(fit$S) == 2L && as.integer(fit$L_f) == 1L &&
                   as.integer(L_s2) == 1L &&
                   as.integer(fit$M_s[[2L]][1L]) == 2L,
                 "Frozen C1 study-2 B2 dimensions changed.")
  K_total <- as.integer(fit$K_total)
  p <- as.integer(fit$p)
  n <- as.integer(fit$n_s[[2L]])
  .stage7i_check(K_total == fit$K + 2L &&
                   length(C2) == n && length(observation$Y[[2L]]) == n &&
                   is.matrix(observation$Z[[2L]]) &&
                   nrow(observation$Z[[2L]]) == n &&
                   ncol(observation$Z[[2L]]) == fit$d,
                 "Frozen observation dimensions changed.")
  C_g <- fit$C_g
  time_g <- fit$time_g
  .stage7i_check(is.matrix(C_g) && ncol(C_g) == K_total &&
                   nrow(C_g) == length(time_g),
                 "Frozen dense-grid basis is invalid.")
  weights <- .stage7i_trap_weights(time_g)
  Y2 <- lapply(observation$Y[[2L]], function(y) do.call(cbind, y))
  N_obs <- sum(vapply(C2, nrow, integer(1L)))
  for (i in seq_len(n)) {
    .stage7i_check(is.matrix(C2[[i]]) && ncol(C2[[i]]) == K_total &&
                     is.matrix(Y2[[i]]) &&
                     identical(dim(Y2[[i]]), c(nrow(C2[[i]]), p)) &&
                     .stage7i_finite(C2[[i]]) && .stage7i_finite(Y2[[i]]),
                   paste0("Frozen study-2 design/response invalid, subject ", i))
    expected <- fit$list_cp_C[[2L]][[i]]
    observed <- crossprod(C2[[i]])
    .stage7i_check(is.matrix(expected) &&
                     identical(dim(expected), dim(observed)) &&
                     max(abs(observed - expected)) <=
                       1e-9 * (1 + max(abs(expected))),
                   paste0("C_i does not match frozen C_i^T C_i, subject ", i))
  }
  coeff_mean <- do.call(cbind, snapshot50$mu_q_nu_mu[[2L]])
  .stage7i_check(identical(dim(coeff_mean), c(K_total, p)),
                 "Study-2 mean coefficients have unexpected shape.")
  coeff_beta <- lapply(seq_len(fit$d), function(r)
    do.call(cbind, lapply(snapshot50$mu_q_nu_beta, `[[`, r)))
  .stage7i_check(all(vapply(coeff_beta, function(x)
    identical(dim(x), c(K_total, p)), logical(1L))),
    "Beta coefficient shapes changed.")
  phi <- snapshot50$mu_q_nu_phi[[1L]]
  zeta <- snapshot50$mu_q_zeta[[2L]][[1L]]
  a <- snapshot50$mu_q_a[, 1L]
  .stage7i_check(is.matrix(phi) && nrow(phi) == K_total &&
                   is.matrix(zeta) && nrow(zeta) == n &&
                   ncol(phi) == ncol(zeta) && length(a) == p,
                 "Shared-state shapes changed.")
  R <- vector("list", n)
  G_R <- matrix(0, p, p)
  for (i in seq_len(n)) {
    Ci <- C2[[i]]
    mu_i <- Ci %*% coeff_mean
    beta_i <- matrix(0, nrow(Ci), p)
    for (r in seq_len(fit$d)) {
      beta_i <- beta_i + observation$Z[[2L]][i, r] *
        (Ci %*% coeff_beta[[r]])
    }
    shared_time <- as.vector(Ci %*% phi %*% zeta[i, ])
    R[[i]] <- Y2[[i]] - mu_i - beta_i - tcrossprod(shared_time, a)
    .stage7i_check(all(is.finite(R[[i]])),
                   paste0("Nonfinite leave-B2-out residual, subject ", i))
    G_R <- G_R + crossprod(R[[i]])
  }
  eig_R <- eigen(.stage7i_sym(G_R), symmetric = TRUE)
  .stage7i_check(is.finite(eig_R$values[1L]) && eig_R$values[1L] > 0,
                 "No leading residual variable signal.")
  v <- eig_R$vectors[, 1L]
  if (v[which.max(abs(v))] < 0) v <- -v
  B_target <- sqrt(eig_R$values[1L] / N_obs)
  tau <- snapshot50$mu_q_recip_sigsq_eps[2L, ]
  .stage7i_check(length(tau) == p && all(is.finite(tau)) &&
                   all(tau > 0), "Invalid frozen study-2 noise precision.")
  sigma_v2 <- sum(v^2 / tau)
  lambda <- sigma_v2 / B_target^2
  .stage7i_check(is.finite(lambda) && lambda > 0,
                 "Spectral time penalty is invalid.")
  Sigma_beta <- fit$list_hyper$Sigma_beta
  .stage7i_check(is.matrix(Sigma_beta) &&
                   identical(dim(Sigma_beta), c(2L, 2L)) &&
                   .stage7i_finite(Sigma_beta),
                 "Frozen linear prior covariance is invalid.")
  inv_beta <- .stage7i_solve(Sigma_beta, diag(2L),
                             "frozen linear prior covariance")
  P <- matrix(0, K_total, K_total)
  P[1:2, 1:2] <- inv_beta
  P[3:K_total, 3:K_total] <- diag(K_total - 2L)
  CTC <- lapply(C2, crossprod)
  H_data <- Reduce(`+`, CTC)
  H <- .stage7i_sym(H_data + lambda * P)
  q <- lapply(R, function(Ri) as.vector(Ri %*% v))
  delta_rhs <- Reduce(`+`, Map(function(Ci, qi) crossprod(Ci, qi), C2, q))
  delta <- as.vector(.stage7i_solve(H, delta_rhs, "pooled time mean"))
  r <- Map(function(qi, Ci) as.vector(qi - Ci %*% delta), q, C2)
  D <- matrix(0, K_total, K_total)
  for (i in seq_len(n)) {
    x <- as.vector(crossprod(C2[[i]], r[[i]]))
    D <- D + tcrossprod(x)
  }
  H_chol <- tryCatch(chol(H), error = function(e) e)
  .stage7i_check(!inherits(H_chol, "error"),
                 "Pooled time metric is not positive definite.")
  whitened <- solve(t(H_chol), D) %*% solve(H_chol)
  eig_time <- eigen(.stage7i_sym(whitened), symmetric = TRUE)
  d <- eig_time$values
  .stage7i_check(all(is.finite(d)) && d[1L] > 0 &&
                   d[2L] / d[1L] > numerical$min_rank_ratio,
                 "Spectral time rank is below the fixed threshold.")
  .stage7i_check((d[2L] - d[3L]) / d[2L] > numerical$min_rank_ratio,
                 "Second spectral time direction is unresolved.")
  Psi_S <- solve(H_chol, eig_time$vectors[, 1:2, drop = FALSE])
  for (m in 1:2) {
    grid_norm <- .stage7i_grid_norm(Psi_S[, m], C_g, weights)
    .stage7i_check(is.finite(grid_norm) &&
                     grid_norm > numerical$min_mean_norm,
                   "A spectral time direction has zero grid energy.")
    Psi_S[, m] <- Psi_S[, m] / grid_norm
    signed <- .stage7i_sign_grid(Psi_S[, m], 1, C_g)
    Psi_S[, m] <- signed$coef
  }
  grid_gram <- crossprod(C_g, weights * C_g)
  time_gram <- crossprod(Psi_S, grid_gram %*% Psi_S)
  .stage7i_check(rcond(time_gram) > numerical$min_grid_rcond,
                 "Spectral time directions are nearly collinear on grid.")
  Xi_S <- t(vapply(seq_len(n), function(i)
    .stage7i_subject_score(C2[[i]], r[[i]], Psi_S,
                           B_target, sigma_v2), numeric(2L)))
  for (round in seq_len(numerical$als_rounds)) {
    for (m in 1:2) {
      other <- 3L - m
      lhs <- lambda * P
      rhs <- numeric(K_total)
      for (i in seq_len(n)) {
        Ci <- C2[[i]]
        xi_m <- Xi_S[i, m]
        partial <- r[[i]] / B_target -
          as.vector(Ci %*% Psi_S[, other]) * Xi_S[i, other]
        lhs <- lhs + xi_m^2 * CTC[[i]]
        rhs <- rhs + xi_m * as.vector(crossprod(Ci, partial))
      }
      Psi_S[, m] <- as.vector(.stage7i_solve(lhs, rhs,
        paste0("spectral time ALS m=", m, " round=", round)))
    }
    Xi_S <- t(vapply(seq_len(n), function(i)
      .stage7i_subject_score(C2[[i]], r[[i]], Psi_S,
                             B_target, sigma_v2), numeric(2L)))
  }
  Xi_S <- sweep(Xi_S, 2L, colMeans(Xi_S), "-")
  score_scale_S <- sqrt(colMeans(Xi_S^2))
  .stage7i_check(all(is.finite(score_scale_S)) &&
                   all(score_scale_S > numerical$min_mean_norm),
                 "S score column is degenerate after centering.")
  Xi_S <- sweep(Xi_S, 2L, score_scale_S, "/")
  Psi_S <- sweep(Psi_S, 2L, score_scale_S, "*")
  for (m in 1:2) {
    signed <- .stage7i_sign_grid(Psi_S[, m], Xi_S[, m], C_g)
    Psi_S[, m] <- signed$coef
    Xi_S[, m] <- signed$score
  }
  final_time_gram <- crossprod(Psi_S, grid_gram %*% Psi_S)
  .stage7i_check(is.finite(rcond(final_time_gram)) &&
                   rcond(final_time_gram) > numerical$min_grid_rcond,
                 "Final S time directions are nearly collinear on grid.")
  G_S <- sqrt(.stage7i_g_energy(C2, Psi_S, Xi_S))
  .stage7i_check(is.finite(G_S) && G_S > numerical$min_mean_norm,
                 "S mean time product is degenerate.")
  b_U <- snapshot50$mu_q_b_specific[[2L]][, 1L]
  Psi_U <- snapshot50$mu_q_nu_psi[[2L]][[1L]]
  Xi_U <- snapshot50$mu_q_xi[[2L]][[1L]]
  .stage7i_check(length(b_U) == p &&
                   identical(dim(Psi_U), c(K_total, 2L)) &&
                   identical(dim(Xi_U), c(n, 2L)) &&
                   all(is.finite(b_U)) && all(is.finite(Psi_U)) &&
                   all(is.finite(Xi_U)),
                 "Original B2 mean-state shapes changed.")
  b_U_norm <- .stage7i_norm(b_U)
  G_U <- sqrt(.stage7i_g_energy(C2, Psi_U, Xi_U))
  old_score_scale <- sqrt(colMeans(Xi_U^2))
  .stage7i_check(is.finite(b_U_norm) &&
                   b_U_norm > numerical$min_mean_norm &&
                   is.finite(G_U) && G_U > numerical$min_mean_norm &&
                   all(is.finite(old_score_scale)) &&
                   all(old_score_scale > numerical$min_mean_norm),
                 "Original B2 directions are numerically zero.")
  b_A <- B_target * b_U / b_U_norm
  Xi_A <- sweep(Xi_U, 2L, old_score_scale, "/")
  Psi_A <- sweep(Psi_U, 2L, old_score_scale * G_S / G_U, "*")
  .stage7i_check(all(is.finite(b_A)) && all(is.finite(Xi_A)) &&
                   all(is.finite(Psi_A)) && all(is.finite(Psi_S)) &&
                   all(is.finite(Xi_S)),
                 "An injected B2 mean contains nonfinite values.")
  temporal_scaling_error <- max(vapply(seq_len(n), function(i) {
    old_g <- as.vector(C2[[i]] %*% Psi_U %*% Xi_U[i, ])
    new_g <- as.vector(C2[[i]] %*% Psi_A %*% Xi_A[i, ])
    max(abs(new_g - (G_S / G_U) * old_g))
  }, numeric(1L)))
  .stage7i_check(temporal_scaling_error <=
                   numerical$match_atol + numerical$match_rtol * G_S,
                 "A does not preserve the original time-product directions.")
  state_U <- snapshot50
  state_A <- .stage7i_install_B2_mean(snapshot50, b_A, Psi_A, Xi_A)
  state_S <- .stage7i_install_B2_mean(snapshot50, B_target * v,
                                      Psi_S, Xi_S)
  energies <- rbind(
    U = .stage7i_state_energy(state_U, C2, p, N_obs),
    A = .stage7i_state_energy(state_A, C2, p, N_obs),
    S = .stage7i_state_energy(state_S, C2, p, N_obs)
  )
  .stage7i_check(all(is.finite(energies)) && all(energies >= 0),
                 "A branch has invalid initial energy diagnostics.")
  fixed_drivers <- list(
    U = .stage7i_fixed_residual_driver(R, C2, Psi_U, Xi_U, tau, v),
    A = .stage7i_fixed_residual_driver(R, C2, Psi_A, Xi_A, tau, v),
    S = .stage7i_fixed_residual_driver(R, C2, Psi_S, Xi_S, tau, v)
  )
  driver_kappa <- G_S / G_U
  h_identity_error <- max(abs(fixed_drivers$A$h -
                                driver_kappa * fixed_drivers$U$h))
  U_identity_error <- abs(fixed_drivers$A$U_mean -
                            driver_kappa^2 * fixed_drivers$U$U_mean)
  normalized_h_identity_error <- abs(fixed_drivers$A$normalized_h -
                                       fixed_drivers$U$normalized_h)
  .stage7i_check(
    is.finite(h_identity_error) &&
      h_identity_error <= 1e-8 +
        1e-10 * max(abs(fixed_drivers$A$h)) &&
      is.finite(U_identity_error) &&
      U_identity_error <= 1e-8 +
        1e-10 * fixed_drivers$A$U_mean &&
      is.finite(normalized_h_identity_error) &&
      normalized_h_identity_error <= 1e-8,
    "A did not preserve the exact fixed-residual driver identity."
  )
  match_close <- function(x, y) {
    abs(x - y) <= numerical$match_atol +
      numerical$match_rtol * max(abs(x), abs(y))
  }
  .stage7i_check(match_close(energies["A", "effective_loading_norm"],
                             energies["S", "effective_loading_norm"]) &&
                   match_close(energies["A", "mean_time_product_norm"],
                               energies["S", "mean_time_product_norm"]) &&
                   match_close(energies["A", "mean_product_rms"],
                               energies["S", "mean_product_rms"]),
                 "A and S do not match the three preregistered amplitudes.")
  .stage7i_check(identical(names(state_U), names(state_A)) &&
                   identical(names(state_U), names(state_S)),
                 "Injected state field schema changed.")
  frozen_fields <- setdiff(names(state_U), c(
    "mu_q_normal_b", "mu_q_b_specific", "term_b_specific",
    "mu_q_nu_psi", "mu_q_xi"))
  .stage7i_check(all(vapply(frozen_fields, function(name)
    identical(state_U[[name]], state_A[[name]]) &&
      identical(state_U[[name]], state_S[[name]]), logical(1L))),
    "A or S altered an unrelated snapshot field.")
  .stage7i_check(identical(state_U$mu_q_gamma_b, state_A$mu_q_gamma_b) &&
                   identical(state_U$mu_q_gamma_b, state_S$mu_q_gamma_b) &&
                   identical(state_U$Sigma_q_normal_b,
                             state_A$Sigma_q_normal_b) &&
                   identical(state_U$Sigma_q_normal_b,
                             state_S$Sigma_q_normal_b) &&
                   identical(state_U$Sigma_q_nu_psi,
                             state_A$Sigma_q_nu_psi) &&
                   identical(state_U$Sigma_q_nu_psi,
                             state_S$Sigma_q_nu_psi) &&
                   identical(state_U$Sigma_q_xi,
                             state_A$Sigma_q_xi) &&
                   identical(state_U$Sigma_q_xi,
                             state_S$Sigma_q_xi),
                 "A or S changed a preserved PPI or covariance field.")
  for (name in c("mu_q_normal_b", "mu_q_b_specific",
                 "term_b_specific", "mu_q_nu_psi", "mu_q_xi")) {
    .stage7i_check(identical(state_U[[name]][[1L]],
                             state_A[[name]][[1L]]) &&
                     identical(state_U[[name]][[1L]],
                               state_S[[name]][[1L]]),
                   paste0("A or S altered the study-1 ", name, " state."))
  }
  scalar_projected_sse <- function(Psi, Xi) {
    sum(vapply(seq_len(n), function(i) {
      error <- r[[i]] - B_target *
        as.vector(C2[[i]] %*% Psi %*% Xi[i, ])
      sum(error^2)
    }, numeric(1L)))
  }
  details <- list(
    protocol = "Stage7I_UAS_20260929_v1",
    numerical = numerical,
    n_subject = n,
    N_obs = N_obs,
    p = p,
    K_total = K_total,
    leading_variable_eigenvalue = eig_R$values[1L],
    second_variable_eigenvalue = eig_R$values[2L],
    variable_direction = v,
    variable_direction_sign_anchor = which.max(abs(v)),
    variable_direction_abs_cosine_to_old_loading =
      abs(sum(v * b_U)) / b_U_norm,
    B_target = B_target,
    projected_noise_variance = sigma_v2,
    lambda = lambda,
    P = P,
    pooled_H_rcond = rcond(H),
    time_eigenvalues = d,
    time_grid_gram_initial = time_gram,
    time_grid_gram_final = final_time_gram,
    pooled_mean_coef = delta,
    pooled_mean_removed_energy =
      sum(vapply(seq_len(n), function(i)
        sum(as.vector(C2[[i]] %*% delta)^2), numeric(1L))) /
      sum(vapply(q, function(qi) sum(qi^2), numeric(1L))),
    projected_centered_energy = sum(vapply(r, function(ri)
      sum(ri^2), numeric(1L))),
    S_projected_centered_sse = scalar_projected_sse(Psi_S, Xi_S),
    initial_energy = energies,
    fixed_residual_driver = fixed_drivers,
    A_driver_scaling_kappa = driver_kappa,
    A_driver_h_identity_max_abs_error = h_identity_error,
    A_driver_U_mean_identity_abs_error = U_identity_error,
    A_driver_normalized_h_identity_abs_error =
      normalized_h_identity_error,
    score_scale_S_before_normalization = score_scale_S,
    score_scale_U_before_normalization = old_score_scale,
    A_temporal_scaling_max_abs_error = temporal_scaling_error,
    amplitudes_matched = TRUE
  )
  list(U = state_U, A = state_A, S = state_S, blind_diagnostics = details)
}
