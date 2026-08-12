# =============================================================================
# orthonormalise_multi.R — Posterior-moment FPCA post-processing
#
# Converts the unconstrained CAVI representation (eta, chi, theta, kappa) to
# the theoretical FPCA representation.  The covariance operator uses complete
# variational second moments of both scores and time-function coefficients.
# No rotation is performed across sparse factor columns.
# =============================================================================

.trap_weights <- function(time_g) {
  time_g <- as.numeric(time_g)
  if (length(time_g) < 2L || any(!is.finite(time_g)) ||
      any(diff(time_g) <= 0)) {
    stop("time_g must contain at least two finite, strictly increasing points.")
  }
  delta <- diff(time_g)
  w <- numeric(length(time_g))
  w[1] <- delta[1] / 2
  w[length(w)] <- delta[length(delta)] / 2
  if (length(w) > 2L) {
    w[2:(length(w) - 1L)] <-
      (delta[1:(length(delta) - 1L)] + delta[2:length(delta)]) / 2
  }
  w
}

.symmetrise <- function(x) (x + t(x)) / 2

.zero_coef_cov <- function(coef_mean) {
  replicate(ncol(coef_mean),
            matrix(0, nrow(coef_mean), nrow(coef_mean)),
            simplify = FALSE)
}

.posterior_factor_fpca <- function(C_g, weights,
                                   coef_mean, coef_cov,
                                   score_mean, score_cov,
                                   pve_cutoff = 0.99,
                                   tol = 1e-10) {
  coef_mean <- as.matrix(coef_mean)
  score_mean <- as.matrix(score_mean)
  M <- ncol(coef_mean)
  if (ncol(score_mean) != M) {
    stop("Score and time-function component dimensions do not agree.")
  }
  N <- nrow(score_mean)
  if (length(score_cov) != N) {
    stop("There must be one score posterior covariance per individual.")
  }
  if (is.null(coef_cov)) coef_cov <- .zero_coef_cov(coef_mean)
  if (length(coef_cov) != M) {
    stop("There must be one coefficient posterior covariance per component.")
  }

  # R = N^{-1} sum_i E(eta_i eta_i^T).
  R_score <- matrix(0, M, M)
  for (i in seq_len(N)) {
    mui <- as.numeric(score_mean[i, ])
    R_score <- R_score + score_cov[[i]] + tcrossprod(mui)
  }
  R_score <- .symmetrise(R_score / N)

  # Coefficient-space representation of the grid covariance kernel:
  # B = E(H) R E(H)^T + sum_m R_mm Var(H_m), K_g = C_g B C_g^T.
  B_coef <- coef_mean %*% R_score %*% t(coef_mean)
  for (m in seq_len(M)) {
    B_coef <- B_coef + R_score[m, m] * coef_cov[[m]]
  }
  B_coef <- .symmetrise(B_coef)

  # Weighted covariance-operator eigendecomposition in coefficient space.
  # This is equivalent to eigen(sqrt(W) K_g sqrt(W)) but costs O((K+2)^3)
  # instead of O(n_g^3).
  G <- crossprod(C_g, weights * C_g)
  G <- .symmetrise(G)
  eg <- eigen(G, symmetric = TRUE)
  g_tol <- tol * max(1, max(abs(eg$values)))
  g_keep <- which(eg$values > g_tol)
  if (length(g_keep) == 0L) stop("The weighted spline Gram matrix has zero rank.")
  Ug <- eg$vectors[, g_keep, drop = FALSE]
  dg <- eg$values[g_keep]
  G_half <- Ug %*% (sqrt(dg) * t(Ug))
  G_inv_half <- Ug %*% ((1 / sqrt(dg)) * t(Ug))

  eh <- eigen(.symmetrise(G_half %*% B_coef %*% G_half),
              symmetric = TRUE)
  e_tol <- tol * max(1, max(abs(eh$values)))
  e_keep <- which(eh$values > e_tol)
  if (length(e_keep) == 0L) {
    stop("Posterior factor covariance has no positive direction.")
  }
  eig_all <- eh$values[e_keep]
  full_integrated_variance <- sum(eig_all)

  # Each fitted factor contains M computational time components.  Posterior
  # uncertainty in their spline coefficients can make E_q(K) have rank > M,
  # but those extra directions quantify parameter uncertainty rather than new
  # population FPCA components.  The reported point estimate is therefore the
  # best rank-M approximation of the complete posterior second moment.
  rank_cap <- min(M, length(eig_all))
  eig_raw <- eig_all[seq_len(rank_cap)]
  coef_eigen <- G_inv_half %*%
    eh$vectors[, e_keep[seq_len(rank_cap)], drop = FALSE]
  functions <- C_g %*% coef_eigen

  # Remove numerical norm drift and impose a deterministic component sign.
  norms <- sqrt(colSums(weights * functions^2))
  functions <- sweep(functions, 2, norms, "/")
  coef_eigen <- sweep(coef_eigen, 2, norms, "/")
  component_sign <- vapply(seq_len(ncol(functions)), function(k) {
    idx <- which.max(abs(functions[, k]))
    if (functions[idx, k] < 0) -1 else 1
  }, numeric(1))
  functions <- sweep(functions, 2, component_sign, "*")
  coef_eigen <- sweep(coef_eigen, 2, component_sign, "*")

  # c^2 is the integrated variance of the retained best rank-M point estimate.
  # Dividing factor scores by c and multiplying the corresponding loading
  # column by c gives unit total integrated variance without counting omitted
  # posterior-uncertainty directions as population components.
  scale_sq <- sum(eig_raw)
  if (!is.finite(scale_sq) || scale_sq <= tol) {
    stop("Posterior factor has non-positive total integrated variance.")
  }
  factor_scale <- sqrt(scale_sq)

  P <- crossprod(C_g, weights * functions)
  mean_map <- t(coef_mean) %*% P
  score_hat_raw <- score_mean %*% mean_map
  cov_hat <- vector("list", N)
  for (i in seq_len(N)) {
    mui <- as.numeric(score_mean[i, ])
    Ri <- score_cov[[i]] + tcrossprod(mui)
    second <- t(mean_map) %*% Ri %*% mean_map
    for (m in seq_len(M)) {
      second <- second + Ri[m, m] *
        crossprod(P, coef_cov[[m]] %*% P)
    }
    cov_i <- .symmetrise(second - tcrossprod(score_hat_raw[i, ]))
    cov_hat[[i]] <- cov_i / scale_sq
  }
  score_hat <- score_hat_raw / factor_scale
  eigenvalues <- pmax(eig_raw / scale_sq, 0)
  pve <- cumsum(eigenvalues) / sum(eigenvalues)
  effective_M <- which(pve >= pve_cutoff)[1]
  if (is.na(effective_M)) effective_M <- length(eigenvalues)
  keep <- seq_len(effective_M)

  list(
    functions = functions[, keep, drop = FALSE],
    scores = score_hat[, keep, drop = FALSE],
    score_cov = lapply(cov_hat, function(x) x[keep, keep, drop = FALSE]),
    eigenvalues = eigenvalues[keep],
    cumulative_pve = 100 * pve[keep],
    effective_M = effective_M,
    factor_scale = factor_scale,
    R_score = R_score,
    integrated_variance = scale_sq,
    full_posterior_integrated_variance = full_integrated_variance,
    omitted_uncertainty_variance =
      max(full_integrated_variance - scale_sq, 0),
    rank_cap = M
  )
}

.loading_sign <- function(x, tol = 1e-12) {
  if (length(x) == 0L || max(abs(x)) <= tol) return(1)
  idx <- which.max(abs(x))
  if (x[idx] < 0) -1 else 1
}

.factor_ppi <- function(gamma) {
  if (is.null(gamma) || ncol(as.matrix(gamma)) == 0L) return(numeric(0))
  gamma <- pmin(pmax(gamma, 0), 1)
  1 - apply(1 - gamma, 2, prod)
}

#' Recover theoretical FPCA parameters from the CAVI representation
#'
#' Legacy input names retain their historical R interface: zeta/xi are the
#' working scores eta/chi, while nu_phi/nu_psi are coefficients of theta/kappa.
#' The returned loadings are scaled, signed and sorted canonical loadings.
#'
#' @param Sigma_q_nu_phi Posterior covariance of shared time coefficients.
#' @param Sigma_q_nu_psi Posterior covariance of specific time coefficients.
#' @param zeta_true Deprecated simulation-only argument.  Truth-based alignment
#'   is deliberately not performed in the model post-processing routine.
#' @noRd
#' @export
orthonormalise_multi <- function(C_g, time_g,
                                  mu_q_nu_mu,
                                  mu_q_nu_beta = NULL,
                                  mu_q_nu_phi, mu_q_nu_psi,
                                  mu_q_zeta, Sigma_q_zeta,
                                  mu_q_xi, Sigma_q_xi,
                                  mu_q_a, mu_q_b_specific,
                                  mu_q_gamma_a, mu_q_gamma_b,
                                  S, n_s, p, L_f, L_s, M_f, M_s,
                                  d = if (is.null(mu_q_nu_beta)) 0L else length(mu_q_nu_beta[[1]]),
                                  response_center = NULL,
                                  response_scale = NULL,
                                  Sigma_q_nu_phi = NULL,
                                  Sigma_q_nu_psi = NULL,
                                  zeta_true = NULL,
                                  pve_cutoff = 0.99) {
  if (!is.null(zeta_true)) {
    warning("zeta_true is ignored here; truth-based matching belongs in simulation evaluation code.")
  }
  weights <- .trap_weights(time_g)

  if (is.null(response_center)) response_center <- rep(0, p)
  if (is.null(response_scale)) response_scale <- rep(1, p)
  if (length(response_center) != p || length(response_scale) != p ||
      any(!is.finite(response_center)) || any(!is.finite(response_scale)) ||
      any(response_scale <= 0)) {
    stop("response_center and response_scale must be finite vectors of length p; scales must be positive.")
  }

  # Original-scale mean functions, indexed as study -> variable.
  list_mu_hat <- vector("list", S)
  for (s in seq_len(S)) {
    list_mu_hat[[s]] <- vector("list", p)
    for (j in seq_len(p)) {
      list_mu_hat[[s]][[j]] <- response_center[j] +
        response_scale[j] * as.vector(C_g %*% mu_q_nu_mu[[s]][[j]])
    }
  }

  # Original-scale covariate coefficient functions, indexed as variable ->
  # covariate.  Covariates themselves are not rescaled here.
  list_beta_hat <- NULL
  if (d > 0L && !is.null(mu_q_nu_beta)) {
    list_beta_hat <- lapply(seq_len(p), function(j) {
      lapply(seq_len(d), function(r) {
        response_scale[j] * as.vector(C_g %*% mu_q_nu_beta[[j]][[r]])
      })
    })
  }

  shared <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    score_mean <- do.call(rbind, lapply(seq_len(S), function(s) {
      mu_q_zeta[[s]][[l]]
    }))
    score_cov <- do.call(c, lapply(seq_len(S), function(s) {
      Sigma_q_zeta[[s]][[l]]
    }))
    shared[[l]] <- .posterior_factor_fpca(
      C_g, weights,
      coef_mean = mu_q_nu_phi[[l]],
      coef_cov = if (is.null(Sigma_q_nu_phi)) NULL else Sigma_q_nu_phi[[l]],
      score_mean = score_mean,
      score_cov = score_cov,
      pve_cutoff = pve_cutoff)
  }

  scale_shared <- vapply(shared, `[[`, numeric(1), "factor_scale")
  loadings_shared <- sweep(mu_q_a, 1, response_scale, "*")
  loadings_shared <- sweep(loadings_shared, 2, scale_shared, "*")
  sign_shared <- vapply(seq_len(L_f), function(l) {
    .loading_sign(loadings_shared[, l])
  }, numeric(1))
  loadings_shared <- sweep(loadings_shared, 2, sign_shared, "*")
  for (l in seq_len(L_f)) shared[[l]]$scores <- shared[[l]]$scores * sign_shared[l]
  strength_shared <- colSums(loadings_shared^2)
  order_shared <- order(-strength_shared, seq_len(L_f))
  shared <- shared[order_shared]
  loadings_shared <- loadings_shared[, order_shared, drop = FALSE]
  gamma_shared <- mu_q_gamma_a[, order_shared, drop = FALSE]
  scale_shared <- scale_shared[order_shared]
  sign_shared <- sign_shared[order_shared]

  specific <- vector("list", S)
  loadings_specific <- vector("list", S)
  gamma_specific <- vector("list", S)
  order_specific <- vector("list", S)
  scale_specific <- vector("list", S)
  sign_specific <- vector("list", S)
  if (.has_specific(L_s) && !is.null(mu_q_nu_psi)) {
    for (s in seq_len(S)) {
      L_ss <- .L_s_at(L_s, s)
      specific[[s]] <- vector("list", L_ss)
      if (L_ss == 0L) {
        loadings_specific[[s]] <- matrix(0, p, 0L)
        gamma_specific[[s]] <- matrix(0, p, 0L)
        order_specific[[s]] <- integer(0)
        scale_specific[[s]] <- numeric(0)
        sign_specific[[s]] <- numeric(0)
        next
      }
      for (l in seq_len(L_ss)) {
        specific[[s]][[l]] <- .posterior_factor_fpca(
          C_g, weights,
          coef_mean = mu_q_nu_psi[[s]][[l]],
          coef_cov = if (is.null(Sigma_q_nu_psi)) NULL else Sigma_q_nu_psi[[s]][[l]],
          score_mean = mu_q_xi[[s]][[l]],
          score_cov = Sigma_q_xi[[s]][[l]],
          pve_cutoff = pve_cutoff)
      }
      scales <- vapply(specific[[s]], `[[`, numeric(1), "factor_scale")
      Bhat <- sweep(mu_q_b_specific[[s]], 1, response_scale, "*")
      Bhat <- sweep(Bhat, 2, scales, "*")
      signs <- vapply(seq_len(L_ss), function(l) .loading_sign(Bhat[, l]),
                      numeric(1))
      Bhat <- sweep(Bhat, 2, signs, "*")
      for (l in seq_len(L_ss)) {
        specific[[s]][[l]]$scores <- specific[[s]][[l]]$scores * signs[l]
      }
      ord <- order(-colSums(Bhat^2), seq_len(L_ss))
      specific[[s]] <- specific[[s]][ord]
      loadings_specific[[s]] <- Bhat[, ord, drop = FALSE]
      gamma_specific[[s]] <- mu_q_gamma_b[[s]][, ord, drop = FALSE]
      order_specific[[s]] <- ord
      scale_specific[[s]] <- scales[ord]
      sign_specific[[s]] <- signs[ord]
    }
  }

  list(
    list_Phi_hat = lapply(shared, `[[`, "functions"),
    list_Zeta_hat = lapply(shared, `[[`, "scores"),
    list_Cov_zeta_hat = lapply(shared, `[[`, "score_cov"),
    list_eigenvalues = lapply(shared, `[[`, "eigenvalues"),
    list_cumulated_pve = lapply(shared, `[[`, "cumulative_pve"),
    list_effective_M = lapply(shared, `[[`, "effective_M"),
    list_full_posterior_integrated_variance =
      lapply(shared, `[[`, "full_posterior_integrated_variance"),
    list_omitted_uncertainty_variance =
      lapply(shared, `[[`, "omitted_uncertainty_variance"),
    list_rank_cap = lapply(shared, `[[`, "rank_cap"),
    list_Phi_hat_spec = lapply(specific, function(x) lapply(x, `[[`, "functions")),
    list_Zeta_hat_spec = lapply(specific, function(x) lapply(x, `[[`, "scores")),
    list_Cov_zeta_hat_spec = lapply(specific, function(x) lapply(x, `[[`, "score_cov")),
    list_eigenvalues_spec = lapply(specific, function(x) lapply(x, `[[`, "eigenvalues")),
    list_pve_spec = lapply(specific, function(x) lapply(x, `[[`, "cumulative_pve")),
    list_effective_M_spec = lapply(specific, function(x) lapply(x, `[[`, "effective_M")),
    list_full_posterior_integrated_variance_spec = lapply(
      specific, function(x) lapply(x, `[[`, "full_posterior_integrated_variance")),
    list_omitted_uncertainty_variance_spec = lapply(
      specific, function(x) lapply(x, `[[`, "omitted_uncertainty_variance")),
    list_rank_cap_spec = lapply(specific, function(x) lapply(x, `[[`, "rank_cap")),
    list_mu_hat = list_mu_hat,
    list_mu_hat_legacy = list_mu_hat[[1]],
    list_beta_hat = list_beta_hat,
    mu_q_a_hat = loadings_shared,
    mu_q_b_specific_hat = loadings_specific,
    mu_q_gamma_a_hat = gamma_shared,
    mu_q_gamma_b_hat = gamma_specific,
    factor_scale_shared = scale_shared,
    factor_scale_specific = scale_specific,
    factor_sign_shared = sign_shared,
    factor_sign_specific = sign_specific,
    factor_order_shared = order_shared,
    factor_order_specific = order_specific,
    factor_ppi_shared = .factor_ppi(gamma_shared),
    factor_ppi_specific = if (.has_specific(L_s)) lapply(gamma_specific, .factor_ppi) else NULL
  )
}
