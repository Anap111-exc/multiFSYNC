# =============================================================================
# orthonormalise_multi.R — Post-processing for multi-study functional factor model
#
# Covers derivations.md §6 (orthonormalisation).
# Produces orthogonal eigenfunctions, rotated scores, and PVEs.
#
# Procedure per shared factor l (M_f[l] basis functions):
#   1. Evaluate eigenfunctions on dense grid: Phi_g = C_g * nu_phi_l   (n_g x M_f[l])
#   2. SVD: Phi_g = U D V^T
#   3. Rotate scores: Zeta_rot = Zeta * V * D
#   4. Eigen-decompose covariance of Zeta_rot: Q Lambda Q^T
#   5. Orthonormal eigenfunctions: Phi_hat = U * Q * sqrt(Lambda)
#   6. Orthonormal scores: Zeta_hat = Zeta_rot * Q * 1/sqrt(Lambda)
#   7. Normalise via trapezoidal integration
#   8. Compute PVE
#
# For specific factors: same procedure per study   (M_s[[s]][l] basis functions).
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Orthonormalise shared and specific factors post-CAVI
#'
#' @param C_g Dense grid design matrix (n_g x (K+2))
#' @param time_g Dense grid time points (length n_g)
#' @param mu_q_nu_mu List[[s]][[j]] of mean function coefficients
#' @param mu_q_nu_phi List[[l]] of shared eigenfunction coefficients
#' @param mu_q_nu_psi List[[s]][[l]] of specific eigenfunction coefficients
#' @param mu_q_zeta Shared factor scores, list[[s]][[l]]
#' @param Sigma_q_zeta Posterior covariance of zeta, list[[s]][[l]][[i]]
#' @param mu_q_xi Specific factor scores, list[[s]][[l]]
#' @param Sigma_q_xi Posterior covariance of xi, list[[s]][[l]][[i]]
#' @param mu_q_a Shared loadings a_{jl}, p x L_f matrix
#' @param mu_q_b_specific Specific loadings b_{sjl}, list[[s]] of p x L_s
#' @param mu_q_gamma_a Shared PPI, p x L_f matrix
#' @param mu_q_gamma_b Specific PPI, list[[s]] of p x L_s
#' @param S, n_s, p, L_f, L_s, M_f, M_s Dimensions; L_f/L_s = #factors;
#'   M_f[l] / M_s[[s]][l] = #basis functions per shared/specific factor
#'
#' @return List with orthonormalised eigenfunctions, scores, PVEs, and factor PPIs.
#'
#' @export
orthonormalise_multi <- function(C_g, time_g,
                                  mu_q_nu_mu,
                                  mu_q_nu_phi, mu_q_nu_psi,
                                  mu_q_zeta, Sigma_q_zeta,
                                  mu_q_xi, Sigma_q_xi,
                                  mu_q_a, mu_q_b_specific,
                                  mu_q_gamma_a, mu_q_gamma_b,
                                  S, n_s, p, L_f, L_s, M_f, M_s) {

  n_g <- length(time_g)
  one_N_all <- lapply(n_s, function(ns) rep(1, ns))

  # ---- Mean function reconstruction ----
  list_mu_hat <- vector("list", p)
  for (j in 1:p) {
    mu_funcs <- lapply(1:S, function(s) {
      sapply(mu_q_nu_mu[[s]], function(nu) as.vector(C_g %*% nu))  # n_g x p
    })
    list_mu_hat[[j]] <- mu_funcs[[1]][, j]
  }

  # ---- Shared factor orthonormalisation ----
  list_Phi_hat <- vector("list", L_f)
  list_Zeta_hat <- vector("list", L_f)
  list_Cov_zeta_hat <- vector("list", L_f)
  list_eigenvalues <- vector("list", L_f)
  list_cumulated_pve <- vector("list", L_f)

  for (l in seq_len(L_f)) {

    # Stack all individuals' scores for factor l
    Zeta_all <- do.call(rbind, lapply(1:S, function(s) mu_q_zeta[[s]][[l]]))
    N_total <- nrow(Zeta_all)

    # Evaluate eigenfunctions on dense grid
    Phi_g <- C_g %*% mu_q_nu_phi[[l]]  # n_g x M_f[l]

    # SVD of eigenfunctions
    phi_svd <- svd(Phi_g)
    U_orth <- phi_svd$u
    if (M_f[l] > 1) {
      D_diag <- diag(phi_svd$d, nrow = M_f[l], ncol = M_f[l])
      V_orth <- phi_svd$v
    } else {
      D_diag <- matrix(phi_svd$d)
      V_orth <- matrix(phi_svd$v, ncol = 1)
    }

    # Rotate scores
    Zeta_rot <- Zeta_all %*% V_orth %*% D_diag

    # Eigen-decompose rotated score covariance
    eigen_Zeta <- eigen(cov(Zeta_rot))
    Q_mat <- eigen_Zeta$vectors

    if (M_f[l] > 1) {
      Lambda_mat <- diag(abs(eigen_Zeta$values) + 1e-10, nrow = M_f[l], ncol = M_f[l])
      Lambda_inv <- diag(1 / (abs(eigen_Zeta$values) + 1e-10), nrow = M_f[l], ncol = M_f[l])
    } else {
      Lambda_mat <- matrix(abs(eigen_Zeta$values) + 1e-10)
      Lambda_inv <- matrix(1 / (abs(eigen_Zeta$values) + 1e-10))
    }

    S_mat <- Q_mat %*% sqrt(Lambda_mat)
    S_inv <- sqrt(Lambda_inv) %*% t(Q_mat)

    # Orthonormal eigenfunctions
    Phi_hat <- U_orth %*% S_mat
    Zeta_hat <- Zeta_rot %*% t(S_inv)

    # Normalise via trapezoidal integration
    norm_const <- rep(NA, M_f[l])
    for (m in 1:M_f[l]) {
      norm_const[m] <- sqrt(trapint(time_g, Phi_hat[, m]^2))
      if (norm_const[m] > 1e-10) {
        Phi_hat[, m] <- Phi_hat[, m] / norm_const[m]
        Zeta_hat[, m] <- Zeta_hat[, m] * norm_const[m]
      }
    }

    # Transform posterior covariance of zeta
    scale_mat <- if (M_f[l] > 1) diag(norm_const, nrow = M_f[l], ncol = M_f[l]) else matrix(norm_const)
    Cov_zeta_hat <- vector("list", N_total)
    idx <- 1
    for (s in 1:S) {
      for (i in 1:n_s[s]) {
        mat_transform <- S_inv %*% t(D_diag) %*% t(V_orth)
        Cov_zeta_dot <- tcrossprod(mat_transform %*% Sigma_q_zeta[[s]][[l]][[i]], mat_transform)
        Cov_zeta_hat[[idx]] <- tcrossprod(scale_mat %*% Cov_zeta_dot, scale_mat)
        idx <- idx + 1
      }
    }

    # PVE
    eigenvalues <- apply(Zeta_hat, 2, var)
    cumulated_pve <- cumsum(eigenvalues) / sum(eigenvalues) * 100

    list_Phi_hat[[l]] <- Phi_hat
    list_Zeta_hat[[l]] <- Zeta_hat
    list_Cov_zeta_hat[[l]] <- Cov_zeta_hat
    list_eigenvalues[[l]] <- eigenvalues
    list_cumulated_pve[[l]] <- cumulated_pve
  }

  # ---- Specific factor orthonormalisation ----
  list_Phi_hat_spec <- vector("list", S)
  list_Zeta_hat_spec <- vector("list", S)
  list_pve_spec <- vector("list", S)

  if (L_s > 0 && !is.null(mu_q_nu_psi)) {
    for (s in 1:S) {
      list_Phi_hat_spec[[s]] <- vector("list", L_s)
      list_Zeta_hat_spec[[s]] <- vector("list", L_s)
      list_pve_spec[[s]] <- vector("list", L_s)

      for (l in seq_len(L_s)) {
        Xi_s <- mu_q_xi[[s]][[l]]
        Phi_psi_g <- C_g %*% mu_q_nu_psi[[s]][[l]]  # n_g x M_s[[s]][l]

        psi_svd <- svd(Phi_psi_g)
        U_orth <- psi_svd$u
        if (M_s[[s]][l] > 1) {
          D_diag <- diag(psi_svd$d, nrow = M_s[[s]][l], ncol = M_s[[s]][l])
          V_orth <- psi_svd$v
        } else {
          D_diag <- matrix(psi_svd$d)
          V_orth <- matrix(psi_svd$v, ncol = 1)
        }

        Xi_rot <- Xi_s %*% V_orth %*% D_diag
        eigen_Xi <- eigen(cov(Xi_rot))
        Q_mat <- eigen_Xi$vectors

        if (M_s[[s]][l] > 1) {
          Lambda_mat <- diag(abs(eigen_Xi$values) + 1e-10, nrow = M_s[[s]][l], ncol = M_s[[s]][l])
          Lambda_inv <- diag(1 / (abs(eigen_Xi$values) + 1e-10), nrow = M_s[[s]][l], ncol = M_s[[s]][l])
        } else {
          Lambda_mat <- matrix(abs(eigen_Xi$values) + 1e-10)
          Lambda_inv <- matrix(1 / (abs(eigen_Xi$values) + 1e-10))
        }

        S_mat <- Q_mat %*% sqrt(Lambda_mat)
        S_inv <- sqrt(Lambda_inv) %*% t(Q_mat)

        Phi_psi_hat <- U_orth %*% S_mat
        Xi_hat <- Xi_rot %*% t(S_inv)

        norm_const <- rep(NA, M_s[[s]][l])
        for (m in 1:M_s[[s]][l]) {
          norm_const[m] <- sqrt(trapint(time_g, Phi_psi_hat[, m]^2))
          if (norm_const[m] > 1e-10) {
            Phi_psi_hat[, m] <- Phi_psi_hat[, m] / norm_const[m]
            Xi_hat[, m] <- Xi_hat[, m] * norm_const[m]
          }
        }

        list_Phi_hat_spec[[s]][[l]] <- Phi_psi_hat
        list_Zeta_hat_spec[[s]][[l]] <- Xi_hat
        eigenvalues <- apply(Xi_hat, 2, var)
        list_pve_spec[[s]][[l]] <- cumsum(eigenvalues) / sum(eigenvalues) * 100
      }
    }
  }

  # ---- Factor PPIs ----
  factor_ppi_shared <- 1 - exp(colSums(log1p(-mu_q_gamma_a)))
  factor_ppi_specific <- NULL
  if (L_s > 0 && !is.null(mu_q_gamma_b)) {
    factor_ppi_specific <- lapply(mu_q_gamma_b, function(gb) {
      1 - exp(colSums(log1p(-gb)))
    })
  }

  create_named_list(
    list_Phi_hat, list_Zeta_hat, list_Cov_zeta_hat,
    list_eigenvalues, list_cumulated_pve,
    list_Phi_hat_spec, list_Zeta_hat_spec, list_pve_spec,
    list_mu_hat,
    factor_ppi_shared, factor_ppi_specific
  )
}
