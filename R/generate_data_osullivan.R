# =============================================================================
# generate_data_osullivan.R
#
# Data generation using O'Sullivan spline basis (SAME basis as inference)
# with decreasing-variance eigenfunction coefficients for natural FPCA ordering.
#
# Two modes:
#   1. Random (default): nu_phi ~ N(0, σ²/m), uses use_explicit=FALSE
#   2. Explicit: sinusoidal functions projected to O'Sullivan basis, use_explicit=TRUE
# =============================================================================

#' Project a function evaluated on a grid onto O'Sullivan spline basis
#'
#' @param f_vals Function values on t_grid (length n_g)
#' @param K Number of O'Sullivan internal knots
#' @param t_grid Dense grid of time points (length n_g)
#' @param ridge Ridge regularization for projection matrix (default 1e-6)
#' @return Coefficient vector of length K+2
#' @keywords internal
project_to_osullivan <- function(f_vals, K, t_grid, ridge = 1e-6) {
  Xg <- cbind(1, t_grid)
  Zg <- ZOSull(t_grid, range.x = c(0, 1),
               intKnots = quantile(t_grid, seq(0, 1, length = K)[-c(1, K)]))
  Cg <- cbind(Xg, Zg)
  Omega <- OmegaOSull(0, 1, intKnots = quantile(t_grid, seq(0, 1, length = K)[-c(1, K)]))
  P <- solve(crossprod(Cg) + ridge * Omega, t(Cg))
  as.vector(P %*% f_vals)
}

#' Generate explicit sinusoidal eigenfunctions projected to O'Sullivan basis
#'
#' For each factor l, generates M_f[l] orthonormal-like sinusoidal functions
#' with decreasing amplitude (1/√m), then projects them onto O'Sullivan basis.
#' Different factors use different base frequencies (factor l → frequency 2πl).
#'
#' @param K Number of O'Sullivan internal knots
#' @param L_f Number of shared factors
#' @param M_f Vector of FPCA components per factor
#' @param L_s Number of specific factors
#' @param M_s List of S vectors of components per specific factor
#' @param S Number of studies
#' @param n_g Dense grid size (default 200)
#' @return List with nu_phi (shared) and nu_psi (specific) as O'Sullivan coefficients
#' @export
generate_explicit_eigenfunctions <- function(K, L_f, M_f, L_s = 0, M_s = NULL, S = 1, n_g = 200) {
  t_g <- seq(0, 1, length.out = n_g)
  K_total <- K + 2

  # Build spline Gram matrix G = C^T C / n_g for spline-space orthogonalization
  Xg <- cbind(1, t_g)
  Zg <- ZOSull(t_g, range.x = c(0, 1),
               intKnots = quantile(t_g, seq(0, 1, length = K)[-c(1, K)]))
  Cg <- cbind(Xg, Zg)
  G_mat <- crossprod(Cg) / n_g

  # Orthogonalize coefficients in the spline metric: Nu^T G Nu = I
  orthogonalize_basis <- function(Nu_mat, G) {
    chol_G <- chol(G)
    transformed_Nu <- chol_G %*% Nu_mat
    svd_res <- svd(transformed_Nu)
    solve(chol_G, svd_res$u %*% t(svd_res$v))
  }

  nu_phi <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    raw_nu <- matrix(0, K_total, M_l)
    for (m in 1:M_l) {
      # Genuine orthogonal Fourier harmonics (sin/cos with integer multiples)
      if (m %% 2 == 1) {
        f <- sin(2 * pi * l * ((m + 1) / 2) * t_g)
      } else {
        f <- cos(2 * pi * l * (m / 2) * t_g)
      }
      raw_nu[, m] <- project_to_osullivan(f, K, t_g)
    }
    # Orthogonalize in spline space, then apply decreasing amplitude
    ortho_nu <- orthogonalize_basis(raw_nu, G_mat)
    for (m in 1:M_l) ortho_nu[, m] <- ortho_nu[, m] * (1 / sqrt(m))
    nu_phi[[l]] <- ortho_nu
  }

  nu_psi <- NULL
  if (L_s > 0) {
    nu_psi <- vector("list", S)
    for (s in 1:S) {
      nu_psi[[s]] <- vector("list", L_s)
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        raw_nu <- matrix(0, K_total, M_sl)
        base_freq <- L_f + l
        for (m in 1:M_sl) {
          if (m %% 2 == 1) {
            f <- sin(2 * pi * base_freq * ((m + 1) / 2) * t_g)
          } else {
            f <- cos(2 * pi * base_freq * (m / 2) * t_g)
          }
          raw_nu[, m] <- project_to_osullivan(f, K, t_g)
        }
        ortho_nu <- orthogonalize_basis(raw_nu, G_mat)
        for (m in 1:M_sl) ortho_nu[, m] <- ortho_nu[, m] * (1 / sqrt(m))
        nu_psi[[s]][[l]] <- ortho_nu
      }
    }
  }

  list(nu_phi = nu_phi, nu_psi = nu_psi)
}

#' Generate explicit sinusoidal mean functions projected to O'Sullivan basis
#'
#' μ_j(t) = sin(2π(t + φ_j)) + 0.5·sin(4πt + 2φ_j), projected to O'Sullivan basis
#'
#' @param K Number of O'Sullivan internal knots
#' @param p Number of variables
#' @param n_g Dense grid size
#' @return List of p coefficient vectors (K+2 each)
#' @keywords internal
generate_explicit_mean_functions <- function(K, p, n_g = 200) {
  t_g <- seq(0, 1, length.out = n_g)
  phases <- runif(p, 0, 2 * pi)
  lapply(1:p, function(j) {
    f <- sin(2 * pi * (t_g + phases[j])) + 0.5 * sin(4 * pi * t_g + 2 * phases[j])
    project_to_osullivan(f, K, t_g)
  })
}
# Design principles:
#   1. O'Sullivan basis for BOTH generation AND inference → zero basis mismatch
#   2. Decreasing ν_φ variance (σ²_m = 1/m²) → interpretable FPCA ordering
#   3. Unit score variance → matches Bayesian prior N(0,1)
#   4. Eigenvalue scale absorbed into ν_φ variance → consistent with model identifiability
#
# Model:
#   y_{sij}(t) = μ_{sj}(t) + Σ_r z_{sir} β_{jr}(t)
#              + Σ_l a_{jl} f^{(l)}_{si}(t) + Σ_l b_{sjl} g^{(l)}_{si}(t) + ε
#
#   f^{(l)}_{si}(t) = Σ_m ζ^{(l)}_{sim} φ^{(l)}_m(t),  φ^{(l)}_m(t) = C(t) ν^{(l)}_{φ,m}
#   g^{(l)}_{si}(t) = Σ_m ξ^{(sl)}_{sim} ψ^{(sl)}_m(t),  ψ^{(sl)}_m(t) = C(t) ν^{(sl)}_{ψ,m}
#
#   ν^{(l)}_{φ,m} ∼ N(0, (σ_φ/m)² I_K)     ← decreasing variance
#   ν^{(sl)}_{ψ,m} ∼ N(0, (σ_ψ/m)² I_K)     ← decreasing variance
#   ζ, ξ ∼ N(0, 1)                           ← unit variance (matches prior)
# =============================================================================

#' @export
simulate_multi_study_osullivan <- function(
    S = 2,
    n_s = c(30, 35),
    p = 10,
    d = 0,
    L_f = 3,
    L_s = 2,
    M_f = NULL,
    M_s = NULL,
    K = NULL,
    n_obs = 50,
    common_grid = TRUE,
    sigma_eps = 0.1,
    bool_sparse_loadings = TRUE,
    prop_sparse = 0.5,
    sd_nu_mu = 0.5,
    sd_nu_beta = 0.3,
    sd_nu_phi = 1.0,
    sd_nu_psi = 1.0,
    seed = NULL,
    use_explicit = FALSE
) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Dimension validation ----
  stopifnot(length(n_s) == S)
  stopifnot(L_f >= 0)
  stopifnot(L_s >= 0)
  d_use <- max(d, 0)

  if (is.null(M_f)) {
    M_f <- if (L_f > 0) rep(2, L_f) else integer(0)
  }
  stopifnot(length(M_f) == max(0, L_f))

  if (is.null(M_s)) {
    if (L_s > 0) {
      M_s <- replicate(S, rep(2, L_s), simplify = FALSE)
    } else {
      M_s <- replicate(S, integer(0), simplify = FALSE)
    }
  }
  stopifnot(is.list(M_s), length(M_s) == S)

  if (is.null(K)) {
    K <- max(7, min(floor(n_obs / 4), 40))
  }
  K_total <- K + 2

  # ---- Step 1: Observation time grid ----
  time_obs <- vector("list", S)
  if (common_grid) {
    t_grid <- seq(0, 1, length.out = n_obs)
    for (s in 1:S) {
      time_obs[[s]] <- lapply(1:n_s[s], function(i) t_grid)
    }
  } else {
    for (s in 1:S) {
      time_obs[[s]] <- lapply(1:n_s[s], function(i) sort(runif(n_obs, 0, 1)))
    }
  }

  # ---- Step 2: O'Sullivan spline design matrix C ----
  # SAME construction as used in inference (get_grid_objects)
  C <- vector("list", S)
  for (s in 1:S) {
    C[[s]] <- vector("list", n_s[s])
    for (i in 1:n_s[s]) {
      X <- cbind(1, time_obs[[s]][[i]])
      Z <- ZOSull(time_obs[[s]][[i]], range.x = c(0, 1),
                  intKnots = quantile(time_obs[[s]][[i]],
                                      seq(0, 1, length = K)[-c(1, K)]))
      C[[s]][[i]] <- cbind(X, Z)
    }
  }

  # ---- Step 3: Generate true parameters ----

  # 3a: Mean function coefficients (per study, per variable)
  nu_mu_true <- vector("list", S)
  if (use_explicit) {
    cat("  Using explicit sinusoidal mean functions...\n")
    mu_explicit <- generate_explicit_mean_functions(K, p)
    for (s in 1:S) {
      nu_mu_true[[s]] <- mu_explicit
    }
  } else {
    for (s in 1:S) {
      nu_mu_true[[s]] <- lapply(1:p, function(j) {
        rnorm(K_total, mean = 0, sd = sd_nu_mu)
      })
    }
  }

  # 3b: Regression coefficients (shared across studies)
  nu_beta_true <- NULL
  if (d_use > 0) {
    nu_beta_true <- lapply(1:p, function(j) {
      lapply(1:d_use, function(r) {
        rnorm(K_total, mean = 0, sd = sd_nu_beta)
      })
    })
  }

  # 3c: Shared eigenfunction coefficients
  if (use_explicit) {
    cat("  Using explicit sinusoidal eigenfunctions...\n")
    ef <- generate_explicit_eigenfunctions(K = K, L_f = L_f, M_f = M_f,
           L_s = L_s, M_s = M_s, S = S, n_g = 200)
    nu_phi_true <- ef$nu_phi
    nu_psi_true <- ef$nu_psi
    sigma_phi_used <- lapply(1:L_f, function(l) rep(1/sqrt(1:M_f[l])))
    sigma_psi_used <- if (L_s > 0) lapply(1:S, function(s) lapply(1:L_s, function(l) rep(1/sqrt(1:M_s[[s]][l])))) else NULL
  } else {
    nu_phi_true <- vector("list", L_f)
    sigma_phi_used <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      nu_mat <- matrix(0, K_total, M_l)
      sd_vec <- rep(NA, M_l)
      for (m in 1:M_l) {
        sd_m <- sd_nu_phi / sqrt(m)
        nu_mat[, m] <- rnorm(K_total, mean = 0, sd = sd_m)
        sd_vec[m] <- sd_m
      }
      nu_phi_true[[l]] <- nu_mat
      sigma_phi_used[[l]] <- sd_vec
    }
  }

  # 3d: Study-specific eigenfunction coefficients (random mode only)
  if (!use_explicit) {
    nu_psi_true <- NULL
    sigma_psi_used <- NULL
    if (L_s > 0) {
    nu_psi_true <- vector("list", S)
    sigma_psi_used <- vector("list", S)
    for (s in 1:S) {
      nu_psi_true[[s]] <- vector("list", L_s)
      sigma_psi_used[[s]] <- vector("list", L_s)
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        nu_mat <- matrix(0, K_total, M_sl)
        sd_vec <- rep(NA, M_sl)
        for (m in 1:M_sl) {
          sd_m <- sd_nu_psi / sqrt(m)
          nu_mat[, m] <- rnorm(K_total, mean = 0, sd = sd_m)
          sd_vec[m] <- sd_m
        }
        nu_psi_true[[s]][[l]] <- nu_mat
        sigma_psi_used[[s]][[l]] <- sd_vec
      }
    }
    }
  } # end if(!use_explicit)

  # 3e: FPCA scores — unit variance (matching Bayesian prior)
  zeta_true <- vector("list", S)
  for (s in 1:S) {
    zeta_true[[s]] <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      zeta_true[[s]][[l]] <- matrix(
        rnorm(n_s[s] * M_l, mean = 0, sd = 1),
        nrow = n_s[s], ncol = M_l)
    }
  }

  xi_true <- NULL
  if (L_s > 0) {
    xi_true <- vector("list", S)
    for (s in 1:S) {
      xi_true[[s]] <- vector("list", L_s)
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        xi_true[[s]][[l]] <- matrix(
          rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
          nrow = n_s[s], ncol = M_sl)
      }
    }
  }

  # 3f: Loadings with spike-and-slab
  a_true <- matrix(0, nrow = p, ncol = L_f)
  gamma_a_true <- matrix(0, nrow = p, ncol = L_f)
  for (l in seq_len(L_f)) {
    if (bool_sparse_loadings) {
      n_active <- max(1, round(p * (1 - prop_sparse)))
      active_idx <- sample(1:p, n_active)
      val_mean <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      a_true[active_idx, l] <- rnorm(n_active, mean = val_mean, sd = 0.2)
      gamma_a_true[active_idx, l] <- 1
    } else {
      a_true[, l] <- rnorm(p, mean = 0, sd = 0.5)
      gamma_a_true[, l] <- 1
    }
  }

  b_true <- NULL
  gamma_b_true <- NULL
  if (L_s > 0) {
    b_true <- vector("list", S)
    gamma_b_true <- vector("list", S)
    for (s in 1:S) {
      b_true[[s]] <- matrix(0, nrow = p, ncol = L_s)
      gamma_b_true[[s]] <- matrix(0, nrow = p, ncol = L_s)
      for (l in seq_len(L_s)) {
        if (bool_sparse_loadings) {
          n_active <- max(1, round(p * (1 - prop_sparse)))
          active_idx <- sample(1:p, n_active)
          val_mean <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
          b_true[[s]][active_idx, l] <- rnorm(n_active, mean = val_mean, sd = 0.2)
          gamma_b_true[[s]][active_idx, l] <- 1
        } else {
          b_true[[s]][, l] <- rnorm(p, mean = 0, sd = 0.5)
          gamma_b_true[[s]][, l] <- 1
        }
      }
    }
  }

  # 3g: Covariates
  Z <- NULL
  if (d_use > 0) {
    Z <- vector("list", S)
    for (s in 1:S) {
      Z[[s]] <- matrix(rnorm(n_s[s] * d_use, mean = 0, sd = 1),
                        nrow = n_s[s], ncol = d_use)
    }
  }

  # 3h: Measurement error
  sigma2_eps_true <- matrix(sigma_eps^2, nrow = S, ncol = p)

  # ---- Step 4: Generate observed data Y ----
  Y <- vector("list", S)

  for (s in 1:S) {
    Y[[s]] <- vector("list", n_s[s])
    for (i in 1:n_s[s]) {
      Y[[s]][[i]] <- vector("list", p)
      C_si <- C[[s]][[i]]

      # Shared factor processes: f^{(l)}_{si} = C_si * nu_phi_l * zeta_{si,l}^T
      f_si_list <- vector("list", L_f)
      for (l in seq_len(L_f)) {
        f_si_list[[l]] <- as.vector(C_si %*% nu_phi_true[[l]] %*% zeta_true[[s]][[l]][i, ])
      }

      # Specific factor processes
      g_si_list <- NULL
      if (L_s > 0) {
        g_si_list <- vector("list", L_s)
        for (l in seq_len(L_s)) {
          g_si_list[[l]] <- as.vector(
            C_si %*% nu_psi_true[[s]][[l]] %*% xi_true[[s]][[l]][i, ])
        }
      }

      # Mean functions: μ_{sj} = C_si * nu_mu_{sj}
      mu_si_list <- lapply(1:p, function(j) as.vector(C_si %*% nu_mu_true[[s]][[j]]))

      # Beta contributions
      beta_si_list <- NULL
      if (d_use > 0) {
        beta_si_list <- lapply(1:p, function(j) {
          lapply(1:d_use, function(r) as.vector(C_si %*% nu_beta_true[[j]][[r]]))
        })
      }

      for (j in 1:p) {
        y_sij <- mu_si_list[[j]]

        if (d_use > 0) {
          for (r in 1:d_use) {
            y_sij <- y_sij + Z[[s]][i, r] * beta_si_list[[j]][[r]]
          }
        }

        for (l in seq_len(L_f)) {
          y_sij <- y_sij + a_true[j, l] * f_si_list[[l]]
        }

        if (L_s > 0) {
          for (l in seq_len(L_s)) {
            y_sij <- y_sij + b_true[[s]][j, l] * g_si_list[[l]]
          }
        }

        y_sij <- y_sij + rnorm(length(y_sij), mean = 0, sd = sqrt(sigma2_eps_true[s, j]))
        Y[[s]][[i]][[j]] <- y_sij
      }
    }
  }

  # ---- Step 5: Assemble true parameter list ----
  true_params <- create_named_list(
    nu_mu_true,
    nu_beta_true,
    nu_phi_true,
    nu_psi_true,
    sigma_phi_used,
    sigma_psi_used,
    zeta_true,
    xi_true,
    a_true,
    gamma_a_true,
    b_true,
    gamma_b_true,
    sigma2_eps_true,
    K, L_f, L_s, M_f, M_s,
    S, n_s, p, d_use, n_obs,
    sd_nu_mu, sd_nu_beta, sd_nu_phi, sd_nu_psi
  )

  create_named_list(Y, Z, time_obs, C, true_params)
}
