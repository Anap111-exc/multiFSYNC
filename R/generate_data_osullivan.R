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
project_to_osullivan <- function(f_vals, K, t_grid, ridge = 1e-6,
                                 int_knots = NULL) {
  if (is.null(int_knots)) {
    int_knots <- unname(quantile(
      t_grid, seq(0, 1, length = K)[-c(1, K)]))
  }
  Xg <- cbind(1, t_grid)
  Zg <- ZOSull(t_grid, range.x = c(0, 1),
               intKnots = int_knots)
  Cg <- cbind(Xg, Zg)
  # Cg is already in the O'Sullivan reparameterised coordinates [1, t, Z].
  # The corresponding curvature penalty is zero on the linear null space and
  # identity on the K nonlinear coordinates; OmegaOSull is in the original
  # B-spline coordinates and must not be mixed with Cg here.
  penalty <- diag(c(0, 0, rep(1, K)))
  P <- solve(crossprod(Cg) + ridge * penalty, t(Cg))
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
#' @noRd
#' @export
generate_explicit_eigenfunctions <- function(K, L_f, M_f, L_s = 0,
                                             M_s = NULL, S = 1, n_g = 200,
                                             int_knots = NULL,
                                             specific_time_heterogeneity = TRUE) {
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  t_g <- seq(0, 1, length.out = n_g)
  K_total <- K + 2
  if (is.null(int_knots)) {
    int_knots <- unname(quantile(
      t_g, seq(0, 1, length = K)[-c(1, K)]))
  }

  # Build spline Gram matrix G = C^T C / n_g for spline-space orthogonalization
  Xg <- cbind(1, t_g)
  Zg <- ZOSull(t_g, range.x = c(0, 1),
               intKnots = int_knots)
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
      raw_nu[, m] <- project_to_osullivan(f, K, t_g,
                                          int_knots = int_knots)
    }
    # Orthogonalize in spline space, then apply decreasing eigenvalues whose
    # sum is one, as required by the factor-process scale convention.
    ortho_nu <- orthogonalize_basis(raw_nu, G_mat)
    lambda <- (1 / seq_len(M_l)) / sum(1 / seq_len(M_l))
    for (m in seq_len(M_l)) ortho_nu[, m] <- ortho_nu[, m] * sqrt(lambda[m])
    nu_phi[[l]] <- ortho_nu
  }

  nu_psi <- NULL
  if (has_specific) {
    if (!is.list(M_s) || length(M_s) != S ||
        any(lengths(M_s) != L_s)) {
      stop("M_s must be a length-S list with length(M_s[[s]]) = L_s[s].")
    }
    nu_psi <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      nu_psi[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        raw_nu <- matrix(0, K_total, M_sl)
        study_offset <- if (s > 1L) sum(L_s[seq_len(s - 1L)]) else 0L
        base_freq <- L_f + l +
          if (specific_time_heterogeneity) study_offset else 0L
        for (m in 1:M_sl) {
          if (m %% 2 == 1) {
            f <- sin(2 * pi * base_freq * ((m + 1) / 2) * t_g)
          } else {
            f <- cos(2 * pi * base_freq * (m / 2) * t_g)
          }
          raw_nu[, m] <- project_to_osullivan(f, K, t_g,
                                              int_knots = int_knots)
        }
        ortho_nu <- orthogonalize_basis(raw_nu, G_mat)
        lambda <- (1 / seq_len(M_sl)) / sum(1 / seq_len(M_sl))
        for (m in seq_len(M_sl)) ortho_nu[, m] <- ortho_nu[, m] * sqrt(lambda[m])
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
generate_explicit_mean_functions <- function(K, p, S = 1, n_g = 200,
                                             int_knots = NULL,
                                             mean_structure = c("shared", "study_specific")) {
  mean_structure <- match.arg(mean_structure)
  t_g <- seq(0, 1, length.out = n_g)
  phases <- runif(p, 0, 2 * pi)
  base <- lapply(seq_len(p), function(j) {
    f <- sin(2 * pi * (t_g + phases[j])) + 0.5 * sin(4 * pi * t_g + 2 * phases[j])
    project_to_osullivan(f, K, t_g, int_knots = int_knots)
  })
  lapply(seq_len(S), function(s) {
    lapply(seq_len(p), function(j) {
      ans <- base[[j]]
      if (mean_structure == "study_specific") {
        study_shift <- (s - (S + 1) / 2) / max(S, 1)
        ans[1] <- ans[1] + study_shift * (0.4 + 0.1 * sin(phases[j]))
        ans[2] <- ans[2] + study_shift * 0.25
      }
      ans
    })
  })
}
# Design principles:
#   1. O'Sullivan basis for BOTH generation AND inference → zero basis mismatch
#   2. Decreasing component weights proportional to 1/m and normalized to sum 1
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
#   realized time-function blocks are rescaled to unit total integrated variance
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
    use_explicit = FALSE,
    mean_structure = c("shared", "study_specific"),
    specific_time_heterogeneity = TRUE,
    identified_loadings = FALSE
) {

  if (!is.null(seed)) set.seed(seed)
  mean_structure <- match.arg(mean_structure)

  # ---- Dimension validation ----
  stopifnot(length(n_s) == S)
  stopifnot(L_f >= 0)
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  d_use <- max(d, 0)

  if (is.null(M_f)) {
    M_f <- if (L_f > 0) rep(2, L_f) else integer(0)
  }
  stopifnot(length(M_f) == max(0, L_f))

  if (is.null(M_s)) {
    if (has_specific) {
      M_s <- lapply(seq_len(S), function(s) rep(2L, L_s[[s]]))
    } else {
      M_s <- replicate(S, integer(0), simplify = FALSE)
    }
  }
  stopifnot(is.list(M_s), length(M_s) == S)
  if (any(lengths(M_s) != L_s)) {
    stop("M_s must satisfy length(M_s[[s]]) = L_s[s] for every study.")
  }

  if (is.null(K)) {
    K <- max(7, min(floor(n_obs / 4), 40))
  }
  if (length(K) != 1L || !is.finite(K) || K != as.integer(K) || K < 2L) {
    stop("K must be one integer of at least 2.")
  }
  K <- as.integer(K)
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
  pooled_time <- sort(unique(unlist(time_obs, recursive = TRUE,
                                    use.names = FALSE)))
  int_knots <- unname(quantile(
    pooled_time, seq(0, 1, length = K)[-c(1, K)]))
  for (s in 1:S) {
    C[[s]] <- vector("list", n_s[s])
    for (i in 1:n_s[s]) {
      X <- cbind(1, time_obs[[s]][[i]])
      Z <- ZOSull(time_obs[[s]][[i]], range.x = c(0, 1),
                  intKnots = int_knots)
      C[[s]][[i]] <- cbind(X, Z)
    }
  }

  # ---- Step 3: Generate true parameters ----

  # 3a: Mean function coefficients (per study, per variable)
  nu_mu_true <- vector("list", S)
  if (use_explicit) {
    cat("  Using explicit sinusoidal mean functions...\n")
    nu_mu_true <- generate_explicit_mean_functions(
      K, p, S = S, int_knots = int_knots,
      mean_structure = mean_structure)
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
           L_s = L_s, M_s = M_s, S = S, n_g = 200,
           int_knots = int_knots,
           specific_time_heterogeneity = specific_time_heterogeneity)
    nu_phi_true <- ef$nu_phi
    nu_psi_true <- ef$nu_psi
    sigma_phi_used <- lapply(seq_len(L_f), function(l) {
      lambda <- (1 / seq_len(M_f[l])) / sum(1 / seq_len(M_f[l]))
      sqrt(lambda)
    })
    sigma_psi_used <- if (has_specific) lapply(seq_len(S), function(s) {
      lapply(seq_len(L_s[[s]]), function(l) {
        lambda <- (1 / seq_len(M_s[[s]][l])) /
          sum(1 / seq_len(M_s[[s]][l]))
        sqrt(lambda)
      })
    }) else NULL
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
    if (has_specific) {
    nu_psi_true <- vector("list", S)
    sigma_psi_used <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      nu_psi_true[[s]] <- vector("list", L_ss)
      sigma_psi_used[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
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

  factor_scale_before_normalization_shared <- rep(1, L_f)
  factor_scale_before_normalization_specific <- if (has_specific) {
    lapply(seq_len(S), function(s) rep(1, L_s[[s]]))
  } else NULL
  t_norm <- seq(0, 1, length.out = 1001L)
  C_norm <- cbind(
    1, t_norm,
    ZOSull(t_norm, range.x = c(0, 1), intKnots = int_knots))
  dt_norm <- diff(t_norm)
  w_norm <- c(dt_norm[1L] / 2,
              (dt_norm[-1L] + dt_norm[-length(dt_norm)]) / 2,
              dt_norm[length(dt_norm)] / 2)
  G_norm <- crossprod(C_norm, w_norm * C_norm)
  normalise_one <- function(Nu) {
    total_var <- sum(diag(crossprod(Nu, G_norm %*% Nu)))
    if (!is.finite(total_var) || total_var <= 0) {
      stop("Generated factor process has non-positive integrated variance.")
    }
    list(coef = Nu / sqrt(total_var), scale = sqrt(total_var))
  }
  for (l in seq_len(L_f)) {
    norm_l <- normalise_one(nu_phi_true[[l]])
    nu_phi_true[[l]] <- norm_l$coef
    sigma_phi_used[[l]] <- sigma_phi_used[[l]] / norm_l$scale
    factor_scale_before_normalization_shared[l] <- norm_l$scale
  }
  if (has_specific) {
    for (s in seq_len(S)) {
      for (l in seq_len(L_s[[s]])) {
        norm_sl <- normalise_one(nu_psi_true[[s]][[l]])
        nu_psi_true[[s]][[l]] <- norm_sl$coef
        sigma_psi_used[[s]][[l]] <- sigma_psi_used[[s]][[l]] /
          norm_sl$scale
        factor_scale_before_normalization_specific[[s]][l] <- norm_sl$scale
      }
    }
  }

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
  if (has_specific) {
    xi_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      xi_true[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        xi_true[[s]][[l]] <- matrix(
          rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
          nrow = n_s[s], ncol = M_sl)
      }
    }
  }

  # 3f: Loadings with spike-and-slab
  loading_draw <- if (identified_loadings) {
    generate_identified_loadings(p, L_f, L_s, S,
      sparse = bool_sparse_loadings, prop_sparse = prop_sparse)
  } else NULL
  a_true <- if (identified_loadings) loading_draw$a else matrix(0, p, L_f)
  gamma_a_true <- if (identified_loadings) loading_draw$gamma_a else matrix(0, p, L_f)
  for (l in if (identified_loadings) integer(0) else seq_len(L_f)) {
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

  b_true <- if (identified_loadings) loading_draw$b else NULL
  gamma_b_true <- if (identified_loadings) loading_draw$gamma_b else NULL
  if (has_specific && !identified_loadings) {
    b_true <- vector("list", S)
    gamma_b_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
      gamma_b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
      for (l in seq_len(L_ss)) {
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
    L_ss <- L_s[[s]]
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
      if (L_ss > 0L) {
        g_si_list <- vector("list", L_ss)
        for (l in seq_len(L_ss)) {
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

        if (L_ss > 0L) {
          for (l in seq_len(L_ss)) {
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
    factor_scale_before_normalization_shared,
    factor_scale_before_normalization_specific,
    zeta_true,
    xi_true,
    a_true,
    gamma_a_true,
    b_true,
    gamma_b_true,
    sigma2_eps_true,
    K, L_f, L_s, M_f, M_s,
    S, n_s, p, d_use, n_obs,
    sd_nu_mu, sd_nu_beta, sd_nu_phi, sd_nu_psi,
    int_knots, mean_structure, specific_time_heterogeneity,
    identified_loadings
  )
  true_params$L_s_by_study <- L_s
  true_params$L_s <- .compact_L_s(L_s)

  create_named_list(Y, Z, time_obs, C, true_params)
}
