# =============================================================================
# generate_data_structured.R
#
# Data generation using EXPLICIT functional forms (following Jaoua et al. 2026 §5.1)
# Difference from generate_data.R: eigenfunctions are orthonormalised B-splines,
# mean functions are periodic+phase-shifted, scores have decreasing variance.
#
# Model:
#   y_{sij}(t) = mu_{sj}(t) + Σ_r z_{sir} * beta_{jr}(t)
#             + Σ_l a_{jl} * f^{(l)}_{si}(t)
#             + Σ_l b_{sjl} * g^{(l)}_{si}(t) + eps,  eps ~ N(0, sigma^2_eps I)
#
# with
#   f^{(l)}_{si}(t) = Σ_{m=1}^{M_f[l]} zeta^{(l)}_{sim} * phi^{(l)}_m(t)
#   g^{(l)}_{si}(t) = Σ_{m=1}^{M_s[[s]][l]} xi^{(sl)}_{sim} * psi^{(sl)}_m(t)
# =============================================================================

#' Construct L2-orthonormal B-spline basis functions on [0,1]
#'
#' @param M Number of basis functions
#' @param n_g Dense grid size
#' @param degree B-spline degree (default 3 = cubic)
#' @return List: phi_matrix (n_g x M), t_grid (n_g x 1)
#' @keywords internal
construct_orthonormal_eigenfunctions <- function(M, n_g = 200, degree = 3) {
  t_grid <- seq(0, 1, length.out = n_g)

  # Create B-spline basis with slightly more df to allow flexibility
  df_bs <- max(M + degree, M * 2)
  bs_mat <- splines::bs(t_grid, df = df_bs, degree = degree, intercept = TRUE)

  # Select first M columns and orthonormalize via Gram-Schmidt
  phi <- matrix(0, n_g, M)
  for (m in 1:M) {
    v <- bs_mat[, m]
    for (k in seq_len(m - 1)) {
      ip <- trapint(t_grid, v * phi[, k])
      v <- v - ip * phi[, k]
    }
    norm_v <- sqrt(trapint(t_grid, v^2))
    if (norm_v > 1e-10) {
      phi[, m] <- v / norm_v
    } else {
      phi[, m] <- v
    }
  }
  list(phi = phi, t_grid = t_grid)
}

#' Evaluate an eigenfunction at arbitrary time points by interpolation
#'
#' @param t_pts Time points to evaluate
#' @param phi_dense Dense-grid eigenfunction vector
#' @param t_grid Dense grid
#' @return Vector of interpolated values
#' @keywords internal
eval_eigenfun <- function(t_pts, phi_dense, t_grid) {
  stats::approx(t_grid, phi_dense, xout = t_pts, rule = 2)$y
}

#' Generate periodic mean function with phase shift
#'
#' mu_j(t) = sin(2*pi*t + phase_j) + 0.5 * sin(4*pi*t + 2*phase_j)
#'
#' @param t_pts Time points
#' @param phase_j Phase shift for variable j
#' @param amplitude Scaling factor
#' @return Vector of function values
#' @keywords internal
mean_function_periodic <- function(t_pts, phase_j, amplitude = 1) {
  amplitude * (sin(2 * pi * t_pts + phase_j) +
               0.5 * sin(4 * pi * t_pts + 2 * phase_j))
}

#' Generate regression coefficient function
#'
#' beta_{jr}(t) = amplitude * sin(pi * r * t + phase)
#'
#' @param t_pts Time points
#' @param r Covariate index
#' @param phase Phase shift (random per variable-covariate pair)
#' @param amplitude Scaling factor
#' @return Vector of function values
#' @keywords internal
beta_function_periodic <- function(t_pts, r, phase, amplitude = 0.3) {
  amplitude * sin(pi * r * t_pts + phase)
}

#' Simulate multi-study functional factor model data with structured functions
#'
#' @inheritParams simulate_multi_study_data
#' @param score_var_decay If TRUE, scores have decreasing variance 1/m^2.
#'   If FALSE, unit variance (matches Bayesian prior).
#' @param n_dense Number of grid points for dense eigenfunction construction
#' @param mean_amp Amplitude of mean functions
#' @param beta_amp Amplitude of regression coefficient functions
#' @param bs_degree Degree of B-splines for eigenfunction construction
#'
#' @return A list with elements Y, Z, time_obs, C, true_params.
#'
#' @export
simulate_multi_study_structured <- function(
    S = 2,
    n_s = c(20, 30),
    p = 10,
    d = 2,
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
    score_var_decay = TRUE,
    n_dense = 200,
    mean_amp = 1.0,
    beta_amp = 0.3,
    bs_degree = 3,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Dimension validation ----
  stopifnot(length(n_s) == S)
  stopifnot(L_f >= 0)
  stopifnot(L_s >= 0)
  if (d == 0) d <- 0
  d_use <- max(d, 0)

  if (is.null(M_f)) {
    M_f <- if (L_f > 0) rep(3, L_f) else integer(0)
  }
  stopifnot(length(M_f) == max(0, L_f))

  if (is.null(M_s)) {
    if (L_s > 0) {
      M_s <- replicate(S, rep(3, L_s), simplify = FALSE)
    } else {
      M_s <- replicate(S, integer(0), simplify = FALSE)
    }
  }
  stopifnot(is.list(M_s), length(M_s) == S)
  for (s in 1:S) {
    stopifnot(length(M_s[[s]]) == L_s)
  }

  if (is.null(K)) {
    K <- max(7, min(floor(n_obs / 4), 40))
  }
  K_total <- K + 2

  # ---- Step 1: Construct observation time grid ----
  time_obs <- vector("list", S)
  if (common_grid) {
    t_grid_obs <- seq(0, 1, length.out = n_obs)
    for (s in 1:S) {
      time_obs[[s]] <- lapply(1:n_s[s], function(i) t_grid_obs)
    }
  } else {
    for (s in 1:S) {
      time_obs[[s]] <- lapply(1:n_s[s], function(i) sort(runif(n_obs, 0, 1)))
    }
  }

  # ---- Step 2: Construct spline design matrix C (O'Sullivan) ----
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

  # ---- Step 3: Generate true structured functions ----

  max_M_f <- if (L_f > 0) max(M_f) else 0
  max_M_s <- if (L_s > 0) max(sapply(1:S, function(s)
    if (L_s > 0) max(M_s[[s]]) else 0)) else 0
  max_M <- max(max_M_f, max_M_s)

  # 3a: Construct orthonormal eigenfunctions on dense grid
  if (max_M > 0) {
    ef_obj <- construct_orthonormal_eigenfunctions(max_M, n_dense, bs_degree)
    phi_dense_all <- ef_obj$phi
    t_grid_dense <- ef_obj$t_grid
  }

  # 3b: Assign eigenfunctions to each shared factor
  phi_dense_list <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    phi_dense_list[[l]] <- phi_dense_all[, 1:M_l, drop = FALSE]
  }

  # 3c: Assign eigenfunctions to each study-specific factor
  psi_dense_list <- NULL
  if (L_s > 0) {
    psi_dense_list <- vector("list", S)
    for (s in 1:S) {
      psi_dense_list[[s]] <- vector("list", L_s)
      used_cols <- 0
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        # Use different columns for each study to ensure distinguishability
        start_col <- ((s - 1) * max(L_s, 1) + l)
        cols <- ((start_col - 1) %% max_M) + 1
        cols <- cols:(cols + M_sl - 1)
        cols <- ((cols - 1) %% max_M) + 1
        psi_dense_list[[s]][[l]] <- phi_dense_all[, cols, drop = FALSE]
      }
    }
  }

  # 3d: Generate phase shifts for mean functions and beta functions
  phase_mu <- runif(p * S, 0, 2 * pi)
  phase_mu_mat <- matrix(phase_mu, nrow = S, ncol = p)

  phase_beta <- NULL
  if (d_use > 0) {
    n_phase_beta <- p * d_use
    phase_beta <- runif(n_phase_beta, 0, 2 * pi)
  }

  # 3e: Generate FPCA scores
  zeta_true <- vector("list", S)
  for (s in 1:S) {
    zeta_true[[s]] <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      sd_vec <- if (score_var_decay) 1 / sqrt(1:M_l) else rep(1, M_l)
      zeta_true[[s]][[l]] <- matrix(
        rnorm(n_s[s] * M_l, mean = 0, sd = rep(sd_vec, each = n_s[s])),
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
        sd_vec <- if (score_var_decay) 1 / sqrt(1:M_sl) else rep(1, M_sl)
        xi_true[[s]][[l]] <- matrix(
          rnorm(n_s[s] * M_sl, mean = 0, sd = rep(sd_vec, each = n_s[s])),
          nrow = n_s[s], ncol = M_sl)
      }
    }
  }

  # 3f: Generate loadings with spike-and-slab
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
  mu_true_values <- vector("list", S)
  f_true_values <- vector("list", S)
  g_true_values <- vector("list", S)
  beta_true_values <- vector("list", S)

  for (s in 1:S) {
    Y[[s]] <- vector("list", n_s[s])
    mu_true_values[[s]] <- vector("list", n_s[s])
    f_true_values[[s]] <- vector("list", n_s[s])
    g_true_values[[s]] <- vector("list", n_s[s])
    beta_true_values[[s]] <- vector("list", n_s[s])

    for (i in 1:n_s[s]) {
      Y[[s]][[i]] <- vector("list", p)
      t_i <- time_obs[[s]][[i]]
      n_i <- length(t_i)

      mu_si <- matrix(0, n_i, p)
      f_si <- matrix(0, n_i, L_f)
      g_si <- if (L_s > 0) matrix(0, n_i, L_s) else NULL
      beta_si_contrib <- rep(0, n_i)

      # Mean functions
      for (j in 1:p) {
        mu_si[, j] <- mean_function_periodic(t_i, phase_mu_mat[s, j], mean_amp)
      }

      # Shared factor processes
      for (l in seq_len(L_f)) {
        M_l <- M_f[l]
        zeta_i <- zeta_true[[s]][[l]][i, ]
        phi_val <- matrix(0, n_i, M_l)
        for (m in 1:M_l) {
          phi_val[, m] <- eval_eigenfun(t_i, phi_dense_list[[l]][, m], t_grid_dense)
        }
        f_si[, l] <- as.vector(phi_val %*% zeta_i)
      }

      # Specific factor processes
      if (L_s > 0) {
        for (l in seq_len(L_s)) {
          M_sl <- M_s[[s]][l]
          xi_i <- xi_true[[s]][[l]][i, ]
          psi_val <- matrix(0, n_i, M_sl)
          for (m in 1:M_sl) {
            psi_val[, m] <- eval_eigenfun(t_i, psi_dense_list[[s]][[l]][, m],
                                          t_grid_dense)
          }
          g_si[, l] <- as.vector(psi_val %*% xi_i)
        }
      }

      # Per-variable response
      for (j in 1:p) {
        y_sij <- mu_si[, j]

        # Covariate contributions
        if (d_use > 0) {
          for (r in 1:d_use) {
            phase_idx <- (j - 1) * d_use + r
            beta_val <- beta_function_periodic(t_i, r, phase_beta[phase_idx], beta_amp)
            y_sij <- y_sij + Z[[s]][i, r] * beta_val
            beta_si_contrib <- beta_si_contrib + Z[[s]][i, r] * beta_val
          }
        }

        # Shared factor contributions
        for (l in seq_len(L_f)) {
          y_sij <- y_sij + a_true[j, l] * f_si[, l]
        }

        # Specific factor contributions
        if (L_s > 0) {
          for (l in seq_len(L_s)) {
            y_sij <- y_sij + b_true[[s]][j, l] * g_si[, l]
          }
        }

        # Noise
        y_sij <- y_sij + rnorm(n_i, mean = 0, sd = sqrt(sigma2_eps_true[s, j]))
        Y[[s]][[i]][[j]] <- y_sij
      }

      mu_true_values[[s]][[i]] <- mu_si
      f_true_values[[s]][[i]] <- f_si
      if (L_s > 0) g_true_values[[s]][[i]] <- g_si
      beta_true_values[[s]][[i]] <- beta_si_contrib
    }
  }

  # ---- Step 5: Assemble true parameter list ----
  true_params <- create_named_list(
    a_true,
    gamma_a_true,
    b_true,
    gamma_b_true,
    zeta_true,
    xi_true,
    phi_dense_list,
    psi_dense_list,
    t_grid_dense,
    phase_mu_mat,
    phase_beta,
    sigma2_eps_true,
    mu_true_values,
    f_true_values,
    g_true_values,
    beta_true_values,
    K, L_f, L_s, M_f, M_s,
    S, n_s, p, d_use, n_obs,
    score_var_decay
  )

  create_named_list(Y, Z, time_obs, C, true_params)
}
