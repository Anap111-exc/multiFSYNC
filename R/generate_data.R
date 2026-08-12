# =============================================================================
# generate_data.R — Multi-study functional factor model data simulation
#
# Generates synthetic data following the model in derivations.md §1.2-1.3:
#   y_{sij} = mu_{sj} + sum_r z_{sir} * beta_{jr}
#           + sum_l a_{jl} * f^{(l)}_{si}
#           + sum_l b_{sjl} * g^{(l)}_{si}
#           + epsilon,  epsilon ~ N(0, sigma^2_{eps,sj} I_n)
#
# Symbol convention (aligned with 创新点.pdf):
#   L_f = number of shared factors (scalar)
#   L_s = number of study-specific factors (scalar or one count per study)
#   M_f = vector length L_f: M_f[l] = FPCA truncation for shared factor l
#   M_s = list of S vectors: M_s[[s]][l] = FPCA truncation for specific factor l in study s
#
# Indexing convention: study [[s]] is always the outermost dimension.
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Simulate multi-study functional factor model data
#'
#' @param S Number of studies (>= 1).
#' @param n_s Vector of length S, number of individuals per study.
#' @param p Number of functional variables.
#' @param d Number of scalar covariates (0 for no covariates).
#' @param L_f Number of shared latent factors (scalar).
#' @param L_s Number of study-specific latent factors. A scalar is recycled to
#'   all studies; a length-S vector allows different counts by study.
#' @param M_f Vector of length L_f, number of FPCA components per shared factor.
#'   If NULL, defaults to rep(2, L_f).
#' @param M_s List of S vectors, with length M_s[[s]] = L_s[s].
#'   M_s[[s]][l] is the FPCA truncation for factor l in study s.
#' @param K Number of O'Sullivan spline basis functions (determines design matrix
#'   size K+2). If NULL, defaults to max(7, floor(min(n_obs/4, 40))).
#' @param n_obs Number of observation grid points per curve.
#' @param common_grid Logical. If TRUE, all individuals share the same grid.
#' @param sigma_eps Standard deviation of measurement error (scalar).
#' @param bool_sparse_loadings Logical. If TRUE, use sparse spike-and-slab
#'   pattern for loadings with p_as_spike proportion at exactly zero.
#' @param prop_sparse Proportion of loadings set to zero when
#'   bool_sparse_loadings = TRUE.
#' @param identified_loadings If TRUE, construct loadings satisfying the
#'   population Gram diagonalisation, strict norm ordering and full-rank
#'   conditions used in the thesis simulations.
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with elements:
#'   \item{Y}{List S of lists n_s[s] of lists p, each an n_obs-length vector.}
#'   \item{Z}{List of S matrices (n_s[s] x d), or NULL if d=0.}
#'   \item{time_obs}{List S of lists n_s[s] of time point vectors.}
#'   \item{C}{List S of lists n_s[s] of design matrices (n_obs x (K+2)).}
#'   \item{true_params}{List of true parameter values.}
#'
#' @noRd
#' @export
simulate_multi_study_data <- function(
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
    prop_sparse = 0.7,
    identified_loadings = FALSE,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Dimension validation ----
  stopifnot(length(n_s) == S)
  stopifnot(L_f >= 0)
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  if (d == 0) d <- 0
  d_use <- max(d, 0)

  if (is.null(M_f)) {
    M_f <- if (L_f > 0) rep(2, L_f) else integer(0)
  }
  stopifnot(length(M_f) == max(0, L_f))
  stopifnot(all(M_f >= 0))

  if (is.null(M_s)) {
    if (has_specific) {
      M_s <- lapply(seq_len(S), function(s) rep(2L, L_s[[s]]))
    } else {
      M_s <- replicate(S, integer(0), simplify = FALSE)
    }
  }
  stopifnot(is.list(M_s), length(M_s) == S)
  for (s in 1:S) {
    stopifnot(length(M_s[[s]]) == L_s[[s]])
  }

  if (is.null(K)) {
    K <- max(7, min(floor(n_obs / 4), 40))
  }

  K_total <- K + 2

  # ---- Step 1: Construct observation time grid ----
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

  # ---- Step 2: Construct spline design matrix C ----
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

  # --- 3a: Spline coefficients ---

  # Mean function coefficients nu_{mu,sj}: (K+2) x 1
  nu_mu_true <- vector("list", S)
  for (s in 1:S) {
    nu_mu_true[[s]] <- lapply(1:p, function(j) {
      rnorm(K_total, mean = 0, sd = 0.5)
    })
  }

  # Regression coefficients nu_{beta,jr}: (K+2) x 1 (shared across studies)
  nu_beta_true <- NULL
  if (d_use > 0) {
    nu_beta_true <- lapply(1:p, function(j) {
      lapply(1:d_use, function(r) {
        rnorm(K_total, mean = 0, sd = 0.3)
      })
    })
  }

  # --- 3b: Shared eigenfunction coefficients nu_{phi,ml}: (K+2) x 1 ---
  nu_phi_true <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    nu_phi_true[[l]] <- matrix(rnorm(K_total * M_l, mean = 0, sd = 0.5),
                                nrow = K_total, ncol = M_l)
  }

  # Study-specific eigenfunction coefficients nu_{psi,slm}: (K+2) x 1
  nu_psi_true <- NULL
  if (has_specific) {
    nu_psi_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      nu_psi_true[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        nu_psi_true[[s]][[l]] <- matrix(rnorm(K_total * M_sl, mean = 0, sd = 0.5),
                                          nrow = K_total, ncol = M_sl)
      }
    }
  }

  # --- 3c: FPCA scores ---

  # Shared factor scores zeta^{(l)}_{sim} ~ N(0, 1)
  zeta_true <- vector("list", S)
  for (s in 1:S) {
    zeta_true[[s]] <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      zeta_true[[s]][[l]] <- matrix(rnorm(n_s[s] * M_l, mean = 0, sd = 1),
                                     nrow = n_s[s], ncol = M_l)
    }
  }

  # Specific factor scores xi^{(sl)}_{sim} ~ N(0, 1)
  xi_true <- NULL
  if (has_specific) {
    xi_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      xi_true[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        xi_true[[s]][[l]] <- matrix(rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
                                     nrow = n_s[s], ncol = M_sl)
      }
    }
  }

  # --- 3d: Loadings with optional sparsity ---

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

  # --- 3e: Covariates ---
  Z <- NULL
  if (d_use > 0) {
    Z <- vector("list", S)
    for (s in 1:S) {
      Z[[s]] <- matrix(rnorm(n_s[s] * d_use, mean = 0, sd = 1),
                        nrow = n_s[s], ncol = d_use)
    }
  }

  # --- 3f: Measurement error variance ---
  sigma2_eps_true <- matrix(sigma_eps^2, nrow = S, ncol = p)

  # ---- Step 4: Generate observed data Y ----

  Y <- vector("list", S)

  for (s in 1:S) {
    L_ss <- L_s[[s]]
    Y[[s]] <- vector("list", n_s[s])
    for (i in 1:n_s[s]) {
      Y[[s]][[i]] <- vector("list", p)

      # Precompute shared factor process for this (s,i)
      f_si_list <- vector("list", L_f)
      for (l in seq_len(L_f)) {
        f_si_list[[l]] <- as.vector(
          C[[s]][[i]] %*% nu_phi_true[[l]] %*% zeta_true[[s]][[l]][i, ]
        )
      }

      # Precompute specific factor process for this (s,i) if L_s > 0
      g_si_list <- NULL
      if (L_ss > 0L) {
        g_si_list <- vector("list", L_ss)
        for (l in seq_len(L_ss)) {
          g_si_list[[l]] <- as.vector(
            C[[s]][[i]] %*% nu_psi_true[[s]][[l]] %*% xi_true[[s]][[l]][i, ]
          )
        }
      }

      # Precompute mean for each variable
      mu_si_list <- vector("list", p)
      for (j in 1:p) {
        mu_si_list[[j]] <- as.vector(C[[s]][[i]] %*% nu_mu_true[[s]][[j]])
      }

      # Precompute beta contribution for each (j,r)
      beta_si_list <- NULL
      if (d_use > 0) {
        beta_si_list <- vector("list", p)
        for (j in 1:p) {
          beta_si_list[[j]] <- vector("list", d_use)
          for (r in 1:d_use) {
            beta_si_list[[j]][[r]] <- as.vector(
              C[[s]][[i]] %*% nu_beta_true[[j]][[r]]
            )
          }
        }
      }

      for (j in 1:p) {
        # sum_{r} z_{sir} * beta_{jr}
        y_sij <- mu_si_list[[j]]
        if (d_use > 0) {
          for (r in 1:d_use) {
            y_sij <- y_sij + Z[[s]][i, r] * beta_si_list[[j]][[r]]
          }
        }

        # sum_{l} a_{jl} * f^{(l)}_{si}
        for (l in seq_len(L_f)) {
          y_sij <- y_sij + a_true[j, l] * f_si_list[[l]]
        }

        # sum_{l} b_{sjl} * g^{(l)}_{si}
        if (L_ss > 0L) {
          for (l in seq_len(L_ss)) {
            y_sij <- y_sij + b_true[[s]][j, l] * g_si_list[[l]]
          }
        }

        # Add measurement error
        y_sij <- y_sij + rnorm(n_obs, mean = 0, sd = sqrt(sigma2_eps_true[s, j]))

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
    zeta_true,
    xi_true,
    a_true,
    gamma_a_true,
    b_true,
    gamma_b_true,
    sigma2_eps_true,
    K, L_f, L_s, M_f, M_s,
    S, n_s, p, d_use, n_obs, int_knots, identified_loadings
  )
  true_params$L_s_by_study <- L_s
  true_params$L_s <- .compact_L_s(L_s)

  # ---- Return ----
  create_named_list(Y, Z, time_obs, C, true_params)
}
