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

#' Construct factor-specific L2-orthonormal B-spline functions on [0,1]
#'
#' @param M Number of basis functions
#' @param n_g Dense grid size
#' @param degree Pool of B-spline degrees. Values are cycled across components.
#' @param block_id Positive identifier used to construct a distinct function block.
#' @return List containing the orthonormal functions, grid and construction metadata.
#' @keywords internal
construct_orthonormal_eigenfunctions <- function(M, n_g = 200,
                                                 degree = c(2L, 3L),
                                                 block_id = 1L) {
  if (length(M) != 1L || !is.finite(M) || M != as.integer(M) || M < 1L) {
    stop("M must be a positive integer.")
  }
  if (length(n_g) != 1L || !is.finite(n_g) || n_g != as.integer(n_g) ||
      n_g < max(50L, 10L * M)) {
    stop("n_g must be an integer of at least max(50, 10 * M).")
  }
  degree <- as.integer(degree)
  if (length(degree) < 1L || any(!is.finite(degree)) ||
      any(degree < 1L | degree > 5L)) {
    stop("degree must contain integers between 1 and 5.")
  }
  if (length(block_id) != 1L || !is.finite(block_id) ||
      block_id != as.integer(block_id) || block_id < 1L) {
    stop("block_id must be a positive integer.")
  }

  t_grid <- seq(0, 1, length.out = n_g)

  dt <- diff(t_grid)
  w_grid <- c(dt[1L] / 2,
              (dt[-1L] + dt[-length(dt)]) / 2,
              dt[length(dt)] / 2)

  # Jaoua-style truth: ordinary B-splines with equally spaced knots.  Each
  # block is generated independently and is never projected to the fitted
  # O'Sullivan basis.  Random coefficient mixtures avoid reusing the same
  # first M basis columns for every latent factor.
  for (attempt in seq_len(25L)) {
    raw_fun <- matrix(0, nrow = n_g, ncol = M)
    component_meta <- vector("list", M)
    for (m in seq_len(M)) {
      degree_m <- degree[((block_id + m - 2L) %% length(degree)) + 1L]
      n_internal <- max(3L, M + 2L + ((block_id + m) %% 3L))
      knots_m <- seq(0, 1, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
      bs_mat <- splines::bs(
        t_grid, knots = knots_m, degree = degree_m,
        Boundary.knots = c(0, 1), intercept = TRUE)
      coef_m <- stats::rnorm(ncol(bs_mat))
      raw_fun[, m] <- as.vector(bs_mat %*% coef_m)
      component_meta[[m]] <- list(
        degree = degree_m,
        knots = knots_m,
        coefficients = coef_m)
    }

    gram_raw <- crossprod(raw_fun, w_grid * raw_fun)
    chol_raw <- tryCatch(chol(gram_raw), error = function(e) NULL)
    if (!is.null(chol_raw) && rcond(gram_raw) > 1e-10) {
      phi <- raw_fun %*% backsolve(chol_raw, diag(M))
      for (m in seq_len(M)) {
        pivot <- which.max(abs(phi[, m]))
        if (phi[pivot, m] < 0) {
          phi[, m] <- -phi[, m]
          component_meta[[m]]$coefficients <-
            -component_meta[[m]]$coefficients
        }
      }
      gram_phi <- crossprod(phi, w_grid * phi)
      if (max(abs(gram_phi - diag(M))) < 1e-8) {
        return(list(
          phi = phi,
          t_grid = t_grid,
          w_grid = w_grid,
          gram = gram_phi,
          block_id = as.integer(block_id),
          components = component_meta))
      }
    }
  }
  stop("Unable to construct a numerically full-rank B-spline truth block.")
}

#' Build the scheme-1 FPCA scale for one latent factor
#' @keywords internal
make_scheme1_eigenvalues <- function(M, score_var_decay = TRUE) {
  raw_lambda <- if (score_var_decay) {
    1 / (seq_len(M)^2)
  } else {
    rep(1, M)
  }
  raw_lambda / sum(raw_lambda)
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
#' @noRd
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
#' @param bs_degree Degree pool for B-spline truth construction. A vector such
#'   as `c(2, 3)` alternates quadratic and cubic components.
#' @param identified_loadings If TRUE, construct loadings satisfying the
#'   population Gram conditions used in the thesis simulations.
#' @param n_obs_range Optional two-integer range for subject-specific numbers
#'   of observations. Requires `common_grid = FALSE`. When NULL, every subject
#'   has `n_obs` observations, preserving the legacy interface.
#' @param sparsity_mode Either `"fixed"` for the legacy fixed active proportion
#'   or `"factor_beta"` for factor-specific inclusion probabilities.
#' @param omega_beta Beta shape parameters used by `sparsity_mode =
#'   "factor_beta"`.
#'
#' @return A list with elements Y, Z, time_obs, C, true_params.
#'
#' @noRd
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
    bs_degree = c(2L, 3L),
    identified_loadings = FALSE,
    seed = NULL,
    n_obs_range = NULL,
    sparsity_mode = c("fixed", "factor_beta"),
    omega_beta = c(1, 10)
) {

  if (!is.null(seed)) set.seed(seed)
  sparsity_mode <- match.arg(sparsity_mode)

  # ---- Dimension validation ----
  stopifnot(length(n_s) == S)
  stopifnot(L_f >= 0)
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  if (d == 0) d <- 0
  d_use <- max(d, 0)

  if (is.null(M_f)) {
    M_f <- if (L_f > 0) rep(3, L_f) else integer(0)
  }
  stopifnot(length(M_f) == max(0, L_f))

  if (is.null(M_s)) {
    if (has_specific) {
      M_s <- lapply(seq_len(S), function(s) rep(3L, L_s[[s]]))
    } else {
      M_s <- replicate(S, integer(0), simplify = FALSE)
    }
  }
  stopifnot(is.list(M_s), length(M_s) == S)
  for (s in 1:S) {
    stopifnot(length(M_s[[s]]) == L_s[[s]])
  }

  if (length(n_obs) != 1L || !is.finite(n_obs) ||
      n_obs != as.integer(n_obs) || n_obs < 2L) {
    stop("n_obs must be an integer of at least 2.")
  }
  if (!is.null(n_obs_range)) {
    if (isTRUE(common_grid)) {
      stop("n_obs_range requires common_grid = FALSE.")
    }
    if (length(n_obs_range) != 2L || any(!is.finite(n_obs_range)) ||
        any(n_obs_range != as.integer(n_obs_range)) ||
        n_obs_range[1L] < 2L || n_obs_range[2L] < n_obs_range[1L]) {
      stop("n_obs_range must contain two ordered integers, both at least 2.")
    }
    n_obs_range <- as.integer(n_obs_range)
  }
  if (length(n_dense) != 1L || !is.finite(n_dense) ||
      n_dense != as.integer(n_dense) || n_dense < 50L) {
    stop("n_dense must be an integer of at least 50.")
  }
  bs_degree <- as.integer(bs_degree)
  if (length(bs_degree) < 1L || any(!is.finite(bs_degree)) ||
      any(bs_degree < 1L | bs_degree > 5L)) {
    stop("bs_degree must contain integers between 1 and 5.")
  }
  if (length(omega_beta) != 2L || any(!is.finite(omega_beta)) ||
      any(omega_beta <= 0)) {
    stop("omega_beta must contain two positive finite shape parameters.")
  }

  # ---- Step 1: Construct observation time grid ----
  time_obs <- vector("list", S)
  n_obs_by_subject <- vector("list", S)
  if (common_grid) {
    t_grid_obs <- seq(0, 1, length.out = n_obs)
    for (s in 1:S) {
      time_obs[[s]] <- lapply(1:n_s[s], function(i) t_grid_obs)
      n_obs_by_subject[[s]] <- rep(as.integer(n_obs), n_s[s])
    }
  } else {
    for (s in 1:S) {
      n_obs_by_subject[[s]] <- if (is.null(n_obs_range)) {
        rep(as.integer(n_obs), n_s[s])
      } else {
        sample(seq.int(n_obs_range[1L], n_obs_range[2L]),
               size = n_s[s], replace = TRUE)
      }
      time_obs[[s]] <- lapply(n_obs_by_subject[[s]], function(n_i) {
        sort(runif(n_i, 0, 1))
      })
    }
  }

  if (is.null(K)) {
    typical_n <- stats::median(unlist(n_obs_by_subject, use.names = FALSE))
    K <- max(7L, min(floor(typical_n / 4), 40L))
  }
  if (length(K) != 1L || !is.finite(K) || K != as.integer(K) || K < 2L) {
    stop("K must be an integer of at least 2.")
  }
  K <- as.integer(K)
  K_total <- K + 2L

  # ---- Step 2: Construct spline design matrix C (O'Sullivan) ----
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

  # ---- Step 3: Generate true structured functions ----

  # 3a: Construct a distinct orthonormal truth block for every shared factor.
  # Theoretical functions have unit norm.  Scheme-1 computational functions
  # absorb sqrt(lambda), while the legacy zeta interface stores unit-variance
  # computational scores eta.
  t_grid_dense <- seq(0, 1, length.out = n_dense)
  phi_dense_list <- vector("list", L_f)
  theta_dense_list <- vector("list", L_f)
  lambda_phi_true <- vector("list", L_f)
  phi_basis_meta <- vector("list", L_f)
  block_id <- 0L
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    block_id <- block_id + 1L
    ef_obj <- construct_orthonormal_eigenfunctions(
      M_l, n_g = n_dense, degree = bs_degree, block_id = block_id)
    lambda_l <- make_scheme1_eigenvalues(M_l, score_var_decay)
    phi_dense_list[[l]] <- ef_obj$phi
    theta_dense_list[[l]] <- sweep(
      ef_obj$phi, 2L, sqrt(lambda_l), FUN = "*")
    lambda_phi_true[[l]] <- lambda_l
    phi_basis_meta[[l]] <- ef_obj[c("block_id", "components", "gram")]
  }

  # 3b: Construct separate study-specific truth blocks.  They are not reused
  # across studies or factors, preventing artificial shared/specific aliases.
  psi_dense_list <- NULL
  kappa_dense_list <- NULL
  lambda_psi_true <- NULL
  psi_basis_meta <- NULL
  if (has_specific) {
    psi_dense_list <- vector("list", S)
    kappa_dense_list <- vector("list", S)
    lambda_psi_true <- vector("list", S)
    psi_basis_meta <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      psi_dense_list[[s]] <- vector("list", L_ss)
      kappa_dense_list[[s]] <- vector("list", L_ss)
      lambda_psi_true[[s]] <- vector("list", L_ss)
      psi_basis_meta[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        block_id <- block_id + 1L
        ef_obj <- construct_orthonormal_eigenfunctions(
          M_sl, n_g = n_dense, degree = bs_degree, block_id = block_id)
        lambda_sl <- make_scheme1_eigenvalues(M_sl, score_var_decay)
        psi_dense_list[[s]][[l]] <- ef_obj$phi
        kappa_dense_list[[s]][[l]] <- sweep(
          ef_obj$phi, 2L, sqrt(lambda_sl), FUN = "*")
        lambda_psi_true[[s]][[l]] <- lambda_sl
        psi_basis_meta[[s]][[l]] <-
          ef_obj[c("block_id", "components", "gram")]
      }
    }
  }

  # 3c: Generate phase shifts for mean functions and beta functions
  phase_mu <- runif(p * S, 0, 2 * pi)
  phase_mu_mat <- matrix(phase_mu, nrow = S, ncol = p)

  phase_beta <- NULL
  if (d_use > 0) {
    n_phase_beta <- p * d_use
    phase_beta <- runif(n_phase_beta, 0, 2 * pi)
  }

  # 3d: Generate scheme-1 scores.  The legacy fields zeta_true and xi_true
  # are computational eta/chi scores with unit variance.  The theoretical
  # FPCA scores are stored separately and equal sqrt(lambda) times eta/chi.
  zeta_true <- vector("list", S)
  zeta_fpca_true <- vector("list", S)
  for (s in 1:S) {
    zeta_true[[s]] <- vector("list", L_f)
    zeta_fpca_true[[s]] <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      zeta_true[[s]][[l]] <- matrix(
        rnorm(n_s[s] * M_l, mean = 0, sd = 1),
        nrow = n_s[s], ncol = M_l)
      zeta_fpca_true[[s]][[l]] <- sweep(
        zeta_true[[s]][[l]], 2L, sqrt(lambda_phi_true[[l]]), FUN = "*")
    }
  }
  eta_true <- zeta_true

  xi_true <- NULL
  xi_fpca_true <- NULL
  chi_true <- NULL
  if (has_specific) {
    xi_true <- vector("list", S)
    xi_fpca_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      xi_true[[s]] <- vector("list", L_ss)
      xi_fpca_true[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        xi_true[[s]][[l]] <- matrix(
          rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
          nrow = n_s[s], ncol = M_sl)
        xi_fpca_true[[s]][[l]] <- sweep(
          xi_true[[s]][[l]], 2L,
          sqrt(lambda_psi_true[[s]][[l]]), FUN = "*")
      }
    }
    chi_true <- xi_true
  }

  # 3e: Generate loadings with spike-and-slab
  loading_draw <- if (identified_loadings) {
    generate_identified_loadings(p, L_f, L_s, S,
      sparse = bool_sparse_loadings, prop_sparse = prop_sparse)
  } else NULL
  a_true <- if (identified_loadings) loading_draw$a else matrix(0, p, L_f)
  gamma_a_true <- if (identified_loadings) loading_draw$gamma_a else matrix(0, p, L_f)
  omega_a_true <- if (L_f > 0L) rep(NA_real_, L_f) else numeric(0)

  draw_active_indices <- function(omega) {
    if (!bool_sparse_loadings) return(seq_len(p))
    if (sparsity_mode == "fixed") {
      return(sample(seq_len(p), max(1L, round(p * omega))))
    }
    active <- which(stats::runif(p) < omega)
    if (length(active) == 0L) active <- sample(seq_len(p), 1L)
    active
  }

  for (l in if (identified_loadings) integer(0) else seq_len(L_f)) {
    omega_a_true[l] <- if (!bool_sparse_loadings) {
      1
    } else if (sparsity_mode == "fixed") {
      1 - prop_sparse
    } else {
      stats::rbeta(1L, omega_beta[1L], omega_beta[2L])
    }
    active_idx <- draw_active_indices(omega_a_true[l])
    if (bool_sparse_loadings) {
      val_mean <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      a_true[active_idx, l] <- rnorm(length(active_idx), mean = val_mean, sd = 0.2)
      gamma_a_true[active_idx, l] <- 1
    } else {
      a_true[, l] <- rnorm(p, mean = 0, sd = 0.5)
      gamma_a_true[, l] <- 1
    }
  }
  if (identified_loadings && L_f > 0L) {
    omega_a_true <- colMeans(gamma_a_true)
  }

  b_true <- if (identified_loadings) loading_draw$b else NULL
  gamma_b_true <- if (identified_loadings) loading_draw$gamma_b else NULL
  omega_b_true <- if (has_specific) vector("list", S) else NULL
  if (has_specific && !identified_loadings) {
    b_true <- vector("list", S)
    gamma_b_true <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
      gamma_b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
      omega_b_true[[s]] <- rep(NA_real_, L_ss)
      for (l in seq_len(L_ss)) {
        omega_b_true[[s]][l] <- if (!bool_sparse_loadings) {
          1
        } else if (sparsity_mode == "fixed") {
          1 - prop_sparse
        } else {
          stats::rbeta(1L, omega_beta[1L], omega_beta[2L])
        }
        active_idx <- draw_active_indices(omega_b_true[[s]][l])
        if (bool_sparse_loadings) {
          val_mean <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
          b_true[[s]][active_idx, l] <-
            rnorm(length(active_idx), mean = val_mean, sd = 0.2)
          gamma_b_true[[s]][active_idx, l] <- 1
        } else {
          b_true[[s]][, l] <- rnorm(p, mean = 0, sd = 0.5)
          gamma_b_true[[s]][, l] <- 1
        }
      }
    }
  }
  if (has_specific && identified_loadings) {
    for (s in seq_len(S)) {
      omega_b_true[[s]] <- colMeans(gamma_b_true[[s]])
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
  signal_true_values <- vector("list", S)
  noise_true_values <- vector("list", S)

  for (s in 1:S) {
    L_ss <- L_s[[s]]
    Y[[s]] <- vector("list", n_s[s])
    mu_true_values[[s]] <- vector("list", n_s[s])
    f_true_values[[s]] <- vector("list", n_s[s])
    g_true_values[[s]] <- vector("list", n_s[s])
    beta_true_values[[s]] <- vector("list", n_s[s])
    signal_true_values[[s]] <- vector("list", n_s[s])
    noise_true_values[[s]] <- vector("list", n_s[s])

    for (i in 1:n_s[s]) {
      Y[[s]][[i]] <- vector("list", p)
      t_i <- time_obs[[s]][[i]]
      n_i <- length(t_i)

      mu_si <- matrix(0, n_i, p)
      f_si <- matrix(0, n_i, L_f)
      g_si <- if (L_ss > 0L) matrix(0, n_i, L_ss) else NULL
      beta_si_contrib <- matrix(0, nrow = n_i, ncol = p)
      signal_si <- matrix(0, nrow = n_i, ncol = p)
      noise_si <- matrix(0, nrow = n_i, ncol = p)

      # Mean functions
      for (j in 1:p) {
        mu_si[, j] <- mean_function_periodic(t_i, phase_mu_mat[s, j], mean_amp)
      }

      # Shared factor processes
      for (l in seq_len(L_f)) {
        M_l <- M_f[l]
        eta_i <- zeta_true[[s]][[l]][i, ]
        theta_val <- matrix(0, n_i, M_l)
        for (m in 1:M_l) {
          theta_val[, m] <- eval_eigenfun(
            t_i, theta_dense_list[[l]][, m], t_grid_dense)
        }
        f_si[, l] <- as.vector(theta_val %*% eta_i)
      }

      # Specific factor processes
      if (L_ss > 0L) {
        for (l in seq_len(L_ss)) {
          M_sl <- M_s[[s]][l]
          chi_i <- xi_true[[s]][[l]][i, ]
          kappa_val <- matrix(0, n_i, M_sl)
          for (m in 1:M_sl) {
            kappa_val[, m] <- eval_eigenfun(
              t_i, kappa_dense_list[[s]][[l]][, m], t_grid_dense)
          }
          g_si[, l] <- as.vector(kappa_val %*% chi_i)
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
            beta_si_contrib[, j] <- beta_si_contrib[, j] +
              Z[[s]][i, r] * beta_val
          }
        }

        # Shared factor contributions
        for (l in seq_len(L_f)) {
          y_sij <- y_sij + a_true[j, l] * f_si[, l]
        }

        # Specific factor contributions
        if (L_ss > 0L) {
          for (l in seq_len(L_ss)) {
            y_sij <- y_sij + b_true[[s]][j, l] * g_si[, l]
          }
        }

        signal_si[, j] <- y_sij
        noise_si[, j] <- rnorm(
          n_i, mean = 0, sd = sqrt(sigma2_eps_true[s, j]))
        Y[[s]][[i]][[j]] <- signal_si[, j] + noise_si[, j]
      }

      mu_true_values[[s]][[i]] <- mu_si
      f_true_values[[s]][[i]] <- f_si
      if (L_ss > 0L) g_true_values[[s]][[i]] <- g_si
      beta_true_values[[s]][[i]] <- beta_si_contrib
      signal_true_values[[s]][[i]] <- signal_si
      noise_true_values[[s]][[i]] <- noise_si
    }
  }

  # ---- Step 5: Assemble true parameter list ----
  truth_basis_family <- "ordinary_B_spline_external"
  truth_projected_to_fitted_basis <- FALSE
  true_params <- create_named_list(
    a_true,
    gamma_a_true,
    omega_a_true,
    b_true,
    gamma_b_true,
    omega_b_true,
    zeta_true,
    xi_true,
    eta_true,
    chi_true,
    zeta_fpca_true,
    xi_fpca_true,
    phi_dense_list,
    psi_dense_list,
    theta_dense_list,
    kappa_dense_list,
    lambda_phi_true,
    lambda_psi_true,
    phi_basis_meta,
    psi_basis_meta,
    t_grid_dense,
    phase_mu_mat,
    phase_beta,
    sigma2_eps_true,
    mu_true_values,
    f_true_values,
    g_true_values,
    beta_true_values,
    signal_true_values,
    noise_true_values,
    K, L_f, L_s, M_f, M_s,
    S, n_s, p, d_use, n_obs,
    n_obs_by_subject,
    n_obs_range,
    score_var_decay,
    bs_degree,
    sparsity_mode,
    omega_beta,
    int_knots,
    identified_loadings,
    truth_basis_family,
    truth_projected_to_fitted_basis
  )
  true_params$L_s_by_study <- L_s
  true_params$L_s <- .compact_L_s(L_s)

  create_named_list(Y, Z, time_obs, C, true_params)
}
