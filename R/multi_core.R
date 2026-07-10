# =============================================================================
# multi_core.R — Multi-study functional factor model: main functions
#
# Contains:
#   bayesSYNC_multi()       — User-facing interface
#   bayesSYNC_multi_core()  — Initialization, precomputation, main CAVI loop
#   compute_elbo_multi()    — ELBO computation
#
# Symbol convention (aligned with 创新点.pdf):
#   L_f = number of shared factors (scalar)
#   L_s = number of study-specific factors (scalar)
#   M_f = vector length L_f: M_f[l] = FPCA truncation for shared factor l
#   M_s = list of S vectors: M_s[[s]][l] = FPCA truncation per study-factor
#
# Indexing convention: study [[s]] is ALWAYS the outermost dimension.
# [[s]] → study, [[i]] → individual, [[j]] → variable, [[l]] → factor,
# [[m]] → FPCA component, [[r]] → covariate
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Multi-Study Bayesian Functional Factor Model
#'
#' @param Y List of length S, where Y[[s]] is a list of n_s[s] individuals,
#'   each being a list of p variable vectors (length n_i observations).
#' @param Z Optional list of length S, Z[[s]] is an n_s[s] x d covariate matrix.
#' @param time_obs List of length S, time_obs[[s]][[i]] is the observation
#'   time vector for individual i in study s.
#' @param L_f Number of shared latent factors (scalar).
#' @param L_s Number of study-specific latent factors (scalar, 0 for none).
#' @param M_f Vector of length L_f, FPCA truncation per shared factor.
#' @param M_s List of S vectors, FPCA truncation per study-specific factor.
#' @param K Number of O'Sullivan spline functions (default NULL: auto).
#' @param anneal Annealing parameters: c(type, T_max, grid_size) or NULL.
#' @param list_hyper Hyperparameter list from set_hyper() or NULL for defaults.
#' @param n_g Dense grid size.
#' @param time_g Dense grid vector (overrides n_g).
#' @param tol_abs Absolute convergence tolerance on ELBO.
#' @param tol_rel Relative convergence tolerance on ELBO.
#' @param maxit Maximum number of iterations.
#' @param n_cpus Number of CPU cores for parallel execution.
#' @param verbose Print progress messages.
#' @param seed Random seed for reproducibility.
#' @param bool_scale Whether to standardise variables.
#' @param bool_var_spec_prob Use variable-specific spike-and-slab probabilities.
#'
#' @return A list with fitted model results.
#'
#' @export
bayesSYNC_multi <- function(
    Y, Z = NULL, time_obs = NULL,
    L_f, L_s, M_f, M_s,
    K = NULL,
    anneal = c(1, 1.9, 100),
    list_hyper = NULL,
    n_g = 1000, time_g = NULL,
    tol_abs = 1e-3, tol_rel = 1e-5,
    maxit = 1000, n_cpus = 1,
    verbose = TRUE, seed = NULL,
    bool_scale = TRUE, bool_var_spec_prob = FALSE
) {

  if (!is.null(seed)) {
    cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

  check_annealing(anneal, verbose)

  # ---- Input validation ----
  if (!is.list(Y)) stop("Y must be a list of studies.")
  S <- length(Y)
  stopifnot(S >= 1)

  n_s <- sapply(Y, length)

  p <- length(Y[[1]][[1]])
  for (s in 1:S) {
    for (i in 1:n_s[s]) {
      stopifnot(length(Y[[s]][[i]]) == p)
    }
  }

  check_structure(L_f, "vector", "numeric", 1)
  stopifnot(L_f >= 0)
  check_structure(L_s, "vector", "numeric", 1)
  stopifnot(L_s >= 0)

  if (L_f > 0) {
    check_structure(M_f, "vector", "numeric")
    stopifnot(all(M_f >= 0))
  }
  stopifnot(length(M_f) == max(0, L_f))

  stopifnot(is.list(M_s), length(M_s) == S)
  if (L_s > 0) {
    for (s in 1:S) {
      check_structure(M_s[[s]], "vector", "numeric")
      stopifnot(length(M_s[[s]]) == L_s)
      stopifnot(all(M_s[[s]] >= 0))
    }
  }

  # ---- Covariates ----
  if (!is.null(Z)) {
    d <- ncol(Z[[1]])
    for (s in 1:S) {
      stopifnot(nrow(Z[[s]]) == n_s[s])
      stopifnot(ncol(Z[[s]]) == d)
    }
  } else {
    d <- 0
  }

  # ---- Hyperparameters ----
  if (is.null(list_hyper)) {
    list_hyper <- set_hyper(d_0 = p)
  }

  # ---- Flat observation grid ----
  time_obs_flat <- unlist(lapply(1:S, function(s) time_obs[[s]]), recursive = FALSE)

  if (is.null(K)) {
    n_times <- sapply(time_obs_flat, length)
    K <- max(7, min(floor(min(n_times) / 4), 40))
    if (verbose) cat(paste0("K set to ", K, "\n"))
  } else {
    check_natural(K)
  }

  grid_obj <- get_grid_objects(time_obs_flat, K, n_g = n_g, time_g = time_g,
                                format_univ = TRUE)
  C_flat <- grid_obj$C
  n_g <- grid_obj$n_g
  time_g <- grid_obj$time_g
  C_g <- grid_obj$C_g

  # Restructure C by study
  C <- vector("list", S)
  idx_start <- 1
  for (s in 1:S) {
    C[[s]] <- C_flat[idx_start:(idx_start + n_s[s] - 1)]
    idx_start <- idx_start + n_s[s]
  }

  # ---- Scaling ----
  if (bool_scale) {
    mean_per_subject <- lapply(1:p, function(j) {
      unlist(lapply(1:S, function(s) {
        sapply(1:n_s[s], function(i) mean(Y[[s]][[i]][[j]]))
      }))
    })
    mean_mean_across_subjects <- sapply(mean_per_subject, mean)
    sd_mean_across_subjects <- sapply(mean_per_subject, function(x) {
      sdx <- sd(x); ifelse(sdx < 1e-10, 1, sdx)
    })

    for (s in 1:S) {
      for (i in 1:n_s[s]) {
        Y[[s]][[i]] <- lapply(1:p, function(j) {
          (Y[[s]][[i]][[j]] - mean_mean_across_subjects[j]) /
            sd_mean_across_subjects[j]
        })
      }
    }
  } else {
    mean_mean_across_subjects <- sd_mean_across_subjects <- NULL
  }

  # ---- Variable names ----
  var_names <- names(Y[[1]][[1]])
  if (is.null(var_names)) {
    var_names <- paste0("variable_", 1:p)
    for (s in 1:S) {
      for (i in 1:n_s[s]) {
        names(Y[[s]][[i]]) <- var_names
      }
    }
  }

  # ---- Call core ----
  res <- bayesSYNC_multi_core(
    S = S, n_s = n_s, p = p, d = d,
    L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s,
    K = K, C = C, C_g = C_g,
    Y = Y, Z = Z,
    mean_mean_across_subjects = mean_mean_across_subjects,
    sd_mean_across_subjects = sd_mean_across_subjects,
    anneal = anneal, list_hyper = list_hyper,
    time_obs = time_obs, bool_var_spec_prob = bool_var_spec_prob,
    n_g = n_g, time_g = time_g,
    tol_abs = tol_abs, tol_rel = tol_rel, maxit = maxit,
    n_cpus = n_cpus, verbose = verbose
  )

  res
}


#' Core multi-study CAVI algorithm with annealing
#'
#' @keywords internal
bayesSYNC_multi_core <- function(
    S, n_s, p, d,
    L_f, L_s, M_f, M_s,
    K, C, C_g,
    Y, Z,
    mean_mean_across_subjects,
    sd_mean_across_subjects,
    anneal, list_hyper,
    time_obs, bool_var_spec_prob,
    n_g, time_g,
    tol_abs, tol_rel, maxit,
    n_cpus, verbose
) {

  eps <- .Machine$double.eps^0.5
  eps_elbo <- 1e-11

  # ---- Annealing setup ----
  if (is.null(anneal)) {
    annealing <- FALSE
    c_val <- 1
    i_iter_init <- 1
  } else {
    annealing <- TRUE
    ladder <- get_annealing_ladder_(anneal, verbose)
    c_val <- ladder[1]
    i_iter_init <- anneal[3]
  }

  K_total <- K + 2

  # ---- Hyperparameters ----
  sigma_zeta <- list_hyper$sigma_zeta
  mu_beta <- list_hyper$mu_beta
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A
  c_0 <- list_hyper$c_0
  d_0 <- list_hyper$d_0

  inv_Sigma_beta <- solve(Sigma_beta)

  # ---- Precomputation ----
  list_cp_C <- vector("list", S)
  list_cp_C_Y <- vector("list", S)
  list_cp_Y <- vector("list", S)
  sum_list_cp_C <- vector("list", S)

  for (s in 1:S) {
    list_cp_C[[s]] <- parallel::mclapply(1:n_s[s], function(i) {
      crossprod(C[[s]][[i]])
    }, mc.cores = n_cpus)

    list_cp_C_Y[[s]] <- parallel::mclapply(1:n_s[s], function(i) {
      sapply(1:p, function(j) crossprod(C[[s]][[i]], Y[[s]][[i]][[j]]))
    }, mc.cores = n_cpus)

    list_cp_Y[[s]] <- simplify2array(parallel::mclapply(1:p, function(j) {
      sapply(1:n_s[s], function(i) crossprod(Y[[s]][[i]][[j]]))
    }, mc.cores = n_cpus))

    sum_list_cp_C[[s]] <- Reduce("+", list_cp_C[[s]])
  }

  # ---- Total observation count per study ----
  sum_obs_s <- sapply(1:S, function(s) {
    sum(sapply(time_obs[[s]], length))
  })

  # ===================================================================
  #  INITIALIZATION OF VARIATIONAL PARAMETERS
  # ===================================================================

  # --- 1. Mean function spline coefficients: nu_{mu,sj} ---
  mu_q_recip_sigsq_mu <- matrix(1, nrow = S, ncol = p)
  mu_q_recip_a_mu <- matrix(1, nrow = S, ncol = p)

  mu_q_nu_mu <- vector("list", S)
  Sigma_q_nu_mu <- vector("list", S)
  inv_Sigma_q_nu_mu <- vector("list", S)

  for (s in 1:S) {
    mu_q_nu_mu[[s]] <- lapply(1:p, function(j) rnorm(K_total, mean = 0, sd = 1))
    Sigma_q_nu_mu[[s]] <- lapply(1:p, function(j) diag(K_total))
    inv_Sigma_q_nu_mu[[s]] <- lapply(1:p, function(j) diag(K_total))
  }

  # --- 2. Beta regression coefficients: nu_{beta,jr} ---
  mu_q_nu_beta <- NULL
  Sigma_q_nu_beta <- NULL
  inv_Sigma_q_nu_beta <- NULL
  mu_q_recip_sigsq_beta <- NULL
  mu_q_recip_a_beta <- NULL

  if (d > 0) {
    mu_q_recip_sigsq_beta <- matrix(1, nrow = p, ncol = d)
    mu_q_recip_a_beta <- matrix(1, nrow = p, ncol = d)

    mu_q_nu_beta <- lapply(1:p, function(j) {
      lapply(1:d, function(r) rnorm(K_total, mean = 0, sd = 1))
    })
    inv_Sigma_q_nu_beta <- lapply(1:p, function(j) {
      lapply(1:d, function(r) {
        blkdiag(inv_Sigma_beta, mu_q_recip_sigsq_beta[j, r] * diag(K))
      })
    })
    Sigma_q_nu_beta <- lapply(1:p, function(j) {
      lapply(1:d, function(r) diag(K_total))
    })
  }

  # --- 3. Shared eigenfunction coefficients: nu_{phi,ml} ---
  # Storage is IRREGULAR: mu_q_nu_phi[[l]] is K_total x M_f[l]
  mu_q_recip_sigsq_phi <- if (L_f > 0) lapply(1:L_f, function(l) rep(1, M_f[l])) else list()
  mu_q_recip_a_phi <- if (L_f > 0) lapply(1:L_f, function(l) rep(1, M_f[l])) else list()

  mu_q_nu_phi <- vector("list", L_f)
  inv_Sigma_q_nu_phi <- vector("list", L_f)
  Sigma_q_nu_phi <- vector("list", L_f)

  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    mu_q_nu_phi[[l]] <- matrix(rnorm(K_total * M_l, mean = 0, sd = 1),
                                nrow = K_total, ncol = M_l)
    inv_Sigma_q_nu_phi[[l]] <- lapply(1:M_l, function(m) {
      blkdiag(inv_Sigma_beta, mu_q_recip_sigsq_phi[[l]][m] * diag(K))
    })
    Sigma_q_nu_phi[[l]] <- lapply(1:M_l, function(m) {
      diag(K_total)
    })
  }

  # --- 4. Study-specific eigenfunction coefficients: nu_{psi,slm} ---
  mu_q_nu_psi <- NULL
  inv_Sigma_q_nu_psi <- NULL
  Sigma_q_nu_psi <- NULL
  mu_q_recip_sigsq_psi <- NULL
  mu_q_recip_a_psi <- NULL

  if (L_s > 0) {
    mu_q_recip_sigsq_psi <- lapply(1:S, function(s) {
      lapply(1:L_s, function(l) rep(1, M_s[[s]][l]))
    })
    mu_q_recip_a_psi <- lapply(1:S, function(s) {
      lapply(1:L_s, function(l) rep(1, M_s[[s]][l]))
    })

    mu_q_nu_psi <- vector("list", S)
    inv_Sigma_q_nu_psi <- vector("list", S)
    Sigma_q_nu_psi <- vector("list", S)

    for (s in 1:S) {
      mu_q_nu_psi[[s]] <- vector("list", L_s)
      inv_Sigma_q_nu_psi[[s]] <- vector("list", L_s)
      Sigma_q_nu_psi[[s]] <- vector("list", L_s)

      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        mu_q_nu_psi[[s]][[l]] <- matrix(rnorm(K_total * M_sl, mean = 0, sd = 1),
                                         nrow = K_total, ncol = M_sl)
        inv_Sigma_q_nu_psi[[s]][[l]] <- lapply(1:M_sl, function(m) {
          blkdiag(inv_Sigma_beta,
                  mu_q_recip_sigsq_psi[[s]][[l]][m] * diag(K))
        })
        Sigma_q_nu_psi[[s]][[l]] <- lapply(1:M_sl, function(m) {
          diag(K_total)
        })
      }
    }
  }

  # --- 5. Shared factor scores: zeta^{(l)}_{si} ---
  mu_q_zeta <- vector("list", S)
  Sigma_q_zeta <- vector("list", S)

  for (s in 1:S) {
    mu_q_zeta[[s]] <- vector("list", L_f)
    Sigma_q_zeta[[s]] <- vector("list", L_f)
    for (l in seq_len(L_f)) {
      M_l <- M_f[l]
      mu_q_zeta[[s]][[l]] <- matrix(rnorm(n_s[s] * M_l, mean = 0, sd = 1),
                                     nrow = n_s[s], ncol = M_l)
      Sigma_q_zeta[[s]][[l]] <- lapply(1:n_s[s], function(i) diag(M_l))
    }
  }

  # --- 6. Study-specific factor scores: xi^{(sl)}_{si} ---
  mu_q_xi <- NULL
  Sigma_q_xi <- NULL

  if (L_s > 0) {
    mu_q_xi <- vector("list", S)
    Sigma_q_xi <- vector("list", S)
    for (s in 1:S) {
      mu_q_xi[[s]] <- vector("list", L_s)
      Sigma_q_xi[[s]] <- vector("list", L_s)
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        mu_q_xi[[s]][[l]] <- matrix(rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
                                     nrow = n_s[s], ncol = M_sl)
        Sigma_q_xi[[s]][[l]] <- lapply(1:n_s[s], function(i) diag(M_sl))
      }
    }
  }

  # --- 7. Shared loadings a_{jl} and PPI gamma_{jl} (Spike-and-Slab) ---
  mu_q_normal_a <- matrix(rnorm(p * max(1, L_f)), nrow = p, ncol = max(1, L_f))
  Sigma_q_normal_a <- matrix(1, nrow = p, ncol = max(1, L_f))
  mu_q_gamma_a <- matrix(0.5, nrow = p, ncol = max(1, L_f))
  if (L_f == 0) {
    mu_q_normal_a <- Sigma_q_normal_a <- mu_q_gamma_a <- matrix(NA, nrow = p, ncol = 0)
  }
  mu_q_a <- mu_q_gamma_a * mu_q_normal_a
  term_a <- (Sigma_q_normal_a + mu_q_normal_a^2) * mu_q_gamma_a

  # --- 8. Study-specific loadings b_{sjl} and PPI ---
  mu_q_normal_b <- NULL
  Sigma_q_normal_b <- NULL
  mu_q_gamma_b <- NULL
  mu_q_b_specific <- NULL
  term_b_specific <- NULL

  if (L_s > 0) {
    mu_q_normal_b <- vector("list", S)
    Sigma_q_normal_b <- vector("list", S)
    mu_q_gamma_b <- vector("list", S)
    mu_q_b_specific <- vector("list", S)
    term_b_specific <- vector("list", S)

    for (s in 1:S) {
      mu_q_normal_b[[s]] <- matrix(rnorm(p * L_s), nrow = p, ncol = L_s)
      Sigma_q_normal_b[[s]] <- matrix(1, nrow = p, ncol = L_s)
      mu_q_gamma_b[[s]] <- matrix(0.5, nrow = p, ncol = L_s)
      mu_q_b_specific[[s]] <- mu_q_gamma_b[[s]] * mu_q_normal_b[[s]]
      term_b_specific[[s]] <- (Sigma_q_normal_b[[s]] +
                               mu_q_normal_b[[s]]^2) * mu_q_gamma_b[[s]]
    }
  }

  # --- 9. Variance parameters ---
  mu_q_recip_sigsq_eps <- matrix(1, nrow = S, ncol = p)
  mu_q_recip_a_eps <- matrix(1, nrow = S, ncol = p)

  # --- 10. Omega parameters ---
  c_1_omega_a <- rep(NA, L_f)
  d_1_omega_a <- rep(NA, L_f)
  mu_q_log_omega_a <- rep(NA, L_f)
  mu_q_log_1_omega_a <- rep(NA, L_f)

  c_1_omega_b <- NULL
  d_1_omega_b <- NULL
  mu_q_log_omega_b <- NULL
  mu_q_log_1_omega_b <- NULL

  if (L_s > 0) {
    c_1_omega_b <- vector("list", S)
    d_1_omega_b <- vector("list", S)
    mu_q_log_omega_b <- vector("list", S)
    mu_q_log_1_omega_b <- vector("list", S)
    for (s in 1:S) {
      c_1_omega_b[[s]] <- rep(NA, L_s)
      d_1_omega_b[[s]] <- rep(NA, L_s)
      mu_q_log_omega_b[[s]] <- rep(NA, L_s)
      mu_q_log_1_omega_b[[s]] <- rep(NA, L_s)
    }
  }

  # ---- Initial omega update ----
  if (L_f > 0) {
    cs_mu_q_gamma_a <- colSums(mu_q_gamma_a)
    c_1_omega_a <- c_val * (c_0 + cs_mu_q_gamma_a) - c_val + 1
    d_1_omega_a <- c_val * (d_0 + p - cs_mu_q_gamma_a) - c_val + 1
    dig_a <- digamma(c_val * (c_0 + d_0 + p) - 2 * c_val + 2)
    mu_q_log_omega_a <- digamma(c_1_omega_a) - dig_a
    mu_q_log_1_omega_a <- digamma(d_1_omega_a) - dig_a
  }

  if (L_s > 0) {
    for (s in 1:S) {
      cs_mu_q_gamma_b <- colSums(mu_q_gamma_b[[s]])
      c_1_omega_b[[s]] <- c_val * (c_0 + cs_mu_q_gamma_b) - c_val + 1
      d_1_omega_b[[s]] <- c_val * (d_0 + p - cs_mu_q_gamma_b) - c_val + 1
      dig_b <- digamma(c_val * (c_0 + d_0 + p) - 2 * c_val + 2)
      mu_q_log_omega_b[[s]] <- digamma(c_1_omega_b[[s]]) - dig_b
      mu_q_log_1_omega_b[[s]] <- digamma(d_1_omega_b[[s]]) - dig_b
    }
  }

  # ===================================================================
  #  MAIN CAVI LOOP
  # ===================================================================

  ELBO <- NULL

  for (i_iter in 1:maxit) {

    # ------ Block 1: Update q(nu_mu) ------
    res_mu <- update_nu_mu(
      Y = Y, C = C,
      list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      sum_list_cp_C = sum_list_cp_C,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_sigsq_mu = mu_q_recip_sigsq_mu,
      inv_Sigma_beta = inv_Sigma_beta,
      mu_q_nu_beta = mu_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
      mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
      mu_q_a = mu_q_a, mu_q_b_specific = mu_q_b_specific,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, K = K, K_total = K_total,
      c_val = c_val, n_cpus = n_cpus)
    mu_q_nu_mu <- res_mu$mu_q_nu_mu
    Sigma_q_nu_mu <- res_mu$Sigma_q_nu_mu

    # ------ Block 2: Update q(nu_beta) ------
    if (d > 0 && !is.null(mu_q_nu_beta)) {
      res_beta <- update_nu_beta(
        Y = Y, C = C,
        list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
        mu_q_nu_mu = mu_q_nu_mu,
        mu_q_nu_beta = mu_q_nu_beta, Sigma_q_nu_beta = Sigma_q_nu_beta, Z = Z,
        mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
        mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
        mu_q_a = mu_q_a, mu_q_b_specific = mu_q_b_specific,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        mu_q_recip_sigsq_beta = mu_q_recip_sigsq_beta,
        inv_Sigma_beta = inv_Sigma_beta,
        S = S, n_s = n_s, p = p, d = d,
        L_f = L_f, L_s = L_s, K = K, K_total = K_total,
        c_val = c_val, n_cpus = n_cpus)
      mu_q_nu_beta <- res_beta$mu_q_nu_beta
      Sigma_q_nu_beta <- res_beta$Sigma_q_nu_beta
    }

    # ------ Block 3: Update q(nu_phi) ------
    res_phi <- update_nu_phi(
      Y = Y, C = C,
      list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi,
      mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
      mu_q_a = mu_q_a, term_a = term_a,
      mu_q_b_specific = mu_q_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      inv_Sigma_beta = inv_Sigma_beta,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, K = K, K_total = K_total,
      c_val = c_val, n_cpus = n_cpus)
    mu_q_nu_phi <- res_phi$mu_q_nu_phi
    Sigma_q_nu_phi <- res_phi$Sigma_q_nu_phi
    inv_Sigma_q_nu_phi <- res_phi$inv_Sigma_q_nu_phi

    # ------ Block 4: Update q(nu_psi) ------
    if (L_s > 0) {
      res_psi <- update_nu_psi(
        Y = Y, C = C,
        list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
        mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
        mu_q_nu_beta = mu_q_nu_beta, Z = Z,
        mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
        mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
        mu_q_nu_psi = mu_q_nu_psi,
        mu_q_a = mu_q_a,
        mu_q_b_specific = mu_q_b_specific,
        term_b_specific = term_b_specific,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
        inv_Sigma_beta = inv_Sigma_beta,
        S = S, n_s = n_s, p = p, d = d,
        L_f = L_f, L_s = L_s, M_s = M_s,
        K = K, K_total = K_total,
        c_val = c_val, n_cpus = n_cpus)
      mu_q_nu_psi <- res_psi$mu_q_nu_psi
      Sigma_q_nu_psi <- res_psi$Sigma_q_nu_psi
      inv_Sigma_q_nu_psi <- res_psi$inv_Sigma_q_nu_psi
    }

    # ------ Block 5: Update q(zeta) ------
    res_zeta <- update_zeta(
      Y = Y, C = C,
      list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
      mu_q_nu_mu = mu_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi, Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
      mu_q_a = mu_q_a, term_a = term_a,
      mu_q_b_specific = mu_q_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f,
      c_val = c_val, n_cpus = n_cpus)
    mu_q_zeta <- res_zeta$mu_q_zeta
    Sigma_q_zeta <- res_zeta$Sigma_q_zeta
    tr_qi_shared <- res_zeta$tr_qi_shared

    # ------ Block 6: Update q(xi) ------
    if (L_s > 0) {
      res_xi <- update_xi(
        Y = Y, C = C,
        list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
        mu_q_nu_mu = mu_q_nu_mu,
        mu_q_nu_beta = mu_q_nu_beta, Z = Z,
        mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
        mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
        mu_q_nu_psi = mu_q_nu_psi, Sigma_q_nu_psi = Sigma_q_nu_psi,
        mu_q_a = mu_q_a,
        mu_q_b_specific = mu_q_b_specific,
        term_b_specific = term_b_specific,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        S = S, n_s = n_s, p = p, d = d,
        L_f = L_f, L_s = L_s, M_s = M_s,
        c_val = c_val, n_cpus = n_cpus)
      mu_q_xi <- res_xi$mu_q_xi
      Sigma_q_xi <- res_xi$Sigma_q_xi
      tr_xi_specific <- res_xi$tr_xi_specific
    }

    # ------ Block 7-8: Update all variance parameters ------
    res_var <- update_all_variances(
      Y = Y, C = C, list_cp_C = list_cp_C,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Sigma_q_nu_beta = Sigma_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi, Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
      mu_q_nu_psi = mu_q_nu_psi, Sigma_q_nu_psi = Sigma_q_nu_psi,
      mu_q_a = mu_q_a, term_a = term_a,
      mu_q_b_specific = mu_q_b_specific, term_b_specific = term_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_a_eps = mu_q_recip_a_eps,
      mu_q_recip_sigsq_mu = mu_q_recip_sigsq_mu,
      mu_q_recip_a_mu = mu_q_recip_a_mu,
      mu_q_recip_sigsq_beta = mu_q_recip_sigsq_beta,
      mu_q_recip_a_beta = mu_q_recip_a_beta,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      mu_q_recip_a_phi = mu_q_recip_a_phi,
      mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
      mu_q_recip_a_psi = mu_q_recip_a_psi,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
      c_val = c_val, n_cpus = n_cpus)

    mu_q_recip_sigsq_eps <- res_var$mu_q_recip_sigsq_eps
    mu_q_recip_a_eps <- res_var$mu_q_recip_a_eps
    mu_q_recip_sigsq_mu <- res_var$mu_q_recip_sigsq_mu
    mu_q_recip_a_mu <- res_var$mu_q_recip_a_mu
    mu_q_recip_sigsq_phi <- res_var$mu_q_recip_sigsq_phi
    mu_q_recip_a_phi <- res_var$mu_q_recip_a_phi
    if (!is.null(res_var$mu_q_recip_sigsq_beta)) {
      mu_q_recip_sigsq_beta <- res_var$mu_q_recip_sigsq_beta
      mu_q_recip_a_beta <- res_var$mu_q_recip_a_beta
    }
    if (!is.null(res_var$mu_q_recip_sigsq_psi)) {
      mu_q_recip_sigsq_psi <- res_var$mu_q_recip_sigsq_psi
      mu_q_recip_a_psi <- res_var$mu_q_recip_a_psi
    }
    kappa_q_a <- res_var$kappa_q_a
    mu_q_log_sigsq_eps <- res_var$mu_q_log_sigsq_eps
    mu_q_log_sigsq_mu <- res_var$mu_q_log_sigsq_mu
    lambda_q_sigsq_eps <- res_var$lambda_q_sigsq_eps
    lambda_q_sigsq_mu <- res_var$lambda_q_sigsq_mu
    kappa_q_sigsq_eps <- res_var$kappa_q_sigsq_eps
    kappa_q_sigsq_mu <- res_var$kappa_q_sigsq_mu
    kappa_q_sigsq_phi <- res_var$kappa_q_sigsq_phi
    lambda_q_sigsq_phi <- res_var$lambda_q_sigsq_phi

    # ------ Block 9: Update q(omega) ------
    if (!is.null(mu_q_gamma_b)) {
      res_omega <- update_omega(
        mu_q_gamma_a = mu_q_gamma_a,
        mu_q_gamma_b = mu_q_gamma_b,
        p = p, S = S,
        c_0 = c_0, d_0 = d_0,
        c_val = c_val, bool_var_spec_prob = bool_var_spec_prob)
    } else {
      res_omega <- update_omega(
        mu_q_gamma_a = mu_q_gamma_a, mu_q_gamma_b = NULL,
        p = p, S = S,
        c_0 = c_0, d_0 = d_0,
        c_val = c_val, bool_var_spec_prob = bool_var_spec_prob)
    }
    c_1_omega_a <- res_omega$c_1_omega_a
    d_1_omega_a <- res_omega$d_1_omega_a
    mu_q_log_omega_a <- res_omega$mu_q_log_omega_a
    mu_q_log_1_omega_a <- res_omega$mu_q_log_1_omega_a
    if (!is.null(res_omega$c_1_omega_b)) {
      c_1_omega_b <- res_omega$c_1_omega_b
      d_1_omega_b <- res_omega$d_1_omega_b
      mu_q_log_omega_b <- res_omega$mu_q_log_omega_b
      mu_q_log_1_omega_b <- res_omega$mu_q_log_1_omega_b
    }

    # ------ Block 10: Update q(a_loadings) ------
    res_a <- update_a_loadings(
      Y = Y, C = C,
      list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
      mu_q_nu_mu = mu_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
      mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
      mu_q_a = mu_q_a, term_a = term_a, mu_q_gamma_a = mu_q_gamma_a,
      mu_q_normal_a = mu_q_normal_a, Sigma_q_normal_a = Sigma_q_normal_a,
      mu_q_b_specific = mu_q_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_log_omega_a = mu_q_log_omega_a,
      mu_q_log_1_omega_a = mu_q_log_1_omega_a,
      tr_qi_shared = tr_qi_shared,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, K_total = K_total,
      c_val = c_val, n_cpus = n_cpus)
    mu_q_normal_a <- res_a$mu_q_normal_a
    Sigma_q_normal_a <- res_a$Sigma_q_normal_a
    mu_q_gamma_a <- res_a$mu_q_gamma_a
    mu_q_a <- res_a$mu_q_a
    term_a <- res_a$term_a

    # ------ Block 11: Update q(b_loadings) ------
    if (L_s > 0) {
      res_b <- update_b_loadings(
        Y = Y, C = C,
        list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
        mu_q_nu_mu = mu_q_nu_mu,
        mu_q_nu_beta = mu_q_nu_beta, Z = Z,
        mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
        mu_q_xi = mu_q_xi, mu_q_nu_psi = mu_q_nu_psi,
        mu_q_a = mu_q_a,
        mu_q_b_specific = mu_q_b_specific,
        term_b_specific = term_b_specific,
        mu_q_gamma_b = mu_q_gamma_b,
        mu_q_normal_b = mu_q_normal_b,
        Sigma_q_normal_b = Sigma_q_normal_b,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        mu_q_log_omega_b = mu_q_log_omega_b,
        mu_q_log_1_omega_b = mu_q_log_1_omega_b,
        tr_xi_specific = tr_xi_specific,
        S = S, n_s = n_s, p = p, d = d,
        L_f = L_f, L_s = L_s, K_total = K_total,
        c_val = c_val, n_cpus = n_cpus)
      mu_q_normal_b <- res_b$mu_q_normal_b
      Sigma_q_normal_b <- res_b$Sigma_q_normal_b
      mu_q_gamma_b <- res_b$mu_q_gamma_b
      mu_q_b_specific <- res_b$mu_q_b_specific
      term_b_specific <- res_b$term_b_specific
    }

    # ------ ELBO computation ------
    n_obs_per_indiv <- length(Y[[1]][[1]][[1]])
    ELBO_iter <- compute_elbo_multi(
      Y = Y, C = C, list_cp_C = list_cp_C,
      n_obs_per_indiv = n_obs_per_indiv,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Sigma_q_nu_beta = Sigma_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi, Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
      mu_q_nu_psi = mu_q_nu_psi, Sigma_q_nu_psi = Sigma_q_nu_psi,
      mu_q_a = mu_q_a, mu_q_normal_a = mu_q_normal_a,
      Sigma_q_normal_a = Sigma_q_normal_a, mu_q_gamma_a = mu_q_gamma_a,
      mu_q_b_specific = mu_q_b_specific,
      mu_q_normal_b = mu_q_normal_b, Sigma_q_normal_b = Sigma_q_normal_b,
      mu_q_gamma_b = mu_q_gamma_b,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_a_eps = mu_q_recip_a_eps,
      mu_q_recip_sigsq_mu = mu_q_recip_sigsq_mu,
      mu_q_recip_a_mu = mu_q_recip_a_mu,
      mu_q_recip_sigsq_beta = mu_q_recip_sigsq_beta,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      mu_q_recip_a_phi = mu_q_recip_a_phi,
      mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
      mu_q_log_sigsq_eps = mu_q_log_sigsq_eps,
      mu_q_log_sigsq_mu = mu_q_log_sigsq_mu,
      kappa_q_sigsq_eps = kappa_q_sigsq_eps,
      kappa_q_sigsq_mu = kappa_q_sigsq_mu,
      lambda_q_sigsq_eps = lambda_q_sigsq_eps,
      lambda_q_sigsq_mu = lambda_q_sigsq_mu,
      mu_q_log_omega_a = mu_q_log_omega_a,
      mu_q_log_1_omega_a = mu_q_log_1_omega_a,
      inv_Sigma_beta = inv_Sigma_beta,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
      c_val = c_val)

    ELBO <- c(ELBO, ELBO_iter)

    if (verbose) {
      cat(paste0("ELBO = ", format(ELBO_iter), "\n\n"))
    }

    # Convergence check
    if (i_iter > i_iter_init) {
      l_ELBO <- length(ELBO)
      ELBO_diff <- ELBO[l_ELBO] - ELBO[l_ELBO - 1]

      if (ELBO_diff < -eps) {
        warning(paste0("ELBO not increasing monotonically. Diff: ", ELBO_diff))
      }

      rel_converged <- (abs(max(ELBO[l_ELBO - 1], ELBO[l_ELBO]) /
                          min(ELBO[l_ELBO - 1], ELBO[l_ELBO]) - 1) < tol_rel)
      abs_converged <- (abs(ELBO_diff) < tol_abs)

      if (rel_converged | abs_converged) {
        if (verbose) {
          cat(paste0("Convergence obtained after ", i_iter, " iterations.\n"))
        }
        break
      } else if (i_iter == maxit) {
        warning(paste0("Max iterations reached before convergence."))
      }
    }
  } # end main loop

  # ===================================================================
  #  RESULTS ASSEMBLY
  # ===================================================================

  res <- create_named_list(
    S, n_s, p, d,
    L_f, L_s, M_f, M_s,
    K, K_total,
    mu_q_nu_mu, Sigma_q_nu_mu,
    mu_q_nu_beta, Sigma_q_nu_beta,
    mu_q_nu_phi, Sigma_q_nu_phi,
    mu_q_nu_psi, Sigma_q_nu_psi,
    mu_q_zeta, Sigma_q_zeta,
    mu_q_xi, Sigma_q_xi,
    mu_q_normal_a, Sigma_q_normal_a, mu_q_gamma_a, mu_q_a, term_a,
    mu_q_normal_b, Sigma_q_normal_b, mu_q_gamma_b, mu_q_b_specific, term_b_specific,
    mu_q_recip_sigsq_eps, mu_q_recip_a_eps,
    mu_q_recip_sigsq_mu, mu_q_recip_a_mu,
    mu_q_recip_sigsq_beta, mu_q_recip_a_beta,
    mu_q_recip_sigsq_phi, mu_q_recip_a_phi,
    mu_q_recip_sigsq_psi, mu_q_recip_a_psi,
    c_1_omega_a, d_1_omega_a, mu_q_log_omega_a, mu_q_log_1_omega_a,
    c_1_omega_b, d_1_omega_b, mu_q_log_omega_b, mu_q_log_1_omega_b,
    ELBO, ELBO_iter, i_iter,
    mean_mean_across_subjects, sd_mean_across_subjects,
    time_g, C_g, n_g,
    list_cp_C, list_cp_C_Y, list_cp_Y, sum_list_cp_C,
    inv_Sigma_q_nu_mu, inv_Sigma_q_nu_phi, inv_Sigma_q_nu_psi
  )

  res
}


#' Compute multi-study ELBO for convergence monitoring
#'
#' @keywords internal
compute_elbo_multi <- function(Y, C, list_cp_C,
                                n_obs_per_indiv,
                                mu_q_nu_mu, Sigma_q_nu_mu,
                                mu_q_nu_beta, Sigma_q_nu_beta, Z,
                                mu_q_zeta, Sigma_q_zeta,
                                mu_q_nu_phi, Sigma_q_nu_phi,
                                mu_q_xi, Sigma_q_xi,
                                mu_q_nu_psi, Sigma_q_nu_psi,
                                mu_q_a, mu_q_normal_a,
                                Sigma_q_normal_a, mu_q_gamma_a,
                                mu_q_b_specific,
                                mu_q_normal_b, Sigma_q_normal_b,
                                mu_q_gamma_b,
                                mu_q_recip_sigsq_eps, mu_q_recip_a_eps,
                                mu_q_recip_sigsq_mu, mu_q_recip_a_mu,
                                mu_q_recip_sigsq_beta,
                                mu_q_recip_sigsq_phi, mu_q_recip_a_phi,
                                mu_q_recip_sigsq_psi,
                                mu_q_log_sigsq_eps, mu_q_log_sigsq_mu,
                                kappa_q_sigsq_eps, kappa_q_sigsq_mu,
                                lambda_q_sigsq_eps, lambda_q_sigsq_mu,
                                mu_q_log_omega_a, mu_q_log_1_omega_a,
                                inv_Sigma_beta,
                                S, n_s, p, d,
                                L_f, L_s, M_f, M_s, K,
                                c_val = 1) {

  eps_elbo <- 1e-11
  nobs <- n_obs_per_indiv

  # ---- Data fitting ----
  elbo_y <- 0
  for (s in 1:S) {
    for (j in 1:p) {
      n_obs_sj <- n_s[s] * nobs
      elbo_y <- elbo_y - 0.5 * n_obs_sj * (log(2 * pi) + mu_q_log_sigsq_eps[s, j]) -
        mu_q_recip_sigsq_eps[s, j] * (lambda_q_sigsq_eps[s, j] - mu_q_recip_a_eps[s, j])
    }
  }

  # ---- Mean function: Gaussian entropy + prior + sigma^2 entropy ----
  elbo_mu <- 0
  for (s in 1:S) {
    for (j in 1:p) {
      det_val <- determinant(Sigma_q_nu_mu[[s]][[j]], logarithm = TRUE)
      log_det <- det_val$modulus * det_val$sign
      nu_pen <- mu_q_nu_mu[[s]][[j]][-c(1:2)]
      Sigma_pen <- Sigma_q_nu_mu[[s]][[j]][-c(1:2), -c(1:2)]

      elbo_mu <- elbo_mu +
        0.5 * log_det + 0.5 * K -
        0.5 * as.numeric(crossprod(nu_pen)) * mu_q_recip_sigsq_mu[s, j] -
        0.5 * tr(Sigma_pen) * mu_q_recip_sigsq_mu[s, j] +
        K / 2
    }
  }

  # ---- Sigma^2 entropy ---
  elbo_sigma <- 0
  for (s in 1:S) {
    for (j in 1:p) {
      ke <- kappa_q_sigsq_eps[s, j]
      le <- lambda_q_sigsq_eps[s, j]
      elbo_sigma <- elbo_sigma +
        (ke - 0.5) * mu_q_log_sigsq_eps[s, j] -
        ke * log(le) - lgamma(0.5) + lgamma(ke)

      elbo_sigma <- elbo_sigma -
        mu_q_recip_sigsq_eps[s, j] * (mu_q_recip_a_eps[s, j] - lambda_q_sigsq_eps[s, j])

      km <- kappa_q_sigsq_mu
      lm <- lambda_q_sigsq_mu[s, j]
      elbo_sigma <- elbo_sigma +
        (km - 0.5) * mu_q_log_sigsq_mu[s, j] -
        km * log(lm) - lgamma(0.5) + lgamma(km)
    }
  }

  # ---- Loadings + PPI ----
  elbo_loadings <- 0
  for (l in seq_len(L_f)) {
    for (j in 1:p) {
      ppi <- mu_q_gamma_a[j, l]
      sigma2 <- max(Sigma_q_normal_a[j, l], eps_elbo)
      mu_val <- mu_q_normal_a[j, l]
      elbo_loadings <- elbo_loadings +
        0.5 * ppi * (log(sigma2) + 1) -
        0.5 * ppi * (mu_val^2 + sigma2) +
        ppi * mu_q_log_omega_a[l] +
        (1 - ppi) * mu_q_log_1_omega_a[l] -
        ppi * log(max(ppi, eps_elbo)) -
        (1 - ppi) * log(max(1 - ppi, eps_elbo))
    }
  }
  if (L_s > 0 && !is.null(mu_q_gamma_b)) {
    for (s in 1:S) {
      for (l in seq_len(L_s)) {
        for (j in 1:p) {
          ppi <- mu_q_gamma_b[[s]][j, l]
          sigma2 <- max(Sigma_q_normal_b[[s]][j, l], eps_elbo)
          mu_val <- mu_q_normal_b[[s]][j, l]
          elbo_loadings <- elbo_loadings +
            0.5 * ppi * (log(sigma2) + 1) -
            0.5 * ppi * (mu_val^2 + sigma2)
        }
      }
    }
  }

  # ---- Zeta score entropy ----
  elbo_zeta <- 0
  for (s in 1:S) {
    for (l in seq_len(L_f)) {
      for (i in 1:n_s[s]) {
        det_val <- determinant(Sigma_q_zeta[[s]][[l]][[i]], logarithm=TRUE)
        log_det <- det_val$modulus * det_val$sign
        elbo_zeta <- elbo_zeta +
          0.5 * log_det -
          0.5 * sum(mu_q_zeta[[s]][[l]][i, ]^2) -
          0.5 * tr(Sigma_q_zeta[[s]][[l]][[i]])
      }
    }
  }

  # ---- nu_phi prior + entropy ----
  elbo_phi <- 0
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    for (m in 1:M_l) {
      det_val <- determinant(Sigma_q_nu_phi[[l]][[m]], logarithm=TRUE)
      log_det <- det_val$modulus * det_val$sign
      nu_pen <- mu_q_nu_phi[[l]][-c(1:2), m]
      Sigma_pen <- Sigma_q_nu_phi[[l]][[m]][-c(1:2), -c(1:2)]
      elbo_phi <- elbo_phi +
        0.5 * log_det + 0.5 * K -
        0.5 * sum(nu_pen^2) * mu_q_recip_sigsq_phi[[l]][m] -
        0.5 * tr(Sigma_pen) * mu_q_recip_sigsq_phi[[l]][m] + K/2
    }
  }

  # ---- nu_psi prior + entropy ----
  elbo_psi <- 0
  if (L_s > 0 && !is.null(mu_q_nu_psi)) {
    for (s in 1:S) {
      for (l in seq_len(L_s)) {
        M_sl <- M_s[[s]][l]
        for (m in 1:M_sl) {
          det_val <- determinant(Sigma_q_nu_psi[[s]][[l]][[m]], logarithm=TRUE)
          log_det <- det_val$modulus * det_val$sign
          nu_pen <- mu_q_nu_psi[[s]][[l]][-c(1:2), m]
          Sigma_pen <- Sigma_q_nu_psi[[s]][[l]][[m]][-c(1:2), -c(1:2)]
          elbo_psi <- elbo_psi +
            0.5 * log_det + 0.5 * K -
            0.5 * sum(nu_pen^2) * mu_q_recip_sigsq_psi[[s]][[l]][m] -
            0.5 * tr(Sigma_pen) * mu_q_recip_sigsq_psi[[s]][[l]][m] + K/2
        }
      }
    }
  }

  # ---- Xi score entropy ----
  elbo_xi <- 0
  if (L_s > 0 && !is.null(mu_q_xi)) {
    for (s in 1:S) {
      for (l in seq_len(L_s)) {
        for (i in 1:n_s[s]) {
          det_val <- determinant(Sigma_q_xi[[s]][[l]][[i]], logarithm=TRUE)
          log_det <- det_val$modulus * det_val$sign
          elbo_xi <- elbo_xi +
            0.5 * log_det -
            0.5 * sum(mu_q_xi[[s]][[l]][i, ]^2) -
            0.5 * tr(Sigma_q_xi[[s]][[l]][[i]])
        }
      }
    }
  }

  # ---- nu_beta prior + entropy ----
  elbo_beta <- 0
  if (d > 0 && !is.null(mu_q_nu_beta)) {
    for (j in 1:p) {
      for (r in 1:d) {
        det_val <- determinant(Sigma_q_nu_beta[[j]][[r]], logarithm=TRUE)
        log_det <- det_val$modulus * det_val$sign
        nu_pen <- mu_q_nu_beta[[j]][[r]][-c(1:2)]
        Sigma_pen <- Sigma_q_nu_beta[[j]][[r]][-c(1:2), -c(1:2)]
        elbo_beta <- elbo_beta +
          0.5 * log_det + 0.5 * K -
          0.5 * sum(nu_pen^2) * mu_q_recip_sigsq_beta[j, r] -
          0.5 * tr(Sigma_pen) * mu_q_recip_sigsq_beta[j, r] + K/2
      }
    }
  }

  ELBO <- elbo_y + elbo_mu + elbo_sigma + elbo_loadings +
          elbo_zeta + elbo_phi + elbo_psi + elbo_xi + elbo_beta

  return(ELBO)
}
