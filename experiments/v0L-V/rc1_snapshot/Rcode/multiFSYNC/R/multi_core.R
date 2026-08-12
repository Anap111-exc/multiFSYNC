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
#   L_s = number of study-specific factors (scalar or one count per study)
#   M_f = vector length L_f: M_f[l] = FPCA truncation for shared factor l
#   M_s = list of S vectors: M_s[[s]][l] = FPCA truncation per study-factor
#
# Theory/computation convention (Scheme 1; see ../NOTATION.md):
#   theoretical FPCA:  (zeta, xi, phi, psi, lambda)
#   CAVI computation:  (eta, chi, theta, kappa), with Var(eta)=Var(chi)=1
#   legacy R objects:   zeta -> eta, xi -> chi,
#                       nu_phi -> coefficients of theta,
#                       nu_psi -> coefficients of kappa
# Public object names are kept for backward compatibility.  They must not be
# interpreted as theoretical FPCA scores/eigenfunctions until the
# orthonormalise_multi() post-processing step has been applied.
#
# Indexing convention: study [[s]] is ALWAYS the outermost dimension.
# [[s]] → study, [[i]] → individual, [[j]] → variable, [[l]] → factor,
# [[m]] → FPCA component, [[r]] → covariate
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

# Maximum absolute and scaled changes over nested lists of parameter means.
# This remains a useful fallback and an independent diagnostic even when the
# complete ordinary ELBO is available.
.parameter_change <- function(old, new) {
  if (is.null(old) && is.null(new)) return(c(abs = 0, rel = 0))
  if (is.list(old) || is.list(new)) {
    if (!is.list(old) || !is.list(new) || length(old) != length(new)) {
      return(c(abs = Inf, rel = Inf))
    }
    if (length(old) == 0L) return(c(abs = 0, rel = 0))
    vals <- Map(.parameter_change, old, new)
    return(c(abs = max(vapply(vals, `[[`, numeric(1), "abs")),
             rel = max(vapply(vals, `[[`, numeric(1), "rel"))))
  }
  if (length(old) != length(new)) return(c(abs = Inf, rel = Inf))
  if (length(old) == 0L) return(c(abs = 0, rel = 0))
  abs_change <- max(abs(as.numeric(new) - as.numeric(old)), na.rm = TRUE)
  scale <- 1 + max(abs(as.numeric(old)), na.rm = TRUE)
  c(abs = abs_change, rel = abs_change / scale)
}

# Compact summaries used by the parameter-based stopping rule.  Monitoring
# diagonals and traces captures covariance movement without copying every
# off-diagonal entry of the variational covariance arrays.
.posterior_cov_summary <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) return(lapply(x, .posterior_cov_summary))
  if (is.matrix(x)) {
    return(c(trace = sum(diag(x)), diagonal = diag(x)))
  }
  as.numeric(x)
}

#' Multi-Study Bayesian Functional Factor Model
#'
#' @details The CAVI implementation uses standardised working scores and
#'   unconstrained spline time functions.  For backward compatibility the
#'   corresponding R objects retain the historical names `zeta`, `xi`,
#'   `nu_phi`, and `nu_psi`.  Canonical FPCA scores, eigenfunctions, and
#'   eigenvalues are produced by the post-processing step; see `NOTATION.md`.
#'
#' @param Y List of length S, where Y[[s]] is a list of n_s[s] individuals,
#'   each being a list of p variable vectors (length n_i observations).
#' @param Z Optional list of length S, Z[[s]] is an n_s[s] x d covariate matrix.
#' @param time_obs List of length S, time_obs[[s]][[i]] is the observation
#'   time vector for individual i in study s.
#' @param L_f Number of shared latent factors (scalar).
#' @param L_s Number of study-specific latent factors. A scalar is recycled to
#'   all studies; a length-S vector allows different counts by study.
#' @param M_f Vector of length L_f, FPCA truncation per shared factor.
#' @param M_s List of S vectors, FPCA truncation per study-specific factor.
#' @param K Number of O'Sullivan spline functions (default NULL: auto).
#' @param anneal Annealing parameters: c(type, T_max, grid_size) or NULL.
#' @param list_hyper Hyperparameter list from set_hyper() or NULL for defaults.
#' @param n_g Dense grid size.
#' @param time_g Dense grid vector (overrides n_g).
#' @param tol_abs Absolute convergence tolerance.
#' @param tol_rel Relative convergence tolerance.
#' @param convergence_rule Stopping rule. `"elbo"` uses the complete ordinary
#'   ELBO at T = 1 only when `lambda_orth = 0` and no adaptive Cholesky jitter
#'   is used during the T = 1 phase; `"practical"` additionally requires
#'   short- and long-window stability of fitted means, expected RSS, and PPI;
#'   `"parameters"` is the strict-coordinate fallback otherwise.
#' @param practical_control Optional named list overriding the practical-rule
#'   defaults returned by the internal `.default_practical_control()`. Recent
#'   20-sweep changes use strict fitted/RSS/PPI tolerances; cumulative
#'   60-sweep changes have separately calibrated long-window tolerances.
#' @param lambda_orth Optional time-function separation regulariser.  The
#'   default 0 disables it; positive values are not identification constraints.
#' @param maxit Maximum number of iterations.
#' @param n_cpus Number of CPU cores for parallel execution.
#' @param verbose Print progress messages.
#' @param seed Random seed for reproducibility.
#' @param bool_scale Whether to standardise variables.
#' @param bool_var_spec_prob Use variable-factor-level rather than factor-level
#'   omega probabilities.
#' @param d_0 Optional second Beta-prior shape for omega; defaults to `p`.
#' @param initialization Iteration-zero parameter construction. `"random"`
#'   preserves the original seeded initialization; `"residual_fpca"` uses only
#'   the working-scale observed data to initialize mean/covariate curves and
#'   shared/specific factor means before running the unchanged CAVI updates.
#' @param initialization_control Optional named list overriding the
#'   residual-FPCA controls returned by the internal
#'   `.default_initialization_control()`: `grid_size`, `perturb_sd`, and the
#'   relative `rank_tol`.
#' @param continuation_state Optional validated internal state used to continue
#'   a dimension-reduced model at `T = 1`.  `NULL` preserves the ordinary
#'   initialization path.  This argument is intended for reproducible staged
#'   model-selection experiments and does not alter any CAVI update.
#'
#' @return A fitted object containing working-scale variational parameters,
#'   original-scale scientific summaries, complete ordinary-ELBO diagnostics,
#'   practical stopping diagnostics, and truth-free initialization diagnostics.
#'   The latter include the method used, explained residual fraction, per-block
#'   fallback reasons, construction time, and an independence identifier for
#'   multi-start support counting.
#'
#' @importFrom splines spline.des
#' @importFrom stats median qnorm quantile rnorm runif sd setNames
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
    bool_scale = TRUE, bool_var_spec_prob = FALSE,
    d_0 = NULL,
    convergence_rule = c("parameters", "elbo", "practical"),
    lambda_orth = 0,
    practical_control = NULL,
    initialization = c("random", "residual_fpca"),
    initialization_control = NULL,
    continuation_state = NULL
) {

  scale_trace_control <- getOption(
    .multiFSYNC_scale_trace_option, NULL
  )
  if (!is.null(scale_trace_control)) {
    scale_trace_control <-
      .validate_scale_trace_control(scale_trace_control)
  }
  convergence_rule_requested <- match.arg(convergence_rule)
  practical_control <- .validate_practical_control(practical_control)
  initialization <- match.arg(initialization)
  initialization_control <-
    .validate_initialization_control(initialization_control)
  if (length(lambda_orth) != 1L || !is.finite(lambda_orth) || lambda_orth < 0) {
    stop("lambda_orth must be one finite non-negative number.")
  }
  if (!is.null(continuation_state)) {
    if (!inherits(continuation_state, "multiFSYNC_continuation_state")) {
      stop("continuation_state must be NULL or a validated internal continuation state.")
    }
    if (!is.null(anneal)) {
      stop("continuation_state is supported only with anneal = NULL (T = 1).")
    }
    if (lambda_orth != 0) {
      stop("continuation_state currently requires lambda_orth = 0.")
    }
    if (!identical(initialization, "random")) {
      stop("continuation_state cannot be combined with a separate data-driven initialization.")
    }
  }
  if (lambda_orth > 0 &&
      convergence_rule_requested %in% c("elbo", "practical")) {
    warning("lambda_orth > 0 has no validated global ELBO objective; ",
            "downgrading the requested objective-based convergence rule ",
            "to parameter stopping. ",
            "The ordinary-model ELBO is retained as a diagnostic only.")
    convergence_rule <- "parameters"
  } else {
    convergence_rule <- convergence_rule_requested
  }
  elbo_role <- if (lambda_orth == 0) "objective" else "diagnostic_only"
  if (length(tol_abs) != 1L || !is.finite(tol_abs) || tol_abs < 0) {
    stop("tol_abs must be one finite non-negative number.")
  }
  if (length(tol_rel) != 1L || !is.finite(tol_rel) || tol_rel < 0) {
    stop("tol_rel must be one finite non-negative number.")
  }

  if (length(maxit) != 1L || !is.finite(maxit) || !is_int(maxit) || maxit < 1L) {
    stop("maxit must be one positive integer.")
  }
  maxit <- as.integer(maxit)
  dense_gate_sweeps <- .driver_diagnostic_dense_gate_sweeps()
  if (dense_gate_sweeps >= maxit) {
    stop(
      "The private dense-gate diagnostic requires maxit greater than ",
      "control$dense_gate_sweeps so the original spike-and-slab update is ",
      "released for at least one sweep."
    )
  }
  if (length(n_cpus) != 1L || !is.finite(n_cpus) || !is_int(n_cpus) || n_cpus < 1L) {
    stop("n_cpus must be one positive integer.")
  }
  n_cpus <- as.integer(n_cpus)
  if (.Platform$OS.type == "windows" && n_cpus > 1L) {
    warning("parallel::mclapply does not support n_cpus > 1 on Windows; using n_cpus = 1. Parallelise simulation replicates at the script level instead.")
    n_cpus <- 1L
  }

  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed) ||
        !is_int(seed) || seed < 0 || seed > .Machine$integer.max) {
      stop("seed must be NULL or one non-negative integer.")
    }
    seed <- as.integer(seed)
    if (verbose) cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

  check_annealing(anneal, verbose)
  if (!is.null(anneal) && maxit <= as.integer(anneal[3]) + 2L) {
    stop("With annealing, maxit must be greater than anneal[3] + 2 so that iterations at T = 1 remain available for convergence assessment.")
  }
  if (identical(convergence_rule, "practical")) {
    annealing_sweeps_planned <- if (is.null(anneal)) {
      0L
    } else {
      as.integer(anneal[3]) - 1L
    }
    practical_budget_required <- annealing_sweeps_planned +
      practical_control$max_t1
    if (maxit < practical_budget_required) {
      stop(
        "For convergence_rule='practical', maxit must be at least ",
        practical_budget_required, " (= planned annealing sweeps + ",
        "practical_control$max_t1) so the fit can either converge or be ",
        "classified as a slow case."
      )
    }
    if (any(practical_control$checkpoints >
            practical_control$max_t1)) {
      stop(
        "For convergence_rule='practical', every practical checkpoint ",
        "must be no greater than practical_control$max_t1. Use a separate ",
        "forced reference fit for later checkpoints."
      )
    }
  }

  # ---- Input validation ----
  if (!is.list(Y) || length(Y) == 0L) stop("Y must be a non-empty list of studies.")
  S <- length(Y)
  n_s <- vapply(Y, length, integer(1))
  if (any(n_s < 1L)) stop("Every study must contain at least one individual.")

  p <- length(Y[[1]][[1]])
  if (p < 1L) stop("Every individual must contain at least one functional variable.")
  if (is.null(time_obs) || !is.list(time_obs) || length(time_obs) != S) {
    stop("time_obs must be a list with one element per study.")
  }
  var_names <- names(Y[[1]][[1]])
  if (is.null(var_names)) var_names <- paste0("variable_", seq_len(p))
  for (s in seq_len(S)) {
    if (!is.list(Y[[s]]) || length(time_obs[[s]]) != n_s[s]) {
      stop("Y and time_obs must contain the same individuals in every study.")
    }
    for (i in seq_len(n_s[s])) {
      if (!is.list(Y[[s]][[i]]) || length(Y[[s]][[i]]) != p) {
        stop("Every individual must contain the same p functional variables.")
      }
      ti <- time_obs[[s]][[i]]
      if (!is.numeric(ti) || length(ti) < 2L || any(!is.finite(ti)) ||
          any(ti < 0 | ti > 1)) {
        stop("Each observation-time vector must contain at least two finite values in [0,1].")
      }
      for (j in seq_len(p)) {
        yij <- Y[[s]][[i]][[j]]
        if (!is.numeric(yij) || length(yij) != length(ti) || any(!is.finite(yij))) {
          stop("Each response vector must be finite and aligned with its individual's time_obs vector.")
        }
      }
    }
  }

  check_structure(L_f, "vector", "numeric", 1)
  if (!is_int(L_f) || L_f < 0) stop("L_f must be a non-negative integer.")
  L_f <- as.integer(L_f)
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)

  if (L_f > 0) {
    check_structure(M_f, "vector", "numeric")
    if (length(M_f) != L_f || any(!is_int(M_f)) || any(M_f < 1L)) {
      stop("M_f must contain one positive integer for each shared factor.")
    }
    M_f <- as.integer(M_f)
  } else if (length(M_f) != 0L) {
    stop("M_f must be empty when L_f = 0.")
  }

  if (!is.list(M_s) || length(M_s) != S) {
    stop("M_s must be a list with one element per study.")
  }
  if (has_specific) {
    for (s in seq_len(S)) {
      L_ss <- L_s[[s]]
      if (L_ss == 0L) {
        if (length(M_s[[s]]) != 0L) {
          stop("M_s[[", s, "]] must be empty because L_s[", s, "] = 0.")
        }
        M_s[[s]] <- integer(0)
        next
      }
      check_structure(M_s[[s]], "vector", "numeric")
      if (length(M_s[[s]]) != L_ss || any(!is_int(M_s[[s]])) ||
          any(M_s[[s]] < 1L)) {
        stop("M_s[[", s, "]] must contain one positive integer for each of that study's specific factors.")
      }
      M_s[[s]] <- as.integer(M_s[[s]])
    }
  } else if (any(lengths(M_s) != 0L)) {
    stop("Every M_s[[s]] must be empty when L_s = 0.")
  }

  # ---- Covariates ----
  if (!is.null(Z)) {
    if (!is.list(Z) || length(Z) != S || !is.matrix(Z[[1]])) {
      stop("Z must be a list of numeric matrices, one per study.")
    }
    d <- ncol(Z[[1]])
    if (d < 1L) stop("Use Z = NULL when there are no covariates.")
    for (s in seq_len(S)) {
      if (!is.matrix(Z[[s]]) || !is.numeric(Z[[s]]) ||
          nrow(Z[[s]]) != n_s[s] || ncol(Z[[s]]) != d ||
          any(!is.finite(Z[[s]]))) {
        stop("Each Z[[s]] must be a finite n_s[s] by d numeric matrix.")
      }
    }
    Z_within <- do.call(rbind, lapply(Z, function(z) {
      sweep(z, 2, colMeans(z), "-")
    }))
    if (qr(Z_within)$rank < d) {
      stop("Covariates are not full rank after removing study means; beta(t) is confounded with the study-specific mean functions.")
    }
  } else {
    d <- 0
  }

  # ---- Hyperparameters ----
  if (is.null(list_hyper)) list_hyper <- set_hyper(d_0 = d_0)
  if (!inherits(list_hyper, "hyper")) stop("list_hyper must be created by set_hyper().")
  if (!is.null(d_0)) list_hyper$d_0 <- d_0
  if (is.null(list_hyper$d_0)) list_hyper$d_0 <- p
  check_positive(list_hyper$c_0)
  check_positive(list_hyper$d_0)

  # ---- Flat observation grid ----
  time_obs_flat <- unlist(lapply(1:S, function(s) time_obs[[s]]), recursive = FALSE)

  if (is.null(K)) {
    n_times <- sapply(time_obs_flat, length)
    K <- max(7, min(floor(min(n_times) / 4), 40))
    if (verbose) cat(paste0("K set to ", K, "\n"))
  } else {
    check_natural(K)
    if (K < 2L) stop("K must be at least 2 for the O'Sullivan spline construction.")
  }

  grid_obj <- get_grid_objects(time_obs_flat, K, n_g = n_g, time_g = time_g,
                                format_univ = TRUE)
  C_flat <- grid_obj$C
  n_g <- grid_obj$n_g
  time_g <- grid_obj$time_g
  C_g <- grid_obj$C_g
  if (ncol(C_g) != K + 2L) {
    stop(sprintf("The O'Sullivan design has %d columns, but K + 2 = %d. Check K and the knot construction.",
                 ncol(C_g), K + 2L))
  }

  basis_rank <- qr(C_g)$rank
  if (L_f > 0L && any(M_f > basis_rank)) {
    stop(sprintf("Each M_f[l] must not exceed rank(C_g) = %d.", basis_rank))
  }
  if (has_specific && any(unlist(M_s, use.names = FALSE) > basis_rank)) {
    stop(sprintf("Each M_s[[s]][l] must not exceed rank(C_g) = %d.", basis_rank))
  }

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
      sdx <- sd(x)
      if (!is.finite(sdx) || sdx < 1e-10) 1 else sdx
    })

    for (s in 1:S) {
      for (i in 1:n_s[s]) {
        Y[[s]][[i]] <- lapply(seq_len(p), function(j) {
          (Y[[s]][[i]][[j]] - mean_mean_across_subjects[j]) /
            sd_mean_across_subjects[j]
        })
        names(Y[[s]][[i]]) <- var_names
      }
    }
  } else {
    mean_mean_across_subjects <- sd_mean_across_subjects <- NULL
  }

  # ---- Variable names ----
  for (s in seq_len(S)) {
    for (i in seq_len(n_s[s])) {
      names(Y[[s]][[i]]) <- var_names
    }
  }

  # ---- Spline dimension ----
  K_total <- K + 2

  # ---- Call core ----
  res <- bayesSYNC_multi_core(
    S = S, n_s = n_s, p = p, d = d,
    L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s,
    K = K, C = C, C_g = C_g,
    Y = Y, Z = Z,
    var_names = var_names,
    mean_mean_across_subjects = mean_mean_across_subjects,
    sd_mean_across_subjects = sd_mean_across_subjects,
    anneal = anneal, list_hyper = list_hyper,
    time_obs = time_obs, bool_var_spec_prob = bool_var_spec_prob,
    n_g = n_g, time_g = time_g,
    tol_abs = tol_abs, tol_rel = tol_rel, maxit = maxit,
    convergence_rule = convergence_rule,
    convergence_rule_requested = convergence_rule_requested,
    practical_control = practical_control,
    scale_trace_control = scale_trace_control,
    elbo_role = elbo_role, lambda_orth = lambda_orth,
    initialization = initialization,
    initialization_control = initialization_control,
    continuation_state = continuation_state,
    fit_seed = seed,
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
    Y, Z, var_names,
    mean_mean_across_subjects,
    sd_mean_across_subjects,
    anneal, list_hyper,
    time_obs, bool_var_spec_prob,
    n_g, time_g,
    tol_abs, tol_rel, maxit,
    convergence_rule, convergence_rule_requested, elbo_role, lambda_orth,
    practical_control,
    initialization, initialization_control, continuation_state, fit_seed,
    n_cpus, verbose,
    scale_trace_control = NULL
) {

  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  eps <- .Machine$double.eps^0.5
  eps_elbo <- 1e-11
  .reset_spd_diagnostics()
  practical_objective_size <- sum(vapply(
    Y, function(study) {
      sum(vapply(study, function(individual) {
        sum(lengths(individual))
      }, numeric(1)))
    }, numeric(1)))
  scale_trace_enabled <- !is.null(scale_trace_control)
  scale_trace_rows <- list()
  scale_trace_initial_state_audit <- NULL
  private_block_schedule <- .private_block_schedule_probe()
  variance_first_loading_schedule <-
    identical(private_block_schedule, "variance_omega_loading")
  # ---- Annealing setup ----
  if (is.null(anneal)) {
    annealing <- FALSE
    c_val <- 1
  } else {
    annealing <- TRUE
    ladder <- get_annealing_ladder_(anneal, verbose)
    c_val <- ladder[1]
  }

  K_total <- K + 2

  # ---- Hyperparameters ----
  sigma_zeta <- list_hyper$sigma_zeta
  mu_beta <- list_hyper$mu_beta
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A
  c_0 <- list_hyper$c_0
  d_0 <- list_hyper$d_0

  inv_Sigma_beta <- .inverse_spd(
    Sigma_beta, context = "linear_prior_covariance")

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

  # ---- Observation counts (per study, per variable for irregular grids) ----
  sum_obs_s <- sapply(1:S, function(s) {
    sum(sapply(time_obs[[s]], length))
  })
  total_obs_sj <- matrix(NA, S, p)
  for (s in 1:S) {
    for (j in 1:p) {
      total_obs_sj[s, j] <- sum(sapply(1:n_s[s], function(i) length(Y[[s]][[i]][[j]])))
    }
  }

  # ===================================================================
  #  INITIALIZATION OF VARIATIONAL PARAMETERS
  # ===================================================================

  # --- 1. Mean function spline coefficients: nu_{mu,sj} ---
  mu_q_recip_sigsq_mu <- matrix(1, nrow = S, ncol = p)
  mu_q_recip_a_mu <- matrix(1, nrow = S, ncol = p)
  mu_q_log_a_mu <- lambda_q_a_mu <- NULL

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
  mu_q_log_a_beta <- lambda_q_a_beta <- NULL

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
  # Hierarchical initialization: decreasing variance per component (1/m)
  mu_q_recip_sigsq_phi <- if (L_f > 0) lapply(1:L_f, function(l) seq(1, M_f[l], by = 1)) else list()
  mu_q_recip_a_phi <- if (L_f > 0) lapply(1:L_f, function(l) rep(1, M_f[l])) else list()
  mu_q_log_a_phi <- lambda_q_a_phi <- NULL

  mu_q_nu_phi <- vector("list", L_f)
  inv_Sigma_q_nu_phi <- vector("list", L_f)
  Sigma_q_nu_phi <- vector("list", L_f)

  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    mu_q_nu_phi[[l]] <- matrix(0, nrow = K_total, ncol = M_l)
    for (m in 1:M_l) {
      mu_q_nu_phi[[l]][, m] <- rnorm(K_total, mean = 0, sd = 1 / sqrt(m))
    }
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
  mu_q_log_a_psi <- lambda_q_a_psi <- NULL

  if (has_specific) {
    mu_q_recip_sigsq_psi <- lapply(1:S, function(s) {
      lapply(seq_len(L_s[[s]]), function(l) seq(1, M_s[[s]][l], by = 1))
    })
    mu_q_recip_a_psi <- lapply(1:S, function(s) {
      lapply(seq_len(L_s[[s]]), function(l) rep(1, M_s[[s]][l]))
    })

    mu_q_nu_psi <- vector("list", S)
    inv_Sigma_q_nu_psi <- vector("list", S)
    Sigma_q_nu_psi <- vector("list", S)

    for (s in 1:S) {
      L_ss <- L_s[[s]]
      mu_q_nu_psi[[s]] <- vector("list", L_ss)
      inv_Sigma_q_nu_psi[[s]] <- vector("list", L_ss)
      Sigma_q_nu_psi[[s]] <- vector("list", L_ss)

      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        mu_q_nu_psi[[s]][[l]] <- matrix(0, nrow = K_total, ncol = M_sl)
        for (m in 1:M_sl) {
          mu_q_nu_psi[[s]][[l]][, m] <- rnorm(K_total, mean = 0, sd = 1 / sqrt(m))
        }
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
      Sigma_q_zeta[[s]][[l]] <-
        lapply(1:n_s[s], function(i) diag(M_l))
    }
  }

  # --- 6. Study-specific factor scores: xi^{(sl)}_{si} ---
  mu_q_xi <- NULL
  Sigma_q_xi <- NULL

  if (has_specific) {
    mu_q_xi <- vector("list", S)
    Sigma_q_xi <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      mu_q_xi[[s]] <- vector("list", L_ss)
      Sigma_q_xi[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        mu_q_xi[[s]][[l]] <- matrix(rnorm(n_s[s] * M_sl, mean = 0, sd = 1),
                                     nrow = n_s[s], ncol = M_sl)
        Sigma_q_xi[[s]][[l]] <-
          lapply(1:n_s[s], function(i) diag(M_sl))
      }
    }
  }

  # --- 7. Shared loadings a_{jl} and PPI gamma_{jl} (Spike-and-Slab) ---
  mu_q_normal_a <- matrix(rnorm(p * max(1, L_f)), nrow = p, ncol = max(1, L_f))
  Sigma_q_normal_a <-
    matrix(1, nrow = p, ncol = max(1, L_f))
  initial_gamma <- if (.driver_diagnostic_initial_dense_gate()) 1 else 0.5
  mu_q_gamma_a <- matrix(initial_gamma, nrow = p, ncol = max(1, L_f))
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

  if (has_specific) {
    mu_q_normal_b <- vector("list", S)
    Sigma_q_normal_b <- vector("list", S)
    mu_q_gamma_b <- vector("list", S)
    mu_q_b_specific <- vector("list", S)
    term_b_specific <- vector("list", S)

    for (s in 1:S) {
      L_ss <- L_s[[s]]
      mu_q_normal_b[[s]] <- matrix(rnorm(p * L_ss), nrow = p, ncol = L_ss)
      Sigma_q_normal_b[[s]] <- matrix(1, nrow = p, ncol = L_ss)
      mu_q_gamma_b[[s]] <- matrix(initial_gamma, nrow = p, ncol = L_ss)
      mu_q_b_specific[[s]] <- mu_q_gamma_b[[s]] * mu_q_normal_b[[s]]
      term_b_specific[[s]] <- (Sigma_q_normal_b[[s]] +
                               mu_q_normal_b[[s]]^2) * mu_q_gamma_b[[s]]
    }
  }

  # --- 8b. Optional data-driven coherent warm start -----------------------
  # The original seeded random state is always created first. This guarantees
  # exact backward compatibility for initialization = "random" and provides
  # explicit per-block random fallbacks when the residual spectrum is rank
  # deficient. Only iteration-zero means are replaced; covariance, variance,
  # PPI, omega, annealing, and every CAVI formula remain unchanged.
  initialization_diagnostics <- list(
    requested_method = initialization,
    used_method = "random",
    complete = TRUE,
    fallback_count = 0L,
    grid_size = NA_integer_,
    residualized_sse = NA_real_,
    remaining_sse = NA_real_,
    explained_fraction = NA_real_,
    factor_diagnostics = data.frame(),
    elapsed_seconds = 0
  )
  if (identical(initialization, "residual_fpca")) {
    warm <- .residual_fpca_initialization(
      Y = Y,
      Z = Z,
      time_obs = time_obs,
      C_g = C_g,
      time_g = time_g,
      S = S,
      n_s = n_s,
      p = p,
      d = d,
      L_f = L_f,
      L_s = L_s,
      M_f = M_f,
      M_s = M_s,
      control = initialization_control
    )
    mu_q_nu_mu <- warm$mu_q_nu_mu
    if (d > 0L) mu_q_nu_beta <- warm$mu_q_nu_beta

    study_end <- cumsum(n_s)
    study_start <- c(1L, head(study_end, -1L) + 1L)
    for (l in seq_len(L_f)) {
      warm_l <- warm$shared[[l]]
      if (!isTRUE(warm_l$success)) next
      mu_q_nu_phi[[l]] <- warm_l$time_coef
      for (s in seq_len(S)) {
        rows_s <- study_start[s]:study_end[s]
        mu_q_zeta[[s]][[l]] <-
          warm_l$scores[rows_s, , drop = FALSE]
      }
      mu_q_normal_a[, l] <-
        warm_l$loading / pmax(mu_q_gamma_a[, l], eps)
    }
    if (L_f > 0L) {
      mu_q_a <- mu_q_gamma_a * mu_q_normal_a
      term_a <-
        (Sigma_q_normal_a + mu_q_normal_a^2) * mu_q_gamma_a
    }

    if (has_specific) {
      for (s in seq_len(S)) {
        for (l in seq_len(L_s[[s]])) {
          warm_sl <- warm$specific[[s]][[l]]
          if (!isTRUE(warm_sl$success)) next
          mu_q_nu_psi[[s]][[l]] <- warm_sl$time_coef
          mu_q_xi[[s]][[l]] <- warm_sl$scores
          mu_q_normal_b[[s]][, l] <-
            warm_sl$loading / pmax(mu_q_gamma_b[[s]][, l], eps)
        }
        mu_q_b_specific[[s]] <-
          mu_q_gamma_b[[s]] * mu_q_normal_b[[s]]
        term_b_specific[[s]] <-
          (Sigma_q_normal_b[[s]] + mu_q_normal_b[[s]]^2) *
            mu_q_gamma_b[[s]]
      }
    }
    initialization_diagnostics <- warm$diagnostics
    if (!isTRUE(initialization_diagnostics$complete)) {
      failed <- initialization_diagnostics$factor_diagnostics
      failed <- failed$block[!failed$success]
      warning(
        "residual_fpca initialization used seeded random fallback for: ",
        paste(failed, collapse = ", "),
        ". See initialization_diagnostics."
      )
    }
  }
  successful_warm_factors <- if (
      identical(initialization, "residual_fpca") &&
      nrow(initialization_diagnostics$factor_diagnostics)) {
    sum(initialization_diagnostics$factor_diagnostics$success)
  } else {
    0L
  }
  initialization_seed_dependent <-
    identical(initialization, "random") ||
    initialization_diagnostics$fallback_count > 0L ||
    (initialization_control$perturb_sd > 0 &&
       p > 1L &&
       successful_warm_factors > 0L)
  initialization_independence_id <- if (
      identical(initialization, "residual_fpca") &&
      !initialization_seed_dependent) {
    "residual_fpca_deterministic"
  } else {
    paste0(
      initialization, "_seed_",
      if (is.null(fit_seed)) "unspecified" else fit_seed
    )
  }

  # --- 8c. Private data-free random-scale calibration -------------------
  # This branch is reachable only through the unexported diagnostic wrapper.
  # It uses C_g and random variational means, never Y, Z, or simulation truth.
  random_scale_calibration_diagnostics <- data.frame()
  random_scale_calibration <-
    .driver_diagnostic_random_scale_calibration()
  if (identical(initialization, "random") &&
      !identical(random_scale_calibration, "none")) {
    calibrated <- .calibrate_random_factor_state(
      C_g = C_g, time_g = time_g,
      mu_q_nu_phi = mu_q_nu_phi,
      mu_q_zeta = mu_q_zeta,
      mu_q_normal_a = mu_q_normal_a,
      Sigma_q_normal_a = Sigma_q_normal_a,
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_nu_psi = mu_q_nu_psi,
      mu_q_xi = mu_q_xi,
      mu_q_normal_b = mu_q_normal_b,
      Sigma_q_normal_b = Sigma_q_normal_b,
      mu_q_gamma_b = mu_q_gamma_b,
      S = S, L_f = L_f, L_s = L_s,
      M_f = M_f, M_s = M_s,
      mode = random_scale_calibration
    )
    mu_q_nu_phi <- calibrated$mu_q_nu_phi
    mu_q_zeta <- calibrated$mu_q_zeta
    mu_q_normal_a <- calibrated$mu_q_normal_a
    mu_q_a <- calibrated$mu_q_a
    term_a <- calibrated$term_a
    mu_q_nu_psi <- calibrated$mu_q_nu_psi
    mu_q_xi <- calibrated$mu_q_xi
    mu_q_normal_b <- calibrated$mu_q_normal_b
    mu_q_b_specific <- calibrated$mu_q_b_specific
    term_b_specific <- calibrated$term_b_specific
    random_scale_calibration_diagnostics <- calibrated$diagnostics
  }

  # --- 9. Variance parameters ---
  mu_q_recip_sigsq_eps <- matrix(1, nrow = S, ncol = p)
  mu_q_recip_a_eps <- matrix(1, nrow = S, ncol = p)
  mu_q_log_a_eps <- lambda_q_a_eps <- NULL

  # --- 9b. Optional validated T=1 continuation state -------------------
  continuation_diagnostics <- list(
    used = FALSE,
    source_fit_seed = NA_integer_,
    source_iterations = NA_integer_,
    source_final_elbo = NA_real_,
    shared_reported_indices = integer(),
    shared_working_indices = integer(),
    inherited_specific_upper_bound = L_s,
    temperature = 1
  )
  if (!is.null(continuation_state)) {
    parameter_names <- .continuation_parameter_names()
    template_parameters <- mget(
      parameter_names, envir = environment(), inherits = FALSE
    )
    response_center_now <- if (is.null(mean_mean_across_subjects)) {
      rep(0, p)
    } else {
      mean_mean_across_subjects
    }
    response_scale_now <- if (is.null(sd_mean_across_subjects)) {
      rep(1, p)
    } else {
      sd_mean_across_subjects
    }
    .validate_continuation_state(
      state = continuation_state,
      template_parameters = template_parameters,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
      bool_var_spec_prob = bool_var_spec_prob,
      response_center = response_center_now,
      response_scale = response_scale_now,
      time_g = time_g
    )
    list2env(continuation_state$parameters, envir = environment())
    continuation_diagnostics <- continuation_state$diagnostics
    initialization_diagnostics <- list(
      requested_method = "continuation_state",
      used_method = "continuation_state",
      complete = TRUE,
      fallback_count = 0L,
      grid_size = length(time_g),
      residualized_sse = NA_real_,
      remaining_sse = NA_real_,
      explained_fraction = NA_real_,
      factor_diagnostics = data.frame(),
      elapsed_seconds = 0,
      continuation = continuation_diagnostics
    )
    initialization_independence_id <- paste0(
      "continuation_from_seed_",
      if (is.null(continuation_diagnostics$source_fit_seed)) {
        "unspecified"
      } else {
        continuation_diagnostics$source_fit_seed
      }
    )
    random_scale_calibration_diagnostics <- data.frame()
  }

  # --- 10. Omega parameters ---
  c_1_omega_a <- rep(NA, L_f)
  d_1_omega_a <- rep(NA, L_f)
  mu_q_log_omega_a <- rep(NA, L_f)
  mu_q_log_1_omega_a <- rep(NA, L_f)

  c_1_omega_b <- NULL
  d_1_omega_b <- NULL
  mu_q_log_omega_b <- NULL
  mu_q_log_1_omega_b <- NULL

  if (has_specific) {
    c_1_omega_b <- vector("list", S)
    d_1_omega_b <- vector("list", S)
    mu_q_log_omega_b <- vector("list", S)
    mu_q_log_1_omega_b <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      c_1_omega_b[[s]] <- rep(NA, L_ss)
      d_1_omega_b[[s]] <- rep(NA, L_ss)
      mu_q_log_omega_b[[s]] <- rep(NA, L_ss)
      mu_q_log_1_omega_b[[s]] <- rep(NA, L_ss)
    }
  }

  # ---- Initial omega update ----
  # Use the same implementation as all subsequent iterations so the public
  # switch changes the hierarchy consistently from the first loading update.
  res_omega <- update_omega(
    mu_q_gamma_a = mu_q_gamma_a,
    mu_q_gamma_b = if (has_specific) mu_q_gamma_b else NULL,
    p = p, S = S, c_0 = c_0, d_0 = d_0,
    c_val = c_val, bool_var_spec_prob = bool_var_spec_prob)
  c_1_omega_a <- res_omega$c_1_omega_a
  d_1_omega_a <- res_omega$d_1_omega_a
  mu_q_log_omega_a <- res_omega$mu_q_log_omega_a
  mu_q_log_1_omega_a <- res_omega$mu_q_log_1_omega_a
  c_1_omega_b <- res_omega$c_1_omega_b
  d_1_omega_b <- res_omega$d_1_omega_b
  mu_q_log_omega_b <- res_omega$mu_q_log_omega_b
  mu_q_log_1_omega_b <- res_omega$mu_q_log_1_omega_b

  # Before the first variance update, the algorithm stores only working
  # precision expectations and has not yet created complete IG q-densities.
  # The trace therefore uses a conditional objective with all unavailable
  # variance-only terms fixed/omitted. Its differences are exact for the
  # coefficient, score, and loading blocks that precede the variance update.
  tempered_objective_ready <- FALSE
  compute_current_tempered_objective <- function(expected_rss_now) {
    if (tempered_objective_ready) {
      log_eps <- mu_q_log_sigsq_eps
      log_mu <- mu_q_log_sigsq_mu
      log_beta <- mu_q_log_sigsq_beta
      log_phi <- mu_q_log_sigsq_phi
      log_psi <- mu_q_log_sigsq_psi
      k_eps <- kappa_q_sigsq_eps
      k_mu <- kappa_q_sigsq_mu
      k_beta <- kappa_q_sigsq_beta
      k_phi <- kappa_q_sigsq_phi
      k_psi <- kappa_q_sigsq_psi
      rate_eps <- lambda_q_sigsq_eps
      rate_mu <- lambda_q_sigsq_mu
      rate_beta <- lambda_q_sigsq_beta
      rate_phi <- lambda_q_sigsq_phi
      rate_psi <- lambda_q_sigsq_psi
      k_aux <- kappa_q_a
      rate_aux_eps <- lambda_q_a_eps
      rate_aux_mu <- lambda_q_a_mu
      rate_aux_beta <- lambda_q_a_beta
      rate_aux_phi <- lambda_q_a_phi
      rate_aux_psi <- lambda_q_a_psi
      include_hc <- TRUE
      scope <- "full"
    } else {
      log_eps <- matrix(0, nrow = S, ncol = p)
      log_mu <- matrix(0, nrow = S, ncol = p)
      log_beta <- if (d > 0) matrix(0, nrow = p, ncol = d) else NULL
      log_phi <- lapply(M_f, function(M_l) rep(0, M_l))
      log_psi <- if (has_specific) {
        lapply(M_s, function(study) {
          lapply(study, function(M_sl) rep(0, M_sl))
        })
      } else {
        NULL
      }
      k_eps <- k_mu <- k_beta <- k_phi <- k_psi <- NULL
      rate_eps <- rate_mu <- rate_beta <- rate_phi <- rate_psi <- NULL
      k_aux <- NULL
      rate_aux_eps <- rate_aux_mu <- rate_aux_beta <- NULL
      rate_aux_phi <- rate_aux_psi <- NULL
      include_hc <- FALSE
      scope <- "prevariance_conditional"
    }

    .compute_variational_objective_multi(
      Y = Y, C = C, list_cp_C = list_cp_C,
      total_obs_sj = total_obs_sj,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta, Sigma_q_nu_beta = Sigma_q_nu_beta,
      Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi, Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
      mu_q_nu_psi = mu_q_nu_psi, Sigma_q_nu_psi = Sigma_q_nu_psi,
      mu_q_a = mu_q_a, mu_q_normal_a = mu_q_normal_a,
      Sigma_q_normal_a = Sigma_q_normal_a,
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_b_specific = mu_q_b_specific,
      mu_q_normal_b = mu_q_normal_b,
      Sigma_q_normal_b = Sigma_q_normal_b,
      mu_q_gamma_b = mu_q_gamma_b,
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
      mu_q_log_sigsq_eps = log_eps,
      mu_q_log_sigsq_mu = log_mu,
      mu_q_log_sigsq_beta = log_beta,
      mu_q_log_sigsq_phi = log_phi,
      mu_q_log_sigsq_psi = log_psi,
      kappa_q_sigsq_eps = k_eps,
      kappa_q_sigsq_mu = k_mu,
      kappa_q_sigsq_beta = k_beta,
      kappa_q_sigsq_phi = k_phi,
      kappa_q_sigsq_psi = k_psi,
      lambda_q_sigsq_eps = rate_eps,
      lambda_q_sigsq_mu = rate_mu,
      lambda_q_sigsq_beta = rate_beta,
      lambda_q_sigsq_phi = rate_phi,
      lambda_q_sigsq_psi = rate_psi,
      kappa_q_a = k_aux,
      mu_q_log_omega_a = mu_q_log_omega_a,
      mu_q_log_1_omega_a = mu_q_log_1_omega_a,
      mu_q_log_omega_b = mu_q_log_omega_b,
      mu_q_log_1_omega_b = mu_q_log_1_omega_b,
      c_1_omega_a = c_1_omega_a, d_1_omega_a = d_1_omega_a,
      c_1_omega_b = c_1_omega_b, d_1_omega_b = d_1_omega_b,
      c_0 = c_0, d_0 = d_0,
      inv_Sigma_beta = inv_Sigma_beta,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
      c_val = c_val,
      term_a = term_a, term_b_specific = term_b_specific,
      expected_rss = expected_rss_now,
      lambda_q_a_eps = rate_aux_eps,
      lambda_q_a_mu = rate_aux_mu,
      lambda_q_a_beta = rate_aux_beta,
      lambda_q_a_phi = rate_aux_phi,
      lambda_q_a_psi = rate_aux_psi,
      A = A, bool_var_spec_prob = bool_var_spec_prob,
      include_half_cauchy = include_hc,
      objective_scope = scope)
  }

  record_scale_trace <- function(stage, iteration, t1_sweep, temperature,
                                 expected_rss, cached_fitted = NULL,
                                 elbo_result = NULL,
                                 elbo_reason = NA_character_,
                                 block = "end_of_sweep") {
    if (!scale_trace_enabled) return(invisible(NULL))
    state <- list(
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s,
      mu_q_nu_mu = mu_q_nu_mu,
      Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta,
      Sigma_q_nu_beta = Sigma_q_nu_beta,
      mu_q_nu_phi = mu_q_nu_phi,
      Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_nu_psi = mu_q_nu_psi,
      Sigma_q_nu_psi = Sigma_q_nu_psi,
      mu_q_zeta = mu_q_zeta,
      Sigma_q_zeta = Sigma_q_zeta,
      mu_q_xi = mu_q_xi,
      Sigma_q_xi = Sigma_q_xi,
      mu_q_normal_a = mu_q_normal_a,
      Sigma_q_normal_a = Sigma_q_normal_a,
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_a = mu_q_a,
      term_a = term_a,
      mu_q_normal_b = mu_q_normal_b,
      Sigma_q_normal_b = Sigma_q_normal_b,
      mu_q_gamma_b = mu_q_gamma_b,
      mu_q_b_specific = mu_q_b_specific,
      term_b_specific = term_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_sigsq_mu = mu_q_recip_sigsq_mu,
      mu_q_recip_sigsq_beta = mu_q_recip_sigsq_beta,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
      mu_q_recip_a_eps = mu_q_recip_a_eps,
      mu_q_recip_a_mu = mu_q_recip_a_mu,
      mu_q_recip_a_beta = mu_q_recip_a_beta,
      mu_q_recip_a_phi = mu_q_recip_a_phi,
      mu_q_recip_a_psi = mu_q_recip_a_psi,
      c_1_omega_a = c_1_omega_a,
      d_1_omega_a = d_1_omega_a,
      c_1_omega_b = c_1_omega_b,
      d_1_omega_b = d_1_omega_b
    )
    if (identical(stage, "initial") &&
        is.null(scale_trace_initial_state_audit)) {
      # Private trace-only raw snapshot. Aggregated RMS summaries are not
      # injective; retaining the original iteration-zero state makes future
      # paired initialization audits exact while leaving ordinary fits inert.
      scale_trace_initial_state_audit <<- state
    }
    scale_trace_rows[[length(scale_trace_rows) + 1L]] <<-
      .make_scale_trace_row(
        stage = stage,
        iteration = iteration,
        t1_sweep = t1_sweep,
        temperature = temperature,
        Y = Y, C = C, C_g = C_g, Z = Z, state = state,
        expected_rss = expected_rss,
        cached_fitted = cached_fitted,
        elbo_result = elbo_result,
        elbo_reason = elbo_reason,
        block = block
      )
    invisible(NULL)
  }

  record_scale_trace_block <- function(block, iteration) {
    trace_scale_block <- scale_trace_enabled && annealing &&
      iteration %in% scale_trace_control$block_annealing_sweeps
    trace_rss_probe <- .block_rss_probe_should_trace(
      block = block, iteration = iteration, annealing = annealing
    )
    if (!trace_scale_block && !trace_rss_probe) {
      return(invisible(NULL))
    }
    block_rss <- compute_rss_cache(
      Y, C, list_cp_C,
      mu_q_nu_mu, Sigma_q_nu_mu,
      mu_q_nu_beta, Sigma_q_nu_beta, Z,
      mu_q_zeta, Sigma_q_zeta,
      mu_q_nu_phi, Sigma_q_nu_phi,
      mu_q_xi, Sigma_q_xi,
      mu_q_nu_psi, Sigma_q_nu_psi,
      mu_q_a, term_a,
      mu_q_b_specific, term_b_specific,
      S, n_s, p, L_f, L_s,
      return_fitted = TRUE
    )
    if (trace_rss_probe) {
      .block_rss_probe_record(
        block = block,
        iteration = iteration,
        temperature = 1 / c_val,
        Y = Y,
        fitted_values = block_rss$fitted_values,
        expected_rss_sum = block_rss$expected_rss_sum,
        total_obs_sj = total_obs_sj,
        C_g = C_g,
        mu_q_nu_phi = mu_q_nu_phi,
        Sigma_q_nu_phi = Sigma_q_nu_phi,
        mu_q_zeta = mu_q_zeta,
        Sigma_q_zeta = Sigma_q_zeta,
        mu_q_normal_a = mu_q_normal_a,
        Sigma_q_normal_a = Sigma_q_normal_a,
        mu_q_gamma_a = mu_q_gamma_a,
        mu_q_a = mu_q_a,
        term_a = term_a,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        L_f = L_f
      )
    }
    if (!trace_scale_block) {
      rm(block_rss)
      return(invisible(NULL))
    }
    block_objective <- compute_current_tempered_objective(
      block_rss$expected_rss)
    record_scale_trace(
      stage = "annealing_block", iteration = iteration,
      t1_sweep = NA_integer_, temperature = 1 / c_val,
      expected_rss = block_rss$expected_rss,
      cached_fitted = block_rss$fitted_values,
      elbo_result = block_objective,
      elbo_reason = if (identical(block_objective$scope, "full")) {
        "complete tempered objective"
      } else {
        paste(
          "conditional objective; variance-only terms unavailable before",
          "the first variance update"
        )
      },
      block = block
    )
    invisible(NULL)
  }

  if (scale_trace_enabled && scale_trace_control$include_initial) {
    initial_rss <- compute_rss_cache(
      Y, C, list_cp_C,
      mu_q_nu_mu, Sigma_q_nu_mu,
      mu_q_nu_beta, Sigma_q_nu_beta, Z,
      mu_q_zeta, Sigma_q_zeta,
      mu_q_nu_phi, Sigma_q_nu_phi,
      mu_q_xi, Sigma_q_xi,
      mu_q_nu_psi, Sigma_q_nu_psi,
      mu_q_a, term_a,
      mu_q_b_specific, term_b_specific,
      S, n_s, p, L_f, L_s,
      return_fitted = TRUE
    )
    initial_objective <- compute_current_tempered_objective(
      initial_rss$expected_rss)
    record_scale_trace(
      stage = "initial", iteration = 0L, t1_sweep = 0L,
      temperature = 1 / c_val,
      expected_rss = initial_rss$expected_rss,
      cached_fitted = initial_rss$fitted_values,
      elbo_result = initial_objective,
      elbo_reason = paste(
        "conditional objective; variance-only terms unavailable before",
        "the first variance update"
      ),
      block = "initial"
    )
    rm(initial_rss)
  }

  # ===================================================================
  #  MAIN CAVI LOOP
  # ===================================================================

  ELBO <- NULL
  ELBO_iter <- NULL
  ELBO_components <- list()
  ELBO_finite <- list()
  ELBO_temperature <- numeric()
  ELBO_diagnostics <- .elbo_history_diagnostics(numeric())
  elbo_result <- NULL
  elbo_t1_jitter_count <- 0L
  jitter_downgrade_warned <- FALSE
  parameter_change <- NULL
  stable_iterations <- 0L
  converged <- FALSE
  practical_streak <- 0L
  practical_snapshot_buffer <- list()
  practical_snapshot_now <- NULL
  practical_diagnostics <- data.frame(
    t1_sweep = integer(),
    eligible = logical(),
    objective_pass = logical(),
    short_output_pass = logical(),
    long_output_pass = logical(),
    output_pass = logical(),
    pass = logical(),
    streak = integer(),
    elbo_abs_rate = double(),
    elbo_per_response_rate = double(),
    elbo_rel_rate = double(),
    fitted_nrmse = double(),
    rss_rel = double(),
    ppi_max_abs = double(),
    ppi_quantile_abs = double(),
    factor_ppi_max_abs = double(),
    long_fitted_nrmse = double(),
    long_rss_rel = double(),
    long_ppi_max_abs = double(),
    long_ppi_quantile_abs = double(),
    long_factor_ppi_max_abs = double(),
    monotone = logical(),
    long_monotone = logical(),
    finite = logical(),
    stringsAsFactors = FALSE
  )
  practical_checkpoints <- list()
  practical_final <- NULL
  objective_converged <- FALSE
  practical_converged <- FALSE
  slow_case <- FALSE
  convergence_status <- "running"
  convergence_reason <- ""

  run_variance_coordinate_update <- function(return_fitted) {
    update_all_variances(
      Y = Y, C = C, list_cp_C = list_cp_C,
      mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta,
      Sigma_q_nu_beta = Sigma_q_nu_beta, Z = Z,
      mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
      mu_q_nu_phi = mu_q_nu_phi,
      Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
      mu_q_nu_psi = mu_q_nu_psi,
      Sigma_q_nu_psi = Sigma_q_nu_psi,
      mu_q_a = mu_q_a, term_a = term_a,
      mu_q_b_specific = mu_q_b_specific,
      term_b_specific = term_b_specific,
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
      total_obs_sj = total_obs_sj, A = A,
      c_val = c_val, n_cpus = n_cpus,
      return_fitted = return_fitted
    )
  }

  install_variance_coordinate_result <- function(
      result, persist_fitted_values) {
    mu_q_recip_sigsq_eps <<- result$mu_q_recip_sigsq_eps
    mu_q_recip_a_eps <<- result$mu_q_recip_a_eps
    mu_q_log_a_eps <<- result$mu_q_log_a_eps
    lambda_q_a_eps <<- result$lambda_q_a_eps
    expected_rss <<- result$expected_rss
    expected_rss_sum <<- result$expected_rss_sum
    fitted_values_this_sweep <<- result$fitted_values
    fitted_values_current <<- if (persist_fitted_values) {
      fitted_values_this_sweep
    } else {
      NULL
    }
    mu_q_recip_sigsq_mu <<- result$mu_q_recip_sigsq_mu
    mu_q_recip_a_mu <<- result$mu_q_recip_a_mu
    mu_q_log_a_mu <<- result$mu_q_log_a_mu
    lambda_q_a_mu <<- result$lambda_q_a_mu
    mu_q_recip_sigsq_phi <<- result$mu_q_recip_sigsq_phi
    mu_q_recip_a_phi <<- result$mu_q_recip_a_phi
    mu_q_log_a_phi <<- result$mu_q_log_a_phi
    lambda_q_a_phi <<- result$lambda_q_a_phi
    if (!is.null(result$mu_q_recip_sigsq_beta)) {
      mu_q_recip_sigsq_beta <<- result$mu_q_recip_sigsq_beta
      mu_q_recip_a_beta <<- result$mu_q_recip_a_beta
      mu_q_log_a_beta <<- result$mu_q_log_a_beta
      lambda_q_a_beta <<- result$lambda_q_a_beta
    }
    if (!is.null(result$mu_q_recip_sigsq_psi)) {
      mu_q_recip_sigsq_psi <<- result$mu_q_recip_sigsq_psi
      mu_q_recip_a_psi <<- result$mu_q_recip_a_psi
      mu_q_log_a_psi <<- result$mu_q_log_a_psi
      lambda_q_a_psi <<- result$lambda_q_a_psi
    }
    kappa_q_a <<- result$kappa_q_a
    mu_q_log_sigsq_eps <<- result$mu_q_log_sigsq_eps
    mu_q_log_sigsq_mu <<- result$mu_q_log_sigsq_mu
    lambda_q_sigsq_eps <<- result$lambda_q_sigsq_eps
    lambda_q_sigsq_mu <<- result$lambda_q_sigsq_mu
    kappa_q_sigsq_eps <<- result$kappa_q_sigsq_eps
    kappa_q_sigsq_mu <<- result$kappa_q_sigsq_mu
    kappa_q_sigsq_phi <<- result$kappa_q_sigsq_phi
    lambda_q_sigsq_phi <<- result$lambda_q_sigsq_phi
    mu_q_log_sigsq_phi <<- result$mu_q_log_sigsq_phi
    kappa_q_sigsq_beta <<- result$kappa_q_sigsq_beta
    lambda_q_sigsq_beta <<- result$lambda_q_sigsq_beta
    mu_q_log_sigsq_beta <<- result$mu_q_log_sigsq_beta
    kappa_q_sigsq_psi <<- result$kappa_q_sigsq_psi
    lambda_q_sigsq_psi <<- result$lambda_q_sigsq_psi
    mu_q_log_sigsq_psi <<- result$mu_q_log_sigsq_psi
    tempered_objective_ready <<- TRUE
    invisible(NULL)
  }

  for (i_iter in 1:maxit) {
    .driver_diagnostic_set_context(
      iteration = i_iter, temperature = 1 / c_val, annealing = annealing)
    .phi_driver_decomposition_set_context(
      iteration = i_iter, temperature = 1 / c_val, annealing = annealing)
    jitter_count_before_iteration <-
      .get_spd_diagnostics()$jitter_count

    old_state <- if (!annealing) list(
      mu_q_nu_mu = mu_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta,
      mu_q_nu_phi = mu_q_nu_phi,
      mu_q_nu_psi = mu_q_nu_psi,
      mu_q_zeta = mu_q_zeta,
      mu_q_xi = mu_q_xi,
      mu_q_a = mu_q_a,
      mu_q_b_specific = mu_q_b_specific,
      slab_a = list(mu_q_normal_a, Sigma_q_normal_a),
      slab_b = list(mu_q_normal_b, Sigma_q_normal_b),
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_gamma_b = mu_q_gamma_b,
      omega_a = list(c_1_omega_a, d_1_omega_a),
      omega_b = list(c_1_omega_b, d_1_omega_b),
      smooth_precision = list(
        mu_q_recip_sigsq_mu, mu_q_recip_sigsq_beta,
        mu_q_recip_sigsq_phi, mu_q_recip_sigsq_psi),
      half_cauchy_aux = list(
        mu_q_recip_a_eps, mu_q_recip_a_mu, mu_q_recip_a_beta,
        mu_q_recip_a_phi, mu_q_recip_a_psi),
      noise_precision = mu_q_recip_sigsq_eps,
      covariance_summary = list(
        .posterior_cov_summary(Sigma_q_nu_mu),
        .posterior_cov_summary(Sigma_q_nu_beta),
        .posterior_cov_summary(Sigma_q_nu_phi),
        .posterior_cov_summary(Sigma_q_nu_psi),
        .posterior_cov_summary(Sigma_q_zeta),
        .posterior_cov_summary(Sigma_q_xi))
    ) else NULL

    record_scale_trace_block("before_updates", i_iter)

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
    inv_Sigma_q_nu_mu <- res_mu$inv_Sigma_q_nu_mu
    record_scale_trace_block("after_mu", i_iter)

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
      record_scale_trace_block("after_beta", i_iter)
    }

    # ------ Private path diagnostic: pre-update score coordinates ------
    # The ordinary path skips this branch.  When enabled, it performs exact
    # score-coordinate updates before the function coordinates, followed by
    # the usual score updates below: score -> function -> score -> loading.
    if (.driver_diagnostic_pre_score_active()) {
      .driver_diagnostic_set_phase("pre_score")
      pre_zeta <- update_zeta(
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
        c_val = c_val, n_cpus = n_cpus
      )
      mu_q_zeta <- pre_zeta$mu_q_zeta
      Sigma_q_zeta <- pre_zeta$Sigma_q_zeta

      if (has_specific) {
        pre_xi <- update_xi(
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
          c_val = c_val, n_cpus = n_cpus
        )
        mu_q_xi <- pre_xi$mu_q_xi
        Sigma_q_xi <- pre_xi$Sigma_q_xi
      }
      .driver_diagnostic_set_phase("main")
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
      term_b_specific = term_b_specific,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      inv_Sigma_beta = inv_Sigma_beta,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, K = K, K_total = K_total,
      c_val = c_val, n_cpus = n_cpus,
      lambda_orth = lambda_orth)
    mu_q_nu_phi <- res_phi$mu_q_nu_phi
    Sigma_q_nu_phi <- res_phi$Sigma_q_nu_phi
    inv_Sigma_q_nu_phi <- res_phi$inv_Sigma_q_nu_phi
    record_scale_trace_block("after_phi", i_iter)
    if (.phi_driver_decomposition_stop_after_phi()) {
      return(.phi_driver_decomposition_early_result())
    }

    # ------ Block 4: Update q(nu_psi) ------
    if (has_specific) {
      res_psi <- update_nu_psi(
        Y = Y, C = C,
        list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
        mu_q_nu_mu = mu_q_nu_mu, Sigma_q_nu_mu = Sigma_q_nu_mu,
        mu_q_nu_beta = mu_q_nu_beta, Z = Z,
        mu_q_zeta = mu_q_zeta, mu_q_nu_phi = mu_q_nu_phi,
        mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
        mu_q_nu_psi = mu_q_nu_psi,
        mu_q_a = mu_q_a,
        term_a = term_a,
        mu_q_b_specific = mu_q_b_specific,
        term_b_specific = term_b_specific,
        mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
        mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
        inv_Sigma_beta = inv_Sigma_beta,
        S = S, n_s = n_s, p = p, d = d,
        L_f = L_f, L_s = L_s, M_s = M_s,
        K = K, K_total = K_total,
        c_val = c_val, n_cpus = n_cpus,
        lambda_orth = lambda_orth)
      mu_q_nu_psi <- res_psi$mu_q_nu_psi
      Sigma_q_nu_psi <- res_psi$Sigma_q_nu_psi
      inv_Sigma_q_nu_psi <- res_psi$inv_Sigma_q_nu_psi
      record_scale_trace_block("after_psi", i_iter)
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
    record_scale_trace_block("after_zeta", i_iter)

    # ------ Block 6: Update q(xi) ------
    if (has_specific) {
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
      record_scale_trace_block("after_xi", i_iter)
    }

    pending_t1_sweep <- if (!annealing) length(ELBO) + 1L else NA_integer_
    trace_annealing_sweep <- scale_trace_enabled &&
      annealing && scale_trace_control$include_annealing
    trace_t1_sweep <- scale_trace_enabled &&
      !annealing &&
      pending_t1_sweep %in% scale_trace_control$t1_sweeps
    persist_fitted_values <- !annealing && (
      identical(convergence_rule, "practical") ||
        pending_t1_sweep %in% practical_control$checkpoints
    )
    need_fitted_values <- trace_annealing_sweep ||
      trace_t1_sweep ||
      persist_fitted_values

    # Private paired-schedule experiment only.  Updating the noise variance
    # here lets the loading coordinate use a precision calibrated to the
    # current mean/time-function/score state.  Ordinary fits never enter this
    # branch and retain the original public algorithm.
    if (variance_first_loading_schedule) {
      res_var <- run_variance_coordinate_update(
        return_fitted = need_fitted_values
      )
      install_variance_coordinate_result(
        res_var, persist_fitted_values = persist_fitted_values
      )
      record_scale_trace_block("after_variances", i_iter)

      res_omega <- update_omega(
        mu_q_gamma_a = mu_q_gamma_a,
        mu_q_gamma_b = if (!is.null(mu_q_gamma_b)) {
          mu_q_gamma_b
        } else {
          NULL
        },
        p = p, S = S,
        c_0 = c_0, d_0 = d_0,
        c_val = c_val, bool_var_spec_prob = bool_var_spec_prob
      )
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
      record_scale_trace_block("after_omega", i_iter)
    }

    # ------ Block 7: Update q(a_loadings) ------
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
    record_scale_trace_block("after_shared_loadings", i_iter)

    # ------ Block 8: Update q(b_loadings) ------
    if (has_specific) {
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
        mu_q_normal_b_specific = mu_q_normal_b,
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
      record_scale_trace_block("after_specific_loadings", i_iter)
    }

    # ------ Block 9: Update all variance parameters ------
    if (!variance_first_loading_schedule) {
      res_var <- run_variance_coordinate_update(
        return_fitted = need_fitted_values
      )
    } else {
      # The variance density above was correctly updated against the pre-load
      # state.  The loading update changes the expected residual only, so
      # refresh that cache without applying a second variance coordinate.
      post_loading_rss <- compute_rss_cache(
        Y, C, list_cp_C,
        mu_q_nu_mu, Sigma_q_nu_mu,
        mu_q_nu_beta, Sigma_q_nu_beta, Z,
        mu_q_zeta, Sigma_q_zeta,
        mu_q_nu_phi, Sigma_q_nu_phi,
        mu_q_xi, Sigma_q_xi,
        mu_q_nu_psi, Sigma_q_nu_psi,
        mu_q_a, term_a,
        mu_q_b_specific, term_b_specific,
        S, n_s, p, L_f, L_s,
        return_fitted = need_fitted_values
      )
      res_var$expected_rss <- post_loading_rss$expected_rss
      res_var$expected_rss_sum <- post_loading_rss$expected_rss_sum
      res_var$fitted_values <- post_loading_rss$fitted_values
      rm(post_loading_rss)
    }

    install_variance_coordinate_result(
      res_var, persist_fitted_values = persist_fitted_values
    )
    if (!variance_first_loading_schedule) {
      record_scale_trace_block("after_variances", i_iter)
    }

    # ------ Block 10: Update q(omega) ------
    if (!variance_first_loading_schedule) {
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
    .annealing_light_trace_record(
      iteration = i_iter,
      temperature = 1 / c_val,
      annealing = annealing,
      C_g = C_g,
      expected_rss_sum = expected_rss_sum,
      total_obs_sj = total_obs_sj,
      mu_q_nu_mu = mu_q_nu_mu,
      Sigma_q_nu_mu = Sigma_q_nu_mu,
      mu_q_nu_phi = mu_q_nu_phi,
      Sigma_q_nu_phi = Sigma_q_nu_phi,
      mu_q_zeta = mu_q_zeta,
      Sigma_q_zeta = Sigma_q_zeta,
      mu_q_normal_a = mu_q_normal_a,
      Sigma_q_normal_a = Sigma_q_normal_a,
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_a = mu_q_a,
      term_a = term_a,
      mu_q_recip_sigsq_eps = mu_q_recip_sigsq_eps,
      mu_q_recip_sigsq_mu = mu_q_recip_sigsq_mu,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      mu_q_recip_a_eps = mu_q_recip_a_eps,
      mu_q_recip_a_mu = mu_q_recip_a_mu,
      mu_q_recip_a_phi = mu_q_recip_a_phi,
      c_1_omega_a = c_1_omega_a,
      d_1_omega_a = d_1_omega_a,
      L_f = L_f
    )
    if (!variance_first_loading_schedule) {
      record_scale_trace_block("after_omega", i_iter)
    }
    if (trace_annealing_sweep) {
      tempered_result <- compute_current_tempered_objective(expected_rss)
      record_scale_trace(
        stage = "annealing", iteration = i_iter,
        t1_sweep = NA_integer_, temperature = 1 / c_val,
        expected_rss = expected_rss,
        cached_fitted = fitted_values_this_sweep,
        elbo_result = tempered_result,
        elbo_reason = "complete tempered objective"
      )
    }
    # ------ ELBO computation (skip during annealing) ------
    if (annealing) {
      ELBO_iter <- NULL
    } else {
    # A fallback jitter changes a Gaussian coordinate update.  It is recorded
    # explicitly and prevents the unperturbed ordinary ELBO from being used as
    # a stopping objective for this fit.  Jitter used only before the T = 1
    # phase is harmless because the starting variational state is arbitrary.
    solver_diagnostics_now <- .get_spd_diagnostics()
    new_t1_jitter <- max(
      0L,
      solver_diagnostics_now$jitter_count -
        jitter_count_before_iteration
    )
    if (new_t1_jitter > 0L) {
      elbo_t1_jitter_count <- elbo_t1_jitter_count + new_t1_jitter
      elbo_role <- "diagnostic_only"
      if (convergence_rule %in% c("elbo", "practical")) {
        convergence_rule <- "parameters"
        stable_iterations <- 0L
        practical_streak <- 0L
        if (!jitter_downgrade_warned) {
          warning(
            "Adaptive Cholesky jitter was required during the T = 1 phase; ",
            "downgrading the requested objective-based convergence rule ",
            "to parameter stopping. ",
            "The unperturbed ordinary-model ELBO is retained as a diagnostic."
          )
          jitter_downgrade_warned <- TRUE
        }
      }
    }

    elbo_result <- compute_elbo_multi(
      Y = Y, C = C, list_cp_C = list_cp_C,
      total_obs_sj = total_obs_sj,
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
      mu_q_recip_a_beta = mu_q_recip_a_beta,
      mu_q_recip_sigsq_phi = mu_q_recip_sigsq_phi,
      mu_q_recip_a_phi = mu_q_recip_a_phi,
      mu_q_recip_sigsq_psi = mu_q_recip_sigsq_psi,
      mu_q_recip_a_psi = mu_q_recip_a_psi,
      mu_q_log_sigsq_eps = mu_q_log_sigsq_eps,
      mu_q_log_sigsq_mu = mu_q_log_sigsq_mu,
      mu_q_log_sigsq_beta = mu_q_log_sigsq_beta,
      mu_q_log_sigsq_phi = mu_q_log_sigsq_phi,
      mu_q_log_sigsq_psi = mu_q_log_sigsq_psi,
      kappa_q_sigsq_eps = kappa_q_sigsq_eps,
      kappa_q_sigsq_mu = kappa_q_sigsq_mu,
      kappa_q_sigsq_beta = kappa_q_sigsq_beta,
      kappa_q_sigsq_phi = kappa_q_sigsq_phi,
      kappa_q_sigsq_psi = kappa_q_sigsq_psi,
      lambda_q_sigsq_eps = lambda_q_sigsq_eps,
      lambda_q_sigsq_mu = lambda_q_sigsq_mu,
      lambda_q_sigsq_beta = lambda_q_sigsq_beta,
      lambda_q_sigsq_phi = lambda_q_sigsq_phi,
      lambda_q_sigsq_psi = lambda_q_sigsq_psi,
      kappa_q_a = kappa_q_a,
      mu_q_log_omega_a = mu_q_log_omega_a,
      mu_q_log_1_omega_a = mu_q_log_1_omega_a,
      mu_q_log_omega_b = mu_q_log_omega_b,
      mu_q_log_1_omega_b = mu_q_log_1_omega_b,
      c_1_omega_a = c_1_omega_a, d_1_omega_a = d_1_omega_a,
      c_1_omega_b = c_1_omega_b, d_1_omega_b = d_1_omega_b,
      c_0 = c_0, d_0 = d_0,
      inv_Sigma_beta = inv_Sigma_beta,
      S = S, n_s = n_s, p = p, d = d,
      L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
      c_val = c_val,
      term_a = term_a, term_b_specific = term_b_specific,
      expected_rss = expected_rss,
      expected_rss_sum = expected_rss_sum,
      lambda_q_a_eps = lambda_q_a_eps,
      lambda_q_a_mu = lambda_q_a_mu,
      lambda_q_a_beta = lambda_q_a_beta,
      lambda_q_a_phi = lambda_q_a_phi,
      lambda_q_a_psi = lambda_q_a_psi,
      A = A, bool_var_spec_prob = bool_var_spec_prob,
      elbo_history = ELBO)

    ELBO_iter <- elbo_result$total
    ELBO <- c(ELBO, ELBO_iter)
    ELBO_components[[length(ELBO_components) + 1L]] <-
      elbo_result$components
    ELBO_finite[[length(ELBO_finite) + 1L]] <- elbo_result$finite
    ELBO_temperature <- c(ELBO_temperature, 1)
    ELBO_diagnostics <- elbo_result$diagnostics

    if (trace_t1_sweep) {
      record_scale_trace(
        stage = "t1", iteration = i_iter,
        t1_sweep = length(ELBO), temperature = 1,
        expected_rss = expected_rss,
        cached_fitted = fitted_values_this_sweep,
        elbo_result = elbo_result
      )
    }

    if (verbose) {
      label <- if (elbo_role == "objective") "ELBO" else "ELBO diagnostic"
      cat(paste0(label, " = ", format(ELBO_iter), "\n"))
    }

    new_state <- list(
      mu_q_nu_mu = mu_q_nu_mu,
      mu_q_nu_beta = mu_q_nu_beta,
      mu_q_nu_phi = mu_q_nu_phi,
      mu_q_nu_psi = mu_q_nu_psi,
      mu_q_zeta = mu_q_zeta,
      mu_q_xi = mu_q_xi,
      mu_q_a = mu_q_a,
      mu_q_b_specific = mu_q_b_specific,
      slab_a = list(mu_q_normal_a, Sigma_q_normal_a),
      slab_b = list(mu_q_normal_b, Sigma_q_normal_b),
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_gamma_b = mu_q_gamma_b,
      omega_a = list(c_1_omega_a, d_1_omega_a),
      omega_b = list(c_1_omega_b, d_1_omega_b),
      smooth_precision = list(
        mu_q_recip_sigsq_mu, mu_q_recip_sigsq_beta,
        mu_q_recip_sigsq_phi, mu_q_recip_sigsq_psi),
      half_cauchy_aux = list(
        mu_q_recip_a_eps, mu_q_recip_a_mu, mu_q_recip_a_beta,
        mu_q_recip_a_phi, mu_q_recip_a_psi),
      noise_precision = mu_q_recip_sigsq_eps,
      covariance_summary = list(
        .posterior_cov_summary(Sigma_q_nu_mu),
        .posterior_cov_summary(Sigma_q_nu_beta),
        .posterior_cov_summary(Sigma_q_nu_phi),
        .posterior_cov_summary(Sigma_q_nu_psi),
        .posterior_cov_summary(Sigma_q_zeta),
        .posterior_cov_summary(Sigma_q_xi))
    )
    change_now <- .parameter_change(old_state, new_state)
    parameter_change <- rbind(parameter_change, change_now)

    if (verbose) {
      cat(paste0("Parameter change: abs=", format(change_now[["abs"]]),
                 ", rel=", format(change_now[["rel"]]), "\n\n"))
    }

    t1_sweeps_now <- length(ELBO)
    practical_snapshot_now <- NULL
    checkpoint_now <- t1_sweeps_now %in% practical_control$checkpoints
    if (identical(convergence_rule, "practical") || checkpoint_now) {
      practical_snapshot_now <- .practical_snapshot(
        fitted_values = fitted_values_current,
        expected_rss = expected_rss,
        mu_q_gamma_a = mu_q_gamma_a,
        mu_q_gamma_b = mu_q_gamma_b
      )
      if (checkpoint_now) {
        practical_checkpoints[[as.character(t1_sweeps_now)]] <- list(
          t1_sweep = t1_sweeps_now,
          total_iteration = i_iter,
          ELBO = ELBO_iter,
          data_likelihood =
            unname(elbo_result$components[["data_likelihood"]]),
          snapshot = practical_snapshot_now,
          parameter_change = change_now,
          t1_jitter_count = elbo_t1_jitter_count,
          solver_jitter_count =
            .get_spd_diagnostics()$jitter_count
        )
      }
    }

    practical_window_now <- NULL
    if (convergence_rule == "practical") {
      practical_snapshot_buffer[[as.character(t1_sweeps_now)]] <-
        practical_snapshot_now
      reference_sweep <- t1_sweeps_now - practical_control$window
      reference_snapshot <- practical_snapshot_buffer[[
        as.character(reference_sweep)]]
      long_reference_sweep <-
        t1_sweeps_now - practical_control$long_window
      long_reference_snapshot <- practical_snapshot_buffer[[
        as.character(long_reference_sweep)]]
      practical_window_now <- .practical_window_diagnostic(
        sweep = t1_sweeps_now,
        ELBO = ELBO,
        ELBO_diagnostics = ELBO_diagnostics,
        current_snapshot = practical_snapshot_now,
        reference_snapshot = reference_snapshot,
        long_reference_snapshot = long_reference_snapshot,
        control = practical_control,
        objective_size = practical_objective_size
      )
      practical_streak <- if (isTRUE(practical_window_now$pass)) {
        practical_streak + 1L
      } else {
        0L
      }
      objective_converged <- objective_converged ||
        isTRUE(practical_window_now$objective_pass)
      practical_converged <- practical_streak >=
        practical_control$consecutive
      practical_final <- c(
        practical_window_now,
        list(streak = practical_streak)
      )
      practical_diagnostics <- rbind(
        practical_diagnostics,
        data.frame(
          t1_sweep = t1_sweeps_now,
          eligible = practical_window_now$eligible,
          objective_pass = practical_window_now$objective_pass,
          short_output_pass =
            practical_window_now$short_output_pass,
          long_output_pass =
            practical_window_now$long_output_pass,
          output_pass = practical_window_now$output_pass,
          pass = practical_window_now$pass,
          streak = practical_streak,
          elbo_abs_rate = practical_window_now$elbo_abs_rate,
          elbo_per_response_rate =
            practical_window_now$elbo_per_response_rate,
          elbo_rel_rate = practical_window_now$elbo_rel_rate,
          fitted_nrmse = practical_window_now$fitted_nrmse,
          rss_rel = practical_window_now$rss_rel,
          ppi_max_abs = practical_window_now$ppi_max_abs,
          ppi_quantile_abs =
            practical_window_now$ppi_quantile_abs,
          factor_ppi_max_abs =
            practical_window_now$factor_ppi_max_abs,
          long_fitted_nrmse =
            practical_window_now$long_fitted_nrmse,
          long_rss_rel = practical_window_now$long_rss_rel,
          long_ppi_max_abs =
            practical_window_now$long_ppi_max_abs,
          long_ppi_quantile_abs =
            practical_window_now$long_ppi_quantile_abs,
          long_factor_ppi_max_abs =
            practical_window_now$long_factor_ppi_max_abs,
          monotone = practical_window_now$monotone,
          long_monotone = practical_window_now$long_monotone,
          finite = practical_window_now$finite,
          stringsAsFactors = FALSE
        )
      )
      keep_sweeps <- suppressWarnings(as.integer(
        names(practical_snapshot_buffer)))
      practical_snapshot_buffer <- practical_snapshot_buffer[
        keep_sweeps >= max(1L, long_reference_sweep)]
    }

    # The ordinary ELBO can stop the model-derived CAVI only when
    # lambda_orth = 0. Positive soft regularisation is downgraded to the
    # parameter rule at the public interface.
    converged_now <- FALSE
    if (convergence_rule == "parameters") {
      small_change <- change_now[["abs"]] < tol_abs ||
        change_now[["rel"]] < tol_rel
      stable_iterations <- if (small_change) stable_iterations + 1L else 0L
      converged_now <- stable_iterations >= 2L
    } else if (convergence_rule == "elbo" && length(ELBO) >= 2) {
      l_ELBO <- length(ELBO)
      ELBO_diff <- ELBO[l_ELBO] - ELBO[l_ELBO - 1]
      rel_converged <- (abs(ELBO_diff) /
        (1 + abs(ELBO[l_ELBO - 1])) < tol_rel)
      abs_converged <- (abs(ELBO_diff) < tol_abs)
      allowed_decrease <- utils::tail(
        ELBO_diagnostics$scaled_tolerance, 1L)
      monotone_step <- is.finite(ELBO_diff) &&
        ELBO_diff >= -allowed_decrease
      converged_now <- monotone_step && (rel_converged || abs_converged)
    } else if (convergence_rule == "practical") {
      converged_now <- practical_converged
    }
    # A private dense-gate diagnostic must reach the release sweep before
    # stopping.  This branch is unreachable in ordinary fits.
    if (.driver_diagnostic_dense_gate_active()) converged_now <- FALSE

    if (converged_now) {
      converged <- TRUE
      convergence_status <- paste0("converged_", convergence_rule)
      convergence_reason <- if (convergence_rule == "practical") {
        "finite_window_elbo_and_scientific_outputs_stable"
      } else if (convergence_rule == "elbo") {
        "ordinary_elbo_increment_below_tolerance"
      } else {
        "strict_parameter_change_below_tolerance"
      }
      if (verbose) {
        cat(paste0("Convergence obtained after ", i_iter,
                   " iterations using the ", convergence_rule, " rule.\n"))
      }
      break
    } else if (convergence_rule_requested == "practical" &&
               t1_sweeps_now >= practical_control$max_t1) {
      slow_case <- TRUE
      convergence_status <- "slow_case"
      convergence_reason <- if (convergence_rule == "practical") {
        "practical_rule_not_satisfied_before_max_t1"
      } else {
        "objective_rule_downgraded_and_parameter_rule_not_satisfied_before_max_t1"
      }
      warning(
        "The requested practical fit did not produce an eligible result within ",
        practical_control$max_t1,
        " T = 1 sweeps; marking this fit as a slow case."
      )
      break
    } else if (i_iter == maxit) {
      warning("Max iterations reached before convergence.")
    }
    } # end if(!annealing)
    
    # ---- Temperature annealing step ----
    if (annealing) {
      if (verbose) cat(paste0("Temperature = ", format(1 / c_val, digits = 4), "\n\n"))
      c_val <- ifelse(i_iter < length(ladder), ladder[i_iter + 1], 1)
      if (isTRUE(all.equal(c_val, 1))) {
        annealing <- FALSE
        if (verbose) cat("** Exiting annealing mode. **\n\n")
      }
    }
  } # end main loop

  annealing_completed <- !annealing && isTRUE(all.equal(c_val, 1))
  final_temperature <- 1 / c_val
  t1_sweeps <- length(ELBO)
  annealing_sweeps <- i_iter - t1_sweeps
  if (identical(convergence_status, "running")) {
    if (!annealing_completed) {
      convergence_status <- "annealing_incomplete"
      convergence_reason <- "fit_ended_before_T1"
    } else if (!converged) {
      convergence_status <- "maxit_reached"
      convergence_reason <-
        "total_iteration_limit_reached_before_convergence"
    }
  }
  if (!annealing_completed) {
    warning("The fit ended before reaching T = 1; post-processing results should not be interpreted as a final posterior fit.")
  }

  # ===================================================================
  #  POST-PROCESSING: Orthonormalisation
  # ===================================================================
  
  res_orth <- orthonormalise_multi(C_g = C_g, time_g = time_g,
    mu_q_nu_mu = mu_q_nu_mu,
    mu_q_nu_beta = mu_q_nu_beta,
    mu_q_nu_phi = mu_q_nu_phi, mu_q_nu_psi = mu_q_nu_psi,
    mu_q_zeta = mu_q_zeta, Sigma_q_zeta = Sigma_q_zeta,
    mu_q_xi = mu_q_xi, Sigma_q_xi = Sigma_q_xi,
    Sigma_q_nu_phi = Sigma_q_nu_phi,
    Sigma_q_nu_psi = Sigma_q_nu_psi,
    mu_q_a = mu_q_a, mu_q_b_specific = mu_q_b_specific,
    mu_q_gamma_a = mu_q_gamma_a, mu_q_gamma_b = mu_q_gamma_b,
    S = S, n_s = n_s, p = p, d = d,
    L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s,
    response_center = if (is.null(mean_mean_across_subjects)) rep(0, p) else mean_mean_across_subjects,
    response_scale = if (is.null(sd_mean_across_subjects)) rep(1, p) else sd_mean_across_subjects)

  # Preserve working-scale variational objects for reproducibility, and add
  # explicit original-scale posterior summaries for scientific interpretation.
  response_center <- if (is.null(mean_mean_across_subjects)) rep(0, p) else mean_mean_across_subjects
  response_scale <- if (is.null(sd_mean_across_subjects)) rep(1, p) else sd_mean_across_subjects
  mu_q_nu_mu_original <- lapply(seq_len(S), function(s) {
    lapply(seq_len(p), function(j) {
      ans <- response_scale[j] * mu_q_nu_mu[[s]][[j]]
      ans[1] <- ans[1] + response_center[j]
      ans
    })
  })
  Sigma_q_nu_mu_original <- lapply(seq_len(S), function(s) {
    lapply(seq_len(p), function(j) response_scale[j]^2 * Sigma_q_nu_mu[[s]][[j]])
  })
  mu_q_nu_beta_original <- Sigma_q_nu_beta_original <- NULL
  if (d > 0L) {
    mu_q_nu_beta_original <- lapply(seq_len(p), function(j) {
      lapply(seq_len(d), function(r) response_scale[j] * mu_q_nu_beta[[j]][[r]])
    })
    Sigma_q_nu_beta_original <- lapply(seq_len(p), function(j) {
      lapply(seq_len(d), function(r) response_scale[j]^2 * Sigma_q_nu_beta[[j]][[r]])
    })
  }
  mu_q_a_original <- sweep(mu_q_a, 1, response_scale, "*")
  term_a_original <- sweep(term_a, 1, response_scale^2, "*")
  mu_q_b_specific_original <- if (has_specific) {
    lapply(mu_q_b_specific, function(x) sweep(x, 1, response_scale, "*"))
  } else NULL
  term_b_specific_original <- if (has_specific) {
    lapply(term_b_specific, function(x) sweep(x, 1, response_scale^2, "*"))
  } else NULL
  mu_q_recip_sigsq_eps_original <- sweep(
    mu_q_recip_sigsq_eps, 2, response_scale^2, "/")
  lambda_q_sigsq_eps_original <- sweep(
    lambda_q_sigsq_eps, 2, response_scale^2, "*")
  sigsq_eps_hat_working <- ifelse(
    kappa_q_sigsq_eps > 1,
    lambda_q_sigsq_eps / (kappa_q_sigsq_eps - 1),
    1 / mu_q_recip_sigsq_eps)
  sigsq_eps_hat <- sweep(sigsq_eps_hat_working, 2, response_scale^2, "*")
  linear_solver_diagnostics <- .get_spd_diagnostics()
  elbo_objective_valid <-
    lambda_orth == 0 && elbo_t1_jitter_count == 0L
  elbo_invalid_reasons <- character()
  if (lambda_orth > 0) {
    elbo_invalid_reasons <- c(
      elbo_invalid_reasons,
      "positive_soft_orthogonality_penalty"
    )
  }
  if (elbo_t1_jitter_count > 0L) {
    elbo_invalid_reasons <- c(
      elbo_invalid_reasons,
      "adaptive_cholesky_jitter_at_T1"
    )
  }
  if (!elbo_objective_valid) {
    elbo_role <- "diagnostic_only"
  }
  ELBO_decrease_count <- ELBO_diagnostics$decrease_count
  ELBO_min_difference <- ELBO_diagnostics$minimum_difference
  fit_control <- list(
    anneal = anneal,
    tol_abs = tol_abs,
    tol_rel = tol_rel,
    maxit = maxit,
    convergence_rule_requested = convergence_rule_requested,
    practical_control = practical_control,
    lambda_orth = lambda_orth,
    initialization = initialization,
    initialization_control = initialization_control,
    continuation_used = isTRUE(continuation_diagnostics$used)
  )
  # Retain only the compact, truth-free outputs required to compare
  # independent initialisations.  This is much smaller than retaining every
  # fitted variational object in the multi-start wrapper.
  stability_snapshot <- if (!is.null(practical_snapshot_now)) {
    practical_snapshot_now
  } else {
    .practical_snapshot(
      fitted_values = fitted_values_current,
      expected_rss = expected_rss,
      mu_q_gamma_a = mu_q_gamma_a,
      mu_q_gamma_b = mu_q_gamma_b
    )
  }

  # ===================================================================
  #  RESULTS ASSEMBLY
  # ===================================================================

  res <- create_named_list(
    S, n_s, p, d,
    L_f, L_s, M_f, M_s,
    K, K_total,
    list_hyper, c_0, d_0, n_cpus_used = n_cpus,
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
    mu_q_log_sigsq_eps, mu_q_log_sigsq_mu,
    mu_q_log_sigsq_beta, mu_q_log_sigsq_phi, mu_q_log_sigsq_psi,
    kappa_q_sigsq_eps, lambda_q_sigsq_eps,
    kappa_q_sigsq_mu, lambda_q_sigsq_mu,
    kappa_q_sigsq_beta, lambda_q_sigsq_beta,
    kappa_q_sigsq_phi, lambda_q_sigsq_phi,
    kappa_q_sigsq_psi, lambda_q_sigsq_psi,
    kappa_q_a,
    mu_q_log_a_eps, lambda_q_a_eps,
    mu_q_log_a_mu, lambda_q_a_mu,
    mu_q_log_a_beta, lambda_q_a_beta,
    mu_q_log_a_phi, lambda_q_a_phi,
    mu_q_log_a_psi, lambda_q_a_psi,
    expected_rss, expected_rss_sum,
    c_1_omega_a, d_1_omega_a, mu_q_log_omega_a, mu_q_log_1_omega_a,
    c_1_omega_b, d_1_omega_b, mu_q_log_omega_b, mu_q_log_1_omega_b,
    ELBO, ELBO_iter, ELBO_components, ELBO_finite, ELBO_temperature,
    ELBO_diagnostics, ELBO_decrease_count, ELBO_min_difference,
    elbo_result, elbo_role, elbo_objective_valid,
    elbo_invalid_reasons, elbo_t1_jitter_count,
    parameter_change, convergence_rule, convergence_rule_requested,
    practical_control, practical_diagnostics, practical_final,
    practical_objective_size,
    practical_checkpoints, objective_converged, practical_converged,
    stability_snapshot,
    slow_case, convergence_status, convergence_reason,
    lambda_orth, bool_var_spec_prob, i_iter,
    t1_sweeps, annealing_sweeps, fit_seed, fit_control,
    converged, annealing_completed, final_temperature,
    linear_solver_diagnostics,
    mean_mean_across_subjects, sd_mean_across_subjects,
    response_center, response_scale, var_names,
    time_g, C_g, n_g,
    list_cp_C, list_cp_C_Y, list_cp_Y, sum_list_cp_C,
    inv_Sigma_q_nu_mu, inv_Sigma_q_nu_phi, inv_Sigma_q_nu_psi,

    # Explicit original-scale posterior summaries.  Unsuffixed variational
    # objects above remain on the working scale for backward compatibility.
    mu_q_nu_mu_original, Sigma_q_nu_mu_original,
    mu_q_nu_beta_original, Sigma_q_nu_beta_original,
    mu_q_a_original, term_a_original,
    mu_q_b_specific_original, term_b_specific_original,
    mu_q_recip_sigsq_eps_original, lambda_q_sigsq_eps_original,
    sigsq_eps_hat_working, sigsq_eps_hat,
    
    # Orthonormalised outputs
    list_Phi_hat = res_orth$list_Phi_hat,
    list_Zeta_hat = res_orth$list_Zeta_hat,
    list_Cov_zeta_hat = res_orth$list_Cov_zeta_hat,
    list_eigenvalues = res_orth$list_eigenvalues,
    list_cumulated_pve = res_orth$list_cumulated_pve,
    list_effective_M = res_orth$list_effective_M,
    list_full_posterior_integrated_variance =
      res_orth$list_full_posterior_integrated_variance,
    list_omitted_uncertainty_variance =
      res_orth$list_omitted_uncertainty_variance,
    list_rank_cap = res_orth$list_rank_cap,
    list_Phi_hat_spec = res_orth$list_Phi_hat_spec,
    list_Zeta_hat_spec = res_orth$list_Zeta_hat_spec,
    list_Cov_zeta_hat_spec = res_orth$list_Cov_zeta_hat_spec,
    list_eigenvalues_spec = res_orth$list_eigenvalues_spec,
    list_pve_spec = res_orth$list_pve_spec,
    list_effective_M_spec = res_orth$list_effective_M_spec,
    list_full_posterior_integrated_variance_spec =
      res_orth$list_full_posterior_integrated_variance_spec,
    list_omitted_uncertainty_variance_spec =
      res_orth$list_omitted_uncertainty_variance_spec,
    list_rank_cap_spec = res_orth$list_rank_cap_spec,
    list_mu_hat = res_orth$list_mu_hat,
    list_mu_hat_legacy = res_orth$list_mu_hat_legacy,
    list_beta_hat = res_orth$list_beta_hat,
    mu_q_a_hat = res_orth$mu_q_a_hat,
    mu_q_b_specific_hat = res_orth$mu_q_b_specific_hat,
    mu_q_gamma_a_hat = res_orth$mu_q_gamma_a_hat,
    mu_q_gamma_b_hat = res_orth$mu_q_gamma_b_hat,
    factor_scale_shared = res_orth$factor_scale_shared,
    factor_scale_specific = res_orth$factor_scale_specific,
    factor_sign_shared = res_orth$factor_sign_shared,
    factor_sign_specific = res_orth$factor_sign_specific,
    factor_order_shared = res_orth$factor_order_shared,
    factor_order_specific = res_orth$factor_order_specific,
    factor_ppi_shared = res_orth$factor_ppi_shared,
    factor_ppi_specific = res_orth$factor_ppi_specific,

    # Truth-free audit of iteration-zero construction.
    initialization, initialization_control,
    initialization_diagnostics, initialization_independence_id,
    continuation_diagnostics
  )
  res$L_s_by_study <- L_s
  res$L_s_max <- .max_L_s(L_s)
  res$L_s <- .compact_L_s(L_s)

  if (scale_trace_enabled) {
    scale_trace <- if (length(scale_trace_rows)) {
      result <- do.call(rbind, scale_trace_rows)
      rownames(result) <- NULL
      result
    } else {
      data.frame()
    }
    attr(scale_trace, "trace_control") <- scale_trace_control
    attr(scale_trace, "initial_state_audit") <-
      scale_trace_initial_state_audit
    res$scale_trace <- scale_trace
  }

  res
}


# ---- Complete variational-objective building blocks -----------------------

.elbo_logdet_positive <- function(x) {
  if (!is.matrix(x) || nrow(x) != ncol(x) || any(!is.finite(x))) {
    return(NA_real_)
  }
  det_x <- determinant((x + t(x)) / 2, logarithm = TRUE)
  if (!is.finite(det_x$modulus) || det_x$sign <= 0) return(NA_real_)
  as.numeric(det_x$modulus)
}

.elbo_osullivan_block <- function(mean, covariance,
                                  expected_recip_sigsq,
                                  expected_log_sigsq,
                                  inv_Sigma_0, logdet_Sigma_0,
                                  K, temperature = 1) {
  D <- K + 2L
  if (length(temperature) != 1L || !is.finite(temperature) ||
      temperature <= 0) {
    return(NA_real_)
  }
  if (length(mean) != D ||
      !is.matrix(covariance) ||
      !isTRUE(all(dim(covariance) == c(D, D))) ||
      !is.matrix(inv_Sigma_0) ||
      !isTRUE(all(dim(inv_Sigma_0) == c(2L, 2L)))) {
    return(NA_real_)
  }
  logdet_q <- .elbo_logdet_positive(covariance)
  alpha <- mean[1:2]
  covariance_alpha <- covariance[1:2, 1:2, drop = FALSE]
  u_index <- 2L + seq_len(K)
  u <- mean[u_index]
  covariance_u <- covariance[u_index, u_index, drop = FALSE]

  alpha_second <- as.numeric(crossprod(alpha, inv_Sigma_0 %*% alpha)) +
    tr(inv_Sigma_0 %*% covariance_alpha)
  u_second <- as.numeric(crossprod(u)) + tr(covariance_u)

  0.5 * temperature * logdet_q - 0.5 * logdet_Sigma_0 +
    temperature * D / 2 +
    (temperature - 1) * D / 2 * log(2 * pi) -
    0.5 * alpha_second -
    K / 2 * expected_log_sigsq -
    0.5 * expected_recip_sigsq * u_second
}

.elbo_standard_normal_score <- function(mean, covariance,
                                        temperature = 1) {
  M <- length(mean)
  if (length(temperature) != 1L || !is.finite(temperature) ||
      temperature <= 0) {
    return(NA_real_)
  }
  if (!is.matrix(covariance) ||
      !isTRUE(all(dim(covariance) == c(M, M)))) {
    return(NA_real_)
  }
  0.5 * temperature * .elbo_logdet_positive(covariance) +
    temperature * M / 2 +
    (temperature - 1) * M / 2 * log(2 * pi) -
    0.5 * (as.numeric(crossprod(mean)) + tr(covariance))
}

.elbo_ig_entropy <- function(shape, rate) {
  if (length(shape) != 1L || length(rate) != 1L ||
      !is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
    return(NA_real_)
  }
  log(rate) + lgamma(shape) - (shape + 1) * digamma(shape) + shape
}

.elbo_half_cauchy_pair <- function(shape_sigsq, rate_sigsq,
                                   shape_aux, rate_aux, A,
                                   temperature = 1) {
  vals <- c(shape_sigsq, rate_sigsq, shape_aux, rate_aux, A,
            temperature)
  if (any(!is.finite(vals)) || any(vals <= 0)) return(NA_real_)

  elog_sigsq <- log(rate_sigsq) - digamma(shape_sigsq)
  erecip_sigsq <- shape_sigsq / rate_sigsq
  elog_aux <- log(rate_aux) - digamma(shape_aux)
  erecip_aux <- shape_aux / rate_aux

  expected_log_prior_sigsq_given_aux <-
    -0.5 * elog_aux - lgamma(0.5) -
    1.5 * elog_sigsq - erecip_aux * erecip_sigsq
  expected_log_prior_aux <-
    0.5 * log(1 / A^2) - lgamma(0.5) -
    1.5 * elog_aux - erecip_aux / A^2

  expected_log_prior_sigsq_given_aux + expected_log_prior_aux +
    temperature * .elbo_ig_entropy(shape_sigsq, rate_sigsq) +
    temperature * .elbo_ig_entropy(shape_aux, rate_aux)
}

.elbo_xlogx <- function(x) {
  ifelse(x == 0, 0, x * log(x))
}

.elbo_spike_slab_scalar <- function(ppi, slab_mean, slab_variance,
                                    expected_log_omega,
                                    expected_log_one_minus_omega,
                                    temperature = 1) {
  vals <- c(ppi, slab_mean, slab_variance,
            expected_log_omega, expected_log_one_minus_omega,
            temperature)
  if (any(!is.finite(vals)) || ppi < 0 || ppi > 1 ||
      slab_variance <= 0) {
    return(NA_real_)
  }
  0.5 * ppi * (
    temperature * (log(slab_variance) + 1) +
      (temperature - 1) * log(2 * pi) -
      slab_mean^2 - slab_variance
  ) +
    ppi * expected_log_omega +
    (1 - ppi) * expected_log_one_minus_omega -
    temperature * (.elbo_xlogx(ppi) + .elbo_xlogx(1 - ppi))
}

.elbo_beta_prior_entropy <- function(shape1_q, shape2_q,
                                     shape1_prior, shape2_prior,
                                     temperature = 1) {
  vals <- c(shape1_q, shape2_q, shape1_prior, shape2_prior,
            temperature)
  if (any(!is.finite(vals)) || any(vals <= 0)) return(NA_real_)
  elog_omega <- digamma(shape1_q) - digamma(shape1_q + shape2_q)
  elog_one_minus <- digamma(shape2_q) - digamma(shape1_q + shape2_q)
  (shape1_prior - 1 - temperature * (shape1_q - 1)) * elog_omega +
    (shape2_prior - 1 - temperature * (shape2_q - 1)) *
      elog_one_minus +
    temperature * lbeta(shape1_q, shape2_q) -
    lbeta(shape1_prior, shape2_prior)
}

.elbo_history_diagnostics <- function(history,
                                      tolerance = sqrt(.Machine$double.eps)) {
  history <- as.numeric(history)
  differences <- diff(history)
  if (length(differences)) {
    previous <- utils::head(history, -1L)
    current <- utils::tail(history, -1L)
    scaled_tolerance <- tolerance *
      (1 + pmax(abs(previous), abs(current)))
    decreases <- is.finite(differences) &
      differences < -scaled_tolerance
    finite_differences <- differences[is.finite(differences)]
    minimum_difference <- if (length(finite_differences)) {
      min(finite_differences)
    } else {
      NA_real_
    }
    last_difference <- utils::tail(differences, 1L)
  } else {
    scaled_tolerance <- numeric()
    decreases <- logical()
    minimum_difference <- NA_real_
    last_difference <- NA_real_
  }

  list(
    history_length = length(history),
    differences = differences,
    scaled_tolerance = scaled_tolerance,
    decrease_count = sum(decreases),
    minimum_difference = minimum_difference,
    last_difference = last_difference,
    monotone_within_tolerance = !any(decreases),
    tolerance = tolerance
  )
}

.compute_variational_objective_multi <- function(Y, C, list_cp_C,
                                total_obs_sj = NULL,
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
                                mu_q_recip_sigsq_beta, mu_q_recip_a_beta,
                                mu_q_recip_sigsq_phi, mu_q_recip_a_phi,
                                mu_q_recip_sigsq_psi, mu_q_recip_a_psi,
                                mu_q_log_sigsq_eps, mu_q_log_sigsq_mu,
                                mu_q_log_sigsq_beta, mu_q_log_sigsq_phi,
                                mu_q_log_sigsq_psi,
                                kappa_q_sigsq_eps, kappa_q_sigsq_mu,
                                kappa_q_sigsq_beta, kappa_q_sigsq_phi,
                                kappa_q_sigsq_psi,
                                lambda_q_sigsq_eps, lambda_q_sigsq_mu,
                                lambda_q_sigsq_beta, lambda_q_sigsq_phi,
                                lambda_q_sigsq_psi,
                                kappa_q_a,
                                mu_q_log_omega_a, mu_q_log_1_omega_a,
                                mu_q_log_omega_b, mu_q_log_1_omega_b,
                                c_1_omega_a, d_1_omega_a,
                                c_1_omega_b, d_1_omega_b,
                                c_0, d_0,
                                inv_Sigma_beta,
                                S, n_s, p, d,
                                L_f, L_s, M_f, M_s, K,
                                c_val = 1,
                                term_a = NULL,
                                term_b_specific = NULL,
                                expected_rss = NULL,
                                expected_rss_sum = NULL,
                                lambda_q_a_eps = NULL,
                                lambda_q_a_mu = NULL,
                                lambda_q_a_beta = NULL,
                                lambda_q_a_phi = NULL,
                                lambda_q_a_psi = NULL,
                                 A = 1e5,
                                 bool_var_spec_prob = FALSE,
                                 elbo_history = numeric(),
                                 monotonicity_tolerance =
                                   sqrt(.Machine$double.eps),
                                 include_half_cauchy = TRUE,
                                 objective_scope = "full") {
  L_s <- .normalize_L_s(L_s, S)
  has_specific <- .has_specific(L_s)
  if (length(c_val) != 1L || !is.finite(c_val) || c_val <= 0) {
    stop("c_val must be one finite positive inverse temperature.")
  }
  temperature <- 1 / c_val
  if (length(A) != 1L || !is.finite(A) || A <= 0) {
    stop("A must be one finite positive Half-Cauchy scale.")
  }
  if (length(include_half_cauchy) != 1L ||
      is.na(include_half_cauchy)) {
    stop("include_half_cauchy must be TRUE or FALSE.")
  }
  if (length(objective_scope) != 1L || is.na(objective_scope) ||
      !nzchar(objective_scope)) {
    stop("objective_scope must be one non-empty string.")
  }
  if (length(monotonicity_tolerance) != 1L ||
      !is.finite(monotonicity_tolerance) ||
      monotonicity_tolerance < 0) {
    stop("monotonicity_tolerance must be finite and non-negative.")
  }

  if (is.null(term_a)) {
    term_a <- (Sigma_q_normal_a + mu_q_normal_a^2) * mu_q_gamma_a
  }
  if (has_specific && is.null(term_b_specific)) {
    term_b_specific <- lapply(seq_len(S), function(s) {
      (Sigma_q_normal_b[[s]] + mu_q_normal_b[[s]]^2) *
        mu_q_gamma_b[[s]]
    })
  }

  if (is.null(expected_rss_sum)) {
    if (!is.null(expected_rss)) {
      expected_rss_sum <- do.call(rbind, lapply(expected_rss, colSums))
      rss_source <- "supplied_subject_cache"
    } else {
      rss_result <- compute_rss_cache(
        Y, C, list_cp_C,
        mu_q_nu_mu, Sigma_q_nu_mu,
        mu_q_nu_beta, Sigma_q_nu_beta, Z,
        mu_q_zeta, Sigma_q_zeta,
        mu_q_nu_phi, Sigma_q_nu_phi,
        mu_q_xi, Sigma_q_xi,
        mu_q_nu_psi, Sigma_q_nu_psi,
        mu_q_a, term_a,
        mu_q_b_specific, term_b_specific,
        S, n_s, p, L_f, L_s)
      expected_rss <- rss_result$expected_rss
      expected_rss_sum <- rss_result$expected_rss_sum
      rss_source <- "direct_recalculation"
    }
  } else {
    rss_source <- "variance_update_cache"
  }
  if (!isTRUE(all(dim(expected_rss_sum) == c(S, p)))) {
    stop("expected_rss_sum must be an S by p matrix.")
  }

  if (is.null(total_obs_sj)) {
    total_obs_sj <- matrix(0, nrow = S, ncol = p)
    for (s in seq_len(S)) {
      for (j in seq_len(p)) {
        total_obs_sj[s, j] <- sum(vapply(
          seq_len(n_s[s]),
          function(i) length(Y[[s]][[i]][[j]]),
          integer(1)))
      }
    }
  }

  component_names <- c(
    "data_likelihood",
    "osullivan_mu", "osullivan_beta",
    "osullivan_phi", "osullivan_psi",
    "score_zeta", "score_xi",
    "spike_slab_shared", "spike_slab_specific",
    "omega_shared", "omega_specific",
    "half_cauchy_eps", "half_cauchy_mu",
    "half_cauchy_beta", "half_cauchy_phi",
    "half_cauchy_psi"
  )
  components <- stats::setNames(rep(0, length(component_names)),
                                component_names)

  # Data likelihood with the direct current E[RSS].
  for (s in seq_len(S)) {
    for (j in seq_len(p)) {
      components["data_likelihood"] <-
        components["data_likelihood"] - 0.5 * (
          total_obs_sj[s, j] * log(2 * pi) +
          total_obs_sj[s, j] * mu_q_log_sigsq_eps[s, j] +
          mu_q_recip_sigsq_eps[s, j] * expected_rss_sum[s, j])
    }
  }

  # O'Sullivan linear/nonlinear prior blocks plus full Gaussian entropy.
  logdet_inv_Sigma_beta <- .elbo_logdet_positive(inv_Sigma_beta)
  logdet_Sigma_beta <- -logdet_inv_Sigma_beta

  for (s in seq_len(S)) {
    for (j in seq_len(p)) {
      components["osullivan_mu"] <- components["osullivan_mu"] +
        .elbo_osullivan_block(
          mu_q_nu_mu[[s]][[j]], Sigma_q_nu_mu[[s]][[j]],
          mu_q_recip_sigsq_mu[s, j], mu_q_log_sigsq_mu[s, j],
           inv_Sigma_beta, logdet_Sigma_beta, K,
           temperature = temperature)
    }
  }

  if (d > 0 && !is.null(mu_q_nu_beta)) {
    for (j in seq_len(p)) {
      for (r in seq_len(d)) {
        components["osullivan_beta"] <- components["osullivan_beta"] +
          .elbo_osullivan_block(
            mu_q_nu_beta[[j]][[r]], Sigma_q_nu_beta[[j]][[r]],
            mu_q_recip_sigsq_beta[j, r],
            mu_q_log_sigsq_beta[j, r],
             inv_Sigma_beta, logdet_Sigma_beta, K,
             temperature = temperature)
      }
    }
  }

  for (l in seq_len(L_f)) {
    for (m in seq_len(M_f[l])) {
      components["osullivan_phi"] <- components["osullivan_phi"] +
        .elbo_osullivan_block(
          mu_q_nu_phi[[l]][, m], Sigma_q_nu_phi[[l]][[m]],
          mu_q_recip_sigsq_phi[[l]][m],
          mu_q_log_sigsq_phi[[l]][m],
           inv_Sigma_beta, logdet_Sigma_beta, K,
           temperature = temperature)
    }
  }

  if (has_specific) {
    for (s in seq_len(S)) {
      for (l in seq_len(L_s[[s]])) {
        for (m in seq_len(M_s[[s]][l])) {
          components["osullivan_psi"] <- components["osullivan_psi"] +
            .elbo_osullivan_block(
              mu_q_nu_psi[[s]][[l]][, m],
              Sigma_q_nu_psi[[s]][[l]][[m]],
              mu_q_recip_sigsq_psi[[s]][[l]][m],
              mu_q_log_sigsq_psi[[s]][[l]][m],
               inv_Sigma_beta, logdet_Sigma_beta, K,
               temperature = temperature)
        }
      }
    }
  }

  # Standard-normal computational scores eta (legacy zeta) and chi (legacy xi).
  for (s in seq_len(S)) {
    for (l in seq_len(L_f)) {
      for (i in seq_len(n_s[s])) {
        components["score_zeta"] <- components["score_zeta"] +
          .elbo_standard_normal_score(
            mu_q_zeta[[s]][[l]][i, ],
             Sigma_q_zeta[[s]][[l]][[i]],
             temperature = temperature)
      }
    }
  }
  if (has_specific) {
    for (s in seq_len(S)) {
      for (l in seq_len(L_s[[s]])) {
        for (i in seq_len(n_s[s])) {
          components["score_xi"] <- components["score_xi"] +
            .elbo_standard_normal_score(
              mu_q_xi[[s]][[l]][i, ],
               Sigma_q_xi[[s]][[l]][[i]],
               temperature = temperature)
        }
      }
    }
  }

  # Shared and study-specific spike-and-slab blocks.
  if (L_f > 0) {
    for (l in seq_len(L_f)) {
      for (j in seq_len(p)) {
        elog_omega <- if (is.matrix(mu_q_log_omega_a)) {
          mu_q_log_omega_a[j, l]
        } else {
          mu_q_log_omega_a[l]
        }
        elog_one_minus <- if (is.matrix(mu_q_log_1_omega_a)) {
          mu_q_log_1_omega_a[j, l]
        } else {
          mu_q_log_1_omega_a[l]
        }
        components["spike_slab_shared"] <-
          components["spike_slab_shared"] +
          .elbo_spike_slab_scalar(
            mu_q_gamma_a[j, l],
            mu_q_normal_a[j, l],
            Sigma_q_normal_a[j, l],
             elog_omega, elog_one_minus,
             temperature = temperature)
      }
    }
  }

  if (has_specific) {
    for (s in seq_len(S)) {
      for (l in seq_len(L_s[[s]])) {
        for (j in seq_len(p)) {
          elog_omega <- if (is.matrix(mu_q_log_omega_b[[s]])) {
            mu_q_log_omega_b[[s]][j, l]
          } else {
            mu_q_log_omega_b[[s]][l]
          }
          elog_one_minus <- if (is.matrix(mu_q_log_1_omega_b[[s]])) {
            mu_q_log_1_omega_b[[s]][j, l]
          } else {
            mu_q_log_1_omega_b[[s]][l]
          }
          components["spike_slab_specific"] <-
            components["spike_slab_specific"] +
            .elbo_spike_slab_scalar(
              mu_q_gamma_b[[s]][j, l],
              mu_q_normal_b[[s]][j, l],
              Sigma_q_normal_b[[s]][j, l],
               elog_omega, elog_one_minus,
               temperature = temperature)
        }
      }
    }
  }

  # Beta priors and Beta entropies for both omega hierarchies.
  if (L_f > 0) {
    if (bool_var_spec_prob) {
      if (!is.matrix(c_1_omega_a) ||
          !isTRUE(all(dim(c_1_omega_a) == c(p, L_f)))) {
        stop("Variable-factor shared omega parameters must be p by L_f.")
      }
      for (l in seq_len(L_f)) {
        for (j in seq_len(p)) {
          components["omega_shared"] <- components["omega_shared"] +
            .elbo_beta_prior_entropy(
               c_1_omega_a[j, l], d_1_omega_a[j, l], c_0, d_0,
               temperature = temperature)
        }
      }
    } else {
      if (is.matrix(c_1_omega_a) || length(c_1_omega_a) != L_f) {
        stop("Factor-level shared omega parameters must have length L_f.")
      }
      for (l in seq_len(L_f)) {
        components["omega_shared"] <- components["omega_shared"] +
          .elbo_beta_prior_entropy(
             c_1_omega_a[l], d_1_omega_a[l], c_0, d_0,
             temperature = temperature)
      }
    }
  }

  if (has_specific) {
    for (s in seq_len(S)) {
      L_ss <- L_s[[s]]
      if (bool_var_spec_prob) {
        if (!is.matrix(c_1_omega_b[[s]]) ||
            !isTRUE(all(dim(c_1_omega_b[[s]]) == c(p, L_ss)))) {
          stop("Variable-factor specific omega parameters in study ", s,
               " must be p by L_s[s].")
        }
        for (l in seq_len(L_ss)) {
          for (j in seq_len(p)) {
            components["omega_specific"] <- components["omega_specific"] +
              .elbo_beta_prior_entropy(
                c_1_omega_b[[s]][j, l],
                 d_1_omega_b[[s]][j, l], c_0, d_0,
                 temperature = temperature)
          }
        }
      } else {
        if (is.matrix(c_1_omega_b[[s]]) ||
            length(c_1_omega_b[[s]]) != L_ss) {
          stop("Factor-level specific omega parameters in study ", s,
               " must have length L_s[s].")
        }
        for (l in seq_len(L_ss)) {
          components["omega_specific"] <- components["omega_specific"] +
            .elbo_beta_prior_entropy(
              c_1_omega_b[[s]][l],
               d_1_omega_b[[s]][l], c_0, d_0,
               temperature = temperature)
        }
      }
    }
  }

  # Complete IG--IG Half-Cauchy hierarchies. These terms are omitted only by
  # the explicit pre-variance conditional diagnostic, before the first full
  # IG variational distributions exist.
  if (include_half_cauchy && is.null(lambda_q_a_eps)) {
    lambda_q_a_eps <- kappa_q_a / mu_q_recip_a_eps
  }
  if (include_half_cauchy && is.null(lambda_q_a_mu)) {
    lambda_q_a_mu <- kappa_q_a / mu_q_recip_a_mu
  }
  if (include_half_cauchy && d > 0 && is.null(lambda_q_a_beta)) {
    lambda_q_a_beta <- kappa_q_a / mu_q_recip_a_beta
  }
  if (include_half_cauchy && is.null(lambda_q_a_phi)) {
    lambda_q_a_phi <- lapply(mu_q_recip_a_phi,
                              function(x) kappa_q_a / x)
  }
  if (include_half_cauchy && has_specific && is.null(lambda_q_a_psi)) {
    lambda_q_a_psi <- lapply(mu_q_recip_a_psi, function(study) {
      lapply(study, function(x) kappa_q_a / x)
    })
  }

  if (include_half_cauchy) {
    for (s in seq_len(S)) {
      for (j in seq_len(p)) {
        components["half_cauchy_eps"] <- components["half_cauchy_eps"] +
          .elbo_half_cauchy_pair(
            kappa_q_sigsq_eps[s, j], lambda_q_sigsq_eps[s, j],
            kappa_q_a, lambda_q_a_eps[s, j], A,
            temperature = temperature)
        components["half_cauchy_mu"] <- components["half_cauchy_mu"] +
          .elbo_half_cauchy_pair(
            kappa_q_sigsq_mu, lambda_q_sigsq_mu[s, j],
            kappa_q_a, lambda_q_a_mu[s, j], A,
            temperature = temperature)
      }
    }

    if (d > 0) {
      for (j in seq_len(p)) {
        for (r in seq_len(d)) {
          components["half_cauchy_beta"] <-
            components["half_cauchy_beta"] +
            .elbo_half_cauchy_pair(
              kappa_q_sigsq_beta, lambda_q_sigsq_beta[j, r],
              kappa_q_a, lambda_q_a_beta[j, r], A,
              temperature = temperature)
        }
      }
    }

    for (l in seq_len(L_f)) {
      for (m in seq_len(M_f[l])) {
        components["half_cauchy_phi"] <- components["half_cauchy_phi"] +
          .elbo_half_cauchy_pair(
            kappa_q_sigsq_phi, lambda_q_sigsq_phi[[l]][m],
            kappa_q_a, lambda_q_a_phi[[l]][m], A,
            temperature = temperature)
      }
    }

    if (has_specific) {
      for (s in seq_len(S)) {
        for (l in seq_len(L_s[[s]])) {
          for (m in seq_len(M_s[[s]][l])) {
            components["half_cauchy_psi"] <-
              components["half_cauchy_psi"] +
              .elbo_half_cauchy_pair(
                kappa_q_sigsq_psi,
                lambda_q_sigsq_psi[[s]][[l]][m],
                kappa_q_a, lambda_q_a_psi[[s]][[l]][m], A,
                temperature = temperature)
          }
        }
      }
    }
  }

  total <- sum(components)
  component_finite <- stats::setNames(is.finite(components),
                                      names(components))
  finite_check <- list(
    all = is.finite(total) && all(component_finite),
    total = is.finite(total),
    components = component_finite,
    nonfinite_components = names(component_finite)[!component_finite]
  )
  diagnostics <- .elbo_history_diagnostics(
    c(elbo_history, total), tolerance = monotonicity_tolerance)

  result <- list(
    total = as.numeric(total),
    components = components,
    finite = finite_check,
    diagnostics = diagnostics,
    rss = list(
      source = rss_source,
      expected_rss = expected_rss,
      expected_rss_sum = expected_rss_sum
    ),
    temperature = temperature,
    inverse_temperature = c_val,
    scope = objective_scope,
    omega_hierarchy = if (bool_var_spec_prob) {
      "variable_factor"
    } else {
      "factor"
    }
  )
  class(result) <- c("multiFSYNC_variational_objective", "list")
  result
}

#' Compute the complete ordinary multi-study ELBO
#'
#' @keywords internal
compute_elbo_multi <- function(...) {
  args <- list(...)
  c_now <- if (is.null(args$c_val)) 1 else args$c_val
  if (!isTRUE(all.equal(c_now, 1))) {
    stop("compute_elbo_multi() computes the ordinary T = 1 ELBO only.")
  }
  args$c_val <- 1
  result <- do.call(.compute_variational_objective_multi, args)
  result$scope <- "full"
  result$objective <- "ordinary_elbo"
  class(result) <- c("multiFSYNC_elbo", class(result))
  result
}

#' Compute the complete annealed multi-study variational objective
#'
#' This diagnostic is E_q log p(Y,Theta) - T E_q log q(Theta) under the fixed
#' reference measures in derivations.md. It never replaces the T = 1 stopping
#' objective.
#'
#' @keywords internal
compute_tempered_objective_multi <- function(..., temperature = NULL,
                                              c_val = NULL) {
  if (is.null(temperature) && is.null(c_val)) {
    stop("Supply temperature or c_val.")
  }
  if (!is.null(temperature)) {
    if (length(temperature) != 1L || !is.finite(temperature) ||
        temperature <= 0) {
      stop("temperature must be one finite positive number.")
    }
    implied_c <- 1 / temperature
    if (!is.null(c_val) && !isTRUE(all.equal(c_val, implied_c))) {
      stop("temperature and c_val are inconsistent.")
    }
    c_val <- implied_c
  }
  if (length(c_val) != 1L || !is.finite(c_val) || c_val <= 0) {
    stop("c_val must be one finite positive inverse temperature.")
  }
  args <- list(...)
  args$c_val <- c_val
  result <- do.call(.compute_variational_objective_multi, args)
  result$objective <- "tempered_variational_objective"
  class(result) <- c("multiFSYNC_tempered_objective", class(result))
  result
}
