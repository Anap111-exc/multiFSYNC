# =============================================================================
# utils_multi.R — Utility functions for multiFSYNC
#
# Based on code from bayesSYNC (https://github.com/hruffieux/bayesSYNC)
# Original authors: Salima Jaoua, Daniel Temko, Helene Ruffieux
# License: GPL (>= 3)
#
# This file aggregates utility functions from:
#   - bayesSYNC/R/utils.R
#   - bayesSYNC/R/OSullivan_splines.R
#   - bayesSYNC/R/set_hyper.R
#
# The code has been adapted for multi-study functional factor model use.
# =============================================================================

# ---- from bayesSYNC/R/utils.R ----

create_named_list <- function(...) {
  nms <- names(match.call()[-1])
  if (is.null(nms)) {
    nms <- as.character(match.call()[-1])
  } else {
    no_name <- nms == ""
    if (any(no_name)) nms[no_name] <- as.character(match.call()[-1])[no_name]
  }
  setNames(list(...), nms)
}

all_same <- function(x) length(unique(x)) == 1

# Study-specific factor counts may be supplied either as the historical scalar
# upper dimension or as one count per study.  Internally multiFSYNC always uses
# the latter representation so that every jagged parameter block has an
# unambiguous dimension.
.normalize_L_s <- function(L_s, S, name = "L_s") {
  if (!is.numeric(S) || length(S) != 1L || !is.finite(S) ||
      abs(S - round(S)) > .Machine$double.eps^0.5 || S < 1L) {
    stop("S must be a positive integer before ", name, " is normalised.")
  }
  S <- as.integer(S)
  if (!is.numeric(L_s) || !is.null(dim(L_s)) ||
      !(length(L_s) %in% c(1L, S)) || any(!is.finite(L_s)) ||
      any(abs(L_s - round(L_s)) > .Machine$double.eps^0.5) ||
      any(L_s < 0)) {
    stop(name, " must be a non-negative integer scalar or a length-S vector.")
  }
  if (length(L_s) == 1L) L_s <- rep(L_s, S)
  as.integer(L_s)
}

.L_s_at <- function(L_s, s) {
  as.integer(if (length(L_s) == 1L) L_s[[1L]] else L_s[[s]])
}

.has_specific <- function(L_s) {
  length(L_s) > 0L && any(L_s > 0L)
}

.max_L_s <- function(L_s) {
  if (length(L_s)) as.integer(max(L_s)) else 0L
}

# Preserve the historical scalar field when all studies use the same count;
# L_s_by_study is added separately by public fit/generator return objects.
.compact_L_s <- function(L_s) {
  L_s <- as.integer(L_s)
  if (length(L_s) && all(L_s == L_s[[1L]])) L_s[[1L]] else L_s
}

# Small dependency-free block-diagonal constructor used by the O'Sullivan
# two-block priors.  This replaces pracma::blkdiag and keeps the package's
# runtime dependency set minimal.
blkdiag <- function(...) {
  mats <- list(...)
  if (length(mats) == 0L) return(matrix(numeric(0), 0L, 0L))
  if (any(!vapply(mats, is.matrix, logical(1)))) {
    stop("Every blkdiag argument must be a matrix.")
  }
  nr <- vapply(mats, nrow, integer(1))
  nc <- vapply(mats, ncol, integer(1))
  ans <- matrix(0, sum(nr), sum(nc))
  r0 <- cumsum(c(0L, nr))
  c0 <- cumsum(c(0L, nc))
  for (k in seq_along(mats)) {
    ans[(r0[k] + 1L):r0[k + 1L],
        (c0[k] + 1L):c0[k + 1L]] <- mats[[k]]
  }
  ans
}

check_natural <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}

check_positive <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste0(deparse(substitute(x)), " must be positive (larger than precision zero in R).")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_zero_one <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste0(deparse(substitute(x)), " must lie between 0 and 1.")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_structure <- function(x, struct, type, size = NULL,
                            null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0)
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0)
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
                size_mess, inf_mess, na_mess, null_mess, "."))
  }
}

#' Create a temporal-grid objects for use in FPCA algorithms.
#'
#' @param time_obs Vector or list of vectors containing time of observations.
#' @param K Number of O'Sullivan spline functions.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector.
#' @param int_knots Position of interior knots.
#' @param format_univ Boolean indicating whether the univariate format is used.
#'
#' @return A list containing C, n_g, time_g, C_g.
#'
#' @noRd
#' @export
#'
get_grid_objects <- function(time_obs, K, n_g = 1000, time_g = NULL,
                             int_knots = NULL,
                             format_univ = FALSE) {

  if (!is.list(time_obs) || length(time_obs) == 0L) {
    stop("time_obs must be a non-empty list of observation-time vectors.")
  }
  if (!is.null(K) && (any(!is.finite(K)) || any(K != as.integer(K)) ||
                      any(K < 2L))) {
    stop("Every supplied K must be an integer of at least 2.")
  }
  all_times <- unlist(time_obs, recursive = TRUE, use.names = FALSE)
  if (!is.numeric(all_times) || any(!is.finite(all_times)) ||
      any(all_times < 0 | all_times > 1)) {
    stop("All observation times must be finite and pre-normalised to [0,1].")
  }

  if (is.null(int_knots)) {
    if(format_univ) {
      if (is.null(K)) {
        K <- max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i))/4), 40)), 7)
      }
      unique_time_obs <- sort(unique(Reduce(c, time_obs)))
      int_knots <- quantile(unique_time_obs, seq(0, 1, length=K)[-c(1,K)])
    } else {
      p <- length(time_obs[[1]])
      if (is.null(K)) {
        K <- sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7))
      } else {
        check_structure(K, "vector", "numeric", c(1, p))
        if (length(K)==1) K <- rep(K, p)
      }
      unique_time_obs <- unname(sort(unlist(time_obs)))
      int_knots <- lapply(
        K,
        function(x) quantile(unique_time_obs, seq(0, 1, length = x)[-c(1, x)])
      )
    }
  }

  N <- length(time_obs)

  if (is.null(time_g)) {
    if (is.null(n_g) || length(n_g) != 1L || !is.finite(n_g) ||
        !is_int(n_g) || n_g < 2L) {
      stop("n_g must be an integer of at least 2 when time_g is not supplied.")
    }
    n_g <- as.integer(n_g)
    time_g <- seq(0, 1, length.out = n_g)
  } else {
    if (!is.numeric(time_g) || length(time_g) < 2L ||
        any(!is.finite(time_g)) || any(time_g < 0 | time_g > 1) ||
        any(diff(time_g) <= 0)) {
      stop("time_g must contain at least two finite, unique, strictly increasing points in [0,1].")
    }
    n_g <- length(time_g)
  }

  C <- vector("list", length=N)

  if (format_univ) {
    for(i in 1:N) {
      X <- X_design(time_obs[[i]])
      Z <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
      C[[i]] <- cbind(X, Z)
    }
    X_g <- X_design(time_g)
    Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
    C_g <- cbind(X_g, Z_g)
  } else {
    p <- length(time_obs[[1]])
    for(i in 1:N) {
      C[[i]] <- vector("list", length = p)
      for(j in 1:p) {
        X <- X_design(time_obs[[i]][[j]])
        Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots[[j]])
        C[[i]][[j]] <- cbind(X, Z)
      }
    }
    X_g <- X_design(time_g)
    Z_g <- lapply(
      int_knots,
      function(x) ZOSull(time_g, range.x = c(0, 1), intKnots = x)
    )
    C_g <- lapply(Z_g, function(Z) cbind(X_g, Z))
  }

  create_named_list(C, n_g, time_g, C_g)
}

X_design <- function(x) {
  if(is.list(x)) {
    x <- do.call(cbind, x)
  }
  X <- cbind(1, x)
  return(X)
}

is_int <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

tr <- function(X) {
  if(nrow(X)!=ncol(X)) stop("X must be a square matrix.")
  ans <- sum(diag(X))
  return(ans)
}

# ---- Positive-definite linear algebra -------------------------------------
#
# CAVI precision matrices are theoretically positive definite. Adding a fixed
# ridge before every solve changes the objective even when no stabilisation is
# needed. Try the unmodified matrix first and add scale-adaptive jitter only
# after Cholesky failure.

.multiFSYNC_spd_state <- new.env(parent = emptyenv())

.reset_spd_diagnostics <- function() {
  .multiFSYNC_spd_state$total_calls <- 0L
  .multiFSYNC_spd_state$jitter_count <- 0L
  .multiFSYNC_spd_state$events <- data.frame(
    context = character(),
    jitter = numeric(),
    relative_jitter = numeric(),
    attempts = integer(),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

.reset_spd_diagnostics()

.record_spd_diagnostic <- function(context, jitter, relative_jitter, attempts) {
  .multiFSYNC_spd_state$total_calls <-
    .multiFSYNC_spd_state$total_calls + 1L
  if (jitter > 0) {
    .multiFSYNC_spd_state$jitter_count <-
      .multiFSYNC_spd_state$jitter_count + 1L
    .multiFSYNC_spd_state$events <- rbind(
      .multiFSYNC_spd_state$events,
      data.frame(
        context = as.character(context),
        jitter = as.numeric(jitter),
        relative_jitter = as.numeric(relative_jitter),
        attempts = as.integer(attempts),
        stringsAsFactors = FALSE
      )
    )
  }
  invisible(NULL)
}

.get_spd_diagnostics <- function() {
  events <- .multiFSYNC_spd_state$events
  list(
    total_calls = .multiFSYNC_spd_state$total_calls,
    jitter_count = .multiFSYNC_spd_state$jitter_count,
    used_jitter = .multiFSYNC_spd_state$jitter_count > 0L,
    max_jitter = if (nrow(events)) max(events$jitter) else 0,
    events = events
  )
}

.inverse_spd <- function(precision, context = "unspecified",
                         jitter_relative = c(1e-12, 1e-10, 1e-8,
                                             1e-6, 1e-4)) {
  if (!is.matrix(precision) || nrow(precision) != ncol(precision) ||
      !is.numeric(precision) || any(!is.finite(precision))) {
    stop(context, ": precision must be a finite numeric square matrix.")
  }
  if (!length(jitter_relative) || any(!is.finite(jitter_relative)) ||
      any(jitter_relative <= 0) ||
      is.unsorted(jitter_relative, strictly = TRUE)) {
    stop("jitter_relative must be a strictly increasing positive vector.")
  }

  asymmetry <- max(abs(precision - t(precision)))
  matrix_scale <- max(1, max(abs(precision)))
  if (asymmetry > 1e-10 * matrix_scale) {
    stop(context, ": precision is not numerically symmetric.")
  }
  precision_sym <- (precision + t(precision)) / 2

  chol_factor <- tryCatch(chol(precision_sym), error = function(e) NULL)
  jitter <- 0
  relative_jitter <- 0
  attempts <- 1L

  if (is.null(chol_factor)) {
    diagonal_scale <- max(1, max(abs(diag(precision_sym))))
    for (relative_candidate in jitter_relative) {
      attempts <- attempts + 1L
      jitter_candidate <- diagonal_scale * relative_candidate
      chol_factor <- tryCatch(
        chol(precision_sym + jitter_candidate * diag(nrow(precision_sym))),
        error = function(e) NULL
      )
      if (!is.null(chol_factor)) {
        jitter <- jitter_candidate
        relative_jitter <- relative_candidate
        break
      }
    }
  }

  if (is.null(chol_factor)) {
    stop(context, ": Cholesky factorisation failed after adaptive jitter up to ",
         format(max(jitter_relative)), " times the diagonal scale.")
  }

  inverse <- chol2inv(chol_factor)
  inverse <- (inverse + t(inverse)) / 2
  attr(inverse, "solver_diagnostic") <- list(
    context = as.character(context),
    used_jitter = jitter > 0,
    jitter = jitter,
    relative_jitter = relative_jitter,
    attempts = attempts
  )
  .record_spd_diagnostic(context, jitter, relative_jitter, attempts)
  inverse
}

# Construct loadings that satisfy the population Gram conditions used in the
# thesis simulations.  Sparse loadings use disjoint supports; dense loadings
# use orthonormal columns.  In both cases [A, B_s]'[A, B_s] is diagonal and the
# reference-study column norms are strictly decreasing.
generate_identified_loadings <- function(p, L_f, L_s, S,
                                         sparse = TRUE,
                                         prop_sparse = 0.5) {
  L_s <- .normalize_L_s(L_s, S)
  L_s_max <- .max_L_s(L_s)
  has_specific <- .has_specific(L_s)
  q <- L_f + L_s_max
  if (q == 0L) {
    return(list(a = matrix(0, p, 0L), gamma_a = matrix(0, p, 0L),
                b = NULL, gamma_b = NULL))
  }
  if (p <= q) {
    stop("identified_loadings = TRUE requires p > L_f + max(L_s).")
  }
  target_norms <- seq(1.6, 0.8, length.out = q)
  A <- matrix(0, p, L_f)
  B <- if (has_specific) vector("list", S) else NULL

  if (sparse) {
    groups <- split(seq_len(p), rep(seq_len(q), length.out = p))
    make_on_support <- function(idx, target) {
      desired <- max(1L, min(length(idx),
        as.integer(round((1 - prop_sparse) * p))))
      active <- sort(sample(idx, desired))
      value <- rnorm(desired)
      value <- target * value / sqrt(sum(value^2))
      list(index = active, value = value)
    }
    if (L_f > 0L) {
      for (l in seq_len(L_f)) {
        z <- make_on_support(groups[[l]], target_norms[l])
        A[z$index, l] <- z$value
      }
    }
    if (has_specific) {
      for (s in seq_len(S)) {
        L_ss <- L_s[[s]]
        B[[s]] <- matrix(0, p, L_ss)
        for (l in seq_len(L_ss)) {
          idx <- groups[[L_f + l]]
          z <- make_on_support(idx, target_norms[L_f + l])
          # Study-specific random values on the same disjoint support preserve
          # Gram orthogonality while making B_s genuinely study dependent.
          B[[s]][z$index, l] <- z$value
        }
      }
    }
  } else {
    Q_ref <- qr.Q(qr(matrix(rnorm(p * q), p, q)))
    if (L_f > 0L) {
      A <- sweep(Q_ref[, seq_len(L_f), drop = FALSE], 2,
                 target_norms[seq_len(L_f)], "*")
    }
    if (has_specific) {
      for (s in seq_len(S)) {
        L_ss <- L_s[[s]]
        if (L_ss == 0L) {
          B[[s]] <- matrix(0, p, 0L)
          next
        }
        if (s == 1L) {
          QB <- Q_ref[, L_f + seq_len(L_ss), drop = FALSE]
        } else {
          raw <- matrix(rnorm(p * L_ss), p, L_ss)
          if (L_f > 0L) {
            QA <- sweep(A, 2, sqrt(colSums(A^2)), "/")
            raw <- raw - QA %*% crossprod(QA, raw)
          }
          QB <- qr.Q(qr(raw))[, seq_len(L_ss), drop = FALSE]
        }
        B[[s]] <- sweep(QB, 2, target_norms[L_f + seq_len(L_ss)], "*")
      }
    }
  }

  list(
    a = A,
    gamma_a = 1 * (A != 0),
    b = B,
    gamma_b = if (has_specific) lapply(B, function(x) 1 * (x != 0)) else NULL
  )
}

cprod <- function(x, y) {
  if(missing(y)) {
    if(!is.vector(x)) {
      stop("Use the crossprod function for matrix inner products")
    }
    y <- x
  }
  if(!is.vector(y) & !is.vector(x)) {
    stop("Use the crossprod function for matrix inner products")
  }
  ans <- as.vector(crossprod(x, y))
  return(ans)
}

normalise <- function(x) {
  ans <- x/sqrt(cprod(x))
  return(ans)
}

E_cprod <- function(mean_1, Cov_21, mean_2, A) {
  if(missing(A)) {
    A <- diag(length(mean_1))
  }
  tr_term <- tr(Cov_21 %*% A)
  cprod_term <- cprod(mean_1, A %*% mean_2)
  ans <- tr_term + cprod_term
  return(ans)
}

E_h <- function(L, mean_1, Cov_21, mean_2, A) {
  if(missing(A)) {
    A <- diag(length(mean_1))
  }
  d_1 <- length(mean_1)
  inds_1 <- matrix(1:d_1, ncol = L)
  ans <- rep(NA, L)
  for(l in 1:L) {
    mean_1_l <- mean_1[inds_1[, l]]
    Cov_21_l <- Cov_21[, inds_1[, l]]
    ans[l] <- E_cprod(mean_1_l, Cov_21_l, mean_2, A)
  }
  return(ans)
}

E_H <- function(L_1, L_2, mean_1, Cov_21, mean_2, A) {
  if(missing(A)) {
    A <- diag(length(mean_1))
  }
  d_1 <- length(mean_1)
  d_2 <- length(mean_2)
  L <- L_1 + L_2
  inds_1 <- matrix(1:d_1, ncol = L_1)
  inds_2 <- matrix(1:d_2, ncol = L_2)
  ans <- matrix(NA, L_1, L_2)
  for(l in 1:L_1) {
    mean_1_l <- mean_1[inds_1[, l]]
    for(k in 1:L_2) {
      mean_2_k <- mean_2[inds_2[, k]]
      Cov_21_kl <- Cov_21[inds_2[, k], inds_1[, l]]
      ans[l, k] <- E_cprod(mean_1_l, Cov_21_kl, mean_2_k, A)
    }
  }
  return(ans)
}

vec <- function(A) {
  return(as.vector(A))
}

vecInverse <- function(a) {
  is.wholenumber <- function(x,tol=sqrt(.Machine$double.eps))
    return(abs(x-round(x))<tol)
  a <- as.vector(a)
  if (!is.wholenumber(sqrt(length(a))))
    stop("input vector must be a perfect square in length")
  dmnVal <- round(sqrt(length(a)))
  A <- matrix(NA,dmnVal,dmnVal)
  for (j in 1:dmnVal)
    A[,j] <- a[((j-1)*dmnVal+1):(j*dmnVal)]
  return(A)
}

#' Function performing trapesoidal integration.
#'
#' @param xgrid Grid.
#' @param fgrid Function on the grid.
#' @return Integration result.
#'
#' @noRd
#' @export
#'
trapint <- function(xgrid,fgrid) {
  ng <- length(xgrid)
  xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
  fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]
  integ <- sum(xvec*fvec)/2
  return(integ)
}

wait <- function() {
  cat("Hit Enter to continue\n")
  ans <- readline()
  invisible()
}

#' Flip the sign of eigenfunctions and corresponding scores
#'
#' @param vec_flip Vector of size L.
#' @param list_Psi_hat List or matrix of eigenfunctions.
#' @param Zeta_hat Matrix of estimated scores.
#' @param zeta_ellipse 95% posterior credible boundaries.
#'
#' @return An object containing the eigenfunctions, scores and credible boundaries.
#'
#' @noRd
#' @export
#'
flip_sign <- function(vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse = NULL) {
  stopifnot(all(vec_flip == 1 | vec_flip == -1))
  list_Psi_hat[,seq_along(vec_flip)] <- sapply(seq_along(vec_flip), function(ll) {
    list_Psi_hat[,ll] * vec_flip[ll]
  })
  Zeta_hat[,seq_along(vec_flip)] <- sweep(Zeta_hat[,seq_along(vec_flip)], 2, vec_flip, "*")
  if (!is.null(zeta_ellipse)) {
    zeta_ellipse <- lapply(zeta_ellipse, function(zeta_ellipse_subj) sweep(zeta_ellipse_subj, 2, vec_flip, "*"))
  }
  create_named_list(list_Psi_hat, Zeta_hat, zeta_ellipse)
}

frobenius_norm <- function(A, B) {
  sqrt(sum((A - B)^2))
}

.ordered_subsets <- function(n, r) {
  if (r < 0L || r > n) stop("r must satisfy 0 <= r <= n.")
  if (r == 0L) return(matrix(integer(0), 1L, 0L))
  recurse <- function(values, k) {
    if (k == 1L) return(matrix(values, ncol = 1L))
    do.call(rbind, lapply(values, function(first) {
      rest <- recurse(values[values != first], k - 1L)
      cbind(first, rest)
    }))
  }
  recurse(seq_len(n), r)
}

# Solve a rectangular minimum-cost one-to-one assignment.  This is the
# Hungarian algorithm with rows assigned to distinct columns; nrow(cost) must
# not exceed ncol(cost).  Keeping it internal avoids an additional package
# dependency in simulation-only evaluation code.
.hungarian_assignment <- function(cost) {
  cost <- as.matrix(cost)
  n <- nrow(cost)
  m <- ncol(cost)
  if (n == 0L) return(integer(0))
  if (n > m) stop("Hungarian assignment requires no more rows than columns.")
  if (any(!is.finite(cost))) {
    finite_cost <- cost[is.finite(cost)]
    replacement <- if (length(finite_cost)) max(finite_cost) + 1e6 else 1e6
    cost[!is.finite(cost)] <- replacement
  }

  # Index 1 is the dummy zero column used by the standard primal-dual form.
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L)
  way <- integer(m + 1L)

  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 1L
    minv <- rep(Inf, m)
    used <- rep(FALSE, m + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- NA_integer_
      for (j in seq_len(m)) {
        jj <- j + 1L
        if (!used[jj]) {
          cur <- cost[i0, j] - u[i0 + 1L] - v[jj]
          if (cur < minv[j]) {
            minv[j] <- cur
            way[jj] <- j0
          }
          if (minv[j] < delta) {
            delta <- minv[j]
            j1 <- jj
          }
        }
      }
      if (!is.finite(delta) || is.na(j1)) stop("Hungarian assignment failed.")
      used_idx <- which(used)
      for (jj in used_idx) {
        u[p[jj] + 1L] <- u[p[jj] + 1L] + delta
        v[jj] <- v[jj] - delta
      }
      free_cols <- which(!used[-1L])
      minv[free_cols] <- minv[free_cols] - delta
      j0 <- j1
      if (p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }

  assignment <- integer(n)
  for (j in seq_len(m)) {
    if (p[j + 1L] > 0L) assignment[p[j + 1L]] <- j
  }
  if (any(assignment == 0L)) stop("Hungarian assignment returned an incomplete match.")
  assignment
}

.safe_correlation <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L ||
      any(!is.finite(x)) || any(!is.finite(y))) return(0)
  sx <- stats::sd(x)
  sy <- stats::sd(y)
  if (is.finite(sx) && is.finite(sy) && sx > 0 && sy > 0) {
    out <- suppressWarnings(stats::cor(x, y))
    if (is.finite(out)) return(out)
  }
  denom <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denom) || denom <= 0) return(0)
  sum(x * y) / denom
}

#' @noRd
#' @export
match_factor_and_sign <- function(B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat,
                                  list_list_Phi_hat, list_cumulated_pve,
                                  list_Cov_zeta_hat, list_h_hat = NULL,
                                  list_var_vec = NULL) {
  if (!is.null(list_h_hat)) {
    stopifnot(!is.null(list_var_vec))
  }
  N <- nrow(Zeta[[1]])
  Q <- ncol(B_hat)
  Q_true <- ncol(B)

  perm_factor <- match_factors(B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat,
                               list_list_Phi_hat, list_cumulated_pve,
                               list_Cov_zeta_hat, list_h_hat, list_var_vec)

  best_perm <- perm_factor$best_perm
  perm_sign_factor <- perm_factor$perm_sign_factor
  perm_list_cumulated_pve <- perm_factor$perm_list_cumulated_pve
  perm_list_Cov_zeta_hat <- perm_factor$perm_list_Cov_zeta_hat

  perm_B_hat <- perm_factor$perm_B_hat
  perm_ppi <- perm_factor$perm_ppi
  perm_factor_ppi <- perm_factor$perm_factor_ppi
  perm_list_Zeta_hat <- perm_factor$perm_list_Zeta_hat
  perm_list_list_Phi_hat <- perm_factor$perm_list_list_Phi_hat

  perm_B_hat_untrimmed <- perm_factor$perm_B_hat_untrimmed
  perm_ppi_untrimmed <- perm_factor$perm_ppi_untrimmed
  perm_list_Zeta_hat_untrimmed <- perm_factor$perm_list_Zeta_hat_untrimmed
  perm_list_list_Phi_hat_untrimmed <- perm_factor$perm_list_list_Phi_hat_untrimmed

  perm_sign <- match_sign_components(Zeta, perm_list_Zeta_hat, perm_list_list_Phi_hat)

  perm_sign_fpca <- perm_sign$perm_sign_fpca
  perm_component <- perm_sign$perm_component

  perm_list_Zeta_hat <- perm_sign$perm_list_Zeta_hat
  perm_list_list_Phi_hat <- perm_sign$perm_list_list_Phi_hat
  for (q in seq_len(Q_true)) {
    perm_list_Zeta_hat_untrimmed[[q]] <- perm_list_Zeta_hat[[q]]
    perm_list_list_Phi_hat_untrimmed[[q]] <- perm_list_list_Phi_hat[[q]]
  }

  # Apply the same component permutation/sign transformation to posterior
  # score covariance matrices used by simulation diagnostics.
  for (q in seq_len(Q_true)) {
    Lq <- ncol(Zeta[[q]])
    idx <- perm_component[q, seq_len(Lq)]
    signs <- perm_sign_fpca[q, seq_len(Lq)]
    perm_list_Cov_zeta_hat[[q]] <- lapply(
      perm_list_Cov_zeta_hat[[q]],
      function(V) {
        V_new <- V[idx, idx, drop = FALSE]
        V_new * tcrossprod(signs)
      })
  }

  perm_list_h_hat <- perm_factor$perm_list_h_hat
  perm_list_h_hat_untrimmed <- perm_factor$perm_list_h_hat_untrimmed
  perm_list_var_vec <- perm_factor$perm_list_var_vec

  if (!is.null(perm_list_h_hat)) {
    perm_list_h_low <- perm_list_h_upp <- lapply(1:N, function(i) vector("list", Q_true))
    perm_list_h_low_untrimmed <-  perm_list_h_upp_untrimmed <- lapply(1:N, function(i) vector("list", Q))
  } else {
    perm_list_h_low <- perm_list_h_upp <- perm_list_h_low_untrimmed <-  perm_list_h_upp_untrimmed <- NULL
  }

  bool_rescale_loadings_and_scores <- T
  if (bool_rescale_loadings_and_scores) {
    norm_col_B <- sqrt(colSums(B^2))
    norm_col_B_hat <- sqrt(colSums(perm_B_hat^2))
    for (q in 1:Q_true) {
      if (!is.finite(norm_col_B[q]) || norm_col_B[q] <= 0 ||
          !is.finite(norm_col_B_hat[q]) || norm_col_B_hat[q] <= 0) {
        stop("Cannot rescale a factor with a zero or non-finite loading norm.")
      }
      perm_B_hat[,q] <- perm_B_hat_untrimmed[,q] <- perm_B_hat[,q] * norm_col_B[q] / norm_col_B_hat[q]
      perm_list_Zeta_hat[[q]] <- perm_list_Zeta_hat_untrimmed[[q]] <- perm_list_Zeta_hat[[q]] * norm_col_B_hat[q] / norm_col_B[q]
      perm_list_Cov_zeta_hat[[q]] <- lapply(perm_list_Cov_zeta_hat[[q]], function(perm_list_Cov_zeta_hat_q_i)  perm_list_Cov_zeta_hat_q_i * norm_col_B_hat[q]^2 / norm_col_B[q]^2)
      if (!is.null(perm_list_h_hat)) {
        for (i in 1:N) {
          perm_list_h_hat[[i]][[q]] <- perm_list_h_hat_untrimmed[[i]][[q]] <- perm_list_h_hat[[i]][[q]] * norm_col_B_hat[q] / norm_col_B[q]
          sd_i_q <- sqrt(perm_list_var_vec[[q]][[i]]) * norm_col_B_hat[q] / norm_col_B[q]
          perm_list_h_low[[i]][[q]] <- perm_list_h_low_untrimmed[[i]][[q]] <- perm_list_h_hat[[i]][[q]] + qnorm(0.025) * sd_i_q
          perm_list_h_upp[[i]][[q]] <- perm_list_h_upp_untrimmed[[i]][[q]] <- perm_list_h_hat[[i]][[q]] + qnorm(0.975) * sd_i_q
        }
      }
    }
  }

  create_named_list(best_perm, perm_sign_factor, perm_sign_fpca, perm_component,
                    perm_factor_ppi,
                    perm_list_cumulated_pve, perm_list_Cov_zeta_hat,
                    perm_B_hat, perm_ppi, perm_list_Zeta_hat, perm_list_list_Phi_hat,
                    perm_B_hat_untrimmed, perm_ppi_untrimmed,
                    perm_list_Zeta_hat_untrimmed, perm_list_list_Phi_hat_untrimmed,
                    perm_list_h_hat,
                    perm_list_h_low,
                    perm_list_h_upp,
                    perm_list_h_hat_untrimmed,
                    perm_list_h_low_untrimmed,
                    perm_list_h_upp_untrimmed,
                    perm_list_var_vec)
}

match_factors <- function(B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat,
                          list_list_Phi_hat, list_cumulated_pve, list_Cov_zeta_hat,
                          list_h_hat = NULL, list_var_vec = NULL) {
  N <- nrow(Zeta[[1]])
  Q_true <- ncol(B)
  Q <- ncol(B_hat)

  if (Q_true > Q) {
    stop("The number of estimated factors, Q, must be larger than the number of true factors.")
  } else if (Q_true < Q) {
    warning("Dropping superfluous factors based on the true loadings.")
  }
  abs_correlation <- outer(seq_len(Q_true), seq_len(Q), Vectorize(function(q, q_hat) {
    abs(.safe_correlation(B[, q], B_hat[, q_hat]))
  }))
  best_perm <- .hungarian_assignment(1 - abs_correlation)

  perm_ppi <- ppi[, best_perm, drop = F]
  perm_B_hat <- B_hat[, best_perm, drop = F]

  perm_sign_factor <- rep(NA, Q_true)
  perm_list_Zeta_hat <- perm_list_list_Phi_hat <- vector("list", Q_true)
  perm_list_Zeta_hat_untrimmed <- perm_list_list_Phi_hat_untrimmed <- vector("list", Q)

  if (!is.null(list_h_hat)) {
    perm_list_h_hat <- lapply(1:N, function(i) vector("list", Q_true))
    perm_list_h_hat_untrimmed <- lapply(1:N, function(i) vector("list", Q))
  } else {
    perm_list_h_hat <- perm_list_h_hat_untrimmed <- NULL
  }

  for (q in 1:Q) {
    if (q <= Q_true) {
      loading_cor <- .safe_correlation(B[, q], perm_B_hat[, q])
      perm_sign_factor[q] <- ifelse(loading_cor < 0, -1, 1)
      perm_B_hat[,q] <- perm_sign_factor[q]*perm_B_hat[,q]

      perm_list_Zeta_hat[[q]] <- perm_list_Zeta_hat_untrimmed[[q]] <- perm_sign_factor[q]*list_Zeta_hat[[best_perm[q]]]
      perm_list_list_Phi_hat[[q]] <- perm_list_list_Phi_hat_untrimmed[[q]] <- list_list_Phi_hat[[best_perm[q]]]

      if (!is.null(list_h_hat)) {
        for (i in 1:N) {
          perm_list_h_hat[[i]][[q]] <- perm_list_h_hat_untrimmed[[i]][[q]] <- perm_sign_factor[q]*list_h_hat[[i]][[best_perm[q]]]
        }
      }
    } else {
      perm_list_Zeta_hat_untrimmed[[q]] <- list_Zeta_hat[[setdiff(1:Q, best_perm)[q-Q_true]]]
      perm_list_list_Phi_hat_untrimmed[[q]] <- list_list_Phi_hat[[setdiff(1:Q, best_perm)[q-Q_true]]]
      if (!is.null(list_h_hat)) {
        for (i in 1:N) {
          perm_list_h_hat_untrimmed[[i]][[q]] <- list_h_hat[[i]][[setdiff(1:Q, best_perm)[q-Q_true]]]
        }
      }
    }
  }

  perm_ppi_untrimmed <- cbind(perm_ppi, ppi[, -best_perm, drop = F])
  perm_B_hat_untrimmed <- cbind(perm_B_hat, B_hat[, -best_perm, drop = F])
  perm_factor_ppi <- c(factor_ppi[best_perm], factor_ppi[-best_perm])
  perm_list_cumulated_pve <- c(list_cumulated_pve[best_perm], list_cumulated_pve[-best_perm])
  perm_list_Cov_zeta_hat <- c(list_Cov_zeta_hat[best_perm], list_Cov_zeta_hat[-best_perm])

  if (!is.null(list_var_vec)) {
    perm_list_var_vec <- c(list_var_vec[best_perm], list_var_vec[-best_perm])
  } else {
    perm_list_var_vec <- NULL
  }

  create_named_list(best_perm, perm_sign_factor, perm_factor_ppi,
                    perm_B_hat_untrimmed, perm_ppi_untrimmed,
                    perm_B_hat, perm_ppi,
                    perm_list_Zeta_hat, perm_list_list_Phi_hat,
                    perm_list_Zeta_hat_untrimmed, perm_list_list_Phi_hat_untrimmed,
                    perm_list_cumulated_pve, perm_list_Cov_zeta_hat,
                    perm_list_h_hat,
                    perm_list_h_hat_untrimmed,
                    perm_list_var_vec)
}

match_sign_components <- function(Zeta,
                                  list_Zeta_hat,
                                  list_list_Phi_hat){
  Q_true <- length(Zeta)
  perm_list_Zeta_hat <- list_Zeta_hat
  perm_list_list_Phi_hat <- list_list_Phi_hat

  L_true_vec <- vapply(Zeta, ncol, integer(1))
  L_hat_vec <- vapply(list_Zeta_hat, ncol, integer(1))
  if (any(L_true_vec > L_hat_vec)) {
    stop("Each factor must have at least as many estimated as true FPCA components.")
  }
  max_L_true <- max(L_true_vec, 0L)
  perm_component <- matrix(NA_integer_, nrow = Q_true, ncol = max_L_true)
  perm_sign_fpca <- matrix(NA_real_, nrow = Q_true, ncol = max_L_true)

  for (q in seq_len(Q_true)) {
    L_true <- L_true_vec[q]
    L_hat <- L_hat_vec[q]
    if (L_true == 0L) next
    corr_component <- outer(seq_len(L_true), seq_len(L_hat),
      Vectorize(function(l, l_hat) {
        .safe_correlation(Zeta[[q]][, l], list_Zeta_hat[[q]][, l_hat])
      }))
    assignment <- .hungarian_assignment(1 - abs(corr_component))
    signs <- sign(corr_component[cbind(seq_len(L_true), assignment)])
    signs[!is.finite(signs) | signs == 0] <- 1

    perm_component[q, seq_len(L_true)] <- assignment
    perm_sign_fpca[q, seq_len(L_true)] <- signs
    perm_list_Zeta_hat[[q]] <- sweep(
      list_Zeta_hat[[q]][, assignment, drop = FALSE], 2, signs, "*")
    perm_list_list_Phi_hat[[q]] <- sweep(
      list_list_Phi_hat[[q]][, assignment, drop = FALSE], 2, signs, "*")
  }
  create_named_list(perm_sign_fpca, perm_component,
                    perm_list_Zeta_hat, perm_list_list_Phi_hat)
}

log_one_plus_exp_ <- function(x) {
  m <- x
  m[x < 0] <- 0
  log(exp(x - m) + exp(- m)) + m
}

log1pExp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

get_annealing_ladder_ <- function(anneal, verbose) {
  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if(anneal[1] == 1) {
    type <- "geometric"
    delta_k <- k_m^(1 / (1 - m)) - 1
    ladder <- (1 + delta_k)^(1 - m:1)
  } else if (anneal[1] == 2) {
    type <- "harmonic"
    delta_k <- ( 1 / k_m - 1) / (m - 1)
    ladder <- 1 / (1 + delta_k * (m:1 - 1))
  } else {
    type <- "linear"
    delta_k <- (1 - k_m) / (m - 1)
    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose != 0)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder
}

check_annealing <- function(anneal, verbose) {
  check_structure(anneal, "vector", "numeric", 3, null_ok = TRUE)

  if (!is.null(anneal)) {
    if (verbose) cat("== Checking the annealing schedule ... \n\n")
    check_natural(anneal[c(1, 3)])
    check_positive(anneal[2])

    if (!(anneal[1] %in% 1:3))
      stop(paste0("The annealing spacing scheme must be set to 1 for geometric ",
                  "2 for harmonic or 3 for linear spacing."))

    if (anneal[2] >= 2)
      stop(paste0("Initial annealing temperature must be strictly smaller than 2.\n ",
                  "Please decrease it."))

    if (anneal[2] <= 1)
      stop("Initial annealing temperature must be strictly greater than 1.")

    if (anneal[3] < 2)
      stop("Temperature grid size anneal[3] must be at least 2.")

    if (anneal[3] > 1000)
      stop(paste0("Temperature grid size very large. This may be unnecessarily ",
                  "computationally demanding. Please decrease it."))

    if (verbose) cat("... done. == \n\n")
  }
}

# ---- from bayesSYNC/R/OSullivan_splines.R ----

ZOSull <- function(x,range.x,intKnots,drv=0)
{
  if (!missing(range.x))
  {
    if (length(range.x)!=2) stop("range.x must be of length 2.")
    if (range.x[1]>range.x[2]) stop("range.x[1] exceeds range.x[2].")
    if (range.x[1]>min(x)) stop("range.x[1] must be <= than min(x).")
    if (range.x[2]<max(x)) stop("range.x[2] must be >= than max(x).")
  }

  if (drv>2) stop("splines not smooth enough for more than 2 derivatives")

  if (missing(range.x))
    range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))

  if (missing(intKnots))
  {
    numIntKnots <- min(length(unique(x)),35)
    intKnots <- quantile(unique(x),seq(0,1,length=
                                         (numIntKnots+2))[-c(1,(numIntKnots+2))])
  }
  numIntKnots <- length(intKnots)

  allKnots <- c(rep(range.x[1],4),intKnots,rep(range.x[2],4))
  K <- length(intKnots) ; L <- 3*(K+8)
  xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+
               rep(allKnots,each=3)[-c(1,2,L)])/2
  wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
  Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                    outer.ok=TRUE)$design
  Omega     <- crossprod(Bdd*wts,Bdd)

  svdOmega <- svd(Omega)
  indsZ <- 1:(numIntKnots+2)
  UZ <- svdOmega$u[,indsZ]
  LZ <- t(t(UZ)/sqrt(svdOmega$d[indsZ]))

  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- svdOmega$u[,indsX]
  L <- cbind(UX,LZ)
  stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))
  if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
    print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")

  B <- spline.des(allKnots,x,derivs=rep(drv,length(x)),
                  outer.ok=TRUE)$design

  Z <- B%*%LZ

  attr(Z,"range.x") <- range.x
  attr(Z,"intKnots") <- intKnots

  return(Z)
}

OmegaOSull <- function(a, b, intKnots) {
  allKnots <- c(rep(a,4),intKnots,rep(b,4))
  K <- length(intKnots) ; L <- 3*(K+8)
  xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+
               rep(allKnots,each=3)[-c(1,2,L)])/2
  wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
  Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                    outer.ok=TRUE)$design
  Omega     <- crossprod(Bdd*wts,Bdd)
  return(Omega)
}

# ---- from bayesSYNC/R/set_hyper.R ----

#' Gather model hyperparameters.
#'
#' @param sigma_beta Vector of size 1 or 2 for the standard deviation of the
#'                   spline coefficients for the mean function.
#' @param A Positive real number for the top-level hyperparameter.
#' @param c_0 Beta prior shape1.
#' @param d_0 Beta prior shape2. If NULL, bayesSYNC_multi() sets it to p.
#'
#' @return An object containing the hyperparameter settings.
#'
#' @noRd
#' @export
#'
set_hyper <- function(sigma_beta = 1e5, A = 1e5, c_0 = 1, d_0 = NULL) {

  check_structure(sigma_beta, "vector", "numeric", c(1, 2))
  check_positive(sigma_beta)

  if (length(sigma_beta) == 1) sigma_beta <- rep(sigma_beta, 2)
  Sigma_beta <- sigma_beta^2*diag(2)

  check_structure(A, "vector", "numeric", 1)
  check_positive(A)

  check_structure(c_0, "vector", "numeric", 1)
  check_positive(c_0)
  check_structure(d_0, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(d_0)) check_positive(d_0)

  sigma_zeta <- 1
  mu_beta <- rep(0, 2)

  list_hyper <- create_named_list(sigma_zeta, mu_beta, Sigma_beta, A, c_0, d_0)

  class(list_hyper) <- "hyper"

  list_hyper
}
