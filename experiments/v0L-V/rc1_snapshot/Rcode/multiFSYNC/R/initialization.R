# Data-driven initialization for the multi-study functional factor model.
#
# These routines only construct iteration-zero variational means. They do not
# change a prior, a CAVI update, the annealing ladder, or the ordinary ELBO.
# All calculations use the working-scale response supplied to the core fitter.

.default_initialization_control <- function() {
  list(
    grid_size = 81L,
    perturb_sd = 0.05,
    rank_tol = 1e-8
  )
}

.validate_initialization_control <- function(control = NULL) {
  defaults <- .default_initialization_control()
  if (is.null(control)) return(defaults)
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control)))) {
    stop("initialization_control must be NULL or a fully named list.")
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown initialization_control field(s): ",
      paste(unknown, collapse = ", ")
    )
  }
  resolved <- utils::modifyList(defaults, control)

  grid_size <- resolved$grid_size
  if (length(grid_size) != 1L || !is.finite(grid_size) ||
      !is_int(grid_size) || grid_size < 3L) {
    stop("initialization_control$grid_size must be one integer of at least 3.")
  }
  resolved$grid_size <- as.integer(grid_size)

  for (field in c("perturb_sd", "rank_tol")) {
    value <- resolved[[field]]
    if (length(value) != 1L || !is.finite(value) || value < 0) {
      stop(
        "initialization_control$", field,
        " must be one finite non-negative number."
      )
    }
  }
  if (resolved$rank_tol <= 0 || resolved$rank_tol >= 1) {
    stop("initialization_control$rank_tol must lie strictly between 0 and 1.")
  }
  resolved
}

.initialization_grid_indices <- function(n_g, K_total, grid_size) {
  target <- min(n_g, max(K_total, grid_size))
  index <- unique(as.integer(round(seq(1L, n_g, length.out = target))))
  if (length(index) < K_total) index <- seq_len(n_g)
  index
}

.initialization_interpolate <- function(time, response, grid) {
  ord <- order(time)
  time <- time[ord]
  response <- response[ord]
  if (length(unique(time)) < 2L) {
    return(rep(mean(response), length(grid)))
  }
  as.numeric(stats::approx(
    x = time,
    y = response,
    xout = grid,
    rule = 2,
    ties = mean
  )$y)
}

.initialization_weighted_sse <- function(residual, weights) {
  sum(sweep(residual^2, 2L, weights, "*"))
}

.initialization_canonical_sign <- function(loading) {
  pivot <- which.max(abs(loading))
  if (length(pivot) && loading[pivot] < 0) -loading else loading
}

.initialization_factor_failure <- function(label, reason, before_sse) {
  list(
    success = FALSE,
    label = label,
    reason = reason,
    residual = NULL,
    loading = NULL,
    time_coef = NULL,
    scores = NULL,
    fitted = NULL,
    before_sse = before_sse,
    after_sse = before_sse,
    explained_fraction = 0,
    variable_eigenvalue = NA_real_,
    signal_fraction = NA_real_,
    temporal_rank = 0L
  )
}

# Extract one separable variable direction, followed by temporal FPCA of its
# projected subject curves. The temporal functions are represented in the
# current O'Sullivan span and the computational scores have empirical RMS one.
.initialization_one_factor <- function(
    residual, C_grid, weights, M, control, label) {
  if (length(M) != 1L || !is.finite(M) ||
      !is_int(M) || M < 1L) {
    stop("M must be one positive integer in initialization_one_factor().")
  }
  M <- as.integer(M)
  dims <- dim(residual)
  n_subject <- dims[1L]
  n_grid <- dims[2L]
  p <- dims[3L]
  before_sse <- .initialization_weighted_sse(residual, weights)
  if (!is.finite(before_sse) || before_sse <= 0) {
    return(.initialization_factor_failure(
      label, "zero_or_nonfinite_residual_energy", before_sse))
  }

  variable_covariance <- matrix(0, p, p)
  for (i in seq_len(n_subject)) {
    residual_i <- matrix(
      residual[i, , ], nrow = n_grid, ncol = p)
    variable_covariance <- variable_covariance +
      crossprod(residual_i, weights * residual_i)
  }
  variable_covariance <-
    (variable_covariance + t(variable_covariance)) / 2
  eig <- eigen(variable_covariance, symmetric = TRUE)
  total_eigenvalue <- sum(pmax(eig$values, 0))
  leading_eigenvalue <- eig$values[1L]
  if (!is.finite(total_eigenvalue) || total_eigenvalue <= 0 ||
      !is.finite(leading_eigenvalue) ||
      leading_eigenvalue <=
        control$rank_tol * total_eigenvalue) {
    return(.initialization_factor_failure(
      label, "insufficient_variable_signal", before_sse))
  }

  loading <- eig$vectors[, 1L]
  if (control$perturb_sd > 0) {
    loading <- loading +
      control$perturb_sd * stats::rnorm(p) / sqrt(p)
  }
  loading_norm <- sqrt(sum(loading^2))
  if (!is.finite(loading_norm) ||
      loading_norm <= control$rank_tol) {
    return(.initialization_factor_failure(
      label, "nonfinite_loading_direction", before_sse))
  }
  loading <- .initialization_canonical_sign(loading / loading_norm)

  projected <- matrix(0, n_subject, n_grid)
  for (i in seq_len(n_subject)) {
    residual_i <- matrix(
      residual[i, , ], nrow = n_grid, ncol = p)
    projected[i, ] <- as.vector(residual_i %*% loading)
  }
  projected <- sweep(projected, 2L, colMeans(projected), "-")

  qr_grid <- qr(C_grid, tol = control$rank_tol)
  if (qr_grid$rank < ncol(C_grid)) {
    return(.initialization_factor_failure(
      label, "rank_deficient_spline_grid", before_sse))
  }
  projected_coef <- qr.coef(qr_grid, t(projected))
  if (any(!is.finite(projected_coef))) {
    return(.initialization_factor_failure(
      label, "nonfinite_projected_spline_coefficients", before_sse))
  }
  projected_smooth <- t(C_grid %*% projected_coef)
  projected_smooth <-
    sweep(projected_smooth, 2L, colMeans(projected_smooth), "-")

  weighted_projected <- sweep(
    projected_smooth, 2L, sqrt(weights), "*")
  max_rank <- min(dim(weighted_projected))
  if (max_rank < M) {
    return(.initialization_factor_failure(
      label, "insufficient_subject_time_rank", before_sse))
  }
  sv <- svd(
    weighted_projected,
    nu = min(M, n_subject),
    nv = min(M, n_grid)
  )
  if (!length(sv$d) || !is.finite(sv$d[1L]) || sv$d[1L] <= 0) {
    return(.initialization_factor_failure(
      label, "insufficient_temporal_signal", before_sse))
  }
  temporal_rank <- sum(
    sv$d > control$rank_tol * sv$d[1L])
  if (temporal_rank < M) {
    return(.initialization_factor_failure(
      label, "insufficient_temporal_signal", before_sse))
  }

  raw_scores <- sweep(
    sv$u[, seq_len(M), drop = FALSE],
    2L,
    sv$d[seq_len(M)],
    "*"
  )
  raw_scores <- sweep(raw_scores, 2L, colMeans(raw_scores), "-")
  score_scale <- sqrt(colMeans(raw_scores^2))
  if (any(!is.finite(score_scale)) ||
      any(score_scale <= 0)) {
    return(.initialization_factor_failure(
      label, "degenerate_score_scale", before_sse))
  }
  scores <- sweep(raw_scores, 2L, score_scale, "/")
  time_grid <- sweep(
    sweep(
      sv$v[, seq_len(M), drop = FALSE],
      1L,
      sqrt(weights),
      "/"
    ),
    2L,
    score_scale,
    "*"
  )
  time_coef <- qr.coef(qr_grid, time_grid)
  if (any(!is.finite(time_coef))) {
    return(.initialization_factor_failure(
      label, "nonfinite_time_coefficients", before_sse))
  }
  time_fitted <- C_grid %*% time_coef
  projected_fitted <- scores %*% t(time_fitted)

  fitted <- array(0, dim = dims)
  for (i in seq_len(n_subject)) {
    fitted[i, , ] <- tcrossprod(
      projected_fitted[i, ], loading)
  }
  residual_after <- residual - fitted
  after_sse <- .initialization_weighted_sse(
    residual_after, weights)
  if (!is.finite(after_sse) ||
      after_sse > before_sse * (1 + 1e-8)) {
    return(.initialization_factor_failure(
      label, "factor_reconstruction_increased_sse", before_sse))
  }

  perturbed_eigenvalue <- as.numeric(
    crossprod(loading, variable_covariance %*% loading))
  list(
    success = TRUE,
    label = label,
    reason = "",
    residual = residual_after,
    loading = loading,
    time_coef = time_coef,
    scores = scores,
    fitted = fitted,
    before_sse = before_sse,
    after_sse = after_sse,
    explained_fraction =
      (before_sse - after_sse) / before_sse,
    variable_eigenvalue = perturbed_eigenvalue,
    signal_fraction =
      perturbed_eigenvalue / total_eigenvalue,
    temporal_rank = as.integer(temporal_rank)
  )
}

.residual_fpca_initialization <- function(
    Y, Z, time_obs, C_g, time_g,
    S, n_s, p, d, L_f, L_s, M_f, M_s,
    control) {
  started <- proc.time()[["elapsed"]]
  K_total <- ncol(C_g)
  grid_index <- .initialization_grid_indices(
    length(time_g), K_total, control$grid_size)
  grid <- time_g[grid_index]
  C_grid <- C_g[grid_index, , drop = FALSE]
  weights <- .trap_weights(grid)
  n_subject <- sum(n_s)
  subject_study <- rep(seq_len(S), n_s)

  response_grid <- array(0, dim = c(n_subject, length(grid), p))
  row_now <- 0L
  for (s in seq_len(S)) {
    for (i in seq_len(n_s[s])) {
      row_now <- row_now + 1L
      for (j in seq_len(p)) {
        response_grid[row_now, , j] <-
          .initialization_interpolate(
            time_obs[[s]][[i]],
            Y[[s]][[i]][[j]],
            grid
          )
      }
    }
  }

  study_design <- matrix(0, n_subject, S)
  study_design[cbind(seq_len(n_subject), subject_study)] <- 1
  design <- study_design
  if (d > 0L) {
    design <- cbind(design, do.call(rbind, Z))
  }
  qr_design <- qr(design, tol = control$rank_tol)
  if (qr_design$rank < ncol(design)) {
    stop(
      "residual_fpca initialization failed: study/covariate design ",
      "is rank deficient."
    )
  }
  qr_grid <- qr(C_grid, tol = control$rank_tol)
  if (qr_grid$rank < K_total) {
    stop(
      "residual_fpca initialization failed: selected spline grid ",
      "is rank deficient."
    )
  }

  residual <- array(0, dim = dim(response_grid))
  mean_grid <- array(0, dim = c(S, length(grid), p))
  beta_grid <- if (d > 0L) {
    array(0, dim = c(p, d, length(grid)))
  } else {
    NULL
  }
  for (j in seq_len(p)) {
    response_j <- matrix(
      response_grid[, , j],
      nrow = n_subject,
      ncol = length(grid)
    )
    effect_j <- qr.coef(qr_design, response_j)
    if (any(!is.finite(effect_j))) {
      stop(
        "residual_fpca initialization failed: non-finite mean/covariate fit."
      )
    }
    residual[, , j] <- response_j - design %*% effect_j
    for (s in seq_len(S)) {
      mean_grid[s, , j] <- effect_j[s, ]
    }
    if (d > 0L) {
      for (r in seq_len(d)) {
        beta_grid[j, r, ] <- effect_j[S + r, ]
      }
    }
  }

  mu_q_nu_mu <- lapply(seq_len(S), function(s) {
    lapply(seq_len(p), function(j) {
      as.numeric(qr.coef(qr_grid, mean_grid[s, , j]))
    })
  })
  mu_q_nu_beta <- if (d > 0L) {
    lapply(seq_len(p), function(j) {
      lapply(seq_len(d), function(r) {
        as.numeric(qr.coef(qr_grid, beta_grid[j, r, ]))
      })
    })
  } else {
    NULL
  }

  residualized_sse <- .initialization_weighted_sse(
    residual, weights)
  shared <- vector("list", L_f)
  shared_loadings <- matrix(NA_real_, p, L_f)
  for (l in seq_len(L_f)) {
    factor_result <- .initialization_one_factor(
      residual = residual,
      C_grid = C_grid,
      weights = weights,
      M = M_f[l],
      control = control,
      label = sprintf("shared[%d]", l)
    )
    shared[[l]] <- factor_result
    if (factor_result$success) {
      residual <- factor_result$residual
      shared_loadings[, l] <- factor_result$loading
    }
  }

  has_specific <- .has_specific(L_s)
  specific <- if (has_specific) vector("list", S) else NULL
  specific_loadings <- if (has_specific) {
    lapply(seq_len(S), function(s) {
      matrix(NA_real_, p, .L_s_at(L_s, s))
    })
  } else {
    NULL
  }
  if (has_specific) {
    for (s in seq_len(S)) {
      L_ss <- .L_s_at(L_s, s)
      rows_s <- which(subject_study == s)
      residual_s <- residual[rows_s, , , drop = FALSE]
      specific[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        factor_result <- .initialization_one_factor(
          residual = residual_s,
          C_grid = C_grid,
          weights = weights,
          M = M_s[[s]][l],
          control = control,
          label = sprintf("specific[%d,%d]", s, l)
        )
        specific[[s]][[l]] <- factor_result
        if (factor_result$success) {
          residual_s <- factor_result$residual
          specific_loadings[[s]][, l] <- factor_result$loading
        }
      }
      residual[rows_s, , ] <- residual_s
    }
  }

  factor_results <- c(
    shared,
    if (has_specific) unlist(specific, recursive = FALSE) else list()
  )
  success <- if (length(factor_results)) {
    vapply(factor_results, `[[`, logical(1), "success")
  } else {
    logical()
  }
  factor_diagnostics <- if (length(factor_results)) {
    do.call(rbind, lapply(factor_results, function(x) {
      data.frame(
        block = x$label,
        success = x$success,
        fallback_reason = x$reason,
        before_sse = x$before_sse,
        after_sse = x$after_sse,
        explained_fraction = x$explained_fraction,
        variable_eigenvalue = x$variable_eigenvalue,
        signal_fraction = x$signal_fraction,
        temporal_rank = x$temporal_rank,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      block = character(),
      success = logical(),
      fallback_reason = character(),
      before_sse = double(),
      after_sse = double(),
      explained_fraction = double(),
      variable_eigenvalue = double(),
      signal_fraction = double(),
      temporal_rank = integer(),
      stringsAsFactors = FALSE
    )
  }
  remaining_sse <- .initialization_weighted_sse(residual, weights)

  list(
    mu_q_nu_mu = mu_q_nu_mu,
    mu_q_nu_beta = mu_q_nu_beta,
    shared = shared,
    specific = specific,
    diagnostics = list(
      requested_method = "residual_fpca",
      used_method = if (all(success)) {
        "residual_fpca"
      } else {
        "residual_fpca_partial"
      },
      complete = all(success),
      fallback_count = sum(!success),
      grid_size = length(grid),
      grid_index = grid_index,
      residualized_sse = residualized_sse,
      remaining_sse = remaining_sse,
      explained_fraction = if (residualized_sse > 0) {
        (residualized_sse - remaining_sse) / residualized_sse
      } else {
        0
      },
      shared_loadings = shared_loadings,
      specific_loadings = specific_loadings,
      shared_score_means = lapply(shared, function(x) {
        if (isTRUE(x$success)) colMeans(x$scores) else NULL
      }),
      shared_score_second_moments = lapply(shared, function(x) {
        if (isTRUE(x$success)) colMeans(x$scores^2) else NULL
      }),
      specific_score_means = if (has_specific) {
        lapply(specific, function(study) {
          lapply(study, function(x) {
            if (isTRUE(x$success)) colMeans(x$scores) else NULL
          })
        })
      } else {
        NULL
      },
      specific_score_second_moments = if (has_specific) {
        lapply(specific, function(study) {
          lapply(study, function(x) {
            if (isTRUE(x$success)) colMeans(x$scores^2) else NULL
          })
        })
      } else {
        NULL
      },
      factor_diagnostics = factor_diagnostics,
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
  )
}
