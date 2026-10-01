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

  # Historical structured truth: ordinary B-splines with equally spaced
  # knots. Each block is generated independently and is never projected to the
  # fitted O'Sullivan basis. Random coefficient mixtures avoid reusing the
  # same first M basis columns for every latent factor.
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

#' Construct Jaoua/bayesSYNC orthonormal spline truth blocks
#'
#' This reproduces the construction and common-column sampling convention in
#' `bayesSYNC::construct_splines()` and `generate_bayesSYNC_data()` as bound on
#' 2026-09-09. The upstream generator assumes a common FPCA truncation for all
#' latent factors and supports at most four factors and four components.
#'
#' @param Q Total number of shared and study-specific latent-factor blocks.
#' @param M Common number of FPCA components in every block.
#' @param n_g Dense grid size used only to store truth values.
#' @return A list of factor-specific truth blocks.
#' @keywords internal
construct_jaoua_splinet_blocks <- function(Q, M, n_g = 200) {
  if (!requireNamespace("Splinets", quietly = TRUE)) {
    stop(
      "Package 'Splinets' is required for truth_basis_method = ",
      "'jaoua_splinets'. Install version 1.5.1 or later."
    )
  }
  if (length(Q) != 1L || !is.finite(Q) || Q != as.integer(Q) ||
      Q < 1L || Q > 4L) {
    stop("Jaoua Splinets truth requires a total factor count Q between 1 and 4.")
  }
  if (length(M) != 1L || !is.finite(M) || M != as.integer(M) ||
      M < 2L || M > 4L) {
    stop("Jaoua Splinets truth requires a common component count M between 2 and 4.")
  }
  if (length(n_g) != 1L || !is.finite(n_g) || n_g != as.integer(n_g) ||
      n_g < max(50L, 10L * M)) {
    stop("n_g must be an integer of at least max(50, 10 * M).")
  }

  Q <- as.integer(Q)
  M <- as.integer(M)
  t_grid <- seq(0, 1, length.out = n_g)
  dt <- diff(t_grid)
  w_grid <- c(dt[1L] / 2,
              (dt[-1L] + dt[-length(dt)]) / 2,
              dt[length(dt)] / 2)

  n_int_knots <- M + 2L
  degrees <- c(3L, 2L, 3L, 2L)
  periodic <- c(FALSE, FALSE, TRUE, TRUE)
  knot_offsets <- c(3L, 2L, 4L, 3L)
  spline_systems <- vector("list", Q)
  knots_by_block <- vector("list", Q)

  for (q in seq_len(Q)) {
    knot_count <- n_int_knots + knot_offsets[q]
    knots_q <- seq(0, 1, length.out = knot_count)
    spline_systems[[q]] <- Splinets::splinet(
      knots_q,
      degree = degrees[q],
      periodic = periodic[q],
      norm = TRUE
    )
    knots_by_block[[q]] <- knots_q
  }

  # The upstream generator deliberately uses the first spline system to define
  # a common index pool, then samples an independent subset for every factor.
  first_eval <- Splinets::evspline(spline_systems[[1L]]$os)
  common_pool_size <- ncol(first_eval) - 1L
  if (common_pool_size < M) {
    stop("The Jaoua Splinets construction has fewer available functions than M.")
  }
  selected_indices <- lapply(
    seq_len(Q),
    function(q) sample(seq_len(common_pool_size), M, replace = FALSE)
  )

  source_binding <- list(
    bayesSYNC_repository = "https://github.com/hruffieux/bayesSYNC",
    bayesSYNC_master_commit = "de326142f15c8a087f544c84f18d83511aae50f1",
    bayesSYNC_simulator_commit = "816db2fe27bb6f596012135caaea3c8852eb8da8",
    bayesSYNC_simulator_sha256 =
      "58d981720c0a639a4f941fa4a1f1353255721258b261519fce671dcdd8517639",
    Splinets_version = as.character(utils::packageVersion("Splinets")),
    Splinets_upstream_commit = "a8f79565257da5d292acbd4fc03af62ecd3f899e"
  )

  blocks <- vector("list", Q)
  for (q in seq_len(Q)) {
    evaluated <- Splinets::evspline(spline_systems[[q]]$os, x = t_grid)
    phi <- evaluated[, -1L, drop = FALSE][, selected_indices[[q]], drop = FALSE]
    continuous_gram_all <- Splinets::gramian(spline_systems[[q]]$os)
    continuous_gram <- continuous_gram_all[
      selected_indices[[q]], selected_indices[[q]], drop = FALSE]
    if (max(abs(continuous_gram - diag(M))) > 1e-8) {
      stop(
        "Splinets did not return an orthonormal selected basis for block ", q,
        "; do not use this package/version combination for simulation."
      )
    }
    components <- lapply(seq_len(M), function(m) {
      list(
        degree = degrees[q],
        knots = knots_by_block[[q]],
        periodic = periodic[q],
        selected_index = selected_indices[[q]][m]
      )
    })
    blocks[[q]] <- list(
      phi = phi,
      t_grid = t_grid,
      w_grid = w_grid,
      gram = continuous_gram,
      grid_gram = crossprod(phi, w_grid * phi),
      block_id = q,
      components = components,
      method = "jaoua_splinets",
      degree = degrees[q],
      knots = knots_by_block[[q]],
      periodic = periodic[q],
      selected_indices = selected_indices[[q]],
      common_selection_pool_size = common_pool_size,
      source_binding = source_binding,
      evaluation = list(
        splinet_os = spline_systems[[q]]$os,
        selected_indices = selected_indices[[q]]
      )
    )
  }

  blocks
}

#' Evaluate a Jaoua/bayesSYNC Splinets truth block at arbitrary times
#' @keywords internal
eval_jaoua_splinet_block <- function(t_pts, evaluation) {
  evaluated <- Splinets::evspline(evaluation$splinet_os, x = t_pts)
  evaluated[, -1L, drop = FALSE][
    , evaluation$selected_indices, drop = FALSE]
}

#' Construct exchangeable factor blocks from one orthonormal spline system
#'
#' A common continuous-L2-orthonormal cubic Splinets system is constructed
#' once. Each latent-factor block then receives an independent Gaussian matrix.
#' An optional global second-derivative roughness weight is applied before its
#' thin QR factor is used as a Stiefel coefficient matrix. Splinets forms and
#' evaluates the resulting linear combinations exactly; no dense-grid
#' interpolation is used to generate observations.
#'
#' @param M_by_block Integer vector giving the FPCA size of each factor block.
#' @param n_g Dense grid size used only to store truth values.
#' @param basis_dim Dimension of the common orthonormal spline system.
#' @param smoothness_alpha Non-negative global second-derivative spectral
#'   weight. Zero reproduces the uniform Stiefel version-1 construction.
#' @return A list of independently generated factor-specific truth blocks.
#' @keywords internal
construct_iid_common_splinet_stiefel_blocks <- function(
    M_by_block, n_g = 200, basis_dim = 8L, smoothness_alpha = 0) {
  if (!requireNamespace("Splinets", quietly = TRUE)) {
    stop(
      "Package 'Splinets' is required for truth_basis_method = ",
      "'iid_common_splinet_stiefel'. Install version 1.5.1 or later."
    )
  }
  if (length(M_by_block) < 1L || any(!is.finite(M_by_block)) ||
      any(M_by_block != as.integer(M_by_block)) || any(M_by_block < 1L)) {
    stop("M_by_block must contain positive integers.")
  }
  M_by_block <- as.integer(M_by_block)
  if (length(basis_dim) != 1L || !is.finite(basis_dim) ||
      basis_dim != as.integer(basis_dim) || basis_dim < max(M_by_block)) {
    stop("truth_basis_dim must be one integer at least max(M_f, M_s).")
  }
  basis_dim <- as.integer(basis_dim)
  if (length(n_g) != 1L || !is.finite(n_g) || n_g != as.integer(n_g) ||
      n_g < max(50L, 10L * max(M_by_block))) {
    stop("n_g must be an integer of at least max(50, 10 * max(M_by_block)).")
  }
  if (length(smoothness_alpha) != 1L ||
      !is.finite(smoothness_alpha) || smoothness_alpha < 0) {
    stop("smoothness_alpha must be one non-negative finite number.")
  }
  smoothness_alpha <- as.numeric(smoothness_alpha)

  degree <- 3L
  periodic <- FALSE
  knots <- seq(0, 1, length.out = basis_dim + degree + 1L)
  common_system <- Splinets::splinet(
    knots, degree = degree, periodic = periodic, norm = TRUE)
  common_gram <- Splinets::gramian(common_system$os)
  if (!identical(dim(common_gram), c(basis_dim, basis_dim)) ||
      max(abs(common_gram - diag(basis_dim))) > 1e-8) {
    stop(
      "Splinets did not return the requested continuous-L2-orthonormal ",
      "common spline system."
    )
  }

  second_derivative <- Splinets::deriva(
    Splinets::deriva(common_system$os))
  second_derivative_gram <- Splinets::gramian(second_derivative)
  second_derivative_gram <-
    (second_derivative_gram + t(second_derivative_gram)) / 2
  roughness_eigen <- eigen(second_derivative_gram, symmetric = TRUE)
  roughness_eigenvalues <- pmax(roughness_eigen$values, 0)
  roughness_scale <- stats::median(roughness_eigenvalues)
  if (!is.finite(roughness_scale) || roughness_scale <= 0) {
    stop("Invalid Splinets second-derivative Gram spectrum.")
  }
  spectral_weights <- (1 + smoothness_alpha *
    roughness_eigenvalues / roughness_scale)^(-0.5)
  coefficient_map <- if (smoothness_alpha == 0) {
    diag(basis_dim)
  } else {
    sweep(
      roughness_eigen$vectors, 2L, spectral_weights, FUN = "*") %*%
      t(roughness_eigen$vectors)
  }

  t_grid <- seq(0, 1, length.out = n_g)
  dt <- diff(t_grid)
  w_grid <- c(dt[1L] / 2,
              (dt[-1L] + dt[-length(dt)]) / 2,
              dt[length(dt)] / 2)
  source_binding <- list(
    Splinets_version = as.character(utils::packageVersion("Splinets")),
    construction = paste(
      "Splinets::splinet + Splinets::deriva + Splinets::gramian +",
      "Splinets::lincomb + base::eigen + base::qr")
  )

  blocks <- vector("list", length(M_by_block))
  for (q in seq_along(M_by_block)) {
    M_q <- M_by_block[[q]]
    block <- NULL
    for (attempt in seq_len(10L)) {
      gaussian_q <- matrix(
        stats::rnorm(basis_dim * M_q), basis_dim, M_q)
      z_q <- coefficient_map %*% gaussian_q
      qr_q <- qr(z_q)
      if (qr_q$rank < M_q) next
      Q_q <- qr.Q(qr_q)[, seq_len(M_q), drop = FALSE]
      qr_sign <- sign(diag(qr.R(qr_q))[seq_len(M_q)])
      qr_sign[qr_sign == 0] <- 1
      Q_q <- sweep(Q_q, 2L, qr_sign, FUN = "*")

      combined <- Splinets::lincomb(
        common_system$os, t(Q_q), reduced = TRUE)
      continuous_gram <- Splinets::gramian(combined)
      if (max(abs(continuous_gram - diag(M_q))) > 1e-8) next

      evaluated <- Splinets::evspline(combined, x = t_grid)
      phi <- evaluated[, -1L, drop = FALSE]
      components <- lapply(seq_len(M_q), function(m) {
        list(
          degree = degree,
          knots = knots,
          periodic = periodic,
          stiefel_coefficients = Q_q[, m]
        )
      })
      block <- list(
        phi = phi,
        t_grid = t_grid,
        w_grid = w_grid,
        gram = continuous_gram,
        grid_gram = crossprod(phi, w_grid * phi),
        block_id = as.integer(q),
        components = components,
        method = "iid_common_splinet_stiefel",
        degree = degree,
        knots = knots,
        periodic = periodic,
        common_basis_dim = basis_dim,
        common_basis_gram = common_gram,
        smoothness_alpha = smoothness_alpha,
        second_derivative_gram = second_derivative_gram,
        roughness_eigenvalues = roughness_eigenvalues,
        roughness_scale = roughness_scale,
        spectral_weights = spectral_weights,
        stiefel_coefficients = Q_q,
        source_binding = source_binding,
        evaluation = list(splinet_object = combined)
      )
      break
    }
    if (is.null(block)) {
      stop("Unable to construct a full-rank Stiefel spline truth block.")
    }
    blocks[[q]] <- block
  }

  blocks
}

#' Evaluate an iid-common-Splinets Stiefel truth block at arbitrary times
#' @keywords internal
eval_iid_common_splinet_stiefel_block <- function(t_pts, evaluation) {
  evaluated <- Splinets::evspline(evaluation$splinet_object, x = t_pts)
  evaluated[, -1L, drop = FALSE]
}

#' Evaluate code under an isolated RNG seed and restore the caller RNG stream
#' @keywords internal
.with_preserved_rng_seed <- function(seed, code,
                                     seed_name = "truth_basis_seed") {
  if (length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed) || seed < 0L ||
      seed > .Machine$integer.max) {
    stop(seed_name, " must be one non-negative integer.")
  }
  if (!is.function(code)) {
    stop("code must be a function.")
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    caller_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", caller_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  code()
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

#' Construct independent periodic B-spline mean truths
#'
#' The construction uses one continuous-L2-orthonormal periodic cubic
#' Splinets system. Study-variable coefficient vectors are generated
#' independently inside the exact zero-integral coefficient subspace and are
#' heat-kernel smoothed using the continuous second-derivative Gram matrix.
#' Every function is then scaled to the continuous L2 energy of the historical
#' periodic mean, `amplitude^2 * 5 / 8`.
#'
#' @param S Number of studies.
#' @param p Number of variables.
#' @param n_g Dense grid size used only to store truth values.
#' @param basis_dim Dimension of the periodic orthonormal spline system.
#' @param smoothness_alpha Fixed heat-kernel smoothing parameter.
#' @param amplitude Historical mean amplitude whose L2 energy is matched.
#' @return A list containing dense truth, direct-evaluation objects and
#'   generation metadata.
#' @keywords internal
construct_independent_splinet_means <- function(
    S, p, n_g = 200L, basis_dim = 8L,
    smoothness_alpha = 1.519826404906808, amplitude = 1) {
  if (!requireNamespace("Splinets", quietly = TRUE)) {
    stop(
      "Package 'Splinets' is required for mean_truth_method = ",
      "'independent_splinet'. Install version 1.5.1 or later."
    )
  }
  if (length(S) != 1L || !is.finite(S) || S != as.integer(S) || S < 1L ||
      length(p) != 1L || !is.finite(p) || p != as.integer(p) || p < 1L) {
    stop("S and p must be positive integers.")
  }
  if (length(n_g) != 1L || !is.finite(n_g) || n_g != as.integer(n_g) ||
      n_g < 50L) {
    stop("n_g must be an integer of at least 50.")
  }
  if (length(basis_dim) != 1L || !is.finite(basis_dim) ||
      basis_dim != as.integer(basis_dim) || basis_dim < 7L) {
    stop("mean_basis_dim must be one integer of at least 7.")
  }
  if (length(smoothness_alpha) != 1L ||
      !is.finite(smoothness_alpha) || smoothness_alpha < 0) {
    stop("mean_smoothness_alpha must be one non-negative finite number.")
  }
  if (length(amplitude) != 1L || !is.finite(amplitude) || amplitude <= 0) {
    stop("mean_amp must be one positive finite number for B-spline means.")
  }

  S <- as.integer(S)
  p <- as.integer(p)
  n_g <- as.integer(n_g)
  basis_dim <- as.integer(basis_dim)
  smoothness_alpha <- as.numeric(smoothness_alpha)
  amplitude <- as.numeric(amplitude)
  degree <- 3L
  periodic <- TRUE
  knots <- seq(0, 1, length.out = basis_dim + 1L)
  common_system <- Splinets::splinet(
    knots, degree = degree, periodic = periodic, norm = TRUE)
  common_gram <- Splinets::gramian(common_system$os)
  if (!identical(dim(common_gram), c(basis_dim, basis_dim)) ||
      max(abs(common_gram - diag(basis_dim))) > 1e-8) {
    stop(
      "Splinets did not return the requested continuous-L2-orthonormal ",
      "periodic mean system."
    )
  }

  integrated <- Splinets::integra(common_system$os)
  boundary <- Splinets::evspline(
    integrated, x = c(0, 1))[, -1L, drop = FALSE]
  integral_vector <- as.numeric(boundary[2L, ] - boundary[1L, ])
  if (any(!is.finite(integral_vector)) || sum(integral_vector^2) <= 0) {
    stop("Unable to compute the continuous integral of the mean basis.")
  }
  null_basis <- qr.Q(
    qr(matrix(integral_vector, ncol = 1L)), complete = TRUE)[
      , -1L, drop = FALSE]

  second_derivative <- Splinets::deriva(
    Splinets::deriva(common_system$os))
  second_derivative_gram <- Splinets::gramian(second_derivative)
  second_derivative_gram <-
    (second_derivative_gram + t(second_derivative_gram)) / 2
  restricted_penalty <- crossprod(
    null_basis, second_derivative_gram %*% null_basis)
  restricted_penalty <- (restricted_penalty + t(restricted_penalty)) / 2
  roughness_eigen <- eigen(restricted_penalty, symmetric = TRUE)
  roughness_eigenvalues <- pmax(roughness_eigen$values, 0)
  roughness_scale <- stats::median(roughness_eigenvalues)
  if (!is.finite(roughness_scale) || roughness_scale <= 0) {
    stop("Invalid restricted mean-roughness spectrum.")
  }
  spectral_weights <- exp(
    -0.5 * smoothness_alpha * roughness_eigenvalues / roughness_scale)
  coefficient_map <- null_basis %*%
    sweep(roughness_eigen$vectors, 2L, spectral_weights, FUN = "*") %*%
    t(roughness_eigen$vectors)

  target_energy <- amplitude^2 * 5 / 8
  target_normalized_roughness <- 64 * pi^4
  t_grid <- seq(0, 1, length.out = n_g)
  dense <- vector("list", S)
  evaluation <- vector("list", S)
  basis_meta <- vector("list", S)
  for (s in seq_len(S)) {
    raw_coefficients <- coefficient_map %*% matrix(
      stats::rnorm((basis_dim - 1L) * p), basis_dim - 1L, p)
    raw_energy <- colSums(raw_coefficients^2)
    if (any(!is.finite(raw_energy)) || any(raw_energy <= 1e-14)) {
      stop("Degenerate independent-splinet mean coefficient draw.")
    }
    coefficients <- sweep(
      raw_coefficients, 2L, sqrt(target_energy / raw_energy), FUN = "*")
    combined <- Splinets::lincomb(
      common_system$os, t(coefficients), reduced = TRUE)
    values <- Splinets::evspline(combined, x = t_grid)[
      , -1L, drop = FALSE]
    if (!identical(dim(values), c(n_g, p))) {
      stop("Unexpected direct-evaluation dimensions for B-spline means.")
    }
    continuous_energy <- colSums(coefficients^2)
    continuous_integral <- as.vector(crossprod(
      integral_vector, coefficients))
    normalized_roughness <- colSums(
      coefficients * (second_derivative_gram %*% coefficients)) /
      continuous_energy

    dense[[s]] <- values
    evaluation[[s]] <- list(splinet_object = combined)
    basis_meta[[s]] <- list(
      study = as.integer(s), basis_coefficients = coefficients,
      continuous_energy = continuous_energy,
      continuous_integral = continuous_integral,
      normalized_roughness = normalized_roughness)
  }

  list(
    dense = dense,
    evaluation = evaluation,
    basis_meta = basis_meta,
    common_meta = list(
      method = "independent_splinet",
      degree = degree,
      periodic = periodic,
      knots = knots,
      basis_dim = basis_dim,
      common_basis_gram = common_gram,
      integral_vector = integral_vector,
      null_basis = null_basis,
      second_derivative_gram = second_derivative_gram,
      restricted_roughness_eigenvalues = roughness_eigenvalues,
      roughness_scale = roughness_scale,
      smoothness_alpha = smoothness_alpha,
      spectral_weights = spectral_weights,
      target_energy = target_energy,
      target_normalized_roughness = target_normalized_roughness,
      source_binding = list(
        Splinets_version = as.character(utils::packageVersion("Splinets")),
        construction = paste(
          "Splinets::splinet(periodic=TRUE) + Splinets::integra +",
          "Splinets::deriva + Splinets::gramian + Splinets::lincomb +",
          "base::eigen + base::qr"))
    ),
    t_grid = t_grid)
}

#' Evaluate one study's independent-Splinets means at arbitrary times
#' @keywords internal
eval_independent_splinet_means <- function(t_pts, evaluation) {
  Splinets::evspline(evaluation$splinet_object, x = t_pts)[
    , -1L, drop = FALSE]
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
#' @param sparsity_mode `"fixed"` preserves the legacy fixed active proportion;
#'   `"factor_beta"` preserves the historical factor-specific Beta support with
#'   nonzero values concentrated around a factor-level random mean;
#'   `"jaoua_beta_iid_normal"` uses the Jaoua-style Beta/Bernoulli support law
#'   and independent standard-normal loading values. The latter is extended
#'   blockwise to shared and study-specific factors and requires
#'   `identified_loadings = FALSE` and `bool_sparse_loadings = TRUE`.
#' @param omega_beta Beta shape parameters used by `sparsity_mode =
#'   "factor_beta"`.
#' @param truth_basis_method Method used to construct orthonormal B-spline
#'   truth. `"legacy_dense_gram"` preserves all historical simulations.
#'   `"jaoua_splinets"` reproduces the current upstream bayesSYNC/Splinets
#'   construction and requires 1--4 total factors, a common M in 2--4, and the
#'   `Splinets` package. `"iid_common_splinet_stiefel"` constructs one common
#'   cubic non-periodic continuous-L2-orthonormal Splinets system and generates
#'   every shared or study-specific factor block independently from the same
#'   second-derivative-weighted Stiefel law. The route-B development setting
#'   fixes the global roughness weight at 3. It supports arbitrary factor
#'   counts and factor-specific M.
#' @param truth_basis_dim Dimension of the common orthonormal Splinets system
#'   used by `truth_basis_method = "iid_common_splinet_stiefel"`. It must be at
#'   least the largest factor-specific FPCA truncation. The default 8 is the
#'   version-1 development setting and is ignored by the historical methods.
#' @param truth_basis_seed Optional RNG seed used only for truth-basis
#'   construction. When supplied, the caller RNG stream is restored afterward,
#'   allowing paired generator arms to share observation times, scores,
#'   loadings, mean phases and noise while differing only in the truth basis.
#'   `NULL` preserves the historical RNG sequence.
#' @param truth_loading_seed Optional RNG seed used only for loading generation.
#'   When supplied, the caller RNG stream is restored afterward, allowing
#'   loading-generator sensitivity arms to share observation times, truth
#'   functions, scores, covariates, mean phases and measurement errors while
#'   differing in loadings and their resulting signals. `NULL` preserves the
#'   historical RNG sequence.
#' @param mean_truth_method Method used to generate mean functions.
#'   `"periodic_phase"` preserves the historical analytic periodic means.
#'   `"independent_splinet"` draws independent study-variable functions in a
#'   common continuous-L2-orthonormal periodic cubic Splinets system, centers
#'   them to zero integral and matches the historical continuous L2 energy.
#' @param mean_basis_dim Dimension of the orthonormal periodic Splinets system
#'   used by `mean_truth_method = "independent_splinet"`.
#' @param mean_smoothness_alpha Fixed heat-kernel spectral smoothing parameter
#'   for independent B-spline means. The default was calibrated using only
#'   generator geometry to match the historical normalized roughness.
#' @param mean_truth_seed Optional RNG seed used only for independent B-spline
#'   means. When supplied, the caller RNG stream is restored afterward.
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
    sparsity_mode = c("fixed", "factor_beta", "jaoua_beta_iid_normal"),
    omega_beta = c(1, 10),
    truth_basis_method = c(
      "legacy_dense_gram", "jaoua_splinets",
      "iid_common_splinet_stiefel"),
    truth_basis_dim = 8L,
    truth_basis_seed = NULL,
    truth_loading_seed = NULL,
    mean_truth_method = c("periodic_phase", "independent_splinet"),
    mean_basis_dim = 8L,
    mean_smoothness_alpha = 1.519826404906808,
    mean_truth_seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)
  sparsity_mode <- match.arg(sparsity_mode)
  truth_basis_method <- match.arg(truth_basis_method)
  mean_truth_method <- match.arg(mean_truth_method)
  if (mean_truth_method == "periodic_phase" && !is.null(mean_truth_seed)) {
    stop(
      "mean_truth_seed is only used when mean_truth_method = ",
      "'independent_splinet'."
    )
  }

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

  all_M <- c(as.integer(M_f),
             as.integer(unlist(M_s, recursive = FALSE, use.names = FALSE)))
  total_factors <- length(all_M)
  jaoua_blocks <- NULL
  stiefel_blocks <- NULL
  paired_basis_blocks <- NULL
  if (truth_basis_method == "jaoua_splinets") {
    if (total_factors > 0L && length(unique(all_M)) != 1L) {
      stop(
        "Jaoua Splinets truth requires the same FPCA component count M ",
        "for every shared and study-specific factor."
      )
    }
    if (total_factors > 0L && is.null(truth_basis_seed)) {
      jaoua_blocks <- construct_jaoua_splinet_blocks(
        total_factors, all_M[1L], n_g = n_dense)
    }
  }
  if (truth_basis_method == "iid_common_splinet_stiefel" &&
      total_factors > 0L && is.null(truth_basis_seed)) {
    stiefel_blocks <- construct_iid_common_splinet_stiefel_blocks(
      all_M, n_g = n_dense, basis_dim = truth_basis_dim,
      smoothness_alpha = 3)
  }
  if (!is.null(truth_basis_seed)) {
    paired_basis_blocks <- .with_preserved_rng_seed(
      truth_basis_seed,
      function() {
        if (total_factors < 1L) return(list())
        if (truth_basis_method == "jaoua_splinets") {
          construct_jaoua_splinet_blocks(
            total_factors, all_M[1L], n_g = n_dense)
        } else if (truth_basis_method == "iid_common_splinet_stiefel") {
          construct_iid_common_splinet_stiefel_blocks(
            all_M, n_g = n_dense, basis_dim = truth_basis_dim,
            smoothness_alpha = 3)
        } else {
          lapply(seq_len(total_factors), function(block_id) {
            construct_orthonormal_eigenfunctions(
              all_M[[block_id]], n_g = n_dense, degree = bs_degree,
              block_id = block_id)
          })
        }
      }
    )
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
  if (sparsity_mode == "jaoua_beta_iid_normal" &&
      (isTRUE(identified_loadings) || !isTRUE(bool_sparse_loadings))) {
    stop(
      "sparsity_mode = \"jaoua_beta_iid_normal\" requires ",
      "identified_loadings = FALSE and bool_sparse_loadings = TRUE."
    )
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
  phi_evaluation <- vector("list", L_f)
  block_id <- 0L
  for (l in seq_len(L_f)) {
    M_l <- M_f[l]
    block_id <- block_id + 1L
    ef_obj <- if (!is.null(paired_basis_blocks)) {
      paired_basis_blocks[[block_id]]
    } else if (truth_basis_method == "jaoua_splinets") {
      jaoua_blocks[[block_id]]
    } else if (truth_basis_method == "iid_common_splinet_stiefel") {
      stiefel_blocks[[block_id]]
    } else {
      construct_orthonormal_eigenfunctions(
        M_l, n_g = n_dense, degree = bs_degree, block_id = block_id)
    }
    lambda_l <- make_scheme1_eigenvalues(M_l, score_var_decay)
    phi_dense_list[[l]] <- ef_obj$phi
    theta_dense_list[[l]] <- sweep(
      ef_obj$phi, 2L, sqrt(lambda_l), FUN = "*")
    lambda_phi_true[[l]] <- lambda_l
    phi_basis_meta[[l]] <- ef_obj[intersect(
      c("block_id", "components", "gram", "grid_gram", "method", "degree",
        "knots", "periodic", "selected_indices",
        "common_selection_pool_size", "common_basis_dim",
        "common_basis_gram", "smoothness_alpha",
        "second_derivative_gram", "roughness_eigenvalues",
        "roughness_scale", "spectral_weights", "stiefel_coefficients",
        "source_binding"),
      names(ef_obj))]
    phi_evaluation[[l]] <- ef_obj$evaluation
  }

  # 3b: Construct separate study-specific truth blocks.  They are not reused
  # across studies or factors, preventing artificial shared/specific aliases.
  psi_dense_list <- NULL
  kappa_dense_list <- NULL
  lambda_psi_true <- NULL
  psi_basis_meta <- NULL
  psi_evaluation <- NULL
  if (has_specific) {
    psi_dense_list <- vector("list", S)
    kappa_dense_list <- vector("list", S)
    lambda_psi_true <- vector("list", S)
    psi_basis_meta <- vector("list", S)
    psi_evaluation <- vector("list", S)
    for (s in 1:S) {
      L_ss <- L_s[[s]]
      psi_dense_list[[s]] <- vector("list", L_ss)
      kappa_dense_list[[s]] <- vector("list", L_ss)
      lambda_psi_true[[s]] <- vector("list", L_ss)
      psi_basis_meta[[s]] <- vector("list", L_ss)
      psi_evaluation[[s]] <- vector("list", L_ss)
      for (l in seq_len(L_ss)) {
        M_sl <- M_s[[s]][l]
        block_id <- block_id + 1L
        ef_obj <- if (!is.null(paired_basis_blocks)) {
          paired_basis_blocks[[block_id]]
        } else if (truth_basis_method == "jaoua_splinets") {
          jaoua_blocks[[block_id]]
        } else if (truth_basis_method == "iid_common_splinet_stiefel") {
          stiefel_blocks[[block_id]]
        } else {
          construct_orthonormal_eigenfunctions(
            M_sl, n_g = n_dense, degree = bs_degree, block_id = block_id)
        }
        lambda_sl <- make_scheme1_eigenvalues(M_sl, score_var_decay)
        psi_dense_list[[s]][[l]] <- ef_obj$phi
        kappa_dense_list[[s]][[l]] <- sweep(
          ef_obj$phi, 2L, sqrt(lambda_sl), FUN = "*")
        lambda_psi_true[[s]][[l]] <- lambda_sl
        psi_basis_meta[[s]][[l]] <- ef_obj[intersect(
          c("block_id", "components", "gram", "grid_gram", "method", "degree",
            "knots", "periodic", "selected_indices",
            "common_selection_pool_size", "common_basis_dim",
            "common_basis_gram", "smoothness_alpha",
            "second_derivative_gram", "roughness_eigenvalues",
            "roughness_scale", "spectral_weights", "stiefel_coefficients",
            "source_binding"),
          names(ef_obj))]
        psi_evaluation[[s]][[l]] <- ef_obj$evaluation
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

  mean_truth_object <- NULL
  if (mean_truth_method == "independent_splinet") {
    draw_means <- function() construct_independent_splinet_means(
      S = S, p = p, n_g = n_dense, basis_dim = mean_basis_dim,
      smoothness_alpha = mean_smoothness_alpha, amplitude = mean_amp)
    mean_truth_object <- if (is.null(mean_truth_seed)) {
      draw_means()
    } else {
      .with_preserved_rng_seed(
        mean_truth_seed, draw_means, seed_name = "mean_truth_seed")
    }
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

  # 3e: Generate loadings with spike-and-slab. The optional loading seed keeps
  # all subsequent generator draws paired across loading-sensitivity arms.
  draw_all_loadings <- function() {
    loading_draw <- if (identified_loadings) {
      generate_identified_loadings(
        p, L_f, L_s, S,
        sparse = bool_sparse_loadings, prop_sparse = prop_sparse)
    } else NULL
    a_true <- if (identified_loadings) {
      loading_draw$a
    } else {
      matrix(0, p, L_f)
    }
    gamma_a_true <- if (identified_loadings) {
      loading_draw$gamma_a
    } else {
      matrix(0, p, L_f)
    }
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

    draw_unidentified_column <- function() {
      omega <- if (!bool_sparse_loadings) {
        1
      } else if (sparsity_mode == "fixed") {
        1 - prop_sparse
      } else {
        stats::rbeta(1L, omega_beta[1L], omega_beta[2L])
      }
      if (sparsity_mode == "jaoua_beta_iid_normal") {
        gamma <- stats::rbinom(p, size = 1L, prob = omega)
        if (!any(gamma == 1L)) gamma[sample(seq_len(p), 1L)] <- 1L
        value <- stats::rnorm(p, mean = 0, sd = 1) * gamma
        return(list(value = value, gamma = gamma, omega = omega))
      }

      active_idx <- draw_active_indices(omega)
      value <- numeric(p)
      gamma <- numeric(p)
      if (bool_sparse_loadings) {
        val_mean <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
        value[active_idx] <- rnorm(
          length(active_idx), mean = val_mean, sd = 0.2)
        gamma[active_idx] <- 1
      } else {
        value <- rnorm(p, mean = 0, sd = 0.5)
        gamma[] <- 1
      }
      list(value = value, gamma = gamma, omega = omega)
    }

    for (l in if (identified_loadings) integer(0) else seq_len(L_f)) {
      one <- draw_unidentified_column()
      a_true[, l] <- one$value
      gamma_a_true[, l] <- one$gamma
      omega_a_true[l] <- one$omega
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
      for (s in seq_len(S)) {
        L_ss <- L_s[[s]]
        b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
        gamma_b_true[[s]] <- matrix(0, nrow = p, ncol = L_ss)
        omega_b_true[[s]] <- rep(NA_real_, L_ss)
        for (l in seq_len(L_ss)) {
          one <- draw_unidentified_column()
          b_true[[s]][, l] <- one$value
          gamma_b_true[[s]][, l] <- one$gamma
          omega_b_true[[s]][l] <- one$omega
        }
      }
    }
    if (has_specific && identified_loadings) {
      for (s in seq_len(S)) {
        omega_b_true[[s]] <- colMeans(gamma_b_true[[s]])
      }
    }

    list(
      a_true = a_true,
      gamma_a_true = gamma_a_true,
      omega_a_true = omega_a_true,
      b_true = b_true,
      gamma_b_true = gamma_b_true,
      omega_b_true = omega_b_true
    )
  }
  loading_result <- if (is.null(truth_loading_seed)) {
    draw_all_loadings()
  } else {
    .with_preserved_rng_seed(
      truth_loading_seed, draw_all_loadings,
      seed_name = "truth_loading_seed")
  }
  a_true <- loading_result$a_true
  gamma_a_true <- loading_result$gamma_a_true
  omega_a_true <- loading_result$omega_a_true
  b_true <- loading_result$b_true
  gamma_b_true <- loading_result$gamma_b_true
  omega_b_true <- loading_result$omega_b_true

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

      # Mean functions. The phase draw remains part of both arms so an isolated
      # mean seed leaves all downstream generator randomness exactly paired.
      if (mean_truth_method == "independent_splinet") {
        mu_si <- eval_independent_splinet_means(
          t_i, mean_truth_object$evaluation[[s]])
      } else {
        for (j in 1:p) {
          mu_si[, j] <- mean_function_periodic(
            t_i, phase_mu_mat[s, j], mean_amp)
        }
      }

      # Shared factor processes
      for (l in seq_len(L_f)) {
        M_l <- M_f[l]
        eta_i <- zeta_true[[s]][[l]][i, ]
        if (truth_basis_method == "jaoua_splinets") {
          phi_val <- eval_jaoua_splinet_block(t_i, phi_evaluation[[l]])
          theta_val <- sweep(
            phi_val, 2L, sqrt(lambda_phi_true[[l]]), FUN = "*")
        } else if (truth_basis_method == "iid_common_splinet_stiefel") {
          phi_val <- eval_iid_common_splinet_stiefel_block(
            t_i, phi_evaluation[[l]])
          theta_val <- sweep(
            phi_val, 2L, sqrt(lambda_phi_true[[l]]), FUN = "*")
        } else {
          theta_val <- matrix(0, n_i, M_l)
          for (m in seq_len(M_l)) {
            theta_val[, m] <- eval_eigenfun(
              t_i, theta_dense_list[[l]][, m], t_grid_dense)
          }
        }
        f_si[, l] <- as.vector(theta_val %*% eta_i)
      }

      # Specific factor processes
      if (L_ss > 0L) {
        for (l in seq_len(L_ss)) {
          M_sl <- M_s[[s]][l]
          chi_i <- xi_true[[s]][[l]][i, ]
          if (truth_basis_method == "jaoua_splinets") {
            psi_val <- eval_jaoua_splinet_block(
              t_i, psi_evaluation[[s]][[l]])
            kappa_val <- sweep(
              psi_val, 2L, sqrt(lambda_psi_true[[s]][[l]]), FUN = "*")
          } else if (truth_basis_method == "iid_common_splinet_stiefel") {
            psi_val <- eval_iid_common_splinet_stiefel_block(
              t_i, psi_evaluation[[s]][[l]])
            kappa_val <- sweep(
              psi_val, 2L, sqrt(lambda_psi_true[[s]][[l]]), FUN = "*")
          } else {
            kappa_val <- matrix(0, n_i, M_sl)
            for (m in seq_len(M_sl)) {
              kappa_val[, m] <- eval_eigenfun(
                t_i, kappa_dense_list[[s]][[l]][, m], t_grid_dense)
            }
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
  truth_basis_family <- switch(
    truth_basis_method,
    jaoua_splinets = "jaoua_splinets_orthonormal_B_spline",
    iid_common_splinet_stiefel =
      "iid_common_splinet_stiefel_orthonormal_B_spline",
    legacy_dense_gram = "ordinary_B_spline_external")
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
    truth_basis_method,
    sparsity_mode,
    omega_beta,
    int_knots,
    identified_loadings,
    truth_basis_family,
    truth_projected_to_fitted_basis
  )
  true_params$L_s_by_study <- L_s
  true_params$L_s <- .compact_L_s(L_s)
  if (!is.null(truth_basis_seed)) {
    true_params$truth_basis_seed <- as.integer(truth_basis_seed)
  }
  if (!is.null(truth_loading_seed)) {
    true_params$truth_loading_seed <- as.integer(truth_loading_seed)
  }
  if (truth_basis_method == "iid_common_splinet_stiefel") {
    true_params$truth_basis_dim <- as.integer(truth_basis_dim)
    true_params$truth_basis_smoothness_alpha <- 3
  }
  if (mean_truth_method == "independent_splinet") {
    true_params$mean_truth_method <- mean_truth_method
    if (!is.null(mean_truth_seed)) {
      true_params$mean_truth_seed <- as.integer(mean_truth_seed)
    }
    true_params$mean_basis_dim <- as.integer(mean_basis_dim)
    true_params$mean_smoothness_alpha <- as.numeric(mean_smoothness_alpha)
    true_params$mean_truth_family <-
      "independent_periodic_continuous_L2_orthonormal_B_spline"
    true_params$mean_dense_list <- mean_truth_object$dense
    true_params$mean_basis_meta <- mean_truth_object$basis_meta
    true_params$mean_common_meta <- mean_truth_object$common_meta
    true_params$mean_truth_directly_evaluated <- TRUE
    true_params$mean_periodic_phase_draw_retained_for_pairing <- TRUE
  }

  create_named_list(Y, Z, time_obs, C, true_params)
}
