test_that("Jaoua-style generator is reproducible and uses irregular grids", {
  args <- list(
    S = 2, n_s = c(12, 13), p = 5, d = 1,
    L_f = 2, L_s = 1, M_f = c(2, 3), M_s = list(2, 2),
    K = 5, n_obs = 8, common_grid = FALSE, n_obs_range = c(5, 10),
    sigma_eps = 0.1, sparsity_mode = "factor_beta",
    omega_beta = c(2, 4), n_dense = 160, seed = 73101)

  dat_1 <- do.call(simulate_multi_study_structured, args)
  dat_2 <- do.call(simulate_multi_study_structured, args)

  expect_identical(dat_1, dat_2)
  n_i <- unlist(dat_1$true_params$n_obs_by_subject, use.names = FALSE)
  expect_true(all(n_i >= 5L & n_i <= 10L))
  expect_gt(length(unique(n_i)), 1L)
  expect_true(all(vapply(
    unlist(dat_1$time_obs, recursive = FALSE),
    function(x) all(diff(x) > 0) && all(x >= 0 & x <= 1),
    logical(1))))
  expect_identical(dat_1$true_params$truth_basis_family,
                   "ordinary_B_spline_external")
  expect_false(dat_1$true_params$truth_projected_to_fitted_basis)

  used_degrees <- unlist(lapply(
    dat_1$true_params$phi_basis_meta,
    function(meta) vapply(meta$components, `[[`, integer(1), "degree")))
  expect_true(all(c(2L, 3L) %in% used_degrees))
})

test_that("scheme-1 theoretical and computational truths are equivalent", {
  dat <- simulate_multi_study_structured(
    S = 2, n_s = c(40, 45), p = 4, d = 0,
    L_f = 2, L_s = 1, M_f = c(3, 2), M_s = list(2, 3),
    K = 5, n_obs = 7, common_grid = FALSE, n_obs_range = c(5, 9),
    bool_sparse_loadings = FALSE, n_dense = 180, seed = 73102)
  tp <- dat$true_params

  for (l in seq_len(tp$L_f)) {
    M_l <- tp$M_f[l]
    expected_lambda <- (1 / seq_len(M_l)^2) /
      sum(1 / seq_len(M_l)^2)
    expect_equal(tp$lambda_phi_true[[l]], expected_lambda, tolerance = 1e-14)
    expect_equal(tp$phi_basis_meta[[l]]$gram, diag(M_l), tolerance = 1e-8)
    expect_equal(
      tp$theta_dense_list[[l]],
      sweep(tp$phi_dense_list[[l]], 2L, sqrt(expected_lambda), `*`),
      tolerance = 1e-14)

    for (s in seq_len(tp$S)) {
      expect_identical(tp$zeta_true[[s]][[l]], tp$eta_true[[s]][[l]])
      expect_equal(
        tp$zeta_fpca_true[[s]][[l]],
        sweep(tp$zeta_true[[s]][[l]], 2L, sqrt(expected_lambda), `*`),
        tolerance = 1e-14)
      expect_equal(
        tp$theta_dense_list[[l]] %*% tp$zeta_true[[s]][[l]][1L, ],
        tp$phi_dense_list[[l]] %*% tp$zeta_fpca_true[[s]][[l]][1L, ],
        tolerance = 1e-12)
    }
  }

  for (s in seq_len(tp$S)) {
    for (l in seq_len(tp$L_s)) {
      M_sl <- tp$M_s[[s]][l]
      expected_lambda <- (1 / seq_len(M_sl)^2) /
        sum(1 / seq_len(M_sl)^2)
      expect_equal(tp$lambda_psi_true[[s]][[l]], expected_lambda,
                   tolerance = 1e-14)
      expect_equal(tp$psi_basis_meta[[s]][[l]]$gram, diag(M_sl),
                   tolerance = 1e-8)
      expect_identical(tp$xi_true[[s]][[l]], tp$chi_true[[s]][[l]])
      expect_equal(
        tp$kappa_dense_list[[s]][[l]] %*% tp$xi_true[[s]][[l]][1L, ],
        tp$psi_dense_list[[s]][[l]] %*% tp$xi_fpca_true[[s]][[l]][1L, ],
        tolerance = 1e-12)
    }
  }

  all_eta <- unlist(tp$eta_true, recursive = TRUE, use.names = FALSE)
  expect_lt(abs(stats::var(all_eta) - 1), 0.2)
})

test_that("latent factors use distinct truth blocks", {
  dat <- simulate_multi_study_structured(
    S = 2, n_s = c(6, 6), p = 4, d = 0,
    L_f = 3, L_s = 1, M_f = c(2, 2, 2), M_s = list(2, 2),
    K = 5, n_obs = 6, common_grid = FALSE,
    bool_sparse_loadings = FALSE, n_dense = 140, seed = 73103)
  tp <- dat$true_params

  shared_ids <- vapply(tp$phi_basis_meta, `[[`, integer(1), "block_id")
  specific_ids <- unlist(lapply(
    tp$psi_basis_meta,
    function(x) vapply(x, `[[`, integer(1), "block_id")))
  expect_equal(length(unique(c(shared_ids, specific_ids))),
               length(c(shared_ids, specific_ids)))

  for (l_1 in 1:2) {
    for (l_2 in (l_1 + 1):3) {
      expect_gt(max(abs(tp$phi_dense_list[[l_1]] -
                        tp$phi_dense_list[[l_2]])), 1e-3)
    }
  }
  expect_gt(max(abs(tp$psi_dense_list[[1]][[1]] -
                    tp$psi_dense_list[[2]][[1]])), 1e-3)
})

test_that("stored signal and noise exactly reconstruct every response", {
  dat <- simulate_multi_study_structured(
    S = 2, n_s = c(5, 4), p = 5, d = 1,
    L_f = 2, L_s = 1, M_f = c(2, 1), M_s = list(2, 2),
    K = 5, n_obs = 7, common_grid = FALSE, n_obs_range = c(5, 8),
    sigma_eps = 0.15, sparsity_mode = "factor_beta",
    omega_beta = c(2, 3), n_dense = 140, seed = 73104)
  tp <- dat$true_params

  for (s in seq_len(tp$S)) {
    for (i in seq_len(tp$n_s[s])) {
      y_mat <- do.call(cbind, dat$Y[[s]][[i]])
      expect_equal(
        y_mat,
        tp$signal_true_values[[s]][[i]] + tp$noise_true_values[[s]][[i]],
        tolerance = 0)

      reconstructed_signal <- tp$mu_true_values[[s]][[i]] +
        tp$beta_true_values[[s]][[i]]
      for (l in seq_len(tp$L_f)) {
        reconstructed_signal <- reconstructed_signal +
          outer(tp$f_true_values[[s]][[i]][, l], tp$a_true[, l])
      }
      for (l in seq_len(tp$L_s)) {
        reconstructed_signal <- reconstructed_signal +
          outer(tp$g_true_values[[s]][[i]][, l], tp$b_true[[s]][, l])
      }
      expect_equal(reconstructed_signal,
                   tp$signal_true_values[[s]][[i]], tolerance = 1e-13)
    }
  }
})

test_that("legacy structured-generator arguments remain supported", {
  dat <- simulate_multi_study_structured(
    S = 1, n_s = 5, p = 4, d = 0,
    L_f = 1, L_s = 0, M_f = 2, M_s = list(integer(0)),
    K = 5, n_obs = 9, common_grid = TRUE,
    bool_sparse_loadings = TRUE, prop_sparse = 0.5,
    bs_degree = 3, seed = 73105)

  expect_true(all(lengths(dat$time_obs[[1]]) == 9L))
  expect_equal(sum(dat$true_params$gamma_a_true[, 1]), 2L)
  expect_identical(dat$true_params$sparsity_mode, "fixed")
  expect_null(dat$true_params$xi_true)
  expect_null(dat$true_params$psi_dense_list)
})

test_that("legacy dense-Gram truth remains the default and reproducible", {
  args <- list(
    S = 2, n_s = c(4, 5), p = 3, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 6, common_grid = FALSE,
    bool_sparse_loadings = FALSE, n_dense = 100, seed = 73106)

  implicit <- do.call(simulate_multi_study_structured, args)
  explicit <- do.call(
    simulate_multi_study_structured,
    c(args, list(truth_basis_method = "legacy_dense_gram")))

  expect_identical(implicit, explicit)
  expect_identical(implicit$true_params$truth_basis_method,
                   "legacy_dense_gram")
})

test_that("Jaoua Splinets truth is exact, reproducible and directly evaluated", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 2, n_s = c(4, 5), p = 3, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 8, n_obs = 6, common_grid = FALSE, n_obs_range = c(5, 7),
    bool_sparse_loadings = FALSE, n_dense = 101, seed = 73107,
    truth_basis_method = "jaoua_splinets")

  dat_1 <- do.call(simulate_multi_study_structured, args)
  dat_2 <- do.call(simulate_multi_study_structured, args)
  expect_identical(dat_1, dat_2)

  tp <- dat_1$true_params
  expect_identical(tp$truth_basis_method, "jaoua_splinets")
  expect_identical(tp$truth_basis_family,
                   "jaoua_splinets_orthonormal_B_spline")

  all_meta <- c(tp$phi_basis_meta, unlist(tp$psi_basis_meta, recursive = FALSE))
  expect_equal(length(all_meta), 3L)
  expect_identical(vapply(all_meta, `[[`, integer(1), "block_id"), 1:3)
  for (meta in all_meta) {
    expect_equal(meta$gram, diag(2), tolerance = 1e-12)
    expect_identical(meta$method, "jaoua_splinets")
    expect_identical(meta$common_selection_pool_size, 3L)
    expect_identical(meta$source_binding$Splinets_version,
                     as.character(packageVersion("Splinets")))
  }

  # Reconstruct each exact spline object from stored non-truth metadata and
  # verify that observed factor values did not come from dense-grid interpolation.
  reconstruct_phi <- function(meta, t) {
    so <- Splinets::splinet(
      meta$knots, degree = meta$degree,
      periodic = meta$periodic, norm = TRUE)
    evaluated <- Splinets::evspline(so$os, x = t)[, -1L, drop = FALSE]
    evaluated[, meta$selected_indices, drop = FALSE]
  }
  for (s in seq_len(tp$S)) {
    for (i in seq_len(tp$n_s[s])) {
      t_i <- dat_1$time_obs[[s]][[i]]
      phi_i <- reconstruct_phi(tp$phi_basis_meta[[1L]], t_i)
      theta_i <- sweep(
        phi_i, 2L, sqrt(tp$lambda_phi_true[[1L]]), FUN = "*")
      expect_equal(
        tp$f_true_values[[s]][[i]][, 1L],
        as.vector(theta_i %*% tp$eta_true[[s]][[1L]][i, ]),
        tolerance = 1e-13)

      psi_i <- reconstruct_phi(tp$psi_basis_meta[[s]][[1L]], t_i)
      kappa_i <- sweep(
        psi_i, 2L, sqrt(tp$lambda_psi_true[[s]][[1L]]), FUN = "*")
      expect_equal(
        tp$g_true_values[[s]][[i]][, 1L],
        as.vector(kappa_i %*% tp$chi_true[[s]][[1L]][i, ]),
        tolerance = 1e-13)
    }
  }
})

test_that("Jaoua Splinets support limits fail explicitly", {
  skip_if_not_installed("Splinets")
  base_args <- list(
    S = 2, n_s = c(3, 3), p = 2, d = 0, K = 5,
    n_obs = 5, common_grid = FALSE, bool_sparse_loadings = FALSE,
    n_dense = 80, seed = 73108,
    truth_basis_method = "jaoua_splinets")

  expect_error(
    do.call(simulate_multi_study_structured, c(
      base_args,
      list(L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 3)))),
    "same FPCA component count", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured, c(
      base_args,
      list(L_f = 3, L_s = 1, M_f = c(2, 2, 2), M_s = list(2, 2)))),
    "total factor count Q between 1 and 4", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured, c(
      base_args,
      list(L_f = 1, L_s = 0, M_f = 1,
           M_s = list(integer(0), integer(0))))),
    "common component count M between 2 and 4", fixed = TRUE)
})

test_that("iid common Splinets Stiefel truth is exact and exchangeable", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 3, n_s = c(4, 4, 4), p = 5, d = 0,
    L_f = 2, L_s = c(1, 2, 1), M_f = c(1, 3),
    M_s = list(2, c(1, 4), 3),
    K = 5, n_obs = 6, common_grid = FALSE, n_obs_range = c(5, 7),
    bool_sparse_loadings = FALSE, n_dense = 101, seed = 75101,
    truth_basis_seed = 85101,
    truth_basis_method = "iid_common_splinet_stiefel",
    truth_basis_dim = 8)

  dat_1 <- do.call(simulate_multi_study_structured, args)
  dat_2 <- do.call(simulate_multi_study_structured, args)
  expect_identical(dat_1, dat_2)

  tp <- dat_1$true_params
  expect_identical(tp$truth_basis_method,
                   "iid_common_splinet_stiefel")
  expect_identical(tp$truth_basis_family,
                   "iid_common_splinet_stiefel_orthonormal_B_spline")
  expect_identical(tp$truth_basis_dim, 8L)
  expect_identical(tp$truth_basis_smoothness_alpha, 3)

  all_meta <- c(tp$phi_basis_meta,
                unlist(tp$psi_basis_meta, recursive = FALSE))
  expect_length(all_meta, 6L)
  expect_identical(vapply(all_meta, `[[`, integer(1), "block_id"), 1:6)
  expect_identical(
    vapply(all_meta, `[[`, integer(1), "degree"), rep(3L, 6L))
  expect_true(all(vapply(all_meta, function(x) !x$periodic, logical(1))))
  expect_true(all(vapply(
    all_meta[-1L],
    function(x) identical(x$knots, all_meta[[1L]]$knots), logical(1))))

  expected_M <- c(1L, 3L, 2L, 1L, 4L, 3L)
  for (q in seq_along(all_meta)) {
    meta <- all_meta[[q]]
    expect_equal(meta$gram, diag(expected_M[[q]]), tolerance = 1e-12)
    expect_equal(meta$common_basis_gram, diag(8), tolerance = 1e-12)
    expect_equal(crossprod(meta$stiefel_coefficients),
                 diag(expected_M[[q]]), tolerance = 1e-12)
    expect_identical(meta$common_basis_dim, 8L)
    expect_identical(meta$smoothness_alpha, 3)
    expect_equal(dim(meta$second_derivative_gram), c(8L, 8L))
    expect_length(meta$roughness_eigenvalues, 8L)
    expect_length(meta$spectral_weights, 8L)
    expect_true(all(meta$spectral_weights > 0 &
                    meta$spectral_weights <= 1))
    expect_match(meta$source_binding$construction,
                 "Splinets::splinet")
  }

  # Independent random subspaces share one basis law but are neither forced
  # equal nor forced mutually orthogonal.
  cross_12 <- crossprod(
    all_meta[[1L]]$stiefel_coefficients,
    all_meta[[2L]]$stiefel_coefficients)
  expect_gt(max(abs(cross_12)), 1e-4)
  expect_lt(max(svd(cross_12, nu = 0, nv = 0)$d), 1 - 1e-8)

  # Rebuild a stored exact Splinets linear combination and confirm that the
  # generated process at irregular observation times used direct evaluation.
  rebuild_phi <- function(meta, t) {
    common <- Splinets::splinet(
      meta$knots, degree = meta$degree,
      periodic = meta$periodic, norm = TRUE)
    block <- Splinets::lincomb(
      common$os, t(meta$stiefel_coefficients), reduced = TRUE)
    Splinets::evspline(block, x = t)[, -1L, drop = FALSE]
  }
  t_i <- dat_1$time_obs[[1L]][[1L]]
  phi_i <- rebuild_phi(tp$phi_basis_meta[[1L]], t_i)
  theta_i <- sweep(
    phi_i, 2L, sqrt(tp$lambda_phi_true[[1L]]), FUN = "*")
  expect_equal(
    tp$f_true_values[[1L]][[1L]][, 1L],
    as.vector(theta_i %*% tp$eta_true[[1L]][[1L]][1L, ]),
    tolerance = 1e-13)
})

test_that("iid common Splinets Stiefel basis dimension is validated", {
  skip_if_not_installed("Splinets")
  expect_error(
    simulate_multi_study_structured(
      S = 1, n_s = 3, p = 2, d = 0,
      L_f = 1, L_s = 0, M_f = 4, M_s = list(integer(0)),
      K = 5, n_obs = 5, common_grid = FALSE,
      bool_sparse_loadings = FALSE, n_dense = 80, seed = 75102,
      truth_basis_method = "iid_common_splinet_stiefel",
      truth_basis_dim = 3),
    "truth_basis_dim must be one integer at least", fixed = TRUE)

  expect_error(
    construct_iid_common_splinet_stiefel_blocks(
      M_by_block = c(2L, 2L), n_g = 80, basis_dim = 8,
      smoothness_alpha = -1),
    "smoothness_alpha must be one non-negative finite number", fixed = TRUE)
})

test_that("zero smoothness weight preserves the version-1 uniform draw", {
  skip_if_not_installed("Splinets")
  set.seed(75103)
  expected_gaussian <- matrix(rnorm(8L * 2L), 8L, 2L)
  expected_qr <- qr(expected_gaussian)
  expected_Q <- qr.Q(expected_qr)[, 1:2, drop = FALSE]
  expected_sign <- sign(diag(qr.R(expected_qr))[1:2])
  expected_sign[expected_sign == 0] <- 1
  expected_Q <- sweep(expected_Q, 2L, expected_sign, FUN = "*")

  set.seed(75103)
  block <- construct_iid_common_splinet_stiefel_blocks(
    M_by_block = 2L, n_g = 80, basis_dim = 8,
    smoothness_alpha = 0)[[1L]]
  expect_identical(block$smoothness_alpha, 0)
  expect_equal(block$stiefel_coefficients, expected_Q, tolerance = 0)
  expect_equal(block$spectral_weights, rep(1, 8), tolerance = 0)
})

test_that("route-B spectral weighting reduces sampled spline roughness", {
  skip_if_not_installed("Splinets")
  mean_roughness <- function(blocks) {
    mean(unlist(lapply(blocks, function(block) {
      Q <- block$stiefel_coefficients
      diag(t(Q) %*% block$second_derivative_gram %*% Q)
    })))
  }

  set.seed(75104)
  uniform <- construct_iid_common_splinet_stiefel_blocks(
    M_by_block = rep(2L, 12L), n_g = 101, basis_dim = 8,
    smoothness_alpha = 0)
  set.seed(75104)
  route_b <- construct_iid_common_splinet_stiefel_blocks(
    M_by_block = rep(2L, 12L), n_g = 101, basis_dim = 8,
    smoothness_alpha = 3)

  expect_lt(mean_roughness(route_b), mean_roughness(uniform))
})

test_that("truth basis seed isolates all non-basis generator randomness", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 2, n_s = c(5, 6), p = 7, d = 1,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 6, common_grid = FALSE, n_obs_range = c(5, 7),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE, prop_sparse = 0.5,
    score_var_decay = TRUE, n_dense = 101, mean_amp = 0.6,
    beta_amp = 0.3, identified_loadings = TRUE,
    seed = 73109, truth_basis_seed = 83109)

  legacy <- do.call(
    simulate_multi_study_structured,
    c(args, list(truth_basis_method = "legacy_dense_gram")))
  jaoua <- do.call(
    simulate_multi_study_structured,
    c(args, list(truth_basis_method = "jaoua_splinets")))
  left <- legacy$true_params
  right <- jaoua$true_params

  paired_fields <- c(
    "time_obs", "C", "Z", "phase_mu_mat", "phase_beta",
    "zeta_true", "xi_true", "eta_true", "chi_true",
    "zeta_fpca_true", "xi_fpca_true", "a_true", "b_true",
    "gamma_a_true", "gamma_b_true", "omega_a_true", "omega_b_true",
    "lambda_phi_true", "lambda_psi_true", "mu_true_values",
    "beta_true_values", "noise_true_values", "n_obs_by_subject",
    "int_knots"
  )
  for (field in paired_fields) {
    left_value <- if (field %in% names(legacy)) legacy[[field]] else left[[field]]
    right_value <- if (field %in% names(jaoua)) jaoua[[field]] else right[[field]]
    expect_identical(left_value, right_value, info = field)
  }
  expect_identical(left$truth_basis_seed, 83109L)
  expect_identical(right$truth_basis_seed, 83109L)
  expect_gt(max(abs(left$phi_dense_list[[1L]] -
                    right$phi_dense_list[[1L]])), 1e-6)
  expect_gt(max(abs(left$signal_true_values[[1L]][[1L]] -
                    right$signal_true_values[[1L]][[1L]])), 1e-6)
  expect_gt(max(abs(do.call(cbind, legacy$Y[[1L]][[1L]]) -
                    do.call(cbind, jaoua$Y[[1L]][[1L]]))), 1e-6)
})

test_that("truth basis seed preserves the caller RNG stream", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 2, n_s = c(3, 3), p = 3, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 5, common_grid = FALSE,
    bool_sparse_loadings = FALSE, n_dense = 80,
    seed = 73110, truth_basis_seed = 83110)

  do.call(simulate_multi_study_structured,
          c(args, list(truth_basis_method = "legacy_dense_gram")))
  after_legacy <- runif(5)
  do.call(simulate_multi_study_structured,
          c(args, list(truth_basis_method = "jaoua_splinets")))
  after_jaoua <- runif(5)
  expect_identical(after_legacy, after_jaoua)

  expect_error(
    do.call(simulate_multi_study_structured,
            utils::modifyList(args, list(truth_basis_seed = -1))),
    "truth_basis_seed must be one non-negative integer", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured,
            utils::modifyList(args, list(truth_basis_seed = 1.5))),
    "truth_basis_seed must be one non-negative integer", fixed = TRUE)
})

test_that("Jaoua-style loadings follow the registered blockwise law", {
  loading_seed <- 94101L
  args <- list(
    S = 2, n_s = c(3, 4), p = 13, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 5, common_grid = FALSE,
    bool_sparse_loadings = TRUE, identified_loadings = FALSE,
    sparsity_mode = "jaoua_beta_iid_normal", omega_beta = c(1, 10),
    n_dense = 80, seed = 74101, truth_loading_seed = loading_seed)

  dat <- do.call(simulate_multi_study_structured, args)
  draw_reference_column <- function(p) {
    omega <- stats::rbeta(1L, 1, 10)
    gamma <- stats::rbinom(p, size = 1L, prob = omega)
    if (!any(gamma == 1L)) gamma[sample(seq_len(p), 1L)] <- 1L
    value <- stats::rnorm(p) * gamma
    list(omega = omega, gamma = gamma, value = value)
  }
  set.seed(loading_seed)
  expected <- replicate(3L, draw_reference_column(args$p), simplify = FALSE)
  tp <- dat$true_params

  expect_identical(tp$sparsity_mode, "jaoua_beta_iid_normal")
  expect_identical(tp$truth_loading_seed, loading_seed)
  expect_equal(tp$omega_a_true[[1L]], expected[[1L]]$omega, tolerance = 0)
  expect_identical(as.integer(tp$gamma_a_true[, 1L]),
                   as.integer(expected[[1L]]$gamma))
  expect_equal(tp$a_true[, 1L], expected[[1L]]$value, tolerance = 0)
  for (s in 1:2) {
    expect_equal(tp$omega_b_true[[s]][[1L]],
                 expected[[s + 1L]]$omega, tolerance = 0)
    expect_identical(as.integer(tp$gamma_b_true[[s]][, 1L]),
                     as.integer(expected[[s + 1L]]$gamma))
    expect_equal(tp$b_true[[s]][, 1L],
                 expected[[s + 1L]]$value, tolerance = 0)
  }
  expect_true(all(c(colSums(tp$gamma_a_true),
                    vapply(tp$gamma_b_true, function(x) colSums(x)[[1L]],
                           numeric(1L))) >= 1L))
})

test_that("truth loading seed isolates all non-loading generator randomness", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 2, n_s = c(5, 6), p = 17, d = 1,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 6, common_grid = FALSE, n_obs_range = c(5, 7),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE, prop_sparse = 0.5,
    score_var_decay = TRUE, n_dense = 101, mean_amp = 0.6,
    beta_amp = 0.3, seed = 74102, truth_basis_seed = 84102,
    truth_loading_seed = 94102,
    truth_basis_method = "jaoua_splinets")

  controlled <- do.call(
    simulate_multi_study_structured,
    c(args, list(identified_loadings = TRUE, sparsity_mode = "fixed")))
  jaoua <- do.call(
    simulate_multi_study_structured,
    c(args, list(
      identified_loadings = FALSE,
      sparsity_mode = "jaoua_beta_iid_normal", omega_beta = c(1, 10))))
  left <- controlled$true_params
  right <- jaoua$true_params

  top_fields <- c("time_obs", "C", "Z")
  truth_fields <- c(
    "phase_mu_mat", "phase_beta", "zeta_true", "xi_true", "eta_true",
    "chi_true", "zeta_fpca_true", "xi_fpca_true", "phi_dense_list",
    "psi_dense_list", "theta_dense_list", "kappa_dense_list",
    "lambda_phi_true", "lambda_psi_true", "mu_true_values",
    "f_true_values", "g_true_values", "beta_true_values",
    "noise_true_values", "n_obs_by_subject", "int_knots",
    "sigma2_eps_true")
  for (field in top_fields) {
    expect_identical(controlled[[field]], jaoua[[field]], info = field)
  }
  for (field in truth_fields) {
    expect_identical(left[[field]], right[[field]], info = field)
  }
  expect_identical(left$truth_loading_seed, 94102L)
  expect_identical(right$truth_loading_seed, 94102L)
  expect_gt(max(abs(left$a_true - right$a_true)), 1e-6)
  expect_gt(max(abs(left$signal_true_values[[1L]][[1L]] -
                    right$signal_true_values[[1L]][[1L]])), 1e-6)
  expect_gt(max(abs(do.call(cbind, controlled$Y[[1L]][[1L]]) -
                    do.call(cbind, jaoua$Y[[1L]][[1L]]))), 1e-6)
})

test_that("truth loading seed preserves caller RNG and validates its contract", {
  base_args <- list(
    S = 2, n_s = c(3, 3), p = 9, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 5, common_grid = FALSE,
    bool_sparse_loadings = TRUE, n_dense = 80,
    seed = 74103, truth_loading_seed = 94103)

  do.call(
    simulate_multi_study_structured,
    c(base_args, list(identified_loadings = TRUE, sparsity_mode = "fixed")))
  after_controlled <- runif(5)
  do.call(
    simulate_multi_study_structured,
    c(base_args, list(
      identified_loadings = FALSE,
      sparsity_mode = "jaoua_beta_iid_normal")))
  after_jaoua <- runif(5)
  expect_identical(after_controlled, after_jaoua)

  expect_error(
    do.call(simulate_multi_study_structured,
            utils::modifyList(base_args, list(truth_loading_seed = -1))),
    "truth_loading_seed must be one non-negative integer", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured,
            utils::modifyList(base_args, list(truth_loading_seed = 1.5))),
    "truth_loading_seed must be one non-negative integer", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured, c(
      base_args,
      list(identified_loadings = TRUE,
           sparsity_mode = "jaoua_beta_iid_normal"))),
    "requires identified_loadings = FALSE", fixed = TRUE)
  expect_error(
    do.call(simulate_multi_study_structured, c(
      utils::modifyList(base_args, list(bool_sparse_loadings = FALSE)),
      list(identified_loadings = FALSE,
           sparsity_mode = "jaoua_beta_iid_normal"))),
    "requires identified_loadings = FALSE and bool_sparse_loadings = TRUE",
    fixed = TRUE)
})

test_that("invalid Jaoua-style grid controls fail explicitly", {
  expect_error(
    simulate_multi_study_structured(
      S = 1, n_s = 3, p = 2, d = 0,
      L_f = 1, L_s = 0, M_f = 1, M_s = list(integer(0)),
      common_grid = TRUE, n_obs_range = c(5, 10)),
    "requires common_grid = FALSE", fixed = TRUE)
  expect_error(
    multiFSYNC:::construct_orthonormal_eigenfunctions(2, n_g = 20),
    "n_g must be", fixed = TRUE)
})
