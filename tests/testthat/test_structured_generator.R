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
