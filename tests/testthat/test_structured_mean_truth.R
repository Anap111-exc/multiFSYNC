test_that("historical periodic mean truth remains the exact default", {
  args <- list(
    S = 2, n_s = c(4, 5), p = 3, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 6, common_grid = FALSE,
    bool_sparse_loadings = FALSE, n_dense = 100, seed = 73106)

  implicit <- do.call(simulate_multi_study_structured, args)
  explicit <- do.call(
    simulate_multi_study_structured,
    c(args, list(mean_truth_method = "periodic_phase")))

  expect_identical(implicit, explicit)
  expect_false("mean_truth_method" %in% names(implicit$true_params))
  expect_false("mean_dense_list" %in% names(implicit$true_params))
})

test_that("independent Splinets means obey continuous geometry", {
  skip_if_not_installed("Splinets")
  set.seed(99310000L)
  means <- construct_independent_splinet_means(
    S = 1L, p = 1000L, n_g = 51L, basis_dim = 8L,
    smoothness_alpha = 1.519826404906808, amplitude = 0.6)
  meta <- means$basis_meta[[1L]]

  expect_equal(
    means$common_meta$common_basis_gram, diag(8L), tolerance = 1e-8)
  expect_true(means$common_meta$periodic)
  expect_equal(meta$continuous_energy, rep(0.225, 1000L),
               tolerance = 1e-12)
  expect_lt(max(abs(meta$continuous_integral)), 1e-10)
  expect_lt(abs(
    stats::median(meta$normalized_roughness) /
      means$common_meta$target_normalized_roughness - 1), 0.15)
  expect_equal(
    eval_independent_splinet_means(
      means$t_grid, means$evaluation[[1L]]),
    means$dense[[1L]], tolerance = 0)
})

test_that("mean truth seed creates a fully paired generator bridge", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 2, n_s = c(4, 5), p = 7, d = 1,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 9, n_obs = 6, common_grid = FALSE, n_obs_range = c(5, 7),
    sigma_eps = 0.3, bool_sparse_loadings = FALSE,
    n_dense = 101, mean_amp = 0.6, beta_amp = 0.3,
    seed = 99311001L,
    truth_basis_method = "iid_common_splinet_stiefel",
    truth_basis_seed = 99312001L,
    truth_loading_seed = 99313001L)

  periodic <- do.call(
    simulate_multi_study_structured,
    c(args, list(mean_truth_method = "periodic_phase")))
  splinet <- do.call(
    simulate_multi_study_structured,
    c(args, list(mean_truth_method = "independent_splinet",
                 mean_truth_seed = 99314001L)))
  repeated <- do.call(
    simulate_multi_study_structured,
    c(args, list(mean_truth_method = "independent_splinet",
                 mean_truth_seed = 99314001L)))

  expect_identical(splinet, repeated)
  for (field in c("time_obs", "C", "Z")) {
    expect_identical(periodic[[field]], splinet[[field]], info = field)
  }
  paired_truth <- c(
    "phase_mu_mat", "phase_beta", "zeta_true", "xi_true",
    "eta_true", "chi_true", "zeta_fpca_true", "xi_fpca_true",
    "phi_dense_list", "psi_dense_list", "theta_dense_list",
    "kappa_dense_list", "lambda_phi_true", "lambda_psi_true",
    "a_true", "b_true", "gamma_a_true", "gamma_b_true",
    "omega_a_true", "omega_b_true", "f_true_values", "g_true_values",
    "beta_true_values", "noise_true_values", "n_obs_by_subject",
    "int_knots", "sigma2_eps_true")
  for (field in paired_truth) {
    expect_identical(
      periodic$true_params[[field]], splinet$true_params[[field]],
      info = field)
  }
  expect_identical(splinet$true_params$mean_truth_method,
                   "independent_splinet")
  expect_identical(splinet$true_params$mean_truth_seed, 99314001L)
  expect_gt(max(abs(
    periodic$true_params$mu_true_values[[1L]][[1L]] -
      splinet$true_params$mu_true_values[[1L]][[1L]])), 1e-6)

  set.seed(99314001L)
  direct <- construct_independent_splinet_means(
    S = args$S, p = args$p, n_g = args$n_dense,
    basis_dim = 8L, smoothness_alpha = 1.519826404906808,
    amplitude = args$mean_amp)
  for (s in seq_len(args$S)) for (i in seq_len(args$n_s[[s]])) {
    expected_mu <- eval_independent_splinet_means(
      splinet$time_obs[[s]][[i]], direct$evaluation[[s]])
    expect_equal(
      expected_mu, splinet$true_params$mu_true_values[[s]][[i]],
      tolerance = 0)
    expect_equal(
      do.call(cbind, splinet$Y[[s]][[i]]) -
        do.call(cbind, periodic$Y[[s]][[i]]),
      splinet$true_params$mu_true_values[[s]][[i]] -
        periodic$true_params$mu_true_values[[s]][[i]],
      tolerance = 1e-13)
  }
})

test_that("mean truth seed preserves caller RNG and validates its scope", {
  skip_if_not_installed("Splinets")
  args <- list(
    S = 1, n_s = 3, p = 4, d = 0,
    L_f = 1, L_s = 0, M_f = 2, M_s = list(integer(0)),
    K = 8, n_obs = 5, common_grid = FALSE,
    bool_sparse_loadings = FALSE, n_dense = 80,
    seed = 99316001L, truth_basis_method = "iid_common_splinet_stiefel",
    truth_basis_seed = 99316002L)

  do.call(simulate_multi_study_structured, args)
  after_periodic <- runif(5)
  do.call(
    simulate_multi_study_structured,
    c(args, list(mean_truth_method = "independent_splinet",
                 mean_truth_seed = 99316003L)))
  after_splinet <- runif(5)
  expect_identical(after_periodic, after_splinet)

  expect_error(
    do.call(
      simulate_multi_study_structured,
      c(args, list(mean_truth_method = "periodic_phase",
                   mean_truth_seed = 99316003L))),
    "only used when mean_truth_method", fixed = TRUE)
  expect_error(
    do.call(
      simulate_multi_study_structured,
      c(args, list(mean_truth_method = "independent_splinet",
                   mean_truth_seed = -1L))),
    "mean_truth_seed must be one non-negative integer", fixed = TRUE)
})
