# True d=0 generation. All changes are isolated experiment configuration/calibration.
nc_generator_config <- function(row) {
  S <- as.integer(row$S); lf <- as.integer(row$L_f)
  ls <- if (is.character(row$L_s_by_study)) as.integer(jsonlite::fromJSON(row$L_s_by_study)) else as.integer(row$L_s_by_study)
  assert(S == 2L && lf == 2L && length(ls) == S && all(ls == 2L), 'Unregistered true factor dimension')
  assert(as.integer(row$d_generated) == 0L && as.integer(row$d_fit) == 0L, 'NC d=0 required')
  list(S = S, n_s = rep(as.integer(row$n_each), S), p = as.integer(row$p), d = 0L,
    L_f = lf, L_s = ls, M_f = rep(as.integer(row$M), lf),
    M_s = lapply(ls, function(x) rep(as.integer(row$M), x)), K = as.integer(row$K),
    n_obs = 8L, common_grid = FALSE, n_obs_range = c(6L, 9L), sigma_eps = 0.3,
    bool_sparse_loadings = TRUE, prop_sparse = 0.9, score_var_decay = TRUE,
    n_dense = 201L, mean_amp = 0.6, beta_amp = 0.3, bs_degree = c(2L, 3L),
    identified_loadings = TRUE, sparsity_mode = 'fixed',
    truth_basis_method = 'iid_common_splinet_stiefel', truth_basis_dim = 8L,
    mean_truth_method = 'independent_splinet', mean_basis_dim = 8L,
    mean_smoothness_alpha = 1.519826404906808)
}

nc_calibrate_all_loadings <- function(data) {
  truth <- data$true_params
  calibrate <- function(x) {
    active <- sum(x != 0); norm <- sqrt(sum(x^2))
    assert(active >= 1L && is.finite(norm) && norm > 1e-14, 'Empty or nonfinite loading column')
    3 * sqrt(active / 10) * x / norm
  }
  for (l in seq_len(truth$L_f)) truth$a_true[, l] <- calibrate(truth$a_true[, l])
  for (s in seq_len(truth$S)) for (l in seq_len(ncol(truth$b_true[[s]])))
    truth$b_true[[s]][, l] <- calibrate(truth$b_true[[s]][, l])
  truth$loading_calibration <- 'v0f_constant_per_active_loading_all_columns'
  for (s in seq_len(truth$S)) for (i in seq_len(truth$n_s[s])) {
    signal <- truth$mu_true_values[[s]][[i]] + truth$beta_true_values[[s]][[i]] +
      truth$f_true_values[[s]][[i]] %*% t(truth$a_true) +
      truth$g_true_values[[s]][[i]] %*% t(truth$b_true[[s]])
    truth$signal_true_values[[s]][[i]] <- signal
    observed <- signal + truth$noise_true_values[[s]][[i]]
    data$Y[[s]][[i]] <- lapply(seq_len(truth$p), function(j) observed[, j])
  }
  data$true_params <- truth
  data
}

nc_generate_dataset <- function(row, experiment_id) {
  config <- nc_generator_config(row)
  data <- do.call(multiFSYNC::simulate_multi_study_structured, c(config,
    list(seed = as.integer(row$data_seed), truth_basis_seed = as.integer(row$truth_basis_seed),
      truth_loading_seed = as.integer(row$truth_loading_seed), mean_truth_seed = as.integer(row$mean_truth_seed))))
  assert(is.null(data$Z) && data$true_params$d_use == 0L, 'Generator retained hidden covariates')
  data <- nc_calibrate_all_loadings(data)
  dims <- config[c('S', 'n_s', 'p', 'd', 'L_f', 'L_s', 'M_f', 'M_s', 'K')]
  obs <- list(bundle_class = 'nc_observation_only_v1', experiment_id = experiment_id,
    scenario_id = as.character(row$scenario_id), data_id = as.character(row$data_id),
    Y = data$Y, Z = NULL, time_obs = data$time_obs, dimensions = dims,
    generated_d = 0L, truth_available_to_fit = FALSE, test_response_available_to_fit = FALSE)
  list(data = data, observation = obs, config = config)
}

nc_generation_qc <- function(pair, row) {
  truth <- pair$data$true_params; obs <- pair$observation
  target_support <- max(1L, as.integer(round(0.1 * truth$p)))
  load_rows <- list()
  make_row <- function(x, role, study, factor) {
    active <- sum(x != 0); normsq <- sum(x^2)
    data.frame(data_id = row$data_id, role = role, study = study, factor = factor,
      support = active, squared_norm = normsq, target_squared_norm = 9 * active / 10,
      pass = active == target_support && abs(normsq - 9 * active / 10) < 1e-8)
  }
  for (l in seq_len(truth$L_f)) load_rows[[length(load_rows) + 1L]] <- make_row(truth$a_true[, l], 'shared', 0L, l)
  gram_error <- 0
  for (s in seq_len(truth$S)) {
    for (l in seq_len(ncol(truth$b_true[[s]]))) load_rows[[length(load_rows) + 1L]] <- make_row(truth$b_true[[s]][, l], 'specific', s, l)
    AB <- cbind(truth$a_true, truth$b_true[[s]]); gg <- crossprod(AB); diag(gg) <- 0
    gram_error <- max(gram_error, abs(gg))
  }
  loadings <- do.call(rbind, load_rows)
  basis <- stage7m_truth_basis(truth); recon_error <- 0
  for (s in seq_len(truth$S)) for (i in seq_len(truth$n_s[s])) {
    ans <- stage7m_true_parts(truth, obs, s, i, obs$time_obs[[s]][[i]], basis)
    recon_error <- max(recon_error, abs(ans$signal - truth$signal_true_values[[s]][[i]]))
  }
  lambda <- unlist(c(truth$lambda_phi_true, unlist(truth$lambda_psi_true, recursive = FALSE)))
  lambda_matrix <- matrix(lambda, nrow = 2L)
  checks <- c(generated_d0 = truth$d_use == 0L && is.null(pair$data$Z),
    observation_d0 = obs$dimensions$d == 0L && is.null(obs$Z),
    finite_observations = all_finite(obs$Y), finite_truth = all_finite(truth),
    all_loading_columns_calibrated = all(loadings$pass),
    all_within_study_loading_grams_diagonal = gram_error < 1e-8,
    score_weights = all(abs(lambda_matrix[1, ] - 0.8) < 1e-12) && all(abs(lambda_matrix[2, ] - 0.2) < 1e-12),
    continuous_truth_reconstruction = recon_error < 1e-9,
    subject_counts = identical(as.integer(lengths(obs$Y)), as.integer(truth$n_s)),
    irregular_time_counts_6_9 = all(lengths(unlist(obs$time_obs, recursive = FALSE)) %in% 6:9),
    time_domain_0_1 = all(unlist(obs$time_obs) >= 0 & unlist(obs$time_obs) <= 1),
    actual_generated_basis_columns = all(vapply(unlist(pair$data$C, recursive = FALSE), ncol, integer(1)) == row$K + 2L))
  list(qc = data.frame(data_id = row$data_id, check_id = names(checks), pass = unname(checks)),
    loading_qc = loadings, reconstruction_error = recon_error, gram_error = gram_error)
}
