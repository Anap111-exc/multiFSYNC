make_practical_test_data <- function(seed = 8101L) {
  simulate_multi_study_osullivan(
    S = 1L, n_s = 4L, p = 2L, d = 0L,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, n_obs = 7L, common_grid = TRUE,
    sigma_eps = 0.08, seed = seed, use_explicit = FALSE
  )
}

small_practical_control <- function(slow = FALSE) {
  list(
    min_t1 = 3L,
    window = 2L,
    long_window = 2L,
    consecutive = 2L,
    max_t1 = if (slow) 5L else 10L,
    elbo_abs_rate = if (slow) 0 else 1e100,
    fitted_nrmse = if (slow) 0 else 1e100,
    rss_rel = if (slow) 0 else 1e100,
    ppi_max_abs = if (slow) 0 else 1e100,
    checkpoints = c(3L, 4L)
  )
}
