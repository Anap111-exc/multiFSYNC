testthat::test_that("assignment and boundary guards are deterministic", {
  cost <- matrix(c(1, 2, 9, 2, 9, 1), nrow = 2, byrow = TRUE)
  testthat::expect_identical(
    multiFSYNC:::.hungarian_assignment(cost), c(1L, 3L))
  testthat::expect_error(
    multiFSYNC::get_grid_objects(list(seq(0, 1, length.out = 8)),
                                 K = 1, n_g = 15, format_univ = TRUE),
    "at least 2")
  testthat::expect_error(
    multiFSYNC:::check_annealing(c(1, 1.9, 1), FALSE), "at least 2")
  testthat::expect_error(
    multiFSYNC:::check_annealing(c(1, 0.9, 3), FALSE), "greater than 1")
  testthat::expect_error(
    multiFSYNC:::generate_identified_loadings(2, 1, 1, 2), "p >")
})

testthat::test_that("simulation truth has the declared shape and factor scale", {
  structured <- multiFSYNC::simulate_multi_study_structured(
    S = 2, n_s = c(2, 2), p = 3, d = 1,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 4, n_obs = 8, identified_loadings = TRUE, seed = 5)
  testthat::expect_identical(
    dim(structured$true_params$beta_true_values[[1]][[1]]), c(8L, 3L))

  sim <- multiFSYNC::simulate_multi_study_osullivan(
    S = 2, n_s = c(2, 2), p = 5, d = 0,
    L_f = 1, L_s = 1, M_f = 2, M_s = list(2, 2),
    K = 5, n_obs = 10, use_explicit = TRUE,
    identified_loadings = TRUE, seed = 3)
  t <- seq(0, 1, length.out = 1001)
  knots <- unname(stats::quantile(sort(unique(unlist(sim$time_obs))),
                                  seq(0, 1, length = 5)[-c(1, 5)]))
  Cg <- cbind(1, t, multiFSYNC:::ZOSull(
    t, range.x = c(0, 1), intKnots = knots))
  dt <- diff(t)
  w <- c(dt[1] / 2, (dt[-1] + dt[-length(dt)]) / 2,
         dt[length(dt)] / 2)
  G <- crossprod(Cg, w * Cg)
  ivar <- function(Nu) sum(diag(crossprod(Nu, G %*% Nu)))
  testthat::expect_equal(
    ivar(sim$true_params$nu_phi_true[[1]]), 1, tolerance = 1e-8)
  testthat::expect_equal(
    ivar(sim$true_params$nu_psi_true[[1]][[1]]), 1, tolerance = 1e-8)

  specific_only <- multiFSYNC::simulate_multi_study_osullivan(
    S = 2, n_s = c(2, 2), p = 3, d = 0,
    L_f = 0, L_s = 1, M_f = integer(0), M_s = list(1, 1),
    K = 4, n_obs = 8, use_explicit = TRUE, seed = 4)
  testthat::expect_length(specific_only$true_params$nu_phi_true, 0)
})

testthat::test_that("tiny fit returns refreshed precision and uncertainty diagnostics", {
  sim <- multiFSYNC::simulate_multi_study_osullivan(
    S = 2, n_s = c(2, 2), p = 5, d = 0,
    L_f = 1, L_s = 1, M_f = 1, M_s = list(1, 1),
    K = 5, n_obs = 8, use_explicit = TRUE,
    identified_loadings = TRUE, seed = 7)
  testthat::expect_warning(
    fit <- multiFSYNC::bayesSYNC_multi(
      sim$Y, time_obs = sim$time_obs,
      L_f = 1, L_s = 1, M_f = 1, M_s = list(1, 1),
      K = 5, n_g = 17, anneal = NULL, maxit = 3,
      verbose = FALSE, bool_scale = FALSE, seed = 8),
    "Max iterations")
  testthat::expect_gt(
    max(abs(fit$inv_Sigma_q_nu_mu[[1]][[1]] - diag(7))), 0.1)
  expected <- c(
    "list_full_posterior_integrated_variance",
    "list_omitted_uncertainty_variance", "list_rank_cap",
    "list_full_posterior_integrated_variance_spec",
    "list_omitted_uncertainty_variance_spec", "list_rank_cap_spec")
  testthat::expect_true(all(expected %in% names(fit)))
  testthat::expect_error(
    multiFSYNC::bayesSYNC_multi(
      sim$Y, time_obs = sim$time_obs,
      L_f = 1, L_s = 1, M_f = 1, M_s = list(1, 1),
      K = 5, n_g = 17, anneal = NULL, maxit = 2,
      tol_abs = -1, verbose = FALSE),
    "tol_abs")
})
