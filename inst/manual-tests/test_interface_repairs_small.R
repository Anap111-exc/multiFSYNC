# Lightweight regression tests for interface, post-processing and simulation
# repairs.  Designed to finish in seconds on a small Windows laptop.
library(splines)
library(parallel)

working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(working_dir, "DESCRIPTION"))) {
  helper_file <- file.path(working_dir, "inst", "manual-tests",
                           "helper_project_root.R")
} else {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(script_arg)) stop("Run this test from the multiFSYNC package root.")
  test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]),
                                    winslash = "/", mustWork = TRUE))
  helper_file <- file.path(test_dir, "helper_project_root.R")
}
source(helper_file)
root <- project_root
for (f in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
near <- function(x, y, tol = 1e-7) isTRUE(all.equal(x, y, tolerance = tol))
expect_error <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
pass <- function(x) cat("PASS ", x, "\n", sep = "")

# 1. Full posterior second moments may have rank > M, but reported FPCA rank
# must not exceed M.
t <- seq(0, 1, length.out = 21)
Cg <- cbind(1, t, sin(2 * pi * t))
fp <- .posterior_factor_fpca(
  Cg, .trap_weights(t),
  coef_mean = matrix(c(0.3, -0.2, 0.8), 3, 1),
  coef_cov = list(diag(c(0.2, 0.1, 0.3))),
  score_mean = matrix(c(-0.5, 0.2, 0.7), 3, 1),
  score_cov = replicate(3, matrix(0.25, 1, 1), simplify = FALSE))
assert(ncol(fp$functions) <= 1L && fp$rank_cap == 1L,
       "Posterior FPCA exceeded the input M.")
assert(fp$full_posterior_integrated_variance + 1e-10 >= fp$integrated_variance,
       "Rank-capped variance diagnostic is inconsistent.")
pass("posterior FPCA uses complete moments but reports at most M directions")

# 2. Means are study-by-variable, and canonical scientific outputs are restored
# to the original response scale.
mu <- list(
  list(c(0.1, 0.2, 0), c(-0.2, 0.1, 0.3)),
  list(c(0.4, -0.1, 0.2), c(0.3, 0.2, -0.2)))
beta <- list(list(c(0.2, 0.1, 0)), list(c(-0.1, 0.3, 0.2)))
eta <- list(list(matrix(c(-0.4, 0.2), 2, 1)),
            list(matrix(c(0.3, 0.5), 2, 1)))
eta_cov <- list(
  list(replicate(2, matrix(0.1, 1, 1), simplify = FALSE)),
  list(replicate(2, matrix(0.1, 1, 1), simplify = FALSE)))
pp <- orthonormalise_multi(
  C_g = Cg, time_g = t, mu_q_nu_mu = mu,
  mu_q_nu_beta = beta,
  mu_q_nu_phi = list(matrix(c(0.2, 0.1, 0.7), 3, 1)),
  mu_q_nu_psi = NULL, mu_q_zeta = eta, Sigma_q_zeta = eta_cov,
  mu_q_xi = NULL, Sigma_q_xi = NULL,
  mu_q_a = matrix(c(0.5, 0.3), 2, 1),
  mu_q_b_specific = NULL, mu_q_gamma_a = matrix(0.7, 2, 1),
  mu_q_gamma_b = NULL, S = 2, n_s = c(2, 2), p = 2, d = 1,
  L_f = 1, L_s = 0, M_f = 1, M_s = list(integer(0), integer(0)),
  Sigma_q_nu_phi = list(list(diag(c(0.02, 0.01, 0.03)))),
  response_center = c(10, -3), response_scale = c(2, 4))
assert(length(pp$list_mu_hat) == 2L && all(lengths(pp$list_mu_hat) == 2L),
       "Mean-function output is not indexed by study and variable.")
assert(near(pp$list_mu_hat[[2]][[1]], 10 + 2 * as.vector(Cg %*% mu[[2]][[1]])),
       "Mean function was not restored to the original scale.")
assert(near(pp$list_beta_hat[[2]][[1]], 4 * as.vector(Cg %*% beta[[2]][[1]])),
       "Beta function was not restored to the original scale.")
pass("study-specific means and beta functions return on the original scale")

# 3. Grid and dimension validation.
go <- get_grid_objects(list(seq(0, 1, length.out = 12)), K = 4,
                       time_g = seq(0, 1, length.out = 17),
                       format_univ = TRUE)
assert(go$n_g == 17L, "n_g does not match a supplied time_g.")
assert(expect_error(get_grid_objects(list(c(0, 0.5, 1)), K = 4,
                                     time_g = c(0, 0.5, 0.5, 1),
                                     format_univ = TRUE)),
       "Duplicate time_g values were not rejected.")
pass("custom grid length and [0,1] validation are enforced")

# 4. Identified loading generator and shared O'Sullivan knots.
set.seed(11)
ld <- generate_identified_loadings(12, 2, 2, 2, sparse = TRUE,
                                   prop_sparse = 0.5)
G1 <- crossprod(cbind(ld$a, ld$b[[1]]))
assert(max(abs(G1 - diag(diag(G1)))) < 1e-10,
       "Generated identified loadings do not have diagonal Gram matrix.")
assert(all(diff(diag(G1)) < 0), "Reference loading norms are not strictly decreasing.")
assert(!near(ld$b[[1]], ld$b[[2]]), "Specific loadings are identical across studies.")

sim <- simulate_multi_study_osullivan(
  S = 2, n_s = c(3, 3), p = 6, d = 0, L_f = 1, L_s = 1,
  M_f = 1, M_s = list(1, 1), K = 4, n_obs = 12,
  common_grid = FALSE, use_explicit = TRUE,
  mean_structure = "study_specific", specific_time_heterogeneity = TRUE,
  identified_loadings = TRUE, seed = 12)
knots <- sim$true_params$int_knots
C11 <- cbind(1, sim$time_obs[[1]][[1]], ZOSull(sim$time_obs[[1]][[1]],
  range.x = c(0, 1), intKnots = knots))
assert(near(C11, sim$C[[1]][[1]]), "Simulation did not use pooled common knots.")
assert(!near(sim$true_params$nu_psi_true[[1]][[1]],
             sim$true_params$nu_psi_true[[2]][[1]]),
       "Explicit study-specific time functions are still identical.")
assert(!near(sim$true_params$nu_mu_true[[1]][[1]],
             sim$true_params$nu_mu_true[[2]][[1]]),
       "Study-specific explicit means are still identical.")
pass("formal simulation supports Gram conditions, pooled knots and genuine study heterogeneity")

# 5. Tiny end-to-end interface checks: one subject does not produce NA scaling;
# status/original-scale fields are returned; omega defaults agree.
one <- simulate_multi_study_data(
  S = 1, n_s = 1, p = 1, d = 0, L_f = 0, L_s = 0,
  M_f = integer(0), M_s = list(integer(0)), K = 4, n_obs = 10,
  seed = 13)
fit <- suppressWarnings(bayesSYNC_multi(
  one$Y, time_obs = one$time_obs, L_f = 0, L_s = 0,
  M_f = integer(0), M_s = list(integer(0)), K = 4,
  anneal = NULL, maxit = 2, bool_scale = TRUE, verbose = FALSE,
  list_hyper = set_hyper()))
assert(is.finite(fit$response_scale) && fit$response_scale == 1,
       "Single-subject scaling still returns NA.")
assert(fit$d_0 == fit$p && all(c("converged", "annealing_completed",
       "final_temperature", "sigsq_eps_hat") %in% names(fit)),
       "Fit status/original-scale fields or default d_0 are missing.")
assert(expect_error(bayesSYNC_multi(
  one$Y, time_obs = one$time_obs, L_f = 0, L_s = 0,
  M_f = integer(0), M_s = list(integer(0)), K = 4,
  anneal = c(1, 1.5, 10), maxit = 12, verbose = FALSE)),
  "A fit without post-annealing iterations was not rejected.")
pass("tiny fit handles scaling, hyperparameter defaults and annealing status")

# 6. Both omega hierarchies complete a tiny shared-plus-specific fit after the
# interface and stopping-rule changes.
fit_factor <- suppressWarnings(bayesSYNC_multi(
  sim$Y, time_obs = sim$time_obs, L_f = 1, L_s = 1,
  M_f = 1, M_s = list(1, 1), K = 4, anneal = NULL,
  maxit = 2, bool_scale = FALSE, bool_var_spec_prob = FALSE,
  verbose = FALSE, seed = 21))
fit_variable <- suppressWarnings(bayesSYNC_multi(
  sim$Y, time_obs = sim$time_obs, L_f = 1, L_s = 1,
  M_f = 1, M_s = list(1, 1), K = 4, anneal = NULL,
  maxit = 2, bool_scale = FALSE, bool_var_spec_prob = TRUE,
  verbose = FALSE, seed = 21))
assert(length(fit_factor$c_1_omega_a) == 1L &&
       identical(dim(fit_variable$c_1_omega_a), c(6L, 1L)),
       "Factor- and variable-level omega outputs have incorrect dimensions.")
assert(all(vapply(fit_variable$list_mu_hat, length, integer(1)) == 6L),
       "End-to-end fit lost study-specific mean functions.")
pass("factor- and variable-level omega branches finish a tiny full fit")

cat("ALL PASS interface-repair small-sample tests\n")
