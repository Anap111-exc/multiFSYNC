# Lightweight regression tests for the O'Sullivan K+2 parameterization.
# Run directly with:
#   Rscript tests/test_osullivan_prior_small.R
#
# Deliberately avoids a full CAVI fit and external testthat/pracma dependencies.

test_started <- proc.time()[["elapsed"]]

working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(working_dir, "R", "utils_multi.R"))) {
  # Preferred on Windows: avoids decoding a Chinese --file path from the
  # process command line under a non-UTF-8 startup locale.
  project_root <- working_dir
} else {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Run from the multiFSYNC package root or invoke this file with Rscript.")
  }
  this_file <- normalizePath(
    sub("^--file=", "", file_arg),
    winslash = "/",
    mustWork = TRUE
  )
  project_root <- normalizePath(
    file.path(dirname(this_file), ".."),
    winslash = "/",
    mustWork = TRUE
  )
}
while (!file.exists(file.path(project_root, "DESCRIPTION"))) {
  parent_root <- dirname(project_root)
  if (identical(parent_root, project_root)) stop("Could not locate the package root.")
  project_root <- parent_root
}

suppressPackageStartupMessages(library(splines))

source(file.path(project_root, "R", "utils_multi.R"), local = FALSE)
source(file.path(project_root, "R", "update_mu.R"), local = FALSE)
source(file.path(project_root, "R", "update_variance.R"), local = FALSE)

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_equal <- function(actual, expected, tolerance = 1e-10, message) {
  ok <- isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = TRUE))
  if (!ok) {
    stop(
      paste0(
        message,
        "\nActual: ", paste(capture.output(str(actual)), collapse = " "),
        "\nExpected: ", paste(capture.output(str(expected)), collapse = " ")
      ),
      call. = FALSE
    )
  }
}

report_pass <- function(label) cat(sprintf("PASS  %s\n", label))

# -------------------------------------------------------------------------
# Test 1: O'Sullivan design dimension and the configurable linear prior.
# -------------------------------------------------------------------------
K <- 4L
K_total <- K + 2L
time_grid <- seq(0, 1, length.out = 9L)
internal_knots <- as.numeric(
  stats::quantile(time_grid, seq(0, 1, length.out = K)[-c(1, K)])
)

X <- X_design(time_grid)
Z_os <- ZOSull(time_grid, range.x = c(0, 1), intKnots = internal_knots)
C_os <- cbind(X, Z_os)

assert_equal(dim(X), c(9L, 2L), message = "X must contain [1,t].")
assert_equal(dim(Z_os), c(9L, K), message = "Z must have K nonlinear columns.")
assert_equal(dim(C_os), c(9L, K_total), message = "C must have K+2 columns.")
assert_true(all(is.finite(C_os)), "The small O'Sullivan design contains non-finite values.")

hyper_default <- set_hyper()
hyper_custom <- set_hyper(sigma_beta = c(2, 3))
assert_equal(
  hyper_default$Sigma_beta,
  1e10 * diag(2),
  tolerance = 1e-12,
  message = "Default Sigma_beta must equal 1e10 I_2."
)
assert_equal(
  hyper_custom$Sigma_beta,
  diag(c(4, 9)),
  message = "Custom linear prior variances must be respected."
)
report_pass("K+2 design and configurable two-dimensional linear prior")

# -------------------------------------------------------------------------
# Test 2: the coefficient update uses blockdiag(Sigma_0^-1, tau I_K).
# Data precision is set to zero so the prior precision is observed directly.
# -------------------------------------------------------------------------
inv_Sigma_0 <- diag(c(0.25, 0.5))
tau_mu <- 4
zero_Ktotal <- matrix(0, nrow = K_total, ncol = K_total)

mu_update <- update_nu_mu(
  Y = list(list(list(0))),
  C = list(list(matrix(0, nrow = 1L, ncol = K_total))),
  list_cp_C = list(list(zero_Ktotal)),
  list_cp_C_Y = list(list(matrix(0, nrow = K_total, ncol = 1L))),
  mu_q_nu_mu = list(),
  Sigma_q_nu_mu = list(),
  sum_list_cp_C = list(zero_Ktotal),
  mu_q_recip_sigsq_eps = matrix(0, nrow = 1L, ncol = 1L),
  mu_q_recip_sigsq_mu = matrix(tau_mu, nrow = 1L, ncol = 1L),
  inv_Sigma_beta = inv_Sigma_0,
  mu_q_nu_beta = NULL,
  Z = NULL,
  mu_q_zeta = NULL,
  mu_q_nu_phi = NULL,
  mu_q_xi = NULL,
  mu_q_nu_psi = NULL,
  mu_q_a = NULL,
  mu_q_b_specific = NULL,
  S = 1L,
  n_s = 1L,
  p = 1L,
  d = 0L,
  L_f = 0L,
  L_s = 0L,
  K = K,
  K_total = K_total,
  c_val = 1,
  n_cpus = 1L
)

expected_precision <- blkdiag(inv_Sigma_0, tau_mu * diag(K))
observed_precision <- mu_update$inv_Sigma_q_nu_mu[[1]][[1]]
assert_equal(
  observed_precision,
  expected_precision,
  message = "Mean coefficient prior precision is not the expected 2+K block diagonal matrix."
)
assert_equal(
  observed_precision[1:2, 1:2],
  inv_Sigma_0,
  message = "The first two precision entries must come from Sigma_0^-1."
)
assert_equal(
  observed_precision[3:K_total, 3:K_total],
  tau_mu * diag(K),
  message = "The last K precision entries must equal tau I_K."
)
report_pass("coefficient update uses the correct 2+K block prior precision")

# -------------------------------------------------------------------------
# Test 3: all four smoothing-variance updates ignore the first two entries
# and respond to changes in the last K entries.
# -------------------------------------------------------------------------
nu_base <- c(10, -7, 0.5, -1, 1.5, 0.25)
Sigma_base <- diag(c(20, 30, 0.10, 0.20, 0.30, 0.40))

nu_linear_changed <- nu_base
nu_linear_changed[1:2] <- c(1e6, -1e6)
Sigma_linear_changed <- Sigma_base
Sigma_linear_changed[1:2, 1:2] <- diag(c(1e8, 2e8))

nu_nonlinear_changed <- nu_base
nu_nonlinear_changed[3] <- nu_nonlinear_changed[3] + 2

aux_precision <- 2

get_variance_updates <- function(nu, Sigma) {
  mu_res <- update_sigsq_mu(
    mu_q_nu_mu = list(list(nu)),
    Sigma_q_nu_mu = list(list(Sigma)),
    mu_q_recip_a_mu = matrix(aux_precision, 1L, 1L),
    S = 1L, p = 1L, K = K, c_val = 1, n_cpus = 1L
  )
  beta_res <- update_sigsq_beta(
    mu_q_nu_beta = list(list(nu)),
    Sigma_q_nu_beta = list(list(Sigma)),
    mu_q_recip_a_beta = matrix(aux_precision, 1L, 1L),
    p = 1L, d = 1L, K = K, c_val = 1, n_cpus = 1L
  )
  phi_res <- update_sigsq_phi(
    mu_q_nu_phi = list(matrix(nu, ncol = 1L)),
    Sigma_q_nu_phi = list(list(Sigma)),
    mu_q_recip_a_phi = list(aux_precision),
    L_f = 1L, M_f = 1L, K = K, c_val = 1
  )
  psi_res <- update_sigsq_psi(
    mu_q_nu_psi = list(list(matrix(nu, ncol = 1L))),
    Sigma_q_nu_psi = list(list(list(Sigma))),
    mu_q_recip_a_psi = list(list(aux_precision)),
    S = 1L, L_s = 1L, M_s = list(1L), K = K, c_val = 1
  )
  c(
    mu = mu_res$lambda_q_sigsq_mu[1, 1],
    beta = beta_res$lambda_q_sigsq_beta[1, 1],
    theta = phi_res$lambda_q_sigsq_phi[[1]][1],
    kappa = psi_res$lambda_q_sigsq_psi[[1]][[1]][1]
  )
}

lambda_base <- get_variance_updates(nu_base, Sigma_base)
lambda_linear_changed <- get_variance_updates(nu_linear_changed, Sigma_linear_changed)
lambda_nonlinear_changed <- get_variance_updates(nu_nonlinear_changed, Sigma_base)

expected_lambda <- aux_precision + 0.5 * (
  sum(nu_base[3:K_total]^2) +
    sum(diag(Sigma_base)[3:K_total])
)

assert_equal(
  unname(lambda_base),
  rep(expected_lambda, 4L),
  message = "At least one smoothing-variance update does not use the expected penalized second moment."
)
assert_equal(
  lambda_linear_changed,
  lambda_base,
  message = "Changing only the two linear coefficients must not change smoothing-variance updates."
)
assert_true(
  all(lambda_nonlinear_changed > lambda_base),
  "Changing a nonlinear coefficient must change all smoothing-variance updates."
)
report_pass("mu, beta, theta and kappa variance updates use only the last K coefficients")

elapsed <- proc.time()[["elapsed"]] - test_started
cat(sprintf(
  "ALL PASS  small sample: K=%d, K_total=%d, n=%d, S=p=d=1; elapsed=%.3f seconds\n",
  K, K_total, length(time_grid), elapsed
))
