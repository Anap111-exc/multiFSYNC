# Step 10: End-to-end verification after full symbol rename
library(splines)
library(parallel)

# Source all R files
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
src_dir <- file.path(project_root, "R")
r_files <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

cat("=== All source files loaded ===\n\n")

# Test 1: Simple S=1, L_f=1, L_s=0
cat("Test 1: S=1, L_f=1, L_s=0, M_f=c(2)\n")
set.seed(42)
dat <- simulate_multi_study_data(
  S = 1, n_s = c(20), p = 3, d = 0,
  L_f = 1, L_s = 0,
  M_f = c(2), n_obs = 30, sigma_eps = 0.05
)
cat("  Data generated OK\n")

fit <- bayesSYNC_multi(
  Y = dat$Y,
  time_obs = dat$time_obs,
  L_f = 1, L_s = 0, M_f = c(2), M_s = list(integer(0)),
  K = 8, maxit = 3, verbose = FALSE
)
cat("  Fit completed. Final ELBO:", tail(fit$ELBO, 1), "\n\n")

# Test 2: S=2, L_f=2, L_s=1, with covariates
cat("Test 2: S=2, L_f=2, L_s=1, M_f=c(2,1), M_s=list(c(2), c(2))\n")
set.seed(42)
dat2 <- simulate_multi_study_data(
  S = 2, n_s = c(15, 20), p = 4, d = 2,
  L_f = 2, L_s = 1,
  M_f = c(2, 1), M_s = list(c(2), c(2)),
  n_obs = 25, sigma_eps = 0.05
)
cat("  Data generated OK\n")

fit2 <- bayesSYNC_multi(
  Y = dat2$Y, Z = dat2$Z,
  time_obs = dat2$time_obs,
  L_f = 2, L_s = 1,
  M_f = c(2, 1), M_s = list(c(2), c(2)),
  K = 6, maxit = 3, verbose = FALSE
)
cat("  Fit completed. Final ELBO:", tail(fit2$ELBO, 1), "\n\n")

# Test 3: L_f=0, L_s=0
cat("Test 3: L_f=0, L_s=0\n")
dat0 <- simulate_multi_study_data(
  S = 1, n_s = c(10), p = 3, d = 0,
  L_f = 0, L_s = 0, M_f = integer(0),
  M_s = list(integer(0)), n_obs = 20, sigma_eps = 0.05
)
cat("  Data generated OK\n")

fit0 <- bayesSYNC_multi(
  Y = dat0$Y,
  time_obs = dat0$time_obs,
  L_f = 0, L_s = 0,
  M_f = integer(0), M_s = list(integer(0)),
  K = 6, maxit = 3, verbose = FALSE
)
cat("  Fit completed. Final ELBO:", tail(fit0$ELBO, 1), "\n\n")

cat("=== ALL E2E TESTS PASSED ===\n")
