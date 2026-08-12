library(splines)
library(parallel)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
src_dir <- file.path(project_root, "R")
for (f in list.files(src_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

cat("=== Phase 4 Tests: Mean Function Update ===\n\n")

# ---- Setup: small data ----
set.seed(42)
dat <- simulate_multi_study_data(S=1, n_s=c(20), p=3, d=0, L_f=0, L_s=0,
  K=6, n_obs=30, sigma_eps=0.05)
tp <- dat$true_params

cat("Test 1: update_mu.R compiles ... ")
# source already done, verify function exists
stopifnot(exists("update_nu_mu", mode="function"))
cat("PASS\n")

# ---- Fit with 1 iteration to get initial q(nu_mu) ----
fit <- bayesSYNC_multi(Y=dat$Y, time_obs=dat$time_obs,
  L_f=0, L_s=0, M_f=integer(0), M_s=list(integer(0)),
  K=tp$K, maxit=20, verbose=FALSE, bool_scale=FALSE)

K_total <- fit$K + 2

# ---- Test 2: Sigma_q_nu_mu positive definite ----
cat("Test 2: Sigma_q_nu_mu positive definite ... ")
for (s in 1:1) {
  for (j in 1:3) {
    eig <- eigen(fit$Sigma_q_nu_mu[[s]][[j]], symmetric=TRUE)
    stopifnot(all(eig$values > 0))
  }
}
cat("PASS\n")

# ---- Test 3: mu_q_nu_mu dimension ----
cat("Test 3: mu_q_nu_mu dimension ... ")
for (s in 1:1) {
  for (j in 1:3) {
    stopifnot(length(fit$mu_q_nu_mu[[s]][[j]]) == K_total)
  }
}
cat(sprintf("PASS  (length = %d = K+2)\n", K_total))

# ---- Test 4: Updated mu values are finite and reasonable ----
cat("Test 4: Updated mu values ... ")
mu_vals <- unlist(lapply(fit$mu_q_nu_mu[[1]], function(v) v))
stopifnot(all(is.finite(mu_vals)))
stopifnot(max(abs(mu_vals)) < 100)  # not exploding
cat(sprintf("PASS  (range %.3f - %.3f)\n", min(mu_vals), max(mu_vals)))

# ---- Test 5: ELBO after fit is finite ----
cat("Test 5: ELBO finite ... ")
stopifnot(is.finite(tail(fit$ELBO, 1)))
cat(sprintf("PASS  (ELBO = %.1f)\n", tail(fit$ELBO, 1)))

# ---- Test 6: No NaN/Inf in any mu params ----
cat("Test 6: No NaN/Inf in mu params ... ")
all_mu <- unlist(fit$mu_q_nu_mu)
stopifnot(all(is.finite(all_mu)))
all_sigma <- unlist(lapply(fit$Sigma_q_nu_mu[[1]], diag))
stopifnot(all(is.finite(all_sigma)))
stopifnot(all(all_sigma > 0))
cat("PASS\n")

cat("\n=== ALL Phase 4 TESTS PASSED ===\n")
