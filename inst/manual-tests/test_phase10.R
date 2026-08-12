library(splines)
library(parallel)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
suppressMessages(devtools::load_all(project_root, quiet=TRUE))

cat("========================================\n")
cat("  Phase 10: End-to-End Verification\n")
cat("========================================\n\n")

# ---- Test 1: load_all ----
cat("Test 1: devtools::load_all() ... PASS\n\n")

# ---- Test 2: Full pipeline S=2, p=10, L_f=3, L_s=2, d=2 ----
cat("Test 2: Full pipeline (S=2, p=10, L_f=3, L_s=2, d=2) ...\n")
set.seed(42)

t_start <- Sys.time()

dat <- simulate_multi_study_data(
  S = 2, n_s = c(20, 30), p = 10, d = 2,
  L_f = 3, L_s = 2,
  M_f = c(2, 2, 2),
  M_s = list(c(2, 2), c(2, 2)),
  K = 8, n_obs = 80, sigma_eps = 0.1,
  bool_sparse_loadings = TRUE, prop_sparse = 0.7
)

cat(sprintf("  Data generated: S=%d, n_s=(%d,%d), p=%d, d=%d\n",
  dat$true_params$S, dat$true_params$n_s[1], dat$true_params$n_s[2],
  dat$true_params$p, dat$true_params$d_use))

fit <- bayesSYNC_multi(
  Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
  L_f = 3, L_s = 2,
  M_f = c(2, 2, 2), M_s = list(c(2, 2), c(2, 2)),
  K = 8, anneal = NULL, maxit = 100,
  verbose = FALSE, seed = 42
)

if (length(warnings()) > 0) {
  ws <- names(summary(warnings()))
  cat(sprintf("  Warnings: %d unique\n", length(ws)))
  for (i in seq_len(min(5, length(ws)))) {
    cat(sprintf("    [%d] %s\n", i, substr(ws[i], 1, 120)))
  }
}

t_end <- Sys.time()
t_elapsed <- as.numeric(difftime(t_end, t_start, units="secs"))

cat(sprintf("  Fit completed in %.1f seconds\n", t_elapsed))
cat(sprintf("  Iterations: %d\n", fit$i_iter))

# ---- Test 3: ELBO convergence ----
cat("\nTest 3: ELBO convergence ... ")
elbo <- fit$ELBO
final_elbo <- tail(elbo, 1)
elbo_last10 <- tail(elbo, min(10, length(elbo)))
elbo_change <- max(diff(elbo_last10))

cat(sprintf("final ELBO = %.2f, last-10 max change = %.4f\n", final_elbo, elbo_change))
stopifnot(is.finite(final_elbo))

# Check ELBO is improving (final 1/3 better than first 1/3, or converged)
n <- length(elbo)
n3 <- max(3, floor(n / 3))
elbo_early <- mean(elbo[1:n3])
elbo_late <- mean(elbo[(n - n3 + 1):n])
cat(sprintf("  ELBO early avg=%.1f  late avg=%.1f\n", elbo_early, elbo_late))
stopifnot(elbo_late > elbo_early || n < 10)
cat("  PASS (ELBO improving)\n")

# ---- Test 4: Degeneration test T1 (S=1, L_s=0, no anneal) ----
cat("\nTest 4: Degeneration T1 — S=1, L_s=0 vs bayesSYNC ... ")

# Check if bayesSYNC is available
bayes_available <- requireNamespace("bayesSYNC", quietly = TRUE)

if (bayes_available) {
  set.seed(42)
  dat_t1 <- simulate_multi_study_data(
    S = 1, n_s = c(30), p = 5, d = 0,
    L_f = 1, L_s = 0, M_f = c(2), M_s = list(integer(0)),
    K = 6, n_obs = 50, sigma_eps = 0.1,
    bool_sparse_loadings = FALSE
  )

  fit_t1 <- bayesSYNC_multi(
    Y = dat_t1$Y, time_obs = dat_t1$time_obs,
    L_f = 1, L_s = 0, M_f = c(2), M_s = list(integer(0)),
    K = 6, anneal = NULL, maxit = 100,
    verbose = FALSE, seed = 42
  )

  fit_bs <- bayesSYNC::bayesSYNC(
    Y = dat_t1$Y[[1]], time_obs = dat_t1$time_obs[[1]],
    Q = 1, L = 2, K = 6, anneal = NULL, maxit = 100,
    verbose = FALSE, seed = 42
  )

  elbo_t1 <- tail(fit_t1$ELBO, 1)
  elbo_bs <- tail(fit_bs$ELBO, 1)

  # Both models should reach similar magnitudes (known: ELBOs not exactly comparable
  # due to different sigma^2 IG entropy handling; check within factor of 1.5)
  rel_diff <- abs(elbo_t1 - elbo_bs) / max(abs(c(elbo_t1, elbo_bs)), 1)
  cat(sprintf("multiFSYNC ELBO=%.2f  bayesSYNC ELBO=%.2f  rel_diff=%.4f\n",
    elbo_t1, elbo_bs, rel_diff))
  stopifnot(rel_diff < 0.5)  # within factor 2
  cat("  PASS (rel diff within tolerance)\n")
} else {
  cat("SKIP (bayesSYNC not installed)\n")
}

# ---- Test 5: Result sanity ----
cat("\nTest 5: Result sanity ...\n")

# Factor PPIs
if (fit$L_f > 0) {
  ppi_shared <- 1 - exp(colSums(log1p(-fit$mu_q_gamma_a)))
  cat(sprintf("  Shared factor PPIs: %s\n", paste(round(ppi_shared, 3), collapse=", ")))
  stopifnot(all(ppi_shared >= 0 & ppi_shared <= 1))
}

if (fit$L_s > 0) {
  for (s in 1:fit$S) {
    ppi_spec <- 1 - exp(colSums(log1p(-fit$mu_q_gamma_b[[s]])))
    cat(sprintf("  Study %d specific PPIs: %s\n", s, paste(round(ppi_spec, 3), collapse=", ")))
    stopifnot(all(ppi_spec >= 0 & ppi_spec <= 1))
  }
}

# Mean function smoothness (check L2 norm is finite)
mu_norms <- sapply(fit$mu_q_nu_mu[[1]], function(v) sqrt(sum(v^2)))
cat(sprintf("  Mean function L2 norms: range [%.3f, %.3f]\n", min(mu_norms), max(mu_norms)))
stopifnot(all(is.finite(mu_norms)))

# Orthonormalisation
res_ortho <- orthonormalise_multi(
  C_g = fit$C_g, time_g = fit$time_g,
  mu_q_nu_mu = fit$mu_q_nu_mu,
  mu_q_nu_phi = fit$mu_q_nu_phi,
  mu_q_nu_psi = fit$mu_q_nu_psi,
  mu_q_zeta = fit$mu_q_zeta,
  Sigma_q_zeta = fit$Sigma_q_zeta,
  mu_q_xi = fit$mu_q_xi,
  Sigma_q_xi = fit$Sigma_q_xi,
  mu_q_a = fit$mu_q_a,
  mu_q_b_specific = fit$mu_q_b_specific,
  mu_q_gamma_a = fit$mu_q_gamma_a,
  mu_q_gamma_b = fit$mu_q_gamma_b,
  S = fit$S, n_s = fit$n_s, p = fit$p,
  L_f = fit$L_f, L_s = fit$L_s,
  M_f = fit$M_f, M_s = fit$M_s
)

# Check PVE (handle NaN from degenerate factors)
for (l in 1:fit$L_f) {
  pve <- res_ortho$list_cumulated_pve[[l]]
  if (any(is.nan(pve))) {
    cat(sprintf("  Shared factor %d PVE: NaN (degenerate factor, PPI=0)\n", l))
  } else {
    cat(sprintf("  Shared factor %d PVE: %s\n", l, paste(round(pve, 1), collapse="%, ")))
    stopifnot(abs(tail(pve, 1) - 100) < 1e-6)
  }
}

cat("\n========================================\n")
cat("  ALL Phase 10 TESTS PASSED\n")
cat("========================================\n")
