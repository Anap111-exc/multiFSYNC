# Phase B: Recovery Accuracy Tests
# Measures multiFSYNC recovery of known ground truth under various scenarios.
# Factor models have sign/rotational ambiguity; focus on:
#   - correct dimensions and output structure
#   - PPI selection behavior (not exact loadings)
#   - orthonormalisation validity
#   - signal reconstruction quality
#   - mean function recovery (identifiable)

library(splines)
library(pracma)

# ---- Helper: sign-adjusted column matching ----
safe_cor <- function(x, y) {
  sx <- sd(x, na.rm=TRUE); sy <- sd(y, na.rm=TRUE)
  if (is.na(sx) || is.na(sy) || sx < 1e-10 || sy < 1e-10) return(NA)
  cxy <- cor(x, y)
  if (is.na(cxy)) return(NA) else cxy
}

match_columns <- function(X_est, X_true) {
  k_est <- ncol(X_est); k_true <- ncol(X_true)
  if (k_est == 0 || k_true == 0) return(list(idx = rep(1, k_est), flip = rep(1, k_est)))
  matched_idx <- rep(NA, k_est); sign_flip <- rep(1, k_est)
  available <- rep(TRUE, k_true)
  for (i in seq_len(k_est)) {
    cors <- rep(-1, k_true)
    for (j in which(available)) cors[j] <- abs(safe_cor(X_est[, i], X_true[, j]))
    best <- which.max(cors)
    if (length(best) == 0 || best < 1 || best > k_true) { matched_idx[i] <- 1; next }
    matched_idx[i] <- best; available[best] <- FALSE
    sf <- sign(safe_cor(X_est[, i], X_true[, best]))
    sign_flip[i] <- if (is.na(sf) || sf == 0) 1 else sf
  }
  list(idx = matched_idx, flip = sign_flip)
}

selection_metrics <- function(est_prob, truth_binary) {
  pred <- as.numeric(est_prob > 0.5)
  tp <- sum(pred == 1 & truth_binary == 1)
  fp <- sum(pred == 1 & truth_binary == 0)
  tn <- sum(pred == 0 & truth_binary == 0)
  fn <- sum(pred == 0 & truth_binary == 1)
  tpr <- if (tp + fn > 0) tp / (tp + fn) else 1
  fpr <- if (fp + tn > 0) fp / (fp + tn) else 0
  list(TPR = tpr, FPR = fpr, N_active = sum(truth_binary), N_est_active = sum(pred))
}

# ---- Scenario 1: Degenerate S=1, L_s=0 ----
test_that("S1: S=1 L_s=0 no sparsity — correct structure & PPI valid", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(30), p=6, d=0,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, n_obs=80, common_grid=TRUE,
    sigma_eps=0.05, bool_sparse_loadings=FALSE, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_equal(dim(fit$mu_q_a), c(6, 2))
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))
  expect_true(is.numeric(fit$ELBO_iter) && fit$ELBO_iter > -1e6)

  # Loading correlation (robust to sign ambiguity)
  mm <- match_columns(fit$mu_q_a, dat$true_params$a_true)
  cors <- sapply(seq_len(fit$L_f), function(i) safe_cor(fit$mu_q_a[, i] * mm$flip[i], dat$true_params$a_true[, mm$idx[i]]))
  cat(sprintf("  S1: a cor = %s\n", paste(sprintf("%.3f", cors), collapse=", ")))
})

test_that("S2: S=1 L_s=0 sparse — PPI & MSE report", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(40), p=10, d=0,
    L_f=3, L_s=0, M_f=c(1,1,1), M_s=list(integer(0)), K=8, n_obs=100, common_grid=TRUE,
    sigma_eps=0.05, bool_sparse_loadings=TRUE, prop_sparse=0.7, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=3, L_s=0, M_f=c(1,1,1), M_s=list(integer(0)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_equal(dim(fit$mu_q_a), c(10, 3))
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))

  mm <- match_columns(fit$mu_q_a, dat$true_params$a_true)
  cors <- sapply(seq_len(fit$L_f), function(i) safe_cor(fit$mu_q_a[, i] * mm$flip[i], dat$true_params$a_true[, mm$idx[i]]))
  sel <- selection_metrics(as.vector(fit$mu_q_gamma_a), as.vector(dat$true_params$gamma_a_true))

  cat(sprintf("  S2: a cor = %s\n", paste(sprintf("%.3f", cors), collapse=", ")))
  cat(sprintf("  S2: PPI TPR=%.3f FPR=%.3f (active=%d est=%d)\n",
      sel$TPR, sel$FPR, sel$N_active, sel$N_est_active))
})

# ---- Scenario 3: Multi-study S=2, L_s=1, sparse ----
test_that("S3: S=2 L_f=2 L_s=1 sparse — structure & PPI", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=2, n_s=c(30, 30), p=10, d=0,
    L_f=2, L_s=1, M_f=c(1,1), M_s=list(c(1), c(1)), K=8, n_obs=80, common_grid=TRUE,
    sigma_eps=0.05, bool_sparse_loadings=TRUE, prop_sparse=0.7, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=1, M_f=c(1,1), M_s=list(c(1), c(1)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_equal(dim(fit$mu_q_a), c(10, 2))
  expect_true(all(fit$mu_q_gamma_a >= 0 & fit$mu_q_gamma_a <= 1))

  # Shared PPI
  sel <- selection_metrics(as.vector(fit$mu_q_gamma_a), as.vector(dat$true_params$gamma_a_true))
  cat(sprintf("  S3 shared: PPI TPR=%.3f FPR=%.3f\n", sel$TPR, sel$FPR))

  # Specific loadings exist
  for (s in 1:2) {
    expect_true(!is.null(fit$mu_q_b_specific[[s]]))
    expect_equal(ncol(as.matrix(fit$mu_q_b_specific[[s]])), 1)
  }

  # ELBO finite
  expect_true(is.numeric(fit$ELBO_iter) && fit$ELBO_iter > -1e6)
})

# ---- Scenario 4: Multi-study with covariates ----
test_that("S4: S=2 d=1 covariates — beta structure", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=2, n_s=c(15, 15), p=6, d=1,
    L_f=2, L_s=1, M_f=c(1,1), M_s=list(c(1), c(1)), K=8, n_obs=50, common_grid=TRUE,
    sigma_eps=0.1, bool_sparse_loadings=TRUE, prop_sparse=0.7, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=dat$Z, time_obs=dat$time_obs,
    L_f=2, L_s=1, M_f=c(1,1), M_s=list(c(1), c(1)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_true(!is.null(fit$mu_q_nu_beta))
  expect_equal(length(fit$mu_q_nu_beta), 6)       # p=6
  expect_equal(length(fit$mu_q_nu_beta[[1]]), 1)   # d=1
  expect_true(is.numeric(fit$ELBO_iter) && fit$ELBO_iter > -1e6)
  cat("  S4: covariates OK\n")
})

# ---- Scenario 5: Mean function recovery ----
test_that("S5: Mean function recovery (identifiable)", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(30), p=4, d=0,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=8, n_obs=100, common_grid=TRUE,
    sigma_eps=0.01, bool_sparse_loadings=FALSE, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=1, L_s=0, M_f=c(1), M_s=list(integer(0)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  mse_mu <- mean(sapply(1:4, function(j) mean((fit$mu_q_nu_mu[[1]][[j]] - dat$true_params$nu_mu_true[[1]][[j]])^2)))
  cat(sprintf("  S5: nu_mu MSE = %.6f\n", mse_mu))
  expect_true(mse_mu < 0.5, label = "Mean function MSE < 0.5")
})

# ---- Scenario 6: ELBO structure ----
test_that("S6: ELBO finite & PPI valid", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(30), p=6, d=0,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, n_obs=80, common_grid=TRUE,
    sigma_eps=0.05, bool_sparse_loadings=FALSE, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  expect_true(length(fit$ELBO) >= 5)
  expect_true(is.finite(tail(fit$ELBO, 1)))
  cat(sprintf("  S6: ELBO len=%d final=%.2f\n", length(fit$ELBO), tail(fit$ELBO, 1)))
})

# ---- Scenario 7: Orthonormalised eigenfunctions ----
test_that("S7: Orthonormal eigenfunctions & PVE", {
  set.seed(42)
  dat <- simulate_multi_study_data(S=1, n_s=c(40), p=6, d=0,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, n_obs=80, common_grid=TRUE,
    sigma_eps=0.02, bool_sparse_loadings=FALSE, seed=42)
  fit <- bayesSYNC_multi(Y=dat$Y, Z=NULL, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(1,1), M_s=list(integer(0)), K=8, anneal=NULL, maxit=80,
    n_cpus=1, verbose=FALSE, seed=42, bool_scale=FALSE)

  orth <- orthonormalise_multi(C_g=fit$C_g, time_g=fit$time_g,
    mu_q_nu_mu=fit$mu_q_nu_mu, mu_q_nu_phi=fit$mu_q_nu_phi,
    mu_q_nu_psi=fit$mu_q_nu_psi,
    mu_q_zeta=fit$mu_q_zeta, Sigma_q_zeta=fit$Sigma_q_zeta,
    mu_q_xi=fit$mu_q_xi, Sigma_q_xi=fit$Sigma_q_xi,
    mu_q_a=fit$mu_q_a, mu_q_b_specific=fit$mu_q_b_specific,
    mu_q_gamma_a=fit$mu_q_gamma_a, mu_q_gamma_b=fit$mu_q_gamma_b,
    S=fit$S, n_s=fit$n_s, p=fit$p,
    L_f=fit$L_f, L_s=fit$L_s, M_f=fit$M_f, M_s=fit$M_s)

  # PVE sums to 100
  for (l in seq_len(fit$L_f)) {
    pve <- orth$list_cumulated_pve[[l]]
    expect_true(tail(pve, 1) > 99, label = sprintf("PVE factor %d > 99", l))
  }

  # Eigenfunction correlation with truth
  C_g <- fit$C_g
  for (l in seq_len(fit$L_f)) {
    P_true <- C_g %*% dat$true_params$nu_phi_true[[l]]
    P_hat <- orth$list_Phi_hat[[l]]
    cors <- sapply(seq_len(ncol(P_hat)), function(m) safe_cor(P_hat[, m], P_true[, m]))
    cat(sprintf("  S7 factor %d: eigenfunction cor = %s\n", l,
        paste(sprintf("%.3f", cors), collapse=", ")))
  }

  # Factor PPI
  expect_true(all(orth$factor_ppi_shared >= 0 & orth$factor_ppi_shared <= 1))
  cat(sprintf("  S7: factor PPI = %s\n", paste(sprintf("%.3f", orth$factor_ppi_shared), collapse=", ")))
})
