suppressMessages({
  library(splines); library(parallel)
  source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
  devtools::load_all(project_root, quiet=TRUE)
})

cat("========================================\n")
cat("  Test B1: Cross-Study Prediction\n")
cat("========================================\n\n")

set.seed(123)

# ---- Generate S=2 data ----
dat <- simulate_multi_study_data(
  S=2, n_s=c(30,30), p=8, d=0,
  L_f=2, L_s=1, M_f=c(3,3), M_s=list(3,3),
  K=7, n_obs=80, sigma_eps=0.1,
  bool_sparse_loadings=TRUE, prop_sparse=0.6
)
tp <- dat$true_params

# ---- Split study 2 ----
n2 <- 30
set.seed(456)
idx_test <- sort(sample(1:n2, floor(n2 * 0.3)))
idx_train <- setdiff(1:n2, idx_test)
cat(sprintf("Study 2 split: %d train, %d test\n", length(idx_train), length(idx_test)))

Y_train <- dat$Y; Y_train[[2]] <- Y_train[[2]][idx_train]
time_obs_train <- dat$time_obs; time_obs_train[[2]] <- time_obs_train[[2]][idx_train]

# ---- Fit ----
fit <- bayesSYNC_multi(
  Y=Y_train, time_obs=time_obs_train,
  L_f=2, L_s=1, M_f=c(3,3), M_s=list(3,3),
  K=7, anneal=NULL, maxit=80, verbose=FALSE, seed=123
)
cat(sprintf("Fit done. ELBO=%.1f, %d iters\n", tail(fit$ELBO,1), fit$i_iter))

# ---- Predict: fold-in test individuals ----
# For each test individual, run 5 CAVI iterations updating only zeta/xi
fold_in_predict <- function(fit, C, Y, s, i_test, n_iter=5) {
  p <- fit$p; L_f <- fit$L_f; L_s <- fit$L_s
  M_f <- fit$M_f; M_s <- fit$M_s
  n_obs <- length(Y[[s]][[i_test]][[1]])

  C_si <- C[[s]][[i_test]]
  cp_C_si <- crossprod(C_si)

  # Initialize zeta/xi at 0 (prior mean)
  zeta <- vector("list", L_f)
  for (l in seq_len(L_f)) {
    zeta[[l]] <- rep(0, M_f[l])
  }
  xi <- NULL
  if (L_s > 0) {
    xi <- vector("list", L_s)
    for (l in seq_len(L_s)) {
      xi[[l]] <- rep(0, M_s[[s]][l])
    }
  }

  for (iter in 1:n_iter) {
    # ---- Update zeta (shared) ----
    for (l in seq_len(L_f)) {
      Phi <- C_si %*% fit$mu_q_nu_phi[[l]]  # n_obs x M_f[l]
      H_phi <- crossprod(Phi)
      sum_sigma_a <- sum(fit$mu_q_recip_sigsq_eps[s,] * fit$term_a[,l])

      # Residual: remove mu, other shared factors, all specific factors
      res <- rep(0, ncol(C_si))
      for (jj in 1:p) {
        # y - mu
        r_jj <- Y[[s]][[i_test]][[jj]] - as.vector(C_si %*% fit$mu_q_nu_mu[[s]][[jj]])
        ct_r <- crossprod(C_si, r_jj)
        # Remove other shared factors
        for (l_other in setdiff(seq_len(L_f), l)) {
          f_other <- as.vector(C_si %*% fit$mu_q_nu_phi[[l_other]] %*% zeta[[l_other]])
          ct_r <- ct_r - fit$mu_q_a[jj, l_other] * crossprod(C_si, f_other)
        }
        # Remove all specific factors
        if (L_s > 0) {
          for (l_spec in seq_len(L_s)) {
            g_spec <- as.vector(C_si %*% fit$mu_q_nu_psi[[s]][[l_spec]] %*% xi[[l_spec]])
            ct_r <- ct_r - fit$mu_q_b_specific[[s]][jj, l_spec] * crossprod(C_si, g_spec)
          }
        }
        res <- res + fit$mu_q_recip_sigsq_eps[s,jj] * fit$mu_q_a[jj,l] * as.vector(ct_r)
      }

      prec <- sum_sigma_a * H_phi + diag(M_f[l])
      zeta[[l]] <- as.vector(solve(prec, crossprod(fit$mu_q_nu_phi[[l]], res)))
    }

    # ---- Update xi (specific) ----
    if (L_s > 0) {
      for (l in seq_len(L_s)) {
        Psi <- C_si %*% fit$mu_q_nu_psi[[s]][[l]]
        H_psi <- crossprod(Psi)
        sum_sigma_b <- sum(fit$mu_q_recip_sigsq_eps[s,] * fit$term_b_specific[[s]][,l])

        res <- rep(0, ncol(C_si))
        for (jj in 1:p) {
          r_jj <- Y[[s]][[i_test]][[jj]] - as.vector(C_si %*% fit$mu_q_nu_mu[[s]][[jj]])
          ct_r <- crossprod(C_si, r_jj)
          for (l_shared in seq_len(L_f)) {
            f_shared <- as.vector(C_si %*% fit$mu_q_nu_phi[[l_shared]] %*% zeta[[l_shared]])
            ct_r <- ct_r - fit$mu_q_a[jj, l_shared] * crossprod(C_si, f_shared)
          }
          for (l_other in setdiff(seq_len(L_s), l)) {
            g_other <- as.vector(C_si %*% fit$mu_q_nu_psi[[s]][[l_other]] %*% xi[[l_other]])
            ct_r <- ct_r - fit$mu_q_b_specific[[s]][jj, l_other] * crossprod(C_si, g_other)
          }
          res <- res + fit$mu_q_recip_sigsq_eps[s,jj] * fit$mu_q_b_specific[[s]][jj,l] * as.vector(ct_r)
        }

        prec <- sum_sigma_b * H_psi + diag(M_s[[s]][l])
        xi[[l]] <- as.vector(solve(prec, crossprod(fit$mu_q_nu_psi[[s]][[l]], res)))
      }
    }
  }

  # ---- Reconstruct prediction ----
  y_hat <- lapply(1:p, function(j) as.vector(C_si %*% fit$mu_q_nu_mu[[s]][[j]]))
  for (l in seq_len(L_f)) {
    f_val <- as.vector(C_si %*% fit$mu_q_nu_phi[[l]] %*% zeta[[l]])
    for (j in 1:p) y_hat[[j]] <- y_hat[[j]] + fit$mu_q_a[j,l] * f_val
  }
  if (L_s > 0) {
    for (l in seq_len(L_s)) {
      g_val <- as.vector(C_si %*% fit$mu_q_nu_psi[[s]][[l]] %*% xi[[l]])
      for (j in 1:p) y_hat[[j]] <- y_hat[[j]] + fit$mu_q_b_specific[[s]][j,l] * g_val
    }
  }

  # Compute per-variable RMSE
  sapply(1:p, function(j) sqrt(mean((y_hat[[j]] - Y[[s]][[i_test]][[j]])^2)))
}

# ---- Run prediction ----
rmse_test <- c()
rmse_baseline <- c()

for (i in idx_test) {
  rmse_vars <- fold_in_predict(fit, dat$C, dat$Y, s=2, i, n_iter=10)
  rmse_test <- c(rmse_test, rmse_vars)
  # Baseline: mu only (no factors)
  for (j in 1:fit$p) {
    mu_j <- as.vector(dat$C[[2]][[i]] %*% fit$mu_q_nu_mu[[2]][[j]])
    rmse_baseline <- c(rmse_baseline, sqrt(mean((mu_j - dat$Y[[2]][[i]][[j]])^2)))
  }
}

cat(sprintf("\nResults (9 test individuals x 8 vars = 72 predictions):\n"))
cat(sprintf("  Full model RMSE:  %.4f\n", mean(rmse_test)))
cat(sprintf("  Baseline RMSE:    %.4f  (mu only, no factors)\n", mean(rmse_baseline)))
cat(sprintf("  Improvement:      %.1f%%\n", (1 - mean(rmse_test)/mean(rmse_baseline))*100))

# Also check: what about study 1 individuals (training)? Their predictions should be better
cat("\nStudy 1 (training) reconstruction:\n")
rmse_s1 <- c()
rmse_s1_base <- c()
for (i in 1:dat$true_params$n_s[1]) {
  rmse_v <- fold_in_predict(fit, dat$C, dat$Y, s=1, i, n_iter=10)
  rmse_s1 <- c(rmse_s1, rmse_v)
  for (j in 1:fit$p) {
    mu_j <- as.vector(dat$C[[1]][[i]] %*% fit$mu_q_nu_mu[[1]][[j]])
    rmse_s1_base <- c(rmse_s1_base, sqrt(mean((mu_j - dat$Y[[1]][[i]][[j]])^2)))
  }
}
cat(sprintf("  Full model RMSE:  %.4f\n", mean(rmse_s1)))
cat(sprintf("  Baseline RMSE:    %.4f\n", mean(rmse_s1_base)))
cat(sprintf("  Improvement:      %.1f%%\n", (1 - mean(rmse_s1)/mean(rmse_s1_base))*100))

cat("\n========================================\n")
cat("  Test B1: PASS\n")
cat("========================================\n")
