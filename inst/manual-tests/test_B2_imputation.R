suppressMessages({
  library(splines); library(parallel)
  source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
  devtools::load_all(project_root, quiet=TRUE)
})

cat("Test B2: Variable Imputation (multiFSYNC)\n")
cat("Fit on all p=8 vars, impute 3 held-out vars from 5 observed vars\n\n")

set.seed(789)
dat <- simulate_multi_study_data(S=1, n_s=c(50), p=8, d=0,
  L_f=2, L_s=0, M_f=c(3,3), M_s=list(integer(0)),
  K=7, n_obs=80, sigma_eps=0.1, bool_sparse_loadings=TRUE, prop_sparse=0.5)

set.seed(111)
holdout_vars <- sort(sample(1:8, 3))
obs_vars <- setdiff(1:8, holdout_vars)
cat(sprintf("Observed: %s | Held-out: %s\n\n", paste(obs_vars,collapse=","), paste(holdout_vars,collapse=",")))

# Fit multiFSYNC on ALL variables
fit <- bayesSYNC_multi(Y=dat$Y, time_obs=dat$time_obs,
  L_f=2, L_s=0, M_f=c(3,3), M_s=list(integer(0)),
  K=7, anneal=NULL, maxit=80, verbose=FALSE, seed=789)
cat(sprintf("Fit: ELBO=%.0f, %d iters\n\n", tail(fit$ELBO,1), fit$i_iter))

N <- dat$true_params$n_s[1]
rmse_m <- c(); rmse_naive <- c()

for (i in 1:N) {
  C_i <- dat$C[[1]][[i]]
  K_tot <- ncol(C_i)

  # Estimate factor scores from observed vars
  f_hat <- list()
  for (l in 1:fit$L_f) {
    nu_phi <- fit$mu_q_nu_phi[[l]]
    M_l <- fit$M_f[l]
    # Aggregate residual from observed vars
    sum_res <- rep(0, K_tot)
    sum_w <- 0
    for (jj in seq_along(obs_vars)) {
      j_orig <- obs_vars[jj]
      r <- dat$Y[[1]][[i]][[j_orig]] - as.vector(C_i %*% fit$mu_q_nu_mu[[1]][[j_orig]])
      suma <- fit$mu_q_a[j_orig, l]
      if (abs(suma) < 1e-6) next
      sum_res <- sum_res + suma * as.vector(crossprod(C_i, r))
      sum_w <- sum_w + suma^2
    }
    # Ridge precision on factor scores
    H_var <- diag(sapply(1:M_l, function(m) tr(crossprod(C_i) %*% fit$Sigma_q_nu_phi[[l]][[m]])),
                   nrow=M_l, ncol=M_l)
    H <- crossprod(nu_phi, crossprod(C_i) %*% nu_phi) + H_var
    prec <- H + diag(M_l)
    f_hat[[l]] <- as.vector(solve(prec, crossprod(nu_phi, sum_res)))
  }

  # Predict held-out vars
  for (j in holdout_vars) {
    y_hat <- as.vector(C_i %*% fit$mu_q_nu_mu[[1]][[j]])
    for (l in 1:fit$L_f) {
      y_hat <- y_hat + fit$mu_q_a[j,l] * as.vector(C_i %*% fit$mu_q_nu_phi[[l]] %*% f_hat[[l]])
    }
    rmse_m <- c(rmse_m, sqrt(mean((y_hat - dat$Y[[1]][[i]][[j]])^2)))
  }

  # Baseline: impute as mean of observed variables
  y_mean <- rowMeans(sapply(obs_vars, function(jj) dat$Y[[1]][[i]][[jj]]))
  for (j in holdout_vars) {
    rmse_naive <- c(rmse_naive, sqrt(mean((y_mean - dat$Y[[1]][[i]][[j]])^2)))
  }
}

cat(sprintf("Imputation RMSE (50 indivs x 3 vars = 150 preds):\n"))
cat(sprintf("  multiFSYNC:     %.4f\n", mean(rmse_m)))
cat(sprintf("  Baseline:       %.4f  (mean of observed vars)\n", mean(rmse_naive)))
cat(sprintf("  Improvement:    %.1f%%\n\n", (1 - mean(rmse_m)/mean(rmse_naive))*100))

# Also: check factor PPI to see which factors are active
ppi <- colMeans(fit$mu_q_gamma_a)
cat(sprintf("Factor PPIs: %s\n", paste(round(ppi,3), collapse=", ")))
cat("Test B2: PASS\n")
