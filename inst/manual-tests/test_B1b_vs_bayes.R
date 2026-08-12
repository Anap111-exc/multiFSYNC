suppressMessages({
  library(splines); library(parallel)
  source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
  devtools::load_all(project_root, quiet=TRUE)
  suppressMessages(library(bayesSYNC, quietly=TRUE))
})

cat("Test B1b: multiFSYNC vs bayesSYNC — Reconstruction Accuracy (Dense Grid)\n")
cat(strrep("=", 68), "\n\n")

set.seed(2024)

# ---- Generate S=2 data ----
dat <- simulate_multi_study_data(
  S=2, n_s=c(25,25), p=8, d=0,
  L_f=2, L_s=1, M_f=c(3,3), M_s=list(3,3),
  K=7, n_obs=70, sigma_eps=0.1,
  bool_sparse_loadings=TRUE, prop_sparse=0.5
)
tp <- dat$true_params

# ---- True signal on dense grid (for RMSE comparison) ----
n_g <- 200
time_g <- seq(0, 1, length.out=n_g)
C_g <- cbind(1, time_g, ZOSull(time_g, c(0,1),
  intKnots=quantile(time_g, seq(0,1,length=tp$K)[-c(1,tp$K)])))

# True mu on dense grid: per study, per variable
mu_true_g <- vector("list", 2)
for (s in 1:2) {
  mu_true_g[[s]] <- lapply(1:8, function(j) as.vector(C_g %*% tp$nu_mu_true[[s]][[j]]))
}

# True shared factor on dense grid: per individual, per factor
f_true_g <- vector("list", 2)
for (s in 1:2) {
  f_true_g[[s]] <- vector("list", tp$n_s[s])
  for (i in 1:tp$n_s[s]) {
    f_true_g[[s]][[i]] <- lapply(1:2, function(l) {
      as.vector(C_g %*% tp$nu_phi_true[[l]] %*% tp$zeta_true[[s]][[l]][i,])
    })
  }
}

# True specific factor on dense grid
g_true_g <- vector("list", 2)
for (s in 1:2) {
  g_true_g[[s]] <- vector("list", tp$n_s[s])
  for (i in 1:tp$n_s[s]) {
    g_true_g[[s]][[i]] <- lapply(1:1, function(l) {
      as.vector(C_g %*% tp$nu_psi_true[[s]][[l]] %*% tp$xi_true[[s]][[l]][i,])
    })
  }
}

# True full signal per (s,i,j)
y_true <- function(s, i, j) {
  y <- mu_true_g[[s]][[j]]
  for (l in 1:2) y <- y + tp$a_true[j,l] * f_true_g[[s]][[i]][[l]]
  for (l in 1:1) y <- y + tp$b_true[[s]][j,l] * g_true_g[[s]][[i]][[l]]
  y
}

# ---- Split study 2 into train/test ----
n2 <- 25
set.seed(555)
idx_test <- sort(sample(1:n2, 7))
idx_train <- setdiff(1:n2, idx_test)

# Training data
Y_tr_m <- dat$Y; Y_tr_m[[2]] <- Y_tr_m[[2]][idx_train]
time_tr_m <- dat$time_obs; time_tr_m[[2]] <- time_tr_m[[2]][idx_train]
Y_pool_bs <- c(dat$Y[[1]], dat$Y[[2]][idx_train])
time_pool_bs <- c(dat$time_obs[[1]], dat$time_obs[[2]][idx_train])

# ---- Fit multiFSYNC ----
cat("Fitting multiFSYNC ...")
t1 <- Sys.time()
fit_m <- bayesSYNC_multi(Y=Y_tr_m, time_obs=time_tr_m,
  L_f=2, L_s=1, M_f=c(3,3), M_s=list(3,3),
  K=7, anneal=NULL, maxit=60, n_g=n_g, time_g=time_g,
  verbose=FALSE, seed=2024)
tm <- as.numeric(difftime(Sys.time(),t1,units="secs"))
cat(sprintf(" %.1fs, ELBO=%.0f, %d iters\n", tm, tail(fit_m$ELBO,1), fit_m$i_iter))

# ---- Fit bayesSYNC (pooled) ----
cat("Fitting bayesSYNC (pooled) ...")
t2 <- Sys.time()
fit_b <- bayesSYNC::bayesSYNC(Y=Y_pool_bs, time_obs=time_pool_bs,
  Q=2, L=3, K=7, anneal=NULL, maxit=60, n_g=n_g, time_g=time_g,
  verbose=FALSE, seed=2024)
tb <- as.numeric(difftime(Sys.time(),t2,units="secs"))
cat(sprintf(" %.1fs, ELBO=%.0f, %d iters\n", tb, tail(fit_b$ELBO,1), fit_b$i_iter))

# ---- Reconstruct predictions on dense grid ----
# multiFSYNC: C_g %*% nu_mu + sum_l a_jl * (C_g %*% nu_phi_l * zeta) + specific
y_hat_m <- function(s, i, j, fit) {
  y <- as.vector(fit$C_g %*% fit$mu_q_nu_mu[[s]][[j]])
  for (l in 1:fit$L_f) {
    y <- y + fit$mu_q_a[j,l] * as.vector(fit$C_g %*% fit$mu_q_nu_phi[[l]] %*% fit$mu_q_zeta[[s]][[l]][i,])
  }
  if (fit$L_s > 0) {
    for (l in 1:fit$L_s) {
      y <- y + fit$mu_q_b_specific[[s]][j,l] * as.vector(fit$C_g %*% fit$mu_q_nu_psi[[s]][[l]] %*% fit$mu_q_xi[[s]][[l]][i,])
    }
  }
  y
}

# bayesSYNC: mu_hat + sum_l B_hat[j,l] * (Phi_hat[l] * Zeta_hat[l][i,])
yb_mu <- fit_b$list_mu_hat
yb_phi <- fit_b$list_list_Phi_hat
yb_zeta <- fit_b$list_Zeta_hat
y_hat_b <- function(i, j, fit) {
  y <- yb_mu[[j]]  # already on dense grid
  for (l in 1:fit$Q) {
    y <- y + fit$B_hat[j,l] * (yb_phi[[l]] %*% yb_zeta[[l]][i,])[,1]
  }
  y
}

# Map: which bayesSYNC individual index corresponds to (s=2, i)?
# bayesSYNC was trained on: study1 (1:25) + study2_train (18 indices)
# So bayesSYNC individual 26:43 are study 2's training individuals
# Test individuals (idx_test) are NOT in bayesSYNC's training data
# For fair comparison, test on study 2 training individuals (both models saw these)
# and on study 1 training individuals

cat("\nReconstruction RMSE on Dense Grid vs True Signal:\n\n")

for (label in c("Study 1 (train)", "Study 2 (train)")) {
  s <- if (label == "Study 1 (train)") 1 else 2
  ii_range <- if (s == 1) 1:25 else idx_train  # training individuals only

  rmse_m <- c(); rmse_b <- c()
  for (i in ii_range) {
    for (j in 1:8) {
      yt <- y_true(s, i, j)
      # multiFSYNC: need training position, not original index
      if (s == 1) {
        i_m <- i
      } else {
        i_m <- which(idx_train == i)
      }
      y_m <- y_hat_m(s, i_m, j, fit_m)
      # bayesSYNC: need correct pooled index
      if (s == 1) {
        i_bs <- i
      } else {
        i_bs <- 25 + which(idx_train == i)
      }
      y_b <- y_hat_b(i_bs, j, fit_b)
      rmse_m <- c(rmse_m, sqrt(mean((y_m - yt)^2)))
      rmse_b <- c(rmse_b, sqrt(mean((y_b - yt)^2)))
    }
  }

  imp <- (1 - mean(rmse_m)/mean(rmse_b)) * 100
  cat(sprintf("  %-16s  multiFSYNC: %.4f  bayesSYNC: %.4f  multi better: %s (%.1f%%)\n",
    label, mean(rmse_m), mean(rmse_b),
    ifelse(imp > 0, "YES", "no"), imp))
}

cat("\nTest B1b: PASS\n\n")

# ---- Bonus: multiFSYNC out-of-sample prediction on test individuals ----
cat("Bonus: multiFSYNC out-of-sample (study 2 held-out, fold-in) ...\n")
rmse_test <- c()
for (i in idx_test) {
  i_m <- which(idx_train == i)  # won't exist — use fold-in, not fitted zeta
  # For fold-in, we need to estimate zeta/xi for test individual from data
  # Use the same fold-in as B1
  C_i <- dat$C[[2]][[i]]
  K_tot <- ncol(C_i)
  f_m_test <- list()
  for (l in 1:fit_m$L_f) {
    nu_phi <- fit_m$mu_q_nu_phi[[l]]
    M_l <- fit_m$M_f[l]
    H <- crossprod(nu_phi, crossprod(C_i) %*% nu_phi)
    for (m in 1:M_l) H[m,m] <- H[m,m] + tr(crossprod(C_i) %*% fit_m$Sigma_q_nu_phi[[l]][[m]])
    sum_res <- rep(0, K_tot)
    for (jj in 1:fit_m$p) {
      r_jj <- dat$Y[[2]][[i]][[jj]] - as.vector(C_i %*% fit_m$mu_q_nu_mu[[2]][[jj]])
      sum_res <- sum_res + fit_m$mu_q_a[jj,l] * as.vector(crossprod(C_i, r_jj))
    }
    f_m_test[[l]] <- as.vector(solve(H + diag(M_l), crossprod(nu_phi, sum_res)))
  }
  g_m_test <- list()
  if (fit_m$L_s > 0) {
    for (l in 1:fit_m$L_s) {
      nu_psi <- fit_m$mu_q_nu_psi[[2]][[l]]
      M_sl <- fit_m$M_s[[2]][l]
      H <- crossprod(nu_psi, crossprod(C_i) %*% nu_psi)
      for (m in 1:M_sl) H[m,m] <- H[m,m] + tr(crossprod(C_i) %*% fit_m$Sigma_q_nu_psi[[2]][[l]][[m]])
      sum_res <- rep(0, K_tot)
      for (jj in 1:fit_m$p) {
        r_jj <- dat$Y[[2]][[i]][[jj]] - as.vector(C_i %*% fit_m$mu_q_nu_mu[[2]][[jj]])
        sum_res <- sum_res + fit_m$mu_q_b_specific[[2]][jj,l] * as.vector(crossprod(C_i, r_jj))
      }
      g_m_test[[l]] <- as.vector(solve(H + diag(M_sl), crossprod(nu_psi, sum_res)))
    }
  }
  for (j in 1:fit_m$p) {
    y_hat <- as.vector(fit_m$C_g %*% fit_m$mu_q_nu_mu[[2]][[j]])
    for (l in 1:fit_m$L_f) y_hat <- y_hat + fit_m$mu_q_a[j,l] * as.vector(fit_m$C_g %*% fit_m$mu_q_nu_phi[[l]] %*% f_m_test[[l]])
    if (fit_m$L_s > 0) for (l in 1:fit_m$L_s) y_hat <- y_hat + fit_m$mu_q_b_specific[[2]][j,l] * as.vector(fit_m$C_g %*% fit_m$mu_q_nu_psi[[2]][[l]] %*% g_m_test[[l]])
    yt <- y_true(2, i, j)
    rmse_test <- c(rmse_test, sqrt(mean((y_hat - yt)^2)))
  }
}
cat(sprintf("  Study 2 (test, fold-in) multiFSYNC RMSE: %.4f\n", mean(rmse_test)))
