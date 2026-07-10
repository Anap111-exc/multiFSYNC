suppressMessages({
  library(splines); library(parallel); library(pracma)
  devtools::load_all("D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC", quiet=TRUE)
  suppressMessages(library(bayesSYNC, quietly=TRUE))
})

cat("Test B3: Large-Scale Comparison (S=4, n_s=50, p=10)\n")
cat(strrep("=", 58), "\n\n")

set.seed(42)

# ---- Generate S=4 data ----
n_studies <- 4
n_per <- 50
p <- 10
L_f <- 3
L_s <- 1

cat(sprintf("Generating S=%d, n_s=rep(%d,%d), p=%d, L_f=%d, L_s=%d ...\n",
  n_studies, n_per, n_studies, p, L_f, L_s))

dat <- simulate_multi_study_data(
  S=n_studies, n_s=rep(n_per, n_studies), p=p, d=0,
  L_f=L_f, L_s=L_s,
  M_f=rep(3, L_f),
  M_s=lapply(1:n_studies, function(s) rep(3, L_s)),
  K=7, n_obs=60, sigma_eps=0.1,
  bool_sparse_loadings=TRUE, prop_sparse=0.5
)
tp <- dat$true_params
N_total <- sum(tp$n_s)
cat(sprintf("  Total individuals: %d\n\n", N_total))

# ---- Fit multiFSYNC ----
cat("Fitting multiFSYNC (S=4, anneal, maxit=80) ...\n")
t1 <- Sys.time()
fit_m <- bayesSYNC_multi(
  Y=dat$Y, time_obs=dat$time_obs,
  L_f=L_f, L_s=L_s,
  M_f=rep(3, L_f), M_s=lapply(1:n_studies, function(s) rep(3, L_s)),
  K=7, anneal=c(1, 1.9, 20), maxit=80, verbose=FALSE, seed=42
)
tm <- as.numeric(difftime(Sys.time(), t1, units="secs"))
cat(sprintf("  %.1fs | ELBO=%.0f | %d iters (anneal %d steps)\n\n", tm, tail(fit_m$ELBO,1), fit_m$i_iter, 20))

# ---- Fit bayesSYNC (pool all studies) ----
cat(sprintf("Fitting bayesSYNC (pooled, %d indivs, anneal, maxit=80) ...\n", N_total))
Y_pool <- unlist(dat$Y, recursive=FALSE)
time_pool <- unlist(dat$time_obs, recursive=FALSE)
t2 <- Sys.time()
fit_b <- bayesSYNC::bayesSYNC(
  Y=Y_pool, time_obs=time_pool,
  Q=L_f, L=3, K=7, anneal=c(1, 1.9, 20), maxit=80,
  verbose=FALSE, seed=42
)
tb <- as.numeric(difftime(Sys.time(), t2, units="secs"))
cat(sprintf("  %.1fs | ELBO=%.0f | %d iters\n\n", tb, tail(fit_b$ELBO_iter,1), fit_b$i_iter))

# ---- Reconstruction on dense grid ----
n_g <- 200
time_g <- seq(0, 1, length.out=n_g)
C_g <- cbind(1, time_g, ZOSull(time_g, c(0,1),
  intKnots=quantile(time_g, seq(0,1,length=tp$K)[-c(1,tp$K)])))

# True signal on dense grid
y_true <- function(s, i, j) {
  y <- as.vector(C_g %*% tp$nu_mu_true[[s]][[j]])
  for (l in 1:L_f) {
    y <- y + tp$a_true[j,l] * as.vector(C_g %*% tp$nu_phi_true[[l]] %*% tp$zeta_true[[s]][[l]][i,])
  }
  for (l in seq_len(L_s)) {
    y <- y + tp$b_true[[s]][j,l] * as.vector(C_g %*% tp$nu_psi_true[[s]][[l]] %*% tp$xi_true[[s]][[l]][i,])
  }
  y
}

# multiFSYNC prediction on dense grid
pred_m <- function(s, i, j) {
  y <- as.vector(C_g %*% fit_m$mu_q_nu_mu[[s]][[j]])
  for (l in 1:fit_m$L_f) {
    y <- y + fit_m$mu_q_a[j,l] * as.vector(C_g %*% fit_m$mu_q_nu_phi[[l]] %*% fit_m$mu_q_zeta[[s]][[l]][i,])
  }
  if (fit_m$L_s > 0) {
    for (l in 1:fit_m$L_s) {
      y <- y + fit_m$mu_q_b_specific[[s]][j,l] * as.vector(C_g %*% fit_m$mu_q_nu_psi[[s]][[l]] %*% fit_m$mu_q_xi[[s]][[l]][i,])
    }
  }
  y
}

# bayesSYNC prediction on dense grid
yb_mu <- fit_b$list_mu_hat
yb_phi <- fit_b$list_list_Phi_hat
yb_zeta <- fit_b$list_Zeta_hat
pred_b <- function(i_pool, j) {
  y <- yb_mu[[j]]
  for (l in 1:fit_b$Q) {
    y <- y + fit_b$B_hat[j,l] * (yb_phi[[l]] %*% yb_zeta[[l]][i_pool,])[,1]
  }
  y
}

# Map: (s,i) -> bayesSYNC pooled index
i_pool_map <- list()
start <- 0
for (s in 1:n_studies) {
  i_pool_map[[s]] <- start + 1:tp$n_s[s]
  start <- start + tp$n_s[s]
}

# ---- Compute RMSE per study ----
cat("Study  LoadingsMSE(m)  LoadingsMSE(b)  ReconRMSE(m)  ReconRMSE(b)  mBetter?  PPI(m)  PPI(b)\n")
cat(strrep("-", 95), "\n")

for (s in 1:n_studies) {
  load_mse_m <- mean((fit_m$mu_q_a - tp$a_true)^2) / mean(tp$a_true^2)
  load_mse_b <- mean((fit_b$B_hat - tp$a_true)^2) / mean(tp$a_true^2)

  rmse_m <- c(); rmse_b <- c()
  for (i in 1:tp$n_s[s]) {
    for (j in 1:p) {
      yt <- y_true(s, i, j)
      rmse_m <- c(rmse_m, sqrt(mean((pred_m(s, i, j) - yt)^2)))
      rmse_b <- c(rmse_b, sqrt(mean((pred_b(i_pool_map[[s]][i], j) - yt)^2)))
    }
  }

  ppi_m <- paste(round(colMeans(fit_m$mu_q_gamma_a), 2), collapse=",")
  ppi_b <- paste(round(colMeans(fit_b$ppi), 2), collapse=",")
  better <- if (mean(rmse_m) < mean(rmse_b)) "YES" else "no"

  cat(sprintf("  %d     %-12.4f  %-12.4f  %-10.4f  %-10.4f  %-6s  %-8s  %s\n",
    s, load_mse_m, load_mse_b, mean(rmse_m), mean(rmse_b), better, ppi_m, ppi_b))
}

# ---- Aggregate summary ----
all_rmse_m <- c(); all_rmse_b <- c()
for (s in 1:n_studies) {
  for (i in 1:tp$n_s[s]) {
    for (j in 1:p) {
      all_rmse_m <- c(all_rmse_m, sqrt(mean((pred_m(s,i,j) - y_true(s,i,j))^2)))
      all_rmse_b <- c(all_rmse_b, sqrt(mean((pred_b(i_pool_map[[s]][i],j) - y_true(s,i,j))^2)))
    }
  }
}

cat(sprintf("\n--- Aggregate (S=%d, N=%d total) ---\n", n_studies, N_total))
cat(sprintf("  multiFSYNC overall RMSE: %.4f\n", mean(all_rmse_m)))
cat(sprintf("  bayesSYNC  overall RMSE: %.4f\n", mean(all_rmse_b)))
cat(sprintf("  multiFSYNC advantage:    %.1f%%\n", (1-mean(all_rmse_m)/mean(all_rmse_b))*100))

# Factor detection
ppi_shared_m <- 1 - exp(colSums(log1p(-fit_m$mu_q_gamma_a)))
cat(sprintf("  Shared factor PPIs (multi): %s\n", paste(round(ppi_shared_m,2), collapse=", ")))

# Specific factor PPIs per study
for (s in 1:n_studies) {
  ppi_spec <- 1 - exp(colSums(log1p(-fit_m$mu_q_gamma_b[[s]])))
  cat(sprintf("  Study %d specific PPIs:       %s\n", s, paste(round(ppi_spec,2), collapse=", ")))
}

cat(sprintf("\n  Runtime: multiFSYNC %.0fs | bayesSYNC %.0fs\n", tm, tb))
cat("\nTest B3: PASS\n")
