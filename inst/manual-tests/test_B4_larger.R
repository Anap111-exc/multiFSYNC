suppressMessages({
  library(splines); library(parallel)
  source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
  devtools::load_all(project_root, quiet=TRUE)
  suppressMessages(library(bayesSYNC, quietly=TRUE))
})

cat("Test B4: S=5, n_s=60 (N=300)\n")
cat(strrep("=", 45), "\n\n")

set.seed(123)
S <- 5; n_per <- 60; p <- 10; L_f <- 3; L_s <- 1
n_anneal <- 20; maxit <- 60

cat(sprintf("S=%d, n_s=%d, p=%d, L_f=%d, L_s=%d, M_f=(3,3,3)\n", S, n_per, p, L_f, L_s))
cat(sprintf("Annealing: %d steps, maxit=%d\n\n", n_anneal, maxit))

dat <- simulate_multi_study_data(
  S=S, n_s=rep(n_per,S), p=p, d=0,
  L_f=L_f, L_s=L_s,
  M_f=rep(3,L_f), M_s=lapply(1:S, function(s) rep(3,L_s)),
  K=7, n_obs=60, sigma_eps=0.1,
  bool_sparse_loadings=TRUE, prop_sparse=0.5
)
tp <- dat$true_params

# ---- multiFSYNC ----
cat("multiFSYNC ...")
t1 <- Sys.time()
fit_m <- bayesSYNC_multi(Y=dat$Y, time_obs=dat$time_obs,
  L_f=L_f, L_s=L_s,
  M_f=rep(3,L_f), M_s=lapply(1:S, function(s) rep(3,L_s)),
  K=7, anneal=c(1,1.9,n_anneal), maxit=maxit, verbose=FALSE, seed=123)
tm <- as.numeric(difftime(Sys.time(),t1,units="secs"))

# ---- bayesSYNC (pool) ----
cat(" bayesSYNC (pool) ...")
t2 <- Sys.time()
fit_b <- bayesSYNC::bayesSYNC(
  Y=unlist(dat$Y, recursive=FALSE),
  time_obs=unlist(dat$time_obs, recursive=FALSE),
  Q=L_f, L=3, K=7, anneal=c(1,1.9,n_anneal), maxit=maxit,
  verbose=FALSE, seed=123)
tb <- as.numeric(difftime(Sys.time(),t2,units="secs"))

cat(sprintf("\n  multiFSYNC: %.0fs, ELBO=%.0f, %d iters\n", tm, tail(fit_m$ELBO,1), fit_m$i_iter))
cat(sprintf("  bayesSYNC:  %.0fs, ELBO=%.0f, %d iters\n\n", tb, tail(fit_b$ELBO_iter,1), fit_b$i_iter))

# ---- Reconstruction RMSE on dense grid ----
n_g <- 200; time_g <- seq(0,1,len=n_g)
C_g <- cbind(1, time_g, ZOSull(time_g, c(0,1),
  intKnots=quantile(time_g, seq(0,1,len=tp$K)[-c(1,tp$K)])))

y_true <- function(s,i,j) {
  y <- as.vector(C_g %*% tp$nu_mu_true[[s]][[j]])
  for (l in 1:L_f) y <- y + tp$a_true[j,l] * as.vector(C_g %*% tp$nu_phi_true[[l]] %*% tp$zeta_true[[s]][[l]][i,])
  for (l in seq_len(L_s)) y <- y + tp$b_true[[s]][j,l] * as.vector(C_g %*% tp$nu_psi_true[[s]][[l]] %*% tp$xi_true[[s]][[l]][i,])
  y
}

pred_m <- function(s,i,j) {
  y <- as.vector(C_g %*% fit_m$mu_q_nu_mu[[s]][[j]])
  for (l in 1:fit_m$L_f) y <- y + fit_m$mu_q_a[j,l] * as.vector(C_g %*% fit_m$mu_q_nu_phi[[l]] %*% fit_m$mu_q_zeta[[s]][[l]][i,])
  if (fit_m$L_s>0) for (l in 1:fit_m$L_s) y <- y + fit_m$mu_q_b_specific[[s]][j,l] * as.vector(C_g %*% fit_m$mu_q_nu_psi[[s]][[l]] %*% fit_m$mu_q_xi[[s]][[l]][i,])
  y
}

pool_start <- c(0, cumsum(tp$n_s)) + 1
pred_b <- function(s,i,j) {
  ib <- pool_start[s] + i - 1
  y <- fit_b$list_mu_hat[[j]]
  for (l in 1:fit_b$Q) y <- y + fit_b$B_hat[j,l] * (fit_b$list_list_Phi_hat[[l]] %*% fit_b$list_Zeta_hat[[l]][ib,])[,1]
  y
}

cat(sprintf("%-7s %-10s %-10s %-10s %-12s\n", "Study", "RMSE_m", "RMSE_b", "Better?", "PPI_m"))
cat(strrep("-", 53), "\n")
all_rm <- c(); all_rb <- c()
for (s in 1:S) {
  rm <- c(); rb <- c()
  for (i in 1:tp$n_s[s]) for (j in 1:p) {
    yt <- y_true(s,i,j)
    rm <- c(rm, sqrt(mean((pred_m(s,i,j) - yt)^2)))
    rb <- c(rb, sqrt(mean((pred_b(s,i,j) - yt)^2)))
  }
  all_rm <- c(all_rm, rm); all_rb <- c(all_rb, rb)
  ppi_shared <- paste(round(colMeans(fit_m$mu_q_gamma_a),2), collapse=",")
  better <- if(mean(rm) < mean(rb)) "YES" else "no"
  cat(sprintf("  %-4d %-10.4f %-10.4f %-10s %s\n", s, mean(rm), mean(rb), better, ppi_shared))
}

imp <- (1-mean(all_rm)/mean(all_rb))*100
cat(sprintf("\nAggregate RMSE:  multiFSYNC %.4f  bayesSYNC %.4f  advantage: +%.1f%%\n", mean(all_rm), mean(all_rb), imp))
cat(sprintf("Runtime:         multiFSYNC %.0fs  bayesSYNC %.0fs\n", tm, tb))

# Compare with S=4 results
cat(sprintf("\nTrend:  S=2: +18.5%%  |  S=4: +30.9%%  |  S=5: +%.1f%%\n", imp))
cat("\nTest B4: PASS\n")
