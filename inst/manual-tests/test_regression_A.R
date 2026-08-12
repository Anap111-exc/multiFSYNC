# Regression Test A: S=1 Degeneracy (N=30,50,80)
# Verifies numerical consistency after K_f→L_f, L_f→M_f rename
suppressMessages({
  library(splines); library(parallel)
  source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
  devtools::load_all(project_root, quiet=TRUE)
  suppressMessages(library(bayesSYNC, quietly=TRUE))
})

cat("========================================\n")
cat("  Regression Test A: S=1 Degeneracy\n")
cat("  N=30, 50, 80  |  p=8  |  Q=2  |  L=3\n")
cat("========================================\n\n")

run_one <- function(N, seed) {
  set.seed(seed)

  # Old params: K_f=2, L_f=3  ->  New: L_f=2, M_f=c(3,3)
  dat <- simulate_multi_study_data(
    S=1, n_s=c(N), p=8, d=0,
    L_f=2, L_s=0, M_f=c(3,3), M_s=list(integer(0)),
    K=8, n_obs=80, sigma_eps=0.1,
    bool_sparse_loadings=TRUE, prop_sparse=0.6
  )

  # multiFSYNC
  t1 <- Sys.time()
  fit_m <- bayesSYNC_multi(
    Y=dat$Y, time_obs=dat$time_obs,
    L_f=2, L_s=0, M_f=c(3,3), M_s=list(integer(0)),
    K=8, anneal=NULL, maxit=100, verbose=FALSE, seed=seed
  )
  tm <- as.numeric(difftime(Sys.time(), t1, units="secs"))

  # bayesSYNC reference
  t2 <- Sys.time()
  fit_b <- bayesSYNC::bayesSYNC(
    Y=dat$Y[[1]], time_obs=dat$time_obs[[1]],
    Q=2, L=3, K=8, anneal=NULL, maxit=100, verbose=FALSE, seed=seed
  )
  tb <- as.numeric(difftime(Sys.time(), t2, units="secs"))

  # ---- Metrics ----
  tp <- dat$true_params

  # Loadings MSE (normalized by true norm)
  a_mse <- function(fit_a, true_a) {
    if (ncol(fit_a) == 0) return(1.0)
    fit_a <- as.matrix(fit_a)
    for (l in 1:ncol(true_a)) {
      if (sd(fit_a[,l]) < 1e-10) next  # degenerate factor
      if (cor(fit_a[,l], true_a[,l], use="complete.obs") < 0) fit_a[,l] <- -fit_a[,l]
    }
    mean((fit_a - true_a)^2) / mean(true_a^2)
  }
  load_mse_m <- a_mse(fit_m$mu_q_a, tp$a_true)
  load_mse_b <- a_mse(fit_b$B_hat, tp$a_true)

  # Score correlation (absolute, best matching)
  score_cor <- function(fit_zeta, true_zeta) {
    if (ncol(fit_zeta) == 0) return(0)
    cors <- sapply(1:ncol(true_zeta), function(m) {
      if (sd(fit_zeta[,m]) < 1e-10) return(0)
      abs(cor(fit_zeta[,m], true_zeta[,m]))
    })
    mean(cors)
  }
  # multiFSYNC zeta: study [[1]], factor [[1]], matrix n_s x M_f[1]
  zeta_m <- fit_m$mu_q_zeta[[1]][[1]]
  zeta_true <- tp$zeta_true[[1]][[1]]
  sc_m <- score_cor(zeta_m, zeta_true)

  # bayesSYNC zeta: list_Zeta_hat[[l]] of (N x L) matrices
  zeta_b <- fit_b$list_Zeta_hat[[1]]
  sc_b <- score_cor(zeta_b, zeta_true)

  # PPI
  ppi_m <- colMeans(fit_m$mu_q_gamma_a)
  ppi_b <- colMeans(fit_b$ppi)

  # ELBO
  elbo_m <- tail(fit_m$ELBO, 1)
  elbo_b <- tail(fit_b$ELBO, 1)

  list(
    N = N,
    load_mse_m = load_mse_m, load_mse_b = load_mse_b,
    score_cor_m = sc_m, score_cor_b = sc_b,
    ppi_m = ppi_m, ppi_b = ppi_b,
    elbo_m = elbo_m, elbo_b = elbo_b,
    iter_m = fit_m$i_iter, iter_b = fit_b$i_iter,
    time_m = tm, time_b = tb
  )
}

results <- list()
for (N in c(30, 50, 80)) {
  cat(sprintf("N = %d ...\n", N))
  results[[as.character(N)]] <- run_one(N, seed=42)
}

cat("\n========================================\n")
cat("  Results vs Old Report\n")
cat("========================================\n\n")

cat(sprintf("%-8s %-12s %-12s %-12s %-12s %-10s %-10s\n",
  "N", "loadMSE_m", "loadMSE_b", "scoreCor_m", "scoreCor_b", "PPI_m", "PPI_b"))
cat(strrep("-", 78), "\n")
for (r in results) {
  cat(sprintf("%-8d %-12.3f %-12.3f %-12.3f %-12.3f %-10s %-10s\n",
    r$N,
    r$load_mse_m, r$load_mse_b,
    r$score_cor_m, r$score_cor_b,
    paste(round(r$ppi_m, 2), collapse=","),
    paste(round(r$ppi_b, 2), collapse=",")))
}

cat("\n--- Old Report (multiFSYNC only) ---\n")
cat("N=30: loadMSE=0.798  scoreCor=0.514  PPI=1,1\n")
cat("N=50: loadMSE=0.465  scoreCor=0.202  PPI=1,1\n")
cat("N=80: loadMSE=0.667  scoreCor=0.573  PPI=1,1\n")

cat("\n========================================\n")
cat("  Runtime Summary\n")
cat("========================================\n")
for (r in results) {
  cat(sprintf("N=%d: multiFSYNC %.1fs  bayesSYNC %.1fs  (iter %d | %d)\n",
    r$N, r$time_m, r$time_b, r$iter_m, r$iter_b))
}

cat("\nRegression Test A complete.\n")
