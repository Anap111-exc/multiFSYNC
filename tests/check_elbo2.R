library(splines); library(parallel); library(pracma)
suppressMessages(devtools::load_all("D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC", quiet=TRUE))

set.seed(42)

# Test 1: L_f=0 (pure mean model)
cat("=== Test A: L_f=0, L_s=0 (no factors, no covariates) ===\n")
dat0 <- simulate_multi_study_data(S=1, n_s=c(20), p=3, d=0, L_f=0, L_s=0,
  K=6, n_obs=30, sigma_eps=0.05)

fit0 <- bayesSYNC_multi(Y=dat0$Y, time_obs=dat0$time_obs,
  L_f=0, L_s=0, M_f=integer(0), M_s=list(integer(0)),
  K=6, anneal=NULL, maxit=30, verbose=FALSE, seed=42)

e0 <- fit0$ELBO
d0 <- diff(e0)
cat(sprintf("Iters: %d, Increases: %d, Decreases: %d\n", length(e0), sum(d0>0), sum(d0<0)))
cat("Last 10 diffs:", paste(round(tail(d0, 10), 2), collapse=", "), "\n\n")

# Test 2: L_f=1 no anneal, more iters
cat("=== Test B: L_f=1, L_s=0, 50 iters ===\n")
dat <- simulate_multi_study_data(S=1, n_s=c(20), p=3, d=0, L_f=1, L_s=0,
  M_f=c(2), K=6, n_obs=30, sigma_eps=0.05)

fit <- bayesSYNC_multi(Y=dat$Y, time_obs=dat$time_obs,
  L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)),
  K=6, anneal=NULL, maxit=50, verbose=FALSE, seed=42)

e <- fit$ELBO
d <- diff(e)
cat(sprintf("Iters: %d, Increases: %d, Decreases: %d\n", length(e), sum(d>0), sum(d<0)))
cat("First 5 diffs:", paste(round(head(d, 5), 2), collapse=", "), "\n")
cat("Last  10 diffs:", paste(round(tail(d, 10), 2), collapse=", "), "\n")
cat(sprintf("First ELBO: %.2f  Last ELBO: %.2f\n", e[1], tail(e,1)))
cat(sprintf("Best (iter %d): %.2f\n", which.max(e), max(e)))
