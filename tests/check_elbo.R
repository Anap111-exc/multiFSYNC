library(splines); library(parallel); library(pracma)
suppressMessages(devtools::load_all("D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC", quiet=TRUE))

set.seed(42)
dat <- simulate_multi_study_data(S=1, n_s=c(20), p=3, d=0, L_f=1, L_s=0,
  M_f=c(2), K=6, n_obs=30, sigma_eps=0.05)

fit <- bayesSYNC_multi(Y=dat$Y, time_obs=dat$time_obs,
  L_f=1, L_s=0, M_f=c(2), M_s=list(integer(0)),
  K=6, anneal=NULL, maxit=50, verbose=FALSE, seed=42)

elbo <- fit$ELBO
n <- length(elbo)

cat(sprintf("Total iterations: %d\n", n))
cat(sprintf("First ELBO: %.2f\n", elbo[1]))
cat(sprintf("Last ELBO:  %.2f\n", elbo[n]))
cat(sprintf("Total change: %.2f\n\n", elbo[n] - elbo[1]))

cat("ELBO trajectory:\n")
for (i in 1:n) {
  change <- if (i > 1) elbo[i] - elbo[i-1] else NA
  marker <- if (i > 1 && change < 0) " <== DECREASE" else ""
  if (i <= 5 || i > n-5 || !is.na(marker) && nchar(marker) > 0) {
    cat(sprintf("  iter %3d: %.2f  (diff: %+.2f)%s\n", i, elbo[i], if(i>1) change else 0, marker))
  }
}

ndec <- sum(diff(elbo) < 0)
ninc <- sum(diff(elbo) > 0)
cat(sprintf("\nIncreases: %d, Decreases: %d, Ratio: %.1f%%\n", ninc, ndec, ninc/(ninc+ndec)*100))
cat(sprintf("Max decrease: %.2f, Max increase: %.2f\n", min(diff(elbo)), max(diff(elbo))))
