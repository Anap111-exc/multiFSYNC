library(splines)
library(parallel)
library(pracma)
library(testthat)

src_dir <- "D:/文档/Factor model/多研究函数型因子模型实践/Rcode/multiFSYNC/R"
for (f in list.files(src_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

cat("========================================\n")
cat("  Phase 1.1-2 Tests\n")
cat("========================================\n\n")

# ---- Test 1: simulate_multi_study_data() ----
cat("Test 1: simulate_multi_study_data() ... ")
set.seed(42)
dat <- simulate_multi_study_data(
  S = 2, n_s = c(10, 15), p = 5, d = 2,
  L_f = 2, L_s = 1,
  M_f = c(2, 1),
  M_s = list(2, 2),
  n_obs = 50, sigma_eps = 0.1
)
stopifnot(is.list(dat$Y))
stopifnot(is.list(dat$Z) || is.null(dat$Z))
stopifnot(is.list(dat$true_params))
cat("PASS\n")

# ---- Test 2: Y structure ----
cat("Test 2: Y structure ... ")
S <- length(dat$Y)
n_s <- c(10, 15)
p <- 5
stopifnot(S == 2)
for (s in 1:S) {
  stopifnot(length(dat$Y[[s]]) == n_s[s])
  for (i in 1:n_s[s]) {
    stopifnot(length(dat$Y[[s]][[i]]) == p)
    stopifnot(is.vector(dat$Y[[s]][[i]][[1]]))
  }
}
cat("PASS\n")

# ---- Test 3: Z structure ----
cat("Test 3: Z structure ... ")
stopifnot(length(dat$Z) == S)
for (s in 1:S) {
  stopifnot(dim(dat$Z[[s]]) == c(n_s[s], 2))
}
cat("PASS\n")

# ---- Test 4: C design matrices ----
cat("Test 4: C design matrices ... ")
C <- dat$C
K <- dat$true_params$K
n_obs <- 50
for (s in 1:S) {
  for (i in 1:n_s[s]) {
    stopifnot(dim(C[[s]][[i]]) == c(n_obs, K + 2))
  }
}
cat("PASS\n")

# ---- Test 5: Variational parameter initialization ----
cat("Test 5: Variational parameter initialization ... ")

# Use a minimal fitting call to trigger initialization
fit <- bayesSYNC_multi(
  Y = dat$Y, Z = dat$Z,
  time_obs = dat$time_obs,
  L_f = 2, L_s = 1,
  M_f = c(2, 1), M_s = list(2, 2),
  K = K, maxit = 1, verbose = FALSE
)

K_total <- K + 2

# mu_q_nu_mu: list[[s]][[j]] of length K_total
stopifnot(length(fit$mu_q_nu_mu) == S)
stopifnot(length(fit$mu_q_nu_mu[[1]]) == p)
stopifnot(length(fit$mu_q_nu_mu[[1]][[1]]) == K_total)
cat("mu_q_nu_mu OK; ")

# mu_q_nu_phi: list[[l]] with ncol = M_f[l]
stopifnot(length(fit$mu_q_nu_phi) == 2)
stopifnot(ncol(fit$mu_q_nu_phi[[1]]) == 2)  # M_f[1]=2
stopifnot(ncol(fit$mu_q_nu_phi[[2]]) == 1)  # M_f[2]=1
cat("mu_q_nu_phi OK; ")

# mu_q_nu_psi: list[[s]][[l]]
stopifnot(length(fit$mu_q_nu_psi) == S)
stopifnot(length(fit$mu_q_nu_psi[[1]]) == 1)
stopifnot(ncol(fit$mu_q_nu_psi[[1]][[1]]) == 2)  # M_s[[1]][1]=2
cat("mu_q_nu_psi OK; ")

# mu_q_zeta: list[[s]][[l]], n_s[s] x M_f[l]
stopifnot(length(fit$mu_q_zeta[[1]]) == 2)
stopifnot(dim(fit$mu_q_zeta[[1]][[1]]) == c(10, 2))
stopifnot(dim(fit$mu_q_zeta[[1]][[2]]) == c(10, 1))
cat("mu_q_zeta OK; ")

# mu_q_xi: list[[s]][[l]], n_s[s] x M_s[[s]][l]
stopifnot(dim(fit$mu_q_xi[[1]][[1]]) == c(10, 2))
cat("mu_q_xi OK; ")

# Sigma_q_nu_mu: list[[s]][[j]], K_total x K_total positive definite
stopifnot(dim(fit$Sigma_q_nu_mu[[1]][[1]]) == c(K_total, K_total))
stopifnot(all(eigen(fit$Sigma_q_nu_mu[[1]][[1]])$values > 0))
cat("Sigma_q_nu_mu OK; ")

# Sigma_q_zeta: list[[s]][[l]][[i]], M_f[l] x M_f[l]
stopifnot(dim(fit$Sigma_q_zeta[[1]][[1]][[1]]) == c(2, 2))
stopifnot(dim(fit$Sigma_q_zeta[[1]][[2]][[1]]) == c(1, 1))
cat("Sigma_q_zeta OK; ")

# Loadings: mu_q_a p x L_f
stopifnot(dim(fit$mu_q_a) == c(p, 2))
stopifnot(dim(fit$mu_q_gamma_a) == c(p, 2))
cat("loadings a OK; ")

# Loadings b: list[[s]] p x L_s
stopifnot(dim(fit$mu_q_b_specific[[1]]) == c(p, 1))
cat("loadings b OK; ")

# mu_q_nu_beta (d=2)
stopifnot(length(fit$mu_q_nu_beta) == p)
stopifnot(length(fit$mu_q_nu_beta[[1]]) == 2)
stopifnot(length(fit$mu_q_nu_beta[[1]][[1]]) == K_total)
cat("mu_q_nu_beta OK\n")

# ---- Test 6: Precomputation objects ----
cat("Test 6: Precomputation objects ... ")
stopifnot(length(fit$list_cp_C) == S)
stopifnot(dim(fit$list_cp_C[[1]][[1]]) == c(K_total, K_total))
stopifnot(dim(fit$list_cp_C_Y[[1]][[1]]) == c(K_total, p))
cat("PASS\n")

cat("\n========================================\n")
cat("  ALL Phase 1.1-2 TESTS PASSED\n")
cat("========================================\n")
