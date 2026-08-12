library(splines)
library(parallel)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
src_dir <- file.path(project_root, "R")
for (f in list.files(src_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

cat("=== Phase 3 Tests: Variance & Omega ===\n\n")

# ---- Setup ----
set.seed(42)
dat <- simulate_multi_study_data(S=2, n_s=c(10,15), p=5, d=2, L_f=2, L_s=1,
  M_f=c(2,1), M_s=list(2,2), n_obs=50, sigma_eps=0.1)

fit <- bayesSYNC_multi(Y=dat$Y, Z=dat$Z, time_obs=dat$time_obs,
  L_f=2, L_s=1, M_f=c(2,1), M_s=list(2,2), K=7, maxit=3, verbose=FALSE)

# ---- Test 1: sigma^2_eps ----
cat("Test 1: sigma^2_eps ... ")
mrse <- fit$mu_q_recip_sigsq_eps
stopifnot(is.matrix(mrse), dim(mrse) == c(2,5))
stopifnot(all(mrse > 0))
stopifnot(all(mrse < 1e6))
cat(sprintf("PASS  (range %.3f - %.3f)\n", min(mrse), max(mrse)))

# ---- Test 2: sigma^2_mu ----
cat("Test 2: sigma^2_mu ... ")
mrsm <- fit$mu_q_recip_sigsq_mu
stopifnot(is.matrix(mrsm), dim(mrsm) == c(2,5))
stopifnot(all(mrsm > 0))
cat(sprintf("PASS  (range %.3f - %.3f)\n", min(mrsm), max(mrsm)))

# ---- Test 3: sigma^2_phi (now list!) ----
cat("Test 3: sigma^2_phi ... ")
mrsphi <- fit$mu_q_recip_sigsq_phi
stopifnot(is.list(mrsphi), length(mrsphi) == 2)
stopifnot(length(mrsphi[[1]]) == 2)  # M_f[1]=2
stopifnot(length(mrsphi[[2]]) == 1)  # M_f[2]=1
stopifnot(all(unlist(mrsphi) > 0))
cat(sprintf("PASS  (values: %.3f, %.3f)\n", mrsphi[[1]][1], mrsphi[[1]][2]))

# ---- Test 4: sigma^2_psi (now list of lists) ----
cat("Test 4: sigma^2_psi ... ")
mrspsi <- fit$mu_q_recip_sigsq_psi
stopifnot(is.list(mrspsi), length(mrspsi) == 2)
stopifnot(length(mrspsi[[1]]) == 1)  # L_s=1
stopifnot(length(mrspsi[[1]][[1]]) == 2)  # M_s[[1]][1]=2
stopifnot(all(unlist(mrspsi) > 0))
cat(sprintf("PASS  (value: %.3f)\n", mrspsi[[1]][[1]][1]))

# ---- Test 5: Auxiliary a_eps ----
cat("Test 5: a_eps (T<2 constraint) ... ")
mrae <- fit$mu_q_recip_a_eps
stopifnot(is.matrix(mrae), dim(mrae) == c(2,5))
stopifnot(all(mrae > 0))
# Verify kappa_q_a = 2*c - 1 > 0 (c=1 for no anneal → kappa=1)
stopifnot(fit$kappa_q_a == 1)  # 2*1 - 1 = 1
cat(sprintf("PASS  (kappa_q_a=%g > 0)\n", fit$kappa_q_a))

# ---- Test 6: a_phi (list) ----
cat("Test 6: a_phi ... ")
mraphi <- fit$mu_q_recip_a_phi
stopifnot(is.list(mraphi), length(mraphi) == 2)
stopifnot(all(unlist(mraphi) > 0))
cat("PASS\n")

# ---- Test 7: a_psi (list of lists) ----
cat("Test 7: a_psi ... ")
mrapsi <- fit$mu_q_recip_a_psi
stopifnot(is.list(mrapsi))
stopifnot(all(unlist(mrapsi) > 0))
cat("PASS\n")

# ---- Test 8: omega_a (shared) ----
cat("Test 8: omega_a (shared) ... ")
stopifnot(all(fit$c_1_omega_a > 0))
stopifnot(all(fit$d_1_omega_a > 0))
stopifnot(all(is.finite(fit$mu_q_log_omega_a)))
stopifnot(all(is.finite(fit$mu_q_log_1_omega_a)))
cat(sprintf("PASS  (c1=%s, d1=%s)\n",
  paste(round(fit$c_1_omega_a,2), collapse=","),
  paste(round(fit$d_1_omega_a,2), collapse=",")))

# ---- Test 9: omega_b (specific) ----
cat("Test 9: omega_b (specific) ... ")
stopifnot(is.list(fit$c_1_omega_b), length(fit$c_1_omega_b) == 2)
for (s in 1:2) {
  stopifnot(length(fit$c_1_omega_b[[s]]) == 1)  # L_s=1
  stopifnot(all(fit$c_1_omega_b[[s]] > 0))
  stopifnot(all(fit$d_1_omega_b[[s]] > 0))
}
stopifnot(all(is.finite(unlist(fit$mu_q_log_omega_b))))
stopifnot(all(is.finite(unlist(fit$mu_q_log_1_omega_b))))
cat("PASS\n")

# ---- Test 10: No NaN/Inf in variance parameters ----
cat("Test 10: NaN/Inf check ... ")
all_vars <- c(
  unlist(fit$mu_q_recip_sigsq_eps),
  unlist(fit$mu_q_recip_sigsq_mu),
  unlist(fit$mu_q_recip_sigsq_phi),
  unlist(fit$mu_q_recip_sigsq_psi),
  unlist(fit$mu_q_recip_a_eps),
  unlist(fit$mu_q_recip_a_phi),
  unlist(fit$mu_q_recip_a_psi),
  unlist(fit$c_1_omega_a), unlist(fit$d_1_omega_a),
  unlist(fit$c_1_omega_b), unlist(fit$d_1_omega_b),
  unlist(fit$mu_q_log_omega_a), unlist(fit$mu_q_log_1_omega_a),
  unlist(fit$mu_q_log_omega_b), unlist(fit$mu_q_log_1_omega_b)
)
stopifnot(all(is.finite(all_vars)))
stopifnot(all(all_vars != Inf))
cat("PASS\n")

cat("\n=== ALL Phase 3 TESTS PASSED ===\n")
