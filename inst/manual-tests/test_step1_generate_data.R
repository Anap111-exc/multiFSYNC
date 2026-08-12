library(splines)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "helper_project_root.R"))
source(file.path(project_root, "R", "utils_multi.R"))
source(file.path(project_root, "R", "generate_data.R"))

set.seed(42)

# Test 1: mixed M_f, M_s
dat <- simulate_multi_study_data(
  S = 2, n_s = c(20, 30), L_f = 3, L_s = 2,
  M_f = c(2, 3, 1),
  M_s = list(c(2, 2), c(3, 1)),
  n_obs = 30, sigma_eps = 0.05
)
tp <- dat$true_params

cat("=== Test 1: mixed M_f/M_s ===\n")
cat("L_f =", tp$L_f, " L_s =", tp$L_s, "\n")
cat("M_f =", tp$M_f, "\n")
cat("M_s[[1]] =", tp$M_s[[1]], " M_s[[2]] =", tp$M_s[[2]], "\n\n")

# Check nu_phi: expected K+2 x M_f[l]
for (l in 1:3) {
  d <- dim(tp$nu_phi_true[[l]])
  stopifnot(d[1] == tp$K + 2)
  stopifnot(d[2] == tp$M_f[l])
}
cat("nu_phi_true dimensions: OK\n")

# Check nu_psi
for (s in 1:2) for (l in 1:2) {
  d <- dim(tp$nu_psi_true[[s]][[l]])
  stopifnot(d[1] == tp$K + 2)
  stopifnot(d[2] == tp$M_s[[s]][l])
}
cat("nu_psi_true dimensions: OK\n")

# Check zeta
for (s in 1:2) for (l in 1:3) {
  d <- dim(tp$zeta_true[[s]][[l]])
  stopifnot(d[1] == tp$n_s[s])
  stopifnot(d[2] == tp$M_f[l])
}
cat("zeta_true dimensions: OK\n")

# Check xi
for (s in 1:2) for (l in 1:2) {
  d <- dim(tp$xi_true[[s]][[l]])
  stopifnot(d[1] == tp$n_s[s])
  stopifnot(d[2] == tp$M_s[[s]][l])
}
cat("xi_true dimensions: OK\n")

# Check Y
stopifnot(length(dat$Y[[1]][[1]][[1]]) == 30)
cat("Y output: OK\n\n")

# Test 2: L_f = 0, L_s = 0 (no factors)
dat0 <- simulate_multi_study_data(S = 1, n_s = c(15), L_f = 0, L_s = 0, n_obs = 20)
stopifnot(dat0$true_params$L_f == 0)
stopifnot(dat0$true_params$L_s == 0)
stopifnot(length(dat0$Y[[1]][[1]][[1]]) == 20)
cat("Test 2 (L_f=0, L_s=0): PASS\n")

# Test 3: default M_f, M_s
dat_def <- simulate_multi_study_data(S = 2, n_s = c(10, 15), L_f = 3, L_s = 2, n_obs = 15)
stopifnot(all(dat_def$true_params$M_f == c(2, 2, 2)))
stopifnot(all(dat_def$true_params$M_s[[1]] == c(2, 2)))
cat("Test 3 (defaults): PASS\n")

# Test 4: L_s = 0 (no specific factors)
dat_nospec <- simulate_multi_study_data(S = 2, n_s = c(10, 15), L_f = 2, L_s = 0, n_obs = 20)
stopifnot(dat_nospec$true_params$L_s == 0)
stopifnot(is.null(dat_nospec$true_params$nu_psi_true))
stopifnot(is.null(dat_nospec$true_params$xi_true))
stopifnot(is.null(dat_nospec$true_params$b_true))
cat("Test 4 (L_s=0): PASS\n")

cat("\n=== ALL TESTS PASSED ===\n")
