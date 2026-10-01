# A small synthetic state-contract test for the private Stage 7I initializer.
# It does not read Stage 7C observations, any simulation truth, or run CAVI.
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1L])
source(file.path(dirname(script_path),
                     "build_stage7i_uas_states_20260929_v1.R"))
set.seed(20260929)
K_total <- 11L
p <- 20L
n <- 12L
M <- 2L
C2 <- lapply(seq_len(n), function(i)
  matrix(rnorm((7L + i %% 3L) * K_total), ncol = K_total))
C_g <- matrix(rnorm(51L * K_total), nrow = 51L)
time_g <- seq(0, 1, length.out = 51L)
Psi_true <- matrix(rnorm(K_total * M), K_total, M)
Xi_true <- matrix(rnorm(n * M), n, M)
b_true <- c(3, -2, 1.5, rep(0, p - 3L))
Y2 <- lapply(seq_len(n), function(i) {
  signal <- tcrossprod(as.vector(C2[[i]] %*% Psi_true %*%
                                  Xi_true[i, ]), b_true)
  noise <- matrix(rnorm(length(signal), sd = 0.05), nrow(signal))
  y <- signal + noise
  lapply(seq_len(p), function(j) y[, j])
})
observation <- list(
  Y = list(NULL, Y2),
  Z = list(matrix(0, n, 1L), matrix(0, n, 1L))
)
fit <- list(
  S = 2L, L_f = 1L, L_s = c(1L, 1L),
  M_s = list(2L, 2L), K = 9L, K_total = K_total,
  p = p, n_s = c(n, n), d = 1L,
  C_g = C_g, time_g = time_g,
  list_hyper = list(Sigma_beta = diag(1e10, 2L)),
  list_cp_C = list(NULL, lapply(C2, crossprod))
)
old_m <- matrix(rnorm(p, sd = 0.001), p, 1L)
old_rho <- matrix(0.22, p, 1L)
old_V <- matrix(1, p, 1L)
snapshot <- list(
  iteration = 50L, annealing = TRUE, temperature = 1.3,
  mu_q_nu_mu = list(NULL, lapply(seq_len(p), function(j)
    numeric(K_total))),
  mu_q_nu_beta = lapply(seq_len(p), function(j)
    list(numeric(K_total))),
  mu_q_nu_phi = list(matrix(0, K_total, M)),
  mu_q_zeta = list(list(matrix(0, n, M)),
                   list(matrix(0, n, M))),
  mu_q_a = matrix(0, p, 1L),
  mu_q_gamma_b = list(old_rho, old_rho),
  mu_q_normal_b = list(old_m, old_m),
  Sigma_q_normal_b = list(old_V, old_V),
  mu_q_b_specific = list(old_rho * old_m, old_rho * old_m),
  term_b_specific = list(old_rho * (old_V + old_m^2),
                         old_rho * (old_V + old_m^2)),
  mu_q_nu_psi = list(list(matrix(rnorm(K_total * M), K_total, M)),
                     list(matrix(rnorm(K_total * M), K_total, M))),
  mu_q_xi = list(list(matrix(rnorm(n * M, sd = 1e-4), n, M)),
                 list(matrix(rnorm(n * M, sd = 1e-4), n, M))),
  Sigma_q_nu_psi = list(list(replicate(M, diag(K_total), FALSE)),
                        list(replicate(M, diag(K_total), FALSE))),
  Sigma_q_xi = list(list(replicate(n, diag(M), FALSE)),
                    list(replicate(n, diag(M), FALSE))),
  mu_q_recip_sigsq_eps = matrix(100, 2L, p)
)
out <- stage7i_build_uas_states(observation, fit, snapshot, C2)
stopifnot(identical(out$U, snapshot))
stopifnot(isTRUE(out$blind_diagnostics$amplitudes_matched))
stopifnot(all(is.finite(out$blind_diagnostics$initial_energy)))
stopifnot(identical(out$A$mu_q_gamma_b, snapshot$mu_q_gamma_b))
stopifnot(identical(out$S$Sigma_q_nu_psi, snapshot$Sigma_q_nu_psi))
cat("STAGE7I_SYNTHETIC_STATE_CONTRACT_PASS\n")
