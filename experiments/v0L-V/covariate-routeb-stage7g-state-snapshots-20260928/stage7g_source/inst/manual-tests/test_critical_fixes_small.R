options(warn = 1)

working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(working_dir, "DESCRIPTION")) &&
    file.exists(file.path(working_dir, "R", "update_variance.R"))) {
  root <- working_dir
} else {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot locate the package root.")
  this_file <- normalizePath(sub("^--file=", "", file_arg),
                             winslash = "/", mustWork = TRUE)
  root <- dirname(this_file)
  while (!file.exists(file.path(root, "DESCRIPTION"))) {
    parent_root <- dirname(root)
    if (identical(parent_root, root)) stop("Could not locate the package root.")
    root <- parent_root
  }
}

create_named_list <- function(...) {
  nms <- names(match.call()[-1])
  if (is.null(nms)) {
    nms <- as.character(match.call()[-1])
  } else {
    no_name <- nms == ""
    if (any(no_name)) nms[no_name] <- as.character(match.call()[-1])[no_name]
  }
  setNames(list(...), nms)
}
tr <- function(x) sum(diag(x))
blkdiag <- function(...) {
  mats <- list(...)
  nr <- sum(vapply(mats, nrow, integer(1)))
  nc <- sum(vapply(mats, ncol, integer(1)))
  out <- matrix(0, nr, nc)
  ir <- ic <- 1L
  for (a in mats) {
    rr <- ir:(ir + nrow(a) - 1L)
    cc <- ic:(ic + ncol(a) - 1L)
    out[rr, cc] <- a
    ir <- max(rr) + 1L
    ic <- max(cc) + 1L
  }
  out
}

source(file.path(root, "R", "update_variance.R"))
source(file.path(root, "R", "update_loadings.R"))
source(file.path(root, "R", "update_eigenfunctions.R"))
source(file.path(root, "R", "orthonormalise_multi.R"))
source(file.path(root, "R", "update_omega.R"))
source(file.path(root, "R", "multi_core.R"))

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
near <- function(x, y, tol = 1e-8) isTRUE(all.equal(x, y, tolerance = tol))
pass <- function(message) cat("PASS", message, "\n")

set.seed(17)
K <- 1L
K_total <- 3L
S <- n <- 1L
p <- 2L
Lf <- Ls <- 1L
C1 <- matrix(c(1, 0, 0.4,
               1, 1, 0.7), nrow = 2, byrow = TRUE)
cp1 <- crossprod(C1)
C <- list(list(C1))
list_cp_C <- list(list(cp1))
Y <- list(list(list(c(0.2, -0.1), c(-0.3, 0.4))))
list_cp_C_Y <- list(list(sapply(Y[[1]][[1]], function(y) crossprod(C1, y))))
mu_mu <- list(list(rep(0, K_total), rep(0, K_total)))
Sigma_mu <- list(list(matrix(0, K_total, K_total),
                      matrix(0, K_total, K_total)))
mu_phi <- list(matrix(c(0.3, -0.2, 0.5), K_total, 1))
Sigma_phi_zero <- list(list(matrix(0, K_total, K_total)))
Sigma_phi_pos <- list(list(diag(c(0.05, 0.02, 0.08))))
mu_psi <- list(list(matrix(c(-0.1, 0.4, 0.2), K_total, 1)))
Sigma_psi_zero <- list(list(list(matrix(0, K_total, K_total))))
Sigma_psi_pos <- list(list(list(diag(c(0.03, 0.04, 0.02)))))
mu_eta <- list(list(matrix(0.6, 1, 1)))
Sigma_eta <- list(list(list(matrix(0.4, 1, 1))))
mu_chi <- list(list(matrix(-0.3, 1, 1)))
Sigma_chi <- list(list(list(matrix(0.25, 1, 1))))
mu_a <- matrix(c(0.4, -0.2), p, Lf)
term_a <- matrix(c(0.5, 0.3), p, Lf)
mu_b <- list(matrix(c(-0.3, 0.5), p, Ls))
term_b <- list(matrix(c(0.4, 0.6), p, Ls))

# 1) RSS interface and complete posterior moments.
mu_beta <- list(list(rep(0, K_total)), list(rep(0, K_total)))
Sigma_beta_zero <- list(list(matrix(0, K_total, K_total)),
                        list(matrix(0, K_total, K_total)))
Sigma_beta_pos <- list(list(diag(c(0.2, 0.1, 0.3))),
                       list(diag(c(0.1, 0.2, 0.1))))
Z <- list(matrix(2, 1, 1))

rss_args <- list(
  Y = Y, C = C, list_cp_C = list_cp_C,
  mu_q_nu_mu = mu_mu, Sigma_q_nu_mu = Sigma_mu,
  mu_q_nu_beta = mu_beta, Z = Z,
  mu_q_zeta = mu_eta, Sigma_q_zeta = Sigma_eta,
  mu_q_nu_phi = mu_phi,
  mu_q_xi = mu_chi, Sigma_q_xi = Sigma_chi,
  mu_q_nu_psi = mu_psi,
  mu_q_a = mu_a, term_a = term_a,
  mu_q_b_specific = mu_b, term_b_specific = term_b,
  mu_q_recip_a_eps = matrix(1, S, p),
  S = S, n_s = n, p = p, L_f = Lf, L_s = Ls,
  total_obs_sj = matrix(2, S, p), c_val = 1, n_cpus = 1)

r0 <- do.call(update_sigsq_eps, c(rss_args, list(
  Sigma_q_nu_beta = Sigma_beta_zero,
  Sigma_q_nu_phi = Sigma_phi_zero,
  Sigma_q_nu_psi = Sigma_psi_zero)))
r_beta <- do.call(update_sigsq_eps, c(rss_args, list(
  Sigma_q_nu_beta = Sigma_beta_pos,
  Sigma_q_nu_phi = Sigma_phi_zero,
  Sigma_q_nu_psi = Sigma_psi_zero)))
r_full <- do.call(update_sigsq_eps, c(rss_args, list(
  Sigma_q_nu_beta = Sigma_beta_pos,
  Sigma_q_nu_phi = Sigma_phi_pos,
  Sigma_q_nu_psi = Sigma_psi_pos)))
assert(all(r_beta$lambda_q_sigsq_eps > r0$lambda_q_sigsq_eps),
       "Beta posterior covariance was not added to RSS.")
assert(all(r_full$lambda_q_sigsq_eps > r_beta$lambda_q_sigsq_eps),
       "Time-function posterior covariance was not added to RSS.")
pass("RSS receives beta covariance and uses theta/kappa covariance")

# 2) Loading residuals must not recycle p-vectors against K_total-vectors.
warnings <- character()
capture_warnings <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})

ra <- capture_warnings(update_a_loadings(
  Y, C, list_cp_C, list_cp_C_Y,
  mu_mu, NULL, NULL,
  mu_eta, mu_phi, mu_chi, mu_psi,
  mu_a, term_a, matrix(0.5, p, Lf),
  matrix(0, p, Lf), matrix(1, p, Lf), mu_b,
  matrix(1, S, p), 0, 0,
  array(1, dim = c(Lf, S, n)),
  S, n, p, 0, Lf, Ls, K_total))
rb <- capture_warnings(update_b_loadings(
  Y, C, list_cp_C, list_cp_C_Y,
  mu_mu, NULL, NULL,
  mu_eta, mu_phi, mu_chi, mu_psi,
  mu_a, mu_b, term_b, list(matrix(0.5, p, Ls)),
  list(matrix(0, p, Ls)), list(matrix(1, p, Ls)),
  matrix(1, S, p), list(0), list(0),
  array(1, dim = c(S, Ls, n)),
  S, n, p, 0, Lf, Ls, K_total))
assert(length(warnings) == 0L, paste("Unexpected warning:", warnings[1]))
assert(all(is.finite(ra$mu_q_a)) && all(is.finite(rb$mu_q_b_specific[[1]])),
       "Loading update returned non-finite values.")
pass("shared/specific loading residuals use K_total by p outer products")

# 3) lambda_orth defaults to zero and is safe for L_f=1, L_s=1.
phi_args <- list(
  Y = Y, C = C, list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
  mu_q_nu_mu = mu_mu, Sigma_q_nu_mu = Sigma_mu,
  mu_q_nu_beta = NULL, Z = NULL,
  mu_q_zeta = mu_eta, Sigma_q_zeta = Sigma_eta,
  mu_q_nu_phi = mu_phi, mu_q_xi = mu_chi, mu_q_nu_psi = mu_psi,
  mu_q_a = mu_a, term_a = term_a,
  mu_q_b_specific = mu_b, term_b_specific = term_b,
  mu_q_recip_sigsq_eps = matrix(1, S, p),
  mu_q_recip_sigsq_phi = list(1), inv_Sigma_beta = diag(2),
  S = S, n_s = n, p = p, d = 0, L_f = Lf, L_s = Ls,
  M_f = 1, K = K, K_total = K_total, c_val = 1, n_cpus = 1)
phi_default <- do.call(update_nu_phi, phi_args)
phi_zero <- do.call(update_nu_phi, c(phi_args, list(lambda_orth = 0)))
phi_pos <- do.call(update_nu_phi, c(phi_args, list(lambda_orth = 0.1)))
assert(near(phi_default$Sigma_q_nu_phi, phi_zero$Sigma_q_nu_phi),
       "Default orthogonal penalty is not zero.")
assert(!near(phi_zero$Sigma_q_nu_phi, phi_pos$Sigma_q_nu_phi),
       "Positive optional orthogonal penalty had no effect.")
psi_args <- list(
  Y = Y, C = C, list_cp_C = list_cp_C, list_cp_C_Y = list_cp_C_Y,
  mu_q_nu_mu = mu_mu, Sigma_q_nu_mu = Sigma_mu,
  mu_q_nu_beta = NULL, Z = NULL,
  mu_q_zeta = mu_eta, mu_q_nu_phi = mu_phi,
  mu_q_xi = mu_chi, Sigma_q_xi = Sigma_chi,
  mu_q_nu_psi = mu_psi, mu_q_a = mu_a, term_a = term_a,
  mu_q_b_specific = mu_b, term_b_specific = term_b,
  mu_q_recip_sigsq_eps = matrix(1, S, p),
  mu_q_recip_sigsq_psi = list(list(1)), inv_Sigma_beta = diag(2),
  S = S, n_s = n, p = p, d = 0, L_f = Lf, L_s = Ls,
  M_s = list(1), K = K, K_total = K_total, c_val = 1, n_cpus = 1)
psi_default <- do.call(update_nu_psi, psi_args)
psi_zero <- do.call(update_nu_psi, c(psi_args, list(lambda_orth = 0)))
psi_pos <- do.call(update_nu_psi, c(psi_args, list(lambda_orth = 0.1)))
assert(near(psi_default$Sigma_q_nu_psi, psi_zero$Sigma_q_nu_psi),
       "Specific default orthogonal penalty is not zero.")
assert(!near(psi_zero$Sigma_q_nu_psi, psi_pos$Sigma_q_nu_psi),
       "Positive specific orthogonal penalty had no effect.")
pass("shared/specific soft regulariser is explicit, defaults to zero, and is branch-safe")

# 4) Posterior-moment FPCA, scaling, sorting and loading-based sign.
time_g <- seq(0, 1, length.out = 17)
C_g <- cbind(1, time_g, sin(2 * pi * time_g))
mu_phi2 <- list(
  matrix(c(0.2, 0.1, 0.8), K_total, 1),
  matrix(c(-0.1, 0.5, 0.2), K_total, 1))
Sigma_phi2 <- list(
  list(diag(c(0.03, 0.02, 0.04))),
  list(diag(c(0.02, 0.01, 0.03))))
eta2 <- list(list(matrix(c(-0.6, 0.1, 0.8), 3, 1),
                    matrix(c(0.2, -0.3, 0.4), 3, 1)))
eta_cov2 <- list(list(
  replicate(3, matrix(0.2, 1, 1), simplify = FALSE),
  replicate(3, matrix(0.3, 1, 1), simplify = FALSE)))
A2 <- matrix(c(-0.05, 0.02, 0.01,
               -1.2, 0.3, 0.1), 3, 2)
gamma2 <- matrix(c(0.5, 0.4, 0.3, 0.9, 0.8, 0.7), 3, 2)

pp <- orthonormalise_multi(
  C_g, time_g,
  mu_q_nu_mu = list(list(rep(0, K_total), rep(0, K_total), rep(0, K_total))),
  mu_q_nu_phi = mu_phi2, mu_q_nu_psi = NULL,
  mu_q_zeta = eta2, Sigma_q_zeta = eta_cov2,
  mu_q_xi = NULL, Sigma_q_xi = NULL,
  mu_q_a = A2, mu_q_b_specific = NULL,
  mu_q_gamma_a = gamma2, mu_q_gamma_b = NULL,
  S = 1, n_s = 3, p = 3, L_f = 2, L_s = 0,
  M_f = c(1, 1), M_s = list(numeric()),
  Sigma_q_nu_phi = Sigma_phi2, Sigma_q_nu_psi = NULL)

w <- .trap_weights(time_g)
for (Phi in pp$list_Phi_hat) {
  assert(near(crossprod(Phi, w * Phi), diag(ncol(Phi)), 1e-6),
         "Weighted FPCA functions are not orthonormal.")
}
assert(all(diff(colSums(pp$mu_q_a_hat^2)) <= 1e-12),
       "Shared factors were not sorted by loading column norm.")
assert(all(vapply(seq_len(ncol(pp$mu_q_a_hat)), function(l) {
  x <- pp$mu_q_a_hat[, l]
  x[which.max(abs(x))] >= 0
}, logical(1))), "Loading sign convention was not applied.")
assert(all(vapply(pp$list_eigenvalues, sum, numeric(1)) > 0.98),
       "Retained normalized eigenvalues do not explain at least 99 percent up to tolerance.")

pp_zero_cov <- orthonormalise_multi(
  C_g, time_g,
  mu_q_nu_mu = list(list(rep(0, K_total), rep(0, K_total), rep(0, K_total))),
  mu_q_nu_phi = mu_phi2, mu_q_nu_psi = NULL,
  mu_q_zeta = eta2, Sigma_q_zeta = eta_cov2,
  mu_q_xi = NULL, Sigma_q_xi = NULL,
  mu_q_a = A2, mu_q_b_specific = NULL,
  mu_q_gamma_a = gamma2, mu_q_gamma_b = NULL,
  S = 1, n_s = 3, p = 3, L_f = 2, L_s = 0,
  M_f = c(1, 1), M_s = list(numeric()),
  Sigma_q_nu_phi = lapply(mu_phi2, .zero_coef_cov),
  Sigma_q_nu_psi = NULL)
assert(!near(sort(pp$factor_scale_shared), sort(pp_zero_cov$factor_scale_shared)),
       "Time-function covariance did not affect posterior integrated variance.")
pass("posterior-moment FPCA uses weighted eigendecomposition, scale, order and loading sign")

pp_spec <- orthonormalise_multi(
  C_g, time_g,
  mu_q_nu_mu = list(list(rep(0, K_total), rep(0, K_total), rep(0, K_total))),
  mu_q_nu_phi = mu_phi2[1],
  mu_q_nu_psi = list(list(matrix(c(0.1, -0.4, 0.7), K_total, 1))),
  mu_q_zeta = list(list(matrix(c(-0.6, 0.1, 0.8), 3, 1))),
  Sigma_q_zeta = list(list(replicate(3, matrix(0.2, 1, 1), simplify = FALSE))),
  mu_q_xi = list(list(matrix(c(0.4, -0.2, 0.6), 3, 1))),
  Sigma_q_xi = list(list(replicate(3, matrix(0.15, 1, 1), simplify = FALSE))),
  mu_q_a = matrix(c(0.4, -0.2, 0.1), 3, 1),
  mu_q_b_specific = list(matrix(c(-0.8, 0.1, 0.2), 3, 1)),
  mu_q_gamma_a = matrix(0.5, 3, 1),
  mu_q_gamma_b = list(matrix(0.6, 3, 1)),
  S = 1, n_s = 3, p = 3, L_f = 1, L_s = 1,
  M_f = 1, M_s = list(1),
  Sigma_q_nu_phi = Sigma_phi2[1],
  Sigma_q_nu_psi = list(list(list(diag(c(0.02, 0.03, 0.01))))))
assert(near(crossprod(pp_spec$list_Phi_hat_spec[[1]][[1]],
                      w * pp_spec$list_Phi_hat_spec[[1]][[1]]),
            diag(ncol(pp_spec$list_Phi_hat_spec[[1]][[1]])), 1e-6),
       "Specific posterior-moment FPCA is not weighted-orthonormal.")
bx <- pp_spec$mu_q_b_specific_hat[[1]][, 1]
assert(bx[which.max(abs(bx))] > 0,
       "Specific loading sign convention was not applied.")
pass("specific-factor posterior-moment FPCA uses the same scale and sign rules")

# 6) Both omega hierarchies must have the documented dimensions and TRUE must
# use variable rows rather than factor indices or vector recycling.
gamma_shared <- matrix(c(0.1, 0.8, 0.3, 0.9,
                         0.7, 0.2, 0.6, 0.4), nrow = 4, ncol = 2)
gamma_specific <- list(gamma_shared, 1 - gamma_shared)
omega_factor <- update_omega(
  mu_q_gamma_a = gamma_shared, mu_q_gamma_b = gamma_specific,
  p = 4, S = 2, c_0 = 1, d_0 = 4,
  bool_var_spec_prob = FALSE)
omega_variable <- update_omega(
  mu_q_gamma_a = gamma_shared, mu_q_gamma_b = gamma_specific,
  p = 4, S = 2, c_0 = 1, d_0 = 4,
  bool_var_spec_prob = TRUE)
assert(identical(dim(omega_variable$c_1_omega_a), c(4L, 2L)) &&
         length(omega_factor$c_1_omega_a) == 2L,
       "Shared omega hierarchy returned the wrong dimensions.")
assert(all(vapply(omega_variable$c_1_omega_b,
                  function(x) identical(dim(x), c(4L, 2L)), logical(1))) &&
         all(vapply(omega_factor$c_1_omega_b,
                    function(x) length(x) == 2L, logical(1))),
       "Specific omega hierarchy returned the wrong dimensions.")
assert(near(omega_variable$c_1_omega_a, 1 + gamma_shared) &&
         near(omega_variable$d_1_omega_a, 5 - gamma_shared),
       "Variable-level Beta update does not match its conjugate formula.")

shared_prior_rows <- cbind(c(-5, -1, 1, 5), c(2, 2, 2, 2))
selected <- .loading_log_prior_odds(
  shared_prior_rows, matrix(0, 4, 2), factor_index = 1,
  p = 4, label = "test")
assert(near(selected, c(-5, -1, 1, 5)),
       "Variable-level loading prior was not selected by [j,l].")

ra_variable <- update_a_loadings(
  Y, C, list_cp_C, list_cp_C_Y,
  mu_mu, NULL, NULL,
  mu_eta, mu_phi, mu_chi, mu_psi,
  mu_a, term_a, matrix(0.5, p, Lf),
  matrix(0, p, Lf), matrix(1, p, Lf), mu_b,
  matrix(1, S, p), matrix(c(-6, 6), p, Lf), matrix(0, p, Lf),
  array(1, dim = c(Lf, S, n)),
  S, n, p, 0, Lf, Ls, K_total)
assert(ra_variable$mu_q_gamma_a[2, 1] > ra_variable$mu_q_gamma_a[1, 1],
       "TRUE shared loading update ignored variable-specific omega odds.")

rb_variable <- update_b_loadings(
  Y, C, list_cp_C, list_cp_C_Y,
  mu_mu, NULL, NULL,
  mu_eta, mu_phi, mu_chi, mu_psi,
  mu_a, mu_b, term_b, list(matrix(0.5, p, Ls)),
  list(matrix(0, p, Ls)), list(matrix(1, p, Ls)),
  matrix(1, S, p), list(matrix(c(-6, 6), p, Ls)),
  list(matrix(0, p, Ls)), array(1, dim = c(S, Ls, n)),
  S, n, p, 0, Lf, Ls, K_total)
assert(rb_variable$mu_q_gamma_b[[1]][2, 1] >
         rb_variable$mu_q_gamma_b[[1]][1, 1],
       "TRUE specific loading update ignored variable-specific omega odds.")
pass("factor- and variable-level omega branches use correct dimensions and [j,l] indexing")

assert(identical(eval(formals(bayesSYNC_multi)$lambda_orth), 0),
       "Public lambda_orth default is not zero.")
pc0 <- .parameter_change(list(matrix(1:4, 2)), list(matrix(1:4, 2)))
pc1 <- .parameter_change(list(matrix(1:4, 2)), list(matrix(c(2:5), 2)))
assert(near(pc0, c(abs = 0, rel = 0)) && pc1[["abs"]] > 0,
       "Parameter-change convergence metric is not functioning.")
pass("default stopping support uses a tested parameter-change metric")

cat("ALL PASS elapsed small-sample critical-fix tests\n")
