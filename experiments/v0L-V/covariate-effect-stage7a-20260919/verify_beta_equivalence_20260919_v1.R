#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
if (length(args) != 2L || !args[[1L]] %in% c("baseline", "compare")) {
  stop("Usage: Rscript verify_beta_equivalence_20260919_v1.R baseline|compare output.rds")
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                           mustWork = TRUE)
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)

set.seed(70701L)
S <- 2L
n_s <- c(4L, 4L)
p <- 6L
time_obs <- lapply(seq_len(S), function(s) {
  lapply(seq_len(n_s[[s]]), function(i) sort(runif(7L)))
})
Z <- lapply(seq_len(S), function(s) {
  matrix(seq(-1, 1, length.out = n_s[[s]]), ncol = 1L)
})
Y <- lapply(seq_len(S), function(s) {
  lapply(seq_len(n_s[[s]]), function(i) {
    t <- time_obs[[s]][[i]]
    lapply(seq_len(p), function(j) {
      0.15 * sin(2 * pi * t + j / 5) +
        0.22 * Z[[s]][i, 1L] * sin(pi * t + j / 7) +
        0.08 * cos(3 * pi * t + s / 3) +
        rnorm(length(t), sd = 0.03)
    })
  })
})

fit <- multiFSYNC::bayesSYNC_multi_pre_score(
  Y = Y, Z = Z, time_obs = time_obs,
  L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
  K = 3L, n_g = 31L, anneal = c(1, 1.2, 3), maxit = 7L,
  n_cpus = 1L, verbose = FALSE, seed = 70702L,
  bool_scale = FALSE, bool_var_spec_prob = FALSE,
  tol_abs = 0, tol_rel = 0, convergence_rule = "parameters",
  lambda_orth = 0, initialization = "random",
  pre_score_sweeps = 1L, trace_sweeps = 1L,
  function_initialization = "gram_unit_energy"
)
result <- list(
  mu_q_nu_beta = fit$mu_q_nu_beta,
  Sigma_q_nu_beta = fit$Sigma_q_nu_beta,
  mu_q_nu_mu = fit$mu_q_nu_mu,
  ELBO = fit$ELBO,
  ELBO_components = fit$ELBO_components,
  t1_sweeps = fit$t1_sweeps,
  annealing_sweeps = fit$annealing_sweeps
)
stopifnot(length(result$ELBO) > 0L, all(is.finite(result$ELBO)),
          !is.null(result$mu_q_nu_beta),
          !is.null(result$Sigma_q_nu_beta))

if (args[[1L]] == "baseline") {
  if (file.exists(args[[2L]])) stop("Refusing to overwrite the baseline.")
  dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, args[[2L]], version = 3)
  cat("BASELINE_PASS ordinary_t1_sweeps=", result$t1_sweeps,
      " baseline=", args[[2L]], "\n", sep = "")
} else {
  reference <- readRDS(args[[2L]])
  for (name in names(result)) {
    difference <- all.equal(result[[name]], reference[[name]],
                            tolerance = 1e-12, check.attributes = TRUE)
    if (!isTRUE(difference)) {
      stop(name, " differs from the pre-change baseline: ",
           paste(difference, collapse = "; "))
    }
    cat("EQUIVALENT ", name, "\n", sep = "")
  }
  cat("EQUIVALENCE_PASS\n")
}
