arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript run_transition_micro_smoke_20260909.R <ASCII package root>")
}
package_root <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for this development smoke.")
}
if (!requireNamespace("Splinets", quietly = TRUE)) {
  stop("Package 'Splinets' is required for this development smoke.")
}
pkgload::load_all(package_root, quiet = TRUE, export_all = FALSE,
                  helpers = FALSE)

started <- proc.time()[["elapsed"]]
data <- simulate_multi_study_structured(
  S = 2,
  n_s = c(4, 4),
  p = 20,
  d = 0,
  L_f = 1,
  L_s = c(1, 1),
  M_f = 2,
  M_s = list(2, 2),
  K = 5,
  n_obs = 8,
  common_grid = FALSE,
  n_obs_range = c(6, 9),
  sigma_eps = 0.3,
  mean_amp = 0.6,
  score_var_decay = TRUE,
  n_dense = 201,
  identified_loadings = TRUE,
  seed = 9060901,
  truth_basis_method = "jaoua_splinets"
)
elapsed <- proc.time()[["elapsed"]] - started
truth <- data$true_params
metadata <- c(
  truth$phi_basis_meta,
  unlist(truth$psi_basis_meta, recursive = FALSE)
)
continuous_error <- max(vapply(metadata, function(value) {
  max(abs(value$gram - diag(nrow(value$gram))))
}, numeric(1)))
grid_error <- max(vapply(metadata, function(value) {
  max(abs(value$grid_gram - diag(nrow(value$grid_gram))))
}, numeric(1)))
response_error <- 0
for (s in seq_len(truth$S)) {
  for (i in seq_len(truth$n_s[s])) {
    observed <- do.call(cbind, data$Y[[s]][[i]])
    response_error <- max(
      response_error,
      abs(observed - truth$signal_true_values[[s]][[i]] -
            truth$noise_true_values[[s]][[i]])
    )
  }
}

stopifnot(
  identical(truth$truth_basis_method, "jaoua_splinets"),
  identical(truth$truth_basis_family,
            "jaoua_splinets_orthonormal_B_spline"),
  length(metadata) == 3L,
  continuous_error < 1e-12,
  response_error < 1e-13,
  all(unlist(truth$n_obs_by_subject, use.names = FALSE) >= 6L),
  all(unlist(truth$n_obs_by_subject, use.names = FALSE) <= 9L)
)

cat("GENERATOR_SMOKE_PASS\n")
cat("formal_data_generated=FALSE\n")
cat("fit_started=FALSE\n")
cat("factor_blocks=", length(metadata), "\n", sep = "")
cat("continuous_gram_max_abs_error=", format(continuous_error, digits = 17),
    "\n", sep = "")
cat("dense_grid_gram_max_abs_error=", format(grid_error, digits = 17),
    "\n", sep = "")
cat("response_identity_max_abs_error=", format(response_error, digits = 17),
    "\n", sep = "")
cat("Splinets_version=", as.character(packageVersion("Splinets")), "\n",
    sep = "")
cat("elapsed_seconds=", format(elapsed, digits = 6), "\n", sep = "")
