# Non-formal micro smoke for the supported G12 public interface.
# This script generates one tiny data set and runs two three-sweep fits. It is
# not a v0L-V data generation, formal fit, continuation, or scientific audit.

suppressPackageStartupMessages(library(multiFSYNC))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

data_seed <- 9861L
fit_seed <- 9862L
data <- multiFSYNC::simulate_multi_study_osullivan(
  S = 2L, n_s = c(3L, 3L), p = 4L, d = 0L,
  L_f = 1L, L_s = 1L,
  M_f = 1L, M_s = list(1L, 1L),
  K = 3L, n_obs = 7L, common_grid = TRUE,
  sigma_eps = 0.1, seed = data_seed, use_explicit = FALSE
)

fit_arguments <- list(
  Y = data$Y, Z = NULL, time_obs = data$time_obs,
  L_f = 1L, L_s = 1L, M_f = 1L, M_s = list(1L, 1L),
  K = 3L, n_g = 21L, anneal = NULL, maxit = 3L,
  tol_abs = 0, tol_rel = 0, convergence_rule = "parameters",
  bool_scale = FALSE, verbose = FALSE, seed = fit_seed,
  pre_score_sweeps = 1L, trace_sweeps = 1:2
)

current <- suppressWarnings(do.call(
  multiFSYNC::bayesSYNC_multi_pre_score,
  c(fit_arguments, list(function_initialization = "current_1_over_m"))
))
gram <- suppressWarnings(do.call(
  multiFSYNC::bayesSYNC_multi_pre_score,
  c(fit_arguments, list(function_initialization = "gram_unit_energy"))
))

assert(
  identical(current$pre_score_interface$function_initialization,
            "current_1_over_m"),
  "Current initialization provenance is missing."
)
assert(
  identical(current$driver_diagnostic_control$random_scale_calibration,
            "none"),
  "Current initialization did not preserve the inert calibration."
)
assert(
  identical(gram$pre_score_interface$function_initialization,
            "gram_unit_energy"),
  "G12 initialization provenance is missing."
)
assert(
  identical(gram$driver_diagnostic_control$random_scale_calibration,
            "function"),
  "G12 public option did not map to function-only Gram calibration."
)

diagnostics <- gram$random_scale_calibration_diagnostics
assert(nrow(diagnostics) == 3L, "G12 did not calibrate all three function blocks.")
assert(all(diagnostics$success), "A G12 function calibration failed.")
assert(
  all(diagnostics$block %in%
        c("shared_function_l2", "specific_function_l2")),
  "G12 calibrated a non-function state block."
)
assert(
  isTRUE(all.equal(diagnostics$after, rep(1, 3L), tolerance = 1e-12)),
  "G12 function energy is not one."
)
assert(
  all(is.finite(current$ELBO)) && all(is.finite(gram$ELBO)),
  "A micro fit returned a non-finite ELBO."
)
assert(
  any(current$driver_trace$phase == "pre_score") &&
    any(gram$driver_trace$phase == "pre_score"),
  "The bounded pre-score phase was not executed."
)

cat("G12_PUBLIC_INTERFACE_SMOKE_PASS\n")
cat("formal_experiment=FALSE\n")
cat("data_seed=", data_seed, "\n", sep = "")
cat("fit_seed=", fit_seed, "\n", sep = "")
cat("fits=2\n")
cat("continuation_started=FALSE\n")
