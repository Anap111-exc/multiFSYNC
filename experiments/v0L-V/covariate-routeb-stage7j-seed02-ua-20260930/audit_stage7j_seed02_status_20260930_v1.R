#!/usr/bin/env Rscript

# Read-only recomputation of Stage 7C objective and practical-stop semantics.
# Usage: Rscript SCRIPT RESULTS_PARENT
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected RESULTS_PARENT.", call. = FALSE)
parent <- normalizePath(args[[1L]], mustWork = TRUE)
stopifnot(requireNamespace("digest", quietly = TRUE))
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
out <- file.path(parent, "stage7j_seed02_status_semantics_20260930_v1")
assert(!file.exists(out) && !dir.exists(out), "Output already exists.")

freeze <- file.path(parent, "stage7j_seed02_blind_freeze_20260930_v1")
marker <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
lines <- readLines(marker, warn = FALSE)
assert(identical(lines[[1L]], "STAGE7J_BLIND_ENDPOINTS_FROZEN") &&
         "branches=U,A" %in% lines, "Blind freeze changed.")
manifest_path <- file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv")
assert(paste0("manifest_sha256=", sha(manifest_path)) %in% lines,
       "Blind manifest changed.")
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
assert(identical(manifest$branch, c("U", "A")), "Branch set changed.")

rows <- lapply(seq_len(nrow(manifest)), function(i) {
  path <- manifest$fit_path[[i]]
  assert(identical(sha(path), manifest$fit_sha256[[i]]),
         paste("Fit changed", manifest$branch[[i]]))
  fit <- readRDS(path)
  finite_elbo <- length(fit$ELBO) > 0L && all(is.finite(fit$ELBO))
  finite_components <- length(fit$ELBO_finite) == length(fit$ELBO) &&
    all(vapply(fit$ELBO_finite, function(value)
      is.list(value) && isTRUE(value$all), logical(1L)))
  checks <- c(
    annealing_completed = isTRUE(fit$annealing_completed),
    final_temperature_one = isTRUE(all.equal(
      as.numeric(fit$final_temperature), 1)),
    ordinary_elbo_present_and_finite = finite_elbo,
    ordinary_elbo_components_finite = finite_components,
    ordinary_objective_valid = isTRUE(fit$elbo_objective_valid),
    no_T1_elbo_decrease = identical(as.integer(fit$ELBO_decrease_count), 0L),
    no_T1_solver_jitter = identical(as.integer(fit$elbo_t1_jitter_count), 0L))
  strict <- isTRUE(fit$practical_converged) && isTRUE(fit$converged) &&
    identical(as.character(fit$convergence_status), "converged_practical")
  data.frame(branch = manifest$branch[[i]], fit_id = manifest$fit_id[[i]],
             fit_sha256 = manifest$fit_sha256[[i]],
             objective_eligible = all(checks),
             invalid_reasons = paste(names(checks)[!checks], collapse = ";"),
             annealing_completed = checks[["annealing_completed"]],
             final_temperature_one = checks[["final_temperature_one"]],
             ordinary_elbo_present_and_finite =
               checks[["ordinary_elbo_present_and_finite"]],
             ordinary_elbo_components_finite =
               checks[["ordinary_elbo_components_finite"]],
             ordinary_objective_valid = checks[["ordinary_objective_valid"]],
             no_T1_elbo_decrease = checks[["no_T1_elbo_decrease"]],
             no_T1_solver_jitter = checks[["no_T1_solver_jitter"]],
             strict_practical_converged = strict,
             practical_converged = isTRUE(fit$practical_converged),
             converged = isTRUE(fit$converged),
             convergence_status = as.character(fit$convergence_status),
             ordinary_sweeps = fit$t1_sweeps,
             final_ordinary_elbo = tail(fit$ELBO, 1L),
             stringsAsFactors = FALSE)
})
status <- do.call(rbind, rows)
baseline_path <- file.path(parent, "covariate_routeb_stage7c_20260927_v1",
                           "scientific_evaluation_20260928_v1",
                           "ALL_START_SCIENCE.csv")
assert(identical(sha(baseline_path),
                 "fcb9a4db7a8ae47062f8c73435e737b69448695b245e4050e9a8b2ef74554081"),
       "Stage 7C baseline table changed.")
baseline <- utils::read.csv(baseline_path, stringsAsFactors = FALSE)
baseline <- baseline[baseline$fit_id == "covrb7c_02__C1__start11", , drop = FALSE]
assert(nrow(baseline) == 1L &&
         identical(status$objective_eligible[[1L]],
                   baseline$objective_eligible[[1L]]) &&
         identical(status$strict_practical_converged[[1L]],
                   baseline$strict_practical_converged[[1L]]),
       "U status semantics differ from Stage 7C baseline.")
assert(dir.create(out, mode = "0700"), "Cannot create output.")
path <- file.path(out, "OBJECTIVE_AND_STOP_STATUS.csv")
utils::write.csv(status, path, row.names = FALSE, quote = TRUE, na = "")
writeLines(c("STAGE7J_STATUS_SEMANTICS_COMPLETE",
             "branches=U,A", "u_stage7c_status_parity=PASS",
             "truth_used_for_status=FALSE", "fit_started=FALSE",
             paste0("table_sha256=", sha(path))),
           file.path(out, "STATUS_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7J_STATUS_SEMANTICS_COMPLETE\n")
