#!/usr/bin/env Rscript

# Read-only Stage 7I trajectory reduction. Uses only blind states and data.
# Usage: Rscript SCRIPT RESULTS_PARENT PACKAGE_LIB INITIALIZER_HELPER

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected RESULTS_PARENT PACKAGE_LIB INITIALIZER_HELPER.", call. = FALSE)
parent <- normalizePath(args[[1L]], mustWork = TRUE)
lib <- normalizePath(args[[2L]], mustWork = TRUE)
helper <- normalizePath(args[[3L]], mustWork = TRUE)
.libPaths(unique(c(lib, .libPaths())))
stopifnot(requireNamespace("multiFSYNC", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
assert(identical(sha(helper),
                 "03ae170103e2d70e63889b68553d9d13261233edfa1ec775ce185de39f840c99"),
       "Bound initializer helper changed.")
freeze <- file.path(parent, "stage7i_blind_freeze_20260929_v1")
marker <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
assert(file.exists(marker) &&
         identical(readLines(marker, n = 1L), "STAGE7I_BLIND_ENDPOINTS_FROZEN"),
       "Blind freeze is missing.")
manifest <- utils::read.csv(file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
assert(identical(manifest$branch, c("U", "A", "S")) &&
         all(manifest$status == "COMPLETE"), "Blind branch set changed.")
observation_path <- file.path(parent,
  "covariate_routeb_stage7c_20260927_v1", "data", "covrb7c_01__C1",
  "observation_bundle.rds")
assert(identical(sha(observation_path),
                 "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6"),
       "Bound observation bundle changed.")
observation <- readRDS(observation_path)
fit <- readRDS(file.path(parent, manifest$fit_path[[1L]]))
grid <- multiFSYNC:::get_grid_objects(
  unlist(observation$time_obs, recursive = FALSE), fit$K,
  n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
n1 <- fit$n_s[[1L]]
n2 <- fit$n_s[[2L]]
C2 <- grid$C[(n1 + 1L):(n1 + n2)]
N_obs <- sum(vapply(C2, nrow, integer(1L)))
p <- fit$p
assert(p == 500L && length(C2) == 30L, "Observation dimensions changed.")
rm(fit, grid)
invisible(gc())
helper_env <- new.env(parent = globalenv())
sys.source(helper, envir = helper_env)
state_energy <- helper_env$.stage7i_state_energy
as_root <- file.path(parent, "stage7i_as_20260929_v1")
out <- file.path(parent, "stage7i_trajectory_readout_20260929_v1")
assert(!file.exists(out) && !dir.exists(out), "Trajectory output exists.")
assert(dir.create(out, mode = "0700"), "Cannot create trajectory output.")

all_rows <- vector("list", 3L)
for (b in seq_len(3L)) {
  branch <- manifest$branch[[b]]
  snapshot_path <- file.path(parent, manifest$snapshot_path[[b]])
  assert(identical(sha(snapshot_path), manifest$snapshot_sha256[[b]]),
         paste0("Snapshot changed: ", branch))
  states <- readRDS(snapshot_path)
  initial <- readRDS(file.path(as_root, paste0("STATE_", branch,
                                             "_SWEEP50.rds")))
  states <- c(list(initial), states)
  assert(length(states) == 450L &&
           identical(vapply(states, `[[`, integer(1L), "iteration"),
                     as.integer(50:499)),
         paste0("State index changed: ", branch))
  branch_root <- if (branch == "U")
    file.path(parent, "stage7i_u_20260929_v1", "u_full_v1") else
    file.path(as_root, paste0(branch, "_full_v1"))
  driver <- utils::read.csv(file.path(branch_root, "DRIVER_TRACE.csv"),
                            stringsAsFactors = FALSE)
  driver <- driver[driver$block == "specific_loading" &
                     driver$study == 2L & driver$factor == 1L, , drop = FALSE]
  assert(nrow(driver) == 449L &&
           identical(as.integer(driver$iteration), as.integer(51:499)),
         paste0("B2 loading driver index changed: ", branch))
  rows <- lapply(seq_along(states), function(j) {
    state <- states[[j]]
    energy <- state_energy(state, C2, p, N_obs)
    rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
    rho_a <- as.numeric(state$mu_q_gamma_a)
    rho_b1 <- as.numeric(state$mu_q_gamma_b[[1L]])
    driver_row <- if (j == 1L) NULL else driver[j - 1L, , drop = FALSE]
    data.frame(
      branch = branch, iteration = state$iteration,
      annealing = isTRUE(state$annealing),
      temperature = state$temperature,
      b2_loading_norm = energy[["effective_loading_norm"]],
      b2_mean_time_product_norm = energy[["mean_time_product_norm"]],
      b2_mean_product_rms = energy[["mean_product_rms"]],
      b2_full_posterior_second_rms =
        energy[["full_posterior_product_second_rms"]],
      b2_score_rms_1 = energy[["score_rms_1"]],
      b2_score_rms_2 = energy[["score_rms_2"]],
      b2_ppi_mean = mean(rho), b2_ppi_ge_half = sum(rho >= 0.5),
      b1_ppi_ge_half = sum(rho_b1 >= 0.5),
      shared_ppi_ge_half = sum(rho_a >= 0.5),
      b2_loading_data_driver_norm = if (is.null(driver_row)) NA_real_ else
        driver_row$data_driver_norm,
      b2_loading_data_coherence = if (is.null(driver_row)) NA_real_ else
        driver_row$data_coherence,
      b2_loading_linear_driver_norm = if (is.null(driver_row)) NA_real_ else
        driver_row$linear_driver_norm,
      b2_loading_posterior_mean_norm = if (is.null(driver_row)) NA_real_ else
        driver_row$posterior_mean_norm,
      b2_loading_precision_trace = if (is.null(driver_row)) NA_real_ else
        driver_row$precision_trace,
      b2_loading_ppi_mean_after_update = if (is.null(driver_row)) NA_real_ else
        driver_row$ppi_mean,
      stringsAsFactors = FALSE)
  })
  all_rows[[b]] <- do.call(rbind, rows)
  rm(states, rows, driver)
  invisible(gc())
  cat("STAGE7I_TRAJECTORY_BRANCH_COMPLETE branch=", branch, "\n", sep = "")
}
trajectory <- do.call(rbind, all_rows)
utils::write.csv(trajectory, file.path(out, "STAGE7I_B2_TRAJECTORY.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
key_sweeps <- c(50:70, 75L, 80L, 90L, 100L, 110L, 125L, 150L,
                200L, 300L, 400L, 499L)
utils::write.csv(trajectory[trajectory$iteration %in% key_sweeps, ],
                 file.path(out, "STAGE7I_B2_KEY_SWEEPS.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
summary <- do.call(rbind, lapply(split(trajectory, trajectory$branch), function(x) {
  first <- function(flag) if (any(flag, na.rm = TRUE))
    min(x$iteration[which(flag)]) else NA_integer_
  data.frame(branch = x$branch[[1L]],
    b2_ppi_ge_half_first = first(x$b2_ppi_ge_half >= 1L),
    b2_ppi_40_first = first(x$b2_ppi_ge_half >= 40L),
    b2_mean_rms_005_first = first(x$b2_mean_product_rms >= 0.05),
    b2_mean_rms_010_first = first(x$b2_mean_product_rms >= 0.10),
    b2_ppi_ge_half_final = tail(x$b2_ppi_ge_half, 1L),
    b2_mean_product_rms_final = tail(x$b2_mean_product_rms, 1L),
    b2_full_posterior_second_rms_final =
      tail(x$b2_full_posterior_second_rms, 1L),
    stringsAsFactors = FALSE)
}))
summary <- summary[match(c("U", "A", "S"), summary$branch), , drop = FALSE]
utils::write.csv(summary, file.path(out, "STAGE7I_B2_TRAJECTORY_SUMMARY.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
writeLines(c("STAGE7I_BLIND_TRAJECTORY_READOUT_COMPLETE",
             "branches=U,A,S", "truth_opened=FALSE",
             paste0("freeze_marker_sha256=", sha(marker)),
             paste0("trajectory_sha256=", sha(file.path(
               out, "STAGE7I_B2_TRAJECTORY.csv")))),
           file.path(out, "TRAJECTORY_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7I_BLIND_TRAJECTORY_READOUT_COMPLETE\n")
