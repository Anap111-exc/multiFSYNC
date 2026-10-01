#!/usr/bin/env Rscript

# Read-only reduction of saved Stage 7J states. No fit and no truth access.
# Usage: Rscript SCRIPT RESULTS_PARENT PACKAGE_LIB INITIALIZER_HELPER
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("Expected RESULTS_PARENT PACKAGE_LIB INITIALIZER_HELPER.", call. = FALSE)
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
       "Initializer helper changed.")

freeze <- file.path(parent, "stage7j_seed02_blind_freeze_20260930_v1")
marker <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
lines <- readLines(marker, warn = FALSE)
assert(identical(lines[[1L]], "STAGE7J_BLIND_ENDPOINTS_FROZEN") &&
         all(c("branches=U,A", "fit_id=covrb7c_02__C1__start11") %in% lines),
       "Blind freeze marker changed.")
manifest_path <- file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv")
marker_hash <- lines[startsWith(lines, "manifest_sha256=")]
assert(length(marker_hash) == 1L &&
         identical(sha(manifest_path), sub("^manifest_sha256=", "", marker_hash)),
       "Blind manifest hash changed.")
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
assert(identical(manifest$branch, c("U", "A")) &&
         all(manifest$fit_id == "covrb7c_02__C1__start11"),
       "Branch manifest changed.")
for (i in seq_len(nrow(manifest))) for (role in c("initial", "fit", "snapshots")) {
  path <- manifest[[paste0(role, "_path")]][[i]]
  expected <- manifest[[paste0(role, "_sha256")]][[i]]
  assert(identical(sha(path), expected),
         paste("Frozen input changed", manifest$branch[[i]], role))
}

observation_path <- file.path(parent, "covariate_routeb_stage7c_20260927_v1",
                              "data", "covrb7c_02__C1", "observation_bundle.rds")
assert(identical(sha(observation_path),
                 "e9bff47094381f3005fcb7d29e56d0854bb436fafe4ebc5cd9d4fecc4bc6ae32"),
       "Bound observation bundle changed.")
observation <- readRDS(observation_path)
fit <- readRDS(manifest$fit_path[[1L]])
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

output <- file.path(parent, "stage7j_seed02_blind_trajectory_readout_20260930_v1")
assert(!file.exists(output) && !dir.exists(output), "Readout output exists.")
assert(dir.create(output, mode = "0700"), "Cannot create readout output.")

a_states <- readRDS(manifest$snapshots_path[[2L]])
a_iterations <- vapply(a_states, `[[`, integer(1L), "iteration")
expected <- as.integer(c(51:70, 75, 80, 90, 100, 120, 150, 200,
                         250, 300, 350, 400, 450, 499))
assert(identical(a_iterations, expected), "A snapshot schedule changed.")
rm(a_states)
invisible(gc())

all_rows <- vector("list", 2L)
for (b in seq_len(2L)) {
  branch <- manifest$branch[[b]]
  states <- readRDS(manifest$snapshots_path[[b]])
  iterations <- vapply(states, `[[`, integer(1L), "iteration")
  if (branch == "U") {
    assert(identical(iterations, as.integer(51:499)),
           "U full snapshot schedule changed.")
    states <- states[match(expected, iterations)]
  } else {
    assert(identical(iterations, expected), "A snapshot schedule changed.")
  }
  initial <- readRDS(manifest$initial_path[[b]])
  assert(identical(initial$iteration, 50L), "Initial sweep changed.")
  states <- c(list(initial), states)
  branch_root <- if (branch == "U")
    file.path(parent, "stage7j_seed02_u_20260930_v1", "u_full_v1") else
    file.path(parent, "stage7j_seed02_a_20260930_v1", "A_full_v1")
  driver <- utils::read.csv(file.path(branch_root, "DRIVER_TRACE.csv"),
                            stringsAsFactors = FALSE)
  driver <- driver[driver$block == "specific_loading" &
                     driver$study == 2L & driver$factor == 1L, , drop = FALSE]
  assert(!anyDuplicated(driver$iteration) &&
           all(expected %in% driver$iteration),
         paste("B2 loading driver index changed", branch))
  rows <- lapply(states, function(state) {
    energy <- state_energy(state, C2, p, N_obs)
    rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
    rho_a <- as.numeric(state$mu_q_gamma_a)
    rho_b1 <- as.numeric(state$mu_q_gamma_b[[1L]])
    driver_row <- driver[driver$iteration == state$iteration, , drop = FALSE]
    assert(nrow(driver_row) <= 1L, "Duplicate driver row.")
    data.frame(
      branch = branch, iteration = state$iteration,
      annealing = isTRUE(state$annealing), temperature = state$temperature,
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
      b2_loading_data_driver_norm = if (nrow(driver_row))
        driver_row$data_driver_norm else NA_real_,
      b2_loading_data_coherence = if (nrow(driver_row))
        driver_row$data_coherence else NA_real_,
      b2_loading_posterior_mean_norm = if (nrow(driver_row))
        driver_row$posterior_mean_norm else NA_real_,
      b2_loading_precision_trace = if (nrow(driver_row))
        driver_row$precision_trace else NA_real_,
      stringsAsFactors = FALSE)
  })
  all_rows[[b]] <- do.call(rbind, rows)
  rm(states, rows, driver)
  invisible(gc())
  cat("STAGE7J_TRAJECTORY_BRANCH_COMPLETE branch=", branch, "\n", sep = "")
}
trajectory <- do.call(rbind, all_rows)
path <- file.path(output, "STAGE7J_B2_TRAJECTORY_SELECTED.csv")
utils::write.csv(trajectory, path, row.names = FALSE, quote = TRUE, na = "")
writeLines(c("STAGE7J_BLIND_TRAJECTORY_READOUT_COMPLETE",
             "branches=U,A", "truth_not_used_in_this_readout=TRUE",
             paste0("freeze_marker_sha256=", sha(marker)),
             paste0("trajectory_sha256=", sha(path))),
           file.path(output, "TRAJECTORY_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7J_BLIND_TRAJECTORY_READOUT_COMPLETE\n")
