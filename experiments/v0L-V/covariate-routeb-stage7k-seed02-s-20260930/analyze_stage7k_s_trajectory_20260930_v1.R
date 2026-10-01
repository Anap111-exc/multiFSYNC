#!/usr/bin/env Rscript

# Read-only reduction of frozen Stage 7K S snapshots; no fit or truth read.
# Usage: Rscript SCRIPT RESULTS_PARENT PACKAGE_LIB
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Expected RESULTS_PARENT PACKAGE_LIB.", call. = FALSE)
parent <- normalizePath(args[[1L]], mustWork = TRUE)
lib <- normalizePath(args[[2L]], mustWork = TRUE)
.libPaths(unique(c(lib, .libPaths())))
stopifnot(requireNamespace("multiFSYNC", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                           check.names = FALSE)
check <- function(path, expected, label) {
  assert(file.exists(path) && identical(sha(path), expected),
         paste(label, "is missing or differs from its SHA256 binding."))
}
script <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
freeze <- file.path(parent, "stage7k_seed02_blind_freeze_20260930_v1")
marker_path <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
marker <- readLines(marker_path, warn = FALSE)
assert(identical(marker[[1L]], "STAGE7K_BLIND_ENDPOINTS_FROZEN") &&
         all(c("branches=U,A,S", "fit_id=covrb7c_02__C1__start11") %in% marker),
       "Stage 7K freeze marker changed.")
bound_sha <- function(prefix) {
  line <- marker[startsWith(marker, prefix)]
  assert(length(line) == 1L, paste("Missing", prefix))
  sub(paste0("^", prefix), "", line)
}
manifest_path <- file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv")
check(manifest_path, bound_sha("manifest_sha256="), "Stage 7K manifest")
manifest <- read_csv(manifest_path)
assert(identical(manifest$branch, c("U", "A", "S")),
       "Stage 7K branch order changed.")
s <- manifest[manifest$branch == "S", , drop = FALSE]
assert(nrow(s) == 1L, "S branch is not unique.")
for (role in c("initial", "fit", "snapshots"))
  check(s[[paste0(role, "_path")]][[1L]],
        s[[paste0(role, "_sha256")]][[1L]], paste("S", role))
assert(identical(s$initial_sha256[[1L]],
  "72a19fd46ed9f096966eca6c5b6f9bd673a6010e51a356ddcaa7e5156ad6d292"),
  "S state is not the originally frozen Stage 7J S state.")

u_root <- file.path(parent, "stage7j_seed02_u_20260930_v1")
u_binding <- read_csv(file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv"))
u_paths <- stats::setNames(u_binding$path, u_binding$role)
u_hashes <- stats::setNames(u_binding$sha256, u_binding$role)
observation_path <- u_paths[["observation"]]
check(observation_path, u_hashes[["observation"]], "observation bundle")
observation <- readRDS(observation_path)
assert(identical(observation$data_id, "covrb7c_02__C1") &&
         isFALSE(observation$truth_available_to_fit),
       "Observation-only bundle changed.")
a_binding <- read_csv(file.path(parent, "stage7j_seed02_a_20260930_v1",
                                     "SOURCE_BINDING.csv"))
a_paths <- stats::setNames(a_binding$path, a_binding$role)
a_hashes <- stats::setNames(a_binding$sha256, a_binding$role)
helper_path <- a_paths[["initializer"]]
check(helper_path, a_hashes[["initializer"]], "Stage 7J helper")
assert(identical(sha(helper_path),
  "03ae170103e2d70e63889b68553d9d13261233edfa1ec775ce185de39f840c99"),
  "Stage 7J helper source changed.")

fit <- readRDS(s$fit_path[[1L]])
assert(fit$i_iter == 499L && fit$t1_sweeps == 400L && fit$p == 500L,
       "S fit dimensions or endpoint changed.")
grid <- multiFSYNC:::get_grid_objects(
  unlist(observation$time_obs, recursive = FALSE), fit$K,
  n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
n1 <- fit$n_s[[1L]]
n2 <- fit$n_s[[2L]]
C2 <- grid$C[(n1 + 1L):(n1 + n2)]
N_obs <- sum(vapply(C2, nrow, integer(1L)))
assert(length(C2) == 30L, "Study 2 subject count changed.")
rm(fit, grid)
invisible(gc(verbose = FALSE))
helper <- new.env(parent = globalenv())
sys.source(helper_path, envir = helper)
state_energy <- helper$.stage7i_state_energy

states <- readRDS(s$snapshots_path[[1L]])
selected <- as.integer(c(51:70, 75, 80, 90, 100, 120, 150, 200,
                         250, 300, 350, 400, 450, 499))
assert(identical(vapply(states, `[[`, integer(1L), "iteration"), selected),
       "S selected snapshot schedule changed.")
initial <- readRDS(s$initial_path[[1L]])
assert(identical(initial$iteration, 50L) && isTRUE(initial$annealing),
       "S initial state changed.")
states <- c(list(initial), states)
driver_path <- file.path(parent, "stage7k_seed02_s_20260930_v1",
                         "S_full_v1", "DRIVER_TRACE.csv")
binding <- read_csv(file.path(parent, "stage7k_seed02_s_20260930_v1",
                                   "S_full_v1", "BRANCH_OUTPUT_BINDING.csv"))
driver_sha <- binding$sha256[binding$role == "driver_trace"]
assert(length(driver_sha) == 1L, "S driver binding missing.")
check(driver_path, driver_sha, "S driver trace")
driver <- read_csv(driver_path)
driver <- driver[driver$block == "specific_loading" &
                   driver$study == 2L & driver$factor == 1L, , drop = FALSE]
assert(!anyDuplicated(driver$iteration) &&
         identical(driver$iteration, selected), "S B2 driver index changed.")

rows <- lapply(states, function(state) {
  energy <- state_energy(state, C2, 500L, N_obs)
  rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
  rho_b1 <- as.numeric(state$mu_q_gamma_b[[1L]])
  rho_a <- as.numeric(state$mu_q_gamma_a)
  d <- driver[driver$iteration == state$iteration, , drop = FALSE]
  data.frame(branch = "S", iteration = state$iteration,
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
    b2_loading_data_driver_norm = if (nrow(d)) d$data_driver_norm else NA_real_,
    b2_loading_data_coherence = if (nrow(d)) d$data_coherence else NA_real_,
    b2_loading_posterior_mean_norm =
      if (nrow(d)) d$posterior_mean_norm else NA_real_,
    b2_loading_precision_trace = if (nrow(d)) d$precision_trace else NA_real_,
    stringsAsFactors = FALSE)
})
trajectory <- do.call(rbind, rows)
output <- file.path(parent, "stage7k_seed02_trajectory_readout_20260930_v1")
assert(!file.exists(output) && !dir.exists(output),
       "Refusing existing trajectory readout output.")
assert(dir.create(output, mode = "0700"), "Cannot create trajectory output.")
path <- file.path(output, "STAGE7K_S_B2_TRAJECTORY_SELECTED.csv")
utils::write.csv(trajectory, path, row.names = FALSE, quote = TRUE, na = "")
writeLines(c("STAGE7K_S_TRAJECTORY_READOUT_COMPLETE",
             "branch=S", "truth_read_by_readout=FALSE",
             "executed_after_truth_evaluation=TRUE",
             paste0("script_sha256=", sha(script)),
             paste0("freeze_marker_sha256=", sha(marker_path)),
             paste0("trajectory_sha256=", sha(path))),
           file.path(output, "READOUT_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7K_S_TRAJECTORY_READOUT_COMPLETE rows=", nrow(trajectory), "\n", sep = "")
