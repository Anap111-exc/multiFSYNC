#!/usr/bin/env Rscript

# Freeze seed02 U/A outputs before any truth evaluation. No truth is opened.
# Usage: Rscript SCRIPT U_ROOT A_ROOT OUTPUT_ROOT
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected U_ROOT A_ROOT OUTPUT_ROOT.")
u_root <- normalizePath(args[[1L]], mustWork = TRUE)
a_root <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (file.exists(output) || dir.exists(output)) stop("Refusing existing freeze output.")
script_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_option) != 1L) stop("Cannot identify freeze source.")
this_script <- normalizePath(sub("^--file=", "", script_option))

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
sha256 <- function(path) {
  assert(file.exists(path) && !dir.exists(path), paste("Missing file", path))
  x <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  assert(length(x) == 1L && grepl("^[0-9a-f]{64}[[:space:]]", x[[1L]]),
         paste("sha256sum failed", path))
  substr(x[[1L]], 1L, 64L)
}
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                           check.names = FALSE)
read_first <- function(path) readLines(path, warn = FALSE)[[1L]]
assert(identical(read_first(file.path(u_root, "u_parity_v1", "PARITY_PASS.txt")),
                 "STAGE7J_U_PARITY_PASS"), "U full parity gate did not pass.")
assert(identical(read_first(file.path(a_root, "INJECTION_COMPLETE.txt")),
                 "STAGE7J_A_INJECTION_COMPLETE"), "A injection is incomplete.")
assert(identical(read_first(file.path(a_root, "A_FIT_COMPLETE.txt")),
                 "STAGE7J_A_FIT_COMPLETE"), "A fit is incomplete.")
assert(identical(read_first(file.path(a_root, "A_full_v1", "FIT_COMPLETE.txt")),
                 "STAGE7J_A_FIT_COMPLETE"), "A branch contract is incomplete.")
u_binding <- read_csv(file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv"))
assert(identical(names(u_binding), c("role", "path", "sha256")) &&
         !anyDuplicated(u_binding$role), "U binding schema changed.")
u_path <- stats::setNames(u_binding$path, u_binding$role)
assert(all(vapply(seq_len(nrow(u_binding)), function(i)
  identical(sha256(u_binding$path[[i]]), u_binding$sha256[[i]]),
  logical(1L))), "U source binding changed.")
assert(identical(read_first(file.path(u_root, "PREPARE_COMPLETE.txt")),
                 "STAGE7J_U_PREPARE_COMPLETE") &&
         "fit_id=covrb7c_02__C1__start11" %in%
           readLines(file.path(u_root, "PREPARE_COMPLETE.txt"), warn = FALSE),
       "U prepare marker does not select seed02/C1/start11.")
a_binding <- read_csv(file.path(a_root, "SOURCE_BINDING.csv"))
assert(identical(names(a_binding), c("role", "path", "sha256")) &&
         !anyDuplicated(a_binding$role) &&
         all(vapply(seq_len(nrow(a_binding)), function(i)
           identical(sha256(a_binding$path[[i]]), a_binding$sha256[[i]]),
           logical(1L))), "A source binding changed.")
a_path <- stats::setNames(a_binding$path, a_binding$role)
injected <- read_csv(file.path(a_root, "INJECTED_STATE_BINDING.csv"))
assert(identical(names(injected), c("role", "path", "sha256")) &&
         all(vapply(seq_len(nrow(injected)), function(i)
           identical(sha256(injected$path[[i]]), injected$sha256[[i]]),
           logical(1L))), "Injected-state binding changed.")
branch_binding <- read_csv(file.path(a_root, "A_full_v1",
                                          "BRANCH_OUTPUT_BINDING.csv"))
assert(identical(names(branch_binding), c("role", "path", "sha256")) &&
         all(vapply(seq_len(nrow(branch_binding)), function(i)
           identical(sha256(branch_binding$path[[i]]),
                     branch_binding$sha256[[i]]), logical(1L))),
       "A branch output binding changed.")

library_root <- dirname(dirname(u_path[["stage7i_description"]]))
.libPaths(c(library_root, .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
assert(identical(normalizePath(system.file(package = "multiFSYNC")),
                 normalizePath(dirname(u_path[["stage7i_description"]]))),
       "Loaded package differs from bound Stage 7I build.")
observation <- readRDS(u_path[["observation"]])
assert(identical(observation$data_id, "covrb7c_02__C1") &&
         isFALSE(observation$truth_available_to_fit),
       "Observation-only seed02 input changed.")
reference_fit <- readRDS(u_path[["stage7g_fit"]])
grid <- multiFSYNC:::get_grid_objects(
  unlist(observation$time_obs, recursive = FALSE), reference_fit$K,
  n_g = reference_fit$n_g, time_g = reference_fit$time_g,
  format_univ = TRUE)
n1 <- reference_fit$n_s[[1L]]
n2 <- reference_fit$n_s[[2L]]
C2 <- grid$C[(n1 + 1L):(n1 + n2)]
Nobs <- sum(vapply(C2, nrow, integer(1L)))
p <- reference_fit$p
helper <- new.env(parent = globalenv())
sys.source(a_path[["initializer"]], envir = helper)
assert(identical(sha256(a_path[["initializer"]]),
                 "03ae170103e2d70e63889b68553d9d13261233edfa1ec775ce185de39f840c99"),
       "Initializer helper differs from Stage 7I.")
rm(reference_fit, grid)
invisible(gc(verbose = FALSE))

branch_info <- list(
  U = list(root = file.path(u_root, "u_full_v1"),
           initial = injected$path[injected$role == "U"],
           snapshots = "RAW_STATE_SNAPSHOTS_51_499.rds"),
  A = list(root = file.path(a_root, "A_full_v1"),
           initial = injected$path[injected$role == "A"],
           snapshots = "RAW_STATE_SNAPSHOTS_SELECTED.rds"))
summary_rows <- list()
manifest_rows <- list()
for (branch in names(branch_info)) {
  info <- branch_info[[branch]]
  fit_path <- file.path(info$root, "fit.rds")
  snapshot_path <- file.path(info$root, info$snapshots)
  fit <- readRDS(fit_path)
  states <- readRDS(snapshot_path)
  iterations <- vapply(states, `[[`, integer(1L), "iteration")
  assert(length(states) > 0L && tail(iterations, 1L) == 499L &&
           fit$i_iter == 499L && fit$t1_sweeps == 400L &&
           length(fit$ELBO) == 400L && all(is.finite(fit$ELBO)),
         paste(branch, "endpoint contract changed."))
  if (branch == "U")
    assert(identical(iterations, as.integer(51:499)),
           "U full snapshot schedule changed.")
  state <- states[[length(states)]]
  energy <- helper$.stage7i_state_energy(state, C2, p, Nobs)
  assert(all(is.finite(energy)), paste(branch, "nonfinite endpoint energy."))
  rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
  summary_rows[[branch]] <- data.frame(
    branch = branch, fit_id = "covrb7c_02__C1__start11",
    final_ordinary_elbo = tail(fit$ELBO, 1L),
    b2_ppi_ge_half = sum(rho >= 0.5), b2_ppi_mean = mean(rho),
    b2_loading_norm = unname(energy[["effective_loading_norm"]]),
    b2_mean_time_norm = unname(energy[["mean_time_product_norm"]]),
    b2_mean_product_rms = unname(energy[["mean_product_rms"]]),
    b2_full_posterior_product_second_rms =
      unname(energy[["full_posterior_product_second_rms"]]),
    saved_snapshots = length(states), resumed_sweeps = 449L,
    stringsAsFactors = FALSE)
  manifest_rows[[branch]] <- data.frame(
    branch = branch, fit_id = "covrb7c_02__C1__start11",
    initial_path = normalizePath(info$initial),
    initial_sha256 = sha256(info$initial),
    fit_path = normalizePath(fit_path), fit_sha256 = sha256(fit_path),
    snapshots_path = normalizePath(snapshot_path),
    snapshots_sha256 = sha256(snapshot_path),
    stringsAsFactors = FALSE)
  rm(fit, states, state)
  invisible(gc(verbose = FALSE))
}

paths <- c(
  freeze_script = this_script,
  u_binding = file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv"),
  u_parity_pass = file.path(u_root, "u_parity_v1", "PARITY_PASS.txt"),
  u_state_parity = file.path(u_root, "u_parity_v1", "ALL_STATE_PARITY.csv"),
  u_elbo_parity = file.path(u_root, "u_parity_v1", "ORDINARY_ELBO_PARITY.csv"),
  a_source_binding = file.path(a_root, "SOURCE_BINDING.csv"),
  a_injected_binding = file.path(a_root, "INJECTED_STATE_BINDING.csv"),
  a_injection_marker = file.path(a_root, "INJECTION_COMPLETE.txt"),
  a_output_binding = file.path(a_root, "A_full_v1",
                               "BRANCH_OUTPUT_BINDING.csv"),
  a_fit_complete = file.path(a_root, "A_full_v1", "FIT_COMPLETE.txt"),
  u_fit = manifest_rows$U$fit_path,
  u_snapshots = manifest_rows$U$snapshots_path,
  a_fit = manifest_rows$A$fit_path,
  a_snapshots = manifest_rows$A$snapshots_path)
hash_table <- data.frame(role = names(paths),
                         path = unname(vapply(paths, normalizePath,
                                              character(1L), mustWork = TRUE)),
                         sha256 = unname(vapply(paths, sha256,
                                                 character(1L))),
                         stringsAsFactors = FALSE)
assert(dir.create(output, recursive = TRUE, mode = "0700"),
       "Cannot create blind freeze output.")
write_csv_new <- function(value, name) {
  path <- file.path(output, name)
  assert(!file.exists(path), paste("Refusing overwrite", path))
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
}
write_csv_new(do.call(rbind, summary_rows), "BLIND_ENDPOINT_SUMMARY.csv")
write_csv_new(do.call(rbind, manifest_rows), "BLIND_ENDPOINT_MANIFEST.csv")
write_csv_new(hash_table, "BLIND_FILE_HASHES.csv")
writeLines(c("STAGE7J_BLIND_ENDPOINTS_FROZEN", "branches=U,A",
             "fit_id=covrb7c_02__C1__start11",
             "truth_bundle_opened=FALSE", "formal_candidate_pool=FALSE",
             paste0("manifest_sha256=", sha256(file.path(
               output, "BLIND_ENDPOINT_MANIFEST.csv"))),
             paste0("file_hashes_sha256=", sha256(file.path(
               output, "BLIND_FILE_HASHES.csv"))),
             paste0("summary_sha256=", sha256(file.path(
               output, "BLIND_ENDPOINT_SUMMARY.csv")))),
           file.path(output, "BLIND_ENDPOINTS_FROZEN.txt"), useBytes = TRUE)
cat("STAGE7J_BLIND_ENDPOINTS_FROZEN\n")
