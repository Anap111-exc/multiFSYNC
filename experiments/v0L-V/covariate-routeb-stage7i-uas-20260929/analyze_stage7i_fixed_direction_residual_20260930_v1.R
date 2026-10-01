#!/usr/bin/env Rscript

# Read-only, truth-blind reconstruction of the study-2 leave-B2-out residual.
# The variable direction is fixed at the registered sweep-50 spectral direction.
# Usage: Rscript SCRIPT RESULTS_PARENT PACKAGE_LIB OUTPUT_DIR

options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected RESULTS_PARENT PACKAGE_LIB OUTPUT_DIR.")
parent <- normalizePath(args[[1L]], mustWork = TRUE)
package_lib <- normalizePath(args[[2L]], mustWork = TRUE)
output <- args[[3L]]
if (file.exists(output) || dir.exists(output)) stop("Refusing existing output.")
.libPaths(unique(c(package_lib, .libPaths())))
stopifnot(requireNamespace("digest", quietly = TRUE),
          requireNamespace("multiFSYNC", quietly = TRUE))
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                      check.names = FALSE)
freeze <- file.path(parent, "stage7i_blind_freeze_20260929_v1")
marker <- readLines(file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt"), warn = FALSE)
assert(identical(marker[[1L]], "STAGE7I_BLIND_ENDPOINTS_FROZEN") &&
         "truth_bundle_opened=FALSE" %in% marker, "Blind freeze marker changed.")
marker_field <- function(prefix) {
  hit <- marker[startsWith(marker, prefix)]
  assert(length(hit) == 1L, paste("Missing marker field", prefix))
  substring(hit, nchar(prefix) + 1L)
}
hashes_path <- file.path(freeze, "BLIND_FILE_HASHES.csv")
assert(identical(sha256(hashes_path), marker_field("file_hashes_sha256=")),
       "Blind hash table changed.")
hashes <- csv(hashes_path)
assert(identical(names(hashes), c("path", "sha256")), "Hash schema changed.")
bound <- function(relative) {
  assert(!grepl("(^/|^[A-Za-z]:|\\.\\.)", relative), "Unsafe relative path.")
  hit <- hashes[hashes$path == relative, , drop = FALSE]
  assert(nrow(hit) == 1L, paste("Unbound input", relative))
  path <- file.path(parent, relative)
  assert(identical(sha256(path), hit$sha256[[1L]]),
         paste("Bound input mismatch", relative))
  path
}
manifest_path <- file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv")
assert(identical(sha256(manifest_path), marker_field("manifest_sha256=")),
       "Blind manifest changed.")
manifest <- csv(manifest_path)
assert(identical(manifest$branch, c("U", "A", "S")), "Branch order changed.")
observation_path <- file.path(parent, "covariate_routeb_stage7c_20260927_v1",
                              "data", "covrb7c_01__C1", "observation_bundle.rds")
assert(identical(sha256(observation_path),
                 "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6"),
       "Frozen observation bundle changed.")
observation <- readRDS(observation_path)
fit <- readRDS(bound(manifest$fit_path[[1L]]))
grid <- multiFSYNC:::get_grid_objects(
  unlist(observation$time_obs, recursive = FALSE), fit$K,
  n_g = fit$n_g, time_g = fit$time_g, format_univ = TRUE)
n1 <- fit$n_s[[1L]]
n2 <- fit$n_s[[2L]]
C2 <- grid$C[(n1 + 1L):(n1 + n2)]
Y2 <- lapply(observation$Y[[2L]], function(y) do.call(cbind, y))
Z2 <- observation$Z[[2L]]
p <- fit$p
K_total <- fit$K_total
d <- fit$d
assert(n2 == 30L && p == 500L && all(vapply(seq_len(n2), function(i)
  identical(dim(Y2[[i]]), c(nrow(C2[[i]]), p)), logical(1L))),
  "Frozen study-2 dimensions changed.")
rm(grid, fit)
invisible(gc())

as_root <- "stage7i_as_20260929_v1"
diagnostics <- readRDS(bound(file.path(as_root,
  "BLIND_INITIALIZATION_DIAGNOSTICS.rds")))
v <- as.numeric(diagnostics$variable_direction)
assert(length(v) == p && abs(sum(v^2) - 1) < 1e-10,
       "Registered variable direction is invalid.")
initial_summary <- csv(bound(file.path(as_root,
  "BLIND_INITIALIZATION_SUMMARY.csv")))
branch_roots <- c(U = "stage7i_u_20260929_v1/u_full_v1",
                  A = "stage7i_as_20260929_v1/A_full_v1",
                  S = "stage7i_as_20260929_v1/S_full_v1")

extract <- function(state) list(
  mean = do.call(cbind, state$mu_q_nu_mu[[2L]]),
  beta = lapply(seq_len(d), function(r)
    do.call(cbind, lapply(state$mu_q_nu_beta, `[[`, r))),
  phi = state$mu_q_nu_phi[[1L]],
  zeta = state$mu_q_zeta[[2L]][[1L]],
  a = as.numeric(state$mu_q_a[, 1L]),
  psi = state$mu_q_nu_psi[[2L]][[1L]],
  xi = state$mu_q_xi[[2L]][[1L]],
  tau = as.numeric(state$mu_q_recip_sigsq_eps[2L, ]))

# Only the study-2 mean, beta and shared component are removed. Study-1 B1
# cannot enter this residual. Individual signed contributions are summed before
# taking any square or norm.
measure <- function(branch, iteration, stage, m, tau) {
  assert(identical(dim(m$mean), c(K_total, p)) &&
           identical(dim(m$psi), c(K_total, 2L)) &&
           identical(dim(m$xi), c(n2, 2L)) &&
           length(tau) == p && all(is.finite(tau)) && all(tau > 0),
         "Mixed-state dimensions or noise precision changed.")
  contrib <- matrix(0, n2, p)
  residual_sse <- projected_residual_sse <- g_sq <- 0
  for (i in seq_len(n2)) {
    Ci <- C2[[i]]
    mu_i <- Ci %*% m$mean
    beta_i <- matrix(0, nrow(Ci), p)
    for (r in seq_len(d))
      beta_i <- beta_i + Z2[i, r] * (Ci %*% m$beta[[r]])
    shared_i <- as.vector(Ci %*% m$phi %*% m$zeta[i, ])
    R_i <- Y2[[i]] - mu_i - beta_i - tcrossprod(shared_i, m$a)
    g_i <- as.vector(Ci %*% m$psi %*% m$xi[i, ])
    contrib[i, ] <- tau * as.vector(crossprod(g_i, R_i))
    residual_sse <- residual_sse + sum(R_i^2)
    projected_residual_sse <- projected_residual_sse + sum((R_i %*% v)^2)
    g_sq <- g_sq + sum(g_i^2)
  }
  h <- colSums(contrib)
  v_by_subject <- as.vector(contrib %*% v)
  h_norm <- sqrt(sum(h^2))
  data.frame(branch = branch, iteration = iteration, stage = stage,
             tau_rule = "sweep_previous",
             residual_sse = residual_sse,
             fixed_v_residual_sse = projected_residual_sse,
             b2_g_norm = sqrt(g_sq), loading_data_driver_norm = h_norm,
             fixed_v_signed_driver = sum(v_by_subject),
             fixed_v_driver_abs_cosine = abs(sum(v_by_subject)) /
               max(h_norm, .Machine$double.eps),
             fixed_v_subject_cancellation = abs(sum(v_by_subject)) /
               max(sum(abs(v_by_subject)), .Machine$double.eps),
             global_subject_cancellation = h_norm /
               max(sum(sqrt(rowSums(contrib^2))), .Machine$double.eps))
}
mix <- function(old, new, fields) {
  ans <- old
  for (field in fields) ans[[field]] <- new[[field]]
  ans
}
rows <- list()
bindings <- list()
for (branch in names(branch_roots)) {
  row <- manifest[manifest$branch == branch, , drop = FALSE]
  snapshot_path <- bound(row$snapshot_path[[1L]])
  initial_path <- bound(file.path(as_root,
                                   paste0("STATE_", branch, "_SWEEP50.rds")))
  driver_path <- bound(file.path(branch_roots[[branch]], "DRIVER_TRACE.csv"))
  trace <- csv(driver_path)
  trace <- trace[!is.na(trace$study) & !is.na(trace$factor) &
                 trace$study == 2L & trace$factor == 1L &
                 trace$block == "specific_loading" &
                 trace$iteration %in% 51:55, , drop = FALSE]
  assert(nrow(trace) == 5L, "Specific loading driver trace changed.")
  all_states <- readRDS(snapshot_path)
  assert(length(all_states) == 449L &&
           identical(vapply(all_states[1:5], `[[`, integer(1L), "iteration"),
                     as.integer(51:55)), "Snapshot order changed.")
  states <- c(list(readRDS(initial_path)), all_states[1:5])
  rm(all_states)
  invisible(gc())
  initial <- extract(states[[1L]])
  tau50 <- initial$tau
  initial_row <- measure(branch, 50L, "injected_end", initial, tau50)
  expected <- initial_summary[
    initial_summary$branch == branch, "fixed_residual_driver_norm"]
  assert(length(expected) == 1L &&
           abs(initial_row$loading_data_driver_norm - expected) <
             1e-8 + 1e-8 * abs(expected),
         paste("Sweep-50 fixed driver parity failed", branch))
  out <- list(initial_row)
  for (step in 1:5) {
    old <- extract(states[[step]])
    new <- extract(states[[step + 1L]])
    iteration <- 50L + step
    tau <- old$tau
    m0 <- old
    m1 <- mix(old, new, c("mean", "beta", "phi"))
    m2 <- mix(m1, new, "psi")
    m3 <- mix(m2, new, "zeta")
    m4 <- mix(m3, new, "xi")
    m5 <- mix(m4, new, "a")
    stages <- list(before_updates = m0, after_mu_beta_phi = m1,
                   after_psi = m2, after_zeta = m3,
                   after_xi = m4, before_specific_loading = m5)
    stage_rows <- lapply(names(stages), function(name)
      measure(branch, iteration, name, stages[[name]], tau))
    actual <- stage_rows[[length(stage_rows)]]$loading_data_driver_norm
    reference <- trace$data_driver_norm[trace$iteration == iteration]
    assert(length(reference) == 1L &&
             abs(actual - reference) < 1e-8 + 1e-8 * abs(reference),
           paste("Specific loading driver parity failed", branch, iteration,
                 "actual", actual, "reference", reference))
    out <- c(out, stage_rows)
  }
  rows[[branch]] <- do.call(rbind, out)
  bindings[[branch]] <- data.frame(branch = branch,
    initial_sha256 = sha256(initial_path),
    snapshots_sha256 = sha256(snapshot_path),
    driver_trace_sha256 = sha256(driver_path))
  rm(states)
  invisible(gc())
  cat("FIXED_DIRECTION_RESIDUAL_BRANCH_COMPLETE branch=", branch, "\n", sep = "")
}
assert(dir.create(output, mode = "0700"), "Cannot create output directory.")
output_table <- file.path(output, "B2_FIXED_DIRECTION_RESIDUAL_50_55.csv")
utils::write.csv(do.call(rbind, rows), output_table,
                 row.names = FALSE, quote = TRUE, na = "")
utils::write.csv(do.call(rbind, bindings),
                 file.path(output, "RESIDUAL_INPUT_SHA256.csv"),
                 row.names = FALSE, quote = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                      value = TRUE)[[1L]])
writeLines(c("status=COMPLETE", "branches=U,A,S", "sweeps=50:55",
             "fixed_direction=sweep50_registered_v", "truth_opened=FALSE",
             "new_fit_started=FALSE", "cavi_update_replayed=FALSE",
             paste0("script_sha256=", sha256(script_file)),
             paste0("table_sha256=", sha256(output_table))),
           file.path(output, "FIXED_DIRECTION_RESIDUAL_COMPLETE.txt"),
           useBytes = TRUE)
cat("STAGE7I_FIXED_DIRECTION_RESIDUAL_COMPLETE\n")
