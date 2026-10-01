#!/usr/bin/env Rscript

# Post-freeze truth diagnostic of saved initial and early loading directions.
# No fit, endpoint choice, or algorithm change.
# Usage: Rscript SCRIPT RESULTS_PARENT
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected RESULTS_PARENT.", call. = FALSE)
parent <- normalizePath(args[[1L]], mustWork = TRUE)
stopifnot(requireNamespace("digest", quietly = TRUE))
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
out <- file.path(parent, "stage7j_seed02_initial_loading_truth_diagnostic_20260930_v1")
assert(!file.exists(out) && !dir.exists(out), "Output already exists.")

freeze <- file.path(parent, "stage7j_seed02_blind_freeze_20260930_v1")
marker <- file.path(freeze, "BLIND_ENDPOINTS_FROZEN.txt")
assert(identical(readLines(marker, n = 1L), "STAGE7J_BLIND_ENDPOINTS_FROZEN"),
       "Seed02 blind freeze is missing.")
manifest <- utils::read.csv(file.path(freeze, "BLIND_ENDPOINT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
assert(identical(manifest$branch, c("U", "A")), "Seed02 branches changed.")
for (j in seq_len(nrow(manifest))) for (role in c("initial", "snapshots")) {
  path <- manifest[[paste0(role, "_path")]][[j]]
  assert(identical(sha(path), manifest[[paste0(role, "_sha256")]][[j]]),
         paste("Seed02 bound state changed", role, manifest$branch[[j]]))
}
s_binding <- utils::read.csv(file.path(parent, "stage7j_seed02_a_20260930_v1",
                                          "INJECTED_STATE_BINDING.csv"),
                             stringsAsFactors = FALSE)
s_row <- s_binding[s_binding$role == "S", , drop = FALSE]
assert(nrow(s_row) == 1L &&
         identical(sha(s_row$path[[1L]]), s_row$sha256[[1L]]),
       "Seed02 S construction state changed.")

truth_paths <- c(
  seed01 = file.path(parent, "covariate_routeb_stage7c_20260927_v1",
                     "data", "covrb7c_01__C1", "sealed_truth", "truth_bundle.rds"),
  seed02 = file.path(parent, "covariate_routeb_stage7c_20260927_v1",
                     "data", "covrb7c_02__C1", "sealed_truth", "truth_bundle.rds"))
expected_truth <- c(
  seed01 = "2101ae33f95597b8756437e673a93e690a9d5e7872199307fe6735ed38b46934",
  seed02 = "b6e843ee9d0fc6b346d618670ef16db513dd7eb10b9ef272bacf4d743c4aa86d")
assert(all(vapply(names(truth_paths), function(s)
  identical(sha(truth_paths[[s]]), expected_truth[[s]]), logical(1L))),
  "Sealed truth hash changed.")

# Only after the frozen inputs and sealed truth hashes pass.
targets <- lapply(names(truth_paths), function(s) {
  bundle <- readRDS(truth_paths[[s]])
  expected_id <- if (s == "seed01") "covrb7c_01__C1" else "covrb7c_02__C1"
  assert(identical(bundle$data_id, expected_id), "Truth identity changed.")
  truth <- bundle$data$true_params
  c(shared = list(as.numeric(truth$a_true[, 1L])),
    B2 = list(as.numeric(truth$b_true[[2L]][, 1L])))
})
names(targets) <- names(truth_paths)
cos_abs <- function(x, y) {
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx <= sqrt(.Machine$double.eps) || ny <= 0) return(NA_real_)
  abs(sum(x * y) / (nx * ny))
}

old_root <- file.path(parent, "stage7i_as_20260929_v1")
old_initial <- stats::setNames(lapply(c("U", "A", "S"), function(b)
  file.path(old_root, paste0("STATE_", b, "_SWEEP50.rds"))), c("U", "A", "S"))
new_initial <- c(U = manifest$initial_path[[1L]],
                 A = manifest$initial_path[[2L]],
                 S = s_row$path[[1L]])
all_rows <- list()
input_rows <- list()
index <- 0L
for (seed in c("seed01", "seed02")) {
  t_shared <- targets[[seed]][["shared"]]
  t_b2 <- targets[[seed]][["B2"]]
  assert(length(t_shared) == 500L && length(t_b2) == 500L,
         "Truth loading length changed.")
  paths <- if (seed == "seed01") old_initial else new_initial
  for (branch in names(paths)) {
    path <- paths[[branch]]
    input_rows[[length(input_rows) + 1L]] <- data.frame(
      seed = seed, branch = branch, role = "initial", path = path,
      sha256 = sha(path))
    state <- readRDS(path)
    b <- as.numeric(state$mu_q_b_specific[[2L]][, 1L])
    rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
    index <- index + 1L
    all_rows[[index]] <- data.frame(
      seed = seed, branch = branch, iteration = state$iteration,
      effective_loading_norm = sqrt(sum(b^2)),
      abs_cos_true_B2 = cos_abs(b, t_b2),
      abs_cos_true_shared = cos_abs(b, t_shared),
      true_shared_B2_abs_cos = cos_abs(t_shared, t_b2),
      b2_ppi_mean = mean(rho), b2_ppi_ge_half = sum(rho >= 0.5))
  }
}
states_path <- manifest$snapshots_path[[2L]]
input_rows[[length(input_rows) + 1L]] <- data.frame(
  seed = "seed02", branch = "A", role = "selected_snapshots",
  path = states_path, sha256 = sha(states_path))
states <- readRDS(states_path)
t_shared <- targets[["seed02"]][["shared"]]
t_b2 <- targets[["seed02"]][["B2"]]
for (state in states) {
  if (!(state$iteration %in% c(51:60, 65L, 70L, 80L, 90L, 499L))) next
  b <- as.numeric(state$mu_q_b_specific[[2L]][, 1L])
  rho <- as.numeric(state$mu_q_gamma_b[[2L]][, 1L])
  index <- index + 1L
  all_rows[[index]] <- data.frame(
    seed = "seed02", branch = "A", iteration = state$iteration,
    effective_loading_norm = sqrt(sum(b^2)),
    abs_cos_true_B2 = cos_abs(b, t_b2),
    abs_cos_true_shared = cos_abs(b, t_shared),
    true_shared_B2_abs_cos = cos_abs(t_shared, t_b2),
    b2_ppi_mean = mean(rho), b2_ppi_ge_half = sum(rho >= 0.5))
}
assert(dir.create(out, mode = "0700"), "Cannot create output.")
utils::write.csv(do.call(rbind, all_rows),
                 file.path(out, "LOADING_DIRECTION_TRUTH_DIAGNOSTIC.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
utils::write.csv(do.call(rbind, input_rows),
                 file.path(out, "INPUT_SHA256.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
writeLines(c("STAGE7J_INITIAL_LOADING_TRUTH_DIAGNOSTIC_COMPLETE",
             "fit_or_selection_started=FALSE", "diagnostic_only=TRUE",
             paste0("freeze_marker_sha256=", sha(marker)),
             paste0("table_sha256=", sha(file.path(
               out, "LOADING_DIRECTION_TRUTH_DIAGNOSTIC.csv")))),
           file.path(out, "EVALUATION_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7J_INITIAL_LOADING_TRUTH_DIAGNOSTIC_COMPLETE\n")
