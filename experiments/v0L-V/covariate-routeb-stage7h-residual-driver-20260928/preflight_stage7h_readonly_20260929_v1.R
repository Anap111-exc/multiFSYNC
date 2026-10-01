#!/usr/bin/env Rscript

# Stage 7H read-only preflight. It opens no truth and calls no fit routine.
options(warn = 2, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript preflight_stage7h_readonly_20260929_v1.R ",
       "STAGE7C_ROOT STAGE7G_ROOT OUTPUT_ROOT", call. = FALSE)
}
stage7c <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g <- normalizePath(args[[2L]], mustWork = TRUE)
out <- args[[3L]]
if (dir.exists(out)) stop("Refusing existing output directory: ", out)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Cannot create output directory: ", out)

sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(line) != 1L) stop("sha256sum failed: ", path)
  strsplit(line, "[[:space:]]+")[[1L]][[1L]]
}
max_abs <- function(a, b) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b))
    return(Inf)
  max(abs(as.numeric(a) - as.numeric(b)))
}
stop_if_bad <- function(ok, detail) if (!isTRUE(ok)) stop(detail, call. = FALSE)

.libPaths(c("/root/multiFSYNC-stage7g-lib-20260928-v1", .libPaths()))
suppressPackageStartupMessages(library(multiFSYNC))
stop_if_bad(as.character(utils::packageVersion("multiFSYNC")) ==
              "0.3.0.9000" &&
              identical(normalizePath(system.file(package = "multiFSYNC")),
                        "/root/multiFSYNC-stage7g-lib-20260928-v1/multiFSYNC"),
            "Frozen Stage 7G package is not loaded.")

manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
stop_if_bad(nrow(manifest) == 3L && nrow(binding) == 3L &&
              identical(as.character(manifest$fit_id),
                        as.character(binding$fit_id)),
            "Manifest or input binding changed.")

rows <- vector("list", 3L)
for (k in seq_len(3L)) {
  fit_id <- manifest$fit_id[[k]]
  data_id <- manifest$data_id[[k]]
  obs_path <- file.path(stage7c, "data", data_id,
                        "observation_bundle.rds")
  fit_path <- file.path(stage7g, "fits", fit_id, "fit.rds")
  snapshot_path <- file.path(stage7g, "fits", fit_id,
                             "RAW_STATE_SNAPSHOTS.rds")
  stop_if_bad(file.exists(obs_path) && file.exists(fit_path) &&
                file.exists(snapshot_path),
              paste0("Missing bound input: ", fit_id))
  obs_sha <- sha256(obs_path)
  fit_sha <- sha256(fit_path)
  snapshot_sha <- sha256(snapshot_path)
  stop_if_bad(identical(obs_sha, manifest$observation_sha256[[k]]) &&
                identical(fit_sha, binding$fit_sha256[[k]]) &&
                identical(snapshot_sha, binding$snapshot_sha256[[k]]),
              paste0("Input SHA256 mismatch: ", fit_id))

  obs <- readRDS(obs_path)
  fit <- readRDS(fit_path)
  stop_if_bad(identical(obs$data_id, data_id) &&
                isFALSE(obs$truth_available_to_fit) &&
                identical(as.integer(fit$n_s), c(30L, 30L)) &&
                identical(as.integer(fit$p), 500L) &&
                identical(as.integer(fit$K), 9L) &&
                identical(as.integer(fit$d), 2L),
              paste0("Observation/fit contract mismatch: ", fit_id))
  time_flat <- unlist(obs$time_obs, recursive = FALSE)
  basis <- multiFSYNC:::get_grid_objects(
    time_flat, fit$K, n_g = fit$n_g, time_g = fit$time_g,
    format_univ = TRUE)
  max_grid_error <- max_abs(basis$C_g, fit$C_g)
  max_cp_error <- max_cpy_error <- 0
  index <- 0L
  for (s in 1:2) for (i in seq_len(fit$n_s[[s]])) {
    index <- index + 1L
    C_i <- basis$C[[index]]
    Y_i <- do.call(cbind, obs$Y[[s]][[i]])
    max_cp_error <- max(max_cp_error,
                        max_abs(crossprod(C_i), fit$list_cp_C[[s]][[i]]))
    max_cpy_error <- max(max_cpy_error,
                         max_abs(crossprod(C_i, Y_i),
                                 fit$list_cp_C_Y[[s]][[i]]))
  }
  stop_if_bad(is.finite(max_grid_error) && is.finite(max_cp_error) &&
                is.finite(max_cpy_error) &&
                max(c(max_grid_error, max_cp_error, max_cpy_error)) <= 1e-7,
              paste0("Observation basis/sufficient-statistic mismatch: ",
                     fit_id, " errors=", paste(signif(c(max_grid_error,
                     max_cp_error, max_cpy_error), 7L), collapse = ",")))
  rows[[k]] <- data.frame(
    fit_id = fit_id, data_id = data_id, observation_sha256 = obs_sha,
    fit_sha256 = fit_sha, snapshot_sha256 = snapshot_sha,
    max_grid_error = max_grid_error,
    max_cp_error = max_cp_error,
    max_cpy_error = max_cpy_error,
    n_subjects_study2 = fit$n_s[[2L]], p = fit$p,
    stringsAsFactors = FALSE)
  rm(obs, fit, basis)
  invisible(gc())
}
utils::write.csv(do.call(rbind, rows),
                 file.path(out, "PREFLIGHT_BINDING_AND_DESIGN.csv"),
                 row.names = FALSE)
writeLines("STAGE7H_PREFLIGHT_PASS",
           file.path(out, "PREFLIGHT_PASS.txt"), useBytes = TRUE)
cat("STAGE7H_PREFLIGHT_PASS\n")
