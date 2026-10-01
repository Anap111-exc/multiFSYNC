#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script <- normalizePath(sub("^--file=", "", argument[[1L]]),
                        winslash = "/", mustWork = TRUE)
root <- dirname(script)
repo <- normalizePath(file.path(root, "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required.")
}

package_sources <- list.files(file.path(repo, "R"), pattern = "[.]R$",
                              full.names = TRUE)
fixed <- c(
  file.path(repo, "DESCRIPTION"),
  file.path(repo, "NAMESPACE"),
  file.path(root, "PROTOCOL_JAOUA_SPLINETS_PAIRED_CONFIRMATION_V1_20260909.md"),
  file.path(root, "DATA_MANIFEST.csv"),
  file.path(root, "FIT_MANIFEST.csv"),
  file.path(root, "MANIFEST_QC.csv"),
  file.path(root, "build_paired_manifests_20260909_v1.R"),
  file.path(root, "build_source_binding_20260909_v1.R"),
  file.path(root, "run_paired_confirmation_20260909_v1.R"),
  file.path(root, "launch_server_20260909_v1.sh"),
  file.path(root, "server_autostop_20260909_v1.sh"),
  file.path(repo, "experiments", "v0L-V",
            "g12-dimension-overspec-stage6a-20260831",
            "stage6a_common_20260831_v1.R")
)
patterns <- c(
  "^v0lv_candidate2_runtime[.]R$", "^r_route_v3_candidate2_runner[.]R$",
  "^v0lv_candidate2_controller[.]R$", "^v0f_data_calibration[.]R$"
)
snapshot <- file.path(repo, "experiments", "v0L-V", "rc1_snapshot")
bound_runtime <- unlist(lapply(patterns, function(pattern) {
  hit <- list.files(snapshot, pattern = pattern, recursive = TRUE,
                    full.names = TRUE)
  if (length(hit) != 1L) {
    stop("Expected one runtime file for ", pattern, "; found ", length(hit), ".")
  }
  hit
}), use.names = FALSE)
files <- unique(normalizePath(c(package_sources, fixed, bound_runtime),
                              winslash = "/", mustWork = TRUE))
prefix <- paste0(repo, "/")
if (!all(startsWith(files, prefix))) stop("A source file is outside the repo.")
relative <- substring(files, nchar(prefix) + 1L)
binding <- data.frame(
  relative_path = relative,
  sha256 = vapply(files, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  bytes = unname(file.info(files)$size),
  stringsAsFactors = FALSE
)
binding <- binding[order(binding$relative_path), , drop = FALSE]
write.csv(binding, file.path(root, "SOURCE_BINDING.csv"),
          row.names = FALSE, quote = TRUE, na = "")
cat("SOURCE_BINDING_PASS rows=", nrow(binding), "\n", sep = "")
