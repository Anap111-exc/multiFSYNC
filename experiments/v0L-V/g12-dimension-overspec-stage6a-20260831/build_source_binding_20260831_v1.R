#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
root <- dirname(script_file)
if (!dir.exists(root)) {
  root <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
}
repo <- file.path(root, "..", "..", "..")
snapshot <- file.path(repo, "experiments", "v0L-V", "rc1_snapshot")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")

write_csv_lf <- function(value, path) {
  lines <- capture.output(utils::write.csv(
    value, row.names = FALSE, quote = TRUE, na = ""
  ))
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(lines, connection, sep = "\n", useBytes = TRUE)
  invisible(path)
}

find_one <- function(pattern) {
  hits <- list.files(snapshot, pattern = pattern, recursive = TRUE,
                     full.names = TRUE)
  if (length(hits) != 1L) {
    stop("Expected one snapshot file matching ", pattern, "; found ",
         length(hits), ".")
  }
  hits[[1L]]
}

package_paths <- c(
  list.files(file.path(repo, "R"), pattern = "[.]R$", full.names = TRUE),
  file.path(repo, c("DESCRIPTION", "NAMESPACE", ".gitattributes"))
)
snapshot_paths <- c(
  find_one("^v0f_data_calibration[.]R$"),
  find_one("^v0lv_candidate2_runtime[.]R$"),
  find_one("^r_route_v3_candidate2_runner[.]R$"),
  find_one("^v0lv_candidate2_controller[.]R$"),
  find_one("^v0lv_candidate2_evaluation[.]R$"),
  find_one("^v0g_evaluation[.]R$"),
  find_one("^v0h_evaluation[.]R$"),
  find_one("^v0i_truth_diagnostic[.]R$"),
  find_one("^v0j_evaluation[.]R$")
)
experiment_names <- c(
  "README.md",
  "PROTOCOL_G12_DIMENSION_OVERSPEC_STAGE6A_V1_20260831.md",
  "build_manifests_20260831_v1.R", "build_source_binding_20260831_v1.R",
  "DATA_MANIFEST.csv", "FIT_CONFIGS.csv", "FIT_MANIFEST.csv",
  "METRIC_MANIFEST.csv", "MANIFEST_QC.csv",
  "stage6a_common_20260831_v1.R", "stage6a_runner_20260831_v1.R",
  "test_stage6a_contracts_20260831_v1.R",
  "prepare_server_20260831_v1.sh", "launch_server_20260831_v1.sh",
  "status_server_20260831_v1.sh"
)
experiment_paths <- file.path(
  repo, "experiments", "v0L-V", basename(root), experiment_names
)
paths <- c(package_paths, snapshot_paths, experiment_paths)
if (!all(file.exists(paths))) {
  stop("Missing source-binding file(s): ",
       paste(paths[!file.exists(paths)], collapse = ";"))
}

repo_lexical <- gsub("\\\\", "/", repo)
paths_lexical <- gsub("\\\\", "/", paths)
prefix <- paste0(repo_lexical, "/")
if (!all(startsWith(paths_lexical, prefix))) {
  stop("A source-binding path lies outside the repository root.")
}
relative <- substring(paths_lexical, nchar(prefix) + 1L)
git_index_source <- function(relative_path, fallback_path) {
  temporary <- tempfile("stage6a_git_index_")
  error_file <- tempfile("stage6a_git_index_error_")
  status <- suppressWarnings(system2(
    "git", c("-C", shQuote(repo), "show",
             shQuote(paste0(":", relative_path))),
    stdout = temporary, stderr = error_file
  ))
  unlink(error_file)
  if (identical(status, 0L) && file.exists(temporary)) {
    return(temporary)
  }
  unlink(temporary)
  fallback_path
}
bound_sources <- Map(git_index_source, relative, paths)
on.exit(unlink(unlist(bound_sources)[startsWith(
  basename(unlist(bound_sources)), "stage6a_git_index_"
)]), add = TRUE)
roles <- c(
  rep("current_multifsync_runtime", length(package_paths)),
  rep("frozen_candidate2_evaluation_and_control", length(snapshot_paths)),
  rep("stage6a_protocol_and_runner", length(experiment_paths))
)
binding <- data.frame(
  protocol_id = "G12_DIMENSION_OVERSPEC_STAGE6A_V1_20260831",
  model_parent_commit = "bede48a77993de85cec66604b98209cddf54e61b",
  source_role = roles,
  relative_path = relative,
  bytes = vapply(bound_sources, function(path) {
    as.numeric(file.info(path)$size)
  }, numeric(1L)),
  sha256 = vapply(bound_sources, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  stringsAsFactors = FALSE
)
if (anyDuplicated(binding$relative_path) || any(!nzchar(binding$sha256))) {
  stop("Source binding contains duplicate paths or empty hashes.")
}
write_csv_lf(binding, file.path(root, "SOURCE_BINDING.csv"))
cat("SOURCE_BINDING_COMPLETE rows=", nrow(binding), "\n", sep = "")
