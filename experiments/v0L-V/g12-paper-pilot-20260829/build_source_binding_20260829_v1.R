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

package_paths <- file.path(repo, c(
  "R/generate_data_structured.R",
  "R/diagnostic_driver_trace.R",
  "R/pre_score_interface.R",
  "R/g12_stopping_control.R",
  "R/multi_core.R",
  "DESCRIPTION", "NAMESPACE", ".gitattributes"
))
snapshot_paths <- c(
  find_one("^v0f_data_calibration[.]R$"),
  find_one("^v0lv_candidate2_runtime[.]R$"),
  find_one("^r_route_v3_candidate2_runner[.]R$"),
  find_one("^v0lv_candidate2_controller[.]R$"),
  find_one("^v0lv_candidate2_evaluation[.]R$")
)
bayes_paths <- file.path(snapshot, "Rcode", "bayesSYNC_ref", "R", c(
  "bayesSYNC_package.R", "utils.R", "OSullivan_splines.R",
  "set_hyper.R", "bayesSYNC.R"
))
pilot_names <- c(
  "README.md", "LOCAL_VALIDATION_20260830.md",
  "PROTOCOL_G12_PAPER_PIPELINE_PILOT_V1_20260829.md",
  "build_manifests_20260829_v1.R", "build_source_binding_20260829_v1.R",
  "DATA_MANIFEST.csv",
  "FIT_MANIFEST.csv", "METRIC_MANIFEST.csv", "MANIFEST_QC.csv",
  "paper_pilot_common_20260829_v1.R",
  "paper_pilot_runner_20260829_v1.R",
  "prepare_server_20260829_v1.sh", "launch_server_20260829_v1.sh",
  "status_server_20260829_v1.sh"
)
pilot_paths <- file.path(root, pilot_names)
paths <- c(package_paths, snapshot_paths, bayes_paths, pilot_paths)
if (!all(file.exists(paths))) {
  stop("Missing source-binding file(s): ",
       paste(paths[!file.exists(paths)], collapse = ";"))
}

repo_normalized <- normalizePath(repo, winslash = "/", mustWork = TRUE)
paths_normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
prefix <- paste0(repo_normalized, "/")
if (!all(startsWith(paths_normalized, prefix))) {
  stop("A source-binding path lies outside the repository root.")
}
relative <- substring(paths_normalized, nchar(prefix) + 1L)
roles <- c(
  rep("current_multifsync_runtime", length(package_paths)),
  rep("frozen_candidate2_evaluation_and_control", length(snapshot_paths)),
  rep("frozen_pooled_bayessync_reference", length(bayes_paths)),
  rep("paper_pilot_protocol_and_runner", length(pilot_paths))
)
binding <- data.frame(
  protocol_id = "G12_PAPER_PIPELINE_PILOT_V1_20260829",
  model_parent_commit = "6b10c796c0ad2d53435c9a134d9476313bf94586",
  source_role = roles,
  relative_path = relative,
  bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1L)),
  stringsAsFactors = FALSE
)
if (anyDuplicated(binding$relative_path) || any(!nzchar(binding$sha256))) {
  stop("Source binding contains duplicate paths or empty hashes.")
}
write_csv_lf(binding, file.path(root, "SOURCE_BINDING.csv"))
cat("SOURCE_BINDING_COMPLETE rows=", nrow(binding), "\n", sep = "")
