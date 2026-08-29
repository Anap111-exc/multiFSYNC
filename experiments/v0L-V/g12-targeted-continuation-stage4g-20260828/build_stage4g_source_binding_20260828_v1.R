#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
))
stage4e_root <- Sys.getenv(
  "STAGE4G_STAGE4E_ROOT",
  unset = "/root/v0lv-g12-focused-reachability-stage4e-20260827-v1"
)
development_repo <- Sys.getenv(
  "STAGE4G_DEVELOPMENT_REPO",
  unset = "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5"
)
output <- file.path(root, "SOURCE_BINDING.csv")
if (file.exists(output)) {
  stop("Refusing to overwrite SOURCE_BINDING.csv.", call. = FALSE)
}
assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
sha256 <- function(path) {
  assert(file.exists(path), paste0("Missing source: ", path))
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ))
  }
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  assert(is.null(status) || status == 0L, "sha256sum failed.")
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
parent_path <- file.path(stage4e_root, "SOURCE_BINDING.csv")
assert(file.exists(parent_path) && dir.exists(development_repo),
       "Stage-4E binding or development source archive is missing.")
parent <- utils::read.csv(
  parent_path, stringsAsFactors = FALSE, check.names = FALSE
)
expected_commit <- "72d9a53f0ef5e9d3e5f49d9cf837958207b3c980"
assert(nrow(parent) == 1L &&
         identical(parent$dev_commit[[1L]], expected_commit) &&
         identical(parent$package_code_commit[[1L]], expected_commit),
       "Stage-4E parent source binding commit changed.")
paths <- c(
  scale_trace_sha256 = "R/diagnostic_scale_trace.R",
  multi_core_sha256 = "R/multi_core.R",
  practical_stopping_sha256 = "R/convergence_practical.R",
  g12_stopping_control_sha256 = "R/g12_stopping_control.R",
  description_sha256 = "DESCRIPTION",
  namespace_sha256 = "NAMESPACE"
)
for (field in names(paths)) {
  assert(identical(
    sha256(file.path(development_repo, paths[[field]])),
    parent[[field]][[1L]]
  ), paste0("Stage-4E parent source hash changed: ", paths[[field]]))
}
binding <- data.frame(
  experiment_id = "G12_TARGETED_CONTINUATION_STAGE4G_V1_20260828",
  parent_experiment_id = parent$experiment_id[[1L]],
  parent_source_binding_sha256 = sha256(parent_path),
  dev_commit_recorded = expected_commit,
  source_archive_has_git_metadata = dir.exists(file.path(
    development_repo, ".git"
  )),
  package_version = parent$package_version[[1L]],
  continuation_state_sha256 = sha256(file.path(
    development_repo, "R", "continuation_state.R"
  )),
  scale_trace_sha256 = parent$scale_trace_sha256[[1L]],
  multi_core_sha256 = parent$multi_core_sha256[[1L]],
  practical_stopping_sha256 = parent$practical_stopping_sha256[[1L]],
  g12_stopping_control_sha256 =
    parent$g12_stopping_control_sha256[[1L]],
  description_sha256 = parent$description_sha256[[1L]],
  namespace_sha256 = parent$namespace_sha256[[1L]],
  model_prior_cavi_elbo_changed = FALSE,
  formal_v0lv_result = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(
  binding, output, row.names = FALSE, quote = TRUE, na = ""
)
cat("STAGE4G_SOURCE_BINDING_REGISTERED git_metadata=",
    binding$source_archive_has_git_metadata[[1L]], "\n", sep = "")
