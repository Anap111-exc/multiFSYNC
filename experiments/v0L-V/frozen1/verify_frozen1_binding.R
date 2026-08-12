# Read-only verifier for the v0L-V frozen1 binding.

script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run this verifier with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]),
                winslash = "/", mustWork = TRUE)
}

assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- sub("/+$", "", normalizePath(root, winslash = "/", mustWork = TRUE))
  prefix <- paste0(root, "/")
  assert(startsWith(path, prefix), paste("Path is outside repository:", path))
  substring(path, nchar(prefix) + 1L)
}

arguments <- commandArgs(trailingOnly = TRUE)
unknown <- setdiff(arguments, "--require-runtime")
assert(!length(unknown), paste("Unknown argument(s):", paste(unknown, collapse = ", ")))
require_runtime <- "--require-runtime" %in% arguments

frozen_root <- dirname(script_path())
repo_root <- normalizePath(file.path(frozen_root, "..", "..", ".."),
                           winslash = "/", mustWork = TRUE)
snapshot_root <- file.path(repo_root, "experiments", "v0L-V", "rc1_snapshot")
candidate_base <- file.path(
  snapshot_root, "对比实验", "v0L阶段A收口与v2冻结准备_20260812"
)
candidate_binding <- file.path(
  candidate_base, "protocol_manifests_20260812_v3_candidate2"
)
binding_dir <- file.path(
  frozen_root, "protocol_manifests_20260813_v3_frozen1"
)
runtime_path <- file.path(candidate_base, "v0lv_candidate2_runtime.R")
source(runtime_path, local = FALSE)

status <- read_csv(file.path(binding_dir, "BINDING_STATUS.csv"))
assert(nrow(status) == 1L, "Frozen binding status must have one row.")
assert(identical(status$status[[1L]], "frozen"), "Binding is not frozen.")
assert(!isTRUE(status$formal_execution_authorized[[1L]]),
       "Frozen1 must not authorize formal execution.")
assert(isTRUE(status$frozen_bundle_created[[1L]]),
       "Frozen bundle marker is false.")
assert(!any(status[c("formal_data_generated", "formal_fits_started",
                     "formal_continuations_started")]),
       "A formal-start flag is true.")
assert(!isTRUE(status$requires_separate_user_freeze_authorization[[1L]]) &&
         isTRUE(status$requires_separate_user_start_authorization[[1L]]),
       "Freeze/start authorization flags are inconsistent.")

data_manifest <- read_csv(file.path(binding_dir, "V0LV_DATA_SEEDS_CANDIDATE2.csv"))
fit_manifest <- read_csv(file.path(binding_dir, "V0LV_FIT_MANIFEST_CANDIDATE2.csv"))
audit_manifest <- read_csv(file.path(binding_dir, "V0LV_AUDIT_MANIFEST_CANDIDATE2.csv"))
assert(nrow(data_manifest) == 10L && nrow(fit_manifest) == 210L &&
         nrow(audit_manifest) == 72L, "Manifest counts are not 10/210/72.")
method_counts <- table(factor(fit_manifest$path_id,
                              levels = c("B0", "A", "R", "pooled_bayesSYNC")))
assert(identical(as.integer(method_counts), c(40L, 40L, 120L, 10L)),
       "Method counts are not 40/40/120/10.")
assert(!any(data_manifest$formal_generation_started) &&
         !any(fit_manifest$formal_fit_started) &&
         !any(audit_manifest$formal_continuation_started),
       "A manifest formal-start flag is true.")
assert(!any(fit_manifest$truth_available_to_fit) &&
         !any(fit_manifest$truth_available_to_stopping) &&
         !any(fit_manifest$truth_available_to_selection) &&
         !any(audit_manifest$truth_available_to_continuation) &&
         !any(audit_manifest$truth_available_to_counterfactual_selection),
       "Truth isolation is not intact.")

source_index <- read_csv(file.path(binding_dir, "FULL_RUNTIME_SOURCE_BINDINGS.csv"))
assert(nrow(source_index) == 50L && !anyDuplicated(source_index$relative_path),
       "Runtime source inventory is not 50 unique paths.")
source_paths <- file.path(snapshot_root, source_index$relative_path)
assert(all(file.exists(source_paths) & !dir.exists(source_paths)),
       "A runtime source file is missing.")
source_sizes <- unname(file.info(source_paths)$size)
source_hashes <- vapply(source_paths, c2_sha256, character(1))
assert(all(source_sizes == source_index$bytes), "Runtime source byte mismatch.")
assert(all(source_hashes == tolower(source_index$sha256)),
       "Runtime source SHA-256 mismatch.")

scientific_files <- c(
  "FULL_RUNTIME_SOURCE_BINDINGS.csv", "OBJECTIVE_VALIDITY_CONTRACTS.csv",
  "PARENT_HISTORY_HASHES.csv", "POOLED_EVALUATION_CONTRACT.csv",
  "R_OFFLINE_RACING_CONTRACT.csv", "SEED_COLLISION_AUDIT.csv",
  "V0LV_ANALYSIS_PLAN_CANDIDATE2.csv",
  "V0LV_AUDIT_MANIFEST_CANDIDATE2.csv", "V0LV_CONFIG_CANDIDATE2.csv",
  "V0LV_DATA_SEEDS_CANDIDATE2.csv", "V0LV_FIT_MANIFEST_CANDIDATE2.csv",
  "V0LV_SEED_REGISTRY_CANDIDATE2.csv",
  "V0LV_SUCCESS_CRITERIA_CANDIDATE2.csv"
)
for (name in scientific_files) {
  assert(identical(c2_sha256(file.path(binding_dir, name)),
                   c2_sha256(file.path(candidate_binding, name))),
         paste("Frozen scientific file differs from candidate2:", name))
}
candidate_status <- read_csv(file.path(candidate_binding, "BINDING_STATUS.csv"))
assert(identical(candidate_status$status[[1L]], "candidate") &&
         !isTRUE(candidate_status$formal_execution_authorized[[1L]]) &&
         !isTRUE(candidate_status$frozen_bundle_created[[1L]]),
       "Parent candidate status was altered.")

acceptance_dir <- file.path(repo_root, "experiments", "v0L-V",
                            "acceptance", "linux")
smoke <- read_csv(file.path(acceptance_dir, "SMOKE_STATUS.csv"))
assert(nrow(smoke) == 12L && all(smoke$status == "PASS") &&
         !any(smoke$formal_experiment), "Linux smoke is not 12/12 non-formal PASS.")
smoke_complete <- readLines(file.path(acceptance_dir, "SMOKE_COMPLETE.txt"),
                            warn = FALSE)
required_smoke_lines <- c(
  "status=PASS", "formal_experiment=FALSE", "formal_data_generated=FALSE",
  "formal_fits_started=FALSE", "formal_continuations_started=FALSE"
)
assert(all(required_smoke_lines %in% smoke_complete),
       "Linux smoke completion marker is inconsistent.")

artifact_index_path <- file.path(binding_dir, "FROZEN1_ARTIFACT_HASHES.csv")
artifact_index <- read_csv(artifact_index_path)
assert(nrow(artifact_index) > 0L && !anyDuplicated(artifact_index$relative_path),
       "Frozen artifact index is empty or duplicated.")
indexed_paths <- file.path(repo_root, artifact_index$relative_path)
assert(all(file.exists(indexed_paths) & !dir.exists(indexed_paths)),
       "A frozen indexed artifact is missing.")
assert(all(unname(file.info(indexed_paths)$size) == artifact_index$bytes),
       "Frozen artifact byte mismatch.")
assert(all(vapply(indexed_paths, c2_sha256, character(1)) ==
             tolower(artifact_index$sha256)), "Frozen artifact hash mismatch.")
expected_indexed <- list.files(binding_dir, full.names = TRUE)
expected_indexed <- expected_indexed[!basename(expected_indexed) %in% c(
  "FROZEN1_ARTIFACT_HASHES.csv", "FROZEN1_COMPLETE.txt"
)]
assert(setequal(artifact_index$relative_path,
                vapply(expected_indexed, relative_path, character(1),
                       root = repo_root)), "Frozen artifact inventory mismatch.")
complete <- readLines(file.path(binding_dir, "FROZEN1_COMPLETE.txt"), warn = FALSE)
assert(all(c("status=frozen", "formal_execution_authorized=FALSE",
             "formal_data_generated=FALSE", "formal_fits_started=FALSE",
             "formal_continuations_started=FALSE") %in% complete),
       "Frozen completion marker is inconsistent.")
assert(paste0("artifact_index_sha256=", c2_sha256(artifact_index_path)) %in% complete,
       "Frozen completion marker artifact-index hash mismatch.")

runtime_result <- "SKIPPED"
if (require_runtime) {
  expected_environment <- read_csv(file.path(binding_dir, "ENVIRONMENT_BINDING.csv"))
  actual_environment <- c2_current_environment_binding(unique(
    expected_environment$package[expected_environment$binding_type == "package"]
  ))
  environment_key <- function(value) {
    paste(value$binding_type, value$package, value$key, sep = "\r")
  }
  expected_key <- environment_key(expected_environment)
  actual_key <- environment_key(actual_environment)
  assert(!anyDuplicated(expected_key) && !anyDuplicated(actual_key) &&
           setequal(expected_key, actual_key), "Runtime environment keys differ.")
  actual_environment <- actual_environment[match(expected_key, actual_key), ]
  assert(all(as.character(expected_environment$value) ==
               as.character(actual_environment$value)),
         "Runtime environment values differ from frozen binding.")

  assert(requireNamespace("multiFSYNC", quietly = TRUE),
         "The frozen multiFSYNC package is not installed.")
  namespace_expected <- read_csv(file.path(
    binding_dir, "MULTIFSYNC_NAMESPACE_FUNCTION_BINDINGS.csv"
  ))
  namespace <- asNamespace("multiFSYNC")
  namespace_names <- sort(ls(namespace, all.names = TRUE)[vapply(
    ls(namespace, all.names = TRUE),
    function(name) is.function(get(name, envir = namespace, inherits = FALSE)),
    logical(1)
  )])
  assert(identical(sort(namespace_expected$function_name), namespace_names),
         "Installed namespace function inventory differs.")
  namespace_hashes <- vapply(namespace_names, function(name) {
    c2_function_signature_sha256(get(name, envir = namespace, inherits = FALSE))
  }, character(1))
  expected_hashes <- setNames(namespace_expected$sha256,
                              namespace_expected$function_name)[namespace_names]
  assert(all(namespace_hashes == expected_hashes),
         "Installed namespace function hashes differ.")

  runtime_expected <- read_csv(file.path(binding_dir, "RUNTIME_FUNCTION_BINDINGS.csv"))
  runtime_actual <- c(
    "multiFSYNC::bayesSYNC_multi" = c2_function_signature_sha256(
      getExportedValue("multiFSYNC", "bayesSYNC_multi")
    ),
    "multiFSYNC::bayesSYNC_multi_pre_score" = c2_function_signature_sha256(
      getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score")
    )
  )
  assert(setequal(runtime_expected$function_name, names(runtime_actual)) &&
           all(runtime_actual[runtime_expected$function_name] == runtime_expected$sha256),
         "Installed public runtime function hashes differ.")
  runtime_result <- "PASS"
}

cat(
  "frozen1_static_checks=PASS\n",
  "runtime_source_sha256=50/50\n",
  "runtime_source_bytes=50/50\n",
  "manifest_rows=10/210/72\n",
  "linux_smoke=12/12_PASS\n",
  "runtime_checks=", runtime_result, "\n",
  "status=frozen\n",
  "formal_execution_authorized=FALSE\n",
  "formal_work_started=FALSE\n",
  sep = ""
)
