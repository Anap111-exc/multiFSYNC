# Read-only final verification for candidate2 deliverables. It creates only a
# versioned verification report and never generates formal data or runs a fit.

candidate2_verification_script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run final verification with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/",
                mustWork = TRUE)
}

script_path <- candidate2_verification_script_path()
base_dir <- dirname(script_path)
project_root <- normalizePath(file.path(base_dir, "..", ".."),
                              winslash = "/", mustWork = TRUE)
local_library <- normalizePath(file.path(project_root, "_r_test_lib"),
                               winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))
source(file.path(base_dir, "v0lv_candidate2_runtime.R"), local = FALSE)

manifest_dir <- file.path(base_dir, "protocol_manifests_20260812_v3_candidate2")
test_dir <- file.path(base_dir, "candidate2_test_results_20260812_v7")
smoke_dir <- file.path(base_dir, "candidate2_e2e_smoke_20260812_v8")
package_test_dir <- file.path(base_dir, "candidate2_full_package_tests_20260812_v1")
output_dir <- file.path(base_dir, "candidate2_final_verification_20260812_v1")
if (dir.exists(output_dir)) stop("Final verification output already exists.")
dir.create(output_dir, recursive = TRUE)

parse_paths <- c(
  list.files(file.path(project_root, "Rcode", "multiFSYNC", "R"),
             pattern = "[.]R$", full.names = TRUE),
  list.files(base_dir, pattern = "candidate2.*[.]R$|[.]R$", full.names = TRUE)
)
parse_paths <- sort(unique(parse_paths))
parse_results <- do.call(rbind, lapply(parse_paths, function(path) {
  condition <- tryCatch({ parse(file = path, encoding = "UTF-8"); NULL },
                        error = identity)
  data.frame(
    relative_path = c2_relative_path(path, project_root),
    parse_passed = is.null(condition),
    error_class = if (is.null(condition)) "" else
      paste(class(condition), collapse = ";"),
    error_message = if (is.null(condition)) "" else conditionMessage(condition),
    stringsAsFactors = FALSE
  )
}))
c2_atomic_write_csv(parse_results,
                    file.path(output_dir, "CANDIDATE2_PARSE_RESULTS.csv"))

read_csv <- function(name) utils::read.csv(
  file.path(manifest_dir, name), stringsAsFactors = FALSE,
  check.names = FALSE
)
data_manifest <- read_csv("V0LV_DATA_SEEDS_CANDIDATE2.csv")
fit_manifest <- read_csv("V0LV_FIT_MANIFEST_CANDIDATE2.csv")
audit_manifest <- read_csv("V0LV_AUDIT_MANIFEST_CANDIDATE2.csv")
binding_status <- read_csv("BINDING_STATUS.csv")
integrity <- read_csv("CANDIDATE_INTEGRITY_CHECKS.csv")
source_index <- read_csv("FULL_RUNTIME_SOURCE_BINDINGS.csv")
source_verification <- c2_verify_hash_index(source_index, project_root)
parent_index <- read_csv("PARENT_HISTORY_HASHES.csv")
parent_verification <- c2_verify_hash_index(parent_index, project_root)

tests <- utils::read.csv(file.path(test_dir, "TEST_RESULTS.csv"),
                         stringsAsFactors = FALSE)
smoke <- utils::read.csv(file.path(smoke_dir, "SMOKE_STATUS.csv"),
                         stringsAsFactors = FALSE)
package_log <- readLines(file.path(
  package_test_dir, "FULL_TESTTHAT_INSTALLED_ASCII.log"
), warn = FALSE)
install_log <- readLines(file.path(
  package_test_dir, "R_CMD_INSTALL_ASCII.log"
), warn = FALSE)

environment_expected <- read_csv("ENVIRONMENT_BINDING.csv")
environment_actual <- c2_current_environment_binding(
  environment_expected$package[environment_expected$binding_type == "package"]
)
environment_key <- function(value) paste(
  value$binding_type, value$package, value$key, sep = "\r"
)
expected_key <- environment_key(environment_expected)
actual_key <- environment_key(environment_actual)
environment_match <- !anyDuplicated(expected_key) &&
  !anyDuplicated(actual_key) && setequal(expected_key, actual_key)
if (environment_match) {
  environment_actual <- environment_actual[match(expected_key, actual_key), ]
  environment_match <- all(
    as.character(environment_expected$value) ==
      as.character(environment_actual$value)
  )
}

if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
  stop("Installed multiFSYNC is unavailable for namespace verification.")
}
namespace_expected <- read_csv("MULTIFSYNC_NAMESPACE_FUNCTION_BINDINGS.csv")
namespace <- asNamespace("multiFSYNC")
namespace_names <- sort(ls(namespace, all.names = TRUE)[vapply(
  ls(namespace, all.names = TRUE),
  function(name) is.function(get(name, envir = namespace, inherits = FALSE)),
  logical(1)
)])
namespace_inventory_match <- identical(
  sort(as.character(namespace_expected$function_name)), namespace_names
)
namespace_hash_match <- FALSE
if (namespace_inventory_match) {
  namespace_actual_hash <- vapply(namespace_names, function(name) {
    c2_function_signature_sha256(
      get(name, envir = namespace, inherits = FALSE)
    )
  }, character(1))
  namespace_expected_hash <- setNames(
    namespace_expected$sha256, namespace_expected$function_name
  )[namespace_names]
  namespace_hash_match <- all(namespace_actual_hash == namespace_expected_hash)
}

formal_rejection <- tryCatch({
  c2_verify_formal_binding(manifest_dir, project_root)
  list(rejected = FALSE, message = "")
}, error = function(condition) list(
  rejected = TRUE, message = conditionMessage(condition)
))
formal_result_dir <- file.path(
  base_dir, "v0lv_candidate2_formal_results_NOT_STARTED"
)

method_counts <- table(factor(
  fit_manifest$path_id,
  levels = c("B0", "A", "R", "pooled_bayesSYNC")
))
checks <- data.frame(
  check_id = c(
    "all_candidate2_and_package_R_parse", "manifest_integrity_checks_pass",
    "formal_data_count", "formal_fit_count_and_methods", "audit_record_count",
    "all_formal_start_flags_false", "binding_is_candidate_not_frozen",
    "formal_preflight_rejects_candidate", "no_formal_result_directory",
    "source_hashes_match", "parent_history_hashes_match",
    "environment_binding_matches", "namespace_inventory_matches",
    "namespace_function_hashes_match", "candidate2_regression_tests_pass",
    "candidate2_e2e_smoke_pass", "full_package_install_pass",
    "full_package_testthat_pass"
  ),
  passed = c(
    all(parse_results$parse_passed), all(integrity$passed),
    nrow(data_manifest) == 10L,
    nrow(fit_manifest) == 210L &&
      identical(as.integer(method_counts), c(40L, 40L, 120L, 10L)),
    nrow(audit_manifest) == 72L,
    !any(data_manifest$formal_generation_started) &&
      !any(fit_manifest$formal_fit_started) &&
      !any(audit_manifest$formal_continuation_started),
    binding_status$status[[1L]] == "candidate" &&
      !binding_status$formal_execution_authorized[[1L]] &&
      !binding_status$frozen_bundle_created[[1L]],
    formal_rejection$rejected && grepl(
      "authorized frozen binding", formal_rejection$message, fixed = TRUE
    ),
    !dir.exists(formal_result_dir), all(source_verification$matches),
    all(parent_verification$matches), environment_match,
    namespace_inventory_match, namespace_hash_match,
    nrow(tests) == 17L && all(tests$status == "PASS") &&
      all(!tests$formal_experiment),
    nrow(smoke) == 12L && all(smoke$status == "PASS") &&
      all(!smoke$formal_experiment),
    any(grepl("DONE [(]multiFSYNC[)]", install_log, fixed = FALSE)),
    any(grepl("STAGE_A_ROUND2_FULL_PACKAGE_TESTS_PASS", package_log,
              fixed = TRUE))
  ),
  detail = c(
    paste0(sum(parse_results$parse_passed), "/", nrow(parse_results)),
    paste0(sum(integrity$passed), "/", nrow(integrity)),
    nrow(data_manifest),
    paste(names(method_counts), as.integer(method_counts), collapse = ";"),
    nrow(audit_manifest), "generation=0;fit=0;continuation=0",
    paste0("status=", binding_status$status[[1L]],
           ";authorized=", binding_status$formal_execution_authorized[[1L]]),
    formal_rejection$message, formal_result_dir,
    paste0(sum(source_verification$matches), "/", nrow(source_verification)),
    paste0(sum(parent_verification$matches), "/", nrow(parent_verification)),
    "R/platform/BLAS/LAPACK/key package versions",
    paste0(length(namespace_names), " functions"),
    paste0(length(namespace_names), " function signatures"),
    paste0(sum(tests$status == "PASS"), "/", nrow(tests)),
    paste0(sum(smoke$status == "PASS"), "/", nrow(smoke)),
    "R CMD INSTALL --install-tests exit 0 and DONE marker",
    "installed-package test_dir success marker"
  ), stringsAsFactors = FALSE
)
if (any(!checks$passed)) {
  c2_atomic_write_csv(checks,
                      file.path(output_dir, "CANDIDATE2_FINAL_CHECKS.csv"))
  stop("Final candidate2 checks failed: ", paste(
    checks$check_id[!checks$passed], collapse = ", "
  ))
}
checks_write <- c2_atomic_write_csv(
  checks, file.path(output_dir, "CANDIDATE2_FINAL_CHECKS.csv")
)
c2_atomic_write_csv(source_verification, file.path(
  output_dir, "CANDIDATE2_SOURCE_HASH_VERIFICATION.csv"
))
c2_atomic_write_csv(parent_verification, file.path(
  output_dir, "CANDIDATE2_PARENT_HASH_VERIFICATION.csv"
))

deliverables <- c(
  file.path(base_dir, c(
    "PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V3_CANDIDATE2_20260812.md",
    "PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812.md",
    "v0lv_candidate2_runtime.R", "r_route_v3_candidate2_runner.R",
    "v0lv_candidate2_controller.R", "v0lv_candidate2_evaluation.R",
    "run_v0lv_candidate2_formal.R", "run_v0lv_candidate2_e2e_smoke.R",
    "test_v0lv_candidate2.R", "build_protocol_manifests_v3_candidate2.R",
    "run_candidate2_final_verification.R"
  )),
  list.files(manifest_dir, full.names = TRUE),
  list.files(test_dir, full.names = TRUE),
  list.files(smoke_dir, recursive = TRUE, full.names = TRUE),
  list.files(package_test_dir, full.names = TRUE)
)
deliverables <- deliverables[file.exists(deliverables) & !dir.exists(deliverables)]
deliverable_hashes <- c2_hash_table(
  deliverables, project_root, role = "candidate2_deliverable_or_evidence"
)
hash_write <- c2_atomic_write_csv(
  deliverable_hashes, file.path(output_dir, "CANDIDATE2_DELIVERABLE_HASHES.csv")
)
c2_atomic_write_lines(c(
  "V0LV_CANDIDATE2_FINAL_VERIFICATION_COMPLETE",
  paste0("checks=", nrow(checks)),
  paste0("passed=", sum(checks$passed)), "failed=0",
  "status=candidate_waiting_for_user_review", "formal_experiment=FALSE",
  "formal_data_generated=FALSE", "formal_fits_started=FALSE",
  "formal_continuations_started=FALSE", "frozen_bundle_created=FALSE",
  paste0("checks_sha256=", checks_write$sha256),
  paste0("deliverable_hashes_sha256=", hash_write$sha256),
  capture.output(sessionInfo())
), file.path(output_dir, "CANDIDATE2_FINAL_VERIFICATION_COMPLETE.txt"))
print(checks, row.names = FALSE)
