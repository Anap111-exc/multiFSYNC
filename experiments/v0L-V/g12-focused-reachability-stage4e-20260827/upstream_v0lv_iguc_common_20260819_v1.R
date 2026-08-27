options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

uc_abort <- function(...) stop(paste0(...), call. = FALSE)
uc_assert <- function(value, message) if (!isTRUE(value)) uc_abort(message)
uc_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
uc_iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}
uc_script_root <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  uc_assert(length(argument) == 1L, "Run this script with Rscript.")
  dirname(normalizePath(sub("^--file=", "", argument[[1L]]),
                        winslash = "/", mustWork = TRUE))
}
uc_parse_cli <- function(arguments) {
  result <- list()
  for (argument in arguments) {
    uc_assert(startsWith(argument, "--") && grepl("=", argument, fixed = TRUE),
              "Arguments must use --key=value syntax.")
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    uc_assert(nzchar(key) && !key %in% names(result),
              "Duplicate or empty CLI key.")
    result[[key]] <- paste(pieces[-1L], collapse = "=")
  }
  result
}
`%||%` <- function(left, right) if (is.null(left)) right else left

UC_EXPERIMENT_ID <-
  "V0LV_INIT_GEOMETRY_UNSEEN_CONFIRMATION_V1_20260819"
UC_PARENT_COMMIT <- "aee98b79a80a6535a0ad7a7c5187b29db2f66176"
UC_PARENT_TAG <- "v0.3.0-v0L-V-frozen1"
UC_DEV_BRANCH <- "dev/v0lv-init-geometry-20260818"
UC_DEV_COMMIT <- "90f4fdcd47d573abba4024a87cbf667878ccc0f7"
UC_DATA_IDS <- sprintf("iguc_%02d", 1:10)
UC_METHOD_IDS <- c(
  "current_1_over_m", "jaoua_iid_coefficients", "gram_unit_energy"
)
UC_CALIBRATION_MODES <- c(
  current_1_over_m = "none",
  jaoua_iid_coefficients = "coefficient_iid",
  gram_unit_energy = "function"
)

UC_ROOT <- uc_script_root()
UC_DEV_REPO <- uc_env(
  "V0LV_DEV_REPO", "/root/multiFSYNC-init-geometry-dev-20260818-v1"
)
UC_FROZEN_REPO <- uc_env("V0LV_FROZEN_REPO", "/root/multiFSYNC-git-rc1")
UC_PARENT_FORMAL_ROOT <- uc_env(
  "V0LV_PARENT_FORMAL_ROOT", "/root/v0lv-formal-v0.3.0-frozen1"
)
UC_DEV_LIB <- uc_env("V0LV_DEV_LIB", "/root/v0lv-init-geometry-dev-lib")
UC_BINDING_DIR <- file.path(
  UC_PARENT_FORMAL_ROOT,
  "protocol_manifests_20260813_v3_frozen1_authorized1"
)
UC_SNAPSHOT <- file.path(
  UC_FROZEN_REPO, "experiments", "v0L-V", "rc1_snapshot"
)
UC_FROZEN_LIB <- file.path(UC_SNAPSHOT, "_r_test_lib")

uc_git <- function(repo, arguments) {
  result <- system2("git", c("-C", repo, arguments),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) uc_abort(paste(result, collapse = "\n"))
  result
}
uc_find_one <- function(root, pattern) {
  hits <- list.files(root, pattern = pattern, recursive = TRUE,
                     full.names = TRUE)
  uc_assert(length(hits) == 1L,
            paste0("Expected one file matching ", pattern,
                   "; found ", length(hits), "."))
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}
uc_read_csv <- function(path) {
  uc_assert(file.exists(path), paste0("Required CSV is missing: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
uc_sha256 <- function(path) {
  uc_assert(file.exists(path), paste0("Cannot hash missing file: ", path))
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
uc_verify_environment <- function(
    load_runtime = TRUE,
    package_role = c("development", "frozen_generation")) {
  package_role <- match.arg(package_role)
  uc_assert(dir.exists(UC_DEV_REPO) && dir.exists(UC_FROZEN_REPO),
            "Development or frozen Git checkout is missing.")
  uc_assert(dir.exists(UC_PARENT_FORMAL_ROOT) && dir.exists(UC_BINDING_DIR) &&
              dir.exists(UC_SNAPSHOT),
            "Frozen parent result, binding, or snapshot is missing.")
  uc_assert(dir.exists(UC_DEV_LIB) && dir.exists(UC_FROZEN_LIB),
            "Required R library is missing.")

  binding <- uc_read_csv(file.path(UC_ROOT, "SOURCE_BINDING.csv"))
  uc_assert(nrow(binding) == 1L &&
              all(c("dev_commit", "dev_branch", "package_source_sha256") %in%
                    names(binding)), "Source binding is malformed.")
  dev_commit <- uc_git(UC_DEV_REPO, c("rev-parse", "HEAD"))
  dev_branch <- uc_git(UC_DEV_REPO, c("branch", "--show-current"))
  uc_assert(identical(dev_commit, UC_DEV_COMMIT) &&
              identical(dev_commit, binding$dev_commit[[1L]]) &&
              identical(dev_branch, UC_DEV_BRANCH) &&
              identical(dev_branch, binding$dev_branch[[1L]]) &&
              length(uc_git(UC_DEV_REPO, c("status", "--porcelain=v1"))) == 0L,
            "Development Git binding or cleanliness failed.")
  source_path <- file.path(UC_DEV_REPO, "R", "diagnostic_driver_trace.R")
  uc_assert(identical(uc_sha256(source_path),
                      binding$package_source_sha256[[1L]]),
            "Development initializer source hash mismatch.")

  uc_assert(identical(uc_git(UC_FROZEN_REPO, c("rev-parse", "HEAD")),
                      UC_PARENT_COMMIT) &&
              identical(uc_git(UC_FROZEN_REPO,
                               c("describe", "--tags", "--exact-match", "HEAD")),
                        UC_PARENT_TAG) &&
              length(uc_git(UC_FROZEN_REPO,
                            c("status", "--porcelain=v1"))) == 0L,
            "Frozen Git binding or cleanliness failed.")

  selected_library <- if (package_role == "development") UC_DEV_LIB else
    UC_FROZEN_LIB
  .libPaths(unique(c(selected_library, .libPaths())))
  uc_assert(requireNamespace("multiFSYNC", quietly = TRUE),
            paste0(package_role, " multiFSYNC package is unavailable."))
  package_path <- normalizePath(find.package("multiFSYNC"), winslash = "/")
  expected_library <- normalizePath(selected_library, winslash = "/")
  uc_assert(startsWith(package_path, paste0(expected_library, "/")),
            paste0("multiFSYNC was not loaded from ", package_role, " library."))

  if (isTRUE(load_runtime)) {
    source(uc_find_one(UC_SNAPSHOT, "^v0lv_candidate2_runtime[.]R$"),
           local = FALSE)
    source(uc_find_one(UC_SNAPSHOT, "^r_route_v3_candidate2_runner[.]R$"),
           local = FALSE)
    source(uc_find_one(UC_SNAPSHOT, "^v0lv_candidate2_controller[.]R$"),
           local = FALSE)
    if (package_role == "frozen_generation") {
      c2_verify_formal_binding(UC_BINDING_DIR, UC_SNAPSHOT)
    } else {
      status <- uc_read_csv(file.path(UC_BINDING_DIR, "BINDING_STATUS.csv"))
      index <- uc_read_csv(file.path(
        UC_BINDING_DIR, "FULL_RUNTIME_SOURCE_BINDINGS.csv"
      ))
      verification <- c2_verify_hash_index(index, UC_SNAPSHOT)
      uc_assert(nrow(status) == 1L && status$status[[1L]] == "frozen" &&
                  isTRUE(status$formal_execution_authorized[[1L]]) &&
                  all(verification$matches),
                "Frozen parent file-level binding failed.")
    }
  }
  invisible(binding)
}

uc_validate_data_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "data_id", "data_index", "data_seed",
    "new_unseen_at_registration", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_selection",
    "formal_v0lv_result"
  )
  uc_assert(is.data.frame(manifest) && nrow(manifest) == 10L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
              identical(manifest$data_id, UC_DATA_IDS) &&
              !anyDuplicated(manifest$data_seed) &&
              all(manifest$new_unseen_at_registration) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_selection) &&
              !any(manifest$formal_v0lv_result),
            "Data manifest contract failed.")
  invisible(TRUE)
}

uc_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "method_id",
    "calibration_mode", "seed_index", "fit_seed", "paired_seed_id",
    "n_cpus", "pre_score_sweeps", "planned_annealing_sweeps",
    "maximum_ordinary_T1_sweeps", "objective_eligibility_required",
    "outer_parallel_fit_limit", "actual_racing", "continuation",
    "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_within_method_selection", "formal_v0lv_result"
  )
  uc_assert(is.data.frame(manifest) && nrow(manifest) == 360L &&
              all(required %in% names(manifest)) &&
              all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
              setequal(manifest$data_id, UC_DATA_IDS) &&
              setequal(manifest$method_id, UC_METHOD_IDS) &&
              !anyDuplicated(manifest$fit_id), "Fit manifest identity failed.")
  expected_modes <- unname(UC_CALIBRATION_MODES[manifest$method_id])
  uc_assert(identical(as.character(manifest$calibration_mode), expected_modes) &&
              all(table(manifest$data_id, manifest$method_id) == 12L) &&
              all(manifest$n_cpus == 1L) &&
              all(manifest$pre_score_sweeps == 1L) &&
              all(manifest$planned_annealing_sweeps == 99L) &&
              all(manifest$maximum_ordinary_T1_sweeps == 200L) &&
              all(manifest$objective_eligibility_required) &&
              all(manifest$outer_parallel_fit_limit == 4L) &&
              !any(manifest$actual_racing) && !any(manifest$continuation) &&
              !any(manifest$truth_available_to_fit) &&
              !any(manifest$truth_available_to_stopping) &&
              !any(manifest$truth_available_to_within_method_selection) &&
              !any(manifest$formal_v0lv_result),
            "Fit route, budget, or isolation contract failed.")
  paired <- aggregate(fit_seed ~ data_id + seed_index, manifest,
                      function(value) length(unique(value)))
  uc_assert(all(paired$fit_seed == 1L), "Fit methods are not seed-paired.")
  invisible(TRUE)
}

uc_load_observation <- function(data_id) {
  path <- file.path(UC_ROOT, "data", data_id, "observation_bundle.rds")
  uc_assert(file.exists(path), paste0("Observation bundle missing: ", data_id))
  observation <- readRDS(path)
  uc_assert(identical(observation$bundle_class,
                      "v0lv_observation_only_candidate2") &&
              identical(observation$data_id, data_id),
            "Observation bundle class or identity mismatch.")
  forbidden <- candidate2_forbidden_observation_names(observation)
  uc_assert(!length(forbidden), "Observation bundle leaks truth fields.")
  observation
}

uc_terminal_summary <- function(record) {
  data.frame(
    fit_id = record$fit_id, data_id = record$data_id,
    method_id = record$method_id,
    calibration_mode = record$calibration_mode,
    seed_index = record$seed_index, fit_seed = record$fit_seed,
    paired_seed_id = record$paired_seed_id,
    terminal_status = record$terminal_status,
    actual_T1_sweeps = record$actual_T1_sweeps,
    strict_practical_converged = record$strict_practical_converged,
    objective_eligible = record$objective_eligible,
    final_elbo = record$final_elbo,
    warning_count = if (is.data.frame(record$warnings)) nrow(record$warnings) else 0L,
    elapsed_seconds = record$elapsed_seconds,
    peak_memory_bytes = record$peak_memory_bytes,
    truth_used_for_fit_or_selection = record$truth_used_for_fit_or_selection,
    formal_v0lv_result = record$formal_v0lv_result,
    stringsAsFactors = FALSE
  )
}

uc_select_winners <- function(terminals) {
  uc_assert(is.data.frame(terminals) && nrow(terminals) == 360L &&
              !any(terminals$truth_used_for_fit_or_selection),
            "Truth-free selection input is invalid.")
  groups <- split(terminals, interaction(
    terminals$data_id, terminals$method_id, drop = TRUE, lex.order = TRUE
  ))
  winners <- lapply(groups, function(group) {
    eligible <- group$objective_eligible & is.finite(group$final_elbo) &
      group$terminal_status %in% c(
        "strict_practical_converged", "max_budget_reached"
      )
    uc_assert(any(eligible), paste0(
      "No eligible endpoint for ", group$data_id[[1L]], "/",
      group$method_id[[1L]], "."
    ))
    candidate <- which(eligible)
    ordered <- candidate[order(
      -group$final_elbo[candidate], group$fit_seed[candidate],
      group$seed_index[candidate]
    )]
    winner <- group[ordered[[1L]], , drop = FALSE]
    winner$selection_rule <-
      "maximum_valid_ordinary_T1_ELBO_within_method_data"
    winner
  })
  result <- do.call(rbind, winners)
  rownames(result) <- NULL
  result[order(result$data_id, result$method_id), , drop = FALSE]
}
