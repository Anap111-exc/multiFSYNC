#!/usr/bin/env Rscript

# Post-freeze scientific evaluation of the three Stage 7I diagnostic branches.
# Usage: Rscript SCRIPT RESULTS_PARENT SUPPORT_REPO PACKAGE_LIB STAGE7C_EVALUATOR
# The sealed truth is first deserialized only after the blind freeze and every
# source / endpoint binding have passed. No fitting or endpoint selection occurs.

options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Expected RESULTS_PARENT SUPPORT_REPO PACKAGE_LIB STAGE7C_EVALUATOR.", call. = FALSE)
results_parent <- normalizePath(args[[1L]], mustWork = TRUE)
support_repo <- normalizePath(args[[2L]], mustWork = TRUE)
package_lib <- normalizePath(args[[3L]], mustWork = TRUE)
stage7c_evaluator <- normalizePath(args[[4L]], mustWork = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
script_file <- normalizePath(script_file, mustWork = TRUE)
stopifnot(requireNamespace("digest", quietly = TRUE))
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
check_hash <- function(path, expected, label) {
  assert(file.exists(path), paste0("Missing ", label, ": ", path))
  actual <- sha256(path)
  assert(identical(tolower(actual), tolower(expected)),
         paste0("SHA256 mismatch for ", label, ": ", path))
  actual
}
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE,
                                            check.names = FALSE)
new_csv <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}
new_lines <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  writeLines(x, path, useBytes = TRUE)
}
new_rds <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  saveRDS(x, path, version = 3L)
}

freeze_root <- file.path(results_parent, "stage7i_blind_freeze_20260929_v1")
freeze_marker <- file.path(freeze_root, "BLIND_ENDPOINTS_FROZEN.txt")
manifest_path <- file.path(freeze_root, "BLIND_ENDPOINT_MANIFEST.csv")
file_hashes_path <- file.path(freeze_root, "BLIND_FILE_HASHES.csv")
u_parity_path <- file.path(results_parent, "stage7i_u_20260929_v1",
                           "u_parity_v1", "PARITY_PASS.txt")
assert(file.exists(freeze_marker), "Blind freeze marker is missing.")
marker <- readLines(freeze_marker, warn = FALSE)
assert(identical(marker[[1L]], "STAGE7I_BLIND_ENDPOINTS_FROZEN") &&
         all(c("branch_order=U,A,S", "manifest_rows=3",
               "truth_used_for_initialization=FALSE",
               "truth_used_for_endpoint_selection=FALSE",
               "formal_candidate_pool=FALSE",
               "diagnostic_case_selected_after_prior_truth_review=TRUE",
               "truth_bundle_opened=FALSE") %in% marker),
       "Blind freeze contract is incomplete.")
marker_value <- function(prefix) {
  hit <- marker[startsWith(marker, prefix)]
  assert(length(hit) == 1L, paste0("Missing freeze field ", prefix))
  substring(hit, nchar(prefix) + 1L)
}
check_hash(manifest_path, marker_value("manifest_sha256="), "blind manifest")
check_hash(file_hashes_path, marker_value("file_hashes_sha256="), "blind file hashes")
check_hash(u_parity_path, marker_value("u_parity_pass_sha256="), "U parity")
manifest <- read_csv(manifest_path)
assert(identical(as.character(manifest$branch), c("U", "A", "S")) &&
         nrow(manifest) == 3L && all(manifest$status == "COMPLETE") &&
         all(manifest$output_eligible) &&
         all(manifest$data_id == "covrb7c_01__C1") &&
         all(manifest$source_fit_id == "covrb7c_01__C1__start06") &&
         all(is.finite(manifest$final_ordinary_elbo)),
       "Three-branch manifest contract failed.")
file_hashes <- read_csv(file_hashes_path)
assert(identical(names(file_hashes), c("path", "sha256")) &&
         nrow(file_hashes) == 53L && !anyDuplicated(file_hashes$path),
       "Blind file-hash table changed.")
for (i in seq_len(nrow(file_hashes))) {
  rel <- file_hashes$path[[i]]
  assert(!grepl("(^/|^[A-Za-z]:|\\.\\.)", rel), "Unsafe relative path in blind binding.")
  check_hash(file.path(results_parent, rel), file_hashes$sha256[[i]], rel)
}
for (i in seq_len(nrow(manifest))) {
  branch <- manifest$branch[[i]]
  for (role in c("fit", "snapshot", "blind_metrics")) {
    rel <- manifest[[paste0(role, "_path")]][[i]]
    check_hash(file.path(results_parent, rel),
               manifest[[paste0(role, "_sha256")]][[i]],
               paste(branch, role))
  }
}

expected_stage7c <- "1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8"
support_v0lv <- file.path(support_repo, "experiments", "v0L-V")
stage6a <- file.path(support_v0lv, "g12-dimension-overspec-stage6a-20260831",
                     "stage6a_common_20260831_v1.R")
stage6h <- file.path(support_v0lv, "g12-reachability-causal-audit-stage6h-20260913",
                     "stage6h_evaluation_helpers_20260913_v1.R")
check_hash(stage7c_evaluator, expected_stage7c, "Stage 7C evaluator")
check_hash(stage6a, "bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6", "Stage 6A helper")
check_hash(stage6h, "5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221", "Stage 6H helper")
PP_SNAPSHOT <- file.path(support_v0lv, "rc1_snapshot")
assert(dir.exists(PP_SNAPSHOT), "Frozen rc1 snapshot is missing.")
snapshot_hashes <- c(
  v0lv_candidate2_runtime.R = "6f0eb7d9c1c01ec3d83074eee43b167cec96c53e7d681ed4ec852f4f791f76bb",
  r_route_v3_candidate2_runner.R = "e7d5808bfdd6af756c840eb84b02d87fc3883823742bd93675d1d4d5451003b6",
  v0lv_candidate2_controller.R = "ae88153e32bf15455e03e98ade5f64a91a06a47c2947adc38f3da2472a22c5ab",
  v0lv_candidate2_evaluation.R = "5ef50d51572e946e158a6b67ee28601c4b8d8847f9830fa11d881e682df24d1e",
  v0g_evaluation.R = "6b765191657f88db46a027a95acd3917c602b35cb6cee8947d917652f9c6b271",
  v0h_evaluation.R = "45639c2d9f14e4df0187d519a62837ef3b3be090ed315123cd97c725f8576dd4",
  v0i_truth_diagnostic.R = "127615551ef1e4dbe14b8c486e3e3ad27ba45eea77f764c7ed7c2fba8318b194",
  v0j_evaluation.R = "3c91af1e07b0ba1aac4539a2c1e28d8abcfe2f5ee239a2a2daed7fcb2714a833")
source_paths <- vapply(names(snapshot_hashes), function(name) {
  hits <- list.files(PP_SNAPSHOT, pattern = paste0("^", name, "$"),
                     recursive = TRUE, full.names = TRUE)
  assert(length(hits) == 1L, paste0("Expected one frozen source: ", name))
  check_hash(hits[[1L]], snapshot_hashes[[name]], name)
  hits[[1L]]
}, character(1L))
stage7c_root <- file.path(results_parent, "covariate_routeb_stage7c_20260927_v1")
baseline_root <- file.path(stage7c_root, "scientific_evaluation_20260928_v1")
baseline_fit_root <- file.path(baseline_root, "covrb7c_01__C1__start06")
baseline_science <- file.path(baseline_root, "ALL_START_SCIENCE.csv")
baseline_summary <- file.path(baseline_fit_root, "summary.csv")
baseline_role <- file.path(baseline_fit_root, "process_contribution_by_role.csv")
check_hash(baseline_science, "fcb9a4db7a8ae47062f8c73435e737b69448695b245e4050e9a8b2ef74554081", "Stage 7C science")
check_hash(baseline_summary, "d9aa10a44e6423865b0b035f7c3805df77d460cd3ab1d5122d72f112891da71f", "Stage 7C summary")
check_hash(baseline_role, "5daabb2b29316ee3b5d486b653606b9ce5893b098e7ea25c2e47d8a1c6f51bb3", "Stage 7C role table")
truth_path <- file.path(stage7c_root, "data", "covrb7c_01__C1",
                        "sealed_truth", "truth_bundle.rds")
check_hash(truth_path, "2101ae33f95597b8756437e673a93e690a9d5e7872199307fe6735ed38b46934", "sealed truth")

evaluation_root <- file.path(results_parent, "stage7i_scientific_evaluation_20260929_v4")
assert(!file.exists(evaluation_root) && !dir.exists(evaluation_root),
       "Evaluation output already exists.")
assert(dir.create(evaluation_root, mode = "0700"), "Cannot create evaluation output.")
new_lines(c("status=POST_FREEZE_TRUTH_UNSEAL_AUTHORIZED_BY_USER",
            "authorization_date=2026-09-29",
            "truth_used_for_initialization=FALSE",
            "truth_used_for_endpoint_selection=FALSE",
            "diagnostic_case_selected_after_prior_truth_review=TRUE",
            "formal_candidate_pool=FALSE",
            paste0("freeze_marker_sha256=", sha256(freeze_marker)),
            paste0("evaluator_sha256=", sha256(script_file))),
          file.path(evaluation_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"))
inputs <- c(evaluator = script_file, freeze_marker = freeze_marker,
            blind_manifest = manifest_path, blind_file_hashes = file_hashes_path,
            u_parity = u_parity_path, stage7c_evaluator = stage7c_evaluator,
            stage6a = stage6a, stage6h = stage6h,
            stats::setNames(source_paths, paste0("rc1:", names(source_paths))),
            baseline_science = baseline_science,
            baseline_summary = baseline_summary, baseline_role = baseline_role,
            sealed_truth = truth_path)
new_csv(data.frame(role = names(inputs), path = unname(inputs),
                   sha256 = vapply(inputs, sha256, character(1L)),
                   stringsAsFactors = FALSE),
        file.path(evaluation_root, "EVALUATION_INPUT_SHA256.csv"))

# First sealed-truth deserialization. Every blind and source check is above.
sealed <- readRDS(truth_path)
assert(identical(sealed$bundle_class, "covariate_routeb_pilot_sealed_truth") &&
         identical(sealed$experiment_id, "COVARIATE_ROUTEB_STAGE7C_V1_20260927") &&
         identical(sealed$data_id, "covrb7c_01__C1") &&
         identical(sealed$data$true_params$truth_basis_method,
                   "iid_common_splinet_stiefel") &&
         identical(sealed$data$true_params$mean_truth_method,
                   "independent_splinet") &&
         identical(sealed$data$true_params$beta_truth_method,
                   "independent_cubic_bspline_d2_matched_mean_centered_energy") &&
         is.list(sealed$data$Z),
       "Sealed truth identity changed.")

Sys.setenv(STAGE6A_R_LIB = package_lib, STAGE6A_USE_LOAD_ALL = "false")
source(stage6a, local = FALSE)
source(stage6h, local = FALSE)
PP_SNAPSHOT <- file.path(support_v0lv, "rc1_snapshot")
pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
functions <- c("trapezoid_weights", "interp", "mean_design", "true_mean_dense",
               "estimated_mean_dense", "mean_recovery",
               "install_dense_residual_override", "new_error_accumulator",
               "accumulate_error", "finish_rmse", "finish_nrmse",
               "new_projection_accumulator", "accumulate_projection",
               "finish_projection", "covariate_and_signal_metrics",
               "add_recovery_flags")
definitions <- parse(stage7c_evaluator, keep.source = FALSE)
found <- character()
for (expression in definitions) {
  if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) ||
      !is.symbol(expression[[2L]]) || !is.call(expression[[3L]]) ||
      !identical(expression[[3L]][[1L]], as.name("function"))) next
  name <- as.character(expression[[2L]])
  if (name %in% functions) {
    eval(expression, envir = .GlobalEnv)
    found <- c(found, name)
  }
}
assert(setequal(found, functions) && !anyDuplicated(found),
       "Could not extract the exact Stage 7C scientific functions.")
install_dense_residual_override(evaluator_environment)

evaluate_one <- function(branch) {
  row <- manifest[manifest$branch == branch, , drop = FALSE]
  branch_root <- file.path(evaluation_root, branch)
  assert(dir.create(branch_root, mode = "0700"), "Cannot create branch output.")
  warnings <- character()
  phase <- "load_fit"
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers({
    fit <- readRDS(file.path(results_parent, row$fit_path[[1L]]))
    assert(length(fit$ELBO) == 400L &&
             abs(tail(fit$ELBO, 1L) - row$final_ordinary_elbo[[1L]]) <= 1e-8,
           paste0(branch, " fit differs from the blind endpoint."))
    phase <- "frozen_factor_evaluator"
    historical <- s6h_evaluate_fit(fit, sealed, evaluator_environment)
    historical$evaluation_contract$stage7i_private_resume <- TRUE
    historical$evaluation_contract$source_sweep <- 50L
    phase <- "mean_evaluator"
    mean_result <- mean_recovery(fit, sealed$data)
    phase <- "covariate_and_full_signal_evaluator"
    covariate <- covariate_and_signal_metrics(fit, sealed,
                                              evaluator_environment, "C1")
    list(historical = historical, mean = mean_result,
         covariate = covariate, fit_seed = fit$fit_seed)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  new_lines(warnings, file.path(branch_root, "EVALUATION_WARNINGS.txt"))
  if (inherits(result, "error")) {
    new_lines(c("status=ERROR", paste0("branch=", branch),
                paste0("phase=", phase), paste0("message=", conditionMessage(result))),
              file.path(branch_root, "EVALUATION_ERROR.txt"))
    return(list(ok = FALSE, error = conditionMessage(result),
                phase = phase, elapsed_seconds = elapsed))
  }
  new_rds(result[c("historical", "mean", "covariate")],
          file.path(branch_root, "evaluation.rds"))
  for (name in names(result$historical)) {
    if (is.data.frame(result$historical[[name]])) {
      new_csv(result$historical[[name]],
              file.path(branch_root, paste0(name, ".csv")))
    }
  }
  new_csv(result$covariate$beta_errors,
          file.path(branch_root, "BETA_COEFFICIENT_ERRORS.csv"))
  new_csv(result$mean$by_study,
          file.path(branch_root, "mean_recovery_by_study.csv"))
  new_csv(result$mean$summary,
          file.path(branch_root, "mean_recovery_summary.csv"))
  legacy <- result$historical$summary
  legacy <- legacy[, !names(legacy) %in% c(
    "signal_nrmse_all_candidates", "signal_nrmse_ppi_selected",
    "signal_nrmse_mean_only"), drop = FALSE]
  names(legacy)[grepl("^dense_signal_", names(legacy))] <- paste0(
    "legacy_residual_",
    names(legacy)[grepl("^dense_signal_", names(legacy))])
  activity <- result$historical$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  base <- data.frame(
    branch = branch, fit_id = row$source_fit_id[[1L]],
    data_id = row$data_id[[1L]], logical_data_id = "covrb7c_01",
    beta_arm = "C1", start_index = 6L,
    fit_seed = result$fit_seed,
    final_ordinary_elbo = row$final_ordinary_elbo[[1L]],
    objective_eligible = row$output_eligible[[1L]],
    strict_practical_converged = read_csv(file.path(
      results_parent, row$blind_metrics_path[[1L]]))$practical_converged[[1L]],
    fit_elapsed_seconds = NA_real_,
    truth_free_winner = identical(branch, "U"),
    truth_used_for_fit_stopping_or_selection = FALSE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE)
  compact <- cbind(base, result$covariate$summary, legacy,
                   result$historical$process_contribution_summary,
                   result$historical$loading_scale_summary,
                   result$mean$summary, activity)
  new_lines(c("status=COMPLETE", paste0("branch=", branch),
              "stage7i_private_resume=TRUE", "source_sweep=50",
              "truth_used_for_fit_stopping_or_selection=FALSE",
              "new_fit_or_continuation_started=FALSE"),
            file.path(branch_root, "EVALUATION_COMPLETE.txt"))
  list(ok = TRUE, compact = compact,
       role = result$historical$process_contribution_by_role,
       elapsed_seconds = elapsed)
}

compare_tables <- function(actual, reference, label, excluded = character()) {
  assert(nrow(actual) == nrow(reference), paste0(label, " row count differs."))
  common <- intersect(setdiff(names(actual), excluded), names(reference))
  assert(length(common) > 0L, paste0(label, " has no common fields."))
  rows <- lapply(common, function(name) {
    x <- actual[[name]]
    y <- reference[[name]]
    if (is.numeric(x) && is.numeric(y)) {
      both_na <- is.na(x) & is.na(y)
      finite <- is.finite(x) & is.finite(y)
      difference <- abs(x - y)
      limit <- 1e-8 + 1e-8 * abs(y)
      ok <- all(both_na | (finite & difference <= limit))
      max_abs <- if (any(finite)) max(difference[finite]) else NA_real_
    } else {
      x <- as.character(x)
      y <- as.character(y)
      x[is.na(x) | x == ""] <- "<missing>"
      y[is.na(y) | y == ""] <- "<missing>"
      ok <- identical(x, y)
      max_abs <- NA_real_
    }
    data.frame(table = label, field = name, pass = ok,
               max_abs = max_abs, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

u <- evaluate_one("U")
assert(u$ok, paste0("U scientific evaluation failed: ", u$error))
u_science <- add_recovery_flags(u$compact)
baseline <- read_csv(baseline_science)
baseline <- baseline[baseline$fit_id == "covrb7c_01__C1__start06", , drop = FALSE]
assert(nrow(baseline) == 1L, "Stage 7C U baseline is not unique.")
parity <- rbind(
  compare_tables(read_csv(file.path(evaluation_root, "U", "summary.csv")),
                 read_csv(baseline_summary), "factor_summary"),
  compare_tables(u$role, read_csv(baseline_role), "role_contribution"),
  compare_tables(u_science, baseline, "science_row",
                 excluded = c("fit_elapsed_seconds")))
new_csv(parity, file.path(evaluation_root, "U_BASELINE_SCIENCE_PARITY.csv"))
if (!all(parity$pass)) {
  new_lines(c("status=U_SCIENCE_PARITY_FAILED",
              paste0("failed_fields=", paste(parity$field[!parity$pass],
                                              collapse = ","))),
            file.path(evaluation_root, "EVALUATION_ERROR.txt"))
  stop("U science differs from the frozen Stage 7C baseline; A/S not interpreted.",
       call. = FALSE)
}

branches <- c("A", "S")
other <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(branches, evaluate_one, mc.cores = 2L,
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
} else lapply(branches, evaluate_one)
names(other) <- branches
results <- c(list(U = u), other)
status <- data.frame(branch = names(results),
                     ok = vapply(results, `[[`, logical(1L), "ok"),
                     elapsed_seconds = vapply(results, `[[`, numeric(1L),
                                              "elapsed_seconds"),
                     error = vapply(results, function(x) x$error %||% "",
                                    character(1L)),
                     stringsAsFactors = FALSE)
new_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
assert(all(status$ok), "Branch evaluation errors retained; no selective rerun.")
science <- add_recovery_flags(do.call(rbind, lapply(results, `[[`, "compact")))
science$stage7i_private_resume <- TRUE
science$source_sweep <- 50L
science$full_structural_recovery <- science$core_truth_recovered &
  science$exact_thresholded_factor_structure &
  science$full_true_feature_component_coverage
assert(all(science$truth_global_direction_count == 3L) &&
         all(science$truth_study_role_count == 4L) &&
         all(science$feature_components_total == 8L),
       "Truth structure count differs from the registered Stage 7C case.")
new_csv(science, file.path(evaluation_root, "STAGE7I_ENDPOINT_SCIENCE.csv"))
roles <- do.call(rbind, lapply(names(results), function(branch) {
  cbind(branch = branch, results[[branch]]$role)
}))
new_csv(roles, file.path(evaluation_root, "STAGE7I_ROLE_DAMAGE.csv"))
metric_names <- intersect(c("final_ordinary_elbo", "obs_full_signal_nrmse_selected",
  "dense_full_signal_nrmse_selected", "obs_residual_signal_nrmse_selected",
  "dense_residual_signal_nrmse_selected", "beta_relative_ise_when_identifiable",
  "mean_dense_nrmse", "complete_contribution_nrmse_mean_missing_as_zero",
  "R_A", "R_B1", "R_B2", "L_B1_to_A", "L_B2_to_A", "L_A_to_B1",
  "L_A_to_B2"), names(science))
pairs <- list(c("A", "U"), c("S", "U"), c("S", "A"))
paired <- do.call(rbind, lapply(pairs, function(pair) {
  left <- science[science$branch == pair[[1L]], , drop = FALSE]
  right <- science[science$branch == pair[[2L]], , drop = FALSE]
  row <- data.frame(comparison = paste(pair, collapse = "-"),
                    left_core = left$core_truth_recovered,
                    right_core = right$core_truth_recovered,
                    left_exact = left$exact_thresholded_factor_structure,
                    right_exact = right$exact_thresholded_factor_structure,
                    left_full_coverage = left$full_true_feature_component_coverage,
                    right_full_coverage = right$full_true_feature_component_coverage)
  for (name in metric_names) {
    if (is.numeric(left[[name]]) && is.numeric(right[[name]])) {
      row[[paste0("delta_", name)]] <- left[[name]] - right[[name]]
    }
  }
  row
}))
new_csv(paired, file.path(evaluation_root, "STAGE7I_PAIRED_DIFFERENCES.csv"))
new_lines(c("status=STAGE7I_ENDPOINT_SCIENCE_COMPLETE",
            "branches=U,A,S", "u_baseline_science_parity=PASS",
            "factor_ppi_threshold=0.5", "truth_used_for_initialization=FALSE",
            "truth_used_for_endpoint_selection=FALSE",
            "formal_candidate_pool=FALSE",
            "diagnostic_case_selected_after_prior_truth_review=TRUE",
            "new_fit_or_continuation_started=FALSE",
            paste0("evaluator_sha256=", sha256(script_file)),
            paste0("science_sha256=", sha256(file.path(
              evaluation_root, "STAGE7I_ENDPOINT_SCIENCE.csv")))),
          file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
cat("STAGE7I_ENDPOINT_SCIENCE_COMPLETE branches=U,A,S\n")
