#!/usr/bin/env Rscript

# Freeze the three Stage 7I diagnostic endpoints before any truth is opened.
# Usage: Rscript freeze_stage7i_blind_endpoints_20260929_v1.R RESULTS_PARENT
# The U and A/S roots must have the fixed names used by the Stage 7I supervisor.
options(warn = 1, stringsAsFactors = FALSE)

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) {
  stop("Usage: Rscript freeze_stage7i_blind_endpoints_20260929_v1.R RESULTS_PARENT",
       call. = FALSE)
}
assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
results_parent <- normalizePath(argv[[1L]], mustWork = TRUE)
u_root <- file.path(results_parent, "stage7i_u_20260929_v1")
as_root <- file.path(results_parent, "stage7i_as_20260929_v1")
freeze_root <- file.path(results_parent, "stage7i_blind_freeze_20260929_v1")
file_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
assert(length(file_option) == 1L, "Cannot identify freeze script source.")
this_script <- normalizePath(sub("^--file=", "", file_option),
                             mustWork = TRUE)
assert(!file.exists(freeze_root) && !dir.exists(freeze_root),
       paste0("Refusing existing freeze output: ", freeze_root))

reject_truth_path <- function(path) {
  normalized <- chartr("\\", "/", as.character(path))
  existing <- file.exists(path)
  if (any(existing)) {
    normalized <- c(normalized,
                    normalizePath(path[existing], winslash = "/",
                                  mustWork = TRUE))
  }
  assert(all(!grepl("(^|/)sealed_truth(/|$)|(^|/)truth_bundle[.]rds$",
                    normalized, ignore.case = TRUE)),
         "Truth files and sealed-truth directories are forbidden in blind freezing.")
  invisible(TRUE)
}
reject_truth_path(results_parent)
sha256 <- function(path) {
  reject_truth_path(path)
  assert(file.exists(path) && !dir.exists(path),
         paste0("Missing bound file: ", path))
  lines <- system2("sha256sum", shQuote(path),
                   stdout = TRUE, stderr = TRUE)
  assert(length(lines) == 1L &&
           grepl("^[0-9a-f]{64}[[:space:]]", lines[[1L]]),
         paste0("sha256sum failed: ", path))
  substr(lines[[1L]], 1L, 64L)
}
new_csv <- function(x, path) {
  assert(!file.exists(path) && !dir.exists(path),
         paste0("Refusing overwrite: ", path))
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
new_lines <- function(x, path) {
  assert(!file.exists(path) && !dir.exists(path),
         paste0("Refusing overwrite: ", path))
  writeLines(x, path, useBytes = TRUE)
}
first_line <- function(path, expected) {
  assert(file.exists(path) && !dir.exists(path),
         paste0("Missing completion marker: ", path))
  lines <- readLines(path, warn = FALSE)
  assert(length(lines) >= 1L && identical(lines[[1L]], expected),
         paste0("Invalid completion marker: ", path))
  lines
}
relative_path <- function(path) {
  root <- paste0(normalizePath(results_parent, winslash = "/"), "/")
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  assert(startsWith(resolved, root),
         paste0("Output path escaped the Stage 7I results parent: ", path))
  substring(resolved, nchar(root) + 1L)
}
check_binding <- function(path, required_roles) {
  reject_truth_path(path)
  tab <- utils::read.csv(path, stringsAsFactors = FALSE,
                         check.names = FALSE)
  assert(identical(names(tab), c("role", "path", "sha256")) &&
           nrow(tab) > 0L && !anyDuplicated(tab$role) &&
           all(grepl("^[0-9a-f]{64}$", tab$sha256)),
         paste0("Invalid binding schema: ", path))
  assert(setequal(tab$role, required_roles),
         paste0("Binding roles changed: ", path))
  reject_truth_path(tab$path)
  assert(all(vapply(seq_len(nrow(tab)), function(i) {
    identical(sha256(tab$path[[i]]), tab$sha256[[i]])
  }, logical(1L))), paste0("Bound file changed: ", path))
  tab
}
check_parity <- function(path, expected_rows = NULL) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  assert("pass" %in% names(x) && nrow(x) > 0L &&
           all(!is.na(x$pass)) && all(as.logical(x$pass)) &&
           (is.null(expected_rows) || nrow(x) == expected_rows),
         paste0("U parity table failed or changed: ", path))
  x
}
compare_state_value <- function(actual, reference, zero_exact = FALSE) {
  failed <- list(pass = FALSE, max_abs = Inf)
  if (is.null(actual) || is.null(reference)) {
    return(if (is.null(actual) && is.null(reference))
      list(pass = TRUE, max_abs = 0) else failed)
  }
  if (is.list(actual) || is.list(reference)) {
    if (!is.list(actual) || !is.list(reference) ||
        length(actual) != length(reference) ||
        !identical(names(actual), names(reference))) return(failed)
    if (!length(actual)) return(list(pass = TRUE, max_abs = 0))
    parts <- Map(function(a, b) compare_state_value(a, b, zero_exact),
                 actual, reference)
    return(list(pass = all(vapply(parts, `[[`, logical(1L), "pass")),
                max_abs = max(vapply(parts, `[[`, numeric(1L), "max_abs"))))
  }
  if (!identical(dim(actual), dim(reference)) ||
      !identical(dimnames(actual), dimnames(reference)) ||
      !identical(names(actual), names(reference)) ||
      length(actual) != length(reference)) return(failed)
  if (is.numeric(actual) && is.numeric(reference)) {
    if (!identical(is.na(actual), is.na(reference)) ||
        any(!is.finite(actual[!is.na(actual)])) ||
        any(!is.finite(reference[!is.na(reference)]))) return(failed)
    valid <- !is.na(actual)
    difference <- abs(as.numeric(actual[valid]) -
                      as.numeric(reference[valid]))
    allowance <- 1e-8 + 1e-10 * abs(as.numeric(reference[valid]))
    zero_ok <- !zero_exact ||
      identical(as.logical(actual[valid] == 0),
                as.logical(reference[valid] == 0))
    return(list(pass = all(difference <= allowance) && zero_ok,
                max_abs = max(difference, 0)))
  }
  if (identical(actual, reference)) list(pass = TRUE, max_abs = 0) else failed
}
recheck_u_snapshots <- function(actual_path, reference_path) {
  reject_truth_path(c(actual_path, reference_path))
  actual_sha_before <- sha256(actual_path)
  reference_sha_before <- sha256(reference_path)
  actual <- readRDS(actual_path)
  reference <- readRDS(reference_path)
  assert(length(actual) == 449L && length(reference) == 500L &&
           identical(vapply(actual, `[[`, integer(1L), "iteration"),
                     as.integer(51:499)) &&
           identical(vapply(reference, `[[`, integer(1L), "iteration"),
                     as.integer(0:499)),
         "U or Stage 7G raw snapshot index is incomplete.")
  rows <- vector("list", 449L * length(reference[[1L]]))
  at <- 0L
  for (iteration in 51:499) {
    current <- actual[[iteration - 50L]]
    target <- reference[[iteration + 1L]]
    assert(identical(names(current), names(target)),
           paste0("U snapshot fields changed at sweep ", iteration))
    for (field in names(target)) {
      result <- compare_state_value(
        current[[field]], target[[field]],
        field %in% c("mu_q_gamma_a", "mu_q_gamma_b"))
      at <- at + 1L
      rows[[at]] <- data.frame(
        iteration = as.integer(iteration), field = field,
        max_abs_difference = result$max_abs, pass = result$pass,
        stringsAsFactors = FALSE)
      assert(result$pass,
             paste0("U raw snapshot recheck failed at sweep ",
                    iteration, ", field ", field))
    }
  }
  assert(identical(sha256(actual_path), actual_sha_before) &&
           identical(sha256(reference_path), reference_sha_before),
         "U or Stage 7G snapshot file changed during recheck.")
  rm(actual, reference)
  invisible(gc(verbose = FALSE))
  do.call(rbind, rows[seq_len(at)])
}

assert(dir.exists(u_root) && dir.exists(as_root),
       "The fixed U and A/S output roots are required.")
u_pass_path <- file.path(u_root, "u_parity_v1", "PARITY_PASS.txt")
u_pass <- first_line(u_pass_path, "STAGE7I_U_PARITY_PASS")
u_complete_path <- file.path(u_root, "u_full_v1", "FIT_COMPLETE.txt")
first_line(u_complete_path, "STAGE7I_U_FIT_COMPLETE")
as_complete_path <- file.path(as_root, "AS_FITS_COMPLETE.txt")
as_complete <- first_line(as_complete_path, "STAGE7I_AS_FITS_COMPLETE")
assert("branches=A,S" %in% as_complete,
       "A/S completion marker does not cover both branches.")

u_binding_path <- file.path(u_root, "INPUT_AND_SOURCE_BINDING.csv")
u_roles <- c("runner", "stage7i_archive", "stage7i_description",
             "stage7i_namespace", "stage7i_rdb", "stage7i_rdx",
             "stage7g_runner", "stage7g_archive", "stage7g_manifest",
             "stage7g_input_binding", "stage7h_blind_manifest",
             "observation", "stage7g_fit", "stage7g_snapshots")
u_binding <- check_binding(u_binding_path, u_roles)
u_paths <- stats::setNames(u_binding$path, u_binding$role)
fixed_input_hashes <- c(
  stage7i_archive =
    "b6c2210089061d99e31825e358d9c03b6d867bcd7fdde21e8432e784e4436178",
  stage7g_runner =
    "668788f6cb159badedb15c0630d0754019e2de18aaaa79cf61bff60eb3acd935",
  stage7g_archive =
    "b3e9050f1de7154baff28b7dd986254df0b634bf85d6f50f464d00dbd4c8f17f",
  stage7g_manifest =
    "20a68304cccb3f4176c046b4c5c4f4deb4ca81ffd6d10c99ccc5dfdd20ad4e7a",
  stage7g_input_binding =
    "1d0e6ae8badf33edbb9b2578e615d395f00f5501a7aad2fce3db0f456ce62a54",
  stage7h_blind_manifest =
    "d67bef6e692fa87296f3d0b0e19a6a954d05a0085bd3f3ec4b08350b3f0f3f98",
  observation =
    "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6",
  stage7g_fit =
    "b2166588d932c4432261e026745f1d1d8185dc0b00f741125d6a779de5d7d761",
  stage7g_snapshots =
    "e85ce3fdd8b2770674232e11cdab5453f61b30271977a41cb08d54b6383621c2"
)
assert(all(vapply(names(fixed_input_hashes), function(role) {
  identical(u_binding$sha256[u_binding$role == role],
            unname(fixed_input_hashes[[role]]))
}, logical(1L))), "Frozen input hashes no longer match Stage 7I.")
reject_truth_path(u_binding$path)

parity_dir <- file.path(u_root, "u_parity_v1")
parity_paths <- file.path(parity_dir, c(
  "ALL_STATE_PARITY.csv", "ORDINARY_ELBO_PARITY.csv",
  "ENDPOINT_PARITY.csv", "CONFIGURATION_PARITY.csv",
  "DRIVER_TRACE_PARITY.csv"))
state_parity <- check_parity(parity_paths[[1L]])
assert(identical(sort(unique(as.integer(state_parity$iteration))),
                 as.integer(51:499)),
       "U parity does not span sweeps 51:499.")
check_parity(parity_paths[[2L]], 400L)
for (path in parity_paths[3:5]) check_parity(path)
u_snapshot_recheck <- recheck_u_snapshots(
  file.path(u_root, "u_full_v1", "RAW_STATE_SNAPSHOTS_51_499.rds"),
  u_paths[["stage7g_snapshots"]])

as_source_path <- file.path(as_root, "SOURCE_BINDING.csv")
as_injected_path <- file.path(as_root, "INJECTED_STATE_BINDING.csv")
as_source_roles <- c("as_runner", "initializer", "protocol",
                     "u_binding", "u_parity_pass", "u_state_parity",
                     "u_elbo_parity", "u_endpoint_parity",
                     "u_config_parity", "u_trace_parity", "u_fit")
as_injected_roles <- c("U", "A", "S", "diagnostics", "summary",
                       "driver_coordinates", "driver_subjects")
as_source <- check_binding(as_source_path, as_source_roles)
as_injected <- check_binding(as_injected_path, as_injected_roles)
assert(identical(as_source$sha256[as_source$role == "u_binding"],
                 sha256(u_binding_path)) &&
         identical(as_source$sha256[as_source$role == "u_parity_pass"],
                   sha256(u_pass_path)) &&
         identical(as_source$sha256[as_source$role == "u_fit"],
                   sha256(file.path(u_root, "u_full_v1", "fit.rds"))),
       "A/S source binding is not tied to this U parity and fit.")
assert(setequal(as_injected$role, as_injected_roles),
       "Injected-state binding is incomplete.")
first_line(file.path(as_root, "INJECTION_COMPLETE.txt"),
           "STAGE7I_AS_INJECTION_COMPLETE")

branches <- c("U", "A", "S")
branch_roots <- c(U = file.path(u_root, "u_full_v1"),
                  A = file.path(as_root, "A_full_v1"),
                  S = file.path(as_root, "S_full_v1"))
branch_files <- vector("list", length(branches))
names(branch_files) <- branches
reference_config <- NULL
reference_index <- NULL
manifest <- vector("list", length(branches))
metrics <- vector("list", length(branches))

for (branch in branches) {
  branch_root <- branch_roots[[branch]]
  fit_path <- file.path(branch_root, "fit.rds")
  snapshot_path <- file.path(branch_root,
                             "RAW_STATE_SNAPSHOTS_51_499.rds")
  index_path <- file.path(branch_root, "RAW_STATE_INDEX.csv")
  elbo_path <- file.path(branch_root, "ORDINARY_ELBO.csv")
  driver_path <- file.path(branch_root, "DRIVER_TRACE.csv")
  warnings_path <- file.path(branch_root, "WARNINGS.txt")
  completion_path <- file.path(branch_root, "FIT_COMPLETE.txt")
  expected_marker <- if (branch == "U")
    "STAGE7I_U_FIT_COMPLETE" else "STAGE7I_AS_FIT_COMPLETE"
  marker <- first_line(completion_path, expected_marker)
  if (branch != "U") {
    assert(paste0("branch=", branch) %in% marker,
           paste0("Wrong branch completion marker: ", branch))
    binding_path <- file.path(branch_root, "BRANCH_OUTPUT_BINDING.csv")
    expected_roles <- c("fit", "snapshots", "state_index",
                        "driver_trace", "elbo", "warnings")
    if (file.exists(file.path(branch_root, "SCALE_TRACE.csv"))) {
      expected_roles <- c(expected_roles, "scale_trace")
    }
    binding <- check_binding(binding_path, expected_roles)
    assert(setequal(binding$role, expected_roles) &&
             paste0("binding_sha256=", sha256(binding_path)) %in% marker,
           paste0("Invalid branch output binding: ", branch))
    expected_paths <- c(fit = fit_path, snapshots = snapshot_path,
                        state_index = index_path, driver_trace = driver_path,
                        elbo = elbo_path, warnings = warnings_path)
    if ("scale_trace" %in% binding$role) {
      expected_paths <- c(expected_paths,
                          scale_trace = file.path(branch_root,
                                                  "SCALE_TRACE.csv"))
    }
    assert(all(vapply(names(expected_paths), function(role) {
      identical(normalizePath(binding$path[binding$role == role]),
                normalizePath(expected_paths[[role]]))
    }, logical(1L))),
    paste0("A/S bound output path changed: ", branch))
    branch_files[[branch]] <- c(unname(expected_paths), binding_path,
                                completion_path)
  } else {
    branch_files[[branch]] <- c(fit_path, snapshot_path, index_path,
                                driver_path, elbo_path, warnings_path,
                                completion_path)
    scale_path <- file.path(branch_root, "SCALE_TRACE.csv")
    if (file.exists(scale_path)) {
      branch_files[[branch]] <- c(branch_files[[branch]], scale_path)
    }
    invisible(vapply(branch_files[[branch]], sha256, character(1L)))
  }

  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  assert(identical(as.integer(index$iteration), as.integer(51:499)) &&
           nrow(index) == 449L &&
           all(index$annealing[1:49]) &&
           !any(index$annealing[50:449]) &&
           all(is.finite(index$temperature)) &&
           all(index$temperature[1:49] > 1) &&
           all(abs(index$temperature[50:449] - 1) <= 1e-12),
         paste0("Invalid fixed temperature schedule: ", branch))
  if (is.null(reference_index)) {
    reference_index <- index
  } else {
    assert(identical(index, reference_index),
           paste0("Branch temperature/index differs from U: ", branch))
  }

  fit <- readRDS(fit_path)
  config <- list(
    dimensions = fit[c("S", "n_s", "p", "d", "L_f", "L_s",
                       "M_f", "M_s", "K", "n_g")],
    seed = fit$fit_seed,
    fit_control = fit$fit_control,
    pre_score_interface = fit$pre_score_interface,
    hyper = fit$list_hyper,
    omega_branch = fit$bool_var_spec_prob,
    time_g = fit$time_g,
    response_center = fit$response_center,
    response_scale = fit$response_scale
  )
  if (is.null(reference_config)) {
    reference_config <- config
  } else {
    assert(identical(config, reference_config),
           paste0("Data, dimensions, budget or stopping rules differ: ",
                  branch))
  }
  assert(identical(as.integer(fit$fit_seed), 99716106L) &&
           fit$i_iter == 499L && fit$annealing_sweeps == 99L &&
           fit$t1_sweeps == 400L && length(fit$ELBO) == 400L &&
           all(is.finite(fit$ELBO)) &&
           isTRUE(fit$elbo_objective_valid),
         paste0("Invalid fixed-budget endpoint: ", branch))
  elbo_table <- utils::read.csv(elbo_path, stringsAsFactors = FALSE)
  assert(nrow(elbo_table) == 400L &&
           identical(as.integer(elbo_table$t1_sweep), as.integer(1:400)) &&
           all(is.finite(elbo_table$ordinary_elbo)) &&
           max(abs(as.numeric(elbo_table$ordinary_elbo) -
                   as.numeric(fit$ELBO))) <= 1e-6,
         paste0("ELBO table differs from fit: ", branch))

  b1 <- as.numeric(fit$mu_q_b_specific[[1L]])
  b2 <- as.numeric(fit$mu_q_b_specific[[2L]])
  rho_a <- as.numeric(fit$mu_q_gamma_a)
  rho_b1 <- as.numeric(fit$mu_q_gamma_b[[1L]])
  rho_b2 <- as.numeric(fit$mu_q_gamma_b[[2L]])
  metrics[[branch]] <- data.frame(
    branch = branch, source_fit_id = "covrb7c_01__C1__start06",
    data_id = "covrb7c_01__C1", source_sweep = 50L,
    total_iteration = as.integer(fit$i_iter),
    annealing_sweeps = as.integer(fit$annealing_sweeps),
    ordinary_sweeps = as.integer(fit$t1_sweeps),
    final_ordinary_elbo = as.numeric(tail(fit$ELBO, 1L)),
    objective_valid = isTRUE(fit$elbo_objective_valid),
    convergence_status = as.character(fit$convergence_status),
    practical_converged = isTRUE(fit$practical_converged),
    t1_jitter_count = as.integer(fit$elbo_t1_jitter_count),
    mean_ppi_shared = mean(rho_a),
    mean_ppi_B1 = mean(rho_b1),
    mean_ppi_B2 = mean(rho_b2),
    effective_loading_norm_shared =
      sqrt(sum(as.numeric(fit$mu_q_a)^2)),
    effective_loading_norm_B1 = sqrt(sum(b1^2)),
    effective_loading_norm_B2 = sqrt(sum(b2^2)),
    truth_used_for_fit_or_selection = FALSE,
    stringsAsFactors = FALSE)
  assert(all(is.finite(unlist(metrics[[branch]][,
    vapply(metrics[[branch]], is.numeric, logical(1L)), drop = FALSE]))),
    paste0("Non-finite blind metric: ", branch))
  metrics_path <- file.path(freeze_root,
                            paste0("BLIND_METRICS_", branch, ".csv"))
  manifest[[match(branch, branches)]] <- data.frame(
    branch = branch,
    source_fit_id = "covrb7c_01__C1__start06",
    data_id = "covrb7c_01__C1",
    fit_path = relative_path(fit_path),
    fit_sha256 = sha256(fit_path),
    snapshot_path = relative_path(snapshot_path),
    snapshot_sha256 = sha256(snapshot_path),
    blind_metrics_path = file.path(basename(freeze_root),
                                   basename(metrics_path)),
    blind_metrics_sha256 = NA_character_,
    status = "COMPLETE",
    final_ordinary_elbo = as.numeric(tail(fit$ELBO, 1L)),
    output_eligible = TRUE,
    stringsAsFactors = FALSE)
  rm(fit)
  invisible(gc(verbose = FALSE))
}

assert(identical(vapply(manifest, `[[`, character(1L), "branch"),
                 branches), "Branch order changed.")
assert(dir.create(freeze_root, recursive = FALSE, showWarnings = FALSE),
       "Cannot create blind freeze output directory.")
u_snapshot_recheck_path <- file.path(freeze_root,
                                    "U_SNAPSHOT_RECHECK.csv")
new_csv(u_snapshot_recheck, u_snapshot_recheck_path)
env_path <- file.path(freeze_root, "R_ENVIRONMENT.txt")
new_lines(c(paste0("R_version=", R.version.string),
            paste0("platform=", R.version$platform),
            paste0("R_home=", R.home()),
            capture.output(utils::sessionInfo())), env_path)
metric_paths <- stats::setNames(
  file.path(freeze_root, paste0("BLIND_METRICS_", branches, ".csv")),
  branches)
for (branch in branches) {
  new_csv(metrics[[branch]], metric_paths[[branch]])
  manifest[[match(branch, branches)]]$blind_metrics_sha256 <-
    sha256(metric_paths[[branch]])
}
manifest <- do.call(rbind, manifest)
assert(nrow(manifest) == 3L &&
         identical(as.character(manifest$branch), branches),
       "Blind manifest must contain exactly U/A/S.")
manifest_path <- file.path(freeze_root, "BLIND_ENDPOINT_MANIFEST.csv")
new_csv(manifest, manifest_path)

source_tables <- list(u = u_binding, as = as_source,
                      injected = as_injected)
source_binding <- do.call(rbind, lapply(names(source_tables), function(group) {
  x <- source_tables[[group]]
  data.frame(role = paste0(group, ":", x$role), path = x$path,
             sha256 = x$sha256, stringsAsFactors = FALSE)
}))
source_binding <- rbind(source_binding, data.frame(
  role = c("freeze_script_source", "r_environment",
           "u_input_binding", "as_source_binding", "as_injected_binding"),
  path = c(this_script, env_path, u_binding_path,
           as_source_path, as_injected_path),
  sha256 = vapply(c(this_script, env_path, u_binding_path,
                    as_source_path, as_injected_path),
                  sha256, character(1L)),
  stringsAsFactors = FALSE))
reject_truth_path(source_binding$path)
source_binding_path <- file.path(freeze_root, "BLIND_SOURCE_BINDING.csv")
new_csv(source_binding, source_binding_path)

u_files <- c(u_binding_path, u_pass_path, u_complete_path, parity_paths,
             branch_files[["U"]],
             file.path(u_root, "PREPARE_COMPLETE.txt"),
             file.path(u_root, "smoke_u_v1", "SMOKE_PASS.txt"))
as_files <- c(as_source_path, as_injected_path,
              file.path(as_root, "INJECTION_COMPLETE.txt"),
              as_complete_path, as_injected$path,
              branch_files[["A"]], branch_files[["S"]])
all_files <- unique(c(u_files, as_files, metric_paths,
                      u_snapshot_recheck_path,
                      manifest_path, env_path, source_binding_path))
file_hashes <- data.frame(
  path = vapply(all_files, relative_path, character(1L)),
  sha256 = vapply(all_files, sha256, character(1L)),
  stringsAsFactors = FALSE)
assert(!anyDuplicated(file_hashes$path), "Duplicate blind file path.")
file_hashes_path <- file.path(freeze_root, "BLIND_FILE_HASHES.csv")
new_csv(file_hashes, file_hashes_path)
assert(all(vapply(seq_len(nrow(file_hashes)), function(i) {
  identical(sha256(file.path(results_parent, file_hashes$path[[i]])),
            file_hashes$sha256[[i]])
}, logical(1L))), "Blind file hash verification failed before freeze.")

marker <- c(
  "STAGE7I_BLIND_ENDPOINTS_FROZEN",
  "branch_order=U,A,S", "manifest_rows=3",
  paste0("relative_path_root=", results_parent),
  paste0("manifest_sha256=", sha256(manifest_path)),
  paste0("file_hashes_sha256=", sha256(file_hashes_path)),
  paste0("source_binding_sha256=", sha256(source_binding_path)),
  paste0("metrics_U_sha256=", sha256(metric_paths[["U"]])),
  paste0("metrics_A_sha256=", sha256(metric_paths[["A"]])),
  paste0("metrics_S_sha256=", sha256(metric_paths[["S"]])),
  paste0("u_snapshot_recheck_sha256=",
         sha256(u_snapshot_recheck_path)),
  paste0("u_snapshot_recheck_rows=", nrow(u_snapshot_recheck)),
  paste0("u_parity_pass_sha256=", sha256(u_pass_path)),
  paste0("as_fits_complete_sha256=", sha256(as_complete_path)),
  "truth_used_for_initialization=FALSE",
  "truth_used_for_endpoint_selection=FALSE",
  "formal_candidate_pool=FALSE",
  "diagnostic_case_selected_after_prior_truth_review=TRUE",
  "truth_bundle_opened=FALSE")
new_lines(marker, file.path(freeze_root,
                            "BLIND_ENDPOINTS_FROZEN.txt"))
cat("STAGE7I_BLIND_ENDPOINTS_FROZEN rows=3 files=",
    nrow(file_hashes), "\n", sep = "")
