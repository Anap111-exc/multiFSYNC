#!/usr/bin/env Rscript

# Stage 7J: truth-blind seed02 U continuation from the frozen Stage 7G sweep 50.
# This script never opens a truth bundle. It refuses existing output directories.
# Usage: Rscript this.R ACTION STAGE7C_ROOT STAGE7G_ROOT STAGE7G_RUNNER \
#   STAGE7G_SOURCE_ARCHIVE STAGE7H_ROOT STAGE7I_LIBRARY \
#   STAGE7I_SOURCE_ARCHIVE OUTPUT_ROOT
# ACTION is prepare, u, verify, or status. The actual-data U run is done once:
# a separate two-sweep smoke would exceed the registered seed02 sweep budget.
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 9L || !argv[[1L]] %in%
    c("prepare", "u", "verify", "status")) {
  stop("Usage: Rscript run_stage7j_u_resume_20260930_v1.R ",
       "{prepare|u|verify|status} STAGE7C_ROOT STAGE7G_ROOT ",
       "STAGE7G_RUNNER STAGE7G_SOURCE_ARCHIVE STAGE7H_ROOT ",
       "STAGE7I_LIBRARY STAGE7I_SOURCE_ARCHIVE OUTPUT_ROOT",
       call. = FALSE)
}
action <- argv[[1L]]
stage7c_root <- normalizePath(argv[[2L]], mustWork = TRUE)
stage7g_root <- normalizePath(argv[[3L]], mustWork = TRUE)
stage7g_runner <- normalizePath(argv[[4L]], mustWork = TRUE)
stage7g_archive <- normalizePath(argv[[5L]], mustWork = TRUE)
stage7h_root <- normalizePath(argv[[6L]], mustWork = TRUE)
stage7i_library <- normalizePath(argv[[7L]], mustWork = TRUE)
stage7i_archive <- normalizePath(argv[[8L]], mustWork = TRUE)
output_root <- argv[[9L]]
script_option <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_option) != 1L) stop("Cannot identify runner source path.")
this_script <- normalizePath(sub("^--file=", "", script_option),
                             mustWork = TRUE)

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
sha256 <- function(path) {
  assert(file.exists(path) && !dir.exists(path),
         paste0("Missing bound file: ", path))
  result <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  assert(length(result) == 1L &&
           grepl("^[0-9a-f]{64}[[:space:]]", result[[1L]]),
         paste0("sha256sum failed: ", path))
  substr(result[[1L]], 1L, 64L)
}
write_csv_new <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
write_lines_new <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  writeLines(x, path, useBytes = TRUE)
}
save_rds_new <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  saveRDS(x, path, version = 3L)
}
require_new_dir <- function(path) {
  assert(!file.exists(path) && !dir.exists(path),
         paste0("Refusing existing output directory: ", path))
  assert(dir.create(path, recursive = TRUE, showWarnings = FALSE),
         paste0("Cannot create output directory: ", path))
}

fit_id <- "covrb7c_02__C1__start11"
data_id <- "covrb7c_02__C1"
fit_seed <- 99716211L
manifest_path <- file.path(stage7g_root, "FIT_MANIFEST.csv")
input_binding_path <- file.path(stage7g_root, "analysis_v1", "INPUT_BINDING.csv")
observation_path <- file.path(stage7c_root, "data", data_id,
                              "observation_bundle.rds")
stage7g_fit_path <- file.path(stage7g_root, "fits", fit_id, "fit.rds")
stage7g_snapshot_path <- file.path(stage7g_root, "fits", fit_id,
                                   "RAW_STATE_SNAPSHOTS.rds")
stage7h_blind_path <- file.path(stage7h_root, "blind_all_v3",
                               "BLIND_TABLES.sha256")
installed_package <- file.path(stage7i_library, "multiFSYNC")
installed_files <- file.path(installed_package,
                             c("DESCRIPTION", "NAMESPACE",
                               "R/multiFSYNC.rdb", "R/multiFSYNC.rdx"))

fixed_hashes <- c(
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
    "e9bff47094381f3005fcb7d29e56d0854bb436fafe4ebc5cd9d4fecc4bc6ae32",
  stage7g_fit =
    "9f73a1a162ec04495d2caf398e4b6f25956bd51bd68b910b28646aebd471d0f2",
  stage7g_snapshots =
    "d98a8e755716a79ff78b6b3c9b67511cfb38287a926eba6b55434893a4abd858"
)
assert(all(grepl("^[0-9a-f]{64}$", fixed_hashes)),
       "A fixed source hash literal is malformed.")

binding_paths <- c(
  runner = this_script,
  stage7i_archive = stage7i_archive,
  stage7i_description = installed_files[[1L]],
  stage7i_namespace = installed_files[[2L]],
  stage7i_rdb = installed_files[[3L]],
  stage7i_rdx = installed_files[[4L]],
  stage7g_runner = stage7g_runner,
  stage7g_archive = stage7g_archive,
  stage7g_manifest = manifest_path,
  stage7g_input_binding = input_binding_path,
  stage7h_blind_manifest = stage7h_blind_path,
  observation = observation_path,
  stage7g_fit = stage7g_fit_path,
  stage7g_snapshots = stage7g_snapshot_path
)

make_binding <- function() {
  data.frame(role = names(binding_paths),
             path = unname(vapply(binding_paths, normalizePath,
                                  character(1L), mustWork = TRUE)),
             sha256 = unname(vapply(binding_paths, sha256, character(1L))),
             stringsAsFactors = FALSE)
}
check_frozen_contract <- function() {
  current <- make_binding()
  for (role in names(fixed_hashes)) {
    actual <- current$sha256[current$role == role]
    assert(length(actual) == 1L && identical(actual, fixed_hashes[[role]]),
           paste0("Frozen SHA256 mismatch: ", role))
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  assert(nrow(manifest) == 3L && nrow(row) == 1L &&
           identical(row$data_id[[1L]], data_id) &&
           identical(as.integer(row$fit_seed[[1L]]), fit_seed) &&
           identical(row$observation_sha256[[1L]],
                     fixed_hashes[["observation"]]) &&
           identical(row$stage7f_fit_sha256[[1L]],
                     fixed_hashes[["stage7g_fit"]]),
         "Stage 7G manifest does not bind the selected C1 winner.")
  input <- utils::read.csv(input_binding_path, stringsAsFactors = FALSE)
  row <- input[input$fit_id == fit_id, , drop = FALSE]
  assert(nrow(row) == 1L &&
           identical(row$fit_sha256[[1L]], fixed_hashes[["stage7g_fit"]]) &&
           identical(row$snapshot_sha256[[1L]],
                     fixed_hashes[["stage7g_snapshots"]]),
         "Stage 7G input binding changed.")
  current
}
check_prepare_binding <- function() {
  assert(file.exists(file.path(output_root, "PREPARE_COMPLETE.txt")),
         "Prepare marker missing.")
  file <- file.path(output_root, "INPUT_AND_SOURCE_BINDING.csv")
  assert(file.exists(file), "Input/source binding table missing.")
  previous <- utils::read.csv(file, stringsAsFactors = FALSE)
  current <- check_frozen_contract()
  assert(identical(previous, current),
         "Stage 7I input/source paths or hashes changed after prepare.")
  current
}
load_stage7i_package <- function() {
  .libPaths(c(stage7i_library, .libPaths()))
  suppressPackageStartupMessages(library(multiFSYNC))
  assert(identical(normalizePath(system.file(package = "multiFSYNC")),
                   normalizePath(installed_package)) &&
           identical(as.character(utils::packageVersion("multiFSYNC")),
                     "0.3.0.9000"),
         "Loaded package is not the bound Stage 7I build.")
}
read_observation <- function() {
  x <- readRDS(observation_path)
  assert(identical(x$bundle_class,
                   "covariate_routeb_pilot_observation_only") &&
           isFALSE(x$truth_available_to_fit) &&
           identical(x$data_id, data_id) &&
           x$dimensions$p == 500L && x$dimensions$d == 2L &&
           x$dimensions$K == 9L,
         "Observation-only input contract failed.")
  x
}

# The private driver trace is diagnostic only. Capture every resumed sweep;
# compare the overlap with Stage 7G by its row keys, excluding sequence.
trace_sweeps <- as.integer(51:499)
scale_trace_control <- list(
  include_initial = FALSE, include_annealing = TRUE,
  t1_sweeps = as.integer(c(1:20, 50, 100, 200, 300, 380, 400)),
  block_annealing_sweeps = as.integer(51:60)
)
fit_arguments <- function(observation) {
  control <- multiFSYNC::g12_stopping_control()
  control$checkpoints <- c(380L, 400L)
  list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L, M_s = list(2L, 2L),
    K = observation$dimensions$K, n_g = 51L,
    anneal = c(1, 1.9, 100), maxit = 499L,
    n_cpus = 1L, verbose = FALSE, seed = fit_seed,
    tol_abs = 1e-3, tol_rel = 1e-5,
    bool_scale = FALSE, bool_var_spec_prob = FALSE,
    d_0 = as.integer(observation$dimensions$p),
    convergence_rule = "practical", lambda_orth = 0,
    practical_control = control, initialization = "random",
    initialization_control = list(grid_size = 81L,
                                  perturb_sd = 0.05, rank_tol = 1e-8),
    continuation_state = NULL, pre_score_sweeps = 1L,
    trace_sweeps = trace_sweeps,
    function_initialization = "gram_unit_energy"
  )
}
reference_snapshots <- function() {
  x <- readRDS(stage7g_snapshot_path)
  assert(is.list(x) && length(x) == 500L &&
           identical(vapply(x, `[[`, integer(1L), "iteration"),
                     as.integer(0:499)),
         "Frozen Stage 7G snapshots must cover 0:499 exactly.")
  x
}
new_snapshot_sink <- function(sweeps) {
  x <- new.env(parent = emptyenv())
  x$schema <- "stage7g_snapshot_v1"
  x$sweeps <- as.integer(sweeps)
  x$records <- list()
  x
}
run_u <- function(last_iteration, snapshot_sweeps) {
  source <- reference_snapshots()
  resume <- new.env(parent = emptyenv())
  resume$schema <- "stage7i_sweep50_resume_v1"
  resume$state <- source[[51L]]
  rm(source)
  invisible(gc(verbose = FALSE))
  resume$stop_after_iteration <- as.integer(last_iteration)
  sink <- new_snapshot_sink(snapshot_sweeps)
  old_options <- options(multiFSYNC.__stage7i_resume_20260929__ = resume,
                         multiFSYNC.__stage7g_snapshot_20260928__ = sink,
                         multiFSYNC.__private_scale_trace_8d72c46f__ =
                           scale_trace_control)
  on.exit(options(old_options), add = TRUE)
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(
    do.call(multiFSYNC::bayesSYNC_multi_pre_score,
            fit_arguments(read_observation())),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))
  ), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  list(fit = fit, snapshots = sink$records, warnings = warnings,
       elapsed = elapsed, applied = isTRUE(resume$applied),
       first_temperature = resume$first_temperature)
}

compare_value <- function(actual, reference, zero_exact = FALSE) {
  failed <- list(pass = FALSE, max_abs = Inf, max_excess = Inf,
                 zero_pattern_equal = FALSE)
  if (is.null(actual) || is.null(reference)) {
    if (is.null(actual) && is.null(reference)) {
      return(list(pass = TRUE, max_abs = 0, max_excess = 0,
                  zero_pattern_equal = TRUE))
    }
    return(failed)
  }
  if (is.list(actual) || is.list(reference)) {
    if (!is.list(actual) || !is.list(reference) ||
        length(actual) != length(reference) ||
        !identical(names(actual), names(reference))) return(failed)
    if (!length(actual)) {
      return(list(pass = TRUE, max_abs = 0, max_excess = 0,
                  zero_pattern_equal = TRUE))
    }
    parts <- Map(function(a, b) compare_value(a, b, zero_exact),
                 actual, reference)
    return(list(pass = all(vapply(parts, `[[`, logical(1L), "pass")),
                max_abs = max(vapply(parts, `[[`, numeric(1L), "max_abs")),
                max_excess = max(vapply(parts, `[[`, numeric(1L),
                                        "max_excess")),
                zero_pattern_equal = all(vapply(parts, `[[`, logical(1L),
                                                "zero_pattern_equal"))))
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
    excess <- max(difference - allowance, 0)
    return(list(pass = excess == 0 && zero_ok,
                max_abs = max(difference, 0), max_excess = excess,
                zero_pattern_equal = zero_ok))
  }
  if (identical(actual, reference)) {
    list(pass = TRUE, max_abs = 0, max_excess = 0,
         zero_pattern_equal = TRUE)
  } else failed
}
compare_snapshots <- function(actual, reference, iterations) {
  assert(length(actual) == length(iterations) &&
           identical(vapply(actual, `[[`, integer(1L), "iteration"),
                     as.integer(iterations)),
         "Resume snapshot index is incomplete or out of order.")
  out <- vector("list", length(iterations) * length(reference[[1L]]))
  at <- 0L
  for (k in seq_along(iterations)) {
    iteration <- iterations[[k]]
    current <- actual[[k]]
    target <- reference[[iteration + 1L]]
    assert(identical(names(current), names(target)),
           paste0("Snapshot field names differ at sweep ", iteration))
    for (field in names(target)) {
      result <- compare_value(current[[field]], target[[field]],
                              field %in% c("mu_q_gamma_a", "mu_q_gamma_b"))
      at <- at + 1L
      out[[at]] <- data.frame(
        iteration = iteration, field = field,
        max_abs_difference = result$max_abs,
        max_tolerance_excess = result$max_excess,
        ppi_zero_pattern_equal = result$zero_pattern_equal,
        pass = result$pass
      )
    }
  }
  do.call(rbind, out[seq_len(at)])
}
write_state_index <- function(snapshots, path) {
  write_csv_new(data.frame(
    iteration = vapply(snapshots, `[[`, integer(1L), "iteration"),
    temperature = vapply(snapshots, `[[`, numeric(1L), "temperature"),
    annealing = vapply(snapshots, `[[`, logical(1L), "annealing")
  ), path)
}

if (action == "prepare") {
  source_roots <- c(stage7c_root, stage7g_root, stage7h_root,
                    stage7i_library)
  proposed <- normalizePath(output_root, mustWork = FALSE)
  assert(!any(vapply(source_roots, function(root)
    identical(proposed, root) || startsWith(proposed,
                                            paste0(root, .Platform$file.sep)),
    logical(1L))), "Output root overlaps a frozen input or installed library.")
  binding <- check_frozen_contract()
  load_stage7i_package()
  read_observation()
  source <- reference_snapshots()
  assert(identical(source[[51L]]$iteration, 50L),
         "Source snapshot 50 is missing.")
  require_new_dir(output_root)
  write_csv_new(binding, file.path(output_root,
                                   "INPUT_AND_SOURCE_BINDING.csv"))
  write_lines_new(c(
    "STAGE7J_U_PREPARE_COMPLETE", paste0("fit_id=", fit_id),
    "source_iteration=50", "requested_sweeps=51:499",
    "truth_opened=FALSE", "formal_candidate_pool=FALSE",
    paste0("R_version=", getRversion()),
    paste0("stage7i_package_version=",
           utils::packageVersion("multiFSYNC"))
  ), file.path(output_root, "PREPARE_COMPLETE.txt"))
  cat("STAGE7J_U_PREPARE_COMPLETE\n")
}

if (action == "u") {
  check_prepare_binding()
  load_stage7i_package()
  u_dir <- file.path(output_root, "u_full_v1")
  require_new_dir(u_dir)
  result <- run_u(499L, 51:499)
  write_lines_new(result$warnings, file.path(u_dir, "WARNINGS.txt"))
  save_rds_new(result$snapshots,
               file.path(u_dir, "RAW_STATE_SNAPSHOTS_51_499.rds"))
  write_state_index(result$snapshots,
                    file.path(u_dir, "RAW_STATE_INDEX.csv"))
  if (inherits(result$fit, "error")) {
    write_lines_new(c("STAGE7J_U_FIT_ERROR",
                      conditionMessage(result$fit)),
                    file.path(u_dir, "FIT_ERROR.txt"))
    stop("Full U continuation failed; retained partial snapshots.")
  }
  fit <- result$fit
  save_rds_new(fit, file.path(u_dir, "fit.rds"))
  write_csv_new(fit$driver_trace,
                file.path(u_dir, "DRIVER_TRACE.csv"))
  if (!is.null(fit$scale_trace)) {
    write_csv_new(fit$scale_trace, file.path(u_dir, "SCALE_TRACE.csv"))
  }
  write_csv_new(data.frame(t1_sweep = seq_along(fit$ELBO),
                           ordinary_elbo = as.numeric(fit$ELBO)),
                file.path(u_dir, "ORDINARY_ELBO.csv"))
  valid <- result$applied && fit$i_iter == 499L &&
    fit$annealing_sweeps == 99L && fit$t1_sweeps == 400L &&
    length(fit$ELBO) == 400L && all(is.finite(fit$ELBO)) &&
    isTRUE(fit$elbo_objective_valid) &&
    identical(as.integer(fit$fit_seed), fit_seed) &&
    identical(vapply(result$snapshots, `[[`, integer(1L), "iteration"),
              as.integer(51:499))
  if (!valid) {
    write_lines_new("STAGE7J_U_FIT_CONTRACT_FAILED",
                    file.path(u_dir, "FIT_CONTRACT_FAILED.txt"))
    stop("Full U continuation violated the fixed sweep/output contract.")
  }
  write_lines_new(c(
    "STAGE7J_U_FIT_COMPLETE", "source_iteration=50",
    "resumed_sweeps=449", "ordinary_sweeps=400",
    "endpoint_parity_not_yet_verified=TRUE",
    "truth_opened=FALSE", "formal_candidate_pool=FALSE",
    paste0("final_ordinary_elbo=",
           format(tail(fit$ELBO, 1L), digits = 16L)),
    paste0("elapsed_seconds=", format(result$elapsed, digits = 10L))
  ), file.path(u_dir, "FIT_COMPLETE.txt"))
  cat("STAGE7J_U_FIT_COMPLETE\n")
}

if (action == "verify") {
  check_prepare_binding()
  load_stage7i_package()
  u_dir <- file.path(output_root, "u_full_v1")
  assert(file.exists(file.path(u_dir, "FIT_COMPLETE.txt")),
         "Full U completion marker missing.")
  parity_dir <- file.path(output_root, "u_parity_v1")
  require_new_dir(parity_dir)
  reference <- reference_snapshots()
  current <- readRDS(file.path(u_dir, "RAW_STATE_SNAPSHOTS_51_499.rds"))
  snapshots <- compare_snapshots(current, reference, 51:499)
  write_csv_new(snapshots, file.path(parity_dir, "ALL_STATE_PARITY.csv"))

  original_fit <- readRDS(stage7g_fit_path)
  resumed_fit <- readRDS(file.path(u_dir, "fit.rds"))
  assert(length(original_fit$ELBO) == 400L &&
           length(resumed_fit$ELBO) == 400L,
         "Ordinary ELBO path must have exactly 400 sweeps.")
  elbo <- data.frame(
    t1_sweep = as.integer(1:400),
    stage7g_elbo = as.numeric(original_fit$ELBO),
    u_elbo = as.numeric(resumed_fit$ELBO),
    abs_difference = abs(as.numeric(resumed_fit$ELBO) -
                         as.numeric(original_fit$ELBO))
  )
  elbo$pass <- is.finite(elbo$abs_difference) &
    elbo$abs_difference <= 1e-6
  write_csv_new(elbo, file.path(parity_dir, "ORDINARY_ELBO_PARITY.csv"))

  endpoint_fields <- c(
    "mu_q_gamma_a", "mu_q_gamma_b", "mu_q_normal_a",
    "mu_q_normal_b", "Sigma_q_normal_a", "Sigma_q_normal_b",
    "mu_q_a_original", "mu_q_b_specific_original",
    "mu_q_zeta", "mu_q_xi", "mu_q_nu_phi", "mu_q_nu_psi",
    "Sigma_q_zeta", "Sigma_q_xi", "Sigma_q_nu_phi",
    "Sigma_q_nu_psi", "mu_q_recip_sigsq_eps"
  )
  endpoint <- do.call(rbind, lapply(endpoint_fields, function(field) {
    value <- compare_value(resumed_fit[[field]], original_fit[[field]],
                           field %in% c("mu_q_gamma_a", "mu_q_gamma_b"))
    data.frame(field = field, max_abs_difference = value$max_abs,
               max_tolerance_excess = value$max_excess,
               ppi_zero_pattern_equal = value$zero_pattern_equal,
               pass = value$pass)
  }))
  write_csv_new(endpoint, file.path(parity_dir, "ENDPOINT_PARITY.csv"))
  configuration <- data.frame(
    check = c("fit_seed", "sweeps", "fit_control",
              "initialization_control", "pre_score_sweeps",
              "function_initialization", "final_objective_role"),
    pass = c(
      identical(resumed_fit$fit_seed, original_fit$fit_seed),
      identical(c(resumed_fit$i_iter, resumed_fit$annealing_sweeps,
                  resumed_fit$t1_sweeps),
                c(original_fit$i_iter, original_fit$annealing_sweeps,
                  original_fit$t1_sweeps)),
      identical(resumed_fit$fit_control, original_fit$fit_control),
      identical(resumed_fit$initialization_control,
                original_fit$initialization_control),
      identical(resumed_fit$pre_score_interface$pre_score_sweeps,
                original_fit$pre_score_interface$pre_score_sweeps),
      identical(resumed_fit$pre_score_interface$function_initialization,
                original_fit$pre_score_interface$function_initialization),
      isTRUE(resumed_fit$elbo_objective_valid) &&
        isTRUE(original_fit$elbo_objective_valid)
    )
  )
  write_csv_new(configuration,
                file.path(parity_dir, "CONFIGURATION_PARITY.csv"))

  original_trace <- original_fit$driver_trace
  original_trace <- original_trace[
    original_trace$iteration >= 51L, , drop = FALSE]
  resumed_trace <- resumed_fit$driver_trace
  resumed_trace <- resumed_trace[
    resumed_trace$iteration %in% unique(original_trace$iteration),
    , drop = FALSE]
  original_trace$sequence <- NULL
  resumed_trace$sequence <- NULL
  rownames(original_trace) <- NULL
  rownames(resumed_trace) <- NULL
  trace_keys <- identical(dim(original_trace), dim(resumed_trace)) &&
    identical(names(original_trace), names(resumed_trace)) &&
    identical(resumed_trace$iteration, original_trace$iteration)
  trace_result <- if (trace_keys) {
    compare_value(resumed_trace, original_trace)
  } else list(pass = FALSE, max_abs = Inf, max_excess = Inf,
              zero_pattern_equal = FALSE)
  trace <- data.frame(rows = nrow(resumed_trace),
                      keys_equal = trace_keys,
                      max_abs_difference = trace_result$max_abs,
                      max_tolerance_excess = trace_result$max_excess,
                      pass = trace_result$pass)
  write_csv_new(trace, file.path(parity_dir, "DRIVER_TRACE_PARITY.csv"))
  pass <- all(snapshots$pass) && all(elbo$pass) &&
    all(endpoint$pass) && all(configuration$pass) && trace$pass[[1L]]
  if (!pass) {
    first_bad <- snapshots[!snapshots$pass, , drop = FALSE]
    write_lines_new(c(
      "STAGE7J_U_PARITY_FAILED", "A_and_S_forbidden=TRUE",
      "no_automatic_rerun=TRUE",
      paste0("first_bad_state_sweep=",
             if (nrow(first_bad)) first_bad$iteration[[1L]] else "none"),
      paste0("first_bad_state_field=",
             if (nrow(first_bad)) first_bad$field[[1L]] else "none")
    ), file.path(parity_dir, "PARITY_FAILED.txt"))
    stop("U parity failed; A/S must not start.")
  }
  write_lines_new(c(
    "STAGE7J_U_PARITY_PASS", "snapshots_compared=449",
    paste0("state_field_checks=", nrow(snapshots)),
    paste0("max_state_abs_difference=",
           format(max(snapshots$max_abs_difference), digits = 16L)),
    paste0("max_elbo_abs_difference=",
           format(max(elbo$abs_difference), digits = 16L)),
    "state_tolerance=1e-8+1e-10*abs(reference)",
    "elbo_absolute_tolerance=1e-6", "ppi_zero_pattern_exact=TRUE",
    "truth_opened=FALSE", "formal_candidate_pool=FALSE"
  ), file.path(parity_dir, "PARITY_PASS.txt"))
  cat("STAGE7J_U_PARITY_PASS\n")
}

if (action == "status") {
  cat("prepare=", file.exists(file.path(output_root, "PREPARE_COMPLETE.txt")),
      " u_complete=", file.exists(file.path(output_root, "u_full_v1",
                                         "FIT_COMPLETE.txt")),
      " parity_pass=", file.exists(file.path(output_root, "u_parity_v1",
                                          "PARITY_PASS.txt")),
      " parity_failed=", file.exists(file.path(output_root, "u_parity_v1",
                                            "PARITY_FAILED.txt")), "\n", sep = "")
}
