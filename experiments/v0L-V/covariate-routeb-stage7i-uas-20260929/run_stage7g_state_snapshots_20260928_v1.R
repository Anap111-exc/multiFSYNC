#!/usr/bin/env Rscript

# Stage 7G: replay three frozen Stage 7C C1 winners with raw state snapshots.
# No truth bundle is opened here. Fit inputs and update equations are unchanged.
options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !args[[1L]] %in%
    c("selftest", "prepare", "smoke", "supervisor", "verify", "status")) {
  stop("Usage: Rscript run_stage7g_state_snapshots_20260928_v1.R ",
       "{selftest|prepare|smoke|supervisor|verify|status} ",
       "[STAGE7C_ROOT OUTPUT_ROOT]")
}
action <- args[[1L]]
if (action != "selftest" && length(args) != 3L) {
  stop("STAGE7C_ROOT and OUTPUT_ROOT are required.")
}
stage7c_root <- if (action == "selftest") NULL else
  normalizePath(args[[2L]], mustWork = TRUE)
output_root <- if (action == "selftest") NULL else args[[3L]]

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
write_csv_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
write_lines_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  writeLines(value, path, useBytes = TRUE)
}
save_rds_new <- function(value, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  saveRDS(value, path, version = 3L)
}
sha256 <- function(path) {
  assert(file.exists(path), paste0("Missing source: ", path))
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  assert(length(line) == 1L, paste0("sha256sum failed: ", path))
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

manifest <- data.frame(
  fit_id = c("covrb7c_01__C1__start06", "covrb7c_02__C1__start11",
             "covrb7c_03__C1__start06"),
  data_id = sprintf("covrb7c_%02d__C1", 1:3),
  fit_seed = c(99716106L, 99716211L, 99716306L),
  observation_sha256 = c(
    "2467e559b055ea9efdac3e8352df1cebce3bcd6940c2260dad20030af552ddd6",
    "e9bff47094381f3005fcb7d29e56d0854bb436fafe4ebc5cd9d4fecc4bc6ae32",
    "b681fe2a8b827e6cb9e090b9ae62079ef9bb106ec9547242e5a1698e8a9b1197"),
  reference_fit_sha256 = c(
    "54a0997a5a0f61f905334a59cd41761ae46d149b3dc88c50d487d5ef35b87faa",
    "4d3ba871ccaae77da9edefdf04baa4dae5ab83a92f44d8ea038fda4bbc9f08ff",
    "558f8cc509f247faf7d6200e9a0af5bd1b3f3b765c22b2c8d6295bf4d9d202a6"),
  stage7f_fit_sha256 = c(
    "b2166588d932c4432261e026745f1d1d8185dc0b00f741125d6a779de5d7d761",
    "9f73a1a162ec04495d2caf398e4b6f25956bd51bd68b910b28646aebd471d0f2",
    "fd4ea90f1d1c4cecdf0df45c51409de840ca648120493a6ee4301c42a9bf0e93"),
  stringsAsFactors = FALSE
)
trace_sweeps <- sort(unique(as.integer(c(
  1:20, seq(25, 90, 5), 91:110, seq(115, 200, 5),
  250, 300, 350, 380, 400, 450, 499
))))
old_trace_sweeps <- c(1L, 380L, 400L)
snapshot_sweeps <- as.integer(0:499)
stage7f_root <- "/root/multifsync-results/stage7f_20260928_v1"

observation_path <- function(row) file.path(
  stage7c_root, "data", row$data_id[[1L]], "observation_bundle.rds")
reference_fit_path <- function(row) file.path(
  stage7c_root, "fits", row$fit_id[[1L]], "fit.rds")
stage7f_fit_path <- function(row) file.path(
  stage7f_root, "fits", row$fit_id[[1L]], "fit.rds")
replay_fit_dir <- function(row) file.path(
  output_root, "fits", row$fit_id[[1L]])

check_inputs <- function(check_package = TRUE) {
  assert(file.exists(file.path(stage7c_root,
                               "FITS_AND_SELECTION_COMPLETE.txt")),
         "Stage 7C completion marker missing.")
  frozen <- file.path(stage7c_root, "truth_free_selection",
                      "TRUTH_FREE_WINNERS.csv")
  assert(identical(sha256(frozen),
                   "1e1e7791b3e975287ab195c87a801ab0fa24f97d26558a85fccc3b11af70a877"),
         "Frozen winner table hash changed.")
  winners <- utils::read.csv(frozen, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    winner <- winners[winners$fit_id == row$fit_id[[1L]], , drop = FALSE]
    assert(nrow(winner) == 1L &&
             identical(as.character(winner$data_id), row$data_id[[1L]]) &&
             identical(as.character(winner$beta_arm), "C1") &&
             identical(as.integer(winner$fit_seed), row$fit_seed[[1L]]) &&
             isTRUE(winner$objective_eligible[[1L]]),
           paste0("Not the frozen eligible C1 winner: ", row$fit_id[[1L]]))
    assert(identical(sha256(observation_path(row)),
                     row$observation_sha256[[1L]]),
           paste0("Observation hash changed: ", row$data_id[[1L]]))
    assert(identical(sha256(reference_fit_path(row)),
                     row$reference_fit_sha256[[1L]]),
           paste0("Reference fit hash changed: ", row$fit_id[[1L]]))
    assert(identical(sha256(stage7f_fit_path(row)),
                     row$stage7f_fit_sha256[[1L]]),
           paste0("Stage 7F fit hash changed: ", row$fit_id[[1L]]))
  }
  if (check_package) {
    suppressPackageStartupMessages(library(multiFSYNC))
    archive <- "/root/multiFSYNC-stage7g-20260928-v1.tar.gz"
    assert(identical(sha256(archive),
                     "b3e9050f1de7154baff28b7dd986254df0b634bf85d6f50f464d00dbd4c8f17f"),
           "Stage 7G source archive hash changed.")
    assert(as.character(utils::packageVersion("multiFSYNC")) == "0.3.0.9000" &&
             normalizePath(system.file(package = "multiFSYNC")) ==
             "/root/multiFSYNC-stage7g-lib-20260928-v1/multiFSYNC",
           "Loaded package is not the bound Stage 7G build.")
  }
  invisible(winners)
}

fit_arguments <- function(observation, fit_seed, smoke = FALSE) {
  control <- multiFSYNC::g12_stopping_control()
  control$checkpoints <- c(380L, 400L)
  result <- list(
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
  if (smoke) {
    result["anneal"] <- list(NULL)
    result$maxit <- 2L
    result$tol_abs <- 0
    result$tol_rel <- 0
    result$convergence_rule <- "parameters"
    result["practical_control"] <- list(NULL)
    result$trace_sweeps <- 1:2
  }
  result
}

read_observation <- function(row) {
  observation <- readRDS(observation_path(row))
  assert(identical(observation$bundle_class,
                   "covariate_routeb_pilot_observation_only") &&
           isFALSE(observation$truth_available_to_fit) &&
           identical(observation$data_id, row$data_id[[1L]]) &&
           observation$dimensions$p == 500L &&
           observation$dimensions$d == 2L &&
           observation$dimensions$K == 9L,
         "Observation-only input contract failed.")
  observation
}

do_selftest <- function() {
  assert(nrow(manifest) == 3L && !anyDuplicated(manifest$fit_id),
         "Manifest IDs invalid.")
  assert(identical(trace_sweeps, sort(unique(trace_sweeps))) &&
           all(old_trace_sweeps %in% trace_sweeps) &&
           min(trace_sweeps) == 1L && max(trace_sweeps) == 499L &&
           all(91:110 %in% trace_sweeps),
         "Trace schedule invalid.")
  assert(identical(snapshot_sweeps, sort(unique(snapshot_sweeps))) &&
           identical(snapshot_sweeps, as.integer(0:499)) &&
           all(trace_sweeps %in% snapshot_sweeps) &&
           max(snapshot_sweeps) == 499L,
         "Snapshot schedule invalid.")
  assert(max_abs_difference(c(NA_real_, 1), c(NA_real_, 1)) == 0 &&
           is.infinite(max_abs_difference(c(NA_real_, 1), c(1, NA_real_))) &&
           max_abs_difference(list(a = matrix(1:4, 2)),
                              list(a = matrix(1:4, 2))) == 0,
         "Parity comparison selftest failed.")
  cat("STAGE7G_SELFTEST_PASS fits=3 trace_sweeps=", length(trace_sweeps),
      " snapshot_sweeps=", length(snapshot_sweeps),
      "\n", sep = "")
}

do_prepare <- function() {
  check_inputs()
  assert(!dir.exists(output_root) && !file.exists(output_root),
         "Stage 7G output already exists; refusing overwrite.")
  dir.create(output_root, recursive = TRUE)
  write_csv_new(manifest, file.path(output_root, "FIT_MANIFEST.csv"))
  write_csv_new(data.frame(iteration = trace_sweeps),
                file.path(output_root, "TRACE_SCHEDULE.csv"))
  write_csv_new(data.frame(iteration = snapshot_sweeps),
                file.path(output_root, "SNAPSHOT_SCHEDULE.csv"))
  write_lines_new(c(
    "STAGE7G_PREPARE_COMPLETE", "registered_fits=3",
    "full_fit_arguments_identical_to_stage7f=TRUE",
    "only_change=private_raw_state_snapshot_hook",
    "truth_used_for_fit_or_selection=FALSE",
    "formal_candidate_pool=FALSE",
    paste0("R_version=", getRversion()),
    paste0("package_version=", utils::packageVersion("multiFSYNC"))
  ), file.path(output_root, "PREPARE_COMPLETE.txt"))
  cat("STAGE7G_PREPARE_COMPLETE fits=3\n")
}

do_smoke <- function() {
  check_inputs()
  assert(file.exists(file.path(output_root, "PREPARE_COMPLETE.txt")),
         "Prepare marker missing.")
  assert(!file.exists(file.path(output_root, "MICRO_SMOKE_PASS.txt")),
         "Micro smoke already completed.")
  row <- manifest[3L, , drop = FALSE]
  snapshot_env <- new.env(parent = emptyenv())
  snapshot_env$schema <- "stage7g_snapshot_v1"
  snapshot_env$sweeps <- as.integer(0:2)
  snapshot_env$records <- list()
  old_option <- options(multiFSYNC.__stage7g_snapshot_20260928__ = snapshot_env)
  on.exit(options(old_option), add = TRUE)
  fit <- do.call(multiFSYNC::bayesSYNC_multi_pre_score,
                 fit_arguments(read_observation(row), row$fit_seed[[1L]], TRUE))
  assert(length(fit$ELBO) == 2L && all(is.finite(fit$ELBO)) &&
           fit$t1_sweeps == 2L && fit$annealing_sweeps == 0L &&
           setequal(unique(fit$driver_trace$iteration), 1:2) &&
           identical(vapply(snapshot_env$records, `[[`, integer(1L),
                            "iteration"), 0:2) &&
           max_abs_difference(tail(snapshot_env$records, 1L)[[1L]]$mu_q_xi,
                              fit$mu_q_xi) == 0 &&
           max_abs_difference(tail(snapshot_env$records, 1L)[[1L]]$mu_q_gamma_b,
                              fit$mu_q_gamma_b) == 0,
         "Two-sweep snapshot smoke failed.")
  write_lines_new(c("STAGE7G_MICRO_SMOKE_PASS", "ordinary_sweeps=2",
                    "trace_iterations=1,2", "snapshots=0,1,2",
                    "full_diagnostic_fits=0"),
                  file.path(output_root, "MICRO_SMOKE_PASS.txt"))
  cat("STAGE7G_MICRO_SMOKE_PASS ordinary_sweeps=2\n")
}

do_fit <- function(row) {
  fit_dir <- replay_fit_dir(row)
  assert(!dir.exists(fit_dir) && !file.exists(fit_dir),
         paste0("Refusing replay rerun: ", row$fit_id[[1L]]))
  dir.create(fit_dir, recursive = TRUE)
  snapshot_env <- new.env(parent = emptyenv())
  snapshot_env$schema <- "stage7g_snapshot_v1"
  snapshot_env$sweeps <- snapshot_sweeps
  snapshot_env$records <- list()
  old_option <- options(multiFSYNC.__stage7g_snapshot_20260928__ = snapshot_env)
  on.exit(options(old_option), add = TRUE)
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  outcome <- tryCatch(withCallingHandlers(
    do.call(multiFSYNC::bayesSYNC_multi_pre_score,
            fit_arguments(read_observation(row), row$fit_seed[[1L]])),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))
  ), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  write_lines_new(warnings, file.path(fit_dir, "WARNINGS.txt"))
  save_rds_new(snapshot_env$records,
               file.path(fit_dir, "RAW_STATE_SNAPSHOTS.rds"))
  if (inherits(outcome, "error")) {
    write_lines_new(c("STAGE7G_FIT_ERROR", conditionMessage(outcome)),
                    file.path(fit_dir, "FIT_ERROR.txt"))
    return(FALSE)
  }
  fit <- outcome
  valid <- fit$annealing_sweeps == 99L && fit$t1_sweeps == 400L &&
    fit$n_g == 51L && fit$n_cpus_used == 1L &&
    length(fit$ELBO) == 400L && all(is.finite(fit$ELBO)) &&
    isTRUE(fit$elbo_objective_valid) &&
    identical(as.integer(fit$fit_seed), row$fit_seed[[1L]]) &&
    identical(fit$pre_score_interface$function_initialization,
              "gram_unit_energy") &&
    fit$pre_score_interface$pre_score_sweeps == 1L &&
    identical(as.integer(fit$pre_score_interface$trace_sweeps),
              trace_sweeps) &&
    setequal(unique(fit$driver_trace$iteration), trace_sweeps) &&
    identical(vapply(snapshot_env$records, `[[`, integer(1L),
                     "iteration"), snapshot_sweeps) &&
    max_abs_difference(tail(snapshot_env$records, 1L)[[1L]]$mu_q_xi,
                       fit$mu_q_xi) == 0 &&
    max_abs_difference(tail(snapshot_env$records, 1L)[[1L]]$mu_q_gamma_b,
                       fit$mu_q_gamma_b) == 0
  if (!valid) {
    write_lines_new(c("STAGE7G_FIT_ERROR",
                      "endpoint_or_trace_contract_failed"),
                    file.path(fit_dir, "FIT_ERROR.txt"))
    return(FALSE)
  }
  save_rds_new(fit, file.path(fit_dir, "fit.rds"))
  write_csv_new(fit$driver_trace, file.path(fit_dir, "DRIVER_TRACE.csv"))
  write_csv_new(data.frame(t1_sweep = seq_along(fit$ELBO),
                           ordinary_elbo = as.numeric(fit$ELBO)),
                file.path(fit_dir, "ORDINARY_ELBO.csv"))
  write_csv_new(data.frame(
    iteration = vapply(snapshot_env$records, `[[`, integer(1L), "iteration"),
    temperature = vapply(snapshot_env$records, `[[`, numeric(1L), "temperature"),
    annealing = vapply(snapshot_env$records, `[[`, logical(1L), "annealing")
  ), file.path(fit_dir, "RAW_STATE_INDEX.csv"))
  write_lines_new(c(
    "STAGE7G_FIT_COMPLETE", paste0("fit_id=", row$fit_id[[1L]]),
    paste0("final_ordinary_elbo=",
           format(tail(fit$ELBO, 1L), digits = 16L)),
    paste0("elapsed_seconds=", format(elapsed, digits = 10L)),
    "annealing_sweeps=99", "ordinary_sweeps=400",
    paste0("raw_state_snapshots=", length(snapshot_env$records)),
    "truth_used_for_fit_or_selection=FALSE", "formal_candidate_pool=FALSE"
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("STAGE7G_FIT_COMPLETE fit_id=", row$fit_id[[1L]], "\n", sep = "")
  TRUE
}

do_supervisor <- function() {
  check_inputs()
  assert(file.exists(file.path(output_root, "MICRO_SMOKE_PASS.txt")),
         "Micro smoke marker missing.")
  assert(!dir.exists(file.path(output_root, "fits")) &&
           !file.exists(file.path(output_root, "FITS_STARTED.txt")),
         "Diagnostic fits already started; no automatic rerun.")
  dir.create(file.path(output_root, "fits"))
  dir.create(file.path(output_root, "supervisor_logs"))
  write_lines_new(c("STAGE7G_THREE_FITS_STARTED", "registered=3"),
                  file.path(output_root, "FITS_STARTED.txt"))
  worker <- function(i) {
    row <- manifest[i, , drop = FALSE]
    log_path <- file.path(output_root, "supervisor_logs",
                          paste0(row$fit_id[[1L]], ".log"))
    sink(log_path)
    on.exit(sink(), add = TRUE)
    tryCatch(do_fit(row), error = function(e) {
      cat("WORKER_ERROR ", conditionMessage(e), "\n", sep = "")
      FALSE
    })
  }
  success <- unlist(parallel::mclapply(seq_len(nrow(manifest)), worker,
                                      mc.cores = 3L, mc.preschedule = FALSE,
                                      mc.set.seed = FALSE))
  status <- data.frame(fit_id = manifest$fit_id,
                       complete = as.logical(success))
  write_csv_new(status, file.path(output_root, "SUPERVISOR_STATUS.csv"))
  if (length(success) != 3L || !all(success)) {
    write_lines_new(c("STAGE7G_FITS_INCOMPLETE",
                      "no_automatic_rerun=TRUE"),
                    file.path(output_root, "FITS_INCOMPLETE.txt"))
    stop("One or more replay fits failed; retained for diagnosis.")
  }
  write_lines_new(c("STAGE7G_FITS_COMPLETE", "fits=3",
                    "endpoint_parity_not_yet_verified=TRUE"),
                  file.path(output_root, "FITS_COMPLETE.txt"))
  cat("STAGE7G_FITS_COMPLETE fits=3\n")
}

max_abs_difference <- function(x, y) {
  if (is.null(x) || is.null(y)) {
    return(if (is.null(x) && is.null(y)) 0 else Inf)
  }
  if (is.list(x) || is.list(y)) {
    if (!is.list(x) || !is.list(y) || length(x) != length(y) ||
        !identical(names(x), names(y))) return(Inf)
    if (!length(x)) return(0)
    return(max(vapply(seq_along(x), function(i)
      max_abs_difference(x[[i]], y[[i]]), numeric(1L))))
  }
  if (!is.numeric(x) || !is.numeric(y) ||
      !identical(dim(x), dim(y)) || length(x) != length(y) ||
      !identical(is.na(x), is.na(y))) return(Inf)
  if (!length(x) || all(is.na(x))) return(0)
  if (any(!is.finite(x[!is.na(x)])) ||
      any(!is.finite(y[!is.na(y)]))) return(Inf)
  max(abs(as.numeric(x[!is.na(x)]) - as.numeric(y[!is.na(y)])))
}

do_verify <- function() {
  check_inputs()
  assert(file.exists(file.path(output_root, "FITS_COMPLETE.txt")),
         "Three-fit completion marker missing.")
  parity_dir <- file.path(output_root, "endpoint_parity_v1")
  assert(!dir.exists(parity_dir), "Parity output already exists.")
  dir.create(parity_dir)
  fields <- c(
    "ELBO", "mu_q_gamma_a", "mu_q_gamma_b",
    "mu_q_normal_a", "mu_q_normal_b", "Sigma_q_normal_a",
    "Sigma_q_normal_b", "mu_q_a_original", "mu_q_b_specific_original",
    "mu_q_zeta", "mu_q_xi", "mu_q_nu_phi", "mu_q_nu_psi",
    "Sigma_q_zeta", "Sigma_q_xi", "Sigma_q_nu_phi", "Sigma_q_nu_psi",
    "mu_q_recip_sigsq_eps"
  )
  comparison <- list()
  exact_checks <- list()
  old_trace_comparison <- list()
  stage7f_trace_comparison <- list()
  snapshot_final_comparison <- list()
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    reference <- readRDS(reference_fit_path(row))
    stage7f_reference <- readRDS(stage7f_fit_path(row))
    replay <- readRDS(file.path(replay_fit_dir(row), "fit.rds"))
    snapshots <- readRDS(file.path(replay_fit_dir(row),
                                   "RAW_STATE_SNAPSHOTS.rds"))
    fit_id <- row$fit_id[[1L]]
    exact <- identical(as.integer(replay$fit_seed),
                       as.integer(reference$fit_seed)) &&
      identical(as.integer(c(replay$annealing_sweeps, replay$t1_sweeps,
                             replay$i_iter)),
                as.integer(c(reference$annealing_sweeps,
                             reference$t1_sweeps, reference$i_iter))) &&
      identical(replay$fit_control, reference$fit_control) &&
      identical(replay$initialization_control,
                reference$initialization_control) &&
      identical(replay$pre_score_interface$pre_score_sweeps,
                reference$pre_score_interface$pre_score_sweeps) &&
      identical(replay$pre_score_interface$function_initialization,
                reference$pre_score_interface$function_initialization) &&
      identical(replay$pre_score_interface$random_scale_calibration,
                reference$pre_score_interface$random_scale_calibration)
    exact_checks[[i]] <- data.frame(fit_id = fit_id,
                                    configuration_and_sweeps_equal = exact)
    dense_x <- replay$driver_trace
    dense_y <- stage7f_reference$driver_trace
    dense_keys_equal <- identical(dim(dense_x), dim(dense_y)) &&
      identical(names(dense_x), names(dense_y)) &&
      all(vapply(names(dense_x), function(name) {
        if (is.numeric(dense_x[[name]]) && is.numeric(dense_y[[name]])) {
          TRUE
        } else {
          identical(dense_x[[name]], dense_y[[name]])
        }
      }, logical(1L)))
    dense_numeric <- if (dense_keys_equal) names(dense_x)[vapply(
      dense_x, is.numeric, logical(1L))] else character()
    dense_diff <- if (dense_keys_equal && length(dense_numeric)) max(vapply(
      dense_numeric, function(name)
        max_abs_difference(dense_x[[name]], dense_y[[name]]),
      numeric(1L))) else Inf
    stage7f_trace_comparison[[i]] <- data.frame(
      fit_id = fit_id, rows = nrow(dense_x),
      column_keys_equal = dense_keys_equal,
      max_abs_metric_difference = dense_diff,
      pass = dense_keys_equal && is.finite(dense_diff) && dense_diff == 0)
    snapshot_index_ok <- identical(vapply(snapshots, `[[`, integer(1L),
                                          "iteration"), snapshot_sweeps)
    final_snapshot <- tail(snapshots, 1L)[[1L]]
    snapshot_fields <- setdiff(names(final_snapshot),
                               c("iteration", "temperature", "annealing"))
    snapshot_fields <- snapshot_fields[snapshot_fields %in% names(replay)]
    for (field in snapshot_fields) {
      difference <- max_abs_difference(final_snapshot[[field]],
                                       replay[[field]])
      snapshot_final_comparison[[length(snapshot_final_comparison) + 1L]] <-
        data.frame(fit_id = fit_id, field = field,
                   snapshot_index_ok = snapshot_index_ok,
                   max_abs_difference = difference,
                   pass = snapshot_index_ok && is.finite(difference) &&
                     difference == 0)
    }
    for (field in fields) {
      difference <- max_abs_difference(replay[[field]], reference[[field]])
      comparison[[length(comparison) + 1L]] <- data.frame(
        fit_id = fit_id, field = field,
        max_abs_difference = difference, tolerance = 1e-6,
        pass = is.finite(difference) && difference <= 1e-6)
    }
    for (iteration in old_trace_sweeps) {
      x <- replay$driver_trace[
        replay$driver_trace$iteration == iteration, , drop = FALSE]
      y <- reference$driver_trace[
        reference$driver_trace$iteration == iteration, , drop = FALSE]
      metrics <- intersect(names(x), names(y))
      metrics <- setdiff(metrics, c("sequence", "iteration", "phase",
                                    "block", "study", "factor", "component",
                                    "annealing", "gate_active"))
      key <- c("phase", "block", "study", "factor", "component",
               "annealing", "gate_active")
      key_equal <- nrow(x) == nrow(y) &&
        identical(lapply(key, function(name) x[[name]]),
                  lapply(key, function(name) y[[name]]))
      metric_diff <- if (key_equal) {
        max(vapply(metrics, function(name)
          max_abs_difference(x[[name]], y[[name]]), numeric(1L)))
      } else Inf
      old_trace_comparison[[length(old_trace_comparison) + 1L]] <- data.frame(
        fit_id = fit_id, iteration = iteration,
        row_keys_equal = key_equal,
        max_abs_metric_difference = metric_diff,
        pass = key_equal && is.finite(metric_diff) && metric_diff <= 1e-6)
    }
  }
  comparison <- do.call(rbind, comparison)
  exact_checks <- do.call(rbind, exact_checks)
  old_trace_comparison <- do.call(rbind, old_trace_comparison)
  stage7f_trace_comparison <- do.call(rbind, stage7f_trace_comparison)
  snapshot_final_comparison <- do.call(rbind, snapshot_final_comparison)
  write_csv_new(comparison, file.path(parity_dir, "PARAMETER_PARITY.csv"))
  write_csv_new(exact_checks, file.path(parity_dir, "CONFIGURATION_PARITY.csv"))
  write_csv_new(old_trace_comparison,
                file.path(parity_dir, "ORIGINAL_TRACE_PARITY.csv"))
  write_csv_new(stage7f_trace_comparison,
                file.path(parity_dir, "STAGE7F_DENSE_TRACE_PARITY.csv"))
  write_csv_new(snapshot_final_comparison,
                file.path(parity_dir, "SNAPSHOT_FINAL_STATE_PARITY.csv"))
  pass <- all(comparison$pass) &&
    all(exact_checks$configuration_and_sweeps_equal) &&
    all(old_trace_comparison$pass) &&
    all(stage7f_trace_comparison$pass) &&
    all(snapshot_final_comparison$pass)
  if (!pass) {
    write_lines_new(c("STAGE7G_ENDPOINT_PARITY_FAILED",
                      "mechanism_interpretation_stopped=TRUE",
                      "no_automatic_rerun=TRUE"),
                    file.path(parity_dir, "PARITY_FAILED.txt"))
    cat("STAGE7G_ENDPOINT_PARITY_FAILED\n")
    stop("Endpoint parity failed; mechanism interpretation stopped.")
  }
  trace <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, , drop = FALSE]
    x <- readRDS(file.path(replay_fit_dir(row), "fit.rds"))$driver_trace
    cbind(fit_id = row$fit_id[[1L]], data_id = row$data_id[[1L]], x)
  }))
  write_csv_new(trace, file.path(parity_dir, "TRACE_COMBINED.csv"))
  write_lines_new(c("STAGE7G_ENDPOINT_PARITY_PASS",
                    "max_absolute_tolerance=1e-6",
                    "truth_used_for_fit_or_selection=FALSE"),
                  file.path(parity_dir, "PARITY_PASS.txt"))
  cat("STAGE7G_ENDPOINT_PARITY_PASS\n")
  invisible(TRUE)
}

do_status <- function() {
  fit_dir <- file.path(output_root, "fits")
  complete <- vapply(manifest$fit_id, function(id)
    file.exists(file.path(fit_dir, id, "FIT_COMPLETE.txt")), logical(1L))
  failed <- vapply(manifest$fit_id, function(id)
    file.exists(file.path(fit_dir, id, "FIT_ERROR.txt")), logical(1L))
  cat("complete=", sum(complete), " failed=", sum(failed),
      " parity_pass=", file.exists(file.path(output_root,
                                             "endpoint_parity_v1",
                                             "PARITY_PASS.txt")),
      " parity_failed=", file.exists(file.path(output_root,
                                               "endpoint_parity_v1",
                                               "PARITY_FAILED.txt")),
      "\n", sep = "")
}

if (action == "selftest") do_selftest()
if (action == "prepare") do_prepare()
if (action == "smoke") do_smoke()
if (action == "supervisor") do_supervisor()
if (action == "verify") do_verify()
if (action == "status") do_status()
