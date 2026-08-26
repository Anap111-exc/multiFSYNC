#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

s4c_verify_environment(load_runtime = TRUE, package_role = "development")
manifest <- uc_read_csv(file.path(UC_ROOT, "DATA_MANIFEST.csv"))
s4c_validate_data_manifest(manifest)
row <- manifest[1L, , drop = FALSE]
config <- s4c_base_config(row, smoke = TRUE)
generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
call <- config
call$seed <- 82726991L
generated <- do.call(generator, call)

spec <- s4c_fit_control()
spec$practical_control$checkpoints <- c(380L, 400L)
fit_function <- getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score")
p <- length(generated$Y[[1L]][[1L]])
captured <- c2_capture_conditions(
  s4c_with_private_traces(
    list(
      include_initial = TRUE, include_annealing = FALSE,
      t1_sweeps = c(380L, 400L), block_annealing_sweeps = integer()
    ),
    do.call(fit_function, list(
      Y = generated$Y, Z = generated$Z, time_obs = generated$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 3L,
      anneal = spec$anneal, list_hyper = NULL,
      n_g = spec$n_g, time_g = NULL,
      tol_abs = spec$tol_abs, tol_rel = spec$tol_rel,
      maxit = spec$maxit, n_cpus = 1L, verbose = FALSE,
      seed = 82726992L, bool_scale = FALSE,
      bool_var_spec_prob = FALSE, d_0 = as.integer(p),
      convergence_rule = spec$convergence_rule,
      lambda_orth = 0, practical_control = spec$practical_control,
      initialization = "random",
      initialization_control = spec$initialization_control,
      continuation_state = NULL, pre_score_sweeps = 1L,
      trace_sweeps = 1L, function_initialization = "gram_unit_energy"
    ))
  ), phase = "stage4C_nonformal_micro_smoke"
)
if (!is.null(captured$error)) stop(captured$error$message, call. = FALSE)
fit <- captured$value$fit
directions <- captured$value$direction_snapshots
checks <- c(
  fit_returned = is.list(fit),
  g12_used = identical(
    fit$pre_score_interface$function_initialization, "gram_unit_energy"
  ),
  one_prescore = fit$pre_score_interface$pre_score_sweeps == 1L,
  jaoua_annealing_99 = fit$annealing_sweeps == 99L,
  single_cpu = fit$n_cpus_used == 1L,
  public_fixed400_control = identical(
    fit$practical_control, spec$practical_control
  ),
  fixed_t1_400 = fit$t1_sweeps == 400L,
  direction_names = identical(names(directions), c("0", "380", "400")),
  direction_finite = all(vapply(directions, function(item) {
    all(is.finite(c(
      item$shared_loading,
      unlist(item$specific_loading, use.names = FALSE),
      item$shared_contribution, item$specific_contribution
    )))
  }, logical(1L))),
  objective_eligible = isTRUE(
    r_route_v3_candidate2_objective_status(fit)$eligible
  ),
  truth_not_evaluated = TRUE,
  not_formal_v0lv_result = TRUE
)
status <- data.frame(
  check_id = names(checks), passed = as.logical(checks),
  stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(UC_ROOT, "NONFORMAL_SMOKE_STATUS.csv"),
                 row.names = FALSE, quote = TRUE)
if (!all(checks)) stop("Stage-4C nonformal micro smoke failed.", call. = FALSE)
writeLines(c(
  "status=PASS", paste0("completed_utc=", uc_iso_time()),
  "checks=12/12", "ordinary_T1_sweeps=400",
  "full_size_fit_started=FALSE",
  "truth_evaluated=FALSE", "formal_v0lv_result=FALSE"
), file.path(UC_ROOT, "NONFORMAL_SMOKE_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE4C_NONFORMAL_SMOKE_PASS checks=12/12 fixed_t1=400 truth=0\n")
