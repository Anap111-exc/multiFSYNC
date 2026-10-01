#!/usr/bin/env Rscript

# Stage 7E post-unseal mechanism diagnostic. Six truth-initialized fits are
# isolated from Stage 7C candidate selection and formal Monte Carlo results.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript run_stage7e_oracle_20260928_v1.R ACTION STAGE7C_ROOT OUTPUT_ROOT [FIT_ID]")
}
action <- args[[1L]]
stage7c_root <- normalizePath(args[[2L]], mustWork = TRUE)
output_root <- args[[3L]]
if (!action %in% c("preflight", "prepare", "smoke", "fit", "supervisor",
                   "evaluate", "status")) stop("Unknown action.")
suppressPackageStartupMessages(library(multiFSYNC))

manifest <- do.call(rbind, lapply(seq_len(3L), function(index) {
  do.call(rbind, lapply(c("exact", "near"), function(mode) {
    data.frame(
      fit_id = sprintf("s7e_%02d_C1_%s", index, mode),
      data_id = sprintf("covrb7c_%02d__C1", index),
      data_index = index, mode = mode,
      template_seed = 99718000L + index,
      fit_seed = 99718100L + index,
      perturb_seed = 99718200L + index,
      stringsAsFactors = FALSE
    )
  }))
}))
stopifnot(nrow(manifest) == 6L, !anyDuplicated(manifest$fit_id))
trace_sweeps <- c(1:20, 25L, 30L, 40L, 50L, 75L, 100L,
                  150L, 200L, 300L, 380L, 400L)

observation_path <- function(data_id) file.path(
  stage7c_root, "data", data_id, "observation_bundle.rds"
)
truth_path <- function(data_id) file.path(
  stage7c_root, "data", data_id, "sealed_truth", "truth_bundle.rds"
)
sha256 <- function(path) {
  result <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(result) != 1L) stop("sha256sum failed: ", path)
  strsplit(result[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
write_csv_new <- function(value, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  write.csv(value, path, row.names = FALSE, na = "")
}
write_lines_new <- function(value, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  writeLines(value, path)
}
save_rds_new <- function(value, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  saveRDS(value, path, version = 3)
}

preflight <- function() {
  expected_tar <- "93b5fb86fc9cb687154502d00a0920979f51aa7f593966eabb462b2d73fcf924"
  tar_path <- "/root/multiFSYNC-stage7c-20260927-v1.tar.gz"
  if (!file.exists(tar_path) || sha256(tar_path) != expected_tar) {
    stop("Stage 7C bound package source archive changed.")
  }
  if (as.character(packageVersion("multiFSYNC")) != "0.3.0.9000" ||
      normalizePath(system.file(package = "multiFSYNC")) !=
      "/root/multiFSYNC-stage7c-lib-20260927-v1/multiFSYNC") {
    stop("The installed package is not the Stage 7C isolated package.")
  }
  if (!file.exists(file.path(stage7c_root,
                             "FITS_AND_SELECTION_COMPLETE.txt"))) {
    stop("Stage 7C completed source is missing.")
  }
  for (data_id in unique(manifest$data_id)) {
    if (!file.exists(observation_path(data_id)) ||
        !file.exists(truth_path(data_id))) {
      stop("Stage 7C observation or sealed truth is missing: ", data_id)
    }
  }
  invisible(TRUE)
}

project_grid_matrix <- function(value, truth_grid, fit_grid, C_g) {
  value <- as.matrix(value)
  target <- vapply(seq_len(ncol(value)), function(index) {
    approx(truth_grid, value[, index], xout = fit_grid,
           rule = 2, ties = mean)$y
  }, numeric(length(fit_grid)))
  if (!is.matrix(target)) target <- matrix(target, ncol = ncol(value))
  coefficients <- qr.solve(C_g, target)
  fitted <- C_g %*% coefficients
  error <- sqrt(sum((fitted - target)^2) /
                  max(sum(target^2), .Machine$double.eps))
  list(coefficients = coefficients, relative_rmse = error)
}

near_loading <- function(target, cosine, seed) {
  target <- as.numeric(target)
  if (cosine == 1) return(target)
  set.seed(as.integer(seed))
  active <- which(abs(target) > 0)
  if (length(active) < 2L) stop("Near-truth loading needs two active values.")
  unit <- target[active] / sqrt(sum(target[active]^2))
  orthogonal <- rnorm(length(active))
  orthogonal <- orthogonal - sum(orthogonal * unit) * unit
  orthogonal <- orthogonal / sqrt(sum(orthogonal^2))
  result <- numeric(length(target))
  result[active] <- sqrt(sum(target[active]^2)) *
    (cosine * unit + sqrt(1 - cosine^2) * orthogonal)
  result
}

set_loading_moments <- function(parameters, shared, specific) {
  update_one <- function(target) {
    probability <- ifelse(abs(target) > 0, 0.99, 0.01)
    normal_mean <- target / probability
    list(
      probability = probability, normal_mean = normal_mean,
      normal_variance = rep(1e-4, length(target)),
      second = probability * (1e-4 + normal_mean^2)
    )
  }
  value <- update_one(shared)
  parameters$mu_q_normal_a[, 1L] <- value$normal_mean
  parameters$Sigma_q_normal_a[, 1L] <- value$normal_variance
  parameters$mu_q_gamma_a[, 1L] <- value$probability
  parameters$mu_q_a[, 1L] <- shared
  parameters$term_a[, 1L] <- value$second
  for (study in seq_along(specific)) {
    value <- update_one(specific[[study]])
    parameters$mu_q_normal_b[[study]][, 1L] <- value$normal_mean
    parameters$Sigma_q_normal_b[[study]][, 1L] <- value$normal_variance
    parameters$mu_q_gamma_b[[study]][, 1L] <- value$probability
    parameters$mu_q_b_specific[[study]][, 1L] <- specific[[study]]
    parameters$term_b_specific[[study]][, 1L] <- value$second
  }
  parameters
}

small_covariance <- function(existing, dimension, variance) {
  lapply(existing, function(value) {
    result <- diag(variance, dimension)
    diagnostic <- attr(value, "solver_diagnostic", exact = TRUE)
    if (!is.null(diagnostic)) attr(result, "solver_diagnostic") <- diagnostic
    result
  })
}

make_oracle_state <- function(template, truth, mode, perturb_seed) {
  maker <- getFromNamespace(".make_reduced_continuation_state", "multiFSYNC")
  state <- maker(template)
  parameters <- state$parameters
  truth <- truth$data$true_params
  fit_grid <- template$time_g
  C_g <- template$C_g
  truth_grid <- truth$t_grid_dense
  mean_projection <- lapply(truth$mean_dense_list, function(value) {
    project_grid_matrix(value, truth_grid, fit_grid, C_g)
  })
  beta_projection <- lapply(seq_len(2L), function(covariate) {
    project_grid_matrix(truth$beta_true_grid[, , covariate],
                        truth_grid, fit_grid, C_g)
  })
  shared_projection <- project_grid_matrix(
    truth$theta_dense_list[[1L]], truth_grid, fit_grid, C_g
  )
  specific_projection <- lapply(seq_len(2L), function(study) {
    project_grid_matrix(truth$kappa_dense_list[[study]][[1L]],
                        truth_grid, fit_grid, C_g)
  })

  parameters$mu_q_nu_mu <- lapply(mean_projection, function(value) {
    lapply(seq_len(truth$p), function(variable) {
      value$coefficients[, variable]
    })
  })
  parameters$mu_q_nu_beta <- lapply(seq_len(truth$p), function(variable) {
    lapply(seq_len(2L), function(covariate) {
      beta_projection[[covariate]]$coefficients[, variable]
    })
  })
  parameters$mu_q_nu_phi[[1L]] <- shared_projection$coefficients
  parameters$mu_q_zeta <- lapply(seq_len(2L), function(study) {
    list(as.matrix(truth$eta_true[[study]][[1L]]))
  })
  for (study in seq_len(2L)) {
    parameters$mu_q_nu_psi[[study]][[1L]] <-
      specific_projection[[study]]$coefficients
    parameters$mu_q_xi[[study]][[1L]] <-
      as.matrix(truth$chi_true[[study]][[1L]])
  }
  K_total <- nrow(parameters$mu_q_nu_phi[[1L]])
  parameters$Sigma_q_nu_phi[[1L]] <- small_covariance(
    parameters$Sigma_q_nu_phi[[1L]], K_total, 1e-3
  )
  for (study in seq_len(2L)) {
    parameters$Sigma_q_zeta[[study]][[1L]] <- lapply(
      parameters$Sigma_q_zeta[[study]][[1L]],
      function(value) diag(0.05, 2L)
    )
    parameters$Sigma_q_nu_psi[[study]][[1L]] <- small_covariance(
      parameters$Sigma_q_nu_psi[[study]][[1L]], K_total, 1e-3
    )
    parameters$Sigma_q_xi[[study]][[1L]] <- lapply(
      parameters$Sigma_q_xi[[study]][[1L]],
      function(value) diag(0.05, 2L)
    )
  }
  cosine <- if (mode == "exact") 1 else 0.995
  shared_loading <- near_loading(truth$a_true[, 1L], cosine,
                                 perturb_seed + 1L)
  specific_loading <- lapply(seq_len(2L), function(study) {
    near_loading(truth$b_true[[study]][, 1L], cosine,
                 perturb_seed + 1L + study)
  })
  parameters <- set_loading_moments(
    parameters, shared_loading, specific_loading
  )
  stopifnot(identical(dim(parameters$mu_q_recip_sigsq_eps),
                      dim(truth$sigma2_eps_true)))
  parameters$mu_q_recip_sigsq_eps <- 1 / truth$sigma2_eps_true
  state$parameters <- parameters
  state$diagnostics$stage7e_oracle_mode <- mode
  state$diagnostics$truth_used_for_initialization <- TRUE
  state$diagnostics$truth_used_for_stopping_or_selection <- FALSE
  validator <- getFromNamespace(".validate_continuation_state", "multiFSYNC")
  L_s <- if (length(state$L_s) == 1L) rep(state$L_s, 2L) else state$L_s
  validator(
    state = state, template_parameters = maker(template)$parameters,
    S = 2L, n_s = state$n_s, p = state$p, d = state$d,
    L_f = state$L_f, L_s = as.integer(L_s), M_f = state$M_f,
    M_s = state$M_s, K = state$K,
    bool_var_spec_prob = state$bool_var_spec_prob,
    response_center = state$response_center,
    response_scale = state$response_scale, time_g = state$time_g
  )
  diagnostics <- data.frame(
    mode = mode, target_loading_cosine = cosine,
    mean_projection_max_relative_rmse = max(vapply(
      mean_projection, `[[`, numeric(1L), "relative_rmse")),
    beta_projection_max_relative_rmse = max(vapply(
      beta_projection, `[[`, numeric(1L), "relative_rmse")),
    shared_projection_relative_rmse =
      shared_projection$relative_rmse,
    specific_1_projection_relative_rmse =
      specific_projection[[1L]]$relative_rmse,
    specific_2_projection_relative_rmse =
      specific_projection[[2L]]$relative_rmse,
    shared_loading_abs_cosine = abs(sum(
      shared_loading * truth$a_true[, 1L])) /
      sqrt(sum(shared_loading^2) * sum(truth$a_true[, 1L]^2)),
    specific_1_loading_abs_cosine = abs(sum(
      specific_loading[[1L]] * truth$b_true[[1L]][, 1L])) /
      sqrt(sum(specific_loading[[1L]]^2) *
             sum(truth$b_true[[1L]][, 1L]^2)),
    specific_2_loading_abs_cosine = abs(sum(
      specific_loading[[2L]] * truth$b_true[[2L]][, 1L])) /
      sqrt(sum(specific_loading[[2L]]^2) *
             sum(truth$b_true[[2L]][, 1L]^2))
  )
  list(state = state, diagnostics = diagnostics)
}

fit_arguments <- function(observation, state, seed, sweeps,
                          smoke = FALSE) {
  control <- if (smoke) NULL else multiFSYNC::g12_stopping_control()
  if (!smoke) control$checkpoints <- c(380L, 400L)
  list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = 9L, anneal = NULL,
    list_hyper = NULL, n_g = 51L, time_g = NULL,
    tol_abs = if (smoke) 0 else 1e-3,
    tol_rel = if (smoke) 0 else 1e-5,
    maxit = as.integer(sweeps), n_cpus = 1L,
    verbose = FALSE, seed = as.integer(seed),
    bool_scale = FALSE, bool_var_spec_prob = FALSE,
    d_0 = 500L,
    convergence_rule = if (smoke) "parameters" else "practical",
    practical_control = control,
    lambda_orth = 0, initialization = "random",
    continuation_state = state,
    pre_score_sweeps = 0L,
    trace_sweeps = if (smoke) 1:2 else trace_sweeps,
    function_initialization = "current_1_over_m"
  )
}

do_prepare <- function() {
  preflight()
  if (dir.exists(output_root) || file.exists(output_root)) {
    stop("Stage 7E oracle output already exists; refusing overwrite.")
  }
  dir.create(output_root, recursive = TRUE)
  dir.create(file.path(output_root, "states"))
  write_csv_new(manifest, file.path(output_root, "FIT_MANIFEST.csv"))
  bindings <- list()
  checks <- list()
  for (index in seq_len(3L)) {
    rows <- manifest[manifest$data_index == index, , drop = FALSE]
    data_id <- rows$data_id[[1L]]
    obs_path <- observation_path(data_id)
    tr_path <- truth_path(data_id)
    bindings[[index]] <- data.frame(
      data_id = data_id, observation_sha256 = sha256(obs_path),
      truth_sha256 = sha256(tr_path)
    )
    observation <- readRDS(obs_path)
    truth <- readRDS(tr_path)
    stopifnot(identical(observation$data_id, data_id),
              identical(truth$data_id, data_id),
              identical(truth$bundle_class,
                        "covariate_routeb_pilot_sealed_truth"),
              length(observation$Y) == 2L,
              observation$dimensions$p == 500L,
              observation$dimensions$d == 2L,
              truth$data$true_params$K == 9L)
    template <- multiFSYNC::bayesSYNC_multi(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 9L, anneal = NULL,
      n_g = 51L, tol_abs = 0, tol_rel = 0,
      maxit = 1L, n_cpus = 1L, verbose = FALSE,
      seed = rows$template_seed[[1L]], bool_scale = FALSE,
      bool_var_spec_prob = FALSE, d_0 = 500L,
      convergence_rule = "parameters", lambda_orth = 0,
      initialization = "random", continuation_state = NULL
    )
    reference_id <- c("covrb7c_01__C1__start06",
                      "covrb7c_02__C1__start11",
                      "covrb7c_03__C1__start06")[[index]]
    reference <- readRDS(file.path(stage7c_root, "fits",
                                   reference_id, "fit.rds"))
    design_difference <- max(abs(template$C_g - reference$C_g))
    stopifnot(design_difference < 1e-10,
              identical(template$time_g, reference$time_g))
    for (row_index in seq_len(nrow(rows))) {
      row <- rows[row_index, , drop = FALSE]
      result <- make_oracle_state(
        template, truth, row$mode[[1L]], row$perturb_seed[[1L]]
      )
      state_path <- file.path(output_root, "states",
                              paste0(row$fit_id[[1L]], ".rds"))
      save_rds_new(result$state, state_path)
      checks[[length(checks) + 1L]] <- cbind(
        fit_id = row$fit_id[[1L]], data_id = data_id,
        design_max_abs_difference = design_difference,
        result$diagnostics
      )
    }
  }
  write_csv_new(do.call(rbind, bindings),
                file.path(output_root, "INPUT_BINDING.csv"))
  write_csv_new(do.call(rbind, checks),
                file.path(output_root, "STATE_CHECKS.csv"))
  write_lines_new(c(
    "STAGE7E_PREPARE_COMPLETE", "states=6", "full_diagnostic_fits=0",
    "truth_used_for_initialization=TRUE",
    "truth_used_for_stopping_or_selection=FALSE",
    paste0("R_version=", getRversion()),
    paste0("package_version=", packageVersion("multiFSYNC"))
  ), file.path(output_root, "PREPARE_COMPLETE.txt"))
  cat("STAGE7E_PREPARE_COMPLETE states=6\n")
}

do_smoke <- function() {
  preflight()
  if (!file.exists(file.path(output_root, "PREPARE_COMPLETE.txt"))) {
    stop("Prepare marker missing.")
  }
  row <- manifest[manifest$fit_id == "s7e_03_C1_exact", , drop = FALSE]
  state <- readRDS(file.path(output_root, "states",
                             paste0(row$fit_id, ".rds")))
  observation <- readRDS(observation_path(row$data_id[[1L]]))
  fit <- do.call(multiFSYNC::bayesSYNC_multi_pre_score,
                 fit_arguments(observation, state, row$fit_seed[[1L]],
                               2L, smoke = TRUE))
  stopifnot(length(fit$ELBO) == 2L, all(is.finite(fit$ELBO)),
            fit$t1_sweeps == 2L, fit$annealing_sweeps == 0L,
            setequal(unique(fit$driver_trace$iteration), 1:2))
  write_lines_new(c("STAGE7E_MICRO_SMOKE_PASS", "ordinary_sweeps=2",
                    "truth_used_for_stopping_or_selection=FALSE"),
                  file.path(output_root, "MICRO_SMOKE_PASS.txt"))
  cat("STAGE7E_MICRO_SMOKE_PASS ordinary_sweeps=2\n")
}

do_fit <- function(fit_id) {
  preflight()
  if (!file.exists(file.path(output_root, "MICRO_SMOKE_PASS.txt"))) {
    stop("Micro smoke marker missing.")
  }
  row <- manifest[manifest$fit_id == fit_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown diagnostic fit ID.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir) || file.exists(fit_dir)) {
    stop("Refusing to overwrite existing diagnostic fit directory: ", fit_id)
  }
  dir.create(fit_dir, recursive = TRUE)
  started <- Sys.time()
  warnings <- character()
  outcome <- tryCatch(
    withCallingHandlers({
      observation <- readRDS(observation_path(row$data_id[[1L]]))
      state <- readRDS(file.path(output_root, "states",
                                 paste0(fit_id, ".rds")))
      do.call(multiFSYNC::bayesSYNC_multi_pre_score,
              fit_arguments(observation, state, row$fit_seed[[1L]], 400L))
    }, warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }), error = function(condition) condition
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  write_lines_new(warnings, file.path(fit_dir, "WARNINGS.txt"))
  if (inherits(outcome, "error")) {
    write_lines_new(c("STAGE7E_FIT_ERROR", paste0("fit_id=", fit_id),
                      paste0("message=", conditionMessage(outcome))),
                    file.path(fit_dir, "FIT_ERROR.txt"))
    return(FALSE)
  }
  fit <- outcome
  valid <- length(fit$ELBO) == 400L &&
    all(is.finite(fit$ELBO)) &&
    identical(as.integer(fit$t1_sweeps), 400L) &&
    identical(as.integer(fit$annealing_sweeps), 0L) &&
    setequal(unique(fit$driver_trace$iteration), trace_sweeps)
  if (!valid) {
    write_lines_new(c("STAGE7E_FIT_ERROR", paste0("fit_id=", fit_id),
                      "message=endpoint_or_trace_validation_failed"),
                    file.path(fit_dir, "FIT_ERROR.txt"))
    return(FALSE)
  }
  save_rds_new(fit, file.path(fit_dir, "fit.rds"))
  write_csv_new(fit$driver_trace,
                file.path(fit_dir, "DRIVER_TRACE.csv"))
  write_csv_new(data.frame(t1_sweep = seq_along(fit$ELBO),
                           ordinary_elbo = as.numeric(fit$ELBO)),
                file.path(fit_dir, "ORDINARY_ELBO.csv"))
  write_lines_new(c(
    "STAGE7E_FIT_COMPLETE", paste0("fit_id=", fit_id),
    paste0("final_ordinary_elbo=", format(tail(fit$ELBO, 1L), digits = 16L)),
    paste0("elapsed_seconds=", format(elapsed, digits = 10L)),
    "ordinary_sweeps=400", "annealing_sweeps=0",
    "truth_used_for_initialization=TRUE",
    "truth_used_for_stopping_or_selection=FALSE",
    "formal_candidate_pool=FALSE"
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("STAGE7E_FIT_COMPLETE fit_id=", fit_id, "\n", sep = "")
  TRUE
}

do_supervisor <- function() {
  preflight()
  if (!file.exists(file.path(output_root, "MICRO_SMOKE_PASS.txt"))) {
    stop("Micro smoke marker missing.")
  }
  if (dir.exists(file.path(output_root, "fits")) ||
      file.exists(file.path(output_root, "FITS_COMPLETE.txt"))) {
    stop("Diagnostic fits have already been started; refusing rerun.")
  }
  dir.create(file.path(output_root, "fits"))
  dir.create(file.path(output_root, "supervisor_logs"))
  write_lines_new(c("STAGE7E_SIX_FITS_STARTED", "registered=6"),
                  file.path(output_root, "FITS_STARTED.txt"))
  worker <- function(fit_id) {
    log_path <- file.path(output_root, "supervisor_logs",
                          paste0(fit_id, ".log"))
    sink(log_path)
    on.exit(sink(), add = TRUE)
    tryCatch(do_fit(fit_id), error = function(condition) {
      cat("WORKER_ERROR ", conditionMessage(condition), "\n")
      FALSE
    })
  }
  success <- unlist(parallel::mclapply(
    manifest$fit_id, worker, mc.cores = 6L,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  ))
  status <- data.frame(fit_id = manifest$fit_id,
                       complete = as.logical(success))
  write_csv_new(status, file.path(output_root, "SUPERVISOR_STATUS.csv"))
  if (length(success) != 6L || !all(success)) {
    write_lines_new(c("STAGE7E_FITS_INCOMPLETE", "no_automatic_rerun=TRUE"),
                    file.path(output_root, "FITS_INCOMPLETE.txt"))
    stop("At least one diagnostic fit failed; retained for diagnosis.")
  }
  write_lines_new(c("STAGE7E_FITS_COMPLETE", "fits=6",
                    "formal_candidate_pool=FALSE"),
                  file.path(output_root, "FITS_COMPLETE.txt"))
  cat("STAGE7E_FITS_COMPLETE fits=6\n")
}

abs_cosine <- function(x, y) {
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  abs(sum(x * y) / denominator)
}

do_evaluate <- function() {
  preflight()
  if (!file.exists(file.path(output_root, "FITS_COMPLETE.txt"))) {
    stop("Six-fit completion marker missing.")
  }
  evaluation_dir <- file.path(output_root, "evaluation_v1")
  if (dir.exists(evaluation_dir)) stop("Refusing evaluation overwrite.")
  dir.create(evaluation_dir)
  winners <- read.csv(file.path(stage7c_root, "truth_free_selection",
                                "TRUTH_FREE_WINNERS.csv"))
  factors <- list()
  endpoints <- list()
  for (row_index in seq_len(nrow(manifest))) {
    row <- manifest[row_index, , drop = FALSE]
    fit_id <- row$fit_id[[1L]]
    fit <- readRDS(file.path(output_root, "fits", fit_id, "fit.rds"))
    truth <- readRDS(truth_path(row$data_id[[1L]]))$data$true_params
    selected <- c(fit$factor_ppi_shared[[1L]],
                  fit$factor_ppi_specific[[1L]][[1L]],
                  fit$factor_ppi_specific[[2L]][[1L]]) >= 0.5
    estimates <- list(fit$mu_q_a_hat[, 1L],
                      fit$mu_q_b_specific_hat[[1L]][, 1L],
                      fit$mu_q_b_specific_hat[[2L]][, 1L])
    raw <- list(fit$mu_q_a[, 1L],
                fit$mu_q_b_specific[[1L]][, 1L],
                fit$mu_q_b_specific[[2L]][, 1L])
    truths <- list(truth$a_true[, 1L], truth$b_true[[1L]][, 1L],
                   truth$b_true[[2L]][, 1L])
    gamma <- list(fit$mu_q_gamma_a[, 1L],
                  fit$mu_q_gamma_b[[1L]][, 1L],
                  fit$mu_q_gamma_b[[2L]][, 1L])
    cosine <- vapply(seq_len(3L), function(index) {
      abs_cosine(estimates[[index]], truths[[index]])
    }, numeric(1L))
    for (index in seq_len(3L)) {
      support <- abs(truths[[index]]) > 0
      called <- gamma[[index]] >= 0.5
      factors[[length(factors) + 1L]] <- data.frame(
        fit_id = fit_id, data_id = row$data_id[[1L]],
        mode = row$mode[[1L]],
        block = c("shared", "specific_1", "specific_2")[[index]],
        factor_selected = selected[[index]],
        raw_loading_abs_cosine = abs_cosine(raw[[index]], truths[[index]]),
        canonical_loading_abs_cosine = cosine[[index]],
        true_support = sum(support),
        estimated_support = sum(called),
        support_true_positive = sum(support & called),
        support_false_positive = sum(!support & called),
        loading_ppi_mean = mean(gamma[[index]])
      )
    }
    winner <- winners[winners$data_id == row$data_id[[1L]], , drop = FALSE]
    stopifnot(nrow(winner) == 1L)
    endpoint_elbo <- tail(fit$ELBO, 1L)
    endpoints[[row_index]] <- data.frame(
      fit_id = fit_id, data_id = row$data_id[[1L]],
      mode = row$mode[[1L]],
      final_ordinary_elbo = endpoint_elbo,
      stage7c_winner_fit_id = winner$fit_id[[1L]],
      stage7c_winner_ordinary_elbo = winner$final_ordinary_elbo[[1L]],
      elbo_minus_stage7c_winner =
        endpoint_elbo - winner$final_ordinary_elbo[[1L]],
      all_three_factors_selected = all(selected),
      all_three_loading_cosine_ge_0_9 =
        all(is.finite(cosine) & cosine >= 0.9),
      diagnostic_three_direction_core =
        all(selected) && all(is.finite(cosine) & cosine >= 0.9),
      minimum_loading_cosine = min(cosine, na.rm = TRUE),
      ordinary_sweeps = fit$t1_sweeps,
      strict_practical_converged =
        isTRUE(fit$practical_convergence$strict_converged)
    )
  }
  write_csv_new(do.call(rbind, factors),
                file.path(evaluation_dir, "ORACLE_FACTOR_DIAGNOSTICS.csv"))
  write_csv_new(do.call(rbind, endpoints),
                file.path(evaluation_dir, "ORACLE_ENDPOINT_COMPARISON.csv"))
  write_lines_new(c("STAGE7E_EVALUATION_COMPLETE", "endpoints=6",
                    "truth_used_for_scientific_diagnosis=TRUE",
                    "formal_candidate_pool=FALSE"),
                  file.path(evaluation_dir, "EVALUATION_COMPLETE.txt"))
  cat("STAGE7E_EVALUATION_COMPLETE endpoints=6\n")
}

if (action == "preflight") {
  preflight()
  cat("STAGE7E_PREFLIGHT_PASS\n")
} else if (action == "prepare") {
  do_prepare()
} else if (action == "smoke") {
  do_smoke()
} else if (action == "fit") {
  if (length(args) != 4L) stop("fit action requires FIT_ID.")
  if (!do_fit(args[[4L]])) stop("Diagnostic fit failed.")
} else if (action == "supervisor") {
  do_supervisor()
} else if (action == "evaluate") {
  do_evaluate()
} else if (action == "status") {
  paths <- c("PREPARE_COMPLETE.txt", "MICRO_SMOKE_PASS.txt",
             "FITS_STARTED.txt", "FITS_COMPLETE.txt",
             "FITS_INCOMPLETE.txt", "evaluation_v1/EVALUATION_COMPLETE.txt")
  print(data.frame(marker = paths,
                   exists = file.exists(file.path(output_root, paths))))
}
