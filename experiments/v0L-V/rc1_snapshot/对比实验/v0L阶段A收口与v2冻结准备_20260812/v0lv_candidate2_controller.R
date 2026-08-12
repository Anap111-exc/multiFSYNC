# Candidate2 end-to-end controllers. Formal calls require a separate frozen
# binding and explicit authorization. Candidate-stage tests use mode="smoke".

V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID <-
  "V0LV_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812"
V0LV_CANDIDATE2_MANIFEST_SCHEMA_ID <-
  "V0LV_MANIFEST_SCHEMA_V3_CANDIDATE2_20260812"
V0LV_CANDIDATE2_FACTOR_PPI_RULE <- ">=0.5"

candidate2_generator_config <- function(mode = c("formal", "smoke")) {
  mode <- match.arg(mode)
  if (mode == "formal") {
    return(list(
      S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 5L, n_obs = 8L,
      common_grid = FALSE, n_obs_range = c(6L, 9L),
      sigma_eps = 0.3, bool_sparse_loadings = TRUE,
      prop_sparse = 0.9, score_var_decay = TRUE,
      n_dense = 201L, mean_amp = 0.6, beta_amp = 0,
      bs_degree = c(2L, 3L), identified_loadings = TRUE,
      sparsity_mode = "fixed",
      loading_calibration = list(
        regime = "constant_per_active_loading", anchor_p = 100L,
        active_fraction = 0.1, anchor_norm = 3
      ),
      specific_direction_abs_cosine_max = 0.8,
      formal_experiment = TRUE
    ))
  }
  list(
    S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = 3L, n_obs = 6L,
    common_grid = FALSE, n_obs_range = c(4L, 6L),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE,
    prop_sparse = 0.8, score_var_decay = TRUE,
    n_dense = 51L, mean_amp = 0.6, beta_amp = 0,
    bs_degree = c(2L, 3L), identified_loadings = TRUE,
    sparsity_mode = "fixed",
    loading_calibration = NULL,
    specific_direction_abs_cosine_max = 0.95,
    formal_experiment = FALSE
  )
}

candidate2_trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2L || any(diff(grid) <= 0)) stop("Grid is invalid.")
  delta <- diff(grid)
  c(delta[[1L]] / 2, (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}

candidate2_abs_cosine <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 1e-14) return(NA_real_)
  abs(sum(x * y) / denominator)
}

candidate2_forbidden_observation_names <- function(value) {
  found <- character()
  walk <- function(object) {
    object_names <- names(object)
    if (length(object_names)) {
      bad <- grepl(
        "true|truth|gamma|support|structure|nrmse|ise|recall|purity|match",
        object_names, ignore.case = TRUE
      )
      found <<- c(found, object_names[bad])
    }
    if (is.list(object)) invisible(lapply(object, walk))
    invisible(NULL)
  }
  walk(value)
  unique(found)
}

candidate2_truth_assertions <- function(data, config) {
  truth <- data$true_params
  checks <- list()
  add <- function(id, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      assertion_id = id, passed = isTRUE(passed), detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  add("dimensions_S_n_p_d",
      identical(as.integer(c(truth$S, truth$n_s, truth$p, truth$d_use)),
                as.integer(c(config$S, config$n_s, config$p, config$d))),
      "S,n_s,p,d equal generator contract")
  add("factor_counts_Q_L",
      identical(as.integer(truth$L_f), as.integer(config$L_f)) &&
        identical(as.integer(truth$L_s_by_study), as.integer(config$L_s)),
      "shared/study-specific factor counts equal 1/1/1")
  add("fpca_counts_M",
      identical(as.integer(truth$M_f), as.integer(config$M_f)) &&
        identical(lapply(truth$M_s, as.integer),
                  lapply(config$M_s, as.integer)),
      "factor-internal FPCA truncations match")
  add("spline_K", identical(as.integer(truth$K), as.integer(config$K)),
      "O'Sullivan nonlinear block K matches")
  observed_counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  expected_range <- if (config$formal_experiment) c(6L, 9L) else
    config$n_obs_range
  add("observation_counts", all(observed_counts >= expected_range[[1L]] &
                                  observed_counts <= expected_range[[2L]]),
      paste0("all subject counts in ", paste(expected_range, collapse = "-")))
  support_a <- which(truth$a_true[, 1L] != 0)
  support_b <- lapply(truth$b_true, function(value) which(value[, 1L] != 0))
  expected_active <- as.integer(round((1 - config$prop_sparse) * config$p))
  add("shared_active_count", length(support_a) == expected_active,
      paste0("shared active=", length(support_a)))
  add("specific_active_counts",
      all(lengths(support_b) == expected_active),
      paste0("specific active=", paste(lengths(support_b), collapse = ",")))
  add("shared_specific_support_disjoint",
      all(vapply(support_b, function(value) {
        !length(intersect(support_a, value))
      }, logical(1))),
      "supp(A) is disjoint from each supp(B_s)")
  specific_overlap <- length(intersect(support_b[[1L]], support_b[[2L]]))
  specific_cosine <- candidate2_abs_cosine(
    truth$b_true[[1L]][, 1L], truth$b_true[[2L]][, 1L]
  )
  add("specific_overlap_allowed", TRUE,
      paste0("B1/B2 support overlap allowed; observed overlap=", specific_overlap))
  add("specific_direction_cosine_bound",
      is.finite(specific_cosine) &&
        specific_cosine <= config$specific_direction_abs_cosine_max,
      paste0("abs cosine=", format(specific_cosine, digits = 8),
             "; upper=", config$specific_direction_abs_cosine_max))
  target_norm <- if (config$formal_experiment) 3 * sqrt(5) else NA_real_
  loading_norms <- c(
    sqrt(sum(truth$a_true[, 1L]^2)),
    vapply(truth$b_true, function(value) sqrt(sum(value[, 1L]^2)), numeric(1))
  )
  add("loading_norms", if (config$formal_experiment) {
    max(abs(loading_norms - target_norm)) <= 1e-10
  } else all(is.finite(loading_norms) & loading_norms > 0),
  paste0("norms=", paste(format(loading_norms, digits = 10), collapse = ",")))
  gram_ok <- all(vapply(seq_len(truth$S), function(study) {
    block <- cbind(truth$a_true, truth$b_true[[study]])
    gram <- crossprod(block)
    max(abs(gram - diag(diag(gram)))) <= 1e-10
  }, logical(1)))
  add("loading_gram", gram_ok, "[A,B_s]'[A,B_s] is diagonal")
  weights <- candidate2_trapezoid_weights(truth$t_grid_dense)
  function_blocks <- c(
    truth$phi_dense_list,
    unlist(truth$psi_dense_list, recursive = FALSE)
  )
  function_gram_error <- vapply(function_blocks, function(value) {
    gram <- crossprod(as.matrix(value), weights * as.matrix(value))
    max(abs(gram - diag(ncol(value))))
  }, numeric(1))
  add("eigenfunction_discrete_gram",
      all(function_gram_error <= 1e-8),
      paste0("max error=", format(max(function_gram_error), digits = 8)))
  eigenvalues <- c(
    unlist(truth$lambda_phi_true, use.names = FALSE),
    unlist(truth$lambda_psi_true, recursive = TRUE, use.names = FALSE)
  )
  expected_eigenvalues <- rep(c(0.8, 0.2), 3L)
  add("fpca_eigenvalues", length(eigenvalues) == 6L &&
        max(abs(eigenvalues - expected_eigenvalues)) <= 1e-12,
      paste0("lambda=", paste(eigenvalues, collapse = ",")))
  variance_checks <- list(); variance_index <- 0L
  for (study in seq_len(truth$S)) {
    for (factor in seq_len(truth$L_f)) {
      variance_index <- variance_index + 1L
      variance_checks[[variance_index]] <- list(
        computational = apply(truth$zeta_true[[study]][[factor]], 2L, var),
        fpca = apply(truth$zeta_fpca_true[[study]][[factor]], 2L, var),
        lambda = truth$lambda_phi_true[[factor]]
      )
    }
    for (factor in seq_len(truth$L_s_by_study[[study]])) {
      variance_index <- variance_index + 1L
      variance_checks[[variance_index]] <- list(
        computational = apply(truth$xi_true[[study]][[factor]], 2L, var),
        fpca = apply(truth$xi_fpca_true[[study]][[factor]], 2L, var),
        lambda = truth$lambda_psi_true[[study]][[factor]]
      )
    }
  }
  variance_error <- vapply(variance_checks, function(value) {
    max(abs(value$fpca - value$lambda * value$computational))
  }, numeric(1))
  realized_variances <- unlist(lapply(
    variance_checks, function(value) value$fpca
  ), use.names = FALSE)
  add("fpca_component_variances",
      all(is.finite(realized_variances) & realized_variances > 0) &&
        max(variance_error) <= 1e-12,
      paste0("positive realized FPCA variances; max scaling error=",
             format(max(variance_error), digits = 8)))
  add("measurement_error_sd",
      max(abs(sqrt(truth$sigma2_eps_true) - config$sigma_eps)) <= 1e-12,
      paste0("sigma_eps=", config$sigma_eps))
  mean_error <- 0
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      time <- data$time_obs[[study]][[subject]]
      expected_mean <- vapply(seq_len(truth$p), function(variable) {
        config$mean_amp * (
          sin(2 * pi * time + truth$phase_mu_mat[study, variable]) +
            0.5 * sin(4 * pi * time +
                        2 * truth$phase_mu_mat[study, variable])
        )
      }, numeric(length(time)))
      mean_error <- max(mean_error, abs(
        expected_mean - truth$mu_true_values[[study]][[subject]]
      ))
    }
  }
  add("mean_amplitude", mean_error <= 1e-12,
      paste0("mean_amp=", config$mean_amp, "; max error=", mean_error))
  observation_signal_error <- 0
  observation_y_error <- 0
  dense_scheme_error <- 0
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      shared <- truth$f_true_values[[study]][[subject]] %*% t(truth$a_true)
      specific <- truth$g_true_values[[study]][[subject]] %*%
        t(truth$b_true[[study]])
      signal <- truth$mu_true_values[[study]][[subject]] +
        truth$beta_true_values[[study]][[subject]] + shared + specific
      observation_signal_error <- max(
        observation_signal_error,
        abs(signal - truth$signal_true_values[[study]][[subject]])
      )
      observed <- do.call(cbind, data$Y[[study]][[subject]])
      observation_y_error <- max(
        observation_y_error,
        abs(observed - signal - truth$noise_true_values[[study]][[subject]])
      )
      shared_theta <- truth$theta_dense_list[[1L]] %*%
        truth$eta_true[[study]][[1L]][subject, ]
      shared_phi <- truth$phi_dense_list[[1L]] %*%
        truth$zeta_fpca_true[[study]][[1L]][subject, ]
      specific_kappa <- truth$kappa_dense_list[[study]][[1L]] %*%
        truth$chi_true[[study]][[1L]][subject, ]
      specific_psi <- truth$psi_dense_list[[study]][[1L]] %*%
        truth$xi_fpca_true[[study]][[1L]][subject, ]
      dense_scheme_error <- max(
        dense_scheme_error, abs(shared_theta - shared_phi),
        abs(specific_kappa - specific_psi)
      )
    }
  }
  add("observation_reconstruction_identity",
      max(observation_signal_error, observation_y_error) <= 1e-12,
      paste0("signal error=", observation_signal_error,
             "; Y error=", observation_y_error))
  add("dense_reconstruction_identity", dense_scheme_error <= 1e-12,
      paste0("scheme-1 dense equivalence error=", dense_scheme_error))
  observation_preview <- list(
    Y = data$Y, time_obs = data$time_obs, Z = data$Z,
    dimensions = list(S = config$S, n_s = config$n_s, p = config$p,
                      d = config$d, L_f = config$L_f, L_s = config$L_s,
                      M_f = config$M_f, M_s = config$M_s, K = config$K)
  )
  forbidden <- candidate2_forbidden_observation_names(observation_preview)
  add("observation_bundle_truth_free", !length(forbidden),
      paste0("forbidden names=", paste(forbidden, collapse = ",")))
  result <- do.call(rbind, checks)
  rownames(result) <- NULL
  result
}

candidate2_resolve_generator <- function(
    mode, generator_function = NULL, calibration_function = NULL,
    binding_dir = NULL, project_root = NULL) {
  if (mode == "formal" &&
      (!is.null(generator_function) || !is.null(calibration_function))) {
    stop("Formal generation forbids generator/calibration injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (is.null(generator_function)) {
    if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
      stop("Installed multiFSYNC is required for generation.")
    }
    generator_function <- getExportedValue(
      "multiFSYNC", "simulate_multi_study_structured"
    )
  }
  if (mode == "formal" && is.null(calibration_function)) {
    calibration_env <- new.env(parent = globalenv())
    sys.source(file.path(
      project_root, "对比实验", "v0F维数信息盆地桥接_20260801",
      "v0f_data_calibration.R"
    ), envir = calibration_env)
    calibration_function <- calibration_env$calibrate_loading_information_v0f
  }
  list(generator = generator_function, calibration = calibration_function)
}

candidate2_generate_and_seal <- function(
    data_id, data_seed, output_dir, mode = c("formal", "smoke"),
    generator_function = NULL, calibration_function = NULL,
    binding_dir = NULL, project_root = NULL) {
  mode <- match.arg(mode)
  data_id <- c2_scalar_string(data_id, "data_id")
  data_seed <- c2_integer(data_seed, "data_seed")
  if (dir.exists(output_dir)) stop("Refusing to overwrite data directory.")
  config <- candidate2_generator_config(mode)
  functions <- candidate2_resolve_generator(
    mode, generator_function, calibration_function,
    binding_dir, project_root
  )
  call <- config[setdiff(names(config), c(
    "loading_calibration", "specific_direction_abs_cosine_max",
    "formal_experiment"
  ))]
  call$seed <- data_seed
  captured <- c2_capture_conditions(
    do.call(functions$generator, call), phase = "data_generation"
  )
  if (!is.null(captured$error)) {
    stop("Data generation failed: ", captured$error$message)
  }
  data <- captured$value
  if (!is.null(config$loading_calibration)) {
    if (!is.function(functions$calibration)) {
      stop("Loading calibration function is unavailable.")
    }
    data <- do.call(functions$calibration,
                    c(list(data = data), config$loading_calibration))
  }
  assertions <- candidate2_truth_assertions(data, config)
  if (!all(assertions$passed)) {
    stop("Truth assertions failed: ", paste(
      assertions$assertion_id[!assertions$passed], collapse = ", "
    ))
  }
  observation <- list(
    bundle_class = "v0lv_observation_only_candidate2",
    data_id = data_id, data_seed = data_seed,
    formal_experiment = isTRUE(config$formal_experiment),
    Y = data$Y, time_obs = data$time_obs, Z = data$Z,
    dimensions = config[c("S", "n_s", "p", "d", "L_f", "L_s",
                          "M_f", "M_s", "K")]
  )
  forbidden <- candidate2_forbidden_observation_names(observation)
  if (length(forbidden)) stop("Observation bundle leaks truth fields.")
  sealed <- list(
    bundle_class = "v0lv_sealed_truth_candidate2",
    data_id = data_id, data_seed = data_seed,
    formal_experiment = isTRUE(config$formal_experiment),
    data = data,
    sealed_generation_provenance = list(
      mean_amp = config$mean_amp, generator_config = config,
      generation_warnings = captured$warnings,
      truth_assertions = assertions,
      generated_at_utc = captured$ended_at_utc
    )
  )
  dir.create(output_dir, recursive = TRUE)
  observation_write <- c2_atomic_save_rds(
    observation, file.path(output_dir, "observation_bundle.rds")
  )
  sealed_dir <- file.path(output_dir, "sealed_truth")
  dir.create(sealed_dir)
  truth_write <- c2_atomic_save_rds(
    sealed, file.path(sealed_dir, "truth_bundle.rds")
  )
  c2_atomic_write_csv(assertions, file.path(sealed_dir, "truth_assertions.csv"))
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_DATA_SEAL_COMPLETE",
    paste0("formal_experiment=", toupper(config$formal_experiment)),
    paste0("observation_sha256=", observation_write$sha256),
    paste0("sealed_truth_sha256=", truth_write$sha256),
    "truth_unseal_authorized=FALSE"
  ), file.path(output_dir, "DATA_SEAL_COMPLETE.txt"))
  list(observation = observation, sealed_path = truth_write$path,
       assertions = assertions)
}

candidate2_method_spec <- function(
    method = c("B0", "A"), mode = c("formal", "smoke")) {
  method <- match.arg(method); mode <- match.arg(mode)
  if (mode == "formal") {
    anneal <- if (method == "B0") NULL else c(1, 1.9, 100)
    annealing <- if (method == "B0") 0L else 99L
    return(list(
      route = method, anneal = anneal, annealing_sweeps = annealing,
      maxit = 200L + annealing, t1_budget = 200L, n_g = 51L,
      convergence_rule = "practical", tol_abs = 1e-3, tol_rel = 1e-5,
      practical_control = r_route_v3_candidate2_practical_control(),
      formal_experiment = TRUE
    ))
  }
  anneal <- if (method == "B0") NULL else c(1, 1.2, 3)
  annealing <- if (method == "B0") 0L else 2L
  list(
    route = paste0(method, "_SMOKE_NOT_FORMAL"), anneal = anneal,
    annealing_sweeps = annealing, maxit = 4L + annealing,
    t1_budget = 4L, n_g = 31L, convergence_rule = "parameters",
    tol_abs = 0, tol_rel = 0, practical_control = NULL,
    formal_experiment = FALSE
  )
}

candidate2_resolve_multi_fit <- function(
    mode, fit_function = NULL, binding_dir = NULL, project_root = NULL) {
  if (mode == "formal" && !is.null(fit_function)) {
    stop("Formal B0/A execution forbids fit_function injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (is.null(fit_function)) {
    if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
      stop("Installed multiFSYNC is required.")
    }
    fit_function <- getExportedValue("multiFSYNC", "bayesSYNC_multi")
  }
  if (!is.function(fit_function)) stop("fit_function must be a function.")
  if (mode == "formal") {
    index <- utils::read.csv(file.path(
      binding_dir, "RUNTIME_FUNCTION_BINDINGS.csv"
    ), stringsAsFactors = FALSE)
    row <- index[index$function_name == "multiFSYNC::bayesSYNC_multi",
                 , drop = FALSE]
    if (nrow(row) != 1L ||
        c2_function_signature_sha256(fit_function) != row$sha256[[1L]]) {
      stop("Public bayesSYNC_multi runtime hash mismatch.")
    }
  }
  fit_function
}

candidate2_run_multi_main <- function(
    observation, method = c("B0", "A"), fit_id, fit_seed, seed_index,
    initialization_independence_id, mode = c("formal", "smoke"),
    fit_function = NULL, binding_dir = NULL, project_root = NULL,
    hashes = list(input = NA_character_, config = NA_character_,
                  source = NA_character_)) {
  method <- match.arg(method); mode <- match.arg(mode)
  fit_id <- c2_scalar_string(fit_id, "fit_id")
  fit_seed <- c2_integer(fit_seed, "fit_seed")
  seed_index <- c2_integer(seed_index, "seed_index", 1L)
  initialization_independence_id <- c2_scalar_string(
    initialization_independence_id, "initialization_independence_id"
  )
  if (isTRUE(observation$formal_experiment) != (mode == "formal")) {
    stop("Observation/formal execution mode mismatch.")
  }
  spec <- candidate2_method_spec(method, mode)
  fit_function <- candidate2_resolve_multi_fit(
    mode, fit_function, binding_dir, project_root
  )
  dimension <- observation$dimensions
  call <- list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = dimension$L_f, L_s = dimension$L_s,
    M_f = dimension$M_f, M_s = dimension$M_s, K = dimension$K,
    anneal = spec$anneal, list_hyper = NULL,
    n_g = spec$n_g, time_g = NULL,
    tol_abs = spec$tol_abs, tol_rel = spec$tol_rel,
    maxit = spec$maxit, n_cpus = 1L, verbose = FALSE,
    seed = fit_seed, bool_scale = FALSE, bool_var_spec_prob = FALSE,
    d_0 = as.integer(dimension$p),
    convergence_rule = spec$convergence_rule, lambda_orth = 0,
    practical_control = spec$practical_control,
    initialization = "random", initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ), continuation_state = NULL
  )
  captured <- c2_capture_conditions(
    do.call(fit_function, call), phase = paste0(method, "_main_fit")
  )
  fit <- captured$value
  if (!is.null(fit)) captured$warnings <-
    c2_annotate_fit_warnings(captured$warnings, fit)
  if (is.null(captured$error) && !is.null(fit)) {
    strict <- r_route_v3_candidate2_strict_status(fit)
    t1 <- as.integer(fit$t1_sweeps)
    valid <- identical(as.integer(fit$annealing_sweeps),
                       spec$annealing_sweeps) && t1 >= 1L &&
      t1 <= spec$t1_budget &&
      (mode == "smoke" || t1 == spec$t1_budget ||
         strict$strict_practical_converged)
    if (!valid) {
      captured$error <- list(
        class = "candidate2_route_validation_error",
        message = "Returned B0/A fit violates endpoint semantics.",
        call = "candidate2_run_multi_main", trace_summary = "route validation"
      )
      fit <- NULL
    }
  }
  objective <- if (is.null(fit)) list(
    eligible = FALSE, checks = logical(), invalid_reasons = "terminal_error"
  ) else r_route_v3_candidate2_objective_status(fit)
  strict <- if (is.null(fit)) list(
    strict_practical_converged = FALSE,
    first_strict_convergence_sweep = NA_integer_
  ) else r_route_v3_candidate2_strict_status(fit)
  terminal_status <- if (!is.null(captured$error)) "error" else
    if (strict$strict_practical_converged) "strict_practical_converged" else
      "max_budget_reached"
  record <- list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
    runner_version = "v0lv_candidate2_controller_1.0.0",
    fit_id = fit_id, data_id = observation$data_id,
    method = method, route = spec$route,
    fit_seed = fit_seed, seed_index = seed_index,
    initialization_independence_id = initialization_independence_id,
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = terminal_status,
    terminal_reason = if (terminal_status == "strict_practical_converged") {
      "complete_strict_practical_gate_satisfied"
    } else if (terminal_status == "max_budget_reached") {
      "maximum_registered_T1_budget_without_strict_gate"
    } else "fit_or_route_validation_error",
    error = captured$error %||% list(
      class = "", message = "", call = "", trace_summary = ""
    ), warnings = captured$warnings,
    actual_annealing_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$annealing_sweeps),
    actual_T1_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$t1_sweeps),
    first_strict_convergence_sweep = strict$first_strict_convergence_sweep,
    strict_practical_converged = strict$strict_practical_converged,
    convergence_diagnostics = if (is.null(fit)) list() else list(
      convergence_status = fit$convergence_status,
      convergence_reason = fit$convergence_reason,
      practical_final = fit$practical_final,
      practical_diagnostics = fit$practical_diagnostics
    ),
    auxiliary_output_stability = if (is.null(fit)) {
      r_route_v3_candidate2_output_stability(list())
    } else r_route_v3_candidate2_output_stability(fit),
    objective_eligible = objective$eligible,
    objective_checks = objective$checks,
    objective_invalid_reasons = objective$invalid_reasons,
    final_elbo = if (is.null(fit) || !length(fit$ELBO)) NA_real_ else
      as.numeric(tail(fit$ELBO, 1L)),
    hashes = hashes, formal_experiment = mode == "formal",
    truth_used_for_fit_or_selection = FALSE,
    main_endpoint_policy = "strict_practical_stop_else_max_200_T1"
  )
  c2_validate_terminal_record(record)
  list(fit = fit, terminal_record = record)
}

candidate2_load_bayes_reference <- function(project_root) {
  dependency_packages <- c(
    "ellipse", "lattice", "gtools", "magic", "MASS", "matrixcalc",
    "matrixStats", "pracma", "splines"
  )
  invisible(lapply(dependency_packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing bayesSYNC dependency: ", package)
    }
  }))
  environment <- new.env(parent = globalenv())
  # The historical reference uses several dependency exports without namespace
  # qualifiers. Bind those exports inside the isolated reference environment
  # instead of attaching packages to the process-wide search path.
  for (package in dependency_packages) {
    namespace <- asNamespace(package)
    for (name in getNamespaceExports(package)) {
      if (!exists(name, envir = environment, inherits = TRUE)) {
        assign(name, getExportedValue(package, name), envir = environment)
      }
    }
  }
  reference_dir <- file.path(project_root, "Rcode", "bayesSYNC_ref", "R")
  files <- c("bayesSYNC_package.R", "utils.R", "OSullivan_splines.R",
             "set_hyper.R", "bayesSYNC.R")
  for (file in files) {
    sys.source(file.path(reference_dir, file), envir = environment,
               keep.source = TRUE)
  }
  environment
}

candidate2_run_pooled_main <- function(
    observation, fit_id, fit_seed, seed_index = 1L,
    initialization_independence_id = NULL,
    mode = c("formal", "smoke"),
    fit_function = NULL, binding_dir = NULL, project_root = NULL,
    hashes = list(input = NA_character_, config = NA_character_,
                  source = NA_character_)) {
  mode <- match.arg(mode)
  if (mode == "formal" && !is.null(fit_function)) {
    stop("Formal pooled execution forbids fit_function injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (is.null(fit_function)) {
    fit_function <- candidate2_load_bayes_reference(project_root)$bayesSYNC
  }
  if (!is.function(fit_function)) stop("Pooled fit function is unavailable.")
  fit_id <- c2_scalar_string(fit_id, "fit_id")
  fit_seed <- c2_integer(fit_seed, "fit_seed")
  seed_index <- c2_integer(seed_index, "seed_index", 1L)
  if (is.null(initialization_independence_id)) {
    initialization_independence_id <- paste(
      observation$data_id, "pooled_bayesSYNC", fit_seed, sep = "__"
    )
  }
  initialization_independence_id <- c2_scalar_string(
    initialization_independence_id, "initialization_independence_id"
  )
  pooled_Y <- unname(unlist(observation$Y, recursive = FALSE))
  pooled_time <- unname(unlist(observation$time_obs, recursive = FALSE))
  if (length(pooled_Y) != sum(observation$dimensions$n_s) ||
      length(pooled_time) != length(pooled_Y)) {
    stop("Pooled subject concatenation failed.")
  }
  spec <- if (mode == "formal") list(
    anneal = c(1, 1.9, 100), annealing_sweeps = 99L,
    maxit = 299L, n_g = 51L, K = 5L, formal = TRUE
  ) else list(
    anneal = c(1, 1.2, 3), annealing_sweeps = 2L,
    maxit = 6L, n_g = 31L, K = 3L, formal = FALSE
  )
  captured <- c2_capture_conditions(fit_function(
    time_obs = pooled_time, Y = pooled_Y,
    Q = 3L, L = 2L, K = spec$K,
    anneal = spec$anneal, list_hyper = NULL,
    n_g = spec$n_g, time_g = NULL,
    tol_abs = if (mode == "formal") 1e-3 else 0,
    tol_rel = if (mode == "formal") 1e-5 else 0,
    maxit = spec$maxit, n_cpus = 1L, verbose = FALSE,
    seed = fit_seed, bool_scale = FALSE, bool_var_spec_prob = FALSE,
    show_factor_ppi_progress = FALSE
  ), phase = "pooled_bayesSYNC_main_fit")
  fit <- captured$value
  valid <- is.null(captured$error) && is.list(fit) && fit$Q == 3L &&
    fit$L == 2L && length(fit$ELBO_iter) > 0L &&
    all(is.finite(fit$ELBO_iter)) && all(is.finite(fit$B_hat)) &&
    all(is.finite(fit$factor_ppi))
  objective_reasons <- if (valid) character() else
    if (!is.null(captured$error)) "terminal_error" else "pooled_output_invalid"
  record <- list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
    runner_version = "v0lv_candidate2_pooled_controller_1.0.0",
    fit_id = fit_id, data_id = observation$data_id,
    method = "pooled_bayesSYNC", route = "pooled_two_studies_fit_once_Q3_L2",
    fit_seed = fit_seed, seed_index = seed_index,
    initialization_independence_id = initialization_independence_id,
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = if (valid) "max_budget_reached" else "error",
    terminal_reason = if (valid) {
      "pooled_reference_completed_or_reached_registered_budget"
    } else "pooled_fit_or_output_validation_error",
    error = captured$error %||% list(
      class = "", message = "", call = "", trace_summary = ""
    ), warnings = captured$warnings,
    actual_annealing_sweeps = if (valid) spec$annealing_sweeps else NA_integer_,
    actual_T1_sweeps = if (valid) max(0L, as.integer(fit$i_iter) -
                                       spec$annealing_sweeps) else NA_integer_,
    first_strict_convergence_sweep = NA_integer_,
    strict_practical_converged = FALSE,
    convergence_diagnostics = if (valid) list(
      i_iter = fit$i_iter, final_ELBO_iter = tail(fit$ELBO_iter, 1L)
    ) else list(),
    auxiliary_output_stability = list(available = FALSE),
    objective_eligible = valid,
    objective_checks = c(
      finite_within_model_objective = valid,
      cross_model_ELBO_comparison_forbidden = TRUE
    ),
    objective_invalid_reasons = objective_reasons,
    final_elbo = if (valid) as.numeric(tail(fit$ELBO_iter, 1L)) else NA_real_,
    hashes = hashes, formal_experiment = mode == "formal",
    truth_used_for_fit_or_selection = FALSE,
    studies_combined_once = TRUE, study_labels_passed_to_model = FALSE,
    cross_model_ELBO_comparison_forbidden = TRUE
  )
  c2_validate_terminal_record(record)
  list(fit = if (valid) fit else NULL, terminal_record = record)
}

candidate2_run_R_audit_horizon <- function(
    observation, source_fit, source_record, target_cumulative_T1,
    audit_id, output_dir, mode = c("formal", "smoke"),
    fit_function = NULL, state_function = NULL,
    binding_dir = NULL, project_root = NULL,
    hashes = list(input = NA_character_, config = NA_character_,
                  source = NA_character_)) {
  mode <- match.arg(mode)
  target <- c2_integer(target_cumulative_T1, "target_cumulative_T1", 1L)
  audit_id <- c2_scalar_string(audit_id, "audit_id")
  if (mode == "formal" && (!is.null(fit_function) || !is.null(state_function))) {
    stop("Formal audit execution forbids function injection.")
  }
  if (mode == "formal") c2_verify_formal_binding(binding_dir, project_root)
  if (source_record$terminal_status == "error" || is.null(source_fit)) {
    now <- c2_iso_time()
    record <- list(
      terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
      protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
      runner_version = "v0lv_candidate2_audit_controller_1.0.0",
      fit_id = audit_id, data_id = observation$data_id, method = "R_audit",
      route = "fixed_horizon_continuation", fit_seed = source_record$fit_seed,
      seed_index = source_record$seed_index,
      initialization_independence_id =
        source_record$initialization_independence_id,
      started_at_utc = now, ended_at_utc = now, elapsed_seconds = 0,
      peak_memory_bytes = c2_peak_memory_bytes(gc()),
      terminal_status = "audit_unavailable_due_to_main_error",
      terminal_reason = "main_fit_has_no_recoverable_complete_state",
      error = list(class = "main_fit_error", message =
        source_record$error$message %||% "main fit unavailable",
        call = "", trace_summary = ""),
      warnings = c2_empty_warnings(), actual_annealing_sweeps = 0L,
      actual_T1_sweeps = NA_integer_, first_strict_convergence_sweep =
        NA_integer_, strict_practical_converged = FALSE,
      convergence_diagnostics = list(target_cumulative_T1 = target),
      auxiliary_output_stability = list(available = FALSE),
      objective_eligible = FALSE, objective_checks = logical(),
      objective_invalid_reasons = "main_error_no_audit_state",
      final_elbo = NA_real_, hashes = hashes,
      formal_experiment = mode == "formal",
      truth_used_for_fit_or_selection = FALSE
    )
    return(c2_write_terminal_bundle(record, NULL, output_dir))
  }
  source_cumulative <- as.integer(
    source_record$cumulative_T1_sweeps %||% source_record$actual_T1_sweeps
  )
  if (!is.finite(source_cumulative) || source_cumulative > target) {
    stop("Audit source cumulative T1 is invalid for target.")
  }
  if (source_cumulative == target) {
    now <- c2_iso_time()
    record <- source_record
    record$protocol_id <- V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID
    record$runner_version <- "v0lv_candidate2_audit_controller_1.0.0"
    record$fit_id <- audit_id; record$method <- "R_audit"
    record$route <- "fixed_horizon_anchor_reuse"
    record$started_at_utc <- now; record$ended_at_utc <- now
    record$elapsed_seconds <- 0
    record$terminal_status <- "audit_horizon_reached"
    record$terminal_reason <- "verified_main_or_prior_state_reused_at_target"
    record$actual_annealing_sweeps <- 0L
    record$actual_T1_sweeps <- 0L
    record$cumulative_T1_sweeps <- target
    record$strict_practical_converged <- FALSE
    record$first_strict_convergence_sweep <- NA_integer_
    record$hashes <- hashes
    record$formal_experiment <- mode == "formal"
    return(c2_write_terminal_bundle(record, source_fit, output_dir))
  }
  if (is.null(fit_function)) {
    fit_function <- candidate2_resolve_multi_fit(
      mode, NULL, binding_dir, project_root
    )
  }
  if (is.null(state_function)) {
    if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
      stop("multiFSYNC namespace is required for continuation state.")
    }
    state_function <- getFromNamespace(
      ".make_reduced_continuation_state", "multiFSYNC"
    )
  }
  state <- state_function(source_fit, shared_hat_indices = seq_len(source_fit$L_f))
  additional <- target - source_cumulative
  dimension <- observation$dimensions
  captured <- c2_capture_conditions(do.call(fit_function, list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = state$L_f, L_s = state$L_s,
    M_f = state$M_f, M_s = state$M_s, K = state$K,
    anneal = NULL, list_hyper = source_fit$list_hyper,
    n_g = length(state$time_g), time_g = state$time_g,
    tol_abs = 0, tol_rel = 0, maxit = additional,
    n_cpus = 1L, verbose = FALSE, seed = source_record$fit_seed,
    bool_scale = FALSE,
    bool_var_spec_prob = source_fit$bool_var_spec_prob,
    d_0 = source_fit$d_0, convergence_rule = "parameters",
    lambda_orth = 0, practical_control = NULL,
    initialization = "random",
    initialization_control = source_fit$initialization_control,
    continuation_state = state
  )), phase = "fixed_horizon_audit_continuation")
  fit <- captured$value
  if (!is.null(fit)) captured$warnings <-
    c2_annotate_fit_warnings(captured$warnings, fit)
  reached <- is.null(captured$error) && is.list(fit) &&
    as.integer(fit$t1_sweeps) == additional && length(fit$ELBO) == additional
  objective <- if (reached) r_route_v3_candidate2_objective_status(fit) else
    list(eligible = FALSE, checks = logical(),
         invalid_reasons = if (is.null(captured$error))
           "audit_target_not_reached" else "terminal_error")
  record <- list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID,
    runner_version = "v0lv_candidate2_audit_controller_1.0.0",
    fit_id = audit_id, data_id = observation$data_id, method = "R_audit",
    route = "fixed_horizon_continuation", fit_seed = source_record$fit_seed,
    seed_index = source_record$seed_index,
    initialization_independence_id =
      source_record$initialization_independence_id,
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = if (reached) "audit_horizon_reached" else "error",
    terminal_reason = if (reached) "fixed_cumulative_T1_horizon_reached" else
      "audit_continuation_error_or_shortfall",
    error = captured$error %||% if (reached) list(
      class = "", message = "", call = "", trace_summary = ""
    ) else list(class = "audit_target_shortfall",
                message = "fixed horizon not reached", call = "",
                trace_summary = ""),
    warnings = captured$warnings, actual_annealing_sweeps = 0L,
    actual_T1_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$t1_sweeps),
    cumulative_T1_sweeps = if (reached) target else NA_integer_,
    first_strict_convergence_sweep = NA_integer_,
    strict_practical_converged = FALSE,
    convergence_diagnostics = if (is.null(fit)) list() else list(
      fixed_horizon = target, additional_T1 = additional,
      intermediate_strict_gate_recorded_only = TRUE,
      convergence_status = fit$convergence_status
    ),
    auxiliary_output_stability = list(available = FALSE),
    objective_eligible = reached && objective$eligible,
    objective_checks = objective$checks,
    objective_invalid_reasons = objective$invalid_reasons,
    final_elbo = if (reached) as.numeric(tail(fit$ELBO, 1L)) else NA_real_,
    hashes = hashes, formal_experiment = mode == "formal",
    truth_used_for_fit_or_selection = FALSE,
    main_winner_reopened = FALSE, main_fit_overwritten = FALSE
  )
  c2_write_terminal_bundle(record, if (reached) fit else NULL, output_dir)
}

candidate2_freeze_truth_free_selection <- function(
    r_manifest, terminal_records, output_dir,
    expected_starts = 12L, formal_experiment = FALSE) {
  if (dir.exists(output_dir)) stop("Selection directory already exists.")
  if (!is.list(terminal_records) || !length(terminal_records)) {
    stop("terminal_records must be a non-empty list.")
  }
  summaries <- do.call(rbind, lapply(terminal_records, c2_terminal_summary))
  summaries$formal_experiment <- as.logical(summaries$formal_experiment)
  if (any(summaries$formal_experiment != isTRUE(formal_experiment))) {
    stop("Selection execution-mode mismatch.")
  }
  winner <- select_R_route_v3_candidate2_winner(
    r_manifest, summaries, expected_starts,
    formal_experiment = formal_experiment
  )
  basin <- r_route_v3_candidate2_best_basin(summaries)
  dir.create(output_dir, recursive = TRUE)
  summary_write <- c2_atomic_write_csv(
    summaries, file.path(output_dir, "TRUTH_FREE_ENDPOINTS.csv")
  )
  winner_write <- c2_atomic_write_csv(
    winner, file.path(output_dir, "TRUTH_FREE_SELECTION.csv")
  )
  basin_write <- c2_atomic_write_csv(
    basin, file.path(output_dir, "TRUTH_FREE_ELBO_BASINS.csv")
  )
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_TRUTH_FREE_SELECTION_FROZEN",
    paste0("formal_experiment=", toupper(isTRUE(formal_experiment))),
    paste0("endpoints_sha256=", summary_write$sha256),
    paste0("selection_sha256=", winner_write$sha256),
    paste0("basins_sha256=", basin_write$sha256),
    "selection_used_truth=FALSE",
    "winner_rule=maximum_valid_ordinary_T1_ELBO"
  ), file.path(output_dir, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  list(endpoints = summaries, winner = winner, basin = basin,
       selection_sha256 = winner_write$sha256)
}

candidate2_freeze_generic_truth_free_selection <- function(
    manifest, terminal_records, method = c("B0", "A"), output_dir,
    expected_starts = 4L, formal_experiment = FALSE) {
  method <- match.arg(method)
  if (dir.exists(output_dir)) stop("Selection directory already exists.")
  rows <- manifest[manifest$path_id == method, , drop = FALSE]
  if (nrow(rows) != expected_starts || anyDuplicated(rows$fit_id)) {
    stop("Generic selection manifest has the wrong registered start count.")
  }
  summaries <- do.call(rbind, lapply(terminal_records, c2_terminal_summary))
  manifest_key <- paste(
    rows$fit_id, rows$data_id, rows$fit_seed, rows$seed_index,
    rows$initialization_independence_id, sep = "\r"
  )
  terminal_key <- paste(
    summaries$fit_id, summaries$data_id, summaries$fit_seed,
    summaries$seed_index, summaries$initialization_independence_id,
    sep = "\r"
  )
  if (anyDuplicated(terminal_key) || !setequal(manifest_key, terminal_key)) {
    stop("Generic terminal records do not exactly join to the manifest.")
  }
  summaries <- summaries[match(manifest_key, terminal_key), , drop = FALSE]
  if (any(summaries$method != method) ||
      any(summaries$formal_experiment != isTRUE(formal_experiment)) ||
      any(summaries$truth_used_for_fit_or_selection)) {
    stop("Generic truth-free selection contract was violated.")
  }
  eligible <- summaries$objective_eligible &
    summaries$terminal_status %in% c(
      "strict_practical_converged", "max_budget_reached"
    ) & is.finite(summaries$final_elbo)
  if (!any(eligible)) stop("No objective-eligible endpoint for method ", method)
  candidate <- which(eligible)
  ordered <- candidate[order(
    -summaries$final_elbo[candidate], summaries$fit_seed[candidate],
    summaries$seed_index[candidate]
  )]
  winner <- summaries[ordered[[1L]], , drop = FALSE]
  winner$selection_rule <- "maximum_valid_ordinary_T1_ELBO"
  winner$selected_endpoint_unfinished <-
    !winner$strict_practical_converged[[1L]]
  basin <- r_route_v3_candidate2_best_basin(summaries)
  dir.create(output_dir, recursive = TRUE)
  endpoints_write <- c2_atomic_write_csv(
    summaries, file.path(output_dir, "TRUTH_FREE_ENDPOINTS.csv")
  )
  winner_write <- c2_atomic_write_csv(
    winner, file.path(output_dir, "TRUTH_FREE_SELECTION.csv")
  )
  basin_write <- c2_atomic_write_csv(
    basin, file.path(output_dir, "TRUTH_FREE_ELBO_BASINS.csv")
  )
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_TRUTH_FREE_SELECTION_FROZEN",
    paste0("method=", method),
    paste0("formal_experiment=", toupper(isTRUE(formal_experiment))),
    paste0("endpoints_sha256=", endpoints_write$sha256),
    paste0("selection_sha256=", winner_write$sha256),
    paste0("basins_sha256=", basin_write$sha256),
    "selection_used_truth=FALSE",
    "winner_rule=maximum_valid_ordinary_T1_ELBO"
  ), file.path(output_dir, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  list(endpoints = summaries, winner = winner, basin = basin)
}

candidate2_freeze_combined_truth_free_selection <- function(
    component_selection_dirs, pooled_terminal_record, output_dir,
    formal_experiment = FALSE) {
  if (dir.exists(output_dir)) stop("Combined selection directory exists.")
  required_names <- c("B0", "A", "R")
  if (!is.list(component_selection_dirs) ||
      !identical(sort(names(component_selection_dirs)), sort(required_names))) {
    stop("Combined selection requires named B0, A and R directories.")
  }
  selections <- lapply(required_names, function(method) {
    directory <- component_selection_dirs[[method]]
    marker <- file.path(directory, "TRUTH_FREE_SELECTION_FROZEN.txt")
    path <- file.path(directory, "TRUTH_FREE_SELECTION.csv")
    if (!file.exists(marker) || !file.exists(path)) {
      stop("Component selection is incomplete: ", method)
    }
    marker_text <- readLines(marker, warn = FALSE)
    if (!paste0("selection_sha256=", c2_sha256(path)) %in% marker_text) {
      stop("Component selection marker hash mismatch: ", method)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
  })
  pooled <- c2_terminal_summary(pooled_terminal_record)
  if (!pooled$objective_eligible[[1L]] || pooled$terminal_status == "error" ||
      pooled$formal_experiment[[1L]] != isTRUE(formal_experiment) ||
      pooled$truth_used_for_fit_or_selection[[1L]]) {
    stop("Pooled endpoint is not valid for truth-free selection freeze.")
  }
  pooled$selection_rule <- "single_preregistered_endpoint"
  pooled$selected_endpoint_unfinished <- FALSE
  combined <- do.call(rbind, lapply(c(selections, list(pooled)), function(row) {
    row[, Reduce(intersect, lapply(c(selections, list(pooled)), names)),
        drop = FALSE]
  }))
  dir.create(output_dir, recursive = TRUE)
  selection_write <- c2_atomic_write_csv(
    combined, file.path(output_dir, "TRUTH_FREE_SELECTION.csv")
  )
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_ALL_METHOD_TRUTH_FREE_SELECTION_FROZEN",
    paste0("formal_experiment=", toupper(isTRUE(formal_experiment))),
    paste0("selection_sha256=", selection_write$sha256),
    "methods=B0;A;R;pooled_bayesSYNC", "selection_used_truth=FALSE"
  ), file.path(output_dir, "TRUTH_FREE_SELECTION_FROZEN.txt"))
  invisible(combined)
}

candidate2_freeze_audit_horizon_truth_free <- function(
    audit_manifest, terminal_records, endpoint_fits, target_cumulative_T1,
    output_dir, expected_starts = 12L, formal_experiment = FALSE) {
  target <- c2_integer(target_cumulative_T1, "target_cumulative_T1", 1L)
  if (dir.exists(output_dir)) stop("Audit selection directory already exists.")
  if (!is.data.frame(audit_manifest) ||
      !all(c("audit_id", "data_id", "seed_index", "original_fit_seed",
             "initialization_independence_id", "target_cumulative_T1") %in%
           names(audit_manifest))) {
    stop("Audit manifest lacks exact identity fields.")
  }
  manifest <- audit_manifest[
    as.integer(audit_manifest$target_cumulative_T1) == target, , drop = FALSE
  ]
  if (nrow(manifest) != expected_starts || anyDuplicated(manifest$audit_id)) {
    stop("Audit horizon does not contain the registered start count.")
  }
  summaries <- do.call(rbind, lapply(terminal_records, c2_terminal_summary))
  manifest_key <- paste(
    as.character(manifest$audit_id), as.character(manifest$data_id),
    as.character(manifest$original_fit_seed),
    as.character(manifest$seed_index),
    as.character(manifest$initialization_independence_id), sep = "\r"
  )
  terminal_key <- paste(
    as.character(summaries$fit_id), as.character(summaries$data_id),
    as.character(summaries$fit_seed), as.character(summaries$seed_index),
    as.character(summaries$initialization_independence_id), sep = "\r"
  )
  if (anyDuplicated(terminal_key) || !setequal(manifest_key, terminal_key)) {
    stop("Audit terminal records do not exactly join to the manifest.")
  }
  summaries <- summaries[match(manifest_key, terminal_key), , drop = FALSE]
  if (any(summaries$formal_experiment != isTRUE(formal_experiment)) ||
      any(summaries$truth_used_for_fit_or_selection)) {
    stop("Audit selection mode or truth-isolation contract was violated.")
  }
  eligible <- as.logical(summaries$objective_eligible) &
    summaries$terminal_status == "audit_horizon_reached" &
    is.finite(summaries$final_elbo)
  if (!any(eligible)) stop("No eligible endpoint at the fixed audit horizon.")
  candidate <- which(eligible)
  ordered <- candidate[order(
    -summaries$final_elbo[candidate], summaries$fit_seed[candidate],
    summaries$seed_index[candidate]
  )]
  winner <- summaries[ordered[[1L]], , drop = FALSE]
  winner$selection_rule <- "counterfactual_maximum_valid_T1_ELBO"
  winner$replaces_main_winner <- FALSE
  basin <- r_route_v3_candidate2_best_basin(summaries)
  if (!is.list(endpoint_fits) || is.null(names(endpoint_fits)) ||
      anyDuplicated(names(endpoint_fits)) ||
      !setequal(names(endpoint_fits), summaries$fit_id)) {
    stop("Audit endpoint fits must be named exactly by audit_id.")
  }
  winner_fit <- endpoint_fits[[winner$fit_id[[1L]]]]
  alignment <- do.call(rbind, lapply(summaries$fit_id, function(fit_id) {
    cbind(
      data.frame(fit_id = fit_id, winner_fit_id = winner$fit_id[[1L]],
                 target_cumulative_T1 = target, stringsAsFactors = FALSE),
      r_route_v3_candidate2_endpoint_output_alignment(
        endpoint_fits[[fit_id]], winner_fit
      )
    )
  }))
  dir.create(output_dir, recursive = TRUE)
  endpoints_write <- c2_atomic_write_csv(
    summaries, file.path(output_dir, "AUDIT_TRUTH_FREE_ENDPOINTS.csv")
  )
  winner_write <- c2_atomic_write_csv(
    winner, file.path(output_dir, "AUDIT_COUNTERFACTUAL_WINNER.csv")
  )
  basin_write <- c2_atomic_write_csv(
    basin, file.path(output_dir, "AUDIT_ELBO_BASINS.csv")
  )
  alignment_write <- c2_atomic_write_csv(
    alignment, file.path(output_dir, "AUDIT_ENDPOINT_OUTPUT_ALIGNMENT.csv")
  )
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_AUDIT_TRUTH_FREE_HORIZON_FROZEN",
    paste0("target_cumulative_T1=", target),
    paste0("formal_experiment=", toupper(isTRUE(formal_experiment))),
    paste0("endpoints_sha256=", endpoints_write$sha256),
    paste0("counterfactual_winner_sha256=", winner_write$sha256),
    paste0("basins_sha256=", basin_write$sha256),
    paste0("alignment_sha256=", alignment_write$sha256),
    "truth_used=FALSE", "main_winner_replaced=FALSE"
  ), file.path(output_dir, "AUDIT_TRUTH_FREE_HORIZON_FROZEN.txt"))
  list(endpoints = summaries, winner = winner, basin = basin,
       alignment = alignment)
}

candidate2_authorize_unseal <- function(
    selection_dir, sealed_truth_path, authorization_dir,
    formal_experiment = FALSE, operator_authorized = FALSE) {
  marker <- file.path(selection_dir, "TRUTH_FREE_SELECTION_FROZEN.txt")
  selection <- file.path(selection_dir, "TRUTH_FREE_SELECTION.csv")
  if (!file.exists(marker) || !file.exists(selection)) {
    stop("Truth-free selection is not frozen.")
  }
  marker_text <- readLines(marker, warn = FALSE)
  selection_hash <- c2_sha256(selection)
  if (!paste0("selection_sha256=", selection_hash) %in% marker_text) {
    stop("Selection hash does not match its frozen marker.")
  }
  if (isTRUE(formal_experiment) && !isTRUE(operator_authorized)) {
    stop("Formal truth unsealing requires explicit operator authorization.")
  }
  if (dir.exists(authorization_dir)) stop("Authorization directory exists.")
  dir.create(authorization_dir, recursive = TRUE)
  authorization <- list(
    selection_sha256 = selection_hash,
    sealed_truth_sha256 = c2_sha256(sealed_truth_path),
    formal_experiment = isTRUE(formal_experiment),
    operator_authorized = isTRUE(operator_authorized),
    authorized_at_utc = c2_iso_time()
  )
  c2_atomic_save_rds(
    authorization, file.path(authorization_dir, "UNSEAL_AUTHORIZATION.rds")
  )
  c2_atomic_write_lines(c(
    "V0LV_CANDIDATE2_UNSEAL_AUTHORIZED",
    paste0("selection_sha256=", selection_hash),
    paste0("sealed_truth_sha256=", authorization$sealed_truth_sha256),
    paste0("formal_experiment=", toupper(isTRUE(formal_experiment)))
  ), file.path(authorization_dir, "UNSEAL_AUTHORIZED.txt"))
  invisible(authorization)
}

candidate2_unseal_truth <- function(
    sealed_truth_path, selection_dir, authorization_dir) {
  authorization_path <- file.path(
    authorization_dir, "UNSEAL_AUTHORIZATION.rds"
  )
  if (!file.exists(authorization_path)) {
    stop("Truth bundle is sealed until selection freeze and authorization.")
  }
  authorization <- readRDS(authorization_path)
  selection_hash <- c2_sha256(file.path(
    selection_dir, "TRUTH_FREE_SELECTION.csv"
  ))
  truth_hash <- c2_sha256(sealed_truth_path)
  if (!identical(authorization$selection_sha256, selection_hash) ||
      !identical(authorization$sealed_truth_sha256, truth_hash)) {
    stop("Unseal authorization hashes do not match current artifacts.")
  }
  readRDS(sealed_truth_path)
}
