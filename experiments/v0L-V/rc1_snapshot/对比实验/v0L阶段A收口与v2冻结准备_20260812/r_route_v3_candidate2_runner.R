# Candidate2 runner for the R-only v3 route. The mathematical fit is delegated
# to the frozen public multiFSYNC interface; this file controls identity,
# warnings, endpoint semantics, objective eligibility and truth-free selection.

R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID <-
  "R_ONLY_MULTISTART_V3_CANDIDATE2_20260812"
R_ROUTE_V3_CANDIDATE2_RUNNER_VERSION <-
  "3.0.0-candidate2.1-20260812"
R_ROUTE_V3_CANDIDATE2_ID <- "random__pre1__anneal_default"

r_route_v3_candidate2_practical_control <- function() {
  list(
    min_t1 = 80L, window = 20L, long_window = 60L,
    consecutive = 5L, max_t1 = 200L,
    elbo_abs_rate = 0.03, elbo_per_response_rate = 1e-4,
    fitted_nrmse = 1e-3, rss_rel = 1e-3,
    ppi_max_abs = 1e-2, ppi_gate = "max", ppi_quantile = 0.95,
    ppi_quantile_max_abs = 1e-2, factor_ppi_max_abs = 1e-2,
    long_fitted_nrmse = 3e-3, long_rss_rel = 6e-3,
    long_ppi_max_abs = 1e-2, long_ppi_quantile_max_abs = 1e-2,
    long_factor_ppi_max_abs = 1e-2, checkpoints = integer()
  )
}

r_route_v3_candidate2_spec <- function(mode = c("formal", "smoke")) {
  mode <- match.arg(mode)
  if (mode == "formal") {
    return(list(
      mode = mode, formal_experiment = TRUE,
      route_id = R_ROUTE_V3_CANDIDATE2_ID,
      anneal = c(1, 1.9, 100), annealing_sweeps = 99L,
      t1_budget = 200L, maxit = 299L, n_g = 51L,
      tol_abs = 1e-3, tol_rel = 1e-5,
      convergence_rule = "practical",
      initialization_control = list(
        grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
      ),
      practical_control = r_route_v3_candidate2_practical_control()
    ))
  }
  list(
    mode = mode, formal_experiment = FALSE,
    route_id = "R_V3_CANDIDATE2_SMOKE_NOT_FORMAL",
    anneal = c(1, 1.2, 3), annealing_sweeps = 2L,
    t1_budget = 4L, maxit = 6L, n_g = 31L,
    tol_abs = 0, tol_rel = 0, convergence_rule = "parameters",
    initialization_control = list(
      grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8
    ),
    practical_control = NULL
  )
}

r_route_v3_candidate2_resolve_fit_function <- function(
    mode = c("formal", "smoke"), fit_function = NULL,
    binding_dir = NULL, project_root = NULL) {
  mode <- match.arg(mode)
  if (mode == "formal" && !is.null(fit_function)) {
    stop("Formal R execution forbids fit_function injection.")
  }
  if (mode == "formal") {
    if (is.null(binding_dir) || is.null(project_root)) {
      stop("Formal R execution requires a frozen binding directory.")
    }
    c2_verify_formal_binding(binding_dir, project_root)
  }
  if (is.null(fit_function)) {
    if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
      stop("The installed multiFSYNC namespace is required.")
    }
    fit_function <- getExportedValue(
      "multiFSYNC", "bayesSYNC_multi_pre_score"
    )
  }
  if (!is.function(fit_function)) stop("fit_function must be a function.")
  if (mode == "formal") {
    function_index_path <- file.path(
      binding_dir, "RUNTIME_FUNCTION_BINDINGS.csv"
    )
    if (!file.exists(function_index_path)) {
      stop("Frozen runtime function binding is missing.")
    }
    index <- utils::read.csv(function_index_path, stringsAsFactors = FALSE)
    row <- index[index$function_name ==
                   "multiFSYNC::bayesSYNC_multi_pre_score", , drop = FALSE]
    if (nrow(row) != 1L ||
        c2_function_signature_sha256(fit_function) != row$sha256[[1L]]) {
      stop("Public bayesSYNC_multi_pre_score runtime hash mismatch.")
    }
  }
  fit_function
}

r_route_v3_candidate2_objective_status <- function(fit) {
  finite_elbo <- length(fit$ELBO) > 0L && all(is.finite(fit$ELBO))
  finite_components <- length(fit$ELBO_finite) == length(fit$ELBO) &&
    all(vapply(fit$ELBO_finite, function(value) {
      is.list(value) && isTRUE(value$all)
    }, logical(1)))
  checks <- c(
    annealing_completed = isTRUE(fit$annealing_completed),
    final_temperature_one = isTRUE(all.equal(
      as.numeric(fit$final_temperature), 1
    )),
    ordinary_elbo_present_and_finite = finite_elbo,
    ordinary_elbo_components_finite = finite_components,
    ordinary_objective_valid = isTRUE(fit$elbo_objective_valid),
    no_T1_elbo_decrease = identical(
      as.integer(fit$ELBO_decrease_count), 0L
    ),
    no_T1_solver_jitter = identical(
      as.integer(fit$elbo_t1_jitter_count), 0L
    )
  )
  list(eligible = all(checks), checks = checks,
       invalid_reasons = names(checks)[!checks])
}

r_route_v3_candidate2_strict_status <- function(fit) {
  diagnostics <- fit$practical_diagnostics
  control <- fit$practical_control
  reached <- integer()
  if (is.data.frame(diagnostics) && nrow(diagnostics) && is.list(control) &&
      "streak" %in% names(diagnostics)) {
    reached <- which(as.integer(diagnostics$streak) >=
                       as.integer(control$consecutive))
  }
  strict <- isTRUE(fit$practical_converged) && isTRUE(fit$converged) &&
    identical(as.character(fit$convergence_status), "converged_practical")
  list(
    strict_practical_converged = strict,
    first_strict_convergence_sweep = if (length(reached)) {
      as.integer(diagnostics$t1_sweep[[reached[[1L]]]])
    } else NA_integer_
  )
}

r_route_v3_candidate2_output_stability <- function(fit) {
  diagnostics <- fit$practical_diagnostics
  control <- fit$practical_control
  required <- c(
    "t1_sweep", "objective_pass", "finite", "monotone", "long_monotone",
    "fitted_nrmse", "rss_rel", "ppi_quantile_abs",
    "factor_ppi_max_abs", "long_fitted_nrmse", "long_rss_rel",
    "long_ppi_quantile_abs", "long_factor_ppi_max_abs"
  )
  unavailable <- list(
    available = FALSE, ever_reached_5_consecutive = NA,
    first_reached_sweep = NA_integer_, maximum_streak = NA_integer_,
    endpoint_pass = NA, endpoint_streak = NA_integer_,
    endpoint_stable_5_consecutive = NA
  )
  if (!is.data.frame(diagnostics) || !nrow(diagnostics) ||
      !all(required %in% names(diagnostics)) || !is.list(control)) {
    return(unavailable)
  }
  pass <- vapply(seq_len(nrow(diagnostics)), function(index) {
    row <- diagnostics[index, , drop = FALSE]
    isTRUE(row$objective_pass[[1L]]) && isTRUE(row$finite[[1L]]) &&
      isTRUE(row$monotone[[1L]]) && isTRUE(row$long_monotone[[1L]]) &&
      row$fitted_nrmse[[1L]] <= control$fitted_nrmse &&
      row$rss_rel[[1L]] <= control$rss_rel &&
      row$ppi_quantile_abs[[1L]] <= control$ppi_quantile_max_abs &&
      row$factor_ppi_max_abs[[1L]] <= control$factor_ppi_max_abs &&
      row$long_fitted_nrmse[[1L]] <= control$long_fitted_nrmse &&
      row$long_rss_rel[[1L]] <= control$long_rss_rel &&
      row$long_ppi_quantile_abs[[1L]] <=
        control$long_ppi_quantile_max_abs &&
      row$long_factor_ppi_max_abs[[1L]] <=
        control$long_factor_ppi_max_abs
  }, logical(1))
  streak <- integer(length(pass))
  for (index in seq_along(pass)) {
    streak[[index]] <- if (pass[[index]]) {
      if (index == 1L) 1L else streak[[index - 1L]] + 1L
    } else 0L
  }
  reached <- which(streak >= 5L)
  endpoint_streak <- if (length(streak)) tail(streak, 1L) else 0L
  list(
    available = TRUE,
    ever_reached_5_consecutive = length(reached) > 0L,
    first_reached_sweep = if (length(reached)) {
      as.integer(diagnostics$t1_sweep[[reached[[1L]]]])
    } else NA_integer_,
    maximum_streak = if (length(streak)) max(streak) else 0L,
    endpoint_pass = isTRUE(tail(pass, 1L)),
    endpoint_streak = as.integer(endpoint_streak),
    endpoint_stable_5_consecutive = endpoint_streak >= 5L
  )
}

.r_route_v3_candidate2_relative_l2 <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  sqrt(sum((current - reference)^2)) / (1 + sqrt(sum(reference^2)))
}

.r_route_v3_candidate2_nrmse <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  sqrt(mean((current - reference)^2)) / (1 + sqrt(mean(reference^2)))
}

.r_route_v3_candidate2_max_abs <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  if (!length(current)) 0 else max(abs(current - reference))
}

r_route_v3_candidate2_endpoint_output_alignment <- function(
    current_fit, winner_fit) {
  current <- current_fit$stability_snapshot
  winner <- winner_fit$stability_snapshot
  control <- winner_fit$practical_control
  required <- c("fitted", "rss", "ppi", "factor_ppi")
  if (!is.list(current) || !is.list(winner) || !is.list(control) ||
      !all(required %in% names(current)) ||
      !all(required %in% names(winner))) {
    return(data.frame(
      fitted_nrmse_to_winner = Inf, rss_relative_l2_to_winner = Inf,
      variable_ppi_max_abs_to_winner = Inf,
      factor_ppi_max_abs_to_winner = Inf,
      endpoint_output_aligned = FALSE, stringsAsFactors = FALSE
    ))
  }
  fitted <- .r_route_v3_candidate2_nrmse(current$fitted, winner$fitted)
  rss <- .r_route_v3_candidate2_relative_l2(current$rss, winner$rss)
  variable_ppi <- .r_route_v3_candidate2_max_abs(current$ppi, winner$ppi)
  factor_ppi <- .r_route_v3_candidate2_max_abs(
    current$factor_ppi, winner$factor_ppi
  )
  data.frame(
    fitted_nrmse_to_winner = fitted,
    rss_relative_l2_to_winner = rss,
    variable_ppi_max_abs_to_winner = variable_ppi,
    factor_ppi_max_abs_to_winner = factor_ppi,
    endpoint_output_aligned = all(is.finite(c(
      fitted, rss, variable_ppi, factor_ppi
    ))) && fitted <= 1e-3 && rss <= 1e-3 &&
      variable_ppi <= 1e-2 && factor_ppi <= 1e-2,
    stringsAsFactors = FALSE
  )
}

r_route_v3_candidate2_validate_return <- function(fit, spec, p) {
  driver <- fit$driver_diagnostic_control
  public_interface <- fit$pre_score_interface
  strict <- r_route_v3_candidate2_strict_status(fit)
  t1 <- as.integer(fit$t1_sweeps)
  no_scaling <- is.null(fit$mean_mean_across_subjects) &&
    is.null(fit$sd_mean_across_subjects) && length(fit$response_scale) == p &&
    isTRUE(all.equal(as.numeric(fit$response_scale), rep(1, p)))
  valid <- is.list(driver) && is.list(public_interface) &&
    identical(as.integer(driver$pre_score_sweeps), 1L) &&
    identical(as.integer(driver$dense_gate_sweeps), 0L) &&
    identical(as.character(driver$random_scale_calibration), "none") &&
    identical(as.character(public_interface$interface),
              "bayesSYNC_multi_pre_score") &&
    identical(as.integer(public_interface$pre_score_sweeps), 1L) &&
    identical(as.character(fit$initialization), "random") &&
    identical(as.integer(fit$annealing_sweeps), spec$annealing_sweeps) &&
    identical(as.integer(fit$n_g), spec$n_g) &&
    identical(as.integer(fit$d_0), as.integer(p)) &&
    identical(as.integer(fit$n_cpus_used), 1L) &&
    identical(as.numeric(fit$lambda_orth), 0) && no_scaling &&
    length(t1) == 1L && !is.na(t1) && t1 >= 1L && t1 <= spec$t1_budget
  if (spec$mode == "formal") {
    valid <- valid && identical(as.character(fit$convergence_rule),
                                "practical") &&
      identical(fit$practical_control, spec$practical_control) &&
      (t1 == spec$t1_budget || strict$strict_practical_converged)
  }
  if (!valid) stop("Returned fit violates candidate2 R-route semantics.")
  invisible(TRUE)
}

fit_R_route_v3_candidate2 <- function(
    Y, time_obs, Z = NULL, L_f, L_s, M_f, M_s, K,
    fit_id, data_id, fit_seed, seed_index,
    initialization_independence_id,
    mode = c("formal", "smoke"), fit_function = NULL,
    binding_dir = NULL, project_root = NULL,
    hashes = list(input = NA_character_, config = NA_character_,
                  source = NA_character_), verbose = FALSE) {
  mode <- match.arg(mode)
  fit_id <- c2_scalar_string(fit_id, "fit_id")
  data_id <- c2_scalar_string(data_id, "data_id")
  initialization_independence_id <- c2_scalar_string(
    initialization_independence_id, "initialization_independence_id"
  )
  fit_seed <- c2_integer(fit_seed, "fit_seed")
  seed_index <- c2_integer(seed_index, "seed_index", 1L)
  spec <- r_route_v3_candidate2_spec(mode)
  fit_function <- r_route_v3_candidate2_resolve_fit_function(
    mode, fit_function, binding_dir, project_root
  )
  if (!is.list(Y) || !length(Y) || !is.list(time_obs) || !length(time_obs)) {
    stop("Y and time_obs must be explicit study lists.")
  }
  p <- length(Y[[1L]][[1L]])
  call_args <- list(
    Y = Y, Z = Z, time_obs = time_obs,
    L_f = L_f, L_s = L_s, M_f = M_f, M_s = M_s, K = K,
    anneal = spec$anneal, list_hyper = NULL,
    n_g = spec$n_g, time_g = NULL,
    tol_abs = spec$tol_abs, tol_rel = spec$tol_rel,
    maxit = spec$maxit, n_cpus = 1L, verbose = isTRUE(verbose),
    seed = fit_seed, bool_scale = FALSE, bool_var_spec_prob = FALSE,
    d_0 = as.integer(p), convergence_rule = spec$convergence_rule,
    lambda_orth = 0, practical_control = spec$practical_control,
    initialization = "random",
    initialization_control = spec$initialization_control,
    continuation_state = NULL, pre_score_sweeps = 1L, trace_sweeps = 1L
  )
  captured <- c2_capture_conditions(
    do.call(fit_function, call_args), phase = "main_fit"
  )
  fit <- captured$value
  if (!is.null(fit)) {
    captured$warnings <- c2_annotate_fit_warnings(captured$warnings, fit)
  }
  if (is.null(captured$error)) {
    validation_error <- tryCatch({
      r_route_v3_candidate2_validate_return(fit, spec, p); NULL
    }, error = function(condition) condition)
    if (!is.null(validation_error)) {
      captured$error <- list(
        class = paste(class(validation_error), collapse = ";"),
        message = conditionMessage(validation_error),
        call = c2_condition_call(validation_error),
        trace_summary = "candidate2 route-return validation"
      )
      fit <- NULL
    }
  }
  objective <- if (is.null(fit)) {
    list(eligible = FALSE, checks = logical(),
         invalid_reasons = "terminal_error")
  } else r_route_v3_candidate2_objective_status(fit)
  strict <- if (is.null(fit)) {
    list(strict_practical_converged = FALSE,
         first_strict_convergence_sweep = NA_integer_)
  } else r_route_v3_candidate2_strict_status(fit)
  stability <- if (is.null(fit)) {
    r_route_v3_candidate2_output_stability(list())
  } else r_route_v3_candidate2_output_stability(fit)
  terminal_status <- if (!is.null(captured$error)) {
    "error"
  } else if (strict$strict_practical_converged) {
    "strict_practical_converged"
  } else {
    "max_budget_reached"
  }
  terminal_reason <- switch(
    terminal_status,
    strict_practical_converged = "complete_strict_practical_gate_satisfied",
    max_budget_reached = "maximum_registered_T1_budget_without_strict_gate",
    error = "fit_or_route_validation_error"
  )
  record <- list(
    terminal_schema = V0LV_CANDIDATE2_TERMINAL_SCHEMA,
    protocol_id = R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID,
    runner_version = R_ROUTE_V3_CANDIDATE2_RUNNER_VERSION,
    fit_id = fit_id, data_id = data_id, method = "R",
    route = spec$route_id, fit_seed = fit_seed, seed_index = seed_index,
    initialization_independence_id = initialization_independence_id,
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = terminal_status, terminal_reason = terminal_reason,
    error = captured$error %||% list(
      class = "", message = "", call = "", trace_summary = ""
    ),
    warnings = captured$warnings,
    actual_annealing_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$annealing_sweeps),
    actual_T1_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$t1_sweeps),
    first_strict_convergence_sweep =
      strict$first_strict_convergence_sweep,
    strict_practical_converged = strict$strict_practical_converged,
    convergence_diagnostics = if (is.null(fit)) list() else list(
      convergence_status = fit$convergence_status,
      convergence_reason = fit$convergence_reason,
      practical_final = fit$practical_final,
      practical_diagnostics = fit$practical_diagnostics
    ),
    auxiliary_output_stability = stability,
    objective_eligible = objective$eligible,
    objective_checks = objective$checks,
    objective_invalid_reasons = objective$invalid_reasons,
    final_elbo = if (is.null(fit) || !length(fit$ELBO)) NA_real_ else
      as.numeric(tail(fit$ELBO, 1L)),
    hashes = hashes,
    formal_experiment = spec$formal_experiment,
    truth_used_for_fit_or_selection = FALSE,
    main_endpoint_policy = "strict_practical_stop_else_max_200_T1",
    maximum_T1_sweeps = spec$t1_budget,
    offline_racing_changes_fit = FALSE
  )
  c2_validate_terminal_record(record)
  list(fit = fit, terminal_record = record)
}

validate_R_route_v3_candidate2_manifest <- function(
    manifest, expected_starts = 12L, formal_experiment = TRUE) {
  expected_starts <- c2_integer(expected_starts, "expected_starts", 1L)
  required <- c(
    "fit_id", "data_id", "r_only_protocol_id", "route_id", "fit_seed",
    "seed_index", "initialization_independence_id",
    "terminal_record_required", "actual_racing",
    "truth_available_to_fit", "formal_fit_started"
  )
  if (!is.data.frame(manifest) || !all(required %in% names(manifest))) {
    stop("Candidate2 R manifest lacks required columns.")
  }
  string_columns <- c("fit_id", "data_id", "initialization_independence_id")
  if (any(vapply(manifest[string_columns], function(value) {
    !is.character(value) || anyNA(value) || any(!nzchar(value))
  }, logical(1)))) stop("Manifest identity fields must be non-empty strings.")
  expected_route <- r_route_v3_candidate2_spec(
    if (isTRUE(formal_experiment)) "formal" else "smoke"
  )$route_id
  if (any(manifest$r_only_protocol_id !=
            R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID) ||
      any(manifest$route_id != expected_route) ||
      any(!as.logical(manifest$terminal_record_required)) ||
      any(as.logical(manifest$actual_racing)) ||
      any(as.logical(manifest$truth_available_to_fit)) ||
      any(as.logical(manifest$formal_fit_started))) {
    stop("Candidate2 R manifest route or isolation fields are invalid.")
  }
  if (anyDuplicated(manifest$fit_id) ||
      anyDuplicated(manifest$initialization_independence_id)) {
    stop("Candidate2 R manifest identities are duplicated.")
  }
  by_data <- split(manifest, manifest$data_id)
  for (rows in by_data) {
    if (nrow(rows) != expected_starts || anyDuplicated(rows$fit_seed) ||
        !identical(sort(as.integer(rows$seed_index)), seq_len(expected_starts))) {
      stop("Each data set must contain all independent registered R starts.")
    }
  }
  invisible(TRUE)
}

r_route_v3_candidate2_exact_join <- function(manifest, terminals) {
  required <- c(
    "fit_id", "data_id", "fit_seed", "seed_index",
    "initialization_independence_id"
  )
  if (!all(required %in% names(manifest)) ||
      !all(required %in% names(terminals))) {
    stop("Exact join fields are missing.")
  }
  key <- function(table) paste(
    as.character(table$fit_id), as.character(table$data_id),
    as.character(table$fit_seed), as.character(table$seed_index),
    as.character(table$initialization_independence_id), sep = "\u001f"
  )
  manifest_key <- key(manifest); terminal_key <- key(terminals)
  if (anyDuplicated(manifest_key) || anyDuplicated(terminal_key) ||
      length(manifest_key) != length(terminal_key) ||
      !setequal(manifest_key, terminal_key)) {
    stop("Terminal records do not exactly join to the FIT_MANIFEST identities.")
  }
  terminals[match(manifest_key, terminal_key), , drop = FALSE]
}

select_R_route_v3_candidate2_winner <- function(
    manifest, terminal_summaries, expected_starts = 12L,
    formal_experiment = TRUE) {
  validate_R_route_v3_candidate2_manifest(
    manifest, expected_starts, formal_experiment = formal_experiment
  )
  required_terminal <- c(
    "protocol_id", "fit_id", "data_id", "method", "route", "fit_seed",
    "seed_index", "initialization_independence_id", "terminal_status",
    "final_elbo", "objective_eligible", "strict_practical_converged",
    "formal_experiment", "truth_used_for_fit_or_selection"
  )
  if (!is.data.frame(terminal_summaries) ||
      !all(required_terminal %in% names(terminal_summaries))) {
    stop("Candidate2 terminal summary lacks objective-only fields.")
  }
  rows <- r_route_v3_candidate2_exact_join(manifest, terminal_summaries)
  expected_route <- r_route_v3_candidate2_spec(
    if (isTRUE(formal_experiment)) "formal" else "smoke"
  )$route_id
  if (any(rows$protocol_id != R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID) ||
      any(rows$route != expected_route) ||
      any(rows$method != "R") ||
      any(as.logical(rows$formal_experiment) != isTRUE(formal_experiment)) ||
      any(as.logical(rows$truth_used_for_fit_or_selection))) {
    stop("Winner selection received a wrong route, non-formal or truth-bearing row.")
  }
  selected <- lapply(split(rows, rows$data_id), function(group) {
    eligible <- as.logical(group$objective_eligible) &
      group$terminal_status %in% c(
        "strict_practical_converged", "max_budget_reached"
      ) & is.finite(as.numeric(group$final_elbo))
    if (!any(eligible)) stop("No objective-eligible R endpoint for data set.")
    candidate <- which(eligible)
    ordered <- candidate[order(
      -as.numeric(group$final_elbo[candidate]),
      as.integer(group$fit_seed[candidate]),
      as.integer(group$seed_index[candidate])
    )]
    winner <- group[ordered[[1L]], , drop = FALSE]
    winner$selection_rule <- "maximum_valid_ordinary_T1_ELBO"
    winner$selected_endpoint_unfinished <-
      !isTRUE(winner$strict_practical_converged[[1L]])
    winner
  })
  result <- do.call(rbind, selected); rownames(result) <- NULL
  result
}

r_route_v3_candidate2_best_basin <- function(terminal_summaries) {
  required <- c(
    "data_id", "fit_id", "fit_seed", "seed_index",
    "initialization_independence_id", "final_elbo", "objective_eligible"
  )
  if (!is.data.frame(terminal_summaries) ||
      !all(required %in% names(terminal_summaries))) {
    stop("Terminal table lacks basin fields.")
  }
  result <- terminal_summaries[, required, drop = FALSE]
  result$elbo_gap_from_best <- NA_real_
  result$best_basin_tolerance <- NA_real_
  result$best_elbo_basin_member <- FALSE
  for (data_id in unique(result$data_id)) {
    member <- which(result$data_id == data_id)
    eligible <- member[as.logical(result$objective_eligible[member]) &
      is.finite(as.numeric(result$final_elbo[member]))]
    if (!length(eligible)) next
    best <- max(as.numeric(result$final_elbo[eligible]))
    tolerance <- max(0.05, 1e-4 * abs(best))
    result$elbo_gap_from_best[member] <-
      best - as.numeric(result$final_elbo[member])
    result$best_basin_tolerance[member] <- tolerance
    result$best_elbo_basin_member[member] <-
      as.logical(result$objective_eligible[member]) &
      is.finite(result$elbo_gap_from_best[member]) &
      result$elbo_gap_from_best[member] <= tolerance
  }
  result
}

r_route_v3_candidate2_offline_racing_replay <- function(
    terminal_summaries, elbo_traces, main_winner,
    screen_t1 = 100L, minimum_keep = 4L, elbo_margin = 500) {
  screen_t1 <- c2_integer(screen_t1, "screen_t1", 1L)
  minimum_keep <- c2_integer(minimum_keep, "minimum_keep", 1L)
  if (!is.list(elbo_traces) || is.null(names(elbo_traces)) ||
      anyDuplicated(names(elbo_traces))) stop("ELBO traces must be named by fit_id.")
  required <- c(
    "fit_id", "data_id", "fit_seed", "seed_index", "actual_T1_sweeps",
    "final_elbo", "objective_eligible"
  )
  if (!all(required %in% names(terminal_summaries))) {
    stop("Terminal summaries lack racing fields.")
  }
  rows <- terminal_summaries
  replay_rows <- lapply(seq_len(nrow(rows)), function(index) {
    fit_id <- as.character(rows$fit_id[[index]])
    trace <- elbo_traces[[fit_id]]
    actual <- as.integer(rows$actual_T1_sweeps[[index]])
    if (!is.data.frame(trace) ||
        !all(c("t1_sweep", "elbo") %in% names(trace)) ||
        !is.finite(actual) || nrow(trace) < actual ||
        max(as.integer(trace$t1_sweep)) != actual) {
      stop("Racing trace is incomplete for fit_id ", fit_id, ".")
    }
    used <- min(actual, screen_t1)
    hit <- which(as.integer(trace$t1_sweep) == used)
    if (length(hit) != 1L) stop("Racing trace lacks its real comparison endpoint.")
    data.frame(
      fit_id = fit_id, data_id = as.character(rows$data_id[[index]]),
      fit_seed = as.integer(rows$fit_seed[[index]]),
      seed_index = as.integer(rows$seed_index[[index]]),
      actual_terminal_T1 = actual, comparison_T1 = used,
      completed_before_screen = actual < screen_t1,
      comparison_elbo = as.numeric(trace$elbo[[hit]]),
      objective_eligible = as.logical(rows$objective_eligible[[index]]),
      stringsAsFactors = FALSE
    )
  })
  replay <- do.call(rbind, replay_rows)
  replay$retained <- FALSE
  replay$estimated_T1_sweeps_saved <- 0L
  for (data_id in unique(replay$data_id)) {
    member <- which(replay$data_id == data_id)
    eligible <- member[replay$objective_eligible[member] &
      is.finite(replay$comparison_elbo[member])]
    if (!length(eligible)) next
    ordered <- eligible[order(
      -replay$comparison_elbo[eligible], replay$fit_seed[eligible],
      replay$seed_index[eligible]
    )]
    best <- replay$comparison_elbo[[ordered[[1L]]]]
    completed <- member[replay$completed_before_screen[member]]
    keep <- union(
      completed,
      union(head(ordered, min(minimum_keep, length(ordered))),
            eligible[best - replay$comparison_elbo[eligible] <= elbo_margin])
    )
    replay$retained[keep] <- TRUE
    eliminated <- setdiff(member, keep)
    replay$estimated_T1_sweeps_saved[eliminated] <- pmax(
      replay$actual_terminal_T1[eliminated] - screen_t1, 0L
    )
  }
  winner_ids <- as.character(main_winner$fit_id)
  summary <- do.call(rbind, lapply(unique(replay$data_id), function(data_id) {
    group <- replay[replay$data_id == data_id, , drop = FALSE]
    winner <- main_winner[main_winner$data_id == data_id, , drop = FALSE]
    if (nrow(winner) != 1L) stop("Main winner is missing for racing replay.")
    retained <- group$retained[group$fit_id == winner$fit_id[[1L]]]
    data.frame(
      data_id = data_id, screen_t1 = screen_t1,
      minimum_keep = minimum_keep, elbo_margin = elbo_margin,
      main_winner_fit_id = as.character(winner$fit_id[[1L]]),
      main_winner_retained = length(retained) == 1L && retained,
      false_elimination_of_main_winner =
        !(length(retained) == 1L && retained),
      retained_count = sum(group$retained),
      estimated_T1_sweeps_saved = sum(group$estimated_T1_sweeps_saved),
      replay_uses_only_real_trace_through_actual_terminal = TRUE,
      replay_changes_main_fit_or_winner = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  list(per_fit = replay, per_data = summary)
}
