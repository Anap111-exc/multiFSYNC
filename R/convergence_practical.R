# Practical convergence diagnostics for the T = 1 CAVI phase.
#
# The ordinary ELBO remains the model-derived objective.  This rule adds a
# finite-window stability check for identifiable, scientifically used outputs
# so that slow Half-Cauchy auxiliary coordinates do not by themselves prevent
# a usable fit from stopping.

.default_practical_control <- function() {
  list(
    min_t1 = 80L,
    window = 20L,
    long_window = 60L,
    consecutive = 5L,
    max_t1 = 200L,
    elbo_abs_rate = 0.03,
    elbo_per_response_rate = 1e-4,
    fitted_nrmse = 1e-3,
    rss_rel = 1e-3,
    ppi_max_abs = 1e-2,
    ppi_gate = "max",
    ppi_quantile = 0.95,
    ppi_quantile_max_abs = 1e-2,
    factor_ppi_max_abs = 1e-2,
    long_fitted_nrmse = 3e-3,
    long_rss_rel = 6e-3,
    long_ppi_max_abs = 1e-2,
    long_ppi_quantile_max_abs = 1e-2,
    long_factor_ppi_max_abs = 1e-2,
    checkpoints = integer()
  )
}

.validate_practical_control <- function(control = NULL) {
  defaults <- .default_practical_control()
  if (is.null(control)) return(defaults)
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control)))) {
    stop("practical_control must be NULL or a fully named list.")
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown practical_control field(s): ",
         paste(unknown, collapse = ", "))
  }
  resolved <- utils::modifyList(defaults, control)
  # Before long-window tolerances were exposed, a supplied short tolerance
  # controlled both windows. Preserve that behaviour for existing custom
  # controls while allowing the production defaults to calibrate the two
  # cumulative horizons separately.
  legacy_tolerance_map <- c(
    long_fitted_nrmse = "fitted_nrmse",
    long_rss_rel = "rss_rel",
    long_ppi_max_abs = "ppi_max_abs",
    long_ppi_quantile_max_abs = "ppi_quantile_max_abs",
    long_factor_ppi_max_abs = "factor_ppi_max_abs"
  )
  for (long_field in names(legacy_tolerance_map)) {
    short_field <- legacy_tolerance_map[[long_field]]
    if (!long_field %in% names(control) &&
        short_field %in% names(control)) {
      resolved[[long_field]] <- control[[short_field]]
    }
  }

  integer_fields <- c(
    "min_t1", "window", "long_window", "consecutive", "max_t1")
  for (field in integer_fields) {
    value <- resolved[[field]]
    if (length(value) != 1L || !is.finite(value) ||
        !is_int(value) || value < 1L) {
      stop("practical_control$", field,
           " must be one positive integer.")
    }
    resolved[[field]] <- as.integer(value)
  }
  if (resolved$long_window < resolved$window) {
    stop("practical_control$long_window must be at least window.")
  }
  if (resolved$min_t1 <= resolved$long_window) {
    stop("practical_control$min_t1 must exceed long_window.")
  }
  if (resolved$max_t1 < resolved$min_t1) {
    stop("practical_control$max_t1 must be at least min_t1.")
  }

  tolerance_fields <- c(
    "elbo_abs_rate", "elbo_per_response_rate", "fitted_nrmse",
    "rss_rel", "ppi_max_abs", "ppi_quantile_max_abs",
    "factor_ppi_max_abs",
    "long_fitted_nrmse", "long_rss_rel",
    "long_ppi_max_abs", "long_ppi_quantile_max_abs",
    "long_factor_ppi_max_abs")
  for (field in tolerance_fields) {
    value <- resolved[[field]]
    if (length(value) != 1L || !is.finite(value) || value < 0) {
      stop("practical_control$", field,
           " must be one finite non-negative number.")
    }
  }

  if (length(resolved$ppi_gate) != 1L ||
      !is.character(resolved$ppi_gate) ||
      is.na(resolved$ppi_gate) ||
      !resolved$ppi_gate %in% c("max", "quantile_factor")) {
    stop(
      "practical_control$ppi_gate must be 'max' or ",
      "'quantile_factor'."
    )
  }
  if (length(resolved$ppi_quantile) != 1L ||
      !is.finite(resolved$ppi_quantile) ||
      resolved$ppi_quantile <= 0 ||
      resolved$ppi_quantile > 1) {
    stop(
      "practical_control$ppi_quantile must be one finite number ",
      "in (0, 1]."
    )
  }

  checkpoints <- resolved$checkpoints
  if (is.null(checkpoints) || !length(checkpoints)) {
    resolved$checkpoints <- integer()
  } else {
    if (!is.numeric(checkpoints) || any(!is.finite(checkpoints)) ||
        any(!vapply(checkpoints, is_int, logical(1))) ||
        any(checkpoints < 1L)) {
      stop("practical_control$checkpoints must contain positive integers.")
    }
    resolved$checkpoints <- sort(unique(as.integer(checkpoints)))
  }
  resolved
}

.practical_relative_l2 <- function(current, reference) {
  current <- as.numeric(current)
  reference <- as.numeric(reference)
  if (!length(current) && !length(reference)) return(0)
  if (length(current) != length(reference) ||
      any(!is.finite(current)) || any(!is.finite(reference))) {
    return(Inf)
  }
  sqrt(sum((current - reference)^2)) /
    (1 + sqrt(sum(reference^2)))
}

.practical_fitted_nrmse <- function(current, reference) {
  current <- as.numeric(current)
  reference <- as.numeric(reference)
  if (!length(current) && !length(reference)) return(0)
  if (length(current) != length(reference) ||
      any(!is.finite(current)) || any(!is.finite(reference))) {
    return(Inf)
  }
  sqrt(mean((current - reference)^2)) /
    (1 + sqrt(mean(reference^2)))
}

.practical_max_abs <- function(current, reference) {
  current <- as.numeric(current)
  reference <- as.numeric(reference)
  if (!length(current) && !length(reference)) return(0)
  if (length(current) != length(reference) ||
      any(!is.finite(current)) || any(!is.finite(reference))) {
    return(Inf)
  }
  max(abs(current - reference))
}

.practical_quantile_abs <- function(current, reference, probability) {
  current <- as.numeric(current)
  reference <- as.numeric(reference)
  if (!length(current) && !length(reference)) return(0)
  if (length(current) != length(reference) ||
      any(!is.finite(current)) || any(!is.finite(reference))) {
    return(Inf)
  }
  as.numeric(stats::quantile(
    abs(current - reference),
    probs = probability, names = FALSE, type = 8
  ))
}

.practical_factor_ppi <- function(gamma) {
  if (is.null(gamma) || !length(gamma)) return(numeric())
  gamma <- as.matrix(gamma)
  if (!nrow(gamma) || !ncol(gamma)) return(numeric())
  gamma <- pmin(pmax(gamma, 0), 1)
  1 - apply(1 - gamma, 2, prod)
}

.practical_snapshot <- function(fitted_values, expected_rss,
                                mu_q_gamma_a, mu_q_gamma_b) {
  fitted_values <- as.numeric(fitted_values)
  rss <- as.numeric(unlist(
    expected_rss, recursive = TRUE, use.names = FALSE))
  ppi <- c(
    as.numeric(unlist(
      mu_q_gamma_a, recursive = TRUE, use.names = FALSE)),
    as.numeric(unlist(
      mu_q_gamma_b, recursive = TRUE, use.names = FALSE))
  )
  factor_ppi_a <- .practical_factor_ppi(mu_q_gamma_a)
  factor_ppi_b <- if (is.list(mu_q_gamma_b)) {
    as.numeric(unlist(
      lapply(mu_q_gamma_b, .practical_factor_ppi),
      recursive = TRUE, use.names = FALSE
    ))
  } else {
    .practical_factor_ppi(mu_q_gamma_b)
  }
  list(
    fitted = fitted_values, rss = rss, ppi = ppi,
    factor_ppi = c(factor_ppi_a, factor_ppi_b)
  )
}

.practical_compare <- function(
    current, reference, control, long_window = FALSE) {
  fitted_nrmse <- .practical_fitted_nrmse(
    current$fitted, reference$fitted)
  rss_rel <- .practical_relative_l2(current$rss, reference$rss)
  ppi_max_abs <- .practical_max_abs(current$ppi, reference$ppi)
  ppi_quantile_abs <- .practical_quantile_abs(
    current$ppi, reference$ppi, control$ppi_quantile
  )
  current_factor_ppi <- if (is.null(current$factor_ppi)) {
    numeric()
  } else {
    current$factor_ppi
  }
  reference_factor_ppi <- if (is.null(reference$factor_ppi)) {
    numeric()
  } else {
    reference$factor_ppi
  }
  factor_ppi_max_abs <- .practical_max_abs(
    current_factor_ppi, reference_factor_ppi
  )
  finite <- all(is.finite(c(
    fitted_nrmse, rss_rel, ppi_max_abs,
    ppi_quantile_abs, factor_ppi_max_abs
  )))
  fitted_tolerance <- if (long_window) {
    control$long_fitted_nrmse
  } else {
    control$fitted_nrmse
  }
  rss_tolerance <- if (long_window) {
    control$long_rss_rel
  } else {
    control$rss_rel
  }
  ppi_tolerance <- if (long_window) {
    control$long_ppi_max_abs
  } else {
    control$ppi_max_abs
  }
  ppi_quantile_tolerance <- if (long_window) {
    control$long_ppi_quantile_max_abs
  } else {
    control$ppi_quantile_max_abs
  }
  factor_ppi_tolerance <- if (long_window) {
    control$long_factor_ppi_max_abs
  } else {
    control$factor_ppi_max_abs
  }
  ppi_pass <- if (control$ppi_gate == "max") {
    ppi_max_abs <= ppi_tolerance
  } else {
    ppi_quantile_abs <= ppi_quantile_tolerance &&
      factor_ppi_max_abs <= factor_ppi_tolerance
  }
  pass <- finite &&
    fitted_nrmse <= fitted_tolerance &&
    rss_rel <= rss_tolerance &&
    ppi_pass
  list(
    fitted_nrmse = fitted_nrmse,
    rss_rel = rss_rel,
    ppi_max_abs = ppi_max_abs,
    ppi_quantile_abs = ppi_quantile_abs,
    factor_ppi_max_abs = factor_ppi_max_abs,
    ppi_pass = ppi_pass,
    finite = finite,
    pass = pass
  )
}

.practical_window_diagnostic <- function(
    sweep, ELBO, ELBO_diagnostics, current_snapshot,
    reference_snapshot, long_reference_snapshot, control,
    objective_size = 1) {
  if (length(objective_size) != 1L || !is.finite(objective_size) ||
      objective_size < 1) {
    stop("objective_size must be one finite number of at least one.")
  }
  window <- control$window
  reference_sweep <- sweep - window
  long_window <- control$long_window
  long_reference_sweep <- sweep - long_window
  if (reference_sweep < 1L || long_reference_sweep < 1L ||
      length(ELBO) < sweep || is.null(reference_snapshot) ||
      is.null(long_reference_snapshot)) {
    return(list(
      eligible = FALSE, objective_pass = FALSE,
      short_output_pass = FALSE, long_output_pass = FALSE,
      output_pass = FALSE, pass = FALSE,
      elbo_abs_rate = Inf, elbo_per_response_rate = Inf,
      elbo_rel_rate = Inf,
      fitted_nrmse = Inf, rss_rel = Inf, ppi_max_abs = Inf,
      ppi_quantile_abs = Inf, factor_ppi_max_abs = Inf,
      long_fitted_nrmse = Inf, long_rss_rel = Inf,
      long_ppi_max_abs = Inf,
      long_ppi_quantile_abs = Inf,
      long_factor_ppi_max_abs = Inf,
      monotone = FALSE, long_monotone = FALSE, finite = FALSE
    ))
  }
  current_elbo <- ELBO[sweep]
  reference_elbo <- ELBO[reference_sweep]
  elbo_change <- abs(current_elbo - reference_elbo)
  elbo_abs_rate <- elbo_change / window
  elbo_per_response_rate <- elbo_abs_rate / objective_size
  elbo_rel_rate <- elbo_change /
    (window * (1 + abs(reference_elbo)))

  differences <- diff(ELBO[reference_sweep:sweep])
  tolerance_index <- reference_sweep:(sweep - 1L)
  scaled_tolerance <- ELBO_diagnostics$scaled_tolerance[tolerance_index]
  monotone <- length(differences) == window &&
    length(scaled_tolerance) == window &&
    all(is.finite(differences)) &&
    all(is.finite(scaled_tolerance)) &&
    all(differences >= -scaled_tolerance)

  long_differences <- diff(
    ELBO[long_reference_sweep:sweep])
  long_tolerance_index <- long_reference_sweep:(sweep - 1L)
  long_scaled_tolerance <-
    ELBO_diagnostics$scaled_tolerance[long_tolerance_index]
  long_monotone <- length(long_differences) == long_window &&
    length(long_scaled_tolerance) == long_window &&
    all(is.finite(long_differences)) &&
    all(is.finite(long_scaled_tolerance)) &&
    all(long_differences >= -long_scaled_tolerance)

  short_output_change <- .practical_compare(
    current_snapshot, reference_snapshot, control)
  long_output_change <- .practical_compare(
    current_snapshot, long_reference_snapshot, control,
    long_window = TRUE)
  finite <- all(is.finite(c(
    current_elbo, reference_elbo,
    elbo_abs_rate, elbo_per_response_rate, elbo_rel_rate,
    short_output_change$fitted_nrmse,
    short_output_change$rss_rel,
    short_output_change$ppi_max_abs,
    short_output_change$ppi_quantile_abs,
    short_output_change$factor_ppi_max_abs,
    long_output_change$fitted_nrmse,
    long_output_change$rss_rel,
    long_output_change$ppi_max_abs,
    long_output_change$ppi_quantile_abs,
    long_output_change$factor_ppi_max_abs
  )))
  eligible <- sweep >= control$min_t1
  # The relative rate is retained as a scale diagnostic, but it is not a
  # stopping gate: the complete ELBO can be close to zero because large
  # positive and negative components cancel. The total absolute-rate gate is
  # retained for small problems, while the per-scalar-response rate makes the
  # same practical rule comparable across sample sizes and dimensions.
  objective_pass <- eligible && finite && monotone && long_monotone &&
    (elbo_abs_rate <= control$elbo_abs_rate ||
       elbo_per_response_rate <= control$elbo_per_response_rate)
  short_output_pass <- short_output_change$pass
  long_output_pass <- long_output_change$pass
  output_pass <- short_output_pass && long_output_pass
  pass <- objective_pass && output_pass

  list(
    eligible = eligible,
    objective_pass = objective_pass,
    short_output_pass = short_output_pass,
    long_output_pass = long_output_pass,
    output_pass = output_pass,
    pass = pass,
    elbo_abs_rate = elbo_abs_rate,
    elbo_per_response_rate = elbo_per_response_rate,
    elbo_rel_rate = elbo_rel_rate,
    fitted_nrmse = short_output_change$fitted_nrmse,
    rss_rel = short_output_change$rss_rel,
    ppi_max_abs = short_output_change$ppi_max_abs,
    ppi_quantile_abs = short_output_change$ppi_quantile_abs,
    factor_ppi_max_abs = short_output_change$factor_ppi_max_abs,
    long_fitted_nrmse = long_output_change$fitted_nrmse,
    long_rss_rel = long_output_change$rss_rel,
    long_ppi_max_abs = long_output_change$ppi_max_abs,
    long_ppi_quantile_abs =
      long_output_change$ppi_quantile_abs,
    long_factor_ppi_max_abs =
      long_output_change$factor_ppi_max_abs,
    monotone = monotone,
    long_monotone = long_monotone,
    finite = finite
  )
}
