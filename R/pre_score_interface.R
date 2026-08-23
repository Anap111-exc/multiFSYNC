#' Fit multiFSYNC with a bounded pre-score update
#'
#' This is the supported interface for optimization paths that update the
#' shared and study-specific score coordinates before the ordinary first
#' function update. After the requested bounded pre-score sweeps, the fit uses
#' the unchanged `bayesSYNC_multi()` CAVI schedule and ordinary objective.
#'
#' The interface deliberately does not expose the historical diagnostic dense
#' loading gate or the generic random-state calibration control. It supports
#' only the reviewed current and Gram-unit-energy function initializations;
#' all other diagnostic calibration modes remain private.
#'
#' @param ... Arguments passed to [bayesSYNC_multi()]. A `control` argument is
#'   not accepted; use `pre_score_sweeps` and `trace_sweeps` below.
#' @param pre_score_sweeps Number of initial total sweeps that receive the
#'   additional pre-score coordinate update. The R route uses exactly `1L`.
#' @param trace_sweeps Positive total-sweep indices retained in the compact
#'   driver trace. This affects diagnostics only, not numerical updates.
#' @param function_initialization Random function-mean initialization geometry.
#'   `"current_1_over_m"` preserves the existing coefficient initialization.
#'   `"gram_unit_energy"` rescales each randomly drawn shared and study-specific
#'   function mean to unit discrete Gram energy without changing its direction.
#'
#' @return A `bayesSYNC_multi()` fit with `driver_trace`,
#'   `driver_diagnostic_control`, and `pre_score_interface` provenance fields.
#'
#' @export
bayesSYNC_multi_pre_score <- function(
    ..., pre_score_sweeps = 1L, trace_sweeps = 1L,
    function_initialization = c("current_1_over_m", "gram_unit_energy")) {
  dots <- list(...)
  if ("control" %in% names(dots)) {
    stop(
      "bayesSYNC_multi_pre_score() does not accept a diagnostic control; ",
      "use pre_score_sweeps and trace_sweeps."
    )
  }
  validate_integer <- function(value, name, minimum = 0L) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        !is_int(value) || value < minimum || value > .Machine$integer.max) {
      stop(name, " must be one integer of at least ", minimum, ".")
    }
    as.integer(value)
  }
  pre_score_sweeps <- validate_integer(
    pre_score_sweeps, "pre_score_sweeps", minimum = 0L
  )
  if (!is.numeric(trace_sweeps) || !length(trace_sweeps) ||
      any(!is.finite(trace_sweeps)) ||
      any(!vapply(trace_sweeps, is_int, logical(1))) ||
      any(trace_sweeps < 1L) || any(trace_sweeps > .Machine$integer.max)) {
    stop("trace_sweeps must contain positive integers.")
  }
  trace_sweeps <- sort(unique(as.integer(trace_sweeps)))
  function_initialization <- match.arg(function_initialization)
  random_scale_calibration <- switch(
    function_initialization,
    current_1_over_m = "none",
    gram_unit_energy = "function"
  )
  diagnostic_control <- list(
    trace_sweeps = trace_sweeps,
    dense_gate_sweeps = 0L,
    random_scale_calibration = random_scale_calibration,
    pre_score_sweeps = pre_score_sweeps
  )
  fit <- do.call(
    .bayesSYNC_multi_driver_diagnostic,
    c(dots, list(control = diagnostic_control))
  )
  fit$pre_score_interface <- list(
    interface = "bayesSYNC_multi_pre_score",
    interface_version = "1.1.0",
    pre_score_sweeps = pre_score_sweeps,
    trace_sweeps = trace_sweeps,
    dense_gate_sweeps = 0L,
    function_initialization = function_initialization,
    random_scale_calibration = random_scale_calibration
  )
  fit
}
