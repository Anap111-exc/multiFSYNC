# =============================================================================
# diagnostic_phi_driver_decomposition.R -- private first-sweep decomposition
#
# This hook records exact additive pieces of one shared time-function
# coordinate. It is absent from ordinary fits, draws no random numbers, and
# can request an experiment-only return immediately after the first phi block.
# =============================================================================

.multiFSYNC_phi_driver_decomposition_option <-
  "multiFSYNC.__private_phi_driver_decomposition_d35a729c__"

.validate_phi_driver_decomposition_control <- function(control) {
  if (!is.list(control) || is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control))) {
    stop("Phi-driver decomposition control must be uniquely named.")
  }
  allowed <- c(
    "iteration", "factor", "component", "stop_after_phi"
  )
  unknown <- setdiff(names(control), allowed)
  if (length(unknown)) {
    stop(
      "Unknown phi-driver decomposition field(s): ",
      paste(unknown, collapse = ", "), "."
    )
  }
  required <- c("iteration", "factor", "component")
  if (!all(required %in% names(control))) {
    stop(
      "Phi-driver decomposition requires iteration, factor, and component."
    )
  }
  for (field in required) {
    value <- control[[field]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 ||
        abs(value - round(value)) > sqrt(.Machine$double.eps) ||
        value > .Machine$integer.max) {
      stop("control$", field, " must be one positive integer.")
    }
    control[[field]] <- as.integer(round(value))
  }
  if (!"stop_after_phi" %in% names(control)) {
    control$stop_after_phi <- FALSE
  }
  if (!is.logical(control$stop_after_phi) ||
      length(control$stop_after_phi) != 1L ||
      is.na(control$stop_after_phi)) {
    stop("control$stop_after_phi must be TRUE or FALSE.")
  }
  control[allowed]
}

.new_phi_driver_decomposition_environment <- function(control) {
  environment <- new.env(parent = emptyenv())
  environment$control <-
    .validate_phi_driver_decomposition_control(control)
  environment$context <- list(
    iteration = NA_integer_,
    temperature = NA_real_,
    annealing = NA
  )
  environment$records <- list()
  environment
}

.get_phi_driver_decomposition_environment <- function() {
  environment <- getOption(
    .multiFSYNC_phi_driver_decomposition_option, NULL
  )
  if (is.null(environment)) return(NULL)
  if (!is.environment(environment) ||
      is.null(environment$control) ||
      is.null(environment$context) ||
      is.null(environment$records)) {
    stop("The private phi-driver decomposition option is malformed.")
  }
  environment
}

.phi_driver_decomposition_set_context <- function(
    iteration, temperature, annealing) {
  environment <- .get_phi_driver_decomposition_environment()
  if (is.null(environment)) return(invisible(NULL))
  environment$context <- list(
    iteration = as.integer(iteration),
    temperature = as.numeric(temperature),
    annealing = isTRUE(annealing)
  )
  invisible(NULL)
}

.phi_driver_decomposition_should_trace <- function(factor, component) {
  environment <- .get_phi_driver_decomposition_environment()
  if (is.null(environment)) return(FALSE)
  control <- environment$control
  context <- environment$context
  identical(context$iteration, control$iteration) &&
    identical(as.integer(factor), control$factor) &&
    identical(as.integer(component), control$component)
}

.phi_driver_decomposition_record <- function(record) {
  environment <- .get_phi_driver_decomposition_environment()
  if (is.null(environment)) return(invisible(NULL))
  environment$records[[length(environment$records) + 1L]] <- record
  invisible(NULL)
}

.phi_driver_decomposition_stop_after_phi <- function() {
  environment <- .get_phi_driver_decomposition_environment()
  if (is.null(environment)) return(FALSE)
  isTRUE(environment$control$stop_after_phi) &&
    length(environment$records) > 0L &&
    identical(
      environment$context$iteration,
      environment$control$iteration
    )
}

.phi_driver_decomposition_early_result <- function() {
  environment <- .get_phi_driver_decomposition_environment()
  if (is.null(environment)) {
    stop("No private phi-driver decomposition is active.")
  }
  structure(
    list(
      private_diagnostic = "phi_driver_decomposition",
      stopped_after_phi = TRUE,
      record_count = length(environment$records),
      control = environment$control
    ),
    class = c("multiFSYNC_private_diagnostic", "list")
  )
}
