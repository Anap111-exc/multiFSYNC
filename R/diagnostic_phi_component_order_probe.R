# =============================================================================
# diagnostic_phi_component_order_probe.R -- private within-factor order probe
#
# Ordinary fits update shared time-function components in ascending order.
# The experiment-only option reverses that coordinate order while preserving
# component labels, initial values, conditional formulas, and all outer blocks.
# =============================================================================

.multiFSYNC_phi_component_order_probe_option <-
  "multiFSYNC.__private_phi_component_order_probe_69bc4e12__"

.new_phi_component_order_probe_environment <- function(
    order = "reverse") {
  if (!is.character(order) || length(order) != 1L ||
      is.na(order) || !identical(order, "reverse")) {
    stop("Private phi component order must be 'reverse'.")
  }
  environment <- new.env(parent = emptyenv())
  environment$order <- order
  environment
}

.private_phi_component_order <- function(component_count) {
  if (!is.numeric(component_count) || length(component_count) != 1L ||
      !is.finite(component_count) || component_count < 1L ||
      abs(component_count - round(component_count)) >
        sqrt(.Machine$double.eps)) {
    stop("component_count must be one positive integer.")
  }
  component_count <- as.integer(round(component_count))
  environment <- getOption(
    .multiFSYNC_phi_component_order_probe_option, NULL
  )
  if (is.null(environment)) return(seq_len(component_count))
  if (!is.environment(environment) ||
      !identical(environment$order, "reverse")) {
    stop("The private phi component order probe option is malformed.")
  }
  rev(seq_len(component_count))
}
