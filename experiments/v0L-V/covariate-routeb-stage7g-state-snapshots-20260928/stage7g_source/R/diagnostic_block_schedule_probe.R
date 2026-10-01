# =============================================================================
# diagnostic_block_schedule_probe.R -- private coordinate-schedule experiment
#
# Ordinary fits never set this option and therefore retain the published
# loading -> variance -> omega order.  The experiment-only environment enables
# variance -> omega -> loading without adding a public argument.
# =============================================================================

.multiFSYNC_block_schedule_probe_option <-
  "multiFSYNC.__private_block_schedule_probe_8f3c741a__"

.new_block_schedule_probe_environment <- function(
    schedule = "variance_omega_loading") {
  allowed <- "variance_omega_loading"
  if (!is.character(schedule) || length(schedule) != 1L ||
      is.na(schedule) || !schedule %in% allowed) {
    stop(
      "Private block schedule must be 'variance_omega_loading'."
    )
  }
  environment <- new.env(parent = emptyenv())
  environment$schedule <- schedule
  environment
}

.private_block_schedule_probe <- function() {
  environment <- getOption(
    .multiFSYNC_block_schedule_probe_option, NULL
  )
  if (is.null(environment)) return("loading_variance_omega")
  if (!is.environment(environment) ||
      !identical(environment$schedule, "variance_omega_loading")) {
    stop("The private block schedule probe option is malformed.")
  }
  environment$schedule
}
