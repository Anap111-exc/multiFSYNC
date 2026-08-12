# Data-only calibration helpers for the v0-F dimension-information bridge.

v0f_loading_target_norm <- function(
    active_count,
    regime = c("fixed_loading_norm", "constant_per_active_loading"),
    anchor_active = 10L,
    anchor_norm = 3) {
  regime <- match.arg(regime)
  if (length(active_count) != 1L || !is.finite(active_count) ||
      active_count != as.integer(active_count) || active_count < 1L) {
    stop("active_count must be one positive integer.")
  }
  if (length(anchor_active) != 1L || !is.finite(anchor_active) ||
      anchor_active != as.integer(anchor_active) || anchor_active < 1L) {
    stop("anchor_active must be one positive integer.")
  }
  if (length(anchor_norm) != 1L || !is.finite(anchor_norm) ||
      anchor_norm <= 0) {
    stop("anchor_norm must be one positive finite number.")
  }
  if (regime == "fixed_loading_norm") {
    return(as.numeric(anchor_norm))
  }
  as.numeric(anchor_norm) * sqrt(as.integer(active_count) /
                                   as.integer(anchor_active))
}

v0f_rebuild_observations <- function(data) {
  truth <- data$true_params
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[study])) {
      shared <- truth$f_true_values[[study]][[subject]] %*%
        t(truth$a_true)
      specific <- truth$g_true_values[[study]][[subject]] %*%
        t(truth$b_true[[study]])
      signal <- truth$mu_true_values[[study]][[subject]] +
        truth$beta_true_values[[study]][[subject]] +
        shared + specific
      noise <- truth$noise_true_values[[study]][[subject]]
      truth$signal_true_values[[study]][[subject]] <- signal
      data$Y[[study]][[subject]] <- lapply(
        seq_len(truth$p),
        function(variable) as.numeric(signal[, variable] +
                                        noise[, variable])
      )
    }
  }
  data$true_params <- truth
  data
}

calibrate_loading_information_v0f <- function(
    data,
    regime = c("fixed_loading_norm", "constant_per_active_loading"),
    anchor_p = 100L,
    active_fraction = 0.1,
    anchor_norm = 3) {
  regime <- match.arg(regime)
  truth <- data$true_params
  if (truth$L_f != 1L || truth$L_s != 1L || truth$d_use != 0L) {
    stop("v0-F calibration requires L_f=L_s=1 and no covariates.")
  }
  if (!isTRUE(truth$identified_loadings)) {
    stop("v0-F calibration requires identified_loadings=TRUE.")
  }
  if (length(anchor_p) != 1L || !is.finite(anchor_p) ||
      anchor_p != as.integer(anchor_p) || anchor_p < 1L) {
    stop("anchor_p must be one positive integer.")
  }
  if (length(active_fraction) != 1L || !is.finite(active_fraction) ||
      active_fraction <= 0 || active_fraction > 1) {
    stop("active_fraction must be in (0, 1].")
  }
  anchor_active <- as.integer(round(anchor_p * active_fraction))
  normalize_column <- function(value) {
    value <- as.numeric(value)
    active_count <- sum(value != 0)
    current_norm <- sqrt(sum(value^2))
    if (active_count < 1L || !is.finite(current_norm) ||
        current_norm <= 1e-14) {
      stop("Cannot calibrate a zero/non-finite loading column.")
    }
    target <- v0f_loading_target_norm(
      active_count, regime, anchor_active, anchor_norm
    )
    list(value = target * value / current_norm,
         active_count = active_count, target_norm = target)
  }

  calibrated_a <- normalize_column(truth$a_true[, 1L])
  truth$a_true[, 1L] <- calibrated_a$value
  calibrated_b <- vector("list", truth$S)
  for (study in seq_len(truth$S)) {
    calibrated_b[[study]] <-
      normalize_column(truth$b_true[[study]][, 1L])
    truth$b_true[[study]][, 1L] <- calibrated_b[[study]]$value
  }
  truth$gamma_a_true <- 1L * (truth$a_true != 0)
  truth$gamma_b_true <- lapply(truth$b_true, function(value) {
    1L * (value != 0)
  })
  truth$omega_a_true <- colMeans(truth$gamma_a_true)
  truth$omega_b_true <- lapply(truth$gamma_b_true, colMeans)
  truth$v0f_information_regime <- regime
  truth$v0f_anchor_p <- as.integer(anchor_p)
  truth$v0f_active_fraction <- active_fraction
  truth$v0f_anchor_active <- anchor_active
  truth$v0f_anchor_loading_norm <- anchor_norm
  truth$v0f_target_loading_norm_shared <- calibrated_a$target_norm
  truth$v0f_target_loading_norm_specific <- vapply(
    calibrated_b, `[[`, numeric(1), "target_norm"
  )
  truth$loading_calibration <- paste0("v0f_", regime)
  data$true_params <- truth
  v0f_rebuild_observations(data)
}

