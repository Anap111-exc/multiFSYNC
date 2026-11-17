#' G12 fixed-horizon practical stopping control
#'
#' Constructs the opt-in practical-control profile selected by the G12
#' stopping-semantics development audit.  The profile keeps the ordinary
#' practical diagnostics and quantile-factor PPI gate, but makes sweep 400 the
#' first possible stopping point.  Reaching that horizon is not itself a
#' convergence claim: a fit that does not pass the gate remains a slow case.
#'
#' This helper controls only the T = 1 stopping budget.  It does not activate
#' Gram-unit-energy initialization, pre-scoring, annealing, multi-start fitting,
#' or endpoint selection; callers must request those choices explicitly.
#'
#' @return A fully resolved named list suitable for the `practical_control`
#'   argument of [bayesSYNC_multi()] or [bayesSYNC_multi_pre_score()].
#'
#' @export
g12_stopping_control <- function() {
  .validate_practical_control(list(
    min_t1 = 396L,
    window = 20L,
    long_window = 60L,
    consecutive = 5L,
    max_t1 = 400L,
    elbo_abs_rate = 0.03,
    elbo_per_response_rate = 1e-4,
    fitted_nrmse = 1e-3,
    rss_rel = 1e-3,
    ppi_max_abs = 1e-2,
    ppi_gate = "quantile_factor",
    ppi_quantile = 0.95,
    ppi_quantile_max_abs = 1e-2,
    factor_ppi_max_abs = 1e-2,
    long_fitted_nrmse = 3e-3,
    long_rss_rel = 6e-3,
    long_ppi_max_abs = 1e-2,
    long_ppi_quantile_max_abs = 1e-2,
    long_factor_ppi_max_abs = 1e-2,
    checkpoints = integer()
  ))
}
