\name{bayesSYNC_multi}
\alias{bayesSYNC_multi}
\title{Fit the multiFSYNC Model}
\usage{
bayesSYNC_multi(Y, Z = NULL, time_obs = NULL, L_f, L_s, M_f, M_s,
  K = NULL, anneal = c(1, 1.9, 100), list_hyper = NULL,
  n_g = 1000, time_g = NULL, tol_abs = 1e-3, tol_rel = 1e-5,
  maxit = 1000, n_cpus = 1, verbose = TRUE, seed = NULL,
  bool_scale = TRUE, bool_var_spec_prob = FALSE, d_0 = NULL,
  convergence_rule = c("parameters", "elbo", "practical"),
  lambda_orth = 0, practical_control = NULL,
  initialization = c("random", "residual_fpca"),
  initialization_control = NULL)
}
\arguments{
  \item{Y}{Study-by-individual-by-variable nested list of response vectors.}
  \item{Z}{Optional list of study-specific scalar-covariate matrices.}
  \item{time_obs}{Study-by-individual list of observation-time vectors in [0,1].}
  \item{L_f}{Scalar number of shared factors.}
  \item{L_s}{Number of study-specific factors. A scalar is recycled across
  studies; a length-S vector permits a different count in each study.}
  \item{M_f}{FPCA truncation ranks for the shared factors.}
  \item{M_s}{Length-S list satisfying
  \code{length(M_s[[s]]) = L_s[s]}; its entries are the study-specific FPCA
  truncation ranks.}
  \item{K}{Number of nonlinear O'Sullivan basis functions.}
  \item{anneal}{NULL or c(spacing type, maximum temperature, grid size).}
  \item{list_hyper}{Object returned by \code{set_hyper}.}
  \item{n_g, time_g}{Dense output-grid size or a supplied increasing grid.}
  \item{tol_abs, tol_rel, maxit}{Stopping tolerances and maximum iterations.}
  \item{n_cpus}{Requested internal cores; Windows falls back to one.}
  \item{verbose, seed}{Progress and random seed controls.}
  \item{bool_scale}{Whether to fit on a standardised response scale.}
  \item{bool_var_spec_prob}{Use variable-factor-level rather than factor-level omega.}
  \item{d_0}{Optional Beta prior shape overriding the default p.}
  \item{convergence_rule}{Stopping rule. \code{"parameters"} is retained as
  the legacy-compatible default. \code{"elbo"} uses adjacent complete
  ordinary-model ELBO increments after full \eqn{T=1} sweeps.
  \code{"practical"} combines ELBO safeguards with fitted-value, expected-RSS,
  and PPI stability over short and long windows and marks an unfinished fit as
  a slow case at its finite \eqn{T=1} limit. ELBO-based rules are valid only
  when \code{lambda_orth = 0} and no adaptive Cholesky jitter is used during
  the \eqn{T=1} phase.}
  \item{lambda_orth}{Optional separation regulariser; zero by default.
  Positive values are not strict identification constraints and force the
  ordinary-model ELBO into diagnostic-only mode.}
  \item{practical_control}{Optional fully named list overriding practical-rule
  defaults. The production defaults use a minimum of 80 \eqn{T=1} sweeps,
  20- and 60-sweep windows, five consecutive passes, and a 200-sweep
  slow-case cap. The ELBO objective gate accepts either a total absolute
  change rate no larger than 0.03 or a per-scalar-response rate no larger
  than \eqn{10^{-4}}. Recent
  fitted/RSS/PPI tolerances are 0.001/0.001/0.01; cumulative 60-sweep
  tolerances are 0.003/0.006/0.01.}
  \item{initialization}{Iteration-zero construction. \code{"random"} preserves
  the original seeded path. \code{"residual_fpca"} uses only working-scale
  observed data to initialize mean/covariate functions and shared/specific
  factor means before the unchanged CAVI updates begin.}
  \item{initialization_control}{Optional fully named list overriding the
  residual-FPCA \code{grid_size}, \code{perturb_sd}, and relative
  \code{rank_tol}.}
}
\value{
A fitted object containing working-scale variational parameters, original-scale
scientific summaries, posterior-moment FPCA output, and convergence/annealing
status fields. ELBO output includes the numeric trace, named components,
finite-value and decline diagnostics, current expected RSS, objective-validity
flags, adaptive-Cholesky diagnostics, practical stopping diagnostics, and a
compact truth-free stability snapshot for multi-start comparison.
Initialization output records the method used, explained residual fraction,
per-block fallback reasons, construction time, and an independence identifier
used by the multi-start support audit.
}
\description{Fits the multi-study Bayesian functional factor model by CAVI.}
