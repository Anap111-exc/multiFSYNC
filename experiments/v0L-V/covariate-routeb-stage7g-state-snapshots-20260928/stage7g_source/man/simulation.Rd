\name{simulate_multi_study_data}
\alias{simulate_multi_study_data}
\alias{simulate_multi_study_osullivan}
\alias{simulate_multi_study_structured}
\alias{generate_explicit_eigenfunctions}
\title{Simulation Utilities for multiFSYNC}
\description{
Generate small or formal multi-study continuous functional datasets. The
O'Sullivan generator can use pooled common knots, study-specific mean shifts,
genuinely different study-specific time functions, and loadings constructed to
satisfy the population Gram conditions.
}
\usage{
simulate_multi_study_data(S = 2, n_s = c(20, 30), p = 10, d = 2,
  L_f = 3, L_s = 2, M_f = NULL, M_s = NULL, K = NULL, n_obs = 50,
  common_grid = TRUE, sigma_eps = 0.1, bool_sparse_loadings = TRUE,
  prop_sparse = 0.7, identified_loadings = FALSE, seed = NULL)

simulate_multi_study_osullivan(S = 2, n_s = c(30, 35), p = 10,
  d = 0, L_f = 3, L_s = 2, M_f = NULL, M_s = NULL, K = NULL,
  n_obs = 50, common_grid = TRUE, sigma_eps = 0.1,
  bool_sparse_loadings = TRUE, prop_sparse = 0.5, sd_nu_mu = 0.5,
  sd_nu_beta = 0.3, sd_nu_phi = 1, sd_nu_psi = 1, seed = NULL,
  use_explicit = FALSE, mean_structure = c("shared", "study_specific"),
  specific_time_heterogeneity = TRUE, identified_loadings = FALSE)

simulate_multi_study_structured(S = 2, n_s = c(20, 30), p = 10,
  d = 2, L_f = 3, L_s = 2, M_f = NULL, M_s = NULL, K = NULL,
  n_obs = 50, common_grid = TRUE, sigma_eps = 0.1,
  bool_sparse_loadings = TRUE, prop_sparse = 0.5,
  score_var_decay = TRUE, n_dense = 200, mean_amp = 1,
  beta_amp = 0.3, bs_degree = 3, identified_loadings = FALSE,
  seed = NULL)

generate_explicit_eigenfunctions(K, L_f, M_f, L_s = 0, M_s = NULL,
  S = 1, n_g = 200, int_knots = NULL,
  specific_time_heterogeneity = TRUE)
}
\arguments{
  \item{S, n_s, p, d}{Study, sample, variable and covariate dimensions.}
  \item{K}{O'Sullivan spline dimension.}
  \item{L_f, M_f}{Shared factor count and within-factor FPCA truncations.}
  \item{L_s}{A non-negative integer scalar, recycled across studies, or a
  length-S vector of study-specific factor counts.}
  \item{M_s}{A length-S list satisfying
  \code{length(M_s[[s]]) = L_s[s]}.}
  \item{n_obs, n_g, n_dense}{Observation or construction-grid sizes.}
  \item{common_grid, sigma_eps}{Observation-grid and noise controls.}
  \item{bool_sparse_loadings, prop_sparse, identified_loadings}{Loading controls.}
  \item{sd_nu_mu, sd_nu_beta, sd_nu_phi, sd_nu_psi}{Coefficient scales.}
  \item{seed, use_explicit}{Reproducibility and explicit-function controls.}
  \item{mean_structure, specific_time_heterogeneity}{Study-heterogeneity controls.}
  \item{score_var_decay, mean_amp, beta_amp, bs_degree}{Structured-simulation controls.}
  \item{int_knots}{Optional common O'Sullivan interior knots.}
}
\value{A simulated dataset with true parameters, or explicit spline coefficients.}
