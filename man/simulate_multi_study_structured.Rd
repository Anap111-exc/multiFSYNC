\name{simulate_multi_study_structured}
\alias{simulate_multi_study_structured}
\title{Simulate Multi-Study Functional Factor Data from External Structured Truths}
\usage{
simulate_multi_study_structured(
  S = 2, n_s = c(20, 30), p = 10, d = 2,
  L_f = 3, L_s = 2, M_f = NULL, M_s = NULL,
  K = NULL, n_obs = 50, common_grid = TRUE,
  sigma_eps = 0.1, bool_sparse_loadings = TRUE,
  prop_sparse = 0.5, score_var_decay = TRUE,
  n_dense = 200, mean_amp = 1, beta_amp = 0.3,
  bs_degree = c(2L, 3L), identified_loadings = FALSE,
  seed = NULL, n_obs_range = NULL,
  sparsity_mode = c("fixed", "factor_beta"),
  omega_beta = c(1, 10)
)
}
\arguments{
  \item{S}{Number of studies.}
  \item{n_s}{Numbers of subjects in the studies.}
  \item{p}{Number of functional variables.}
  \item{d}{Number of scalar covariates.}
  \item{L_f}{Scalar number of shared factors.}
  \item{L_s}{Scalar or length-S vector of study-specific factor counts.}
  \item{M_f}{Numbers of FPCA components in the shared factor blocks.}
  \item{M_s}{Length-S list with one FPCA truncation rank for each of
  the \code{L_s[s]} specific factors in study s.}
  \item{K}{O'Sullivan spline dimension used to return the convenience design
  matrices. It does not define the data-generating truth functions.}
  \item{n_obs}{Legacy fixed number of observations per subject.}
  \item{common_grid}{Whether all subjects use a common regular grid.}
  \item{sigma_eps}{Measurement-error standard deviation.}
  \item{bool_sparse_loadings}{Whether the loading matrices are sparse.}
  \item{prop_sparse}{Legacy fixed proportion of zero loadings.}
  \item{score_var_decay}{Whether theoretical FPCA eigenvalues are proportional
  to \eqn{1/m^2}; otherwise they are equal. The values are normalised within
  each factor to sum to one.}
  \item{n_dense}{Dense-grid size for constructing and orthonormalising truth
  functions.}
  \item{mean_amp, beta_amp}{Amplitudes of the analytic mean and covariate
  functions.}
  \item{bs_degree}{Integer degree pool for the ordinary B-spline truth
  functions.}
  \item{identified_loadings}{Whether to use the identified loading generator.}
  \item{seed}{Optional simulation seed.}
  \item{n_obs_range}{Optional two-integer range from which subject-specific
  observation counts are sampled. Requires \code{common_grid = FALSE}.}
  \item{sparsity_mode}{\code{"fixed"} uses \code{prop_sparse};
  \code{"factor_beta"} draws a separate loading-inclusion probability for
  every factor.}
  \item{omega_beta}{Two positive Beta shape parameters used in
  \code{sparsity_mode = "factor_beta"}.}
}
\value{
A list with responses \code{Y}, covariates \code{Z}, observation times
\code{time_obs}, convenience O'Sullivan designs \code{C}, and
\code{true_params}. The legacy truth fields \code{zeta_true} and
\code{xi_true} are the unit-variance computational scores \eqn{\eta} and
\eqn{\chi}. Theoretical FPCA scores, unit-norm eigenfunctions, eigenvalues,
and the scaled computational time functions are stored separately.
}
\description{
Generates multi-study functional factor data from ordinary B-spline truth
functions that are evaluated directly at the observation times. The truth
functions are not projected into the O'Sullivan spline space used for model
fitting. This provides a Jaoua-style external-basis simulation while retaining
the original public argument names.
}
