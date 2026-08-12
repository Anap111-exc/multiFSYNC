\name{multiFSYNC-utilities}
\alias{get_grid_objects}
\alias{set_hyper}
\alias{orthonormalise_multi}
\alias{trapint}
\alias{flip_sign}
\alias{match_factor_and_sign}
\title{multiFSYNC Grid, Hyperparameter and Post-processing Utilities}
\description{
Utilities for O'Sullivan grid construction, hyperparameter specification,
posterior-moment FPCA recovery, numerical integration, and simulation-only
factor matching. Truth-based matching is not used by real-data post-processing.
}
\usage{
get_grid_objects(time_obs, K, n_g = 1000, time_g = NULL,
  int_knots = NULL, format_univ = FALSE)
set_hyper(sigma_beta = 1e5, A = 1e5, c_0 = 1, d_0 = NULL)
orthonormalise_multi(C_g, time_g, mu_q_nu_mu, mu_q_nu_beta = NULL,
  mu_q_nu_phi, mu_q_nu_psi, mu_q_zeta, Sigma_q_zeta, mu_q_xi,
  Sigma_q_xi, mu_q_a, mu_q_b_specific, mu_q_gamma_a,
  mu_q_gamma_b, S, n_s, p, L_f, L_s, M_f, M_s,
  d = if (is.null(mu_q_nu_beta)) 0L else length(mu_q_nu_beta[[1]]),
  response_center = NULL, response_scale = NULL,
  Sigma_q_nu_phi = NULL, Sigma_q_nu_psi = NULL,
  zeta_true = NULL, pve_cutoff = 0.99)
trapint(xgrid, fgrid)
flip_sign(vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse = NULL)
match_factor_and_sign(B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat,
  list_list_Phi_hat, list_cumulated_pve, list_Cov_zeta_hat,
  list_h_hat = NULL, list_var_vec = NULL)
}
\arguments{
  \item{time_obs, K, n_g, time_g, int_knots, format_univ}{Grid-construction inputs.}
  \item{sigma_beta, A, c_0, d_0}{Prior hyperparameters.}
  \item{xgrid, fgrid}{Grid and values for trapezoidal integration.}
  \item{vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse}{Sign-alignment inputs.}
  \item{C_g, mu_q_nu_mu, mu_q_nu_beta, mu_q_nu_phi, mu_q_nu_psi}{Spline-grid posterior inputs.}
  \item{mu_q_zeta, Sigma_q_zeta, mu_q_xi, Sigma_q_xi}{Score posterior inputs.}
  \item{mu_q_a, mu_q_b_specific, mu_q_gamma_a, mu_q_gamma_b}{Loading posterior inputs.}
  \item{S, n_s, p, L_f, L_s, M_f, M_s, d}{Model dimensions.}
  \item{response_center, response_scale}{Response scaling metadata.}
  \item{Sigma_q_nu_phi, Sigma_q_nu_psi, zeta_true, pve_cutoff}{FPCA post-processing controls.}
  \item{B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat}{Simulation matching inputs.}
  \item{list_list_Phi_hat, list_cumulated_pve, list_Cov_zeta_hat, list_h_hat, list_var_vec}{Additional simulation matching inputs.}
}
\value{A function-specific list, matrix, vector, or hyperparameter object.}
