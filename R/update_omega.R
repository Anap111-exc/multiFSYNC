# =============================================================================
# update_omega.R — Omega (sparsity) hyperparameter update
#
# Covers parameter block 4.13:
#   omega_l   — shared factor sparsity probability
#   omega_tilde_sl — study-specific factor sparsity probability
#
# Variational form: q*(omega) = Beta(c_1, d_1)
#
# Key formulas (derivations.md §3.13):
#   c_1 = c*(c_0 + sum_j PPI_jl) - c + 1
#   d_1 = c*(d_0 + p - sum_j PPI_jl) - c + 1
#   E[log omega]   = digamma(c_1) - digamma(c_1 + d_1)
#   E[log(1-omega)] = digamma(d_1) - digamma(c_1 + d_1)
#
# Temperature c appears in c_1, d_1 computation.
#
# Based on bayesSYNC (GPL-3, hruffieux/bayesSYNC).
# =============================================================================

#' Update shared factor sparsity hyperparameter omega_l
#'
#' @param mu_q_gamma_a p x L_f matrix, PPI_{jl} for shared factors
#' @param c_0 Beta prior shape1 (default 1)
#' @param d_0 Beta prior shape2 (default p)
#' @param c_val Temperature constant c = 1/T
#' @param bool_var_spec_prob Use variable-specific probabilities
#' @return List with c_1_omega_a, d_1_omega_a, mu_q_log_omega_a, mu_q_log_1_omega_a
#'
#' @noRd
update_omega_shared <- function(mu_q_gamma_a, p, c_0 = 1, d_0 = NULL,
                                 c_val = 1, bool_var_spec_prob = FALSE) {

  if (is.null(d_0)) d_0 <- p

  L_f <- ncol(mu_q_gamma_a)

  if (bool_var_spec_prob) {

    c_1_omega_a <- c_val * (c_0 + mu_q_gamma_a) - c_val + 1
    d_1_omega_a <- c_val * (d_0 + 1 - mu_q_gamma_a) - c_val + 1
    dig_a <- digamma(c_val * (c_0 + d_0 + 1) - 2 * c_val + 2)

    mu_q_log_omega_a <- digamma(c_1_omega_a) - dig_a
    mu_q_log_1_omega_a <- digamma(d_1_omega_a) - dig_a

  } else {

    cs_mu_q_gamma_a <- colSums(mu_q_gamma_a)
    c_1_omega_a <- c_val * (c_0 + cs_mu_q_gamma_a) - c_val + 1
    d_1_omega_a <- c_val * (d_0 + p - cs_mu_q_gamma_a) - c_val + 1
    dig_a <- digamma(c_val * (c_0 + d_0 + p) - 2 * c_val + 2)

    mu_q_log_omega_a <- digamma(c_1_omega_a) - dig_a
    mu_q_log_1_omega_a <- digamma(d_1_omega_a) - dig_a

  }

  create_named_list(c_1_omega_a, d_1_omega_a, mu_q_log_omega_a, mu_q_log_1_omega_a)
}

#' Update study-specific factor sparsity hyperparameter omega_tilde_{sl}
#'
#' @param mu_q_gamma_b List of S matrices, each p x L_s, PPI for specific factors
#' @param p Number of variables
#' @param S Number of studies
#' @param c_0 Beta prior shape1
#' @param d_0 Beta prior shape2
#' @param c_val Temperature constant
#' @param bool_var_spec_prob Use variable-specific probabilities
#' @return List with c_1_omega_b, d_1_omega_b, mu_q_log_omega_b, mu_q_log_1_omega_b
#'
#' @noRd
update_omega_specific <- function(mu_q_gamma_b, p, S, c_0 = 1, d_0 = NULL,
                                    c_val = 1, bool_var_spec_prob = FALSE) {

  if (is.null(mu_q_gamma_b)) return(NULL)
  if (is.null(d_0)) d_0 <- p

  c_1_omega_b <- vector("list", S)
  d_1_omega_b <- vector("list", S)
  mu_q_log_omega_b <- vector("list", S)
  mu_q_log_1_omega_b <- vector("list", S)

  if (bool_var_spec_prob) {

    dig_b <- digamma(c_val * (c_0 + d_0 + 1) - 2 * c_val + 2)
    for (s in 1:S) {
      c_1_omega_b[[s]] <- c_val * (c_0 + mu_q_gamma_b[[s]]) - c_val + 1
      d_1_omega_b[[s]] <- c_val * (d_0 + 1 - mu_q_gamma_b[[s]]) - c_val + 1
      mu_q_log_omega_b[[s]] <- digamma(c_1_omega_b[[s]]) - dig_b
      mu_q_log_1_omega_b[[s]] <- digamma(d_1_omega_b[[s]]) - dig_b
    }

  } else {

    dig_b <- digamma(c_val * (c_0 + d_0 + p) - 2 * c_val + 2)
    for (s in 1:S) {
      cs_mu_q_gamma_b <- colSums(mu_q_gamma_b[[s]])
      c_1_omega_b[[s]] <- c_val * (c_0 + cs_mu_q_gamma_b) - c_val + 1
      d_1_omega_b[[s]] <- c_val * (d_0 + p - cs_mu_q_gamma_b) - c_val + 1
      mu_q_log_omega_b[[s]] <- digamma(c_1_omega_b[[s]]) - dig_b
      mu_q_log_1_omega_b[[s]] <- digamma(d_1_omega_b[[s]]) - dig_b
    }
  }

  create_named_list(c_1_omega_b, d_1_omega_b, mu_q_log_omega_b, mu_q_log_1_omega_b)
}

#' Update all omega parameters (shared + specific)
#'
#' @param mu_q_gamma_a p x L_f matrix of shared PPIs
#' @param mu_q_gamma_b List of S matrices, specific PPIs
#' @param p Number of variables
#' @param S Number of studies
#' @param c_0 Beta prior shape1
#' @param d_0 Beta prior shape2
#' @param c_val Temperature constant
#' @param bool_var_spec_prob Variable-specific probabilities flag
#' @return List with all updated omega parameters
#'
update_omega <- function(mu_q_gamma_a, mu_q_gamma_b = NULL,
                          p, S, c_0 = 1, d_0 = NULL,
                          c_val = 1, bool_var_spec_prob = FALSE) {

  if (is.null(d_0)) d_0 <- p

  res_shared <- update_omega_shared(mu_q_gamma_a, p, c_0, d_0,
                                     c_val, bool_var_spec_prob)

  res_specific <- update_omega_specific(mu_q_gamma_b, p, S, c_0, d_0,
                                         c_val, bool_var_spec_prob)

  create_named_list(
    c_1_omega_a = res_shared$c_1_omega_a,
    d_1_omega_a = res_shared$d_1_omega_a,
    mu_q_log_omega_a = res_shared$mu_q_log_omega_a,
    mu_q_log_1_omega_a = res_shared$mu_q_log_1_omega_a,
    c_1_omega_b = if (!is.null(res_specific)) res_specific$c_1_omega_b else NULL,
    d_1_omega_b = if (!is.null(res_specific)) res_specific$d_1_omega_b else NULL,
    mu_q_log_omega_b = if (!is.null(res_specific)) res_specific$mu_q_log_omega_b else NULL,
    mu_q_log_1_omega_b = if (!is.null(res_specific)) res_specific$mu_q_log_1_omega_b else NULL
  )
}
