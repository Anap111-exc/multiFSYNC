# Residual-FPCA initialization tests.
#
# These tests deliberately use small deterministic signals.  They verify the
# iteration-zero construction separately from the short CAVI integration so a
# failure identifies either the initializer or the unchanged update sequence.

.make_initialization_quality_data <- function(
    S = 1L, n = 6L, p = 3L, d = 0L,
    shared = FALSE, specific = FALSE,
    noise_sd = 1e-3, seed = 9101L) {
  set.seed(seed)
  stopifnot(d %in% c(0L, 1L))

  time_grid <- seq(0.05, 0.95, length.out = 9L)
  time_obs <- lapply(seq_len(S), function(s) {
    lapply(seq_len(n), function(i) time_grid)
  })
  z <- seq(-1, 1, length.out = n)
  Z <- if (d == 1L) {
    lapply(seq_len(S), function(s) matrix(z, ncol = 1L))
  } else {
    NULL
  }

  subject_design <- if (d == 1L) cbind(1, z) else matrix(1, n, 1L)
  subject_basis <- qr.Q(qr(subject_design), complete = TRUE)
  shared_score <- subject_basis[, ncol(subject_design) + 1L] * sqrt(n)
  specific_score <- subject_basis[, ncol(subject_design) + 2L] * sqrt(n)
  shared_loading <- seq(1, 0.45, length.out = p)
  specific_loading <- lapply(seq_len(S), function(s) {
    direction <- if (s %% 2L) {
      seq(-0.35, 0.95, length.out = p)
    } else {
      seq(0.85, -0.45, length.out = p)
    }
    direction / sqrt(sum(direction^2))
  })

  Y <- lapply(seq_len(S), function(s) {
    lapply(seq_len(n), function(i) {
      lapply(seq_len(p), function(j) {
        tt <- time_obs[[s]][[i]]
        value <- 0.12 * s + 0.04 * j + 0.03 * j * tt
        if (d == 1L) {
          value <- value +
            z[i] * (0.08 + 0.015 * j) * cos(pi * tt)
        }
        if (shared) {
          value <- value +
            0.75 * shared_score[i] * shared_loading[j] *
            sin(2 * pi * tt)
        }
        if (specific) {
          value <- value +
            0.55 * specific_score[i] * specific_loading[[s]][j] *
            cos(2 * pi * tt)
        }
        value + stats::rnorm(length(tt), sd = noise_sd)
      })
    })
  })

  # Intentionally no true_params element: the initializer must use observed
  # working-scale data only.
  list(Y = Y, Z = Z, time_obs = time_obs)
}

.run_initialization_quality_fit <- function(
    dat, L_f, L_s, bool_var_spec_prob = FALSE,
    seed = 9201L, maxit = 4L,
    initialization = "residual_fpca",
    initialization_control = list(
      grid_size = 21L, perturb_sd = 0, rank_tol = 1e-8
    ),
    omit_initialization = FALSE) {
  S <- length(dat$Y)
  args <- list(
    Y = dat$Y,
    Z = dat$Z,
    time_obs = dat$time_obs,
    L_f = L_f,
    L_s = L_s,
    M_f = if (L_f > 0L) rep(1L, L_f) else integer(),
    M_s = lapply(seq_len(S), function(s) {
      if (L_s > 0L) rep(1L, L_s) else integer()
    }),
    K = 2L,
    n_g = 21L,
    anneal = NULL,
    maxit = maxit,
    n_cpus = 1L,
    tol_abs = 0,
    tol_rel = 0,
    convergence_rule = "elbo",
    bool_scale = FALSE,
    bool_var_spec_prob = bool_var_spec_prob,
    verbose = FALSE,
    seed = seed
  )
  if (!omit_initialization) {
    args$initialization <- initialization
    args$initialization_control <- initialization_control
  }
  suppressWarnings(do.call(bayesSYNC_multi, args))
}

.construct_residual_fpca_quality_start <- function(
    dat, L_f = 1L, L_s = 0L, seed = 9301L,
    perturb_sd = 0.05) {
  S <- length(dat$Y)
  n_s <- vapply(dat$Y, length, integer(1))
  p <- length(dat$Y[[1L]][[1L]])
  d <- if (is.null(dat$Z)) 0L else ncol(dat$Z[[1L]])
  time_flat <- unlist(dat$time_obs, recursive = FALSE)
  grid <- multiFSYNC:::get_grid_objects(
    time_flat, K = 2L, n_g = 21L, format_univ = TRUE
  )

  set.seed(seed)
  multiFSYNC:::.residual_fpca_initialization(
    Y = dat$Y,
    Z = dat$Z,
    time_obs = dat$time_obs,
    C_g = grid$C_g,
    time_g = grid$time_g,
    S = S,
    n_s = n_s,
    p = p,
    d = d,
    L_f = L_f,
    L_s = L_s,
    M_f = if (L_f > 0L) rep(1L, L_f) else integer(),
    M_s = lapply(seq_len(S), function(s) {
      if (L_s > 0L) rep(1L, L_s) else integer()
    }),
    control = multiFSYNC:::.validate_initialization_control(list(
      grid_size = 21L,
      perturb_sd = perturb_sd,
      rank_tol = 1e-8
    ))
  )
}

.without_initialization_elapsed <- function(x) {
  x$diagnostics$elapsed_seconds <- NULL
  x
}

test_that("initialization controls are resolved and strictly validated", {
  defaults <- multiFSYNC:::.validate_initialization_control()
  expect_named(defaults, c("grid_size", "perturb_sd", "rank_tol"))
  expect_true(defaults$grid_size >= 3L)
  expect_true(defaults$perturb_sd >= 0)
  expect_true(defaults$rank_tol > 0)

  changed <- multiFSYNC:::.validate_initialization_control(list(
    grid_size = 11L, perturb_sd = 0
  ))
  expect_identical(changed$grid_size, 11L)
  expect_identical(changed$perturb_sd, 0)
  expect_identical(changed$rank_tol, defaults$rank_tol)

  expect_error(
    multiFSYNC:::.validate_initialization_control(1),
    "fully named list"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(11L)),
    "fully named list"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(extra = 1)),
    "Unknown initialization_control"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(grid_size = 3.5)),
    "integer"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(perturb_sd = -0.1)),
    "non-negative"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(rank_tol = 0)),
    "strictly between 0 and 1"
  )
  expect_error(
    multiFSYNC:::.validate_initialization_control(list(rank_tol = 1)),
    "strictly between 0 and 1"
  )
})

test_that("omitted initialization is exactly the seeded random legacy path", {
  dat <- .make_initialization_quality_data(
    S = 1L, shared = TRUE, seed = 9401L
  )
  implicit <- .run_initialization_quality_fit(
    dat, L_f = 1L, L_s = 0L,
    seed = 9402L, maxit = 2L,
    omit_initialization = TRUE
  )
  explicit <- .run_initialization_quality_fit(
    dat, L_f = 1L, L_s = 0L,
    seed = 9402L, maxit = 2L,
    initialization = "random"
  )

  expect_identical(implicit$initialization, "random")
  expect_identical(explicit$initialization, "random")
  expect_identical(implicit$initialization_independence_id,
                   explicit$initialization_independence_id)
  expect_identical(implicit$ELBO, explicit$ELBO)
  expect_identical(implicit$ELBO_components, explicit$ELBO_components)
  expect_identical(implicit$mu_q_nu_mu, explicit$mu_q_nu_mu)
  expect_identical(implicit$mu_q_nu_phi, explicit$mu_q_nu_phi)
  expect_identical(implicit$mu_q_zeta, explicit$mu_q_zeta)
  expect_identical(implicit$mu_q_a, explicit$mu_q_a)
  expect_identical(implicit$mu_q_gamma_a, explicit$mu_q_gamma_a)
  expect_identical(implicit$expected_rss, explicit$expected_rss)
})

test_that("residual-FPCA starts are seeded reproducibly and perturbable", {
  dat <- .make_initialization_quality_data(
    S = 1L, shared = TRUE, seed = 9501L
  )
  first <- .construct_residual_fpca_quality_start(
    dat, seed = 9502L, perturb_sd = 0.12
  )
  repeat_first <- .construct_residual_fpca_quality_start(
    dat, seed = 9502L, perturb_sd = 0.12
  )
  second_seed <- .construct_residual_fpca_quality_start(
    dat, seed = 9503L, perturb_sd = 0.12
  )

  expect_identical(
    .without_initialization_elapsed(first),
    .without_initialization_elapsed(repeat_first)
  )
  expect_true(first$diagnostics$complete)
  expect_true(second_seed$diagnostics$complete)
  expect_gt(
    max(abs(
      first$diagnostics$shared_loadings -
        second_seed$diagnostics$shared_loadings
    )),
    1e-8
  )
})

test_that("warm starts cover small model structures and preserve valid ELBOs", {
  cases <- list(
    pure_mean = list(
      S = 1L, d = 0L, shared = FALSE, specific = FALSE,
      L_f = 0L, L_s = 0L, variable_omega = FALSE, seed = 9601L
    ),
    shared_only = list(
      S = 1L, d = 0L, shared = TRUE, specific = FALSE,
      L_f = 1L, L_s = 0L, variable_omega = FALSE, seed = 9602L
    ),
    shared_specific_covariate_factor_omega = list(
      S = 2L, d = 1L, shared = TRUE, specific = TRUE,
      L_f = 1L, L_s = 1L, variable_omega = FALSE, seed = 9603L
    ),
    shared_specific_covariate_variable_omega = list(
      S = 2L, d = 1L, shared = TRUE, specific = TRUE,
      L_f = 1L, L_s = 1L, variable_omega = TRUE, seed = 9604L
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    dat <- .make_initialization_quality_data(
      S = case$S,
      d = case$d,
      shared = case$shared,
      specific = case$specific,
      seed = case$seed
    )
    expect_false("true_params" %in% names(dat), info = case_name)

    fit <- .run_initialization_quality_fit(
      dat,
      L_f = case$L_f,
      L_s = case$L_s,
      bool_var_spec_prob = case$variable_omega,
      seed = case$seed + 100L,
      maxit = 4L
    )
    diagnostics <- fit$initialization_diagnostics

    expect_identical(fit$initialization, "residual_fpca",
                     info = case_name)
    expect_identical(diagnostics$requested_method, "residual_fpca",
                     info = case_name)
    expect_identical(diagnostics$used_method, "residual_fpca",
                     info = case_name)
    expect_true(diagnostics$complete, info = case_name)
    expect_identical(diagnostics$fallback_count, 0L, info = case_name)
    expect_true(is.finite(diagnostics$explained_fraction),
                info = case_name)
    expect_true(all(is.finite(fit$ELBO)), info = case_name)
    expect_true(all(vapply(
      fit$ELBO_components, function(x) all(is.finite(x)), logical(1)
    )), info = case_name)
    expect_identical(fit$ELBO_decrease_count, 0L, info = case_name)

    if (length(fit$ELBO) > 1L) {
      elbo_difference <- diff(fit$ELBO)
      scaled_tolerance <- fit$ELBO_diagnostics$tolerance *
        (1 + pmax(
          abs(head(fit$ELBO, -1L)),
          abs(tail(fit$ELBO, -1L))
        ))
      expect_true(
        all(elbo_difference >= -scaled_tolerance),
        info = case_name
      )
    }

    expected_hierarchy <- if (case$variable_omega) {
      "variable_factor"
    } else {
      "factor"
    }
    expect_identical(
      fit$elbo_result$omega_hierarchy,
      expected_hierarchy,
      info = case_name
    )

    if (case$L_f > 0L) {
      expect_equal(
        unlist(diagnostics$shared_score_means, use.names = FALSE),
        rep(0, sum(fit$M_f)),
        tolerance = 1e-12,
        info = case_name
      )
      expect_equal(
        unlist(
          diagnostics$shared_score_second_moments,
          use.names = FALSE
        ),
        rep(1, sum(fit$M_f)),
        tolerance = 1e-12,
        info = case_name
      )
    }
    if (case$L_s > 0L) {
      expect_equal(
        unlist(diagnostics$specific_score_means, use.names = FALSE),
        rep(0, sum(unlist(fit$M_s, use.names = FALSE))),
        tolerance = 1e-12,
        info = case_name
      )
      expect_equal(
        unlist(
          diagnostics$specific_score_second_moments,
          use.names = FALSE
        ),
        rep(1, sum(unlist(fit$M_s, use.names = FALSE))),
        tolerance = 1e-12,
        info = case_name
      )
    }
  }
})

test_that("a degenerate residual returns an explicit seeded-fallback reason", {
  residual <- array(0, dim = c(4L, 7L, 2L))
  tt <- seq(0, 1, length.out = 7L)
  result <- multiFSYNC:::.initialization_one_factor(
    residual = residual,
    C_grid = cbind(1, tt, tt^2),
    weights = rep(1 / 7, 7L),
    M = 1L,
    control = multiFSYNC:::.validate_initialization_control(list(
      perturb_sd = 0
    )),
    label = "shared[1]"
  )

  expect_false(result$success)
  expect_identical(
    result$reason,
    "zero_or_nonfinite_residual_energy"
  )
  expect_identical(result$label, "shared[1]")
  expect_identical(result$explained_fraction, 0)
  expect_null(result$residual)
  expect_error(
    multiFSYNC:::.initialization_one_factor(
      residual = residual,
      C_grid = cbind(1, tt, tt^2),
      weights = rep(1 / 7, 7L),
      M = 0L,
      control = multiFSYNC:::.validate_initialization_control(),
      label = "shared[1]"
    ),
    "positive integer"
  )
})

test_that("warm starts without an effective perturb share one source ID", {
  no_factor_dat <- .make_initialization_quality_data(
    S = 1L, p = 2L, shared = FALSE, seed = 9701L
  )
  no_factor_one <- .run_initialization_quality_fit(
    no_factor_dat, L_f = 0L, L_s = 0L,
    seed = 9702L, maxit = 1L,
    initialization_control = list(
      grid_size = 21L, perturb_sd = 0.2, rank_tol = 1e-8
    )
  )
  no_factor_two <- .run_initialization_quality_fit(
    no_factor_dat, L_f = 0L, L_s = 0L,
    seed = 9703L, maxit = 1L,
    initialization_control = list(
      grid_size = 21L, perturb_sd = 0.2, rank_tol = 1e-8
    )
  )
  expect_identical(
    no_factor_one$initialization_independence_id,
    "residual_fpca_deterministic"
  )
  expect_identical(
    no_factor_one$initialization_independence_id,
    no_factor_two$initialization_independence_id
  )

  one_variable_dat <- .make_initialization_quality_data(
    S = 1L, p = 1L, shared = TRUE, seed = 9704L
  )
  one_variable_one <- .run_initialization_quality_fit(
    one_variable_dat, L_f = 1L, L_s = 0L,
    seed = 9705L, maxit = 1L,
    initialization_control = list(
      grid_size = 21L, perturb_sd = 0.2, rank_tol = 1e-8
    )
  )
  one_variable_two <- .run_initialization_quality_fit(
    one_variable_dat, L_f = 1L, L_s = 0L,
    seed = 9706L, maxit = 1L,
    initialization_control = list(
      grid_size = 21L, perturb_sd = 0.2, rank_tol = 1e-8
    )
  )
  expect_true(one_variable_one$initialization_diagnostics$complete)
  expect_identical(
    one_variable_one$initialization_independence_id,
    "residual_fpca_deterministic"
  )
  expect_identical(
    one_variable_one$initialization_independence_id,
    one_variable_two$initialization_independence_id
  )
})

test_that("factor extraction is equivariant to response scale", {
  tt <- seq(0, 1, length.out = 21L)
  C_grid <- cbind(1, tt, sin(2 * pi * tt))
  weights <- multiFSYNC:::.trap_weights(tt)
  scores <- seq(-1, 1, length.out = 6L)
  loading <- c(1, 0.7, -0.25)
  residual <- array(0, dim = c(length(scores), length(tt), length(loading)))
  for (i in seq_along(scores)) {
    residual[i, , ] <- tcrossprod(
      scores[i] * sin(2 * pi * tt), loading
    )
  }
  control <- multiFSYNC:::.validate_initialization_control(list(
    perturb_sd = 0,
    rank_tol = 1e-8
  ))
  unit_scale <- multiFSYNC:::.initialization_one_factor(
    residual, C_grid, weights, M = 1L,
    control = control, label = "unit"
  )
  tiny_scale <- multiFSYNC:::.initialization_one_factor(
    residual * 1e-5, C_grid, weights, M = 1L,
    control = control, label = "tiny"
  )

  expect_true(unit_scale$success)
  expect_true(tiny_scale$success)
  expect_equal(
    tiny_scale$loading, unit_scale$loading,
    tolerance = 1e-10
  )
  expect_equal(
    tiny_scale$explained_fraction,
    unit_scale$explained_fraction,
    tolerance = 1e-10
  )
  expect_equal(
    tiny_scale$fitted / 1e-5,
    unit_scale$fitted,
    tolerance = 1e-9
  )
})
