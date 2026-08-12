test_that("multistart selection uses the highest valid endpoint", {
  summary_table <- data.frame(
    candidate_valid = c(TRUE, TRUE, TRUE, TRUE),
    eligible = c(FALSE, TRUE, TRUE, TRUE),
    final_elbo = c(100, 20, 30, 30)
  )
  expect_equal(
    multiFSYNC:::.select_multistart_summary(summary_table), 1L)
  summary_table$candidate_valid <- FALSE
  expect_true(is.na(
    multiFSYNC:::.select_multistart_summary(summary_table)))
})

mock_multistart_fit <- function(
    seed, elbo, eligible = TRUE,
    fitted = c(0, 0), rss = 1, ppi = 0.5) {
  list(
    fit_seed = seed,
    ELBO = elbo,
    elbo_objective_valid = TRUE,
    ELBO_decrease_count = 0L,
    elbo_t1_jitter_count = 0L,
    practical_converged = eligible,
    converged = eligible,
    slow_case = !eligible,
    i_iter = 10L,
    annealing_sweeps = 2L,
    t1_sweeps = 8L,
    expected_rss_sum = rss,
    stability_snapshot = list(
      fitted = fitted,
      rss = rss,
      ppi = ppi
    ),
    elbo_result = list(
      components = c(
        data_likelihood = elbo - 1,
        score_zeta = 1
      )
    )
  )
}

test_that("initialization controls preserve the legacy argument order", {
  expect_identical(
    tail(names(formals(bayesSYNC_multi)), 7L),
    c(
      "d_0", "convergence_rule", "lambda_orth", "practical_control",
      "initialization", "initialization_control", "continuation_state"
    )
  )
})

test_that("adaptive multistart stops after the first eligible fit", {
  dat <- make_practical_test_data(8401L)
  expect_warning(
    result <- bayesSYNC_multi_multistart(
      Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
      L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      K = 3L, anneal = c(1, 1.25, 2), n_g = 21L,
      maxit = 30L, convergence_rule = "practical",
      practical_control = small_practical_control(FALSE),
      bool_scale = FALSE, verbose = FALSE,
      start_seeds = c(8402L, 8403L, 8404L),
      retry_only = TRUE
    ),
    "not independently resolved"
  )
  expect_false(result$initialization_unresolved)
  expect_equal(result$attempted_starts, 1L)
  expect_equal(result$selected_seed, 8402L)
  expect_true(result$multistart_summary$eligible[1L])
  expect_equal(result$fit$fit_seed, 8402L)
  expect_true(result$objective_competition_unresolved)
  expect_true(result$output_competition_unresolved)
  expect_false(result$initialization_stable)
  expect_equal(result$best_basin_support, 1L)
  expect_equal(result$best_basin_output_support, 1L)
  expect_equal(
    result$fit$t1_sweeps, length(result$fit$ELBO))
  expect_identical(
    result$stage_history$decision,
    "first_eligible_fast_path"
  )
})

test_that("fixed multistart selects the highest valid endpoint ELBO", {
  dat <- make_practical_test_data(8501L)
  result <- bayesSYNC_multi_multistart(
    Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
    L_f = 0L, L_s = 0L,
    M_f = integer(), M_s = list(integer()),
    K = 3L, anneal = c(1, 1.25, 2), n_g = 21L,
    maxit = 30L, convergence_rule = "practical",
    practical_control = small_practical_control(FALSE),
    bool_scale = FALSE, verbose = FALSE,
    start_seeds = c(8502L, 8503L),
    retry_only = FALSE
  )
  candidates <- which(result$multistart_summary$candidate_valid)
  expected <- candidates[which.max(
    result$multistart_summary$final_elbo[candidates])]
  expect_equal(result$attempted_starts, 2L)
  expect_equal(result$selected_start, expected)
  expect_equal(
    result$selected_seed,
    result$multistart_summary$fit_seed[expected])
  expect_equal(sum(result$multistart_summary$selected), 1L)
  expect_true(all(
    result$multistart_summary$output_aligned[
      result$multistart_summary$objective_basin_member
    ]
  ))
  expect_equal(
    result$best_basin_output_support,
    result$best_basin_support
  )
  expect_false(result$output_competition_unresolved)
  expect_true(result$initialization_stable)
  expect_identical(
    result$stage_history$decision,
    "all_requested_starts_completed"
  )
})

test_that("adaptive multistart runs every retry after a failed first start", {
  fake_fit <- function(seed, eligible, elbo) {
    list(
      fit_seed = seed,
      ELBO = elbo,
      elbo_objective_valid = TRUE,
      ELBO_decrease_count = 0L,
      elbo_t1_jitter_count = 0L,
      practical_converged = eligible,
      converged = eligible,
      slow_case = !eligible,
      i_iter = 10L,
      annealing_sweeps = 2L,
      t1_sweeps = 8L,
      expected_rss_sum = 1,
      stability_snapshot = list(
        fitted = c(1, 2),
        rss = 1,
        ppi = 0.5
      )
    )
  }
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      switch(
        as.character(seed),
        "8702" = fake_fit(seed, FALSE, -20),
        "8703" = fake_fit(seed, TRUE, -10),
        "8704" = fake_fit(seed, TRUE, -5),
        stop("Unexpected seed")
      )
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = c(8702L, 8703L, 8704L),
      retry_only = TRUE
    ),
    "not independently resolved"
  )

  expect_equal(result$attempted_starts, 3L)
  expect_equal(result$selected_seed, 8704L)
  expect_equal(result$fit$fit_seed, 8704L)
  expect_false("fits" %in% names(result))
  expect_true(result$objective_competition_unresolved)
  expect_true(result$output_competition_unresolved)
})

test_that("near-best ELBO starts must also agree on scientific outputs", {
  fake_fit <- function(seed, fitted) {
    list(
      fit_seed = seed,
      ELBO = -100,
      elbo_objective_valid = TRUE,
      ELBO_decrease_count = 0L,
      elbo_t1_jitter_count = 0L,
      practical_converged = TRUE,
      converged = TRUE,
      slow_case = FALSE,
      i_iter = 10L,
      annealing_sweeps = 2L,
      t1_sweeps = 8L,
      expected_rss_sum = 1,
      stability_snapshot = list(
        fitted = fitted,
        rss = 1,
        ppi = 0.5
      )
    )
  }
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed == 8801L) {
        fake_fit(seed, c(0, 0))
      } else {
        fake_fit(seed, c(10, 10))
      }
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = c(8801L, 8802L),
      retry_only = FALSE
    ),
    "output competition=TRUE"
  )

  expect_equal(result$best_basin_support, 2L)
  expect_equal(result$best_basin_output_support, 1L)
  expect_false(result$objective_competition_unresolved)
  expect_true(result$output_competition_unresolved)
  expect_false(result$initialization_stable)
})

test_that("all valid slow starts return the best unfinished endpoint", {
  dat <- make_practical_test_data(8601L)
  expect_warning(
    result <- bayesSYNC_multi_multistart(
      Y = dat$Y, Z = dat$Z, time_obs = dat$time_obs,
      L_f = 0L, L_s = 0L,
      M_f = integer(), M_s = list(integer()),
      K = 3L, anneal = NULL, n_g = 21L,
      maxit = 20L, convergence_rule = "practical",
      practical_control = small_practical_control(TRUE),
      bool_scale = FALSE, verbose = FALSE,
      start_seeds = c(8602L, 8603L),
      retry_only = TRUE
    ),
    "selected endpoint unfinished=TRUE"
  )
  expect_false(result$initialization_unresolved)
  expect_false(is.null(result$fit))
  expect_equal(result$attempted_starts, 2L)
  expect_true(all(result$multistart_summary$slow_case))
  expect_true(all(result$multistart_summary$candidate_valid))
  expect_false(any(result$multistart_summary$eligible))
  expect_identical(
    result$selection_reason, "highest_elbo_valid_unfinished_endpoint")
  expect_true(result$selected_endpoint_unfinished)
  expect_false(result$production_ready)
  expect_true("selected_endpoint_unfinished" %in%
                result$unresolved_reasons)
  expect_true(is.finite(result$diagnostic_best_start))
})

test_that("multistart requires explicit valid seeds and a valid objective", {
  expect_error(
    bayesSYNC_multi_multistart(start_seeds = c(1L, 1L)),
    "distinct"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1L, convergence_rule = "parameters"),
    "requires convergence_rule"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1L, lambda_orth = 0.1),
    "requires lambda_orth"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1L, basin_rel_tol = -1),
    "basin_rel_tol"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1L, basin_abs_tol = -1),
    "basin_abs_tol"
  )
})

test_that("default staged gate does not accept apparent support before eight", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      elbo <- if (seed %in% c(1L, 2L)) {
        -100
      } else if (seed %in% c(6L, 7L)) {
        -10
      } else {
        -200 - seed
      }
      mock_multistart_fit(seed, elbo = elbo)
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:10,
    stage_sizes = c(3L, 5L, 10L)
  )

  expect_equal(result$min_exploration_starts, 8L)
  expect_equal(result$attempted_starts, 10L)
  expect_equal(result$selected_seed, 6L)
  expect_true(result$production_ready)
  expect_equal(
    result$stage_history$conditional_initialization_stable,
    c(TRUE, TRUE, TRUE)
  )
  expect_equal(
    result$stage_history$exploration_sufficient,
    c(FALSE, FALSE, TRUE)
  )
  expect_equal(
    result$stage_history$initialization_stable,
    c(FALSE, FALSE, TRUE)
  )
  expect_equal(
    result$stage_history$decision,
    c("expand", "expand", "stop_stable")
  )
  expect_true(all(grepl(
    "insufficient_exploration_starts",
    result$stage_history$unresolved_reasons[1:2],
    fixed = TRUE
  )))
})

test_that("an explicit diagnostic gate can stop after three stable starts", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed > 3L) stop("A stable first stage must not run later seeds.")
      mock_multistart_fit(
        seed,
        elbo = c(`1` = -10, `2` = -10, `3` = -20)[as.character(seed)]
      )
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:10,
    stage_sizes = c(3L, 5L, 10L),
    min_exploration_starts = 3L
  )

  expect_equal(result$min_exploration_starts, 3L)
  expect_equal(result$attempted_starts, 3L)
  expect_equal(result$attempted_start_seeds, 1:3)
  expect_equal(result$completed_stage, 1L)
  expect_identical(result$multistart_stop_reason, "stable_after_3_starts")
  expect_true(result$production_ready)
  expect_identical(result$initialization_status, "stable")
  expect_length(result$unresolved_reasons, 0L)
  expect_equal(result$best_basin_support, 2L)
  expect_equal(result$best_basin_output_support, 2L)
  expect_identical(result$stage_history$decision, "stop_stable")
  expect_equal(unique(result$multistart_summary$stage_id), 1L)
  expect_equal(nrow(result$multistart_elbo_components), 6L)
})

test_that("staged multistart expands from three to five when needed", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed > 5L) stop("The stable second stage must not run later seeds.")
      elbo <- c(
        `1` = -10, `2` = -20, `3` = -30,
        `4` = -10, `5` = -25
      )[as.character(seed)]
      mock_multistart_fit(seed, elbo = elbo)
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:10,
    stage_sizes = c(3L, 5L, 10L),
    min_exploration_starts = 3L
  )

  expect_equal(result$attempted_starts, 5L)
  expect_equal(result$completed_stage, 2L)
  expect_identical(result$multistart_stop_reason, "stable_after_5_starts")
  expect_true(result$production_ready)
  expect_equal(result$selected_seed, 1L)
  expect_equal(result$best_basin_support, 2L)
  expect_equal(result$stage_history$decision, c("expand", "stop_stable"))
  expect_equal(result$stage_history$attempted_starts, c(3L, 5L))
  expect_equal(result$multistart_summary$stage_id, c(1L, 1L, 1L, 2L, 2L))
})

test_that("an eight-to-ten staged schedule is valid", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed > 8L) stop("A stable eight-start stage must stop.")
      mock_multistart_fit(seed, elbo = -10)
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:10,
    stage_sizes = c(8L, 10L)
  )

  expect_equal(result$attempted_starts, 8L)
  expect_equal(result$completed_stage, 1L)
  expect_equal(result$min_exploration_starts, 8L)
  expect_true(result$exploration_sufficient)
  expect_true(result$production_ready)
  expect_identical(result$multistart_stop_reason, "stable_after_8_starts")
})

test_that("staged multistart can resolve only at its ten-start budget", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      elbo <- if (seed %in% c(1L, 9L)) {
        -10
      } else {
        -20 - seed
      }
      mock_multistart_fit(seed, elbo = elbo)
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:10,
    stage_sizes = c(3L, 5L, 10L),
    keep_fits = TRUE
  )

  expect_equal(result$attempted_starts, 10L)
  expect_equal(result$completed_stage, 3L)
  expect_identical(result$multistart_stop_reason, "stable_after_10_starts")
  expect_true(result$production_ready)
  expect_equal(result$selected_seed, 1L)
  expect_equal(result$best_basin_support, 2L)
  expect_equal(
    result$stage_history$decision,
    c("expand", "expand", "stop_stable")
  )
  expect_length(result$fits, 10L)
  expect_equal(vapply(result$fits, `[[`, integer(1), "fit_seed"), 1:10)
})

test_that("a higher unfinished endpoint is selected but never production-ready", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed == 3L) {
        return(mock_multistart_fit(seed, elbo = 0, eligible = FALSE))
      }
      mock_multistart_fit(seed, elbo = -10)
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L)
    ),
    "selected endpoint unfinished=TRUE"
  )

  expect_equal(result$attempted_starts, 5L)
  expect_false(result$production_ready)
  expect_identical(
    result$multistart_stop_reason,
    "budget_exhausted_unresolved"
  )
  expect_equal(result$selected_seed, 3L)
  expect_true(result$selected_endpoint_unfinished)
  expect_false(result$higher_ineligible_competitor)
  expect_length(result$higher_ineligible_start_ids, 0L)
  expect_length(result$higher_ineligible_seeds, 0L)
  expect_true("selected_endpoint_unfinished" %in%
                result$unresolved_reasons)
  expect_false(any(
    result$multistart_summary$higher_than_selected_unfinished))
  expect_equal(
    result$stage_history$decision,
    c("expand", "stop_budget_exhausted")
  )
})

test_that("a higher unreliable finite fit blocks production but is not selected", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      fit <- mock_multistart_fit(seed, elbo = if (seed == 3L) 0 else -10)
      if (seed == 3L) {
        fit$elbo_objective_valid <- FALSE
      }
      fit
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(start_seeds = 1:3),
    "higher unreliable competitor=TRUE"
  )

  expect_equal(result$selected_seed, 1L)
  expect_false(result$production_ready)
  expect_true(result$conditional_initialization_stable == FALSE)
  expect_false(result$higher_ineligible_competitor)
  expect_true(result$higher_unreliable_competitor)
  expect_equal(result$higher_unreliable_start_ids, 3L)
  expect_equal(result$higher_unreliable_seeds, 3L)
  expect_true(
    "higher_unreliable_competitor" %in% result$unresolved_reasons
  )
  expect_false(
    result$multistart_summary$higher_than_selected_unfinished[3L]
  )
  expect_true(
    result$multistart_summary$higher_than_selected_unreliable[3L]
  )
  expect_identical(
    result$stage_history$higher_unreliable_start_ids,
    "3"
  )
  expect_identical(result$stage_history$higher_unreliable_seeds, "3")
})

test_that("a later supported basin can overtake an unfinished competitor", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed == 3L) {
        return(mock_multistart_fit(seed, elbo = 0, eligible = FALSE))
      }
      elbo <- if (seed %in% c(4L, 5L)) 1 else -10
      mock_multistart_fit(seed, elbo = elbo)
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:5,
    stage_sizes = c(3L, 5L)
  )

  expect_equal(result$attempted_starts, 5L)
  expect_true(result$production_ready)
  expect_identical(result$multistart_stop_reason, "stable_after_5_starts")
  expect_equal(result$selected_seed, 4L)
  expect_false(result$higher_ineligible_competitor)
  expect_length(result$higher_ineligible_seeds, 0L)
  expect_equal(
    result$stage_history$higher_ineligible_seeds,
    c("", "")
  )
  expect_equal(
    result$stage_history$decision,
    c("expand", "stop_stable")
  )
})

test_that("staged multistart returns a valid unfinished endpoint at budget", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      mock_multistart_fit(seed, elbo = -seed, eligible = FALSE)
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L)
    ),
    "selected endpoint unfinished=TRUE"
  )

  expect_equal(result$attempted_starts, 5L)
  expect_false(is.null(result$fit))
  expect_false(result$initialization_unresolved)
  expect_true(result$selected_endpoint_unfinished)
  expect_false(result$production_ready)
  expect_identical(
    result$initialization_status,
    "best_valid_unfinished"
  )
  expect_true("selected_endpoint_unfinished" %in%
                result$unresolved_reasons)
  expect_identical(
    result$multistart_stop_reason,
    "budget_exhausted_unresolved"
  )
  expect_equal(
    result$stage_history$decision,
    c("expand", "stop_budget_exhausted")
  )
})

test_that("staged multistart returns no fit when every endpoint is invalid", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      fit <- mock_multistart_fit(seed, elbo = -seed, eligible = FALSE)
      fit$elbo_objective_valid <- FALSE
      fit
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L)
    ),
    "without a valid ordinary-ELBO endpoint"
  )
  expect_null(result$fit)
  expect_true(result$initialization_unresolved)
  expect_false(result$selected_endpoint_unfinished)
  expect_identical(result$initialization_status, "no_valid_endpoint")
  expect_identical(result$unresolved_reasons, "no_valid_endpoint")
  expect_identical(
    result$multistart_stop_reason,
    "budget_exhausted_no_valid_endpoint"
  )
})

test_that("kept fits preserve failed start positions", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      if (seed == 2L) stop("mock failure")
      mock_multistart_fit(seed, elbo = -seed)
    },
    .package = "multiFSYNC"
  )

  expect_warning(
    result <- bayesSYNC_multi_multistart(
      start_seeds = 1:3,
      keep_fits = TRUE
    ),
    "not independently resolved"
  )

  expect_length(result$fits, 3L)
  expect_equal(result$fits[[1L]]$fit_seed, 1L)
  expect_null(result$fits[[2L]])
  expect_equal(result$fits[[3L]]$fit_seed, 3L)
  expect_false(result$multistart_summary$success[2L])
  expect_match(result$multistart_summary$error_text[2L], "mock failure")
})

test_that("staged multistart validates its deterministic schedule", {
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(2L, 5L)
    ),
    "start at 3"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 3L, 5L)
    ),
    "strictly increasing"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 4L)
    ),
    "must equal length"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      retry_only = TRUE,
      stage_sizes = c(3L, 5L)
    ),
    "cannot be combined"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L),
      min_exploration_starts = 0L
    ),
    "positive integer"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L),
      min_exploration_starts = 2.5
    ),
    "positive integer"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:5,
      stage_sizes = c(3L, 5L),
      min_exploration_starts = 6L
    ),
    "no larger"
  )
})

test_that("mixed initialization schedules count independent sources", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed, initialization) {
      fit <- mock_multistart_fit(seed, elbo = -10)
      fit$initialization <- initialization
      fit$initialization_independence_id <- if (
          identical(initialization, "residual_fpca")) {
        "residual_fpca_deterministic"
      } else {
        paste0("random_seed_", seed)
      }
      fit$initialization_diagnostics <- list(
        complete = TRUE,
        fallback_count = 0L,
        explained_fraction = if (
          identical(initialization, "residual_fpca")) 0.8 else NA_real_,
        elapsed_seconds = 0.01
      )
      fit
    },
    .package = "multiFSYNC"
  )

  result <- bayesSYNC_multi_multistart(
    start_seeds = 1:3,
    start_initializations =
      c("residual_fpca", "residual_fpca", "random")
  )

  expect_equal(
    result$attempted_start_initializations,
    c("residual_fpca", "residual_fpca", "random")
  )
  expect_equal(sum(result$multistart_summary$objective_basin_member), 3L)
  expect_equal(result$best_basin_support, 2L)
  expect_equal(result$best_basin_output_support, 2L)
  expect_true(result$production_ready)
  expect_equal(
    length(unique(
      result$multistart_summary$initialization_independence_id)),
    2L
  )
  expect_equal(
    unique(result$multistart_elbo_components$initialization),
    c("residual_fpca", "random")
  )
})

test_that("multistart validates initialization schedules", {
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:2,
      start_initializations = c("random", "unknown")
    ),
    "only 'random' or 'residual_fpca'"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:3,
      start_initializations = c("random", "residual_fpca")
    ),
    "length one or length"
  )
  expect_error(
    bayesSYNC_multi_multistart(
      start_seeds = 1:2,
      initialization = "random",
      start_initializations = "residual_fpca"
    ),
    "either start_initializations or initialization"
  )
})

test_that("an absolute ELBO floor avoids cancellation-driven basin splitting", {
  testthat::local_mocked_bindings(
    bayesSYNC_multi = function(..., seed) {
      mock_multistart_fit(
        seed,
        elbo = if (seed == 1L) 20 else 19.98
      )
    },
    .package = "multiFSYNC"
  )

  resolved <- bayesSYNC_multi_multistart(
    start_seeds = 1:2
  )
  expect_equal(resolved$best_basin_support, 2L)
  expect_true(resolved$production_ready)
  expect_equal(resolved$basin_abs_tol, 0.05)
  expect_equal(
    max(resolved$multistart_summary$basin_absolute_gap),
    0.02,
    tolerance = 1e-12
  )

  expect_warning(
    relative_only <- bayesSYNC_multi_multistart(
      start_seeds = 1:2,
      basin_abs_tol = 0
    ),
    "not independently resolved"
  )
  expect_equal(relative_only$best_basin_support, 1L)
  expect_false(relative_only$production_ready)
})
