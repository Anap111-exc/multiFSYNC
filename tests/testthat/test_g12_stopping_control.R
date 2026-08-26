test_that("G12 stopping control freezes the selected fixed-400 profile", {
  control <- g12_stopping_control()

  expect_identical(
    names(control), names(multiFSYNC:::.default_practical_control())
  )
  expect_identical(control$min_t1, 396L)
  expect_identical(control$window, 20L)
  expect_identical(control$long_window, 60L)
  expect_identical(control$consecutive, 5L)
  expect_identical(control$max_t1, 400L)
  expect_identical(control$ppi_gate, "quantile_factor")
  expect_equal(control$ppi_quantile, 0.95)
  expect_equal(control$ppi_quantile_max_abs, 1e-2)
  expect_equal(control$factor_ppi_max_abs, 1e-2)
  expect_equal(control$long_ppi_quantile_max_abs, 1e-2)
  expect_equal(control$long_factor_ppi_max_abs, 1e-2)
  expect_length(control$checkpoints, 0L)
  expect_identical(
    multiFSYNC:::.validate_practical_control(control), control
  )
})

test_that("G12 profile cannot practically stop before sweep 400", {
  control <- g12_stopping_control()
  snapshot <- list(
    fitted = c(1, 2), rss = c(3, 4),
    ppi = numeric(), factor_ppi = numeric()
  )
  elbo <- -100 + seq_len(400L) * 1e-6
  diagnostics <- multiFSYNC:::.elbo_history_diagnostics(elbo)
  passes <- vapply(391:400, function(sweep) {
    multiFSYNC:::.practical_window_diagnostic(
      sweep = sweep,
      ELBO = elbo,
      ELBO_diagnostics = diagnostics,
      current_snapshot = snapshot,
      reference_snapshot = snapshot,
      long_reference_snapshot = snapshot,
      control = control,
      objective_size = 1
    )$pass
  }, logical(1))

  expect_false(any(passes[1:5]))
  expect_true(all(passes[6:10]))
  streak <- 0L
  first_five <- NA_integer_
  for (index in seq_along(passes)) {
    streak <- if (passes[[index]]) streak + 1L else 0L
    if (streak >= control$consecutive) {
      first_five <- 390L + index
      break
    }
  }
  expect_identical(first_five, 400L)
})
