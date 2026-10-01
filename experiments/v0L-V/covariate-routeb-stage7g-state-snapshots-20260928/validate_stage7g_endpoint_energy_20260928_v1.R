#!/usr/bin/env Rscript

# Independent endpoint check from the already frozen Stage 7F fit objects.
# It compares direct grid reconstruction and covariance-integrated function
# energy with the Stage 7G all-sweep machine table. No new fit or truth read.
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
root <- "stage7g_20260928_v1"
table <- read.csv(file.path(root, "analysis_v1",
                            "STATE_METRICS_ALL_SWEEPS.csv"))
ids <- c("covrb7c_01__C1__start06", "covrb7c_02__C1__start11",
         "covrb7c_03__C1__start06")
source_root <- file.path("..", "covariate-routeb-stage7f-dense-replay-20260928",
                         "stage7f_20260928_v1", "fits")
trap_weights <- function(grid) {
  delta <- diff(grid)
  c(delta[[1L]] / 2,
    (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}
rows <- lapply(ids, function(id) {
  fit <- readRDS(file.path(source_root, id, "fit.rds"))
  target <- table[table$fit_id == id & table$iteration == 499L &
                    table$block == "specific" & table$study == 2L, ]
  stopifnot(nrow(target) == 1L)
  weights <- trap_weights(fit$time_g)
  basis_gram <- crossprod(fit$C_g,
                          sweep(fit$C_g, 1L, weights, "*"))
  coefficients <- fit$mu_q_nu_psi[[2L]][[1L]]
  coefficient_covariances <- fit$Sigma_q_nu_psi[[2L]][[1L]]
  functions <- fit$C_g %*% coefficients
  function_mean_sq <- sum(sweep(functions^2, 1L, weights, "*"))
  function_cov_sq <- sum(vapply(coefficient_covariances, function(value)
    sum(diag(basis_gram %*% value)), numeric(1L)))
  scores <- fit$mu_q_xi[[2L]][[1L]]
  subject_curves <- scores %*% t(functions)
  subject_grid_l2_sq <- rowSums(sweep(subject_curves^2, 2L,
                                     weights, "*"))
  loading <- fit$mu_q_b_specific[[2L]][, 1L]
  direct_mean_signal_rms <- sqrt(mean(loading^2) *
                                   mean(subject_grid_l2_sq))
  data.frame(
    fit_id = id,
    direct_mean_signal_rms = direct_mean_signal_rms,
    table_mean_signal_rms =
      target$mean_signal_l2_rms_over_variables,
    direct_function_mean_l2 = sqrt(function_mean_sq),
    table_function_mean_l2 = target$function_mean_l2,
    direct_function_second_l2 = sqrt(function_mean_sq + function_cov_sq),
    table_function_second_l2 = target$function_second_l2,
    function_covariance_l2_sq = function_cov_sq,
    max_coefficient_covariance_diagonal = max(unlist(lapply(
      coefficient_covariances, diag))),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
stopifnot(all(abs(result$direct_mean_signal_rms -
                    result$table_mean_signal_rms) < 1e-10),
          all(abs(result$direct_function_mean_l2 -
                    result$table_function_mean_l2) < 1e-10),
          all(abs(result$direct_function_second_l2 -
                    result$table_function_second_l2) < 1e-6))
out <- "posthoc_v1/ENDPOINT_ENERGY_INDEPENDENT_CHECK.csv"
if (file.exists(out)) stop("Refusing overwrite: ", out)
write.csv(result, out, row.names = FALSE)
print(result, digits = 7L, row.names = FALSE)
cat("STAGE7G_ENDPOINT_ENERGY_INDEPENDENT_CHECK_PASS\n")
