#!/usr/bin/env Rscript

# Endpoint factor-product second moment under the fitted mean-field factors.
# Uses frozen Stage 7F fit objects, whose raw endpoints exactly match Stage 7G.
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
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
  weights <- trap_weights(fit$time_g)
  basis_gram <- crossprod(fit$C_g,
                          sweep(fit$C_g, 1L, weights, "*"))
  coeff <- fit$mu_q_nu_psi[[2L]][[1L]]
  coeff_cov <- fit$Sigma_q_nu_psi[[2L]][[1L]]
  function_mean <- fit$C_g %*% coeff
  function_gram <- crossprod(function_mean,
                             sweep(function_mean, 1L, weights, "*"))
  function_gram_second <- function_gram
  diag(function_gram_second) <- diag(function_gram_second) +
    vapply(coeff_cov, function(value)
      sum(diag(basis_gram %*% value)), numeric(1L))
  score_mean <- fit$mu_q_xi[[2L]][[1L]]
  score_cov <- fit$Sigma_q_xi[[2L]][[1L]]
  subject_second <- vapply(seq_len(nrow(score_mean)), function(i) {
    as.numeric(crossprod(score_mean[i, ],
                         function_gram_second %*% score_mean[i, ])) +
      sum(diag(score_cov[[i]] %*% function_gram_second))
  }, numeric(1L))
  loading_second <- fit$term_b_specific[[2L]][, 1L]
  stopifnot(all(is.finite(subject_second)),
            all(subject_second >= -1e-6),
            all(loading_second >= 0))
  data.frame(
    fit_id = id,
    mean_loading_second = mean(loading_second),
    integrated_score_function_second = mean(subject_second),
    full_factor_second_rms = sqrt(mean(loading_second) *
                                    mean(pmax(subject_second, 0))),
    function_mean_l2 = sqrt(sum(diag(function_gram))),
    function_second_l2 = sqrt(sum(diag(function_gram_second))),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)
out <- "posthoc_v1/FULL_POSTERIOR_SIGNAL_ENDPOINT.csv"
if (file.exists(out)) stop("Refusing overwrite: ", out)
write.csv(result, out, row.names = FALSE)
print(result, digits = 7L, row.names = FALSE)
cat("STAGE7G_FULL_POSTERIOR_SECOND_SIGNAL_PASS\n")
