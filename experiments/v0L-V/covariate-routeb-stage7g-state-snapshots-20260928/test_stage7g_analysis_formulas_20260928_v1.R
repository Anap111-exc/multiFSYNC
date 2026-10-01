#!/usr/bin/env Rscript

# Truth-free shape and energy sanity check on an already frozen Stage 7F fit.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript TEST ANALYZER STAGE7F_FIT_RDS")
}
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
lines <- readLines(args[[1L]], warn = FALSE)
first <- which(grepl("^trap_weights <- function", lines))
last <- which(grepl("^state_rows <- list\\(\\)", lines))
stopifnot(length(first) == 1L, length(last) == 1L, first < last)
eval(parse(text = paste(lines[first:(last - 1L)], collapse = "\n")))
fit <- readRDS(args[[2L]])
weights <- trap_weights(fit$time_g)
stopifnot(length(weights) == length(fit$time_g),
          abs(sum(weights) - 1) < 1e-12)
for (block in c("shared", "specific")) for (study in 1:2) {
  state <- block_state(fit, block, study)
  scores <- as.matrix(state$scores)
  functions <- fit$C_g %*% as.matrix(state$coefficients)
  gram <- crossprod(functions, sweep(functions, 1L, weights, "*"))
  covariance_energy <- vapply(state$score_covariances,
                              sum_trace_product, numeric(1L), y = gram)
  stopifnot(length(state$loading) == 500L,
            length(state$score_covariances) == 30L,
            nrow(scores) == 30L, ncol(scores) == 2L,
            nrow(functions) == 51L, ncol(functions) == 2L,
            all(is.finite(covariance_energy)),
            min(covariance_energy) > -1e-7,
            is.finite(time_score_l2_rms(scores, functions, weights)))
}
cat("STAGE7G_ANALYSIS_FORMULAS_PASS frozen_fit=1 blocks=4 truth_read=FALSE\n")
