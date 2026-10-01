#!/usr/bin/env Rscript

# Generation-only calibration for the independent B-spline mean truth.
# No model fit, recovery metric, or experimental outcome is read here.
options(warn = 1, stringsAsFactors = FALSE)

stopifnot(requireNamespace("Splinets", quietly = TRUE))

basis_dim <- 8L
degree <- 3L
knots <- seq(0, 1, length.out = basis_dim + 1L)
system <- Splinets::splinet(
  knots, degree = degree, periodic = TRUE, norm = TRUE)
gram <- Splinets::gramian(system$os)
stopifnot(max(abs(gram - diag(basis_dim))) < 1e-8)

second <- Splinets::deriva(Splinets::deriva(system$os))
penalty <- Splinets::gramian(second)
penalty <- (penalty + t(penalty)) / 2
integrated <- Splinets::integra(system$os)
boundary <- Splinets::evspline(integrated, x = c(0, 1))[, -1L, drop = FALSE]
integral_vector <- as.numeric(boundary[2L, ] - boundary[1L, ])
cat("integral_vector_norm2=", format(sum(integral_vector^2), digits = 16),
    "\n", sep = "")
stopifnot(all(is.finite(integral_vector)), sum(integral_vector^2) > 0)

# Work directly in the exact coefficient subspace m'b = 0.  Filtering first
# and projecting afterward can remove the smoothest direction and create an
# artificial roughness floor.
null_basis <- qr.Q(qr(matrix(integral_vector, ncol = 1L)), complete = TRUE)[
  , -1L, drop = FALSE]
restricted_penalty <- crossprod(null_basis, penalty %*% null_basis)
restricted_penalty <- (restricted_penalty + t(restricted_penalty)) / 2
eig <- eigen(restricted_penalty, symmetric = TRUE)
values <- pmax(eig$values, 0)
scale <- stats::median(values)
stopifnot(is.finite(scale), scale > 0)

target <- 64 * pi^4
alpha_grid <- c(0, seq(0.01, 20, length.out = 1000L))
n_draws <- 20000L
set.seed(99310000L)
gaussian <- matrix(
  stats::rnorm((basis_dim - 1L) * n_draws), basis_dim - 1L, n_draws)

calibrate_one <- function(alpha) {
  # Heat-kernel spectral smoothing varies continuously from an isotropic draw
  # at alpha=0 toward the smoothest zero-integral direction.
  weights <- exp(-0.5 * alpha * values / scale)
  map <- null_basis %*%
    sweep(eig$vectors, 2L, weights, FUN = "*") %*% t(eig$vectors)
  coefficients <- map %*% gaussian
  energy <- colSums(coefficients^2)
  roughness <- colSums(coefficients * (penalty %*% coefficients)) / energy
  data.frame(
    alpha = alpha,
    median_normalized_roughness = stats::median(roughness),
    mean_normalized_roughness = mean(roughness),
    target_normalized_roughness = target,
    absolute_log_ratio = abs(log(stats::median(roughness) / target)),
    stringsAsFactors = FALSE)
}

results <- do.call(rbind, lapply(alpha_grid, calibrate_one))
coarse <- results[which.min(results$absolute_log_ratio), , drop = FALSE]
lower <- max(0, coarse$alpha - 0.1)
upper <- coarse$alpha + 0.1
optimum <- stats::optimize(
  function(alpha) calibrate_one(alpha)$absolute_log_ratio,
  interval = c(lower, upper), tol = 1e-10)
best <- calibrate_one(optimum$minimum)

cat("BSPLINE_MEAN_ROUGHNESS_CALIBRATION\n")
cat("basis_dim=", basis_dim, "\n", sep = "")
cat("draws=", n_draws, "\n", sep = "")
cat("target=", format(target, digits = 16), "\n", sep = "")
cat("selected_alpha=", format(best$alpha, digits = 16), "\n", sep = "")
cat("selected_median=",
    format(best$median_normalized_roughness, digits = 16), "\n", sep = "")
cat("selected_ratio=",
    format(best$median_normalized_roughness / target, digits = 16),
    "\n", sep = "")

if (length(commandArgs(trailingOnly = TRUE))) {
  output <- commandArgs(trailingOnly = TRUE)[[1L]]
  if (file.exists(output)) stop("Refusing to overwrite: ", output)
  utils::write.csv(results, output, row.names = FALSE, na = "")
}
