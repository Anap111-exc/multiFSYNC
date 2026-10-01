#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)

repo_root_env <- Sys.getenv("MULTIFSYNC_REPO_ROOT", unset = "")
repo_candidates <- unique(c(
  repo_root_env,
  getwd(),
  file.path(getwd(), "Rcode", "multiFSYNC")))
repo_candidates <- repo_candidates[
  nzchar(repo_candidates) &
    file.exists(file.path(repo_candidates, "DESCRIPTION"))]
if (length(repo_candidates)) {
  repo_root <- normalizePath(repo_candidates[[1L]], mustWork = TRUE)
} else {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot locate package root.")
  script_path <- normalizePath(
    sub("^--file=", "", file_arg), mustWork = TRUE)
  repo_root <- normalizePath(
    file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
}

script_dir <- file.path(
  repo_root, "experiments", "v0L-V",
  "boundary-geometry-no-fit-audit-20260927")
output_dir <- file.path(script_dir, "boundary_geometry_no_fit_audit_20260927_v1")
if (dir.exists(output_dir)) {
  stop("Refusing to overwrite existing output: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("Splinets", "ggplot2", "patchwork", "ragg", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
if (!capabilities("cairo")) stop("R cairo capability is required for PDF export.")

suppressPackageStartupMessages({
  library(splines)
  library(ggplot2)
  library(patchwork)
})
source(file.path(repo_root, "R", "utils_multi.R"), local = .GlobalEnv)
source(file.path(repo_root, "R", "generate_data_structured.R"),
       local = .GlobalEnv)

audit_id <- "BOUNDARY_GEOMETRY_NO_FIT_AUDIT_V1_20260927"
audit_seed <- 20260927L
n_grid <- 2001L
K_fit <- 9L
n_mean_draw <- 20L
mean_p <- 50L
mean_S <- 2L
n_factor_draw <- 500L
M_by_block <- c(2L, 2L, 2L)
mean_basis_dim <- 8L
factor_basis_dim <- 8L
mean_smoothness_alpha <- 1.519826404906808
factor_smoothness_alpha <- 3
mean_amplitude <- 0.6

set.seed(audit_seed)
mean_generator_seeds <- sample.int(.Machine$integer.max, n_mean_draw)
mean_design_seeds <- sample.int(.Machine$integer.max, n_mean_draw)
factor_generator_seeds <- sample.int(.Machine$integer.max, n_factor_draw)
factor_design_seeds <- sample.int(.Machine$integer.max, n_factor_draw)

t_grid <- seq(0, 1, length.out = n_grid)
dt <- diff(t_grid)
w_grid <- c(dt[[1L]] / 2,
            (dt[-1L] + dt[-length(dt)]) / 2,
            dt[[length(dt)]] / 2)
edge_left <- t_grid <= 0.1
edge_right <- t_grid >= 0.9
edge_any <- edge_left | edge_right
interior <- !edge_any

eval_collection <- function(object, derivative = 0L, x = t_grid) {
  evaluated_object <- object
  if (derivative > 0L) {
    for (i in seq_len(derivative)) {
      evaluated_object <- Splinets::deriva(evaluated_object)
    }
  }
  Splinets::evspline(evaluated_object, x = x)[, -1L, drop = FALSE]
}

derivative_object <- function(object, derivative) {
  out <- object
  if (derivative > 0L) {
    for (i in seq_len(derivative)) out <- Splinets::deriva(out)
  }
  out
}

make_osullivan_design <- function(seed) {
  set.seed(seed)
  subject_times <- unlist(lapply(seq_len(60L), function(i) {
    stats::runif(sample(6:9, 1L))
  }), use.names = FALSE)
  pooled <- sort(unique(subject_times))
  int_knots <- unname(stats::quantile(
    pooled, seq(0, 1, length.out = K_fit)[-c(1L, K_fit)],
    type = 7L))
  C_fit <- cbind(
    1, t_grid,
    ZOSull(t_grid, range.x = c(0, 1), intKnots = int_knots))
  list(
    seed = seed,
    int_knots = int_knots,
    C = C_fit,
    C_weighted = sqrt(w_grid) * C_fit,
    n_times = length(subject_times))
}

projection_metrics <- function(values, design) {
  coefficients <- qr.solve(
    design$C_weighted, sqrt(w_grid) * values)
  projected <- design$C %*% coefficients
  residual <- values - projected
  energy <- colSums(w_grid * values^2)
  ise <- colSums(w_grid * residual^2)
  relative_ise <- ise / pmax(energy, .Machine$double.eps)
  list(
    projection_ise = ise,
    projection_relative_ise = relative_ise,
    projection_relative_l2 = sqrt(pmax(relative_ise, 0)))
}

collection_metrics <- function(
    object, family, boundary_condition, draw_id, generator_seed,
    design, role, study = NA_integer_, target_energy) {
  values0 <- eval_collection(object, 0L)
  values1 <- eval_collection(object, 1L)
  values2 <- eval_collection(object, 2L)
  endpoint0 <- eval_collection(object, 0L, c(0, 1))
  endpoint1 <- eval_collection(object, 1L, c(0, 1))
  endpoint2 <- eval_collection(object, 2L, c(0, 1))

  gram0 <- Splinets::gramian(object)
  gram2 <- Splinets::gramian(derivative_object(object, 2L))
  energy_exact <- diag(gram0)
  energy_grid <- colSums(w_grid * values0^2)
  roughness <- diag(gram2)
  normalized_roughness <- roughness / pmax(energy_exact, .Machine$double.eps)
  edge_left_energy <- colSums(
    w_grid[edge_left] * values0[edge_left, , drop = FALSE]^2)
  edge_right_energy <- colSums(
    w_grid[edge_right] * values0[edge_right, , drop = FALSE]^2)
  edge_fraction <- (edge_left_energy + edge_right_energy) /
    pmax(energy_grid, .Machine$double.eps)
  edge_peak <- apply(abs(values0[edge_any, , drop = FALSE]), 2L, max)
  interior_peak <- apply(abs(values0[interior, , drop = FALSE]), 2L, max)
  edge_peak_ratio <- edge_peak / pmax(interior_peak, .Machine$double.eps)
  max_abs0 <- apply(abs(values0), 2L, max)
  max_abs1 <- apply(abs(values1), 2L, max)
  max_abs2 <- apply(abs(values2), 2L, max)

  closure0 <- abs(endpoint0[1L, ] - endpoint0[2L, ])
  closure1 <- abs(endpoint1[1L, ] - endpoint1[2L, ])
  closure2 <- abs(endpoint2[1L, ] - endpoint2[2L, ])
  zero0 <- apply(abs(endpoint0), 2L, max)
  zero1 <- apply(abs(endpoint1), 2L, max)
  zero2 <- apply(abs(endpoint2), 2L, max)
  relative_constraint_error <- if (boundary_condition == "periodic_closure") {
    pmax(
      closure0 / pmax(max_abs0, 1),
      closure1 / pmax(max_abs1, 1),
      closure2 / pmax(max_abs2, 1))
  } else {
    pmax(
      zero0 / pmax(max_abs0, 1),
      zero1 / pmax(max_abs1, 1),
      zero2 / pmax(max_abs2, 1))
  }
  projection <- projection_metrics(values0, design)

  n_component <- ncol(values0)
  data.frame(
    audit_id = audit_id,
    family = family,
    boundary_condition = boundary_condition,
    draw_id = draw_id,
    generator_seed = generator_seed,
    design_seed = design$seed,
    role = role,
    study = study,
    component = seq_len(n_component),
    f_0 = endpoint0[1L, ],
    f_1 = endpoint0[2L, ],
    d1_0 = endpoint1[1L, ],
    d1_1 = endpoint1[2L, ],
    d2_0 = endpoint2[1L, ],
    d2_1 = endpoint2[2L, ],
    closure_error_f = closure0,
    closure_error_d1 = closure1,
    closure_error_d2 = closure2,
    zero_boundary_error_f = zero0,
    zero_boundary_error_d1 = zero1,
    zero_boundary_error_d2 = zero2,
    relative_constraint_error = relative_constraint_error,
    energy_exact = energy_exact,
    energy_grid = energy_grid,
    target_energy = target_energy,
    energy_target_relative_error =
      abs(energy_exact - target_energy) / target_energy,
    roughness = roughness,
    normalized_roughness = normalized_roughness,
    left_edge_energy_fraction =
      edge_left_energy / pmax(energy_grid, .Machine$double.eps),
    right_edge_energy_fraction =
      edge_right_energy / pmax(energy_grid, .Machine$double.eps),
    edge_energy_fraction = edge_fraction,
    edge_peak_ratio = edge_peak_ratio,
    projection_ise = projection$projection_ise,
    projection_relative_ise = projection$projection_relative_ise,
    projection_relative_l2 = projection$projection_relative_l2,
    constraint_fail = relative_constraint_error > 1e-8,
    edge_energy_high_flag = edge_fraction > 0.5,
    edge_energy_low_flag = edge_fraction < 0.01,
    edge_peak_flag = edge_peak_ratio > 3,
    projection_l2_high_flag = projection$projection_relative_l2 > 0.6,
    finite = is.finite(energy_exact) & is.finite(roughness) &
      is.finite(edge_fraction) &
      is.finite(projection$projection_relative_l2),
    stringsAsFactors = FALSE)
}

design_rows <- list()
metric_rows <- list()
design_pos <- 0L
metric_pos <- 0L

for (draw_id in seq_len(n_mean_draw)) {
  design <- make_osullivan_design(mean_design_seeds[[draw_id]])
  design_pos <- design_pos + 1L
  design_rows[[design_pos]] <- data.frame(
    family = "periodic_mean",
    draw_id = draw_id,
    design_seed = design$seed,
    n_times = design$n_times,
    knot_index = seq_along(design$int_knots),
    internal_knot = design$int_knots)

  set.seed(mean_generator_seeds[[draw_id]])
  means <- construct_independent_splinet_means(
    S = mean_S, p = mean_p, n_g = 201L,
    basis_dim = mean_basis_dim,
    smoothness_alpha = mean_smoothness_alpha,
    amplitude = mean_amplitude)
  for (study in seq_len(mean_S)) {
    metric_pos <- metric_pos + 1L
    metric_rows[[metric_pos]] <- collection_metrics(
      object = means$evaluation[[study]]$splinet_object,
      family = "periodic_mean",
      boundary_condition = "periodic_closure",
      draw_id = draw_id,
      generator_seed = mean_generator_seeds[[draw_id]],
      design = design,
      role = "mean",
      study = study,
      target_energy = mean_amplitude^2 * 5 / 8)
  }
}

factor_roles <- c("shared", "study1_specific", "study2_specific")
for (draw_id in seq_len(n_factor_draw)) {
  design <- make_osullivan_design(factor_design_seeds[[draw_id]])
  design_pos <- design_pos + 1L
  design_rows[[design_pos]] <- data.frame(
    family = "nonperiodic_feature",
    draw_id = draw_id,
    design_seed = design$seed,
    n_times = design$n_times,
    knot_index = seq_along(design$int_knots),
    internal_knot = design$int_knots)

  set.seed(factor_generator_seeds[[draw_id]])
  blocks <- construct_iid_common_splinet_stiefel_blocks(
    M_by_block = M_by_block,
    n_g = 201L,
    basis_dim = factor_basis_dim,
    smoothness_alpha = factor_smoothness_alpha)
  for (block_id in seq_along(blocks)) {
    metric_pos <- metric_pos + 1L
    metric_rows[[metric_pos]] <- collection_metrics(
      object = blocks[[block_id]]$evaluation$splinet_object,
      family = "nonperiodic_feature",
      boundary_condition = "zero_boundary",
      draw_id = draw_id,
      generator_seed = factor_generator_seeds[[draw_id]],
      design = design,
      role = factor_roles[[block_id]],
      study = if (block_id == 1L) NA_integer_ else block_id - 1L,
      target_energy = 1)
  }
}

metrics <- do.call(rbind, metric_rows)
design_knots <- do.call(rbind, design_rows)
row.names(metrics) <- NULL
row.names(design_knots) <- NULL
metrics$function_id <- sprintf(
  "%s_d%03d_%s_s%s_c%03d",
  metrics$family, metrics$draw_id, metrics$role,
  ifelse(is.na(metrics$study), "NA", metrics$study), metrics$component)

quantile_value <- function(x, probability) {
  unname(stats::quantile(x, probability, na.rm = TRUE, type = 7L))
}

summarize_group <- function(x) {
  data.frame(
    n_functions = nrow(x),
    energy_min = min(x$energy_exact),
    energy_median = stats::median(x$energy_exact),
    energy_max = max(x$energy_exact),
    energy_target_relative_error_max =
      max(x$energy_target_relative_error),
    roughness_q05 = quantile_value(x$normalized_roughness, 0.05),
    roughness_median = stats::median(x$normalized_roughness),
    roughness_q95 = quantile_value(x$normalized_roughness, 0.95),
    edge_energy_q05 = quantile_value(x$edge_energy_fraction, 0.05),
    edge_energy_median = stats::median(x$edge_energy_fraction),
    edge_energy_q95 = quantile_value(x$edge_energy_fraction, 0.95),
    edge_energy_max = max(x$edge_energy_fraction),
    edge_peak_ratio_q95 = quantile_value(x$edge_peak_ratio, 0.95),
    edge_peak_ratio_max = max(x$edge_peak_ratio),
    relative_constraint_error_max = max(x$relative_constraint_error),
    projection_relative_l2_median =
      stats::median(x$projection_relative_l2),
    projection_relative_l2_q95 =
      quantile_value(x$projection_relative_l2, 0.95),
    projection_relative_l2_max = max(x$projection_relative_l2),
    projection_relative_ise_median =
      stats::median(x$projection_relative_ise),
    projection_relative_ise_q95 =
      quantile_value(x$projection_relative_ise, 0.95),
    constraint_fail_n = sum(x$constraint_fail),
    edge_energy_high_n = sum(x$edge_energy_high_flag),
    edge_energy_low_n = sum(x$edge_energy_low_flag),
    edge_peak_flag_n = sum(x$edge_peak_flag),
    projection_l2_high_n = sum(x$projection_l2_high_flag),
    nonfinite_n = sum(!x$finite),
    stringsAsFactors = FALSE)
}

family_summary <- do.call(rbind, lapply(split(metrics, metrics$family), function(x) {
  cbind(family = x$family[[1L]], summarize_group(x))
}))
role_key <- interaction(metrics$family, metrics$role, drop = TRUE)
role_summary <- do.call(rbind, lapply(split(metrics, role_key), function(x) {
  cbind(family = x$family[[1L]], role = x$role[[1L]], summarize_group(x))
}))
row.names(family_summary) <- NULL
row.names(role_summary) <- NULL

periodic_metrics <- metrics[metrics$family == "periodic_mean", , drop = FALSE]
feature_metrics <- metrics[
  metrics$family == "nonperiodic_feature", , drop = FALSE]
checks <- data.frame(
  check_id = c(
    "expected_function_counts",
    "all_metrics_finite",
    "periodic_endpoint_closure",
    "feature_zero_boundary",
    "mean_energy_target",
    "feature_energy_target",
    "no_model_fit_called"),
  passed = c(
    nrow(periodic_metrics) == n_mean_draw * mean_S * mean_p &&
      nrow(feature_metrics) == n_factor_draw * sum(M_by_block),
    all(metrics$finite),
    max(periodic_metrics$relative_constraint_error) <= 1e-8,
    max(feature_metrics$relative_constraint_error) <= 1e-8,
    max(periodic_metrics$energy_target_relative_error) <= 1e-8,
    max(feature_metrics$energy_target_relative_error) <= 1e-8,
    TRUE),
  detail = c(
    sprintf("periodic_mean=%d;nonperiodic_feature=%d",
            nrow(periodic_metrics), nrow(feature_metrics)),
    sprintf("nonfinite=%d", sum(!metrics$finite)),
    sprintf("max_relative_error=%.3e",
            max(periodic_metrics$relative_constraint_error)),
    sprintf("max_relative_error=%.3e",
            max(feature_metrics$relative_constraint_error)),
    sprintf("max_relative_error=%.3e",
            max(periodic_metrics$energy_target_relative_error)),
    sprintf("max_relative_error=%.3e",
            max(feature_metrics$energy_target_relative_error)),
    "script contains generation, projection, summary and plotting only"),
  stringsAsFactors = FALSE)

select_representatives <- function(x) {
  candidates <- list(
    median_edge = which.min(abs(
      x$edge_energy_fraction - stats::median(x$edge_energy_fraction))),
    max_edge = which.max(x$edge_energy_fraction),
    max_roughness = which.max(x$normalized_roughness),
    max_projection = which.max(x$projection_relative_l2))
  selected <- integer()
  reasons <- character()
  for (reason in names(candidates)) {
    index <- candidates[[reason]]
    if (index %in% selected) {
      ordering_variable <- switch(
        reason,
        max_edge = x$edge_energy_fraction,
        max_roughness = x$normalized_roughness,
        max_projection = x$projection_relative_l2,
        -abs(x$edge_energy_fraction - stats::median(x$edge_energy_fraction)))
      ordered <- order(ordering_variable, decreasing = TRUE)
      index <- ordered[which(!ordered %in% selected)[1L]]
    }
    selected <- c(selected, index)
    reasons <- c(reasons, reason)
  }
  out <- x[selected, , drop = FALSE]
  out$selection_reason <- reasons
  out
}

representatives <- do.call(rbind, lapply(
  split(metrics, metrics$family), select_representatives))
row.names(representatives) <- NULL

regenerate_shape <- function(row) {
  if (row$family == "periodic_mean") {
    set.seed(row$generator_seed)
    object <- construct_independent_splinet_means(
      S = mean_S, p = mean_p, n_g = 201L,
      basis_dim = mean_basis_dim,
      smoothness_alpha = mean_smoothness_alpha,
      amplitude = mean_amplitude)$evaluation[[row$study]]$splinet_object
  } else {
    set.seed(row$generator_seed)
    block_id <- match(row$role, factor_roles)
    object <- construct_iid_common_splinet_stiefel_blocks(
      M_by_block = M_by_block, n_g = 201L,
      basis_dim = factor_basis_dim,
      smoothness_alpha = factor_smoothness_alpha)[[block_id]]$
      evaluation$splinet_object
  }
  plot_grid <- seq(0, 1, length.out = 401L)
  values <- eval_collection(
    object, derivative = 0L, x = plot_grid)[, row$component]
  data.frame(
    function_id = row$function_id,
    family = row$family,
    selection_reason = row$selection_reason,
    time = plot_grid,
    value = values,
    stringsAsFactors = FALSE)
}

shape_rows <- lapply(seq_len(nrow(representatives)), function(i) {
  regenerate_shape(representatives[i, , drop = FALSE])
})
representative_shapes <- do.call(rbind, shape_rows)

source_files <- c(
  file.path(repo_root, "R", "generate_data_structured.R"),
  file.path(repo_root, "R", "utils_multi.R"),
  file.path(script_dir,
            "run_boundary_geometry_no_fit_audit_20260927_v1.R"))
source_binding <- data.frame(
  relative_path = vapply(source_files, function(path) {
    substring(normalizePath(path, winslash = "/"),
              nchar(normalizePath(repo_root, winslash = "/")) + 2L)
  }, character(1)),
  sha256 = vapply(source_files, function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1)),
  bytes = file.info(source_files)$size,
  stringsAsFactors = FALSE)

utils::write.csv(metrics, file.path(output_dir, "FUNCTION_LEVEL_METRICS.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(family_summary, file.path(output_dir, "FAMILY_SUMMARY.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(role_summary, file.path(output_dir, "ROLE_SUMMARY.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(design_knots, file.path(output_dir, "OSULLIVAN_DESIGN_KNOTS.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(checks, file.path(output_dir, "IMPLEMENTATION_CHECKS.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(representatives,
                 file.path(output_dir, "REPRESENTATIVE_SELECTION.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(representative_shapes,
                 file.path(output_dir, "REPRESENTATIVE_SHAPES.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(source_binding, file.path(output_dir, "SOURCE_BINDING.csv"),
                 row.names = FALSE, na = "")

display_labels <- c(
  periodic_mean = "Periodic B-spline mean",
  nonperiodic_feature = "Non-periodic feature")
palette <- c(
  "Periodic B-spline mean" = "#3182BD",
  "Non-periodic feature" = "#D97732")
metrics$display_family <- unname(display_labels[metrics$family])
representative_shapes$display_family <- unname(
  display_labels[representative_shapes$family])
representative_shapes$reason_label <- factor(
  representative_shapes$selection_reason,
  levels = c("median_edge", "max_edge", "max_roughness", "max_projection"),
  labels = c("Median edge energy", "Maximum edge energy",
             "Maximum roughness", "Maximum projection error"))

theme_audit <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.text = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.tag = element_text(size = 9, face = "bold"),
      legend.position = "none")
}

violin_panel <- function(variable, y_label, transform = "identity",
                         reference = NULL) {
  plot_data <- data.frame(
    family = metrics$display_family,
    value = metrics[[variable]])
  if (transform == "log10") {
    # Preserve exact/near-exact zeros in the plotted population by placing
    # them at the smallest positive plotting value rather than dropping rows.
    plot_data$value <- pmax(plot_data$value, .Machine$double.eps)
  }
  p <- ggplot(plot_data, aes(x = family, y = value, fill = family)) +
    geom_violin(trim = TRUE, linewidth = 0.25, alpha = 0.65) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3,
                 fill = "white") +
    scale_fill_manual(values = palette) +
    labs(x = NULL, y = y_label) +
    theme_audit() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  if (transform == "log10") p <- p + scale_y_log10()
  if (!is.null(reference)) {
    p <- p + geom_hline(yintercept = reference, linetype = 2,
                        linewidth = 0.3, colour = "#666666")
  }
  p
}

p_a <- violin_panel(
  "relative_constraint_error", "Relative endpoint-constraint error",
  transform = "log10", reference = 1e-8) +
  ggtitle("Endpoint constraints")
p_b <- violin_panel(
  "edge_energy_fraction", "Energy in t <= 0.1 or t >= 0.9",
  reference = 0.2) +
  ggtitle("Boundary energy")
p_c <- violin_panel(
  "projection_relative_l2", "Relative L2 projection floor",
  transform = "log10", reference = 0.6) +
  ggtitle("O'Sullivan K=9 approximation")
p_d <- violin_panel(
  "normalized_roughness", "Integral(f''(t)^2) / Integral(f(t)^2)",
  transform = "log10") +
  ggtitle("Normalized roughness")

metric_figure <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a")

shape_figure <- ggplot(
  representative_shapes,
  aes(x = time, y = value, colour = display_family)) +
  annotate("rect", xmin = -Inf, xmax = 0.1, ymin = -Inf, ymax = Inf,
           fill = "#D9D9D9", alpha = 0.25) +
  annotate("rect", xmin = 0.9, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D9D9D9", alpha = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#777777") +
  geom_line(linewidth = 0.55) +
  facet_grid(display_family ~ reason_label, scales = "free_y") +
  scale_colour_manual(values = palette) +
  labs(x = "t", y = "Function value") +
  theme_audit() +
  theme(strip.text.y = element_text(angle = 0),
        panel.spacing = grid::unit(3, "mm"))

save_figure <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  # Cairo on Windows cannot open a graphics device directly on this project's
  # non-ASCII absolute path. Render in R's temporary directory, then copy the
  # completed files without changing their bytes or visual content.
  temporary_pdf <- tempfile(pattern = paste0(stem, "_"), fileext = ".pdf")
  temporary_png <- tempfile(pattern = paste0(stem, "_"), fileext = ".png")
  grDevices::cairo_pdf(
    temporary_pdf,
    width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    temporary_png,
    width = width_in, height = height_in, units = "in", res = 300,
    background = "white")
  print(plot)
  grDevices::dev.off()
  copied <- c(
    file.copy(temporary_pdf, file.path(output_dir, paste0(stem, ".pdf")),
              overwrite = FALSE),
    file.copy(temporary_png, file.path(output_dir, paste0(stem, ".png")),
              overwrite = FALSE))
  unlink(c(temporary_pdf, temporary_png))
  if (!all(copied)) stop("Unable to copy completed figure exports.")
}

save_figure(metric_figure, "FIGURE_1_BOUNDARY_METRICS", 183, 145)
save_figure(shape_figure, "FIGURE_2_REPRESENTATIVE_SHAPES", 183, 105)

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "SESSION_INFO.txt"))
writeLines(c(
  paste0("audit_id=", audit_id),
  paste0("completed_at=", format(Sys.time(), tz = "Asia/Shanghai")),
  paste0("periodic_mean_functions=", nrow(periodic_metrics)),
  paste0("nonperiodic_feature_functions=", nrow(feature_metrics)),
  "formal_data_generated=FALSE",
  "model_fit_calls=0",
  "cavi_calls=0",
  "elbo_calls=0",
  "continuation_calls=0"),
  file.path(output_dir, "NO_FIT_BOUNDARY_AUDIT_COMPLETE.txt"))

if (!all(checks$passed)) {
  writeLines("FAIL", file.path(output_dir, "AUDIT_STATUS.txt"))
  stop("Boundary audit failed one or more implementation checks.")
}
writeLines("PASS", file.path(output_dir, "AUDIT_STATUS.txt"))
cat("BOUNDARY_GEOMETRY_NO_FIT_AUDIT_PASS\n")
