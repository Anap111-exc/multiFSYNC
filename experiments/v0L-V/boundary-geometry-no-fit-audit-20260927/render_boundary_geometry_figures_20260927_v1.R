#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
stopifnot(
  requireNamespace("svglite", quietly = TRUE),
  requireNamespace("ragg", quietly = TRUE),
  capabilities("cairo"))

repo_candidates <- unique(c(
  Sys.getenv("MULTIFSYNC_REPO_ROOT", unset = ""),
  getwd(), file.path(getwd(), "Rcode", "multiFSYNC")))
repo_candidates <- repo_candidates[
  nzchar(repo_candidates) &
    file.exists(file.path(repo_candidates, "DESCRIPTION"))]
if (!length(repo_candidates)) stop("Cannot locate package root.")
repo_root <- normalizePath(repo_candidates[[1L]], mustWork = TRUE)
output_dir <- file.path(
  repo_root, "experiments", "v0L-V",
  "boundary-geometry-no-fit-audit-20260927",
  "boundary_geometry_no_fit_audit_20260927_v1")

metrics <- utils::read.csv(
  file.path(output_dir, "FUNCTION_LEVEL_METRICS.csv"),
  check.names = FALSE)
shapes <- utils::read.csv(
  file.path(output_dir, "REPRESENTATIVE_SHAPES.csv"),
  check.names = FALSE)
stopifnot(nrow(metrics) == 5000L, length(unique(shapes$function_id)) == 8L)

display_labels <- c(
  periodic_mean = "Periodic B-spline mean",
  nonperiodic_feature = "Non-periodic feature")
palette <- c(
  "Periodic B-spline mean" = "#3182BD",
  "Non-periodic feature" = "#D97732")
metrics$display_family <- unname(display_labels[metrics$family])
shapes$display_family <- unname(display_labels[shapes$family])
shapes$reason_label <- factor(
  shapes$selection_reason,
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
    plot_data$value <- pmax(plot_data$value, .Machine$double.eps)
  }
  plot <- ggplot(plot_data, aes(x = family, y = value, fill = family)) +
    geom_violin(trim = TRUE, linewidth = 0.25, alpha = 0.65) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3,
                 fill = "white") +
    scale_fill_manual(values = palette) +
    labs(x = NULL, y = y_label) +
    theme_audit() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))
  if (transform == "log10") plot <- plot + scale_y_log10()
  if (!is.null(reference)) {
    plot <- plot + geom_hline(
      yintercept = reference, linetype = 2,
      linewidth = 0.3, colour = "#666666")
  }
  plot
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
figure_1 <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a")

figure_2 <- ggplot(
  shapes, aes(x = time, y = value, colour = display_family)) +
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

export_additional_formats <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  pdf_target <- file.path(output_dir, paste0(stem, ".pdf"))
  png_target <- file.path(output_dir, paste0(stem, ".png"))
  svg_target <- file.path(output_dir, paste0(stem, ".svg"))
  tiff_target <- file.path(output_dir, paste0(stem, ".tiff"))
  if (file.exists(svg_target) || file.exists(tiff_target)) {
    stop("Refusing to overwrite existing vector/raster export for ", stem)
  }
  temporary_pdf <- tempfile(pattern = paste0(stem, "_"), fileext = ".pdf")
  temporary_png <- tempfile(pattern = paste0(stem, "_"), fileext = ".png")
  temporary_svg <- tempfile(pattern = paste0(stem, "_"), fileext = ".svg")
  temporary_tiff <- tempfile(pattern = paste0(stem, "_"), fileext = ".tiff")
  grDevices::cairo_pdf(
    temporary_pdf, width = width_in, height = height_in,
    family = "Arial", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    temporary_png, width = width_in, height = height_in,
    units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()
  svglite::svglite(temporary_svg, width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(
    temporary_tiff, width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw",
    background = "white")
  print(plot)
  grDevices::dev.off()
  copied <- c(
    file.copy(temporary_pdf, pdf_target, overwrite = TRUE),
    file.copy(temporary_png, png_target, overwrite = TRUE),
    file.copy(temporary_svg, svg_target, overwrite = FALSE),
    file.copy(temporary_tiff, tiff_target, overwrite = FALSE))
  unlink(c(temporary_pdf, temporary_png, temporary_svg, temporary_tiff))
  if (!all(copied)) stop("Unable to copy completed figure exports.")
}

export_additional_formats(figure_1, "FIGURE_1_BOUNDARY_METRICS", 183, 145)
export_additional_formats(figure_2, "FIGURE_2_REPRESENTATIVE_SHAPES", 183, 105)
cat("BOUNDARY_GEOMETRY_FIGURE_EXPORT_PASS\n")
