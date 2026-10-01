#!/usr/bin/env Rscript

# Post-freeze loading-driver orientation; truth never enters blind tables.
options(warn = 2, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(
  "Usage: Rscript evaluate_stage7h_driver_truth_20260929_v1.R ",
  "STAGE7C_ROOT STAGE7G_ROOT BLIND_ROOT OUTPUT_ROOT", call. = FALSE)
stage7c <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g <- normalizePath(args[[2L]], mustWork = TRUE)
blind <- normalizePath(args[[3L]], mustWork = TRUE)
out <- args[[4L]]
if (dir.exists(out)) stop("Refusing existing output directory: ", out)
sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (length(line) != 1L) stop("sha256sum failed: ", path)
  strsplit(line, "[[:space:]]+")[[1L]][[1L]]
}
norm2 <- function(x) sqrt(sum(x^2))
hash_lines <- readLines(file.path(blind, "BLIND_TABLES.sha256"),
                        warn = FALSE)
if (length(hash_lines) != 6L ||
    !all(c("mode=all", "point_count=20") %in%
         readLines(file.path(blind, "BLIND_COMPLETE.txt"), warn = FALSE)))
  stop("Blind freeze marker missing.")
for (line in hash_lines) {
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  if (length(fields) != 2L ||
      !identical(sha256(file.path(blind, fields[[2L]])), fields[[1L]]))
    stop("Blind output hash mismatch.")
}
drivers <- utils::read.csv(file.path(blind, "UPDATE_DRIVER_BLIND.csv"),
                           stringsAsFactors = FALSE)
if (nrow(drivers) != 10000L ||
    length(unique(paste(drivers$fit_id, drivers$iteration))) != 20L)
  stop("Blind driver table cardinality changed.")
manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
if (nrow(manifest) != 3L || nrow(binding) != 3L ||
    !identical(as.character(manifest$fit_id),
               as.character(binding$fit_id)))
  stop("Stage 7G input binding changed.")

# Truth boundary: all blind hashes above have passed.
rows <- list()
idx <- 0L
for (k in seq_len(3L)) {
  fit_id <- manifest$fit_id[[k]]
  data_id <- manifest$data_id[[k]]
  truth_path <- file.path(stage7c, "data", data_id,
                          "sealed_truth", "truth_bundle.rds")
  if (!identical(sha256(truth_path), binding$truth_sha256[[k]]))
    stop("Truth binding changed: ", data_id)
  bundle <- readRDS(truth_path)
  if (!identical(bundle$data_id, data_id) ||
      !identical(bundle$bundle_class,
                 "covariate_routeb_pilot_sealed_truth"))
    stop("Truth bundle identity changed.")
  truth <- bundle$data$true_params
  targets <- list(true_shared = as.numeric(truth$a_true[, 1L]),
                  true_B2 = as.numeric(truth$b_true[[2L]][, 1L]))
  if (any(vapply(targets, length, integer(1L)) != 500L))
    stop("Truth loading dimension changed.")
  subset <- drivers[drivers$fit_id == fit_id, , drop = FALSE]
  for (time in sort(unique(subset$iteration))) {
    state <- subset[subset$iteration == time, , drop = FALSE]
    state <- state[order(state$coordinate), , drop = FALSE]
    if (nrow(state) != 500L ||
        any(state$coordinate != seq_len(500L)))
      stop("Driver coordinate index changed: ", fit_id, " t=", time)
    h_without_shared <- state$h_data + state$h_mean + state$h_beta
    vectors <- list(
      h_data = state$h_data,
      h_mean = state$h_mean,
      h_beta = state$h_beta,
      h_shared = state$h_shared,
      h_without_shared = h_without_shared,
      h_total = state$h_total,
      slab_mean = state$slab_mean,
      effective_loading = state$effective_loading)
    for (field in names(vectors)) {
      value <- vectors[[field]]
      value_norm <- norm2(value)
      for (name in names(targets)) {
        for (metric in c("raw_target", "tau_weighted_target")) {
          target <- targets[[name]]
          if (metric == "tau_weighted_target")
            target <- target * state$tau
          target_norm <- norm2(target)
          signed <- if (value_norm > 0 && target_norm > 0)
            sum(value * target) / (value_norm * target_norm) else NA_real_
          idx <- idx + 1L
          rows[[idx]] <- data.frame(
            fit_id = fit_id, data_id = data_id,
            iteration = time, field = field,
            target = name, target_metric = metric,
            signed_cosine = signed,
            abs_cosine = abs(signed),
            vector_norm = value_norm,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  rm(bundle, truth)
}
result <- do.call(rbind, rows)
if (nrow(result) != 640L) stop("Truth driver row count changed: ",
                              nrow(result))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, file.path(out, "DRIVER_TRUTH_ALIGNMENT.csv"),
                 row.names = FALSE)
writeLines(c("STAGE7H_DRIVER_TRUTH_ALIGNMENT_COMPLETE",
             paste0("blind_hash_list_sha256=",
                    sha256(file.path(blind, "BLIND_TABLES.sha256"))),
             paste0("rows=", nrow(result))),
           file.path(out, "DRIVER_ALIGNMENT_COMPLETE.txt"),
           useBytes = TRUE)
cat("STAGE7H_DRIVER_ALIGNMENT_COMPLETE rows=", nrow(result), "\n",
    sep = "")
