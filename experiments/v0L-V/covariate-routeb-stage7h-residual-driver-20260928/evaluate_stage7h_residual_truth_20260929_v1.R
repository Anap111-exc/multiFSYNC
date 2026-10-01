#!/usr/bin/env Rscript

# Truth is opened only after the complete Stage 7H blind tables pass hashes.
options(warn = 2, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(
  "Usage: Rscript evaluate_stage7h_residual_truth_20260929_v1.R ",
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
marker <- readLines(file.path(blind, "BLIND_COMPLETE.txt"), warn = FALSE)
if (!all(c("mode=all", "point_count=20") %in% marker))
  stop("Complete frozen blind table marker missing.")
file_hashes <- readLines(file.path(blind, "BLIND_TABLES.sha256"),
                         warn = FALSE)
if (length(file_hashes) != 6L) stop("Blind hash list malformed.")
for (line in file_hashes) {
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  if (length(fields) != 2L ||
      !identical(sha256(file.path(blind, fields[[2L]])), fields[[1L]]))
    stop("Blind table hash mismatch: ", line)
}
parity <- utils::read.csv(file.path(blind, "RECONSTRUCTION_PARITY.csv"),
                          stringsAsFactors = FALSE)
geometry <- utils::read.csv(file.path(blind, "RESIDUAL_GEOMETRY_BLIND.csv"),
                            stringsAsFactors = FALSE)
if (nrow(parity) != 646L || !all(parity$pass) ||
    nrow(geometry) != 160L ||
    length(unique(paste(geometry$fit_id,
                        geometry$iteration))) != 20L)
  stop("Blind parity or geometry cardinality changed.")
vectors <- readRDS(file.path(blind,
                             "RESIDUAL_EIGENVECTORS_BLIND.rds"))
manifest <- utils::read.csv(file.path(stage7g, "FIT_MANIFEST.csv"),
                            stringsAsFactors = FALSE)
binding <- utils::read.csv(file.path(stage7g, "analysis_v1",
                                     "INPUT_BINDING.csv"),
                           stringsAsFactors = FALSE)
if (nrow(manifest) != 3L || nrow(binding) != 3L ||
    !identical(as.character(manifest$fit_id),
               as.character(binding$fit_id)))
  stop("Stage 7G manifest mismatch.")

# All prior checks are truth-blind and have completed at this boundary.
rows <- list()
input_rows <- list()
index <- 0L
for (k in seq_len(3L)) {
  fit_id <- manifest$fit_id[[k]]
  data_id <- manifest$data_id[[k]]
  truth_path <- file.path(stage7c, "data", data_id,
                          "sealed_truth", "truth_bundle.rds")
  truth_hash <- sha256(truth_path)
  if (!identical(truth_hash, binding$truth_sha256[[k]]))
    stop("Truth bundle hash mismatch: ", data_id)
  bundle <- readRDS(truth_path)
  if (!identical(bundle$bundle_class,
                 "covariate_routeb_pilot_sealed_truth") ||
      !identical(bundle$data_id, data_id))
    stop("Truth bundle contract mismatch: ", data_id)
  truth <- bundle$data$true_params
  targets <- list(true_shared = as.numeric(truth$a_true[, 1L]),
                  true_B2 = as.numeric(truth$b_true[[2L]][, 1L]))
  if (any(vapply(targets, length, integer(1L)) != 500L))
    stop("Truth loading length changed.")
  targets <- lapply(targets, function(x) x / sqrt(sum(x^2)))
  target_cos <- abs(sum(targets$true_shared * targets$true_B2))
  input_rows[[k]] <- data.frame(fit_id = fit_id,
                                data_id = data_id,
                                truth_sha256 = truth_hash,
                                true_shared_B2_abs_cosine = target_cos,
                                stringsAsFactors = FALSE)
  subset <- geometry[geometry$fit_id == fit_id, , drop = FALSE]
  for (r in seq_len(nrow(subset))) {
    row <- subset[r, , drop = FALSE]
    splits <- if (row$stage[[1L]] == "minus_shared")
      c("full", "odd", "even") else "full"
    for (split in splits) {
      domain_key <- if (split == "full") row$domain[[1L]] else
        paste0(row$domain[[1L]], "_", split)
      key <- paste(fit_id, row$iteration[[1L]],
                   row$stage[[1L]], domain_key, sep = "__")
      U <- vectors[[key]]
      if (!is.matrix(U) || !identical(dim(U), c(500L, 5L)) ||
          max(abs(crossprod(U) - diag(5L))) > 1e-6)
        stop("Blind eigenvector contract mismatch: ", key)
      for (name in names(targets)) {
        cproj <- as.vector(crossprod(U, targets[[name]]))
        projections <- cumsum(cproj^2)
        if (any(!is.finite(projections)) ||
            max(projections) > 1 + 1e-6)
          stop("Invalid truth subspace projection: ", key)
        index <- index + 1L
        rows[[index]] <- data.frame(
          fit_id = fit_id, data_id = data_id,
          iteration = row$iteration[[1L]],
          residual_stage = row$stage[[1L]],
          domain = row$domain[[1L]], split = split,
          true_direction = name,
          projection_sq_top1 = projections[[1L]],
          projection_sq_top2 = projections[[2L]],
          projection_sq_top3 = projections[[3L]],
          projection_sq_top5 = projections[[5L]],
          true_shared_B2_abs_cosine = target_cos,
          sigma1 = row$sigma1[[1L]], sigma2 = row$sigma2[[1L]],
          gap1 = row$gap1[[1L]], gap5 = row$gap5[[1L]],
          stringsAsFactors = FALSE)
      }
    }
  }
  rm(bundle, truth)
}

out_rows <- do.call(rbind, rows)
if (nrow(out_rows) != 480L) stop("Truth output cardinality unexpected: ",
                               nrow(out_rows))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Could not create truth output directory.")
utils::write.csv(out_rows, file.path(out, "RESIDUAL_TRUTH_ALIGNMENT.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, input_rows),
                 file.path(out, "TRUTH_INPUT_BINDING.csv"),
                 row.names = FALSE)
writeLines(c("STAGE7H_TRUTH_ALIGNMENT_COMPLETE",
             paste0("blind_hash_list_sha256=",
                    sha256(file.path(blind, "BLIND_TABLES.sha256"))),
             paste0("rows=", nrow(out_rows))),
           file.path(out, "TRUTH_ALIGNMENT_COMPLETE.txt"),
           useBytes = TRUE)
cat("STAGE7H_TRUTH_ALIGNMENT_COMPLETE rows=", nrow(out_rows), "\n",
    sep = "")
