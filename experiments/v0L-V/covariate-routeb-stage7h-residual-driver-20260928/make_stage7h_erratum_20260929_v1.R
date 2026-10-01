#!/usr/bin/env Rscript

# Versioned correction of two Stage 7H presentation issues. No fit or truth
# selection occurs here; all inputs are previously frozen machine tables.
options(warn = 2, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript make_stage7h_erratum_20260929_v1.R STAGE7H_ROOT OUTPUT_DIR",
       call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
out <- args[[2L]]
if (dir.exists(out) || file.exists(out)) {
  stop("Refusing existing output path: ", out, call. = FALSE)
}

truth_path <- file.path(root, "stage7h_20260929_v1", "truth_alignment_v1",
                        "RESIDUAL_TRUTH_ALIGNMENT.csv")
driver_path <- file.path(root, "stage7h_20260929_v1", "blind_all_v3",
                         "UPDATE_DRIVER_BLIND.csv")
block_path <- file.path(root, "stage7h_20260929_v1", "blind_all_v3",
                        "BLOCK_ENERGY_BLIND.csv")
if (!all(file.exists(c(truth_path, driver_path, block_path)))) {
  stop("A frozen Stage 7H input table is missing.", call. = FALSE)
}

# Keep the exact text of every valid field. Split singular values and gaps in
# v1 came from the full sample; those four split-specific values are unknown.
old <- utils::read.csv(truth_path, colClasses = "character",
                       check.names = FALSE, na.strings = character(),
                       fileEncoding = "UTF-8")
if (nrow(old) != 480L) stop("Unexpected truth alignment row count.")
key <- c("fit_id", "iteration", "residual_stage", "domain", "split",
         "true_direction")
spectral <- c("sigma1", "sigma2", "gap1", "gap5")
if (!all(c(key, spectral) %in% names(old)) ||
    anyDuplicated(do.call(paste, c(old[key], sep = "\r")))) {
  stop("Truth alignment schema or key changed.")
}
is_split <- old$split %in% c("odd", "even")
if (sum(is_split) != 160L || sum(old$split == "full") != 320L ||
    any(!old$split %in% c("full", "odd", "even"))) {
  stop("Unexpected split allocation.")
}
fixed <- old
fixed[is_split, spectral] <- NA_character_
fixed$spectrum_scope <- ifelse(is_split, "not_recorded_for_split",
                               "full_sample")
if (!identical(fixed[setdiff(names(old), spectral)],
               old[setdiff(names(old), spectral)]) ||
    !identical(fixed[!is_split, spectral], old[!is_split, spectral]) ||
    !all(is.na(as.matrix(fixed[is_split, spectral])))) {
  stop("Correction changed a valid v1 value.")
}

# U_mean is the sum over subjects of ||C_i psi xi_i||^2 at observed times.
# Normalize the 500-coordinate linear loading driver by sqrt(U_mean).
driver <- utils::read.csv(driver_path, check.names = FALSE)
block <- utils::read.csv(block_path, check.names = FALSE)
driver <- driver[driver$iteration == 50L, , drop = FALSE]
block <- block[block$iteration == 50L, , drop = FALSE]
fits <- c("covrb7c_01__C1__start06", "covrb7c_02__C1__start11",
          "covrb7c_03__C1__start06")
if (nrow(driver) != 1500L || nrow(block) != 3L ||
    !identical(sort(unique(driver$fit_id)), sort(fits))) {
  stop("Sweep 50 driver table changed.")
}
rows <- lapply(fits, function(id) {
  x <- driver[driver$fit_id == id, , drop = FALSE]
  y <- block[block$fit_id == id, , drop = FALSE]
  u <- unique(x$U_mean)
  if (nrow(x) != 500L || nrow(y) != 1L || length(u) != 1L ||
      !is.finite(u) || u <= 0 ||
      !isTRUE(all.equal(u, y$U_mean, tolerance = 1e-12))) {
    stop("Inconsistent U_mean or coordinates for ", id)
  }
  h_norm <- sqrt(sum(x$h_total^2))
  data.frame(fit_id = id, iteration = 50L,
             h_norm = h_norm, U_mean = u,
             normalized_h = h_norm / sqrt(u),
             mean_ppi = mean(x$ppi), stringsAsFactors = FALSE)
})
normalized <- do.call(rbind, rows)
if (any(!is.finite(as.matrix(normalized[, -1L]))) ||
    max(abs(normalized$normalized_h -
            c(83.8931320187, 79.9450610813, 86.4947239580))) > 1e-6) {
  stop("Normalized driver regression values changed.")
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(fixed, file.path(out, "RESIDUAL_TRUTH_ALIGNMENT_SCOPE_CORRECTED.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(normalized,
                 file.path(out, "NORMALIZED_DRIVER_SWEEP50.csv"),
                 row.names = FALSE)

# Round-trip checks guard the output contract as well as the in-memory values.
roundtrip <- utils::read.csv(
  file.path(out, "RESIDUAL_TRUTH_ALIGNMENT_SCOPE_CORRECTED.csv"),
  colClasses = "character", check.names = FALSE, na.strings = "",
  fileEncoding = "UTF-8")
if (nrow(roundtrip) != 480L ||
    !identical(roundtrip[!is_split, spectral], old[!is_split, spectral]) ||
    !all(is.na(as.matrix(roundtrip[is_split, spectral]))) ||
    !identical(roundtrip[setdiff(names(old), spectral)],
               old[setdiff(names(old), spectral)])) {
  stop("Corrected CSV round-trip verification failed.")
}
writeLines(c("STAGE7H_ERRATUM_V2_COMPLETE", "truth_rows=480",
             "split_rows_with_spectrum_unavailable=160",
             "normalized_driver_rows=3",
             "original_v1_tables_preserved=true"),
           file.path(out, "ERRATUM_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7H_ERRATUM_V2_COMPLETE\n")
