#!/usr/bin/env Rscript

# Apply the exact Stage 7C runner convergence label to saved Stage 7E fits.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript SCRIPT STAGE7E_ROOT")
root <- normalizePath(args[[1L]], mustWork = TRUE)
output_root <- file.path(root, "convergence_semantics_v1")
if (file.exists(output_root) || dir.exists(output_root)) {
  stop("Refusing to overwrite convergence_semantics_v1.")
}
endpoints <- read.csv(file.path(root, "evaluation_v2",
                                "ORACLE_ENDPOINT_COMPARISON.csv"))
stopifnot(nrow(endpoints) == 6L,
          file.exists(file.path(root, "FITS_COMPLETE.txt")))
rows <- lapply(seq_len(nrow(endpoints)), function(index) {
  row <- endpoints[index, ]
  fit <- readRDS(file.path(root, "fits", row$fit_id, "fit.rds"))
  stage7c_strict <- isTRUE(fit$practical_converged) &&
    isTRUE(fit$converged) &&
    identical(as.character(fit$convergence_status),
              "converged_practical")
  data.frame(
    fit_id = row$fit_id,
    practical_converged = isTRUE(fit$practical_converged),
    converged = isTRUE(fit$converged),
    convergence_status = as.character(fit$convergence_status),
    stage7c_strict_practical_converged = stage7c_strict,
    endpoint_table_strict = row$strict_practical_converged,
    fields_agree = identical(stage7c_strict,
                             row$strict_practical_converged),
    stringsAsFactors = FALSE)
})
result <- do.call(rbind, rows)
stopifnot(all(result$fields_agree))
dir.create(output_root)
write.csv(result, file.path(output_root, "STRICT_STATUS.csv"),
          row.names = FALSE)
writeLines(c("STAGE7E_STRICT_STATUS_VERIFIED", "fits=6",
             "criterion=practical_converged_and_converged_and_convergence_status_converged_practical",
             "new_fit_started=FALSE"),
           file.path(output_root, "STRICT_STATUS_VERIFIED.txt"))
cat("STAGE7E_STRICT_STATUS_VERIFIED fits=6\n")
