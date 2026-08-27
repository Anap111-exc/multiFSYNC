#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
if (!grepl("^([A-Za-z]:|/)", script_file)) {
  script_file <- file.path(getwd(), script_file)
}
root <- dirname(script_file)

r_files <- list.files(root, pattern = "[.]R$", full.names = TRUE)
parse_ok <- vapply(r_files, function(path) {
  tryCatch({ parse(file = path); TRUE }, error = function(condition) FALSE)
}, logical(1L))
fit <- utils::read.csv(file.path(root, "FIT_MANIFEST.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
old <- utils::read.csv(file.path(root, "OLD_ENDPOINT_MANIFEST.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
data <- utils::read.csv(file.path(root, "TARGET_DATA_BINDING.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
qc <- utils::read.csv(file.path(root, "MANIFEST_QC.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
checks <- c(
  all_r_files_parse = all(parse_ok),
  three_target_data = nrow(data) == 3L,
  eighteen_new_fits = nrow(fit) == 18L,
  six_inherited_endpoints = nrow(old) == 6L,
  eight_total_starts_per_data = all(
    table(c(fit$data_id, old$data_id)) == 8L
  ),
  fit_ids_unique = !anyDuplicated(c(fit$fit_id, old$fit_id)),
  fit_seeds_unique = !anyDuplicated(c(fit$fit_seed, old$fit_seed)),
  manifest_qc_passed = nrow(qc) == 10L && all(qc$passed)
)
result <- data.frame(check_id = names(checks), passed = unname(checks),
                     detail = c(
                       paste0(sum(parse_ok), "/", length(parse_ok)),
                       nrow(data), nrow(fit), nrow(old),
                       paste(table(c(fit$data_id, old$data_id)), collapse = ";"),
                       length(unique(c(fit$fit_id, old$fit_id))),
                       length(unique(c(fit$fit_seed, old$fit_seed))),
                       paste0(sum(qc$passed), "/", nrow(qc))
                     ), stringsAsFactors = FALSE)
utils::write.csv(result, file.path(root, "BUNDLE_VALIDATION.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
if (!all(checks)) stop("Stage-4E bundle validation failed.", call. = FALSE)
cat("STAGE4E_BUNDLE_VALIDATION_PASS checks=8/8 R=",
    length(r_files), "\n", sep = "")
