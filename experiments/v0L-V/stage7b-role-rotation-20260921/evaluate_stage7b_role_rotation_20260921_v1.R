#!/usr/bin/env Rscript

# Post-selection Stage 7B scientific evaluation.  Stage 7B deliberately uses
# the already-audited Stage 7A evaluator without copying or changing any
# scientific formula.  The substitutions below only bind the evaluator to the
# smaller Stage 7B manifest and its sealed-truth identity.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
script_root <- dirname(script_path)
template_path <- normalizePath(file.path(
  script_root, "..", "covariate-effect-stage7a-20260919",
  "evaluate_covariate_stage7a_20260920_v1.R"),
  winslash = "/", mustWork = TRUE)

stopifnot(requireNamespace("digest", quietly = TRUE))
template_sha256 <- tolower(digest::digest(
  file = template_path, algo = "sha256", serialize = FALSE))
stopifnot(identical(
  template_sha256,
  "b6102710a692e8e14a2dbb9cc5a500c7fcad81b5d3ae631d83a06d797f88d97d"
))

code <- readLines(template_path, warn = FALSE, encoding = "UTF-8")
replace_one <- function(old, new) {
  hits <- which(grepl(old, code, fixed = TRUE))
  stopifnot(length(hits) == 1L)
  code[hits] <<- sub(old, new, code[hits], fixed = TRUE)
}

replace_one('limit <- as.integer(cli$limit %||% "144")',
            'limit <- as.integer(cli$limit %||% "24")')
replace_one('"scientific_evaluation_20260920_v1"',
            '"scientific_evaluation_20260921_v1"')
replace_one('workers %in% 1:8 && limit %in% 1:144 &&',
            'workers %in% 1:8 && limit %in% 1:24 &&')
replace_one('"selected_winners=12", "registered_fits=144",',
            '"selected_winners=2", "registered_fits=24",')
replace_one('nrow(manifest) == 144L && nrow(data_manifest) == 12L &&',
            'nrow(manifest) == 24L && nrow(data_manifest) == 2L &&')
replace_one('nrow(endpoints) == 144L && nrow(winners) == 12L &&',
            'nrow(endpoints) == 24L && nrow(winners) == 2L &&')
replace_one('authorization_date=2026-09-20',
            'authorization_date=2026-09-21')
replace_one('identical(x$bundle_class, "stage7a_sealed_truth") &&',
            'identical(x$bundle_class, "stage7b_sealed_truth") &&')
replace_one('"MULTIFSYNC_COVARIATE_EFFECT_STAGE7A_V1_20260919")',
            '"MULTIFSYNC_STAGE7B_ROLE_ROTATION_V1_20260921")')
replace_one('phase = "stage7a_post_selection_scientific_evaluation"',
            'phase = "stage7b_post_selection_scientific_evaluation"')
replace_one('if (limit == 144L) {', 'if (FALSE) {')
replace_one('status=STAGE7A_POST_SELECTION_SCIENTIFIC_EVALUATION_COMPLETE',
            'status=STAGE7B_POST_SELECTION_SCIENTIFIC_EVALUATION_COMPLETE')
replace_one('cat("STAGE7A_EVALUATION_COMPLETE fits="',
            'cat("STAGE7B_EVALUATION_COMPLETE fits="')

# Parsing the transformed, SHA-bound evaluator avoids maintaining a second
# implementation of loading, feature, process, contribution, covariance, and
# covariate metrics.  The fit and winner selection files are read-only here.
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
