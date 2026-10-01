#!/usr/bin/env Rscript

# Add the study-level covariance-operator attribution that is not present in
# the compact Stage 7A/7B science rows.  This reads only completed evaluations.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, all(dir.exists(args[1:2])))
stage7a_output <- args[[1L]]
stage7b_output <- args[[2L]]
stage7b_evaluation_name <- args[[3L]]
stage7a_eval <- file.path(stage7a_output,
                          "scientific_evaluation_20260920_v1")
stage7b_eval <- file.path(stage7b_output, stage7b_evaluation_name)
read_csv <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE)

stage7a_winner <- read_csv(file.path(stage7a_eval, "WINNER_SCIENCE.csv"))
stage7a_winner <- stage7a_winner[
  stage7a_winner$data_id == "s7a_03__sine", , drop = FALSE]
stage7b_winner <- read_csv(file.path(stage7b_eval, "WINNER_SCIENCE.csv"))
winners <- rbind(stage7a_winner, stage7b_winner)
stopifnot(nrow(winners) == 3L, all(winners$truth_free_winner))

manifest <- read_csv(file.path(stage7b_output, "DATA_MANIFEST.csv"))
rotated <- manifest[, c(
  "data_id", "shared_block", "study1_block", "study2_block")]
rotated$assignment <- sub("^s7b_03__", "", rotated$data_id)
maps <- rbind(
  data.frame(data_id = "s7a_03__sine", shared_block = 1L,
             study1_block = 2L, study2_block = 3L,
             assignment = "identity"),
  rotated)

rows <- do.call(rbind, lapply(seq_len(nrow(winners)), function(i) {
  winner <- winners[i, , drop = FALSE]
  data_id <- winner$data_id[[1L]]
  fit_id <- winner$fit_id[[1L]]
  eval_root <- if (data_id == "s7a_03__sine") stage7a_eval else stage7b_eval
  value <- read_csv(file.path(eval_root, fit_id,
                              "covariance_operator_error.csv"))
  map <- maps[maps$data_id == data_id, , drop = FALSE]
  stopifnot(nrow(value) == 2L, nrow(map) == 1L,
            identical(as.integer(value$study), 1:2))
  value$fit_id <- fit_id
  value$data_id <- data_id
  value$assignment <- map$assignment[[1L]]
  value$shared_block <- paste0("q", map$shared_block[[1L]])
  value$specific_block <- paste0("q", c(
    map$study1_block[[1L]], map$study2_block[[1L]]))
  value[, c("fit_id", "data_id", "assignment", "study",
            "shared_block", "specific_block", "truth_operator_hs_norm",
            "estimated_operator_hs_norm",
            "covariance_operator_relative_error")]
}))
stopifnot(nrow(rows) == 6L)
analysis_root <- file.path(stage7b_eval,
                           "descriptive_analysis_20260921_v1")
path <- file.path(analysis_root, "WINNER_STUDY_COVARIANCE_METRICS.csv")
stopifnot(dir.exists(analysis_root), !file.exists(path))
utils::write.csv(rows, path, row.names = FALSE, na = "")
writeLines(c(
  "status=STAGE7B_WINNER_STUDY_COVARIANCE_ATTRIBUTION_COMPLETE",
  "winner_rows=3", "study_rows=6", "post_truth_only=TRUE",
  "new_fit_or_continuation_started=FALSE"),
  file.path(analysis_root, "COVARIANCE_ATTRIBUTION_COMPLETE.txt"),
  useBytes = TRUE)
cat("STAGE7B_COVARIANCE_ATTRIBUTION_COMPLETE rows=6\n")
