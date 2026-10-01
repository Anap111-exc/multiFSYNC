#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
experiment_dir <- dirname(file_arg)
if (!dir.exists(experiment_dir)) stop("Experiment directory does not exist.")

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required to bind evidence with SHA256.")
}

manifest_path <- file.path(experiment_dir, "EVIDENCE_MANIFEST_V1.csv")
status_path <- file.path(
  experiment_dir, "FROZEN_EVIDENCE_STATUS_V1_20260927.txt")
if (file.exists(manifest_path) || file.exists(status_path)) {
  stop("Versioned freeze outputs already exist; refusing to overwrite them.")
}

evaluation_dir <- file.path(
  experiment_dir, "scientific_evaluation_20260927_v1")
summary_dir <- file.path(
  evaluation_dir, "descriptive_summary_20260927_v1")

all_start <- read.csv(
  file.path(evaluation_dir, "ALL_START_SCIENCE.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
winners <- read.csv(
  file.path(evaluation_dir, "WINNER_SCIENCE.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
evaluation_status <- read.csv(
  file.path(evaluation_dir, "EVALUATION_STATUS.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(
  nrow(all_start) == 36L,
  nrow(winners) == 3L,
  nrow(evaluation_status) == 36L,
  all(evaluation_status$ok),
  all(all_start$objective_eligible),
  all(!all_start$truth_used_for_fit_stopping_or_selection),
  all(winners$truth_free_winner),
  identical(sort(unique(all_start$beta_arm)), c("C0", "C1", "C2")),
  all(table(all_start$beta_arm) == 12L),
  all(table(winners$beta_arm) == 1L))

complete_lines <- readLines(
  file.path(evaluation_dir, "EVALUATION_COMPLETE.txt"), warn = FALSE)
required_markers <- c(
  "status=COVARIATE_ROUTEB_SCIENTIFIC_EVALUATION_COMPLETE",
  "fits_evaluated=36",
  "truth_free_winners_evaluated=3",
  paste0("frozen_winners_sha256=",
         "501dcdf11f1ef2ffad15db02fd3176f00f74ee6166fea3936228ee94a1dbee6e"),
  "truth_read_after_truth_free_selection_freeze=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE")
stopifnot(all(required_markers %in% complete_lines))

relative_paths <- c(
  "DGP_AND_BOUNDARY_SPECIFICATION_V1_20260927.md",
  "PROTOCOL_COVARIATE_ROUTEB_PILOT_V1_DRAFT_20260927.md",
  "FIT_MANIFEST_CANDIDATE_V1.csv",
  "SOURCE_BINDING_V1.csv",
  "EXECUTION_SOURCE_BINDING_V1.csv",
  "PILOT_IMPLEMENTATION_AND_EXECUTION_AUTHORIZATION_20260927.txt",
  "TRUTH_UNSEAL_AND_SCIENTIFIC_EVALUATION_AUTHORIZATION_20260927.txt",
  "run_covariate_routeb_pilot_20260927_v1.R",
  "test_covariate_routeb_pilot_20260927_v1.R",
  "launch_covariate_routeb_pilot_20260927_v1.sh",
  "LOCAL_AND_LINUX_ACCEPTANCE_20260927.md",
  "evaluate_covariate_routeb_pilot_20260927_v1.R",
  "summarize_covariate_routeb_pilot_20260927_v1.R",
  "REPORT_COVARIATE_ROUTEB_PILOT_SCIENTIFIC_EVALUATION_20260927.md",
  "EVIDENCE_FREEZE_NOTE_V1_20260927.md",
  "freeze_covariate_routeb_pilot_evidence_20260927_v1.R",
  "covariate-routeb-scientific-evaluation-20260927-v1.tar.gz",
  "scientific_evaluation_20260927_v1/EVALUATION_COMPLETE.txt",
  "scientific_evaluation_20260927_v1/EVALUATION_STATUS.csv",
  "scientific_evaluation_20260927_v1/ALL_START_SCIENCE.csv",
  "scientific_evaluation_20260927_v1/WINNER_SCIENCE.csv",
  "scientific_evaluation_20260927_v1/ARM_SUMMARY.csv",
  "scientific_evaluation_20260927_v1/PAIRED_START_C1_C2.csv",
  "scientific_evaluation_20260927_v1/TRUTH_UNSEAL_AUTHORIZATION.txt",
  paste0("scientific_evaluation_20260927_v1/",
         "descriptive_summary_20260927_v1/",
         c("ALL_START_METRIC_SUMMARY_BY_ARM.csv",
           "ARM_BASIN_AND_CONVERGENCE.csv",
           "BETA_ARM_SUMMARY.csv",
           "DESCRIPTIVE_SUMMARY_COMPLETE.txt",
           "ELBO_SCIENCE_ASSOCIATION.csv",
           "PAIRED_START_C1_C2_DIFFERENCES.csv",
           "PAIRED_START_C1_C2_SUMMARY.csv",
           "RUNTIME_SUMMARY.csv",
           "WINNER_SCIENCE_PRESENTATION.csv")))

roles <- c(
  "dgp_specification", "pilot_protocol", "fit_manifest",
  "source_binding", "execution_source_binding", "execution_authorization",
  "truth_unseal_authorization", "fit_runner", "regression_test",
  "linux_launcher", "acceptance_record", "scientific_evaluator",
  "scientific_summarizer", "scientific_report", "freeze_note",
  "freeze_script", "complete_evaluation_archive",
  "evaluation_complete_marker", "evaluation_status",
  "all_start_science", "truth_free_winner_science", "arm_summary",
  "paired_start_science", "copied_truth_unseal_authorization",
  rep("descriptive_summary", 9L))
stopifnot(length(relative_paths) == length(roles))

absolute_paths <- file.path(experiment_dir, relative_paths)
stopifnot(all(file.exists(absolute_paths)))

manifest <- data.frame(
  evidence_id = "COVARIATE_ROUTEB_PILOT_EVIDENCE_V1_20260927",
  role = roles,
  relative_path = gsub("\\\\", "/", relative_paths),
  bytes = as.numeric(file.info(absolute_paths)$size),
  sha256 = vapply(
    absolute_paths,
    digest::digest,
    character(1L), file = TRUE, algo = "sha256", serialize = FALSE),
  stringsAsFactors = FALSE)
write.csv(manifest, manifest_path, row.names = FALSE, na = "")

archive_row <- manifest[
  manifest$role == "complete_evaluation_archive", , drop = FALSE]
stopifnot(
  nrow(archive_row) == 1L,
  identical(
    archive_row$sha256,
    "33c4109c595e1bbeec1d12ab0b31031926b93608354ae8b8df36e41fe0ff83ed"))

manifest_sha <- digest::digest(
  manifest_path, file = TRUE, algo = "sha256", serialize = FALSE)
status_lines <- c(
  "status=FROZEN_DEVELOPMENT_EVIDENCE",
  "evidence_id=COVARIATE_ROUTEB_PILOT_EVIDENCE_V1_20260927",
  "logical_data_seeds=1",
  "arms=C0;C1;C2",
  "fits_evaluated=36",
  "objective_eligible_fits=36",
  "truth_free_winners=3",
  "evaluation_errors=0",
  "truth_free_selection_preceded_truth_unseal=TRUE",
  "truth_used_for_fit_stopping_or_selection=FALSE",
  "new_fit_or_continuation_started_by_freeze=FALSE",
  "paper_level_monte_carlo=FALSE",
  "raw_fit_objects_duplicated_locally=FALSE",
  paste0("raw_fit_server_directory=",
         "/root/multifsync-results/covariate_routeb_pilot_20260927_v1"),
  paste0("frozen_winners_sha256=",
         "501dcdf11f1ef2ffad15db02fd3176f00f74ee6166fea3936228ee94a1dbee6e"),
  paste0("evaluation_archive_sha256=", archive_row$sha256),
  paste0("evidence_manifest_rows=", nrow(manifest)),
  paste0("evidence_manifest_sha256=", manifest_sha),
  "git_commit_or_tag_created=FALSE",
  "stage7c_execution_authorized=FALSE",
  "stage7c_data_generated=FALSE",
  "stage7c_fits_started=FALSE")
writeLines(status_lines, status_path, useBytes = TRUE)

cat("Frozen evidence manifest rows:", nrow(manifest), "\n")
cat("Manifest SHA256:", manifest_sha, "\n")
cat("Status:", status_path, "\n")
