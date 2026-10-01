#!/usr/bin/env Rscript

# Descriptive post-truth analysis for the fixed-seed Stage 7B role rotation.
# This combines the immutable Stage 7A identity assignment with the two new
# assignments.  It never participates in fitting or endpoint selection.

`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    hit <- regexec("^--([^=]+)=(.*)$", arg)
    value <- regmatches(arg, hit)[[1L]]
    assert(length(value) == 3L, paste0("Malformed argument: ", arg))
    out[[value[[2L]]]] <- value[[3L]]
  }
  out
}
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
stage7a_output <- cli$`stage7a-output-root` %||% ""
stage7b_output <- cli$`stage7b-output-root` %||% ""
stage7b_evaluation_name <- cli$`stage7b-evaluation-name` %||%
  "scientific_evaluation_20260921_v2"
assert(dir.exists(stage7a_output) && dir.exists(stage7b_output),
       "Both Stage 7A and Stage 7B output roots are required.")
assert(grepl("^[A-Za-z0-9_-]+$", stage7b_evaluation_name),
       "Invalid Stage 7B evaluation name.")

stage7a_eval <- file.path(stage7a_output,
                          "scientific_evaluation_20260920_v1")
stage7b_eval <- file.path(stage7b_output, stage7b_evaluation_name)
assert(dir.exists(stage7a_eval) && dir.exists(stage7b_eval),
       "A required scientific-evaluation directory is missing.")

read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing file: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv_new <- function(value, name) {
  path <- file.path(analysis_root, name)
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(value, path, row.names = FALSE, na = "")
}
safe_mean <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) NA_real_ else mean(x[is.finite(x)])
}
safe_median <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) NA_real_ else stats::median(x[is.finite(x)])
}
safe_min <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) NA_real_ else min(x[is.finite(x)])
}
safe_max <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) NA_real_ else max(x[is.finite(x)])
}

stage7a_all_full <- read_csv(file.path(stage7a_eval,
                                       "ALL_START_SCIENCE.csv"))
stage7a_winner_full <- read_csv(file.path(stage7a_eval,
                                          "WINNER_SCIENCE.csv"))
stage7a_all <- stage7a_all_full[
  stage7a_all_full$data_id == "s7a_03__sine", , drop = FALSE]
stage7a_winner <- stage7a_winner_full[
  stage7a_winner_full$data_id == "s7a_03__sine", , drop = FALSE]
stage7b_all <- read_csv(file.path(stage7b_eval, "ALL_START_SCIENCE.csv"))
stage7b_winner <- read_csv(file.path(stage7b_eval, "WINNER_SCIENCE.csv"))
stage7b_status <- read_csv(file.path(stage7b_eval, "EVALUATION_STATUS.csv"))

assert(nrow(stage7a_all) == 12L && nrow(stage7a_winner) == 1L &&
         nrow(stage7b_all) == 24L && nrow(stage7b_winner) == 2L &&
         nrow(stage7b_status) == 24L && all(stage7b_status$ok),
       "Unexpected Stage 7A/7B scientific-evaluation dimensions or errors.")
assert(identical(names(stage7a_all), names(stage7b_all)),
       "Stage 7A and Stage 7B science schemas differ.")
assert(all(stage7a_winner$truth_free_winner) &&
         all(stage7b_winner$truth_free_winner),
       "A combined winner was not selected truth-free.")

stage7b_manifest <- read_csv(file.path(stage7b_output, "DATA_MANIFEST.csv"))
assert(nrow(stage7b_manifest) == 2L &&
         setequal(stage7b_manifest$data_id, unique(stage7b_all$data_id)),
       "Stage 7B data manifest does not bind evaluated data.")
rotated_maps <- stage7b_manifest[, c(
  "data_id", "shared_block", "study1_block", "study2_block")]
rotated_maps$assignment <- sub("^s7b_03__", "", rotated_maps$data_id)
rotated_maps <- rotated_maps[, c(
  "data_id", "assignment", "shared_block", "study1_block", "study2_block")]
maps <- rbind(
  data.frame(data_id = "s7a_03__sine", assignment = "identity",
             shared_block = 1L, study1_block = 2L, study2_block = 3L),
  rotated_maps)
map_blocks <- as.matrix(maps[, c(
  "shared_block", "study1_block", "study2_block")])
assert(!anyNA(maps$assignment) && nrow(maps) == 3L &&
         all(apply(map_blocks, 1L, function(x) identical(
           sort(as.integer(x)), 1:3))) &&
         identical(tabulate(as.integer(map_blocks), nbins = 3L),
                   rep(3L, 3L)),
       "The three role maps are not balanced permutations.")

combined_all <- rbind(stage7a_all, stage7b_all)
combined_winner <- rbind(stage7a_winner, stage7b_winner)
combined_all <- merge(combined_all, maps, by = "data_id", sort = FALSE)
combined_winner <- merge(combined_winner, maps, by = "data_id", sort = FALSE)
combined_all <- combined_all[order(match(combined_all$data_id, maps$data_id),
                                   combined_all$start_index), , drop = FALSE]
combined_winner <- combined_winner[
  order(match(combined_winner$data_id, maps$data_id)), , drop = FALSE]
assert(nrow(combined_all) == 36L && nrow(combined_winner) == 3L &&
         all(table(combined_all$data_id) == 12L),
       "Combined Stage 7B audit is incomplete.")

analysis_root <- file.path(stage7b_eval,
                           "descriptive_analysis_20260921_v1")
assert(!dir.exists(analysis_root),
       "Stage 7B descriptive-analysis directory already exists.")
dir.create(analysis_root, recursive = TRUE, mode = "0700")
write_csv_new(combined_all, "ALL_START_SCIENCE_COMBINED.csv")
write_csv_new(combined_winner, "WINNER_SCIENCE_COMBINED.csv")
write_csv_new(maps, "ROLE_MAPS.csv")

dataset_summary <- do.call(rbind, lapply(split(combined_all,
                                               combined_all$data_id), function(x) {
  winner <- x[x$truth_free_winner, , drop = FALSE]
  assert(nrow(winner) == 1L, "Each assignment must have one frozen winner.")
  ordered_elbo <- sort(x$final_ordinary_elbo, decreasing = TRUE)
  data.frame(
    data_id = x$data_id[[1L]], assignment = x$assignment[[1L]],
    shared_block = x$shared_block[[1L]],
    study1_block = x$study1_block[[1L]],
    study2_block = x$study2_block[[1L]],
    starts = nrow(x), core_basin_hits = sum(x$core_truth_recovered),
    exact_structure_hits = sum(x$exact_thresholded_factor_structure),
    full_feature_coverage_hits = sum(x$full_true_feature_component_coverage),
    strict_practical_hits = sum(x$strict_practical_converged),
    winner_fit_id = winner$fit_id[[1L]],
    winner_start_index = winner$start_index[[1L]],
    winner_core_basin = winner$core_truth_recovered[[1L]],
    winner_exact_structure =
      winner$exact_thresholded_factor_structure[[1L]],
    winner_full_feature_coverage =
      winner$full_true_feature_component_coverage[[1L]],
    winner_strict_practical = winner$strict_practical_converged[[1L]],
    winner_obs_signal_nrmse =
      winner$obs_full_signal_nrmse_selected[[1L]],
    winner_dense_signal_nrmse =
      winner$dense_full_signal_nrmse_selected[[1L]],
    winner_loading_relative_l2 =
      winner$loading_relative_l2_error_mean[[1L]],
    winner_feature_total_ise = winner$feature_total_ise_mean[[1L]],
    winner_feature_projection_floor_ise =
      winner$feature_projection_floor_ise_mean[[1L]],
    winner_feature_excess_ise =
      winner$feature_estimate_to_projection_ise_mean[[1L]],
    winner_factor_process_nrmse =
      winner$factor_process_nrmse_mean_missing_as_zero[[1L]],
    winner_complete_contribution_nrmse =
      winner$complete_contribution_nrmse_mean_missing_as_zero[[1L]],
    winner_covariance_operator_relative_error =
      winner$covariance_operator_relative_error_mean[[1L]],
    elbo_margin_to_second = ordered_elbo[[1L]] - ordered_elbo[[2L]])
}))
dataset_summary <- dataset_summary[
  order(match(dataset_summary$data_id, maps$data_id)), , drop = FALSE]
write_csv_new(dataset_summary, "ASSIGNMENT_SUMMARY.csv")

paired_start <- do.call(rbind, lapply(seq_len(12L), function(start) {
  rows <- combined_all[combined_all$start_index == start, , drop = FALSE]
  rows <- rows[match(maps$data_id, rows$data_id), , drop = FALSE]
  assert(nrow(rows) == 3L && !anyNA(rows$data_id),
         "A paired start is missing an assignment.")
  data.frame(
    start_index = start,
    identity_core = rows$core_truth_recovered[[1L]],
    rot231_core = rows$core_truth_recovered[[2L]],
    rot312_core = rows$core_truth_recovered[[3L]],
    identity_exact = rows$exact_thresholded_factor_structure[[1L]],
    rot231_exact = rows$exact_thresholded_factor_structure[[2L]],
    rot312_exact = rows$exact_thresholded_factor_structure[[3L]],
    core_assignments_recovered = sum(rows$core_truth_recovered),
    exact_assignments_recovered =
      sum(rows$exact_thresholded_factor_structure))
}))
write_csv_new(paired_start, "PAIRED_START_BASIN_PATTERNS.csv")
write_csv_new(data.frame(
  paired_starts = nrow(paired_start),
  core_all_three = sum(paired_start$core_assignments_recovered == 3L),
  core_two = sum(paired_start$core_assignments_recovered == 2L),
  core_one = sum(paired_start$core_assignments_recovered == 1L),
  core_none = sum(paired_start$core_assignments_recovered == 0L),
  exact_all_three = sum(paired_start$exact_assignments_recovered == 3L),
  exact_two = sum(paired_start$exact_assignments_recovered == 2L),
  exact_one = sum(paired_start$exact_assignments_recovered == 1L),
  exact_none = sum(paired_start$exact_assignments_recovered == 0L)),
  "PAIRED_START_BASIN_SUMMARY.csv")

role_key <- function(true_role, study) {
  ifelse(true_role == "shared", "shared", paste0("study", study))
}
detail_root <- function(data_id) {
  if (identical(data_id, "s7a_03__sine")) stage7a_eval else stage7b_eval
}
block_for_role <- function(map, role) {
  as.integer(map[[paste0(role, "_block")]][[1L]])
}
extract_fit_roles <- function(science_row) {
  fit_id <- science_row$fit_id[[1L]]
  data_id <- science_row$data_id[[1L]]
  root <- file.path(detail_root(data_id), fit_id)
  loading <- read_csv(file.path(root, "loading_recovery.csv"))
  feature <- read_csv(file.path(root, "feature_component_errors.csv"))
  kernel <- read_csv(file.path(root, "factor_kernel_errors.csv"))
  process <- read_csv(file.path(root, "process_contribution_by_role.csv"))
  map <- maps[maps$data_id == data_id, , drop = FALSE]
  do.call(rbind, lapply(c("shared", "study1", "study2"), function(role) {
    select_role <- function(value) {
      key <- role_key(value$true_role, value$study)
      value[key == role, , drop = FALSE]
    }
    load_key <- ifelse(loading$true_block == "shared", "shared",
                       paste0("study", loading$true_study))
    lr <- loading[load_key == role, , drop = FALSE]
    ft <- select_role(feature)
    ke <- select_role(kernel)
    pr <- select_role(process)
    assert(nrow(lr) == 1L && nrow(ft) >= 2L && nrow(ke) >= 1L &&
             nrow(pr) >= 1L,
           paste0("Malformed role detail for ", fit_id, ":", role))
    data.frame(
      fit_id = fit_id, data_id = data_id,
      assignment = map$assignment[[1L]], start_index = science_row$start_index,
      truth_free_winner = science_row$truth_free_winner,
      core_truth_recovered = science_row$core_truth_recovered,
      exact_thresholded_factor_structure =
        science_row$exact_thresholded_factor_structure,
      role = role, truth_block = paste0("q", block_for_role(map, role)),
      loading_status = lr$status[[1L]],
      loading_matched = lr$matched[[1L]],
      loading_abs_cosine = lr$loading_abs_cosine[[1L]],
      loading_relative_l2_error = lr$loading_relative_l2_error[[1L]],
      loading_norm_ratio = lr$loading_norm_ratio[[1L]],
      feature_components_total = nrow(ft),
      feature_components_matched = sum(ft$component_matched %in% TRUE),
      feature_total_ise_mean = safe_mean(ft$total_sign_aligned_ise),
      feature_projection_floor_ise_mean =
        safe_mean(ft$projection_floor_ise),
      feature_excess_ise_mean = safe_mean(ft$estimate_to_projection_ise),
      kernel_relative_ise_mean = safe_mean(ke$kernel_relative_ise),
      kernel_relative_hs_error_mean =
        safe_mean(ke$kernel_relative_hs_error),
      kernel_projection_floor_relative_ise_mean =
        safe_mean(ke$kernel_projection_floor_relative_ise),
      kernel_excess_relative_ise_mean =
        safe_mean(ke$kernel_estimate_to_projection_relative_ise),
      factor_process_nrmse =
        safe_mean(pr$factor_process_nrmse_missing_as_zero),
      complete_contribution_nrmse =
        safe_mean(pr$complete_contribution_nrmse_missing_as_zero))
  }))
}
role_metrics <- do.call(rbind, lapply(seq_len(nrow(combined_all)), function(i) {
  extract_fit_roles(combined_all[i, , drop = FALSE])
}))
winner_roles <- role_metrics[role_metrics$truth_free_winner, , drop = FALSE]
assert(nrow(role_metrics) == 108L && nrow(winner_roles) == 9L &&
         all(table(role_metrics$truth_block) == 36L) &&
         all(table(role_metrics$role) == 36L) &&
         all(table(winner_roles$truth_block) == 3L) &&
         all(table(winner_roles$role) == 3L),
       "Role/block balance or detailed metric extraction failed.")
write_csv_new(role_metrics, "ALL_START_ROLE_METRICS.csv")
write_csv_new(winner_roles, "WINNER_ROLE_METRICS.csv")

role_metric_names <- c(
  "loading_abs_cosine", "loading_relative_l2_error", "loading_norm_ratio",
  "feature_total_ise_mean", "feature_projection_floor_ise_mean",
  "feature_excess_ise_mean", "kernel_relative_ise_mean",
  "kernel_relative_hs_error_mean",
  "kernel_projection_floor_relative_ise_mean",
  "kernel_excess_relative_ise_mean", "factor_process_nrmse",
  "complete_contribution_nrmse")
long_summary <- function(value, group) {
  groups <- split(value, value[[group]])
  do.call(rbind, lapply(names(groups), function(label) {
    x <- groups[[label]]
    do.call(rbind, lapply(role_metric_names, function(metric) {
      z <- x[[metric]]
      data.frame(group = group, level = label, metric = metric,
                 rows = nrow(x), finite_rows = sum(is.finite(z)),
                 mean = safe_mean(z), median = safe_median(z),
                 min = safe_min(z), max = safe_max(z))
    }))
  }))
}
write_csv_new(long_summary(role_metrics, "truth_block"),
              "ALL_START_SUMMARY_BY_TRUTH_BLOCK.csv")
write_csv_new(long_summary(role_metrics, "role"),
              "ALL_START_SUMMARY_BY_ROLE.csv")
write_csv_new(long_summary(winner_roles, "truth_block"),
              "WINNER_SUMMARY_BY_TRUTH_BLOCK.csv")
write_csv_new(long_summary(winner_roles, "role"),
              "WINNER_SUMMARY_BY_ROLE.csv")

writeLines(c(
  "status=STAGE7B_POST_TRUTH_DESCRIPTIVE_ANALYSIS_COMPLETE",
  "stage7a_identity_starts=12", "stage7b_rotated_starts=24",
  "combined_starts=36", "truth_free_winners=3",
  "role_metric_rows=108", "winner_role_metric_rows=9",
  "selection_frozen_before_truth_read=TRUE",
  "truth_used_for_fit_or_selection=FALSE",
  "new_fit_or_continuation_started=FALSE",
  "raw_elbo_compared_across_datasets=FALSE"),
  file.path(analysis_root, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE7B_DESCRIPTIVE_ANALYSIS_COMPLETE starts=36 winners=3\n")
