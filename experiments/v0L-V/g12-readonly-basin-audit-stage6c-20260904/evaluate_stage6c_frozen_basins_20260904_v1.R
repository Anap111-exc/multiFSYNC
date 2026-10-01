#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
parse_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  value <- arguments[startsWith(arguments, prefix)]
  if (!length(value)) return("")
  sub(prefix, "", value[[1L]], fixed = TRUE)
}

stage6a_root <- parse_argument("stage6a-root")
stage6b_root <- parse_argument("stage6b-root")
truth_free_root <- parse_argument("truth-free-root")
output_root <- parse_argument("output-root")
tools_path <- "stage6c_basin_tools_20260904_v1.R"
if (!file.exists(tools_path)) {
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                          value = TRUE)[[1L]])
  tools_path <- file.path(dirname(script_path), basename(tools_path))
}
source(tools_path, local = FALSE)

s6c_assert(all(nzchar(c(stage6a_root, stage6b_root, truth_free_root,
                       output_root))),
             "All four root arguments are required.")
s6c_verify_truth_free_freeze(truth_free_root)
s6c_assert(!dir.exists(output_root) ||
             !length(list.files(output_root, all.files = TRUE, no.. = TRUE)),
           "Stage 6C truth-evaluation output directory already contains files.")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  s6c_assert(file.exists(path), paste0("Missing evaluation input: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

stage6a_ids <- sort(list.dirs(file.path(stage6a_root, "evaluation"),
                              recursive = FALSE, full.names = FALSE))
stage6a_ids <- stage6a_ids[nzchar(stage6a_ids)]
s6c_assert(length(stage6a_ids) == 36L,
           "Expected 36 Stage 6A evaluation directories.")
stage6a <- do.call(rbind, lapply(stage6a_ids, function(fit_id) {
  summary <- read_csv(file.path(stage6a_root, "evaluation", fit_id,
                                "summary.csv"))
  record <- readRDS(file.path(stage6a_root, "fits", fit_id,
                              "terminal_record.rds"))
  s6c_assert(nrow(summary) == 1L && identical(record$fit_id, fit_id),
               "Invalid Stage 6A evaluation record.")
  data.frame(
    cohort_id = "s6a_full_data",
    config_group_id = paste("s6a_full_data", record$fit_config_id, sep = "::"),
    fit_id = fit_id,
    fit_config_id = record$fit_config_id,
    seed_index = record$seed_index,
    practical_converged = record$practical_converged,
    summary,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(stage6a) <- NULL

stage6b <- read_csv(file.path(stage6b_root, "evaluation",
                              "ALL_STAGE6B_SCIENTIFIC_RESULTS.csv"))
s6c_assert(nrow(stage6b) == 20L,
           "Expected 20 Stage 6B scientific evaluation rows.")
stage6b$cohort_id <- "s6b_holdout_training"
stage6b$config_group_id <- paste(stage6b$cohort_id,
                                 stage6b$fit_config_id, sep = "::")

common <- intersect(names(stage6a), names(stage6b))
required <- c(
  "cohort_id", "config_group_id", "fit_id", "fit_config_id", "seed_index",
  "practical_converged", "final_elbo", "selected_counts",
  "global_correct", "global_misplaced", "global_missing", "global_extra",
  "global_duplicate", "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate",
  "shared_effective_M", "specific_effective_M_study_1",
  "specific_effective_M_study_2", "signal_nrmse_ppi_selected",
  "dense_signal_nrmse_ppi_selected", "feature_total_ise_mean",
  "covariance_operator_relative_error_mean"
)
s6c_assert(all(required %in% common),
           "Stage 6A and 6B evaluation schemas lack required common fields.")
scientific <- rbind(stage6a[, required, drop = FALSE],
                    stage6b[, required, drop = FALSE])
rownames(scientific) <- NULL

all_components_equal <- function(value, target) {
  components <- suppressWarnings(as.integer(strsplit(
    as.character(value), ";", fixed = TRUE
  )[[1L]]))
  length(components) >= 1L && all(is.finite(components)) &&
    all(components == target)
}
scientific$core_truth_recovered <- with(
  scientific,
  global_correct == 3L & global_misplaced == 0L & global_missing == 0L &
    study_correct == 4L & study_misplaced == 0L & study_missing == 0L
)
scientific$exact_factor_structure_recovered <- with(
  scientific,
  core_truth_recovered & global_extra == 0L & global_duplicate == 0L &
    study_candidate_extra == 0L & study_candidate_duplicate == 0L &
    selected_counts == "1,1,1"
)
scientific$true_effective_M <- vapply(seq_len(nrow(scientific)), function(i) {
  all_components_equal(scientific$shared_effective_M[[i]], 2L) &&
    all_components_equal(scientific$specific_effective_M_study_1[[i]], 2L) &&
    all_components_equal(scientific$specific_effective_M_study_2[[i]], 2L)
}, logical(1L))
scientific$exact_joint_dimension_and_structure <-
  scientific$exact_factor_structure_recovered & scientific$true_effective_M

assignments <- read_csv(file.path(truth_free_root,
                                  "BASIN_ASSIGNMENTS_TRUTH_FREE.csv"))
basin_summary <- read_csv(file.path(truth_free_root,
                                    "BASIN_SUMMARY_TRUTH_FREE.csv"))
s6c_assert(nrow(scientific) == 56L &&
             setequal(scientific$fit_id, unique(assignments$fit_id)),
           "Truth evaluation and frozen endpoints do not join exactly.")

joined <- merge(assignments, scientific,
                by = c("cohort_id", "config_group_id", "fit_id",
                       "fit_config_id", "seed_index"),
                all.x = TRUE, sort = FALSE, suffixes = c("_audit", "_truth"))
s6c_assert(nrow(joined) == nrow(assignments) &&
             !anyNA(joined$core_truth_recovered),
           "Frozen assignment to truth evaluation join failed.")

cluster_groups <- split(
  joined,
  interaction(joined$config_group_id, joined$basin_view,
              joined$similarity_threshold, joined$basin_id, drop = TRUE)
)
basin_truth <- do.call(rbind, lapply(cluster_groups, function(rows) {
  representative <- rows[rows$is_basin_representative, , drop = FALSE]
  s6c_assert(nrow(representative) == 1L,
             "A frozen basin has no unique representative.")
  data.frame(
    cohort_id = rows$cohort_id[[1L]],
    config_group_id = rows$config_group_id[[1L]],
    fit_config_id = rows$fit_config_id[[1L]],
    basin_view = rows$basin_view[[1L]],
    similarity_threshold = rows$similarity_threshold[[1L]],
    basin_id = rows$basin_id[[1L]],
    basin_size = rows$basin_size[[1L]],
    basin_support = rows$basin_support[[1L]],
    is_modal_basin = rows$is_modal_basin[[1L]],
    representative_fit_id = representative$fit_id,
    representative_core_truth = representative$core_truth_recovered,
    representative_exact_factor_structure =
      representative$exact_factor_structure_recovered,
    representative_true_effective_M = representative$true_effective_M,
    representative_exact_joint =
      representative$exact_joint_dimension_and_structure,
    any_core_truth_in_basin = any(rows$core_truth_recovered),
    any_exact_factor_structure_in_basin =
      any(rows$exact_factor_structure_recovered),
    any_exact_joint_in_basin = any(rows$exact_joint_dimension_and_structure),
    core_truth_rate = mean(rows$core_truth_recovered),
    exact_factor_structure_rate =
      mean(rows$exact_factor_structure_recovered),
    exact_joint_rate = mean(rows$exact_joint_dimension_and_structure),
    dense_signal_nrmse_mean = mean(rows$dense_signal_nrmse_ppi_selected),
    dense_signal_nrmse_min = min(rows$dense_signal_nrmse_ppi_selected),
    feature_ise_mean = mean(rows$feature_total_ise_mean),
    covariance_error_mean = mean(
      rows$covariance_operator_relative_error_mean
    ),
    stringsAsFactors = FALSE
  )
}))
rownames(basin_truth) <- NULL

selector_rows <- list()
summary_groups <- split(
  basin_summary,
  interaction(basin_summary$config_group_id, basin_summary$basin_view,
              basin_summary$similarity_threshold, drop = TRUE)
)
for (rows in summary_groups) {
  s6c_assert(nrow(rows) == 1L, "Frozen basin summary key is not unique.")
  candidates <- c(raw_max_elbo = rows$raw_elbo_winner_fit_id,
                  modal_basin_then_elbo = rows$modal_basin_then_elbo_fit_id)
  for (method in names(candidates)) {
    selected <- scientific[scientific$fit_id == candidates[[method]],
                           , drop = FALSE]
    s6c_assert(nrow(selected) == 1L, "Selected endpoint has no truth row.")
    group_scientific <- scientific[
      scientific$config_group_id == rows$config_group_id, , drop = FALSE
    ]
    selector_rows[[length(selector_rows) + 1L]] <- data.frame(
      cohort_id = rows$cohort_id,
      config_group_id = rows$config_group_id,
      fit_config_id = rows$fit_config_id,
      basin_view = rows$basin_view,
      similarity_threshold = rows$similarity_threshold,
      selector = method,
      selected_fit_id = selected$fit_id,
      selected_core_truth = selected$core_truth_recovered,
      selected_exact_factor_structure =
        selected$exact_factor_structure_recovered,
      selected_true_effective_M = selected$true_effective_M,
      selected_exact_joint = selected$exact_joint_dimension_and_structure,
      selected_dense_signal_nrmse = selected$dense_signal_nrmse_ppi_selected,
      selected_feature_ise = selected$feature_total_ise_mean,
      selected_covariance_error =
        selected$covariance_operator_relative_error_mean,
      any_core_truth_available = any(group_scientific$core_truth_recovered),
      any_exact_factor_structure_available =
        any(group_scientific$exact_factor_structure_recovered),
      any_exact_joint_available =
        any(group_scientific$exact_joint_dimension_and_structure),
      stringsAsFactors = FALSE
    )
  }
}
selector_results <- do.call(rbind, selector_rows)
rownames(selector_results) <- NULL

selector_groups <- split(
  selector_results,
  interaction(selector_results$basin_view,
              selector_results$similarity_threshold,
              selector_results$selector, drop = TRUE)
)
selector_summary <- do.call(rbind, lapply(selector_groups, function(rows) {
  data.frame(
    basin_view = rows$basin_view[[1L]],
    similarity_threshold = rows$similarity_threshold[[1L]],
    selector = rows$selector[[1L]],
    strata = nrow(rows),
    strata_with_core_available = sum(rows$any_core_truth_available),
    selected_core_truth = sum(rows$selected_core_truth),
    strata_with_exact_factor_available =
      sum(rows$any_exact_factor_structure_available),
    selected_exact_factor_structure =
      sum(rows$selected_exact_factor_structure),
    strata_with_exact_joint_available = sum(rows$any_exact_joint_available),
    selected_exact_joint = sum(rows$selected_exact_joint),
    dense_signal_nrmse_mean = mean(rows$selected_dense_signal_nrmse),
    feature_ise_mean = mean(rows$selected_feature_ise),
    covariance_error_mean = mean(rows$selected_covariance_error),
    stringsAsFactors = FALSE
  )
}))
rownames(selector_summary) <- NULL
selector_summary <- selector_summary[order(
  selector_summary$basin_view, selector_summary$similarity_threshold,
  selector_summary$selector
), , drop = FALSE]

primary <- selector_results[
  selector_results$basin_view == "contribution" &
    abs(selector_results$similarity_threshold - 0.9) < 1e-12,
  , drop = FALSE
]
primary_wide <- merge(
  primary[primary$selector == "raw_max_elbo", ],
  primary[primary$selector == "modal_basin_then_elbo", ],
  by = c("cohort_id", "config_group_id", "fit_config_id", "basin_view",
         "similarity_threshold"), suffixes = c("_raw_elbo", "_modal_basin"),
  sort = FALSE
)

integrity <- data.frame(
  check_id = c(
    "truth_free_hashes_verified", "scientific_rows_56",
    "assignment_join_complete", "selector_strata_complete",
    "no_fit_or_continuation_started"
  ),
  passed = c(TRUE, nrow(scientific) == 56L,
             nrow(joined) == nrow(assignments),
             all(table(interaction(selector_results$basin_view,
                                   selector_results$similarity_threshold,
                                   selector_results$selector)) == 8L),
             TRUE),
  detail = c(
    "all frozen truth-free CSV SHA256 values match",
    nrow(scientific), nrow(joined),
    "8 strata for every view-threshold-selector combination",
    "evaluation script contains no fitting or continuation call"
  ),
  stringsAsFactors = FALSE
)
s6c_assert(all(integrity$passed), "Stage 6C truth evaluation integrity failed.")

s6c_write_csv(scientific, file.path(output_root,
                                    "SCIENTIFIC_ENDPOINT_RESULTS.csv"))
s6c_write_csv(joined, file.path(output_root,
                                "FROZEN_BASIN_TRUTH_JOIN.csv"))
s6c_write_csv(basin_truth, file.path(output_root,
                                    "BASIN_SCIENTIFIC_SUMMARY.csv"))
s6c_write_csv(selector_results, file.path(output_root,
                                         "SELECTOR_SCIENTIFIC_RESULTS.csv"))
s6c_write_csv(selector_summary, file.path(output_root,
                                         "SELECTOR_COMPARISON_SUMMARY.csv"))
s6c_write_csv(primary_wide, file.path(output_root,
                                     "PRIMARY_0_90_SELECTOR_COMPARISON.csv"))
s6c_write_csv(integrity, file.path(output_root,
                                  "TRUTH_EVALUATION_INTEGRITY.csv"))
s6c_write_lines(c(
  "status=STAGE6C_FROZEN_BASIN_TRUTH_EVALUATION_COMPLETE",
  "scientific_endpoints=56",
  "truth_free_hashes_verified=TRUE",
  "fit_or_continuation_started=FALSE"
), file.path(output_root, "TRUTH_EVALUATION_COMPLETE.txt"))

cat("STAGE6C_TRUTH_EVALUATION_PASS endpoints=", nrow(scientific),
    " selector_rows=", nrow(selector_results), "\n", sep = "")
