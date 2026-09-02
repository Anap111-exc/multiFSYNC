#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
parse_argument <- function(name) {
  prefix <- paste0("--", name, "=")
  value <- arguments[startsWith(arguments, prefix)]
  if (!length(value)) return("")
  sub(prefix, "", value[[1L]], fixed = TRUE)
}
input_root <- parse_argument("input-root")
output_root <- parse_argument("output-root")
if (!nzchar(input_root) || !nzchar(output_root)) {
  stop("--input-root and --output-root are required.", call. = FALSE)
}
if (!file.exists(file.path(input_root, "AUDIT_COMPLETE.txt"))) {
  stop("Input is not a completed Stage 6B-A audit.", call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_table <- function(filename) {
  utils::read.csv(file.path(input_root, filename), stringsAsFactors = FALSE,
                  check.names = FALSE)
}
write_table <- function(value, filename) {
  utils::write.csv(value, file.path(output_root, filename), row.names = FALSE,
                   na = "")
}
bind_groups <- function(value, keys, summarize) {
  groups <- split(value, interaction(value[keys], drop = TRUE))
  result <- do.call(rbind, lapply(groups, summarize))
  rownames(result) <- NULL
  result
}

candidate <- read_table("FIT_CANDIDATE_ACTIVITY.csv")
matches <- read_table("REFERENCE_MATCH_STABILITY.csv")
duplicates <- read_table("WITHIN_ROLE_DUPLICATION.csv")
leakage <- read_table("CROSS_ROLE_LEAKAGE.csv")
loading_subspace <- read_table("LOADING_SUBSPACE_STABILITY.csv")
fpca_subspace <- read_table("FPCA_SUBSPACE_STABILITY.csv")

ppi_energy <- bind_groups(candidate, "fit_config_id", function(group) {
  data.frame(
    fit_config_id = group$fit_config_id[[1L]],
    candidate_rows = nrow(group),
    zero_loading_directions = sum(group$zero_loading_direction),
    zero_loading_rate = mean(group$zero_loading_direction),
    ppi_retained_rows = sum(group$ppi_retained_0_5),
    ppi_retained_share_le_0_01 = sum(
      group$ppi_retained_0_5 & group$weighted_contribution_share <= 0.01
    ),
    ppi_retained_zero_loading = sum(
      group$ppi_retained_0_5 & group$zero_loading_direction
    ),
    median_contribution_share = stats::median(
      group$weighted_contribution_share
    ),
    maximum_contribution_share = max(group$weighted_contribution_share),
    stringsAsFactors = FALSE
  )
})

effective_m_table <- as.data.frame(with(candidate, table(
  fit_config_id, effective_M_99pct
)), stringsAsFactors = FALSE)
names(effective_m_table) <- c("fit_config_id", "effective_M_99pct", "count")
effective_m_table <- effective_m_table[effective_m_table$count > 0L,
                                       , drop = FALSE]

direction_extended <- bind_groups(
  matches, c("fit_config_id", "role_id", "reference_candidate_id"),
  function(group) {
    data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      role_id = group$role_id[[1L]],
      reference_candidate_id = group$reference_candidate_id[[1L]],
      starts = nrow(group),
      loading_abs_cosine_median = stats::median(group$loading_abs_cosine),
      loading_abs_cosine_minimum = min(group$loading_abs_cosine),
      recurrence_ge_0_9 = mean(group$loading_abs_cosine >= 0.9),
      contribution_share_q25 = as.numeric(stats::quantile(
        group$candidate_contribution_share, 0.25
      )),
      contribution_share_median = stats::median(
        group$candidate_contribution_share
      ),
      contribution_share_q75 = as.numeric(stats::quantile(
        group$candidate_contribution_share, 0.75
      )),
      contribution_share_maximum = max(group$candidate_contribution_share),
      complete_contribution_rms_median = stats::median(
        group$candidate_complete_contribution_rms
      ),
      factor_ppi_median = stats::median(group$candidate_factor_ppi),
      stringsAsFactors = FALSE
    )
  }
)

duplication_summary <- bind_groups(
  duplicates, c("fit_config_id", "role_id"), function(group) {
    data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      role_id = group$role_id[[1L]],
      fit_role_units = nrow(group),
      duplicate_rate_ge_0_9 = mean(group$pairs_ge_0_9 > 0L),
      median_maximum_within_role_abs_cosine = stats::median(
        group$maximum_within_role_abs_cosine
      ),
      maximum_within_role_abs_cosine = max(
        group$maximum_within_role_abs_cosine
      ),
      stringsAsFactors = FALSE
    )
  }
)

leakage_summary <- bind_groups(
  leakage, c("fit_config_id", "study"), function(group) {
    data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      study = group$study[[1L]],
      fit_study_units = nrow(group),
      cross_role_leakage_rate_ge_0_9 = mean(group$pairs_ge_0_9 > 0L),
      median_maximum_cross_role_abs_cosine = stats::median(
        group$maximum_cross_role_abs_cosine
      ),
      maximum_cross_role_abs_cosine = max(group$maximum_cross_role_abs_cosine),
      stringsAsFactors = FALSE
    )
  }
)

loading_subspace_summary <- bind_groups(
  loading_subspace, c("fit_config_id", "role_id"), function(group) {
    data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      role_id = group$role_id[[1L]],
      starts = nrow(group),
      median_mean_squared_canonical_correlation = stats::median(
        group$mean_squared_canonical_correlation
      ),
      minimum_mean_squared_canonical_correlation = min(
        group$mean_squared_canonical_correlation
      ),
      median_minimum_canonical_correlation = stats::median(
        group$minimum_canonical_correlation
      ),
      minimum_canonical_correlation = min(group$minimum_canonical_correlation),
      stringsAsFactors = FALSE
    )
  }
)

fpca_subspace_summary <- bind_groups(
  fpca_subspace,
  c("fit_config_id", "role_id", "reference_candidate_id"),
  function(group) {
    data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      role_id = group$role_id[[1L]],
      reference_candidate_id = group$reference_candidate_id[[1L]],
      starts = nrow(group),
      median_mean_squared_canonical_correlation = stats::median(
        group$mean_squared_canonical_correlation
      ),
      minimum_mean_squared_canonical_correlation = min(
        group$mean_squared_canonical_correlation
      ),
      median_minimum_canonical_correlation = stats::median(
        group$minimum_canonical_correlation
      ),
      minimum_canonical_correlation = min(group$minimum_canonical_correlation),
      stringsAsFactors = FALSE
    )
  }
)

write_table(ppi_energy, "PPI_ENERGY_SUMMARY.csv")
write_table(effective_m_table, "EFFECTIVE_M_COUNTS.csv")
write_table(direction_extended, "REFERENCE_DIRECTION_EXTENDED_SUMMARY.csv")
write_table(duplication_summary, "DUPLICATION_ROLE_SUMMARY.csv")
write_table(leakage_summary, "CROSS_ROLE_LEAKAGE_STUDY_SUMMARY.csv")
write_table(loading_subspace_summary, "LOADING_SUBSPACE_SUMMARY.csv")
write_table(fpca_subspace_summary, "FPCA_MATCHED_DIRECTION_SUMMARY.csv")
writeLines(c(
  "status=STAGE6B_A_DESCRIPTIVE_SUMMARY_COMPLETE",
  "truth_fields_read=FALSE",
  "fit_or_continuation_started=FALSE"
), file.path(output_root, "SUMMARY_COMPLETE.txt"), useBytes = TRUE)

cat("STAGE6B_A_DESCRIPTIVE_SUMMARY_PASS output_root=", output_root,
    "\n", sep = "")
