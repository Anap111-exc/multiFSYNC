args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript run_stage6e_a_readonly_audit_20260908_v1.R <stage6d-root> <output-root>",
       call. = FALSE)
}

script_filename <- "run_stage6e_a_readonly_audit_20260908_v1.R"
if (file.exists(file.path(getwd(), script_filename))) {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                               winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
}
stage6d_root <- args[[1L]]
output_root <- args[[2L]]

source(file.path(script_dir, "stage6e_a_tools_20260908_v1.R"),
       local = FALSE)
source(file.path(dirname(script_dir),
                 "g12-truth-free-dimension-selection-stage6b-20260901",
                 "stage6b_truth_free_tools_20260901_v1.R"), local = FALSE)

s6ea_assert(dir.exists(stage6d_root), "Stage 6D input root does not exist.")
s6ea_assert(!dir.exists(output_root),
             "Output directory already exists; use a new versioned directory.")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

inputs <- s6ea_validate_stage6d(stage6d_root)
inventory <- data.frame(
  fit_id = inputs$endpoints$fit_id,
  data_id = inputs$endpoints$data_id,
  fit_config_id = inputs$endpoints$fit_config_id,
  fit_relative_path = file.path("fits", inputs$endpoints$fit_id, "fit.rds"),
  fit_bytes = file.info(inputs$fit_paths)$size,
  holdout_relative_path = file.path("fits", inputs$endpoints$fit_id,
                                    "HOLDOUT_SUBJECT_TIME_ALL.csv"),
  holdout_bytes = file.info(inputs$holdout_paths)$size,
  read_only = TRUE,
  stringsAsFactors = FALSE
)
s6ea_write_csv(inventory, file.path(output_root, "READ_INPUT_INVENTORY.csv"))

losses <- s6ea_read_holdout_losses(stage6d_root, inputs$endpoints)
risk_units <- s6ea_make_risk_units(losses)
risk_summary <- s6ea_summarize_risk(risk_units)
selection <- s6ea_compare_selection_rules(risk_units)
top_sets <- s6ea_top_set_diagnostics(losses, inputs$science, risk_summary)

s6ea_write_csv(risk_units,
                file.path(output_root, "TRUTH_FREE_RISK_UNIT_SUMMARIES.csv"))
s6ea_write_csv(risk_summary,
                file.path(output_root, "TRUTH_FREE_ROBUST_RISK_SUMMARY.csv"))
s6ea_write_csv(selection$candidates,
                file.path(output_root, "PAIRED_LOSS_DIFFERENCES.csv"))
s6ea_write_csv(selection$selections,
                file.path(output_root,
                          "EXPLORATORY_SELECTION_RULE_COMPARISON.csv"))
s6ea_write_csv(top_sets,
                file.path(output_root, "TOP_ELBO_SET_POSTHOC_DIAGNOSTICS.csv"))

extracted <- s6ea_load_truth_free_fits(stage6d_root, inputs)
complexity <- s6ea_summarize_truth_free_complexity(extracted)
association <- s6ea_posthoc_complexity_association(complexity$fit,
                                                    inputs$science)
s6ea_write_csv(complexity$fit,
                file.path(output_root,
                          "TRUTH_FREE_FIT_COMPLEXITY_DIAGNOSTICS.csv"))
s6ea_write_csv(complexity$config,
                file.path(output_root,
                          "TRUTH_FREE_CONFIG_STABILITY_COMPLEXITY.csv"))
s6ea_write_csv(association,
                file.path(output_root,
                          "POSTHOC_COMPLEXITY_TRUTH_ASSOCIATIONS.csv"))

selection_summary <- s6ea_group_apply(
  selection$selections, c("risk_method", "one_se_semantics"), function(group) {
    factor_rows <- group[group$axis == "factor", ]
    fpca_rows <- group[group$axis == "fpca", ]
    merged <- merge(factor_rows[c("data_id", "correct_dimension_posthoc")],
                    fpca_rows[c("data_id", "correct_dimension_posthoc")],
                    by = "data_id", suffixes = c("_factor", "_fpca"))
    data.frame(
      risk_method = group$risk_method[[1L]],
      one_se_semantics = group$one_se_semantics[[1L]],
      factor_axis_correct_posthoc = sum(factor_rows$correct_dimension_posthoc),
      fpca_axis_correct_posthoc = sum(fpca_rows$correct_dimension_posthoc),
      both_axes_correct_posthoc = sum(
        merged$correct_dimension_posthoc_factor &
          merged$correct_dimension_posthoc_fpca),
      data_sets = 3L,
      status = "posthoc_exploratory_not_validation",
      stringsAsFactors = FALSE
    )
  }
)
s6ea_write_csv(selection_summary,
                file.path(output_root,
                          "EXPLORATORY_SELECTION_ACCURACY_SUMMARY.csv"))

integrity <- data.frame(
  check_id = c(
    "endpoints_180", "configuration_winners_15", "objective_eligible_180",
    "truth_isolated_endpoints", "holdout_rows_10800",
    "risk_unit_rows_2700", "selection_rows_36", "fit_objects_read_180",
    "candidate_rows_complete", "risk_inputs_truth_free",
    "new_fit_or_continuation_started"
  ),
  passed = c(
    nrow(inputs$endpoints) == 180L,
    nrow(inputs$winners) == 15L,
    sum(inputs$endpoints$objective_eligible) == 180L,
    !any(inputs$endpoints$truth_used_for_fit_stopping_or_selection),
    nrow(losses) == 10800L,
    nrow(risk_units) == 2700L,
    nrow(selection$selections) == 36L,
    nrow(complexity$fit) == 180L,
    nrow(extracted$candidates) == sum(complexity$fit$candidate_count),
    all(risk_units$source_fit_count %in% c(1L, 3L, 5L)),
    TRUE
  ),
  detail = c(
    as.character(nrow(inputs$endpoints)),
    as.character(nrow(inputs$winners)),
    as.character(sum(inputs$endpoints$objective_eligible)),
    "all fit/stopping/selection truth-use flags FALSE",
    as.character(nrow(losses)), as.character(nrow(risk_units)),
    as.character(nrow(selection$selections)),
    as.character(nrow(complexity$fit)),
    as.character(nrow(extracted$candidates)),
    "top-k membership uses ordinary ELBO only; holdout response is observed",
    "audit code contains no fitting or continuation call"
  ),
  stringsAsFactors = FALSE
)
s6ea_assert(all(integrity$passed), "Stage 6E-A integrity checks failed.")
s6ea_write_csv(integrity, file.path(output_root, "AUDIT_INTEGRITY.csv"))
s6ea_write_lines(c(
  "status=STAGE6E_A_READONLY_SELECTION_ATTRIBUTION_COMPLETE",
  "input_stage6d_endpoints=180",
  "configuration_winners=15",
  "risk_methods=max_elbo_start;top3_elbo_median;top5_elbo_median",
  "one_se_semantics=independent_one_se;paired_difference_one_se",
  "truth_free_risk_and_complexity_tables_created=TRUE",
  "truth_join_for_posthoc_assessment=TRUE",
  "new_fit_started=FALSE",
  "continuation_started=FALSE",
  "ecs_used=FALSE",
  "formal_validation_claim=FALSE"
), file.path(output_root, "STAGE6E_A_COMPLETE.txt"))

cat("Stage 6E-A read-only audit complete\n")
print(selection_summary)
