#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_name <- basename(sub("^--file=", "", script_argument[[1L]]))
audit_root <- "."
if (!file.exists(script_name)) {
  stop("Run from the directory containing the audit script.", call. = FALSE)
}

if (length(commandArgs(trailingOnly = TRUE))) {
  stop("No command-line arguments are allowed.", call. = FALSE)
}

output_dir <- file.path(audit_root, "convergence_semantics_audit_20260824_v1")
if (dir.exists(output_dir)) {
  stop("Versioned audit output already exists.", call. = FALSE)
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing source: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
as_flag <- function(value) {
  if (is.logical(value)) return(value)
  toupper(as.character(value)) == "TRUE"
}
as_number <- function(value) suppressWarnings(as.numeric(value))
safe_mean <- function(value) {
  value <- as_number(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
safe_median <- function(value) {
  value <- as_number(value)
  value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}
safe_quantile <- function(value, probability) {
  value <- as_number(value)
  value <- value[is.finite(value)]
  if (length(value)) {
    unname(stats::quantile(value, probability, names = FALSE, type = 7))
  } else NA_real_
}
bind_rows <- function(values) {
  values <- values[!vapply(values, is.null, logical(1L))]
  if (length(values)) do.call(rbind, values) else data.frame()
}

source_dir <- file.path(audit_root, "sources")
endpoint_paths <- c(
  independent = file.path(source_dir, "G12IC_ENDPOINTS_20260823.csv"),
  cross_scenario = file.path(source_dir, "G12CS_ENDPOINTS_20260824.csv")
)
science_paths <- c(
  independent = file.path(source_dir, "G12IC_SCIENCE_20260823.csv"),
  cross_scenario = file.path(source_dir, "G12CS_SCIENCE_20260824.csv")
)

independent_endpoint <- read_csv(endpoint_paths[["independent"]])
cross_endpoint <- read_csv(endpoint_paths[["cross_scenario"]])
independent_science <- read_csv(science_paths[["independent"]])
cross_science <- read_csv(science_paths[["cross_scenario"]])

if (nrow(independent_endpoint) != 144L || nrow(cross_endpoint) != 96L ||
    nrow(independent_science) != 144L || nrow(cross_science) != 96L) {
  stop("Unexpected source row count.", call. = FALSE)
}

independent_endpoint$experiment_id <- "G12IC_20260822"
independent_endpoint$scenario_id <- "independent_baseline_strong"
cross_endpoint$experiment_id <- "G12CS_20260823"
independent_science$experiment_id <- "G12IC_20260822"
independent_science$scenario_id <- "independent_baseline_strong"
cross_science$experiment_id <- "G12CS_20260823"

required_endpoint <- c(
  "fit_id", "data_id", "scenario_id", "method_id", "seed_index",
  "t1_sweeps", "exact_correct_1_1_1", "strict_practical_converged",
  "original_longest_streak", "quantile_factor_alternative_ever_5",
  "output_only_quantile_factor_ever_5", "endpoint_objective_pass",
  "endpoint_fitted_gate_pass", "endpoint_rss_gate_pass",
  "endpoint_ppi_max_gate_pass", "endpoint_ppi_quantile_factor_gate_pass",
  "minimum_gate_groups_removed_for_5_streak",
  "remove_ppi_max_alone_reaches_5", "endpoint_ppi_max_abs",
  "endpoint_ppi_quantile_abs", "endpoint_factor_ppi_max_abs",
  "endpoint_long_ppi_max_abs", "endpoint_long_ppi_quantile_abs",
  "endpoint_long_factor_ppi_max_abs", "experiment_id"
)
if (!all(required_endpoint %in% names(independent_endpoint)) ||
    !all(required_endpoint %in% names(cross_endpoint))) {
  stop("Endpoint source schema mismatch.", call. = FALSE)
}

endpoints <- rbind(
  independent_endpoint[, required_endpoint],
  cross_endpoint[, required_endpoint]
)
rownames(endpoints) <- NULL

logical_fields <- c(
  "exact_correct_1_1_1", "strict_practical_converged",
  "quantile_factor_alternative_ever_5",
  "output_only_quantile_factor_ever_5", "endpoint_objective_pass",
  "endpoint_fitted_gate_pass", "endpoint_rss_gate_pass",
  "endpoint_ppi_max_gate_pass", "endpoint_ppi_quantile_factor_gate_pass",
  "remove_ppi_max_alone_reaches_5"
)
for (field in logical_fields) endpoints[[field]] <- as_flag(endpoints[[field]])

if (nrow(endpoints) != 240L || anyDuplicated(endpoints$fit_id) ||
    !identical(sort(unique(endpoints$method_id)),
               sort(c("current_1_over_m", "gram_unit_energy"))) ||
    !all(table(endpoints$method_id) == 120L) ||
    length(unique(endpoints$data_id)) != 12L) {
  stop("Combined endpoint identity check failed.", call. = FALSE)
}

required_science <- c(
  "fit_id", "data_id", "scenario_id", "method_id", "objective_eligible",
  "truth_used_for_fit_or_selection", "formal_v0lv_result",
  "continuation_used", "selected_counts", "exact_correct_1_1_1",
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_components_total", "feature_components_matched",
  "feature_total_ise_mean", "feature_estimate_to_projection_ise_mean",
  "loading_relative_l2_error_mean", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean",
  "experiment_id"
)
if (!all(required_science %in% names(independent_science)) ||
    !all(required_science %in% names(cross_science))) {
  stop("Scientific source schema mismatch.", call. = FALSE)
}
science <- rbind(
  independent_science[, required_science],
  cross_science[, required_science]
)
rownames(science) <- NULL
for (field in c("objective_eligible", "truth_used_for_fit_or_selection",
                "formal_v0lv_result", "continuation_used",
                "exact_correct_1_1_1")) {
  science[[field]] <- as_flag(science[[field]])
}

endpoint_index <- match(science$fit_id, endpoints$fit_id)
if (anyNA(endpoint_index) || !identical(science$method_id,
                                        endpoints$method_id[endpoint_index]) ||
    !identical(science$exact_correct_1_1_1,
               endpoints$exact_correct_1_1_1[endpoint_index])) {
  stop("Endpoint/science fit join failed.", call. = FALSE)
}
science$strict_practical_converged <-
  endpoints$strict_practical_converged[endpoint_index]
science$quantile_factor_alternative_ever_5 <-
  endpoints$quantile_factor_alternative_ever_5[endpoint_index]
science$remove_ppi_max_alone_reaches_5 <-
  endpoints$remove_ppi_max_alone_reaches_5[endpoint_index]
science$feature_component_coverage <-
  as_number(science$feature_components_matched) /
  as_number(science$feature_components_total)

summarize_endpoint_group <- function(rows, method_id, scope) {
  incorrect <- !rows$exact_correct_1_1_1
  data.frame(
    scope = scope, method_id = method_id, endpoints = nrow(rows),
    data_sets = length(unique(rows$data_id)),
    exact_correct_n = sum(rows$exact_correct_1_1_1),
    exact_correct_rate = mean(rows$exact_correct_1_1_1),
    strict_n = sum(rows$strict_practical_converged),
    strict_rate = mean(rows$strict_practical_converged),
    correct_and_strict_n = sum(
      rows$exact_correct_1_1_1 & rows$strict_practical_converged
    ),
    incorrect_and_strict_n = sum(
      incorrect & rows$strict_practical_converged
    ),
    endpoint_objective_pass_n = sum(rows$endpoint_objective_pass),
    endpoint_fitted_pass_n = sum(rows$endpoint_fitted_gate_pass),
    endpoint_rss_pass_n = sum(rows$endpoint_rss_gate_pass),
    endpoint_ppi_max_pass_n = sum(rows$endpoint_ppi_max_gate_pass),
    endpoint_quantile_factor_pass_n = sum(
      rows$endpoint_ppi_quantile_factor_gate_pass
    ),
    quantile_factor_alternative_ever_5_n = sum(
      rows$quantile_factor_alternative_ever_5
    ),
    output_only_quantile_factor_ever_5_n = sum(
      rows$output_only_quantile_factor_ever_5
    ),
    remove_ppi_max_alone_reaches_5_n = sum(
      rows$remove_ppi_max_alone_reaches_5
    ),
    minimum_removed_gate_groups_median = safe_median(
      rows$minimum_gate_groups_removed_for_5_streak
    ),
    endpoints_at_200 = sum(as_number(rows$t1_sweeps) == 200L),
    stringsAsFactors = FALSE
  )
}

method_summary <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy"), function(method_id) {
    summarize_endpoint_group(
      endpoints[endpoints$method_id == method_id, , drop = FALSE],
      method_id, "combined_two_experiments"
    )
  }
))

experiment_method_summary <- bind_rows(lapply(
  unique(endpoints$experiment_id), function(experiment_id) {
    bind_rows(lapply(c("current_1_over_m", "gram_unit_energy"),
                     function(method_id) {
      rows <- endpoints[endpoints$experiment_id == experiment_id &
                          endpoints$method_id == method_id, , drop = FALSE]
      summarize_endpoint_group(rows, method_id, experiment_id)
    }))
  }
))

scenario_method_summary <- bind_rows(lapply(
  unique(endpoints$scenario_id), function(scenario_id) {
    bind_rows(lapply(c("current_1_over_m", "gram_unit_energy"),
                     function(method_id) {
      rows <- endpoints[endpoints$scenario_id == scenario_id &
                          endpoints$method_id == method_id, , drop = FALSE]
      if (!nrow(rows)) return(NULL)
      summarize_endpoint_group(rows, method_id, scenario_id)
    }))
  }
))

data_method_summary <- bind_rows(lapply(
  split(endpoints, interaction(endpoints$data_id, endpoints$method_id,
                               drop = TRUE, lex.order = TRUE)),
  function(rows) data.frame(
    experiment_id = rows$experiment_id[[1L]],
    scenario_id = rows$scenario_id[[1L]], data_id = rows$data_id[[1L]],
    method_id = rows$method_id[[1L]], endpoints = nrow(rows),
    exact_correct_n = sum(rows$exact_correct_1_1_1),
    strict_n = sum(rows$strict_practical_converged),
    quantile_factor_alternative_ever_5_n = sum(
      rows$quantile_factor_alternative_ever_5
    ),
    remove_ppi_max_alone_reaches_5_n = sum(
      rows$remove_ppi_max_alone_reaches_5
    ), stringsAsFactors = FALSE
  )
))

gate_fields <- c(
  strict_practical_converged = "strict_practical",
  endpoint_objective_pass = "endpoint_objective",
  endpoint_fitted_gate_pass = "endpoint_fitted",
  endpoint_rss_gate_pass = "endpoint_rss",
  endpoint_ppi_max_gate_pass = "endpoint_ppi_max",
  endpoint_ppi_quantile_factor_gate_pass = "endpoint_quantile_factor",
  remove_ppi_max_alone_reaches_5 = "remove_ppi_max_ever_5",
  quantile_factor_alternative_ever_5 = "quantile_factor_ever_5",
  output_only_quantile_factor_ever_5 = "output_only_quantile_factor_ever_5"
)

gate_by_correctness <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy"), function(method_id) {
    bind_rows(lapply(c(FALSE, TRUE), function(correct) {
      rows <- endpoints[endpoints$method_id == method_id &
                          endpoints$exact_correct_1_1_1 == correct,
                        , drop = FALSE]
      bind_rows(lapply(names(gate_fields), function(field) data.frame(
        method_id = method_id, exact_correct_1_1_1 = correct,
        endpoints = nrow(rows), gate = unname(gate_fields[[field]]),
        pass_n = sum(rows[[field]]), pass_rate = mean(rows[[field]]),
        stringsAsFactors = FALSE
      )))
    }))
  }
))

rule_fields <- c(
  strict_practical = "strict_practical_converged",
  remove_ppi_max_ever_5 = "remove_ppi_max_alone_reaches_5",
  quantile_factor_ever_5 = "quantile_factor_alternative_ever_5",
  output_only_quantile_factor_ever_5 =
    "output_only_quantile_factor_ever_5"
)
rule_semantic_association <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy", "ALL"), function(method_id) {
    rows <- if (method_id == "ALL") endpoints else
      endpoints[endpoints$method_id == method_id, , drop = FALSE]
    bind_rows(lapply(names(rule_fields), function(rule_name) {
      pass <- rows[[rule_fields[[rule_name]]]]
      correct <- rows$exact_correct_1_1_1
      data.frame(
        method_id = method_id, rule = rule_name, endpoints = nrow(rows),
        correct_total = sum(correct), pass_correct = sum(pass & correct),
        pass_rate_among_correct = mean(pass[correct]),
        incorrect_total = sum(!correct), pass_incorrect = sum(pass & !correct),
        pass_rate_among_incorrect = mean(pass[!correct]),
        passed_total = sum(pass),
        correct_fraction_among_passed = if (sum(pass)) {
          mean(correct[pass])
        } else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }
))

minimum_removal_distribution <- as.data.frame(table(
  method_id = endpoints$method_id,
  exact_correct_1_1_1 = endpoints$exact_correct_1_1_1,
  minimum_removed_gate_groups =
    endpoints$minimum_gate_groups_removed_for_5_streak
), stringsAsFactors = FALSE)
minimum_removal_distribution <- minimum_removal_distribution[
  minimum_removal_distribution$Freq > 0L, , drop = FALSE
]

ppi_fields <- c(
  "endpoint_ppi_max_abs", "endpoint_ppi_quantile_abs",
  "endpoint_factor_ppi_max_abs", "endpoint_long_ppi_max_abs",
  "endpoint_long_ppi_quantile_abs", "endpoint_long_factor_ppi_max_abs"
)
ppi_endpoint_distribution <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy"), function(method_id) {
    bind_rows(lapply(c(FALSE, TRUE), function(correct) {
      rows <- endpoints[endpoints$method_id == method_id &
                          endpoints$exact_correct_1_1_1 == correct,
                        , drop = FALSE]
      bind_rows(lapply(ppi_fields, function(metric) data.frame(
        method_id = method_id, exact_correct_1_1_1 = correct,
        metric = metric, n = length(rows[[metric]]),
        minimum = safe_quantile(rows[[metric]], 0),
        q25 = safe_quantile(rows[[metric]], 0.25),
        median = safe_quantile(rows[[metric]], 0.5),
        q75 = safe_quantile(rows[[metric]], 0.75),
        maximum = safe_quantile(rows[[metric]], 1),
        stringsAsFactors = FALSE
      )))
    }))
  }
))

scientific_metrics <- c(
  "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
  "feature_total_ise_mean", "feature_estimate_to_projection_ise_mean",
  "loading_relative_l2_error_mean", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "kernel_relative_ise_mean", "covariance_operator_relative_error_mean",
  "feature_component_coverage"
)
summarize_science <- function(rows, grouping) {
  bind_rows(lapply(scientific_metrics, function(metric) data.frame(
    grouping, metric = metric,
    n_finite = sum(is.finite(as_number(rows[[metric]]))),
    mean = safe_mean(rows[[metric]]), median = safe_median(rows[[metric]]),
    stringsAsFactors = FALSE, check.names = FALSE
  )))
}

science_by_method_strict <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy"), function(method_id) {
    bind_rows(lapply(c(FALSE, TRUE), function(strict) {
      rows <- science[science$method_id == method_id &
                        science$strict_practical_converged == strict,
                      , drop = FALSE]
      if (!nrow(rows)) return(NULL)
      summarize_science(rows, data.frame(
        method_id = method_id, strict_practical_converged = strict,
        endpoints = nrow(rows), stringsAsFactors = FALSE
      ))
    }))
  }
))

science_by_method_correctness <- bind_rows(lapply(
  c("current_1_over_m", "gram_unit_energy"), function(method_id) {
    bind_rows(lapply(c(FALSE, TRUE), function(correct) {
      rows <- science[science$method_id == method_id &
                        science$exact_correct_1_1_1 == correct,
                      , drop = FALSE]
      summarize_science(rows, data.frame(
        method_id = method_id, exact_correct_1_1_1 = correct,
        endpoints = nrow(rows), stringsAsFactors = FALSE
      ))
    }))
  }
))

selected_counts_by_strict <- as.data.frame(table(
  method_id = science$method_id,
  strict_practical_converged = science$strict_practical_converged,
  selected_counts = science$selected_counts,
  exact_correct_1_1_1 = science$exact_correct_1_1_1
), stringsAsFactors = FALSE)
selected_counts_by_strict <- selected_counts_by_strict[
  selected_counts_by_strict$Freq > 0L, , drop = FALSE
]

qc <- data.frame(
  check_id = c(
    "endpoint_rows_240", "science_rows_240", "unique_fit_ids",
    "methods_120_each", "twelve_data_sets", "all_objective_eligible",
    "truth_free_fit_and_selection", "no_continuation", "not_formal_v0lv",
    "endpoint_science_exact_join", "source_experiments_144_96"
  ),
  passed = c(
    nrow(endpoints) == 240L, nrow(science) == 240L,
    !anyDuplicated(endpoints$fit_id) && !anyDuplicated(science$fit_id),
    all(table(endpoints$method_id) == 120L),
    length(unique(endpoints$data_id)) == 12L,
    all(science$objective_eligible),
    !any(science$truth_used_for_fit_or_selection),
    !any(science$continuation_used), !any(science$formal_v0lv_result),
    identical(science$exact_correct_1_1_1,
              endpoints$exact_correct_1_1_1[endpoint_index]),
    identical(as.integer(table(endpoints$experiment_id)[
      c("G12IC_20260822", "G12CS_20260823")]), c(144L, 96L))
  ),
  detail = c(
    "240/240", "240/240", "240/240", "current=120;G12=120", "12",
    "240/240", "0/240", "0/240", "0/240", "240/240", "144;96"
  ), stringsAsFactors = FALSE
)
if (!all(qc$passed)) {
  stop("Audit QC failed: ", paste(qc$check_id[!qc$passed], collapse = ";"),
       call. = FALSE)
}

dir.create(output_dir, recursive = FALSE)
write_table <- function(value, name) utils::write.csv(
  value, file.path(output_dir, name), row.names = FALSE,
  quote = TRUE, na = ""
)
write_table(qc, "AUDIT_QC.csv")
write_table(endpoints, "COMBINED_240_ENDPOINT_ATTRIBUTION.csv")
write_table(method_summary, "METHOD_CONVERGENCE_SUMMARY.csv")
write_table(experiment_method_summary, "EXPERIMENT_METHOD_SUMMARY.csv")
write_table(scenario_method_summary, "SCENARIO_METHOD_SUMMARY.csv")
write_table(data_method_summary, "DATA_METHOD_SUMMARY.csv")
write_table(gate_by_correctness, "GATE_PASS_BY_METHOD_CORRECTNESS.csv")
write_table(rule_semantic_association, "RULE_SEMANTIC_ASSOCIATION.csv")
write_table(minimum_removal_distribution,
            "MINIMUM_GATE_REMOVAL_DISTRIBUTION.csv")
write_table(ppi_endpoint_distribution, "PPI_ENDPOINT_DISTRIBUTION.csv")
write_table(science_by_method_strict,
            "SCIENTIFIC_METRICS_BY_METHOD_STRICT.csv")
write_table(science_by_method_correctness,
            "SCIENTIFIC_METRICS_BY_METHOD_CORRECTNESS.csv")
write_table(selected_counts_by_strict, "SELECTED_COUNTS_BY_STRICT.csv")
capture.output(sessionInfo(), file = file.path(output_dir, "SESSION_INFO.txt"))
writeLines(c(
  "status=PASS", "audit_id=G12_CONVERGENCE_SEMANTICS_OFFLINE_AUDIT_V1_20260824",
  paste0("completed_utc=", format(Sys.time(), tz = "UTC",
                                  format = "%Y-%m-%dT%H:%M:%SZ")),
  "source_endpoints=240", "source_data_sets=12",
  "new_data_generated=0", "new_fits_started=0",
  "continuation_started=FALSE", "package_source_modified=FALSE",
  "historical_results_modified=FALSE"
), file.path(output_dir, "AUDIT_COMPLETE.txt"), useBytes = TRUE)

cat("AUDIT_COMPLETE endpoints=240 data_sets=12 new_fits=0 continuation=0\n")
print(method_summary, row.names = FALSE)
print(rule_semantic_association, row.names = FALSE)
