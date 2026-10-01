#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "validate_stage6g_protocol_20260911_v1.R"
}
if (!file.exists(script_path)) {
  stop("Cannot locate protocol validator: ", script_path)
}
root <- dirname(script_path)

read_registered <- function(name) {
  utils::read.csv(file.path(root, name), stringsAsFactors = FALSE,
                  check.names = FALSE)
}

scenarios <- read_registered("SCENARIOS_STAGE6G_V1.csv")
configs <- read_registered("FIT_CONFIGS_STAGE6G_V1.csv")
splits <- read_registered("HOLDOUT_SPLITS_STAGE6G_V1.csv")
starts <- read_registered("START_SEEDS_STAGE6G_V1.csv")
budget <- read_registered("BUDGET_STAGE6G_V1.csv")

checks <- c(
  three_scenarios = nrow(scenarios) == 3L,
  five_common_configs = nrow(configs) == 5L &&
    all(configs$included_in_all_scenarios),
  two_holdout_splits = nrow(splits) == 2L,
  four_paired_starts = nrow(starts) == 4L &&
    all(starts$paired_across_scenarios_splits_and_dimensions),
  unique_ids = !anyDuplicated(scenarios$scenario_id) &&
    !anyDuplicated(configs$fit_config_id) &&
    !anyDuplicated(splits$split_id) && !anyDuplicated(starts$fit_seed),
  new_seed_sets_internal = !anyDuplicated(scenarios$data_seed) &&
    !anyDuplicated(scenarios$truth_basis_seed) &&
    length(intersect(scenarios$data_seed, scenarios$truth_basis_seed)) == 0L,
  expected_truth_in_grid = all(
    scenarios$expected_truth_fit_config_id %in% configs$fit_config_id),
  jaoua_truth_domain = all(
    scenarios$truth_L_f + scenarios$truth_L_s_1 +
      scenarios$truth_L_s_2 <= 4L) &&
    all(scenarios$truth_M_common >= 2L & scenarios$truth_M_common <= 4L) &&
    all(scenarios$truth_basis_method == "jaoua_splinets"),
  common_fixed_design = length(unique(scenarios$S)) == 1L &&
    length(unique(scenarios$n_1)) == 1L &&
    length(unique(scenarios$n_2)) == 1L &&
    length(unique(scenarios$p)) == 1L &&
    length(unique(scenarios$d)) == 1L &&
    length(unique(scenarios$K)) == 1L,
  no_execution_authorization = !any(scenarios$execution_authorized) &&
    !any(scenarios$formal_paper_mc),
  minimum_five_training_times = all(splits$minimum_training_times >= 5L),
  budget_matches_configs = setequal(budget$fit_config_id,
                                    configs$fit_config_id) &&
    all(budget$registered_endpoints == 24L)
)

p <- unique(scenarios$p)
n1 <- unique(scenarios$n_1)
n2 <- unique(scenarios$n_2)
K <- unique(scenarios$K)
Q <- configs$fit_L_f + configs$fit_L_s_1 + configs$fit_L_s_2
expected_complexity <- p * Q + configs$fit_M_common * (
  (n1 + n2) * configs$fit_L_f + n1 * configs$fit_L_s_1 +
    n2 * configs$fit_L_s_2 + (K + 2) * Q)
checks <- c(checks,
  complexity_proxy_exact = identical(as.numeric(configs$complexity_proxy),
                                     as.numeric(expected_complexity)),
  registered_endpoints_120 = nrow(scenarios) * nrow(configs) *
    nrow(splits) * nrow(starts) == 120L,
  configuration_winners_30 = nrow(scenarios) * nrow(configs) *
    nrow(splits) == 30L,
  budget_total_endpoints_120 = sum(budget$registered_endpoints) == 120L
)

r_files <- list.files(root, pattern = "[.]R$", full.names = TRUE)
parse_ok <- vapply(r_files, function(path) {
  tryCatch({
    parse(file = path)
    TRUE
  }, error = function(e) FALSE)
}, logical(1L))
checks <- c(checks, all_protocol_r_files_parse = all(parse_ok))

result <- data.frame(check_id = names(checks), passed = unname(checks),
                     stringsAsFactors = FALSE)
print(result, row.names = FALSE)
cat(sprintf("registered_endpoints=%d\n",
            nrow(scenarios) * nrow(configs) * nrow(splits) * nrow(starts)))
cat(sprintf("configuration_winners=%d\n",
            nrow(scenarios) * nrow(configs) * nrow(splits)))
cat(sprintf("planning_cpu_hours=%.2f\n", sum(budget$planning_cpu_hours)))

if (!all(checks)) {
  quit(save = "no", status = 1L)
}
cat("STAGE6G_PROTOCOL_VALIDATION_PASS\n")
