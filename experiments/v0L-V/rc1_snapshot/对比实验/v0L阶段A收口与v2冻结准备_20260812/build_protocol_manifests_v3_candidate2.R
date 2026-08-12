# Build the third-round candidate manifest. This script creates metadata only;
# it never generates data, starts a fit, starts continuation, or freezes use.

candidate2_script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run this builder with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/",
                mustWork = TRUE)
}

script_path <- candidate2_script_path()
base_dir <- dirname(script_path)
project_root <- normalizePath(file.path(base_dir, "..", ".."),
                              winslash = "/", mustWork = TRUE)
source(file.path(base_dir, "v0lv_candidate2_runtime.R"), local = FALSE)

MANIFEST_ID <- "V0LV_MANIFEST_SCHEMA_V3_CANDIDATE2_20260812"
TOTAL_ID <- "V0LV_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812"
R_ID <- "R_ONLY_MULTISTART_V3_CANDIDATE2_20260812"
output_dir <- file.path(base_dir, "protocol_manifests_20260812_v3_candidate2")
if (dir.exists(output_dir)) stop("Candidate2 manifest directory already exists.")
dir.create(output_dir, recursive = TRUE)

write_csv <- function(value, name) {
  c2_atomic_write_csv(value, file.path(output_dir, name))
}

data_manifest <- data.frame(
  manifest_schema_id = MANIFEST_ID, total_protocol_id = TOTAL_ID,
  data_id = sprintf("v0lv_%02d", seq_len(10L)), data_index = seq_len(10L),
  data_seed = 82612000L + seq_len(10L), unseen_at_candidate_creation = TRUE,
  formal_generation_started = FALSE, observation_bundle_created = FALSE,
  sealed_truth_bundle_created = FALSE,
  budget_audit_subset = seq_len(10L) %in% c(3L, 8L),
  stringsAsFactors = FALSE
)
write_csv(data_manifest, "V0LV_DATA_SEEDS_CANDIDATE2.csv")

fit_rows <- list(); fit_index <- 0L
add_fit <- function(data_index, path, seed_index, fit_seed) {
  data_id <- sprintf("v0lv_%02d", data_index)
  is_r <- path == "R"; is_pooled <- path == "pooled_bayesSYNC"
  core <- !is_pooled && seed_index <= 4L
  route <- if (is_r) "random__pre1__anneal_default" else path
  fit_index <<- fit_index + 1L
  fit_rows[[fit_index]] <<- data.frame(
    manifest_schema_id = MANIFEST_ID, total_protocol_id = TOTAL_ID,
    r_only_protocol_id = if (is_r) R_ID else "not_applicable",
    fit_id = sprintf("%s_%s_%02d", data_id, path, seed_index),
    data_id = data_id,
    method_family = if (is_pooled) "bayesSYNC_reference" else "multiFSYNC",
    path_id = path, route_id = route, seed_index = as.integer(seed_index),
    fit_seed = as.integer(fit_seed),
    seed_pairing_group = if (core) sprintf("%s_core_%02d", data_id, seed_index)
                         else "none",
    paired_core_seed = core,
    initialization_independence_id = paste(data_id, path, fit_seed, sep = "__"),
    initialization = "random",
    pre_score_sweeps = if (is_r) 1L else 0L,
    anneal = if (path == "B0") "none" else "c(1,1.9,100)",
    planned_annealing_sweeps = if (path == "B0") 0L else 99L,
    main_endpoint_policy = if (is_pooled) {
      "reference_ELBO_stop_or_maxit299"
    } else "strict_practical_stop_or_maximum_200_T1",
    maximum_ordinary_T1_sweeps = 200L,
    terminal_record_required = TRUE,
    selection_pool_size = if (is_r) 12L else if (is_pooled) 1L else 4L,
    selection_rule = if (is_pooled) "single_preregistered_endpoint"
                     else "maximum_eligible_ordinary_T1_ELBO",
    actual_racing = FALSE, offline_racing_replay = is_r,
    input_bundle_class = "observation_only",
    truth_available_to_fit = FALSE, truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE, formal_fit_started = FALSE,
    terminal_record_written = FALSE, stringsAsFactors = FALSE
  )
}

for (data_index in seq_len(10L)) {
  for (path in c("B0", "A")) for (seed_index in seq_len(4L)) {
    add_fit(data_index, path, seed_index,
            82620000L + 100L * data_index + seed_index)
  }
  for (seed_index in seq_len(12L)) {
    add_fit(data_index, "R", seed_index,
            82620000L + 100L * data_index + seed_index)
  }
  add_fit(data_index, "pooled_bayesSYNC", 1L, 82630000L + data_index)
}
fit_manifest <- do.call(rbind, fit_rows)
stopifnot(nrow(fit_manifest) == 210L, !anyDuplicated(fit_manifest$fit_id))
write_csv(fit_manifest, "V0LV_FIT_MANIFEST_CANDIDATE2.csv")

audit_rows <- list(); audit_index <- 0L
for (data_index in c(3L, 8L)) for (seed_index in seq_len(12L)) {
  data_id <- sprintf("v0lv_%02d", data_index)
  fit_seed <- 82620000L + 100L * data_index + seed_index
  initialization_id <- paste(data_id, "R", fit_seed, sep = "__")
  stage <- c("anchor_200", "to_400", "to_800")
  target <- c(200L, 400L, 800L)
  start <- c("main_endpoint_or_verified_main_200",
             "verified_audit_anchor_200", "verified_audit_state_400")
  additional <- c("max(0,200-main_actual_T1)", "200", "400")
  for (part in seq_along(stage)) {
    audit_index <- audit_index + 1L
    audit_rows[[audit_index]] <- data.frame(
      manifest_schema_id = MANIFEST_ID, total_protocol_id = TOTAL_ID,
      r_only_protocol_id = R_ID,
      audit_id = sprintf("%s_R_%02d_%s", data_id, seed_index, stage[[part]]),
      data_id = data_id, path_id = "R",
      route_id = "random__pre1__anneal_default",
      seed_index = seed_index, original_fit_seed = fit_seed,
      continuation_seed = fit_seed,
      initialization_independence_id = initialization_id,
      audit_stage = stage[[part]], start_state_policy = start[[part]],
      additional_T1_sweeps_policy = additional[[part]],
      target_cumulative_T1 = target[[part]],
      fixed_horizon_ignores_intermediate_strict_stop = TRUE,
      sequential_from_verified_state = TRUE,
      main_endpoint_replaced = FALSE, main_selection_reopened = FALSE,
      counterfactual_selection_only = TRUE,
      main_error_policy = "audit_unavailable_due_to_main_error_no_rerun",
      truth_available_to_continuation = FALSE,
      truth_available_to_counterfactual_selection = FALSE,
      formal_continuation_started = FALSE, terminal_record_required = TRUE,
      stringsAsFactors = FALSE
    )
  }
}
audit_manifest <- do.call(rbind, audit_rows)
stopifnot(nrow(audit_manifest) == 72L, !anyDuplicated(audit_manifest$audit_id))
write_csv(audit_manifest, "V0LV_AUDIT_MANIFEST_CANDIDATE2.csv")

config <- data.frame(
  section = c(
    rep("version_binding", 4L), rep("generator", 20L),
    rep("multiFSYNC_model", 14L), rep("R_v3", 12L),
    rep("pooled_bayesSYNC", 9L), rep("evaluation", 8L)
  ),
  key = c(
    "manifest_schema_id", "total_protocol_id", "R_only_protocol_id", "status",
    "function", "S", "n_s", "p", "d", "L_f_true", "L_s_true",
    "M_f_true", "M_s_true", "K", "n_obs", "common_grid", "n_obs_range",
    "sigma_eps", "bool_sparse_loadings", "prop_sparse", "score_var_decay",
    "n_dense", "mean_amp", "beta_amp",
    "L_f_fit", "L_s_fit", "M_f_fit", "M_s_fit", "K", "n_g", "tol_abs",
    "tol_rel", "d_0", "bool_scale", "bool_var_spec_prob", "lambda_orth",
    "n_cpus", "initialization",
    "starts", "pre_score_sweeps", "anneal", "planned_annealing_sweeps",
    "ordinary_T1_maximum", "strict_stop_enabled", "winner_rule",
    "objective_eligibility_separate_from_strict", "actual_racing",
    "offline_racing_screen_T1", "offline_racing_minimum_keep",
    "offline_racing_elbo_margin",
    "combine_studies_once", "study_labels", "Q", "M", "K", "n_g",
    "anneal", "maxit", "cross_model_ELBO_comparison",
    "factor_ppi_rule", "oracle_matching", "bootstrap_cluster", "bootstrap_B",
    "bootstrap_seed", "success_gate_name", "worsening_count_role",
    "oracle_L_M_scope"
  ),
  value = c(
    MANIFEST_ID, TOTAL_ID, R_ID, "candidate_not_frozen_not_authorized",
    "simulate_multi_study_structured", "2", "c(30,30)", "500", "0", "1",
    "c(1,1)", "2", "list(2,2)", "5", "8", "FALSE", "c(6,9)", "0.3",
    "TRUE", "0.9", "TRUE", "201", "0.6", "0",
    "1", "c(1,1)", "2", "list(2,2)", "5", "51", "1e-3", "1e-5",
    "500", "FALSE", "FALSE", "0", "1", "random",
    "12", "1", "c(1,1.9,100)", "99", "200", "TRUE",
    "maximum_eligible_ordinary_T1_ELBO", "TRUE", "FALSE", "100", "4",
    "500",
    "TRUE", "not_passed", "3", "2", "5", "51", "c(1,1.9,100)",
    "299", "forbidden", ">=0.5_auxiliary_only", "one_to_one_A0_B10_B20",
    "data_id", "9999", "82640001", "中位无恶化/优势门槛",
    "descriptive_only_not_independent_gate", "true_L_and_M_auxiliary_only"
  ), stringsAsFactors = FALSE
)
write_csv(config, "V0LV_CONFIG_CANDIDATE2.csv")

objective <- data.frame(
  method_family = c("multiFSYNC", "bayesSYNC_reference"),
  objective_field = c("ordinary_T1_ELBO", "ELBO_iter"),
  eligibility_or_validity = c(
    paste(c("annealing_complete", "ordinary_T1_entered", "finite_ELBO",
            "finite_components", "objective_valid", "no_T1_decrease",
            "no_T1_jitter", "hash_valid"), collapse = ";"),
    paste(c("no_error", "ordinary_T1_entered", "finite_ELBO_iter",
            "finite_primary_outputs", "terminal_and_hash_complete"), collapse = ";")
  ),
  strict_practical_convergence_is_separate = TRUE,
  terminal_error_is_ineligible = TRUE,
  within_model_selection = c("maximum eligible endpoint ELBO",
                             "single preregistered endpoint"),
  cross_model_elbo_comparison = FALSE, stringsAsFactors = FALSE
)
write_csv(objective, "OBJECTIVE_VALIDITY_CONTRACTS.csv")

racing <- data.frame(
  screen_T1 = 100L, minimum_keep = 4L, elbo_margin = 500,
  input_trace = "real_trace_only_through_actual_terminal",
  early_strict_endpoints = "completed_candidates_preserved_with_real_endpoint_ELBO",
  diagnostics_after_actual_terminal_fabricated = FALSE,
  winner_retention_reported = TRUE, false_elimination_reported = TRUE,
  estimated_sweeps_saved_reported = TRUE,
  changes_fit_eligibility_winner_or_success_gate = FALSE,
  stringsAsFactors = FALSE
)
write_csv(racing, "R_OFFLINE_RACING_CONTRACT.csv")

analysis <- data.frame(
  analysis_id = c("A_minus_B0", "R_minus_A", "R_minus_B0",
                  "protocol_level_winners", "R_basin_support",
                  "R_budget_audit", "pooled_reference"),
  unit = c(rep("data_id_then_core_seed_pair", 3L), "data_id", "data_id",
           "audit_data_id_and_R_start", "data_id"),
  primary_summary = c(rep("within_data_median_of_four_paired_differences", 3L),
                      "one_preregistered_winner_per_method_data",
                      "independent_start_support_counts",
                      "fixed_horizon_incumbent_and_counterfactual_trajectories",
                      "single_preregistered_pooled_endpoint"),
  bootstrap_cluster = "data_id", bootstrap_B = 9999L,
  bootstrap_seed = 82640001L, truth_used_before_selection = FALSE,
  stringsAsFactors = FALSE
)
write_csv(analysis, "V0LV_ANALYSIS_PLAN_CANDIDATE2.csv")

success <- data.frame(
  criterion_id = paste0("C", seq_len(8L)),
  level = c("hard_gate", "scientific", "scientific", "scientific",
            "budget", "diagnostic", "production", "baseline_integrity"),
  criterion = c(
    "All data have 12 R terminal rows and at least one eligible R endpoint",
    "At least 8 of 10 R winners recover complete correct 1/1/1 structure",
    "Core-paired R missing plus misplaced is median-no-worse than A and B0",
    "Core-paired reconstruction and feature excess ISE are median-no-worse than A and B0",
    "03 and 08 fixed-horizon audits are complete without replacing main endpoints",
    "Independent best-basin and endpoint-alignment support are fully reported",
    "Production-ready label requires zero unfinished selected R endpoints",
    "All ten pooled terminal records exist; errors and warnings are not imputed"
  ),
  operational_rule = c(
    "12 terminal rows/data;>=1 eligible/data;10 truth-free winners",
    ">=8/10 after authorized unseal",
    "median data-level paired-median difference<=0; worsening count descriptive",
    "each registered metric and contrast median difference<=0; bootstrap and sign test reported",
    "72 audit terminal rows; fixed 200/400/800 or explicit main-error unavailability",
    "counts and identities reported; no unsupported production claim",
    "0 selected_endpoint_unfinished for production-ready label",
    "10 terminal records; validity rate and every warning/error reported"
  ),
  invented_noninferiority_margin = FALSE,
  worsening_count_independent_gate = FALSE,
  selection_uses_criterion = FALSE,
  truth_reveal_required = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
write_csv(success, "V0LV_SUCCESS_CRITERIA_CANDIDATE2.csv")

pooled <- data.frame(
  contract = c("input", "fit", "oracle_match", "primary_candidates",
               "active_study_metrics", "inactive_study_leakage", "N_A",
               "factor_ppi", "ELBO"),
  rule = c(
    "two studies concatenated once; no study labels",
    "Q=3;M=2;single preregistered fit",
    "one-to-one loading-direction match to span[A0,B10,B20]",
    "all three oracle-matched candidates regardless of factor PPI",
    "shared both studies;specific only its true active study",
    "RMS and inactive/(active+inactive) contribution energy",
    "misplaced and block purity are NA with reason; never zero/pass/fail",
    ">=0.5 auxiliary retention rule only",
    "within pooled fit only; never compared to multiFSYNC ELBO"
  ), stringsAsFactors = FALSE
)
write_csv(pooled, "POOLED_EVALUATION_CONTRACT.csv")

seed_rows <- list(); seed_index <- 0L
add_seed <- function(seed, namespace, data_id, path_id, use_id, policy) {
  seed_index <<- seed_index + 1L
  seed_rows[[seed_index]] <<- data.frame(
    seed = as.integer(seed), seed_namespace = namespace, data_id = data_id,
    path_id = path_id, use_id = use_id, reuse_policy = policy,
    formal_use_started = FALSE, stringsAsFactors = FALSE
  )
}
for (row in seq_len(nrow(data_manifest))) add_seed(
  data_manifest$data_seed[[row]], "formal_data_generation",
  data_manifest$data_id[[row]], "generator", data_manifest$data_id[[row]],
  "unique_global_data_seed"
)
for (row in seq_len(nrow(fit_manifest))) {
  fit <- fit_manifest[row, ]
  policy <- if (fit$paired_core_seed) {
    "intentional_B0_A_R_core_pair_within_data"
  } else "unique_within_path_and_data"
  add_seed(fit$fit_seed, "formal_fit", fit$data_id, fit$path_id,
           fit$fit_id, policy)
}
for (row in seq_len(nrow(audit_manifest))) {
  audit <- audit_manifest[row, ]
  add_seed(audit$continuation_seed, "formal_audit_continuation", audit$data_id,
           "R_audit", audit$audit_id,
           "required_reuse_of_original_R_seed_and_complete_state")
}
seed_registry <- do.call(rbind, seed_rows)
write_csv(seed_registry, "V0LV_SEED_REGISTRY_CANDIDATE2.csv")

collision <- do.call(rbind, lapply(split(seed_registry, seed_registry$seed), function(group) {
  if (nrow(group) == 1L) {
    intentional <- TRUE
    collision_class <- "unique"
  } else {
    fit_group <- group[group$seed_namespace == "formal_fit", , drop = FALSE]
    audit_group <- group[
      group$seed_namespace == "formal_audit_continuation", , drop = FALSE
    ]
    generation_group <- group[
      group$seed_namespace == "formal_data_generation", , drop = FALSE
    ]
    fit_paths <- sort(as.character(fit_group$path_id))
    valid_core <- identical(fit_paths, c("A", "B0", "R")) &&
      all(fit_group$reuse_policy ==
            "intentional_B0_A_R_core_pair_within_data")
    valid_extra_r <- identical(fit_paths, "R") &&
      all(fit_group$reuse_policy == "unique_within_path_and_data")
    valid_audit <- nrow(audit_group) %in% c(0L, 3L) &&
      (nrow(audit_group) == 0L || (
        "R" %in% fit_paths &&
          all(audit_group$reuse_policy ==
                "required_reuse_of_original_R_seed_and_complete_state") &&
          all(c("anchor_200", "to_400", "to_800") %in%
                sub("^.*_(anchor_200|to_400|to_800)$", "\\1",
                    audit_group$use_id))
      ))
    intentional <- !nrow(generation_group) &&
      length(unique(group$data_id)) == 1L &&
      (valid_core || valid_extra_r) && valid_audit
    collision_class <- if (valid_core && nrow(audit_group) == 3L) {
      "core_B0_A_R_pair_plus_three_R_audit_stages"
    } else if (valid_core) {
      "core_B0_A_R_pair"
    } else if (valid_extra_r && nrow(audit_group) == 3L) {
      "one_extra_R_start_plus_three_audit_stages"
    } else if (valid_extra_r) {
      "one_extra_R_start"
    } else "unregistered_collision"
  }
  data.frame(seed = group$seed[[1L]], uses = nrow(group),
             data_ids = paste(unique(group$data_id), collapse = ";"),
             use_ids = paste(group$use_id, collapse = ";"),
             collision_class = collision_class,
             intentional_or_unique = intentional,
             stringsAsFactors = FALSE)
}))
if (any(!collision$intentional_or_unique)) stop("Unregistered seed collision.")
write_csv(collision, "SEED_COLLISION_AUDIT.csv")

version_lineage <- data.frame(
  version_line = c("R-only", "v0L-V total", "manifest"),
  parent = c("R-only v2 draft and v3 candidate",
             "v0L-V v1 draft and v2 candidate", "manifest v2 and v3 candidate"),
  candidate2 = c(
    "PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V3_CANDIDATE2_20260812.md",
    "PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812.md",
    "protocol_manifests_20260812_v3_candidate2"
  ), status = "candidate_not_frozen", overwrote_parent = FALSE,
  stringsAsFactors = FALSE
)
write_csv(version_lineage, "VERSION_LINEAGE.csv")

parent_paths <- c(
  file.path(base_dir, c(
    "PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V2_DRAFT_20260812.md",
    "PROTOCOL_V0L_V_UNSEEN_VALIDATION_V1_DRAFT_20260812.md",
    "PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V3_CANDIDATE_20260812.md",
    "PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE_20260812.md"
  )),
  list.files(file.path(base_dir, "protocol_manifests_20260812_v3_candidate"),
             recursive = TRUE, full.names = TRUE)
)
if (any(!file.exists(parent_paths))) stop("A required parent artifact is missing.")
parent_hashes <- c2_hash_table(
  parent_paths, project_root, role = "immutable_parent_or_prior_candidate"
)
write_csv(parent_hashes, "PARENT_HISTORY_HASHES.csv")

source_paths <- c(
  list.files(file.path(project_root, "Rcode", "multiFSYNC", "R"),
             pattern = "[.]R$", full.names = TRUE),
  file.path(project_root, "Rcode", "multiFSYNC", c("DESCRIPTION", "NAMESPACE")),
  list.files(file.path(project_root, "Rcode", "bayesSYNC_ref", "R"),
             pattern = "[.]R$", full.names = TRUE),
  file.path(base_dir, c(
    "v0lv_candidate2_runtime.R", "r_route_v3_candidate2_runner.R",
    "v0lv_candidate2_controller.R", "v0lv_candidate2_evaluation.R",
    "run_v0lv_candidate2_formal.R", "run_v0lv_candidate2_e2e_smoke.R",
    "test_v0lv_candidate2.R", "build_protocol_manifests_v3_candidate2.R",
    "run_candidate2_final_verification.R",
    "run_stage_a_round2_package_tests.R",
    "PROTOCOL_R_ONLY_MULTISTART_VALIDATION_V3_CANDIDATE2_20260812.md",
    "PROTOCOL_V0L_V_UNSEEN_VALIDATION_V2_CANDIDATE2_20260812.md"
  )),
  file.path(project_root, "对比实验", c(
    file.path("v0F维数信息盆地桥接_20260801", "v0f_data_calibration.R"),
    file.path("v0G中性缩减重拟合_20260803", "v0g_evaluation.R"),
    file.path("v0H载荷优先结构评估_20260803", "v0h_evaluation.R"),
    file.path("v0I真值初始化短诊断_20260805", "v0i_truth_diagnostic.R"),
    file.path("v0J无真值可达性析因_20260808", "v0j_evaluation.R")
  ))
)
missing_source <- source_paths[!file.exists(source_paths)]
if (length(missing_source)) {
  stop("Candidate2 source binding input is missing: ",
       paste(missing_source, collapse = ", "))
}
role <- ifelse(grepl("Rcode/multiFSYNC", source_paths, fixed = TRUE),
               "multiFSYNC_complete_package_source",
        ifelse(grepl("Rcode/bayesSYNC_ref", source_paths, fixed = TRUE),
               "bayesSYNC_reference_source", "controller_protocol_evaluator"))
source_index <- c2_hash_table(source_paths, project_root, role = role)
write_csv(source_index, "FULL_RUNTIME_SOURCE_BINDINGS.csv")

key_packages <- c(
  "multiFSYNC", "ellipse", "gtools", "magic", "MASS", "matrixcalc",
  "matrixStats", "pracma", "splines", "testthat", "pkgload", "desc"
)
environment <- c2_current_environment_binding(key_packages)
write_csv(environment, "ENVIRONMENT_BINDING.csv")

if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
  stop("Install the current multiFSYNC source before building namespace bindings.")
}
namespace <- asNamespace("multiFSYNC")
namespace_names <- sort(ls(namespace, all.names = TRUE)[vapply(
  ls(namespace, all.names = TRUE),
  function(name) is.function(get(name, envir = namespace, inherits = FALSE)),
  logical(1)
)])
namespace_bindings <- data.frame(
  function_name = namespace_names,
  sha256 = vapply(namespace_names, function(name) {
    c2_function_signature_sha256(get(name, envir = namespace, inherits = FALSE))
  }, character(1)), stringsAsFactors = FALSE
)
write_csv(namespace_bindings, "MULTIFSYNC_NAMESPACE_FUNCTION_BINDINGS.csv")

runtime_functions <- data.frame(
  function_name = c("multiFSYNC::bayesSYNC_multi",
                    "multiFSYNC::bayesSYNC_multi_pre_score"),
  sha256 = c(
    c2_function_signature_sha256(getExportedValue("multiFSYNC", "bayesSYNC_multi")),
    c2_function_signature_sha256(
      getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score")
    )
  ), stringsAsFactors = FALSE
)
write_csv(runtime_functions, "RUNTIME_FUNCTION_BINDINGS.csv")

binding_status <- data.frame(
  manifest_schema_id = MANIFEST_ID, status = "candidate",
  formal_execution_authorized = FALSE, frozen_bundle_created = FALSE,
  formal_data_generated = FALSE, formal_fits_started = FALSE,
  formal_continuations_started = FALSE,
  requires_separate_user_freeze_authorization = TRUE,
  requires_separate_user_start_authorization = TRUE,
  stringsAsFactors = FALSE
)
write_csv(binding_status, "BINDING_STATUS.csv")

source_verification <- c2_verify_hash_index(source_index, project_root)
integrity_checks <- data.frame(
  check_id = c(
    "ten_data_rows", "two_hundred_ten_fit_rows", "method_row_counts",
    "seventy_two_audit_rows", "unique_fit_and_audit_ids",
    "all_formal_start_flags_false", "truth_isolation_flags_false",
    "seed_collision_contract", "complete_source_inventory",
    "source_hash_self_check", "candidate_not_frozen_or_authorized",
    "parent_history_present"
  ),
  passed = c(
    nrow(data_manifest) == 10L,
    nrow(fit_manifest) == 210L,
    identical(as.integer(table(factor(
      fit_manifest$path_id,
      levels = c("B0", "A", "R", "pooled_bayesSYNC")
    ))), c(40L, 40L, 120L, 10L)),
    nrow(audit_manifest) == 72L,
    !anyDuplicated(fit_manifest$fit_id) &&
      !anyDuplicated(audit_manifest$audit_id),
    !any(fit_manifest$formal_fit_started) &&
      !any(audit_manifest$formal_continuation_started) &&
      !any(data_manifest$formal_generation_started),
    !any(fit_manifest$truth_available_to_fit) &&
      !any(fit_manifest$truth_available_to_stopping) &&
      !any(fit_manifest$truth_available_to_selection) &&
      !any(audit_manifest$truth_available_to_continuation) &&
      !any(audit_manifest$truth_available_to_counterfactual_selection),
    all(collision$intentional_or_unique),
    sum(source_index$role == "multiFSYNC_complete_package_source") == 28L &&
      sum(source_index$role == "bayesSYNC_reference_source") == 5L &&
      !anyDuplicated(source_index$relative_path),
    all(source_verification$matches),
    binding_status$status == "candidate" &&
      !binding_status$formal_execution_authorized &&
      !binding_status$frozen_bundle_created,
    length(parent_paths) > 4L && all(file.exists(parent_paths))
  ),
  detail = c(
    nrow(data_manifest), nrow(fit_manifest),
    paste(names(table(fit_manifest$path_id)),
          as.integer(table(fit_manifest$path_id)), collapse = ";"),
    nrow(audit_manifest), "fit_id and audit_id unique",
    "generation/fit/continuation starts all FALSE",
    "fit/stopping/selection/continuation truth flags all FALSE",
    paste(table(collision$collision_class), collapse = ";"),
    paste0("source rows=", nrow(source_index),
           ";multiFSYNC=28;bayesSYNC_ref=5"),
    paste0("verified=", sum(source_verification$matches), "/",
           nrow(source_verification)),
    "status=candidate;authorized=FALSE;frozen=FALSE",
    paste0("parent files=", length(parent_paths))
  ), stringsAsFactors = FALSE
)
if (any(!integrity_checks$passed)) {
  stop("Candidate2 integrity checks failed: ", paste(
    integrity_checks$check_id[!integrity_checks$passed], collapse = ", "
  ))
}
write_csv(integrity_checks, "CANDIDATE_INTEGRITY_CHECKS.csv")

artifact_paths <- list.files(output_dir, full.names = TRUE)
artifact_paths <- artifact_paths[!dir.exists(artifact_paths)]
artifact_hashes <- c2_hash_table(artifact_paths, project_root,
                                 role = "candidate2_manifest_artifact")
artifact_write <- write_csv(artifact_hashes, "CANDIDATE_ARTIFACT_HASHES.csv")
c2_atomic_write_lines(c(
  "V0LV_PROTOCOL_MANIFEST_V3_CANDIDATE2_COMPLETE",
  paste0("manifest_schema_id=", MANIFEST_ID),
  paste0("total_protocol_id=", TOTAL_ID), paste0("R_only_protocol_id=", R_ID),
  "status=candidate_not_frozen", "formal_execution_authorized=FALSE",
  "formal_data_generated=FALSE", "formal_fits_started=FALSE",
  "formal_continuations_started=FALSE", "frozen_bundle_created=FALSE",
  "formal_fit_rows=210", "formal_audit_rows=72",
  paste0("artifact_index_sha256=", artifact_write$sha256)
), file.path(output_dir, "CANDIDATE_COMPLETE.txt"))

cat("candidate2_manifest_built=", normalizePath(output_dir, winslash = "/"),
    "\nformal_work_started=FALSE\n", sep = "")
