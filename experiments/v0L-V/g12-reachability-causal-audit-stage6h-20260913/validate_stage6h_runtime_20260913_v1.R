#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_root <- dirname(sub("^--file=", "", script_argument[[1L]]))
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(
    sub("^--file=", "", script_argument[[1L]]), winslash = "/",
    mustWork = TRUE
  ))
}
experiments_root <- dirname(script_root)
source(file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)
source(file.path(
  script_root, "stage6h_evaluation_helpers_20260913_v1.R"
), local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6H_REACHABILITY_CAUSAL_AUDIT_V1_20260913"
stage6h_library <- pp_env("STAGE6H_R_LIB")
if (nzchar(stage6h_library)) Sys.setenv(STAGE6A_R_LIB = stage6h_library)
stage6g_output_root <- pp_env(
  "STAGE6H_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
pp_load_package()
pp_source_runtime(TRUE)

observation <- readRDS(file.path(
  stage6g_output_root, "data", "s6g_minimal_L111_M2",
  "observation_bundle.rds"
))
forbidden <- character()
walk <- function(object) {
  object_names <- names(object)
  if (!is.null(object_names)) forbidden <<- c(
    forbidden,
    object_names[grepl("(^|_)(truth|true)(_|$)|signal_true|noise_true",
                       object_names, ignore.case = TRUE)]
  )
  if (is.list(object)) for (item in object) walk(item)
}
walk(observation)
pp_assert(!length(unique(forbidden)),
          "Stage 6H smoke observation contains truth.")
row <- data.frame(
  fit_id = "stage6h_runtime_micro_smoke", data_id = "s6g_minimal_L111_M2",
  scenario_id = "s6g_minimal_L111_M2", fit_config_id = "cand_L111_M2",
  selection_stratum_id = "stage6h_runtime_micro_smoke",
  method_id = "G12", route_id = "random__pre1__jaoua_anneal__G12",
  seed_index = 12L, fit_seed = 84611112L,
  fit_L_f = 1L, fit_L_s_1 = 1L, fit_L_s_2 = 1L,
  fit_M_f = "2", fit_M_s_1 = "2", fit_M_s_2 = "2",
  stringsAsFactors = FALSE
)
smoke <- pp_fit_g12(observation, row, smoke = TRUE)
pp_assert(!is.null(smoke$fit) && isTRUE(smoke$objective$eligible) &&
            identical(smoke$fit$pre_score_interface$function_initialization,
                      "gram_unit_energy") &&
            identical(as.integer(smoke$fit$n_cpus_used), 1L),
          "Stage 6H G12 micro fit failed.")

existing_fit_id <-
  "s6g_minimal_L111_M2__split1_middle__cand_L111_M2__G12_01"
existing_fit <- readRDS(file.path(
  stage6g_output_root, "fits", existing_fit_id, "fit.rds"
))
truth <- readRDS(file.path(
  stage6g_output_root, "data", "s6g_minimal_L111_M2", "sealed_truth",
  "truth_bundle.rds"
))
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
evaluation <- s6h_evaluate_fit(existing_fit, truth, evaluator_environment)
expected <- c(
  "feature_component_errors", "factor_kernel_errors",
  "dense_signal_metrics", "loading_recovery", "threshold_sensitivity",
  "unthresholded_assignment", "covariance_operator_error",
  "ppi_diagnostics", "shared_candidate_metrics", "summary",
  "evaluation_contract", "process_contribution_by_role",
  "process_contribution_summary", "loading_scale_by_block",
  "loading_scale_summary", "candidate_activity",
  "candidate_activity_summary", "fpca_rank_activity",
  "fpca_component_activity"
)
pp_assert(all(expected %in% names(evaluation)) &&
            all(vapply(evaluation[expected], function(value) {
              is.data.frame(value) && nrow(value) > 0L
            }, logical(1L))),
          "Stage 6H evaluator micro smoke failed.")
cat("STAGE6H_RUNTIME_MICRO_SMOKE_PASS fit_route=1 evaluator_tables=19\n")
