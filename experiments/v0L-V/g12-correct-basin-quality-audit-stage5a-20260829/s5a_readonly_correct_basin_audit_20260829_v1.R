#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "Rcode", "multiFSYNC", "DESCRIPTION")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot locate the multiFSYNC project root.", call. = FALSE)
    }
    current <- parent
  }
}
project_root <- find_project_root()

find_named_directory <- function(root, target, maximum_depth = 4L) {
  frontier <- normalizePath(root, winslash = "/", mustWork = TRUE)
  for (depth in seq_len(maximum_depth)) {
    children <- unique(unlist(lapply(frontier, function(path) {
      listed <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      listed[listed != path]
    }), use.names = FALSE))
    hits <- children[basename(children) == target]
    if (length(hits) == 1L) {
      return(normalizePath(hits, winslash = "/", mustWork = TRUE))
    }
    if (length(hits) > 1L) {
      stop(paste0("Multiple directories named ", target, "."), call. = FALSE)
    }
    frontier <- children
  }
  stop(paste0("Cannot locate directory ", target, "."), call. = FALSE)
}

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
independent_root <- env_or(
  "STAGE5A_INDEPENDENT_ROOT",
  find_named_directory(
    project_root,
    "v0lv-g12-independent-confirmation-20260822-v1-full_archive"
  )
)
cross_root <- env_or(
  "STAGE5A_CROSS_ROOT",
  find_named_directory(
    project_root,
    "v0lv-g12-cross-scenario-confirmation-20260823-v1"
  )
)
snapshot_root <- env_or(
  "STAGE5A_SNAPSHOT_ROOT",
  file.path(project_root, "Rcode", "multiFSYNC", "experiments",
            "v0L-V", "rc1_snapshot")
)
result_root <- env_or(
  "STAGE5A_RESULT_ROOT",
  file.path(
    dirname(dirname(independent_root)),
    "v0L-V_G12_stage5A_correct_basin_quality_20260829",
    "g12-correct-basin-quality-audit-stage5a-20260829-v1"
  )
)
output <- file.path(result_root, "analysis_20260829_v1")
check_only <- identical(Sys.getenv("STAGE5A_CHECK_ONLY", unset = "0"), "1")

abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(value, message) if (!isTRUE(value)) abort(message)
read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing input: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(value, name) {
  utils::write.csv(
    value, file.path(output, name), row.names = FALSE,
    quote = TRUE, na = ""
  )
}
safe_mean <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
safe_median <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) stats::median(value) else NA_real_
}
safe_min <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) min(value) else NA_real_
}
safe_max <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) max(value) else NA_real_
}
safe_quantile <- function(value, probability) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) {
    unname(stats::quantile(value, probability, names = FALSE, type = 7L))
  } else NA_real_
}
vector_l2 <- function(value) sqrt(sum(as.numeric(value)^2))
vector_rms <- function(value) sqrt(mean(as.numeric(value)^2))
abs_cosine <- function(left, right) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (length(left) != length(right) || !length(left) ||
      any(!is.finite(c(left, right)))) return(NA_real_)
  denominator <- vector_l2(left) * vector_l2(right)
  if (!is.finite(denominator) || denominator <= 0) return(0)
  abs(sum(left * right) / denominator)
}
nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
  assert(length(estimate) == length(truth) && length(truth) > 0L,
         "NRMSE vectors are not aligned.")
  denominator <- sqrt(mean(truth^2))
  assert(is.finite(denominator) && denominator > 0,
         "Truth NRMSE denominator is invalid.")
  sqrt(mean((estimate - truth)^2)) / denominator
}
iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}
find_one <- function(path, pattern) {
  result <- list.files(
    path, pattern = pattern, recursive = TRUE, full.names = TRUE
  )
  assert(length(result) == 1L,
         paste0("Expected one ", pattern, "; found ", length(result), "."))
  normalizePath(result[[1L]], winslash = "/", mustWork = TRUE)
}

assert(!dir.exists(output), paste0("Refusing to overwrite: ", output))
assert(all(dir.exists(c(independent_root, cross_root, snapshot_root))),
       "A bound Stage-5A input root is missing.")

independent_files <- c(
  science = file.path(independent_root, "evaluation",
                      "ALL_144_SCIENTIFIC_RESULTS.csv"),
  runtime = file.path(independent_root, "analysis_20260823_v1",
                      "ALL_144_RUNTIME_AND_STABILITY.csv"),
  winners = file.path(independent_root, "truth_free_selection",
                      "ALL_TRUTH_FREE_WINNERS.csv"),
  manifest = file.path(independent_root, "FIT_MANIFEST.csv"),
  evaluation_complete = file.path(independent_root, "evaluation",
                                  "EVALUATION_COMPLETE.txt"),
  analysis_complete = file.path(independent_root, "analysis_20260823_v1",
                                "ANALYSIS_COMPLETE.txt")
)
cross_files <- c(
  science = file.path(cross_root, "evaluation",
                      "ALL_96_SCIENTIFIC_RESULTS.csv"),
  runtime = file.path(cross_root, "analysis_20260824_v1",
                      "ALL_96_RUNTIME_AND_STABILITY.csv"),
  winners = file.path(cross_root, "truth_free_selection",
                      "ALL_TRUTH_FREE_WINNERS.csv"),
  manifest = file.path(cross_root, "FIT_MANIFEST.csv"),
  evaluation_complete = file.path(cross_root, "evaluation",
                                  "EVALUATION_COMPLETE.txt"),
  analysis_complete = file.path(cross_root, "analysis_20260824_v1",
                                "ANALYSIS_COMPLETE.txt")
)
assert(all(file.exists(c(independent_files, cross_files))),
       "A registered historical result is incomplete.")

independent_science_all <- read_csv(independent_files[["science"]])
independent_runtime_all <- read_csv(independent_files[["runtime"]])
independent_winners_all <- read_csv(independent_files[["winners"]])
cross_science_all <- read_csv(cross_files[["science"]])
cross_runtime_all <- read_csv(cross_files[["runtime"]])
cross_winners_all <- read_csv(cross_files[["winners"]])

independent <- independent_science_all[
  independent_science_all$method_id == "gram_unit_energy", , drop = FALSE
]
independent$scenario_id <- "baseline_strong"
independent$cohort <- "independent_primary"
cross <- cross_science_all[
  cross_science_all$method_id == "gram_unit_energy" &
    cross_science_all$scenario_id %in%
      c("baseline_strong", "weak_loading_separation"), , drop = FALSE
]
cross$cohort <- "cross_standard_secondary"

append_runtime <- function(science, runtime) {
  match_index <- match(science$fit_id, runtime$fit_id)
  assert(all(!is.na(match_index)), "Runtime join is incomplete.")
  auxiliary <- c(
    "auxiliary_ever_stable_5", "auxiliary_maximum_streak",
    "auxiliary_endpoint_pass", "auxiliary_endpoint_stable_5"
  )
  assert(all(auxiliary %in% names(runtime)),
         "An auxiliary stability field is missing.")
  for (name in auxiliary) science[[name]] <- runtime[[name]][match_index]
  science
}
independent <- append_runtime(independent, independent_runtime_all)
cross <- append_runtime(cross, cross_runtime_all)

independent_winners <- independent_winners_all[
  independent_winners_all$method_id == "gram_unit_energy", , drop = FALSE
]
cross_winners <- cross_winners_all[
  cross_winners_all$method_id == "gram_unit_energy" &
    cross_winners_all$scenario_id %in%
      c("baseline_strong", "weak_loading_separation"), , drop = FALSE
]
winner_ids <- c(independent_winners$fit_id, cross_winners$fit_id)

common_columns <- c(
  "fit_id", "data_id", "scenario_id", "method_id", "seed_index",
  "fit_seed", "terminal_status", "actual_T1_sweeps",
  "strict_practical_converged", "objective_eligible", "final_elbo",
  "warning_count", "elapsed_seconds", "peak_memory_bytes",
  "truth_used_for_fit_or_selection", "formal_v0lv_result",
  "selected_counts", "global_correct", "global_misplaced",
  "global_missing", "global_extra", "global_duplicate",
  "study_correct", "study_misplaced", "study_missing",
  "study_candidate_extra", "study_candidate_duplicate",
  "strict_joint_recovery", "strict_joint_correct_type",
  "strict_joint_misplaced", "matched_function_recall_mean",
  "matched_trajectory_abs_cor_mean", "matched_score_abs_cor_mean",
  "signal_nrmse_ppi_selected", "signal_nrmse_mean_only", "P_A", "P_B1",
  "P_B2", "dense_signal_nrmse_ppi_selected",
  "dense_signal_relative_mise_ppi_selected", "feature_components_total",
  "feature_components_matched", "feature_total_ise_mean",
  "feature_projection_floor_ise_mean",
  "feature_estimate_to_projection_ise_mean", "factor_kernels_total",
  "factor_kernels_matched", "kernel_relative_ise_mean",
  "covariance_operator_relative_error_mean",
  "loading_relative_l2_error_mean", "loading_support_ppi_auc_mean",
  "exact_correct_1_1_1", "auxiliary_ever_stable_5",
  "auxiliary_maximum_streak", "auxiliary_endpoint_pass",
  "auxiliary_endpoint_stable_5", "cohort"
)
assert(all(common_columns %in% names(independent)) &&
         all(common_columns %in% names(cross)),
       "The two cohorts do not share the required frozen evaluation schema.")
endpoints <- rbind(
  independent[, common_columns, drop = FALSE],
  cross[, common_columns, drop = FALSE]
)
endpoints$is_frozen_winner <- endpoints$fit_id %in% winner_ids
endpoints$structure_class <- ifelse(
  endpoints$exact_correct_1_1_1, "exact_identity",
  ifelse(endpoints$selected_counts == "1,1,1",
         "count_correct_identity_wrong", "incomplete_or_other")
)

assert(nrow(independent) == 72L && all(table(independent$data_id) == 12L),
       "Primary G12 cohort is not the frozen 72 endpoints.")
assert(nrow(cross) == 32L && all(table(cross$data_id) == 8L),
       "Secondary standard-density G12 cohort is not 32 endpoints.")
assert(nrow(endpoints) == 104L && length(unique(endpoints$fit_id)) == 104L,
       "Stage-5A endpoint identities are not 104 unique fits.")
assert(sum(endpoints$exact_correct_1_1_1) == 45L,
       "The frozen exact-endpoint count is not 45.")
assert(length(winner_ids) == 10L && sum(endpoints$is_frozen_winner) == 10L &&
         all(endpoints$exact_correct_1_1_1[endpoints$is_frozen_winner]),
       "The ten frozen G12 winners are not intact and exact.")
assert(all(endpoints$objective_eligible) &&
         !any(endpoints$truth_used_for_fit_or_selection) &&
         !any(endpoints$formal_v0lv_result),
       "Historical truth isolation or eligibility has changed.")

primary_fit_paths <- file.path(
  independent_root, "fits", independent$fit_id, "fit.rds"
)
primary_evaluation_paths <- file.path(
  independent_root, "evaluation", independent$fit_id, "evaluation.rds"
)
primary_truth_paths <- file.path(
  independent_root, "data", independent$data_id,
  "sealed_truth", "truth_bundle.rds"
)
assert(all(file.exists(c(primary_fit_paths, primary_evaluation_paths,
                         primary_truth_paths))),
       "A primary-cohort fit/evaluation/truth RDS is missing.")

if (check_only) {
  cat(paste0(
    "STAGE5A_CHECK_ONLY_PASS endpoints=", nrow(endpoints),
    " exact=", sum(endpoints$exact_correct_1_1_1),
    " winners=", sum(endpoints$is_frozen_winner),
    " primary_rds=", length(primary_fit_paths), "\n"
  ))
  quit(save = "no", status = 0L)
}

wrapper <- env_or(
  "STAGE5A_EVALUATOR_WRAPPER",
  find_one(snapshot_root, "^v0lv_candidate2_evaluation[.]R$")
)
assert(file.exists(wrapper), paste0("Missing evaluator wrapper: ", wrapper))
source(
  wrapper, local = globalenv(), encoding = "UTF-8", keep.source = TRUE
)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(snapshot_root)

dir.create(output, recursive = TRUE, mode = "0700")

evaluator_source_paths <- candidate2_evaluation_source_paths(snapshot_root)
input_paths <- c(
  stats::setNames(
    independent_files, paste0("independent_", names(independent_files))
  ),
  stats::setNames(cross_files, paste0("cross_", names(cross_files))),
  evaluator_wrapper = wrapper,
  stats::setNames(
    evaluator_source_paths,
    paste0("evaluator_source_", basename(evaluator_source_paths))
  )
)
input_paths <- input_paths[!duplicated(input_paths)]
path_names <- names(input_paths)
display_paths <- gsub("\\\\", "/", input_paths, fixed = TRUE)
independent_index <- startsWith(path_names, "independent_")
cross_index <- startsWith(path_names, "cross_")
source_index <- startsWith(path_names, "evaluator_source_")
display_paths[independent_index] <- paste0(
  "independent_primary/",
  substring(display_paths[independent_index], nchar(independent_root) + 2L)
)
display_paths[cross_index] <- paste0(
  "cross_standard_secondary/",
  substring(display_paths[cross_index], nchar(cross_root) + 2L)
)
display_paths[path_names == "evaluator_wrapper"] <-
  "rc1_snapshot/evaluator_wrapper/v0lv_candidate2_evaluation.R"
display_paths[source_index] <- paste0(
  "rc1_snapshot/frozen_evaluator_source/", basename(input_paths[source_index])
)
input_binding <- data.frame(
  input_name = path_names,
  logical_path = display_paths,
  bytes = unname(file.info(input_paths)$size),
  sha256 = vapply(
    input_paths,
    function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L)
  ),
  stringsAsFactors = FALSE
)
write_csv(input_binding, "INPUT_BINDING.csv")

truth_cache <- new.env(parent = emptyenv())
process_rows <- vector("list", nrow(independent))
process_endpoint_rows <- vector("list", nrow(independent))
scale_rows <- vector("list", nrow(independent))

for (index in seq_len(nrow(independent))) {
  row <- independent[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  data_id <- row$data_id[[1L]]
  fit <- readRDS(file.path(independent_root, "fits", fit_id, "fit.rds"))
  evaluation <- readRDS(file.path(
    independent_root, "evaluation", fit_id, "evaluation.rds"
  ))
  if (!exists(data_id, envir = truth_cache, inherits = FALSE)) {
    assign(
      data_id,
      readRDS(file.path(independent_root, "data", data_id,
                        "sealed_truth", "truth_bundle.rds")),
      envir = truth_cache
    )
  }
  truth <- get(data_id, envir = truth_cache, inherits = FALSE)
  assert(identical(truth$data_id, data_id) &&
           identical(truth$bundle_class, "v0lv_sealed_truth_candidate2"),
         paste0("Truth identity mismatch for ", fit_id, "."))
  truth_data <- truth$data

  reported_paths <-
    evaluator_environment$v0g_reported_paths_by_study(fit, truth_data)
  true_paths <- evaluator_environment$v0h_truth_paths_by_study(truth_data)
  functional <- evaluation$primary$functional_recovery
  functional <- functional[functional$scope == "ppi_selected", , drop = FALSE]
  assert(nrow(functional) == 4L,
         paste0("Expected four study-role truth rows for ", fit_id, "."))

  endpoint_process <- lapply(seq_len(nrow(functional)), function(item_index) {
    item <- functional[item_index, , drop = FALSE]
    study <- as.integer(item$study[[1L]])
    true_block <- item$true_role[[1L]]
    true_factor <- as.integer(item$true_factor[[1L]])
    true_loading <- if (true_block == "shared") {
      truth_data$true_params$a_true[, true_factor]
    } else {
      truth_data$true_params$b_true[[study]][, true_factor]
    }
    true_path <- true_paths[[true_block]][[study]][, true_factor]
    true_contribution <- as.numeric(outer(true_path, true_loading))
    matched <- isTRUE(item$matched_by_loading[[1L]])
    if (matched) {
      estimated_block <- item$estimated_block[[1L]]
      estimated_factor <- as.integer(item$estimated_factor[[1L]])
      estimated_loading <- if (estimated_block == "shared") {
        fit$mu_q_a_hat[, estimated_factor]
      } else {
        fit$mu_q_b_specific_hat[[study]][, estimated_factor]
      }
      estimated_path <-
        reported_paths[[estimated_block]][[study]][, estimated_factor]
      aligned_path <- item$loading_sign[[1L]] * estimated_path
      process_error <- nrmse(aligned_path, true_path)
      contribution_error <- nrmse(
        as.numeric(outer(estimated_path, estimated_loading)),
        true_contribution
      )
      loading_norm_ratio <-
        vector_l2(estimated_loading) / vector_l2(true_loading)
      loading_direction_cosine <-
        abs_cosine(estimated_loading, true_loading)
      estimated_path_rms <- vector_rms(estimated_path)
      true_path_rms <- vector_rms(true_path)
    } else {
      estimated_block <- NA_character_
      estimated_factor <- NA_integer_
      process_error <- contribution_error <- loading_norm_ratio <-
        loading_direction_cosine <- estimated_path_rms <- NA_real_
      true_path_rms <- vector_rms(true_path)
    }
    data.frame(
      fit_id = fit_id, data_id = data_id,
      seed_index = row$seed_index[[1L]],
      is_frozen_winner = fit_id %in% winner_ids,
      structure_exact = row$exact_correct_1_1_1[[1L]],
      study = study, true_role = true_block,
      true_id = item$true_id[[1L]],
      estimated_block = estimated_block,
      estimated_factor = estimated_factor,
      loading_structure_status = item$loading_structure_status[[1L]],
      matched_by_loading = matched,
      loading_norm_ratio = loading_norm_ratio,
      loading_abs_cosine = loading_direction_cosine,
      true_path_rms = true_path_rms,
      estimated_path_rms = estimated_path_rms,
      factor_process_nrmse_matched = process_error,
      factor_process_nrmse_missing_as_zero =
        if (matched) process_error else 1,
      complete_contribution_nrmse_matched = contribution_error,
      complete_contribution_nrmse_missing_as_zero =
        if (matched) contribution_error else 1,
      stringsAsFactors = FALSE
    )
  })
  endpoint_process <- do.call(rbind, endpoint_process)
  process_rows[[index]] <- endpoint_process
  process_endpoint_rows[[index]] <- data.frame(
    fit_id = fit_id, data_id = data_id,
    seed_index = row$seed_index[[1L]],
    is_frozen_winner = fit_id %in% winner_ids,
    structure_exact = row$exact_correct_1_1_1[[1L]],
    study_role_matched = sum(endpoint_process$matched_by_loading),
    study_role_truth_total = nrow(endpoint_process),
    factor_process_nrmse_mean_matched =
      safe_mean(endpoint_process$factor_process_nrmse_matched),
    factor_process_nrmse_mean_missing_as_zero =
      mean(endpoint_process$factor_process_nrmse_missing_as_zero),
    complete_contribution_nrmse_mean_matched =
      safe_mean(endpoint_process$complete_contribution_nrmse_matched),
    complete_contribution_nrmse_mean_missing_as_zero =
      mean(endpoint_process$complete_contribution_nrmse_missing_as_zero),
    stringsAsFactors = FALSE
  )

  role_rows <- lapply(c("A", "B1", "B2"), function(role) {
    if (role == "A") {
      raw_loading <- fit$mu_q_a_original[, 1L]
      canonical_loading <- fit$mu_q_a_hat[, 1L]
      true_loading <- truth_data$true_params$a_true[, 1L]
      factor_scale <- fit$factor_scale_shared[[1L]]
      factor_ppi <- fit$factor_ppi_shared[[1L]]
    } else {
      study <- if (role == "B1") 1L else 2L
      raw_loading <- fit$mu_q_b_specific_original[[study]][, 1L]
      canonical_loading <- fit$mu_q_b_specific_hat[[study]][, 1L]
      true_loading <- truth_data$true_params$b_true[[study]][, 1L]
      factor_scale <- fit$factor_scale_specific[[study]][[1L]]
      factor_ppi <- fit$factor_ppi_specific[[study]][[1L]]
    }
    data.frame(
      fit_id = fit_id, data_id = data_id,
      seed_index = row$seed_index[[1L]],
      is_frozen_winner = fit_id %in% winner_ids,
      structure_exact = row$exact_correct_1_1_1[[1L]],
      block_role = role, factor_ppi = factor_ppi,
      factor_scale = factor_scale,
      raw_loading_norm_ratio =
        vector_l2(raw_loading) / vector_l2(true_loading),
      canonical_loading_norm_ratio =
        vector_l2(canonical_loading) / vector_l2(true_loading),
      raw_loading_abs_cosine = abs_cosine(raw_loading, true_loading),
      canonical_loading_abs_cosine =
        abs_cosine(canonical_loading, true_loading),
      stringsAsFactors = FALSE
    )
  })
  scale_rows[[index]] <- do.call(rbind, role_rows)

  rm(fit, evaluation, truth, truth_data, endpoint_process)
  invisible(gc(verbose = FALSE))
}

process <- do.call(rbind, process_rows)
process_endpoint <- do.call(rbind, process_endpoint_rows)
loading_scale <- do.call(rbind, scale_rows)
write_csv(process, "PRIMARY_PROCESS_CONTRIBUTION_BY_ROLE.csv")
write_csv(process_endpoint, "PRIMARY_PROCESS_CONTRIBUTION_BY_ENDPOINT.csv")
write_csv(loading_scale, "PRIMARY_LOADING_SCALE_BY_BLOCK.csv")

process_match <- match(endpoints$fit_id, process_endpoint$fit_id)
endpoints$factor_process_nrmse_mean_missing_as_zero <-
  process_endpoint$factor_process_nrmse_mean_missing_as_zero[process_match]
endpoints$complete_contribution_nrmse_mean_missing_as_zero <-
  process_endpoint$complete_contribution_nrmse_mean_missing_as_zero[
    process_match
  ]
endpoints$study_role_matched <-
  process_endpoint$study_role_matched[process_match]
endpoints$study_role_truth_total <-
  process_endpoint$study_role_truth_total[process_match]

endpoint_order <- order(endpoints$cohort, endpoints$data_id,
                        endpoints$seed_index)
endpoints <- endpoints[endpoint_order, , drop = FALSE]
write_csv(endpoints, "ALL_104_STANDARD_DENSITY_G12_ENDPOINTS.csv")
exact_endpoints <- endpoints[endpoints$exact_correct_1_1_1, , drop = FALSE]
write_csv(exact_endpoints, "EXACT_45_ENDPOINT_QUALITY.csv")
write_csv(endpoints[endpoints$is_frozen_winner, , drop = FALSE],
          "FROZEN_10_WINNER_QUALITY.csv")

data_rows <- lapply(split(endpoints, endpoints$data_id), function(group) {
  exact <- group[group$exact_correct_1_1_1, , drop = FALSE]
  incorrect <- group[!group$exact_correct_1_1_1, , drop = FALSE]
  winner <- group[group$is_frozen_winner, , drop = FALSE]
  assert(nrow(winner) == 1L && nrow(exact) >= 1L,
         "Each Stage-5A data set must have one exact frozen winner.")
  data.frame(
    cohort = group$cohort[[1L]], scenario_id = group$scenario_id[[1L]],
    data_id = group$data_id[[1L]], endpoints = nrow(group),
    exact_endpoints = nrow(exact), exact_rate = nrow(exact) / nrow(group),
    count_correct_identity_wrong =
      sum(group$structure_class == "count_correct_identity_wrong"),
    incomplete_or_other =
      sum(group$structure_class == "incomplete_or_other"),
    frozen_winner_fit_id = winner$fit_id[[1L]],
    frozen_winner_exact = winner$exact_correct_1_1_1[[1L]],
    best_exact_elbo = max(exact$final_elbo),
    best_incorrect_elbo = if (nrow(incorrect)) {
      max(incorrect$final_elbo)
    } else NA_real_,
    exact_over_incorrect_elbo_gap = if (nrow(incorrect)) {
      max(exact$final_elbo) - max(incorrect$final_elbo)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
})
per_data_basin <- do.call(rbind, data_rows)
row.names(per_data_basin) <- NULL
write_csv(per_data_basin, "PER_DATA_BASIN_SUPPORT.csv")

metric_specification <- data.frame(
  metric = c(
    "signal_nrmse_ppi_selected", "dense_signal_nrmse_ppi_selected",
    "loading_relative_l2_error_mean", "feature_total_ise_mean",
    "feature_estimate_to_projection_ise_mean", "kernel_relative_ise_mean",
    "covariance_operator_relative_error_mean",
    "matched_function_recall_mean", "matched_trajectory_abs_cor_mean",
    "matched_score_abs_cor_mean",
    "factor_process_nrmse_mean_missing_as_zero",
    "complete_contribution_nrmse_mean_missing_as_zero"
  ),
  better = c(rep("lower", 7L), rep("higher", 3L), rep("lower", 2L)),
  stringsAsFactors = FALSE
)

metric_rows <- list()
association_rows <- list()
row_index <- 0L
association_index <- 0L
for (data_id in unique(exact_endpoints$data_id)) {
  group <- exact_endpoints[exact_endpoints$data_id == data_id, , drop = FALSE]
  winner <- group[group$is_frozen_winner, , drop = FALSE]
  assert(nrow(winner) == 1L, "The exact set lost a frozen winner.")
  for (metric_index in seq_len(nrow(metric_specification))) {
    metric <- metric_specification$metric[[metric_index]]
    better <- metric_specification$better[[metric_index]]
    values <- as.numeric(group[[metric]])
    finite <- is.finite(values)
    winner_value <- as.numeric(winner[[metric]])
    if (sum(finite) >= 1L && is.finite(winner_value)) {
      ranked_values <- if (better == "lower") values else -values
      ranks <- rank(ranked_values[finite], ties.method = "min")
      finite_ids <- group$fit_id[finite]
      winner_rank <- ranks[match(winner$fit_id, finite_ids)]
      median_value <- stats::median(values[finite])
      winner_better_median <- if (better == "lower") {
        winner_value <= median_value
      } else winner_value >= median_value
    } else {
      winner_rank <- median_value <- NA_real_
      winner_better_median <- NA
    }
    row_index <- row_index + 1L
    metric_rows[[row_index]] <- data.frame(
      cohort = group$cohort[[1L]],
      scenario_id = group$scenario_id[[1L]], data_id = data_id,
      metric = metric, better = better,
      exact_endpoints_with_metric = sum(finite),
      exact_min = safe_min(values), exact_median = safe_median(values),
      exact_max = safe_max(values), winner_fit_id = winner$fit_id[[1L]],
      winner_value = winner_value, winner_rank = winner_rank,
      winner_rank_fraction = if (sum(finite) > 1L && is.finite(winner_rank)) {
        (winner_rank - 1) / (sum(finite) - 1)
      } else if (sum(finite) == 1L && is.finite(winner_rank)) 0 else NA_real_,
      winner_better_or_equal_exact_median = winner_better_median,
      stringsAsFactors = FALSE
    )
    if (sum(finite) >= 3L) {
      rho <- suppressWarnings(stats::cor(
        group$final_elbo[finite], values[finite], method = "spearman"
      ))
      association_index <- association_index + 1L
      association_rows[[association_index]] <- data.frame(
        cohort = group$cohort[[1L]],
        scenario_id = group$scenario_id[[1L]], data_id = data_id,
        metric = metric, better = better,
        exact_endpoints = sum(finite), spearman_elbo_metric = rho,
        direction_aligned = if (is.finite(rho)) {
          if (better == "lower") rho <= 0 else rho >= 0
        } else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}
per_data_metric <- do.call(rbind, metric_rows)
elbo_science <- do.call(rbind, association_rows)
write_csv(per_data_metric, "PER_DATA_EXACT_BASIN_METRIC_AUDIT.csv")
write_csv(elbo_science, "WITHIN_DATA_ELBO_SCIENCE_ASSOCIATION.csv")

association_summary <- do.call(rbind, lapply(
  split(elbo_science, elbo_science$metric),
  function(group) data.frame(
    metric = group$metric[[1L]], better = group$better[[1L]],
    data_sets = nrow(group),
    finite_correlations = sum(is.finite(group$spearman_elbo_metric)),
    median_spearman = safe_median(group$spearman_elbo_metric),
    direction_aligned_data_sets = sum(group$direction_aligned, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
))
winner_rank_summary <- do.call(rbind, lapply(
  split(per_data_metric, per_data_metric$metric),
  function(group) {
    available <- is.finite(group$winner_rank_fraction)
    data.frame(
      metric = group$metric[[1L]], better = group$better[[1L]],
      data_sets = sum(available),
      winner_better_or_equal_median = sum(
        group$winner_better_or_equal_exact_median[available], na.rm = TRUE
      ),
      median_winner_rank_fraction =
        safe_median(group$winner_rank_fraction[available]),
      maximum_winner_rank_fraction =
        safe_max(group$winner_rank_fraction[available]),
      stringsAsFactors = FALSE
    )
  }
))
write_csv(association_summary, "ELBO_SCIENCE_ASSOCIATION_SUMMARY.csv")
write_csv(winner_rank_summary, "WINNER_WITHIN_EXACT_RANK_SUMMARY.csv")

quality_rows <- list()
quality_index <- 0L
quality_groups <- list(
  exact = endpoints$exact_correct_1_1_1,
  nonexact = !endpoints$exact_correct_1_1_1,
  frozen_winner = endpoints$is_frozen_winner
)
for (cohort in c("all_cohorts", unique(endpoints$cohort))) {
  cohort_keep <- if (cohort == "all_cohorts") {
    rep(TRUE, nrow(endpoints))
  } else endpoints$cohort == cohort
  for (group_name in names(quality_groups)) {
    keep <- cohort_keep & quality_groups[[group_name]]
    for (metric_index in seq_len(nrow(metric_specification))) {
      metric <- metric_specification$metric[[metric_index]]
      values <- as.numeric(endpoints[[metric]][keep])
      quality_index <- quality_index + 1L
      quality_rows[[quality_index]] <- data.frame(
        cohort = cohort, endpoint_group = group_name,
        endpoints = sum(keep), metric = metric,
        better = metric_specification$better[[metric_index]],
        finite_values = sum(is.finite(values)), mean = safe_mean(values),
        median = safe_median(values), q25 = safe_quantile(values, 0.25),
        q75 = safe_quantile(values, 0.75), min = safe_min(values),
        max = safe_max(values), stringsAsFactors = FALSE
      )
    }
  }
}
quality_summary <- do.call(rbind, quality_rows)
write_csv(quality_summary, "CORRECT_BASIN_QUALITY_SUMMARY.csv")

conditional_convergence_rows <- list()
conditional_convergence_index <- 0L
for (endpoint_pass in c(TRUE, FALSE)) {
  keep <- endpoints$exact_correct_1_1_1 &
    endpoints$auxiliary_endpoint_pass == endpoint_pass
  for (metric_index in seq_len(nrow(metric_specification))) {
    metric <- metric_specification$metric[[metric_index]]
    values <- as.numeric(endpoints[[metric]][keep])
    conditional_convergence_index <- conditional_convergence_index + 1L
    conditional_convergence_rows[[conditional_convergence_index]] <- data.frame(
      structure_group = "exact", auxiliary_endpoint_pass = endpoint_pass,
      endpoints = sum(keep), metric = metric,
      better = metric_specification$better[[metric_index]],
      finite_values = sum(is.finite(values)), mean = safe_mean(values),
      median = safe_median(values), q25 = safe_quantile(values, 0.25),
      q75 = safe_quantile(values, 0.75), min = safe_min(values),
      max = safe_max(values), stringsAsFactors = FALSE
    )
  }
}
conditional_convergence_summary <-
  do.call(rbind, conditional_convergence_rows)
write_csv(
  conditional_convergence_summary,
  "CONVERGENCE_CONDITIONAL_QUALITY_SUMMARY.csv"
)

process$true_component <- ifelse(
  process$true_role == "shared",
  paste0("A_study", process$study), paste0("B", process$study)
)
process_summary_rows <- list()
process_summary_index <- 0L
for (group_name in c("exact", "frozen_winner")) {
  group_keep <- if (group_name == "exact") {
    process$structure_exact
  } else process$is_frozen_winner
  for (component in c("all", sort(unique(process$true_component)))) {
    keep <- group_keep & if (component == "all") {
      rep(TRUE, nrow(process))
    } else process$true_component == component
    process_error <- process$factor_process_nrmse_missing_as_zero[keep]
    contribution_error <-
      process$complete_contribution_nrmse_missing_as_zero[keep]
    process_summary_index <- process_summary_index + 1L
    process_summary_rows[[process_summary_index]] <- data.frame(
      endpoint_group = group_name, true_component = component,
      rows = sum(keep), matched_rows = sum(process$matched_by_loading[keep]),
      matched_rate = safe_mean(process$matched_by_loading[keep]),
      process_nrmse_mean = safe_mean(process_error),
      process_nrmse_median = safe_median(process_error),
      process_nrmse_min = safe_min(process_error),
      process_nrmse_max = safe_max(process_error),
      contribution_nrmse_mean = safe_mean(contribution_error),
      contribution_nrmse_median = safe_median(contribution_error),
      contribution_nrmse_min = safe_min(contribution_error),
      contribution_nrmse_max = safe_max(contribution_error),
      stringsAsFactors = FALSE
    )
  }
}
process_summary <- do.call(rbind, process_summary_rows)
write_csv(process_summary, "PRIMARY_PROCESS_CONTRIBUTION_SUMMARY.csv")

scale_metric_names <- c(
  "factor_ppi", "factor_scale", "raw_loading_norm_ratio",
  "canonical_loading_norm_ratio", "raw_loading_abs_cosine",
  "canonical_loading_abs_cosine"
)
scale_summary_rows <- list()
scale_summary_index <- 0L
for (group_name in c("exact", "frozen_winner")) {
  group_keep <- if (group_name == "exact") {
    loading_scale$structure_exact
  } else loading_scale$is_frozen_winner
  for (role in sort(unique(loading_scale$block_role))) {
    keep <- group_keep & loading_scale$block_role == role
    for (metric in scale_metric_names) {
      values <- as.numeric(loading_scale[[metric]][keep])
      scale_summary_index <- scale_summary_index + 1L
      scale_summary_rows[[scale_summary_index]] <- data.frame(
        endpoint_group = group_name, block_role = role, metric = metric,
        values = sum(is.finite(values)), mean = safe_mean(values),
        median = safe_median(values), q25 = safe_quantile(values, 0.25),
        q75 = safe_quantile(values, 0.75), min = safe_min(values),
        max = safe_max(values), stringsAsFactors = FALSE
      )
    }
  }
}
scale_summary <- do.call(rbind, scale_summary_rows)
write_csv(scale_summary, "PRIMARY_LOADING_SCALE_SUMMARY.csv")

group_definitions <- list(
  all = rep(TRUE, nrow(endpoints)),
  exact = endpoints$exact_correct_1_1_1,
  nonexact = !endpoints$exact_correct_1_1_1,
  frozen_winner = endpoints$is_frozen_winner
)
convergence_rows <- list()
counter <- 0L
for (cohort in c("all_cohorts", unique(endpoints$cohort))) {
  cohort_keep <- if (cohort == "all_cohorts") {
    rep(TRUE, nrow(endpoints))
  } else endpoints$cohort == cohort
  for (group_name in names(group_definitions)) {
    keep <- cohort_keep & group_definitions[[group_name]]
    counter <- counter + 1L
    convergence_rows[[counter]] <- data.frame(
      cohort = cohort, endpoint_group = group_name,
      endpoints = sum(keep),
      objective_eligible = sum(endpoints$objective_eligible[keep]),
      maximum_200_sweeps = sum(endpoints$actual_T1_sweeps[keep] == 200L),
      strict_practical_converged =
        sum(endpoints$strict_practical_converged[keep]),
      auxiliary_endpoint_pass =
        sum(endpoints$auxiliary_endpoint_pass[keep]),
      auxiliary_ever_stable_5 =
        sum(endpoints$auxiliary_ever_stable_5[keep]),
      auxiliary_endpoint_stable_5 =
        sum(endpoints$auxiliary_endpoint_stable_5[keep]),
      stringsAsFactors = FALSE
    )
  }
}
convergence_summary <- do.call(rbind, convergence_rows)
write_csv(convergence_summary, "CONVERGENCE_STRUCTURE_SUMMARY.csv")

audit_summary <- data.frame(
  item = c(
    "standard_density_data_sets", "standard_density_g12_endpoints",
    "primary_endpoints_with_full_rds", "secondary_aggregate_endpoints",
    "exact_identity_endpoints", "frozen_truth_free_winners",
    "exact_frozen_winners", "primary_process_role_rows",
    "primary_loading_scale_rows", "quality_summary_rows",
    "convergence_conditional_quality_rows", "process_summary_rows",
    "loading_scale_summary_rows",
    "new_data_generated", "new_fit_run", "continuation_run",
    "winner_selection_reopened", "formal_v0lv_result"
  ),
  value = c(
    length(unique(endpoints$data_id)), nrow(endpoints), nrow(independent),
    nrow(cross), sum(endpoints$exact_correct_1_1_1),
    sum(endpoints$is_frozen_winner),
    sum(endpoints$is_frozen_winner & endpoints$exact_correct_1_1_1),
    nrow(process), nrow(loading_scale), nrow(quality_summary),
    nrow(conditional_convergence_summary), nrow(process_summary),
    nrow(scale_summary), "FALSE", "FALSE", "FALSE", "FALSE", "FALSE"
  ),
  stringsAsFactors = FALSE
)
write_csv(audit_summary, "AUDIT_SUMMARY.csv")

qc <- data.frame(
  check_id = c(
    "ten_standard_density_data_sets", "one_hundred_four_g12_endpoints",
    "primary_72_complete", "secondary_32_complete",
    "forty_five_exact_endpoints", "ten_frozen_truth_free_winners",
    "all_winners_exact", "all_objectives_eligible",
    "historical_truth_isolation_intact", "primary_process_rows_complete",
    "primary_scale_rows_complete", "summary_tables_complete",
    "no_new_fit_or_continuation", "selection_not_reopened",
    "development_not_formal"
  ),
  passed = c(
    length(unique(endpoints$data_id)) == 10L,
    nrow(endpoints) == 104L,
    nrow(independent) == 72L,
    nrow(cross) == 32L,
    sum(endpoints$exact_correct_1_1_1) == 45L,
    sum(endpoints$is_frozen_winner) == 10L,
    all(endpoints$exact_correct_1_1_1[endpoints$is_frozen_winner]),
    all(endpoints$objective_eligible),
    !any(endpoints$truth_used_for_fit_or_selection),
    nrow(process) == 72L * 4L && nrow(process_endpoint) == 72L,
    nrow(loading_scale) == 72L * 3L,
    nrow(quality_summary) == 108L &&
      nrow(conditional_convergence_summary) == 24L &&
      nrow(process_summary) == 10L && nrow(scale_summary) == 36L,
    TRUE, TRUE, !any(endpoints$formal_v0lv_result)
  ),
  stringsAsFactors = FALSE
)
write_csv(qc, "AUDIT_QC.csv")
assert(all(qc$passed), "Stage-5A read-only audit QC failed.")

writeLines(c(
  "status=PASS", paste0("completed_utc=", iso_time()),
  "standard_density_data_sets=10", "g12_endpoints=104",
  "exact_identity_endpoints=45", "frozen_truth_free_winners=10",
  "primary_process_role_rows=288", "primary_loading_scale_rows=216",
  "quality_summary_rows=108", "process_summary_rows=10",
  "convergence_conditional_quality_rows=24",
  "loading_scale_summary_rows=36",
  "new_data_generated=FALSE", "new_fit_run=FALSE",
  "continuation_run=FALSE", "winner_selection_reopened=FALSE",
  "formal_v0lv_result=FALSE"
), file.path(output, "AUDIT_COMPLETE.txt"), useBytes = TRUE)
utils::capture.output(sessionInfo(), file = file.path(output, "SESSION_INFO.txt"))

cat(paste0(
  "STAGE5A_READONLY_AUDIT_PASS endpoints=", nrow(endpoints),
  " exact=", sum(endpoints$exact_correct_1_1_1),
  " winners=", sum(endpoints$is_frozen_winner), "\n"
))
