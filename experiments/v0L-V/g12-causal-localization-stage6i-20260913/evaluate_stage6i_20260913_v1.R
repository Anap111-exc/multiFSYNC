#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
experiments_root <- dirname(script_root)
source(file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831",
  "stage6a_common_20260831_v1.R"
), local = FALSE)
source(file.path(
  experiments_root, "g12-reachability-causal-audit-stage6h-20260913",
  "stage6h_evaluation_helpers_20260913_v1.R"
), local = FALSE)
source(file.path(script_root, "stage6i_helpers_20260913_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <- "G12_STAGE6I_CAUSAL_LOCALIZATION_V1_20260913"
cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- cli$`output-root` %||% pp_env("STAGE6I_OUTPUT_ROOT")
workers <- as.integer(cli$workers %||% pp_env("STAGE6I_EVAL_WORKERS", "7"))
stage6g_output_root <- pp_env(
  "STAGE6I_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6h_output_root <- pp_env(
  "STAGE6I_STAGE6H_OUTPUT_ROOT",
  "/root/multifsync-results/stage6h_20260913_v1"
)
package_root <- pp_env("STAGE6I_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
stage6i_library <- pp_env("STAGE6I_R_LIB")
if (nzchar(stage6i_library)) Sys.setenv(STAGE6A_R_LIB = stage6i_library)
pp_assert(nzchar(output_root) && dir.exists(output_root) &&
            .Platform$OS.type == "unix" && workers %in% 1:7,
          "Stage 6I evaluation requires a prepared Unix output and 1--7 workers.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

authorization <- readLines(file.path(
  script_root, "STAGE6I_EXECUTION_AUTHORIZATION_20260913.txt"
), warn = FALSE)
pp_assert(all(c(
  "status=AUTHORIZED_POST_UNSEAL_CAUSAL_LOCALIZATION",
  "scientific_evaluation_authorized=TRUE",
  "truth_available_to_oracle_initialization=TRUE",
  "cross_dimension_selection_authorized=FALSE",
  "new_dimension_criterion_authorized=FALSE"
) %in% authorization), "Stage 6I evaluation authorization is incomplete.")
pp_assert(file.exists(file.path(output_root, "STAGE6I_FITS_COMPLETE.txt")),
          "All Stage 6I fits must finish before evaluation.")

manifest <- pp_read_csv(file.path(output_root, "FIT_MANIFEST.csv"))
terminals <- pp_read_csv(file.path(output_root, "ALL_NEW_TERMINALS.csv"))
pp_assert(nrow(manifest) == 10L && nrow(terminals) == 10L &&
            setequal(manifest$fit_id, terminals$fit_id) &&
            all(terminals$objective_eligible),
          "Stage 6I endpoint accounting failed.")

source_binding <- pp_read_csv(file.path(output_root,
                                        "SOURCE_ARTIFACT_BINDING.csv"))
source_info <- file.info(source_binding$path)
source_match <- file.exists(source_binding$path) &
  as.numeric(source_info$size) == as.numeric(source_binding$bytes) &
  vapply(source_binding$path, pp_sha256, character(1L)) == source_binding$sha256
pp_assert(all(source_match), "A read-only Stage 6I source artifact changed.")

truth_path <- source_binding$path[source_binding$role == "truth"]
truth <- readRDS(truth_path)
swapped_truth <- readRDS(file.path(output_root, "inputs", "swapped_truth.rds"))
observation <- readRDS(source_binding$path[source_binding$role == "observation"])
pp_assert(identical(truth$data_id, "s6g_minimal_L111_M2") &&
            identical(swapped_truth$study_permutation, c(2L, 1L)),
          "Stage 6I truth identity failed.")

pp_load_package()
pp_source_runtime(TRUE)
evaluator_environment <- candidate2_load_frozen_multi_evaluator(PP_SNAPSHOT)
evaluation_root <- file.path(output_root, "scientific_evaluation")
pp_assert(!dir.exists(evaluation_root),
          "Stage 6I scientific evaluation already exists.")
dir.create(evaluation_root, recursive = TRUE, showWarnings = FALSE,
           mode = "0700")

s6i_add_recovery_flags <- function(value) {
  value$ppi_factor_counts_exact <- with(
    value, selected_shared == truth_L_f &
      selected_specific_study_1 == truth_L_s_1 &
      selected_specific_study_2 == truth_L_s_2
  )
  value$core_truth_recovered <- with(
    value, global_correct == truth_global_direction_count &
      global_misplaced == 0L & global_missing == 0L &
      study_correct == truth_study_role_count &
      study_misplaced == 0L & study_missing == 0L
  )
  value$exact_thresholded_factor_structure <- with(
    value, core_truth_recovered & ppi_factor_counts_exact &
      global_extra == 0L & global_duplicate == 0L &
      study_candidate_extra == 0L & study_candidate_duplicate == 0L
  )
  value$full_true_feature_component_coverage <- with(
    value, feature_components_matched == feature_components_total
  )
  value$core_with_full_feature_coverage <- with(
    value, core_truth_recovered & full_true_feature_component_coverage
  )
  value
}

evaluate_one <- function(index) {
  row <- manifest[index, , drop = FALSE]
  fit_id <- row$fit_id[[1L]]
  fit_path <- file.path(output_root, "fits", fit_id, "fit.rds")
  result_dir <- file.path(evaluation_root, fit_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  truth_current <- if (row$study_order[[1L]] == "swapped") {
    swapped_truth
  } else truth
  warning_rows <- list()
  started <- Sys.time()
  result <- tryCatch(withCallingHandlers({
    s6h_evaluate_fit(readRDS(fit_path), truth_current,
                     evaluator_environment)
  }, warning = function(condition) {
    warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
      class = paste(class(condition), collapse = ";"),
      message = conditionMessage(condition), stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  }), error = function(condition) condition)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(class = character(), message = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(result_dir, "EVALUATION_WARNINGS.csv"))
  if (inherits(result, "condition")) {
    pp_write_lines(c("status=ERROR", paste0("message=", conditionMessage(result))),
                   file.path(result_dir, "EVALUATION_ERROR.txt"))
    return(list(status = data.frame(
      fit_id = fit_id, ok = FALSE, elapsed_seconds = elapsed,
      error = conditionMessage(result), stringsAsFactors = FALSE
    ), compact = NULL))
  }
  pp_save_rds(result, file.path(result_dir, "evaluation.rds"))
  for (name in names(result)) if (is.data.frame(result[[name]])) {
    pp_write_csv(result[[name]], file.path(result_dir, paste0(name, ".csv")))
  }
  terminal <- terminals[terminals$fit_id == fit_id, , drop = FALSE]
  base <- data.frame(
    fit_id = fit_id, fit_kind = row$fit_kind[[1L]],
    study_order = row$study_order[[1L]],
    oracle_mode = row$oracle_mode[[1L]], replicate = row$replicate[[1L]],
    seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
    truth_used_for_fit = row$truth_available_to_fit[[1L]],
    final_elbo = terminal$final_elbo[[1L]],
    practical_converged = terminal$practical_converged[[1L]],
    convergence_status = terminal$convergence_status[[1L]],
    fit_elapsed_seconds = terminal$elapsed_seconds[[1L]],
    fit_warning_count = terminal$warning_count[[1L]],
    generating_dimension_configuration = TRUE,
    truth_L_f = 1L, truth_L_s_1 = 1L, truth_L_s_2 = 1L,
    truth_M_common = 2L, stringsAsFactors = FALSE
  )
  activity <- result$candidate_activity_summary
  names(activity) <- paste0("activity_", names(activity))
  compact <- cbind(base, result$summary,
                   result$process_contribution_summary,
                   result$loading_scale_summary, activity)
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("truth_used_for_fit=", row$truth_available_to_fit[[1L]]),
    "truth_used_for_stopping_or_selection=FALSE",
    "cross_dimension_selection=FALSE"
  ), file.path(result_dir, "EVALUATION_COMPLETE.txt"))
  list(status = data.frame(
    fit_id = fit_id, ok = TRUE, elapsed_seconds = elapsed, error = "",
    stringsAsFactors = FALSE
  ), compact = compact)
}

results <- parallel::mclapply(
  seq_len(nrow(manifest)), evaluate_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
status <- do.call(rbind, lapply(results, `[[`, "status"))
pp_write_csv(status, file.path(evaluation_root, "EVALUATION_STATUS.csv"))
pp_assert(all(status$ok), "Stage 6I evaluation failed; no selective rerun.")
science <- pp_rbind_fill(lapply(results, `[[`, "compact"))
science <- science[match(manifest$fit_id, science$fit_id), , drop = FALSE]
science <- s6i_add_recovery_flags(science)
pp_write_csv(science, file.path(evaluation_root,
                                "ALL_STAGE6I_NEW_SCIENTIFIC_RESULTS.csv"))

# Import only the four paired original-order Stage 6H full-data endpoints.
stage6h_science <- pp_read_csv(file.path(
  stage6h_output_root, "truth_evaluation_20260913_v1",
  "ALL_STAGE6H_SCIENTIFIC_RESULTS.csv"
))
original <- stage6h_science[
  stage6h_science$condition_id == "full_data_6_9" &
    stage6h_science$seed_index %in% 1:4, , drop = FALSE
]
pp_assert(nrow(original) == 4L, "Paired Stage 6H science import failed.")
original$fit_kind <- "original_g12_read_only"
original$study_order <- "original"
original$oracle_mode <- "none"
original$truth_used_for_fit <- FALSE
original$replicate <- original$seed_index
comparison_science <- pp_rbind_fill(list(original, science))
pp_write_csv(comparison_science, file.path(
  evaluation_root, "ALL_STAGE6I_COMPARISON_SCIENTIFIC_RESULTS.csv"
))

# ----- Continuous versus observed Gram and true signal geometry -----
reference_fit_path <- file.path(
  stage6h_output_root, "fits",
  "s6h_minimal_L111_M2__full_data_6_9__G12_01", "fit.rds"
)
reference_fit <- readRDS(reference_fit_path)
weights <- s6i_trapezoid_weights(reference_fit$time_g)
G_cont <- crossprod(reference_fit$C_g, weights * reference_fit$C_g)
obs_counts <- vapply(observation$time_obs, function(study) {
  sum(vapply(study, length, integer(1L)))
}, integer(1L))
G_obs <- lapply(seq_len(2L), function(study) {
  reference_fit$sum_list_cp_C[[study]] / obs_counts[[study]]
})
G_obs_pooled <- (reference_fit$sum_list_cp_C[[1L]] +
                   reference_fit$sum_list_cp_C[[2L]]) / sum(obs_counts)
gram_names <- c("study1", "study2", "pooled")
gram_values <- c(G_obs, list(G_obs_pooled))
gram_summary <- do.call(rbind, lapply(seq_along(gram_values), function(index) {
  current <- gram_values[[index]]
  eigenvalues <- s6i_generalized_eigenvalues(G_cont, current)
  data.frame(
    gram_id = gram_names[[index]], observations = if (index <= 2L) {
      obs_counts[[index]]
    } else sum(obs_counts),
    rank = qr(current)$rank,
    condition_number = kappa(current, exact = TRUE),
    relative_frobenius_difference =
      sqrt(sum((current - G_cont)^2)) / sqrt(sum(G_cont^2)),
    generalized_eigen_min = min(eigenvalues),
    generalized_eigen_median = stats::median(eigenvalues),
    generalized_eigen_max = max(eigenvalues),
    stringsAsFactors = FALSE
  )
}))
pp_write_csv(gram_summary, file.path(evaluation_root,
                                     "GRAM_DESIGN_SUMMARY.csv"))
eigen_table <- do.call(rbind, lapply(seq_along(gram_values), function(index) {
  values <- s6i_generalized_eigenvalues(G_cont, gram_values[[index]])
  data.frame(gram_id = gram_names[[index]], eigen_index = seq_along(values),
             generalized_eigenvalue = values, stringsAsFactors = FALSE)
}))
pp_write_csv(eigen_table, file.path(evaluation_root,
                                    "GRAM_GENERALIZED_EIGENVALUES.csv"))

true_data <- truth$data
true_parameters <- true_data$true_params
true_feature_rows <- list()
add_true_features <- function(block, study, functions, observed_times) {
  for (component in seq_len(ncol(functions))) {
    dense <- functions[, component]
    observed <- stats::approx(
      true_parameters$t_grid_dense, dense, xout = observed_times, rule = 2
    )$y
    projection <- s6i_project_dense_block(
      matrix(dense, ncol = 1L), true_parameters$t_grid_dense,
      reference_fit$time_g, reference_fit$C_g
    )$coefficients[, 1L]
    total_energy <- as.numeric(crossprod(projection, G_cont %*% projection))
    null <- projection
    null[-c(1L, 2L)] <- 0
    null_energy <- as.numeric(crossprod(null, G_cont %*% null))
    true_feature_rows[[length(true_feature_rows) + 1L]] <<- data.frame(
      block = block, study = study, component = component,
      continuous_energy = sum(weights * stats::approx(
        true_parameters$t_grid_dense, dense,
        xout = reference_fit$time_g, rule = 2
      )$y^2),
      observed_energy = mean(observed^2),
      leverage_ratio = sum(weights * stats::approx(
        true_parameters$t_grid_dense, dense,
        xout = reference_fit$time_g, rule = 2
      )$y^2) / max(mean(observed^2), .Machine$double.eps),
      projection_relative_rmse = s6i_project_dense_block(
        matrix(dense, ncol = 1L), true_parameters$t_grid_dense,
        reference_fit$time_g, reference_fit$C_g
      )$relative_rmse,
      projected_nullspace_energy_fraction = null_energy /
        max(total_energy, .Machine$double.eps),
      stringsAsFactors = FALSE
    )
  }
}
add_true_features(
  "shared", NA_integer_, true_parameters$theta_dense_list[[1L]],
  unlist(observation$time_obs, recursive = TRUE, use.names = FALSE)
)
for (study in 1:2) add_true_features(
  "specific", study, true_parameters$kappa_dense_list[[study]][[1L]],
  unlist(observation$time_obs[[study]], use.names = FALSE)
)
true_features <- do.call(rbind, true_feature_rows)
pp_write_csv(true_features, file.path(evaluation_root,
                                      "TRUE_FEATURE_VISIBILITY.csv"))

true_block_rows <- list()
for (study in 1:2) {
  for (role in c("shared", "specific")) {
    loading <- if (role == "shared") true_parameters$a_true[, 1L] else
      true_parameters$b_true[[study]][, 1L]
    process_observed <- unlist(lapply(
      if (role == "shared") true_parameters$f_true_values[[study]] else
        true_parameters$g_true_values[[study]], function(value) value[, 1L]
    ), use.names = FALSE)
    scores <- if (role == "shared") {
      as.matrix(true_parameters$eta_true[[study]][[1L]])
    } else as.matrix(true_parameters$chi_true[[study]][[1L]])
    functions <- if (role == "shared") {
      true_parameters$theta_dense_list[[1L]]
    } else true_parameters$kappa_dense_list[[study]][[1L]]
    process_dense <- scores %*% t(functions)
    loading_ms <- mean(loading^2)
    contribution_rms <- sqrt(mean(process_observed^2) * loading_ms)
    true_block_rows[[length(true_block_rows) + 1L]] <- data.frame(
      study = study, role = role,
      active_loading_count = sum(abs(loading) > 0),
      loading_l2 = sqrt(sum(loading^2)), loading_rms = sqrt(loading_ms),
      process_observed_rms = sqrt(mean(process_observed^2)),
      process_dense_rms = sqrt(mean(process_dense^2)),
      process_continuous_to_observed_energy =
        mean(process_dense^2) / max(mean(process_observed^2),
                                    .Machine$double.eps),
      complete_contribution_observed_rms = contribution_rms,
      contribution_to_noise_variance_ratio =
        contribution_rms^2 / true_parameters$sigma2_eps_true,
      stringsAsFactors = FALSE
    )
  }
}
true_blocks <- do.call(rbind, true_block_rows)
pp_write_csv(true_blocks, file.path(evaluation_root,
                                    "TRUE_BLOCK_SIGNAL_GEOMETRY.csv"))

# ----- Fitted continuous/observed energy, null-space share, and early PPI -----
fit_registry <- rbind(
  data.frame(
    fit_id = original$fit_id, fit_kind = "original_g12_read_only",
    study_order = "original", fit_path = file.path(
      stage6h_output_root, "fits", original$fit_id, "fit.rds"
    ), stringsAsFactors = FALSE
  ),
  data.frame(
    fit_id = manifest$fit_id, fit_kind = manifest$fit_kind,
    study_order = manifest$study_order, fit_path = file.path(
      output_root, "fits", manifest$fit_id, "fit.rds"
    ), stringsAsFactors = FALSE
  )
)
fitted_rows <- list()
trace_rows <- list()
for (registry_index in seq_len(nrow(fit_registry))) {
  registry <- fit_registry[registry_index, , drop = FALSE]
  fit <- readRDS(registry$fit_path[[1L]])
  fit_weights <- s6i_trapezoid_weights(fit$time_g)
  fit_G_cont <- crossprod(fit$C_g, fit_weights * fit$C_g)
  fit_obs_counts <- vapply(fit$list_cp_C, length, integer(1L))
  # Subject counts above are not time counts; recover time-count normalization
  # from the diagonal scale by using the Stage 6G observation order.
  fit_time_obs <- if (registry$study_order[[1L]] == "swapped") {
    rev(observation$time_obs)
  } else observation$time_obs
  fit_obs_counts <- vapply(fit_time_obs, function(value) {
    sum(vapply(value, length, integer(1L)))
  }, integer(1L))
  fit_G_obs <- lapply(1:2, function(study) {
    fit$sum_list_cp_C[[study]] / fit_obs_counts[[study]]
  })
  fit_G_pool <- (fit$sum_list_cp_C[[1L]] + fit$sum_list_cp_C[[2L]]) /
    sum(fit_obs_counts)
  block_specs <- list(
    list(block = "shared", study = NA_integer_,
         mu = fit$mu_q_nu_phi[[1L]],
         Sigma = fit$Sigma_q_nu_phi[[1L]], gram = fit_G_pool,
         loading = fit$mu_q_a_original[, 1L],
         ppi = fit$mu_q_gamma_a[, 1L]),
    list(block = "specific", study = 1L,
         mu = fit$mu_q_nu_psi[[1L]][[1L]],
         Sigma = fit$Sigma_q_nu_psi[[1L]][[1L]], gram = fit_G_obs[[1L]],
         loading = fit$mu_q_b_specific_original[[1L]][, 1L],
         ppi = fit$mu_q_gamma_b[[1L]][, 1L]),
    list(block = "specific", study = 2L,
         mu = fit$mu_q_nu_psi[[2L]][[1L]],
         Sigma = fit$Sigma_q_nu_psi[[2L]][[1L]], gram = fit_G_obs[[2L]],
         loading = fit$mu_q_b_specific_original[[2L]][, 1L],
         ppi = fit$mu_q_gamma_b[[2L]][, 1L])
  )
  for (spec in block_specs) {
    for (component in seq_len(ncol(spec$mu))) {
      mu <- spec$mu[, component]
      Sigma <- spec$Sigma[[component]]
      E_cont <- s6i_second_moment_energy(mu, Sigma, fit_G_cont)
      E_obs <- s6i_second_moment_energy(mu, Sigma, spec$gram)
      null_mu <- mu
      null_mu[-c(1L, 2L)] <- 0
      null_Sigma <- matrix(0, nrow(Sigma), ncol(Sigma))
      null_Sigma[1:2, 1:2] <- Sigma[1:2, 1:2]
      E_null <- s6i_second_moment_energy(null_mu, null_Sigma, fit_G_cont)
      fitted_rows[[length(fitted_rows) + 1L]] <- data.frame(
        fit_id = registry$fit_id[[1L]], fit_kind = registry$fit_kind[[1L]],
        study_order = registry$study_order[[1L]], block = spec$block,
        study = spec$study, component = component,
        continuous_second_moment_energy = E_cont,
        observed_second_moment_energy = E_obs,
        leverage_ratio = E_cont / max(E_obs, .Machine$double.eps),
        nullspace_continuous_energy_fraction = E_null /
          max(E_cont, .Machine$double.eps),
        loading_l2 = sqrt(sum(spec$loading^2)),
        loading_ppi_mean = mean(spec$ppi), loading_ppi_max = max(spec$ppi),
        factor_ppi = 1 - prod(1 - spec$ppi), stringsAsFactors = FALSE
      )
    }
  }
  trace <- fit$driver_trace
  if (is.data.frame(trace) && nrow(trace)) {
    keep <- trace$block %in% c("shared_loading", "specific_loading") &
      trace$iteration %in% c(1L, 380L, 400L)
    if (any(keep)) {
      current <- trace[keep, c(
        "iteration", "temperature", "phase", "block", "study", "factor",
        "posterior_mean_norm", "ppi_mean", "ppi_max", "gate_active"
      ), drop = FALSE]
      current$fit_id <- registry$fit_id[[1L]]
      current$fit_kind <- registry$fit_kind[[1L]]
      current$study_order <- registry$study_order[[1L]]
      trace_rows[[length(trace_rows) + 1L]] <- current
    }
  }
}
fitted_geometry <- do.call(rbind, fitted_rows)
pp_write_csv(fitted_geometry, file.path(evaluation_root,
                                        "FITTED_FUNCTION_GEOMETRY.csv"))
early_trace <- if (length(trace_rows)) do.call(rbind, trace_rows) else
  data.frame()
pp_write_csv(early_trace, file.path(evaluation_root,
                                    "G12_LOADING_TRACE_CHECKPOINTS.csv"))

scale_by_fit <- aggregate(
  cbind(loading_l2, continuous_second_moment_energy,
        observed_second_moment_energy, leverage_ratio,
        nullspace_continuous_energy_fraction) ~ fit_id + fit_kind +
    study_order + block + study,
  data = fitted_geometry, FUN = mean, na.action = na.pass
)
scale_by_fit$complete_scale_proxy <- scale_by_fit$loading_l2 *
  sqrt(scale_by_fit$continuous_second_moment_energy)
pp_write_csv(scale_by_fit, file.path(evaluation_root,
                                     "SCALE_AND_NULLSPACE_BY_FIT.csv"))

# Static source audit: identical study-loop branches do not exclude order effects.
audit_specs <- data.frame(
  file = c("R/update_xi.R", "R/update_eigenfunctions.R",
           "R/update_loadings.R", "R/multi_core.R",
           "R/pre_score_interface.R"),
  target = c("specific scores", "specific time functions",
             "specific loadings", "initialization and update schedule",
             "pre-score public interface"),
  stringsAsFactors = FALSE
)
code_audit <- do.call(rbind, lapply(seq_len(nrow(audit_specs)), function(index) {
  path <- file.path(package_root, audit_specs$file[[index]])
  lines <- readLines(path, warn = FALSE)
  data.frame(
    file = audit_specs$file[[index]], target = audit_specs$target[[index]],
    sha256 = pp_sha256(path),
    all_study_loop_hits = sum(grepl(
      "for \\(s in 1:S\\)|for \\(s in seq_len\\(S\\)\\)", lines
    )),
    literal_study1_index_hits = sum(grepl("\\[\\[1L?\\]\\]", lines)),
    formula_branch_is_study_generic = TRUE,
    sequential_order_or_rng_can_still_be_position_dependent = TRUE,
    stringsAsFactors = FALSE
  )
}))
pp_write_csv(code_audit, file.path(evaluation_root,
                                   "CODE_STUDY_SYMMETRY_AUDIT.csv"))

# ----- Causal summaries -----
recovery_summary <- do.call(rbind, lapply(
  split(science, interaction(science$fit_kind, science$study_order,
                             science$oracle_mode, drop = TRUE)),
  function(rows) data.frame(
    fit_kind = rows$fit_kind[[1L]], study_order = rows$study_order[[1L]],
    oracle_mode = rows$oracle_mode[[1L]], endpoints = nrow(rows),
    core_hits = sum(rows$core_truth_recovered),
    exact_hits = sum(rows$exact_thresholded_factor_structure),
    full_feature_hits = sum(rows$full_true_feature_component_coverage),
    practical_converged = sum(rows$practical_converged),
    mean_observed_nrmse = mean(rows$signal_nrmse_ppi_selected),
    mean_dense_nrmse = mean(rows$dense_signal_nrmse_ppi_selected),
    mean_covariance_error = mean(rows$covariance_operator_relative_error_mean),
    mean_loading_error = mean(rows$loading_relative_l2_error_mean),
    stringsAsFactors = FALSE
  )
))
pp_write_csv(recovery_summary, file.path(evaluation_root,
                                        "RECOVERY_SUMMARY.csv"))

swap_rows <- science[science$fit_kind == "study_swap_g12", , drop = FALSE]
original_rows <- original[match(swap_rows$seed_index, original$seed_index),
                          , drop = FALSE]
swap_pairs <- data.frame(
  seed_index = swap_rows$seed_index, fit_seed = swap_rows$fit_seed,
  original_order_B10_recovery = original_rows$R_B1,
  original_order_B20_recovery = original_rows$R_B2,
  swapped_position1_original_B20_recovery = swap_rows$R_B1,
  swapped_position2_original_B10_recovery = swap_rows$R_B2,
  original_core = original_rows$core_truth_recovered,
  swapped_core = swap_rows$core_truth_recovered,
  original_elbo = original_rows$final_elbo,
  swapped_elbo = swap_rows$final_elbo,
  stringsAsFactors = FALSE
)
pp_write_csv(swap_pairs, file.path(evaluation_root,
                                   "STUDY_SWAP_PAIRED_RESULTS.csv"))

oracle <- science[science$fit_kind == "truth_neighborhood", , drop = FALSE]
oracle_original_hits <- sum(oracle$core_truth_recovered &
                              oracle$study_order == "original")
oracle_swapped_hits <- sum(oracle$core_truth_recovered &
                             oracle$study_order == "swapped")
oracle_decision <- if (oracle_original_hits >= 2L && oracle_swapped_hits >= 2L) {
  "CORRECT_BASIN_EXISTS_BUT_G12_REACHABILITY_FAILS"
} else if (oracle_original_hits == 0L && oracle_swapped_hits == 0L) {
  "ORACLE_NEIGHBORHOOD_COLLAPSES_CHECK_OBJECTIVE_PRIOR_SCALE_OR_IMPLEMENTATION"
} else {
  "ORACLE_SUCCESS_DEPENDS_ON_STUDY_ORDER_CHECK_POSITION_ASYMMETRY"
}
position1_hits <- sum(swap_rows$R_B1 >= 0.8, na.rm = TRUE)
position2_hits <- sum(swap_rows$R_B2 >= 0.8, na.rm = TRUE)
physical_B10_hits_after_swap <- position2_hits
physical_B20_hits_after_swap <- position1_hits
swap_decision <- if (position1_hits == 0L && position2_hits >= 2L) {
  "FAILURE_FOLLOWS_FIRST_STUDY_POSITION"
} else if (physical_B10_hits_after_swap == 0L &&
           physical_B20_hits_after_swap >= 2L) {
  "FAILURE_FOLLOWS_ORIGINAL_B10_DATA_BLOCK"
} else {
  "STUDY_SWAP_INCONCLUSIVE_OR_GENERAL_COMPETITION"
}

random_original_elbo <- original$final_elbo
core_oracle_original <- oracle[
  oracle$study_order == "original" & oracle$core_truth_recovered, , drop = FALSE
]
elbo_comparison <- data.frame(
  random_original_endpoints = length(random_original_elbo),
  random_original_max_elbo = max(random_original_elbo),
  oracle_original_endpoints = sum(oracle$study_order == "original"),
  oracle_original_core_endpoints = nrow(core_oracle_original),
  oracle_original_core_max_elbo = if (nrow(core_oracle_original)) {
    max(core_oracle_original$final_elbo)
  } else NA_real_,
  core_oracle_minus_random_max_elbo = if (nrow(core_oracle_original)) {
    max(core_oracle_original$final_elbo) - max(random_original_elbo)
  } else NA_real_, stringsAsFactors = FALSE
)
pp_write_csv(elbo_comparison, file.path(evaluation_root,
                                        "ELBO_BASIN_COMPARISON.csv"))

decision <- data.frame(
  oracle_original_core_hits = oracle_original_hits,
  oracle_swapped_core_hits = oracle_swapped_hits,
  swapped_position1_truth_hits = position1_hits,
  swapped_position2_truth_hits = position2_hits,
  original_B10_hits_after_swap = physical_B10_hits_after_swap,
  original_B20_hits_after_swap = physical_B20_hits_after_swap,
  oracle_decision = oracle_decision,
  study_swap_decision = swap_decision,
  continuous_observed_gram_hypothesis_supported =
    max(gram_summary$generalized_eigen_max) >= 2 ||
      max(true_features$leverage_ratio) >= 2,
  nullspace_hypothesis_supported =
    stats::median(fitted_geometry$nullspace_continuous_energy_fraction,
                  na.rm = TRUE) >= 0.5,
  new_dimension_criterion_caused_failure = FALSE,
  stringsAsFactors = FALSE
)
pp_write_csv(decision, file.path(evaluation_root,
                                 "STAGE6I_CAUSAL_DECISION.csv"))

pp_write_lines(c(
  "status=COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "new_endpoints_evaluated=10", "read_only_stage6h_endpoints=4",
  paste0("oracle_decision=", oracle_decision),
  paste0("study_swap_decision=", swap_decision),
  "new_dimension_criterion_used=FALSE",
  "new_dimension_criterion_caused_failure=FALSE",
  "model_prior_cavi_elbo_changed=FALSE", "formal_paper_mc=FALSE"
), file.path(evaluation_root, "EVALUATION_COMPLETE.txt"))
pp_write_lines(c(
  "status=STAGE6I_WORKFLOW_COMPLETE", paste0("completed_utc=", pp_iso_time()),
  "fits=10", "evaluations=10", "geometry_audit=TRUE",
  "formal_paper_mc=FALSE"
), file.path(output_root, "STAGE6I_EVALUATION_COMPLETE.txt"))
cat("STAGE6I_EVALUATION_COMPLETE oracle=", oracle_decision,
    " swap=", swap_decision, "\n", sep = "")
