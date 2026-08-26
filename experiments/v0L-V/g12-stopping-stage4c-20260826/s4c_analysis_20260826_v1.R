#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

s4c_verify_environment(load_runtime = TRUE, package_role = "development")
manifest <- uc_read_csv(file.path(UC_ROOT, "FIT_MANIFEST.csv"))
s4c_validate_fit_manifest(manifest)
uc_assert(file.exists(file.path(UC_ROOT, "FITS_COMPLETE.txt")),
          "Stage-4C fixed-400 fits are not complete.")

analysis_dir <- file.path(UC_ROOT, "analysis_20260826_v1")
uc_assert(!dir.exists(analysis_dir),
          "Analysis directory already exists; refusing to overwrite.")
dir.create(analysis_dir, recursive = TRUE, mode = "0700")
write_csv <- function(value, name) {
  utils::write.csv(value, file.path(analysis_dir, name),
                   row.names = FALSE, quote = TRUE, na = "")
}

nrmse <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  sqrt(mean((current - reference)^2)) /
    (1 + sqrt(mean(reference^2)))
}

relative_distance <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  sqrt(sum((current - reference)^2)) /
    (1 + sqrt(sum(reference^2)))
}

quantile_abs <- function(current, reference, probability = 0.95) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  if (!length(current)) return(0)
  as.numeric(stats::quantile(
    abs(current - reference), probability, names = FALSE, type = 8
  ))
}

max_abs <- function(current, reference) {
  current <- as.numeric(current); reference <- as.numeric(reference)
  if (length(current) != length(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  if (!length(current)) 0 else max(abs(current - reference))
}

subspace_distance <- function(current, reference) {
  current <- as.matrix(current); reference <- as.matrix(reference)
  if (nrow(current) != nrow(reference) ||
      any(!is.finite(c(current, reference)))) return(Inf)
  basis <- function(value) {
    if (!nrow(value) || !ncol(value)) {
      return(matrix(numeric(), nrow = nrow(value), ncol = 0L))
    }
    decomposition <- svd(value, nu = min(dim(value)), nv = 0L)
    tolerance <- sqrt(.Machine$double.eps) * max(dim(value)) *
      max(c(decomposition$d, 1))
    rank <- sum(decomposition$d > tolerance)
    if (!rank) matrix(numeric(), nrow = nrow(value), ncol = 0L) else
      decomposition$u[, seq_len(rank), drop = FALSE]
  }
  left <- basis(current); right <- basis(reference)
  rank_left <- ncol(left); rank_right <- ncol(right)
  max_rank <- max(rank_left, rank_right)
  if (!max_rank) return(0)
  common_rank <- min(rank_left, rank_right)
  singular <- if (common_rank) {
    svd(crossprod(left, right), nu = 0L, nv = 0L)$d
  } else numeric()
  singular <- pmin(1, pmax(0, singular))
  squared <- (abs(rank_left - rank_right) + sum(1 - singular^2)) / max_rank
  sqrt(max(0, min(1, squared)))
}

snapshot_metrics <- function(current, reference) {
  data.frame(
    fitted_nrmse_380_to_400 = nrmse(current$fitted, reference$fitted),
    rss_relative_l2_380_to_400 = relative_distance(
      current$rss, reference$rss
    ),
    variable_ppi_q95_abs_380_to_400 = quantile_abs(
      current$ppi, reference$ppi
    ),
    variable_ppi_max_abs_380_to_400 = max_abs(
      current$ppi, reference$ppi
    ),
    factor_ppi_max_abs_380_to_400 = max_abs(
      current$factor_ppi, reference$factor_ppi
    ), stringsAsFactors = FALSE
  )
}

direction_rows <- function(manifest_row, current, reference) {
  rows <- list()
  add <- function(block, unit, family, value) {
    rows[[length(rows) + 1L]] <<- data.frame(
      fit_id = manifest_row$fit_id[[1L]],
      data_id = manifest_row$data_id[[1L]],
      scenario_id = manifest_row$scenario_id[[1L]],
      seed_index = manifest_row$seed_index[[1L]],
      block = block, unit = unit, distance_family = family,
      distance_380_to_400 = as.numeric(value), stringsAsFactors = FALSE
    )
  }
  add("shared_loading", "shared", "subspace",
      subspace_distance(current$shared_loading, reference$shared_loading))
  for (study in seq_along(current$specific_loading)) {
    add("specific_loading", paste0("study", study), "subspace",
        subspace_distance(current$specific_loading[[study]],
                          reference$specific_loading[[study]]))
  }
  for (factor in seq_along(current$shared_feature)) {
    add("shared_feature", paste0("factor", factor), "subspace",
        subspace_distance(current$shared_feature[[factor]],
                          reference$shared_feature[[factor]]))
  }
  for (study in seq_along(current$specific_feature)) {
    for (factor in seq_along(current$specific_feature[[study]])) {
      add("specific_feature", paste0("study", study, "_factor", factor),
          "subspace", subspace_distance(
            current$specific_feature[[study]][[factor]],
            reference$specific_feature[[study]][[factor]]
          ))
    }
  }
  for (study in seq_along(current$shared_score_mean)) {
    for (factor in seq_along(current$shared_score_mean[[study]])) {
      unit <- paste0("study", study, "_factor", factor)
      add("shared_score_mean", unit, "subspace", subspace_distance(
        current$shared_score_mean[[study]][[factor]],
        reference$shared_score_mean[[study]][[factor]]
      ))
      add("shared_score_kernel", unit, "relative_frobenius",
          relative_distance(
            current$shared_score_kernel[[study]][[factor]],
            reference$shared_score_kernel[[study]][[factor]]
          ))
    }
  }
  for (study in seq_along(current$specific_score_mean)) {
    for (factor in seq_along(current$specific_score_mean[[study]])) {
      unit <- paste0("study", study, "_factor", factor)
      add("specific_score_mean", unit, "subspace", subspace_distance(
        current$specific_score_mean[[study]][[factor]],
        reference$specific_score_mean[[study]][[factor]]
      ))
      add("specific_score_kernel", unit, "relative_frobenius",
          relative_distance(
            current$specific_score_kernel[[study]][[factor]],
            reference$specific_score_kernel[[study]][[factor]]
          ))
    }
  }
  add("shared_contribution", "observed_vector", "relative_l2",
      relative_distance(current$shared_contribution,
                        reference$shared_contribution))
  add("specific_contribution", "observed_vector", "relative_l2",
      relative_distance(current$specific_contribution,
                        reference$specific_contribution))
  do.call(rbind, rows)
}

max_family <- function(frame, pattern) {
  values <- frame$distance_380_to_400[grepl(pattern, frame$block)]
  if (!length(values)) NA_real_ else max(values)
}

thresholds <- c(
  fitted = 0.003, rss = 0.006, loading = 0.001,
  feature = 0.02, score_mean = 0.02,
  score_kernel = 0.05, contribution = 0.01
)
per_fit <- vector("list", nrow(manifest))
all_direction <- vector("list", nrow(manifest))
all_warnings <- list()

for (index in seq_len(nrow(manifest))) {
  row <- manifest[index, , drop = FALSE]
  fit_dir <- file.path(UC_ROOT, "fits", row$fit_id[[1L]])
  record <- readRDS(file.path(fit_dir, "terminal_record.rds"))
  fit <- readRDS(file.path(fit_dir, "fit.rds"))
  directions <- readRDS(file.path(fit_dir, "DIRECTION_SNAPSHOTS.rds"))
  uc_assert(
    identical(record$fit_id, row$fit_id[[1L]]) &&
      identical(record$terminal_status, "fixed_400_complete") &&
      isTRUE(record$objective_eligible) && record$actual_T1_sweeps == 400L &&
      isTRUE(record$fixed_400_endpoint) &&
      isTRUE(record$public_g12_stopping_profile_applied) &&
      !isTRUE(record$truth_used) && !isTRUE(record$continuation_used) &&
      !isTRUE(record$automatic_800_used) &&
      identical(names(directions), as.character(c(0L, S4C_CHECKPOINTS))),
    paste0("Invalid fixed-400 bundle: ", row$fit_id[[1L]])
  )
  current <- fit$practical_checkpoints[["380"]]$snapshot
  reference <- fit$practical_checkpoints[["400"]]$snapshot
  fit_metrics <- snapshot_metrics(current, reference)
  direction <- direction_rows(
    row, directions[["380"]], directions[["400"]]
  )
  all_direction[[index]] <- direction
  loading <- max_family(direction, "loading$")
  feature <- max_family(direction, "feature$")
  score_mean <- max_family(direction, "score_mean$")
  score_kernel <- max_family(direction, "score_kernel$")
  contribution <- max_family(direction, "contribution$")
  safety <- c(
    fitted = fit_metrics$fitted_nrmse_380_to_400 <= thresholds[["fitted"]],
    rss = fit_metrics$rss_relative_l2_380_to_400 <= thresholds[["rss"]],
    loading = loading <= thresholds[["loading"]],
    feature = feature <= thresholds[["feature"]],
    score_mean = score_mean <= thresholds[["score_mean"]],
    score_kernel = score_kernel <= thresholds[["score_kernel"]],
    contribution = contribution <= thresholds[["contribution"]]
  )
  per_fit[[index]] <- data.frame(
    fit_id = row$fit_id[[1L]], data_id = row$data_id[[1L]],
    scenario_id = row$scenario_id[[1L]],
    seed_index = row$seed_index[[1L]], fit_seed = row$fit_seed[[1L]],
    objective_eligible = record$objective_eligible,
    final_elbo = record$final_elbo,
    elbo_gain_380_to_400 = fit$ELBO[[400L]] - fit$ELBO[[380L]],
    practical_converged_at_400 = record$practical_converged,
    slow_case_at_400 = record$slow_case,
    convergence_status = record$convergence_status,
    convergence_reason = record$convergence_reason,
    fit_metrics,
    loading_max_380_to_400 = loading,
    feature_max_380_to_400 = feature,
    score_mean_max_380_to_400 = score_mean,
    score_kernel_max_380_to_400 = score_kernel,
    contribution_max_380_to_400 = contribution,
    fitted_safe = safety[["fitted"]], rss_safe = safety[["rss"]],
    loading_safe = safety[["loading"]],
    feature_safe = safety[["feature"]],
    score_mean_safe = safety[["score_mean"]],
    score_kernel_safe = safety[["score_kernel"]],
    contribution_safe = safety[["contribution"]],
    all_scientific_outputs_locally_stable = all(safety),
    elapsed_seconds = record$elapsed_seconds,
    peak_memory_bytes = record$peak_memory_bytes,
    warning_count = if (is.data.frame(record$warnings))
      nrow(record$warnings) else 0L,
    truth_used = FALSE, continuation_used = FALSE,
    automatic_800_used = FALSE, stringsAsFactors = FALSE
  )
  if (is.data.frame(record$warnings) && nrow(record$warnings)) {
    warning_frame <- record$warnings
    warning_frame$fit_id <- row$fit_id[[1L]]
    all_warnings[[row$fit_id[[1L]]]] <- warning_frame
  }
}

per_fit <- do.call(rbind, per_fit)
direction <- do.call(rbind, all_direction)
rownames(per_fit) <- rownames(direction) <- NULL
warnings <- if (length(all_warnings)) do.call(rbind, all_warnings) else
  data.frame(fit_id = character(), stringsAsFactors = FALSE)
rownames(warnings) <- NULL

winner_rows <- lapply(split(per_fit, per_fit$data_id), function(group) {
  eligible <- group[group$objective_eligible & is.finite(group$final_elbo), ]
  uc_assert(nrow(eligible) == 2L,
            paste0("Two eligible starts are required for ", group$data_id[[1L]]))
  order_index <- order(eligible$final_elbo, decreasing = TRUE)
  winner <- eligible[order_index[[1L]], , drop = FALSE]
  runner_up <- eligible[order_index[[2L]], , drop = FALSE]
  data.frame(
    data_id = winner$data_id, scenario_id = winner$scenario_id,
    selected_fit_id = winner$fit_id,
    selected_seed_index = winner$seed_index,
    selected_fit_seed = winner$fit_seed,
    selected_final_elbo = winner$final_elbo,
    runner_up_fit_id = runner_up$fit_id,
    runner_up_final_elbo = runner_up$final_elbo,
    top_minus_runner_up_elbo = winner$final_elbo - runner_up$final_elbo,
    selection_used_truth = FALSE, stringsAsFactors = FALSE
  )
})
winners <- do.call(rbind, winner_rows)
winners <- winners[match(UC_DATA_IDS, winners$data_id), , drop = FALSE]
rownames(winners) <- NULL

scenario_summary <- do.call(rbind, lapply(
  split(per_fit, per_fit$scenario_id), function(group) data.frame(
    scenario_id = group$scenario_id[[1L]], trajectories = nrow(group),
    practical_converged_at_400 = sum(group$practical_converged_at_400),
    slow_case_at_400 = sum(group$slow_case_at_400),
    locally_stable_380_to_400 =
      sum(group$all_scientific_outputs_locally_stable),
    mean_elapsed_seconds = mean(group$elapsed_seconds),
    max_peak_memory_bytes = max(group$peak_memory_bytes, na.rm = TRUE),
    warning_rows = sum(group$warning_count), stringsAsFactors = FALSE
  )
))
rownames(scenario_summary) <- NULL

metric_summary <- data.frame(
  metric = c(
    "fitted_nrmse", "rss_relative_l2", "loading_subspace",
    "feature_subspace", "score_mean_subspace",
    "score_kernel_relative_frobenius", "contribution_relative_l2"
  ),
  threshold = unname(thresholds),
  median_380_to_400 = c(
    stats::median(per_fit$fitted_nrmse_380_to_400),
    stats::median(per_fit$rss_relative_l2_380_to_400),
    stats::median(per_fit$loading_max_380_to_400),
    stats::median(per_fit$feature_max_380_to_400),
    stats::median(per_fit$score_mean_max_380_to_400),
    stats::median(per_fit$score_kernel_max_380_to_400),
    stats::median(per_fit$contribution_max_380_to_400)
  ),
  maximum_380_to_400 = c(
    max(per_fit$fitted_nrmse_380_to_400),
    max(per_fit$rss_relative_l2_380_to_400),
    max(per_fit$loading_max_380_to_400),
    max(per_fit$feature_max_380_to_400),
    max(per_fit$score_mean_max_380_to_400),
    max(per_fit$score_kernel_max_380_to_400),
    max(per_fit$contribution_max_380_to_400)
  ),
  violating_trajectories = c(
    sum(!per_fit$fitted_safe), sum(!per_fit$rss_safe),
    sum(!per_fit$loading_safe), sum(!per_fit$feature_safe),
    sum(!per_fit$score_mean_safe), sum(!per_fit$score_kernel_safe),
    sum(!per_fit$contribution_safe)
  ), stringsAsFactors = FALSE
)

qc <- data.frame(
  check_id = c(
    "twelve_complete_fits", "all_objective_eligible", "all_fixed_400",
    "public_g12_stopping_profile_used", "six_truth_free_winners",
    "truth_never_used", "continuation_never_used",
    "automatic_800_never_used", "not_formal_v0lv"
  ),
  passed = c(
    nrow(per_fit) == 12L,
    all(per_fit$objective_eligible),
    all(vapply(file.path(
      UC_ROOT, "fits", manifest$fit_id, "terminal_record.rds"
    ), function(path) readRDS(path)$actual_T1_sweeps == 400L, logical(1L))),
    all(vapply(file.path(
      UC_ROOT, "fits", manifest$fit_id, "terminal_record.rds"
    ), function(path) {
      isTRUE(readRDS(path)$public_g12_stopping_profile_applied)
    }, logical(1L))),
    nrow(winners) == 6L && !any(winners$selection_used_truth),
    !any(per_fit$truth_used) && !any(manifest$truth_available_to_analysis),
    !any(per_fit$continuation_used) && !any(manifest$continuation),
    !any(per_fit$automatic_800_used) && !any(manifest$automatic_800),
    !any(manifest$formal_v0lv_result)
  ), stringsAsFactors = FALSE
)
uc_assert(all(qc$passed), "Stage-4C truth-free analysis QC failed.")

write_csv(per_fit, "PER_FIT_FIXED400_STABILITY.csv")
write_csv(direction, "DIRECTION_380_TO_400.csv")
write_csv(metric_summary, "METRIC_STABILITY_SUMMARY.csv")
write_csv(scenario_summary, "SCENARIO_STOPPING_SUMMARY.csv")
write_csv(winners, "TRUTH_FREE_ELBO_WINNERS.csv")
write_csv(warnings, "ALL_WARNINGS.csv")
write_csv(qc, "ANALYSIS_QC.csv")
writeLines(capture.output(sessionInfo()),
           file.path(analysis_dir, "SESSION_INFO.txt"), useBytes = TRUE)
writeLines(c(
  "status=COMPLETE", paste0("completed_utc=", uc_iso_time()),
  "fixed_400_fits=12", "objective_eligible=12",
  paste0("practical_converged_at_400=",
         sum(per_fit$practical_converged_at_400)),
  paste0("locally_stable_380_to_400=",
         sum(per_fit$all_scientific_outputs_locally_stable)),
  "truth_free_elbo_winners=6", "truth_read=FALSE",
  "continuation_used=FALSE", "automatic_800_used=FALSE",
  "formal_v0lv_result=FALSE"
), file.path(analysis_dir, "ANALYSIS_COMPLETE.txt"), useBytes = TRUE)
writeLines(c(
  "status=STAGE4C_FIXED400_INDEPENDENT_VALIDATION_COMPLETE",
  paste0("completed_utc=", uc_iso_time()),
  "fits=12", "analysis_truth_free=TRUE", "truth_unseal_authorized=FALSE",
  "scientific_truth_evaluation_started=FALSE",
  "continuation_started=FALSE", "automatic_800_started=FALSE"
), file.path(UC_ROOT, "STAGE4C_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE4C_ANALYSIS_COMPLETE practical=",
    sum(per_fit$practical_converged_at_400),
    "/12 locally_stable=",
    sum(per_fit$all_scientific_outputs_locally_stable),
    "/12 winners=6 truth=0\n", sep = "")
