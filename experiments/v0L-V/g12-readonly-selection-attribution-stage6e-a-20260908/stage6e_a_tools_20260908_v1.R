s6ea_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

s6ea_read_csv <- function(path) {
  s6ea_assert(file.exists(path), paste0("Missing input CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

s6ea_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, na = "")
  invisible(path)
}

s6ea_write_lines <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(as.character(value)), path, useBytes = TRUE)
  invisible(path)
}

s6ea_stratified_mean_se <- function(loss, study) {
  loss <- as.numeric(loss)
  study <- as.integer(study)
  s6ea_assert(length(loss) == length(study) && length(loss) >= 4L &&
                 all(is.finite(loss)) && all(is.finite(study)),
               "Invalid losses for stratified standard error.")
  split_loss <- split(loss, study)
  sizes <- vapply(split_loss, length, integer(1L))
  s6ea_assert(all(sizes >= 2L),
               "Every study needs at least two held-out units.")
  weights <- sizes / sum(sizes)
  variances <- vapply(split_loss, stats::var, numeric(1L))
  sqrt(sum(weights^2 * variances / sizes))
}

s6ea_group_apply <- function(data, keys, fun) {
  groups <- split(data, interaction(data[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, fun)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

s6ea_axis_metadata <- function(axis = c("factor", "fpca")) {
  axis <- match.arg(axis)
  if (axis == "factor") {
    data.frame(
      axis = axis,
      fit_config_id = c("factor_L1_M2", "factor_L2_M2", "factor_L3_M2"),
      complexity_value = 1:3,
      true_complexity_posthoc = 1L,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      axis = axis,
      fit_config_id = c("factor_L1_M2", "fpca_L1_M3", "fpca_L1_M4"),
      complexity_value = 2:4,
      true_complexity_posthoc = 2L,
      stringsAsFactors = FALSE
    )
  }
}

s6ea_validate_stage6d <- function(stage6d_root) {
  selection_root <- file.path(stage6d_root, "truth_free_selection")
  endpoints <- s6ea_read_csv(file.path(selection_root, "ALL_ENDPOINTS.csv"))
  winners <- s6ea_read_csv(file.path(selection_root, "CONFIG_WINNERS.csv"))
  science <- s6ea_read_csv(file.path(stage6d_root, "analysis",
                                     "START_LEVEL_SCIENTIFIC_RESULTS.csv"))
  dimensions <- s6ea_read_csv(file.path(selection_root,
                                        "DIMENSION_SELECTIONS.csv"))
  s6ea_assert(nrow(endpoints) == 180L && !anyDuplicated(endpoints$fit_id),
               "Stage 6D must contain 180 unique endpoints.")
  s6ea_assert(nrow(winners) == 15L && !anyDuplicated(winners$fit_id) &&
                 all(table(winners$data_id) == 5L),
               "Stage 6D must contain 15 configuration winners.")
  s6ea_assert(nrow(science) == 180L && !anyDuplicated(science$fit_id) &&
                 setequal(science$fit_id, endpoints$fit_id),
               "Stage 6D science table does not match all endpoints.")
  s6ea_assert(nrow(dimensions) == 3L &&
                 identical(sort(as.character(dimensions$data_id)),
                           sprintf("g12s6d_%02d", 1:3)),
               "Stage 6D dimension table is incomplete.")
  s6ea_assert(all(endpoints$objective_eligible) &&
                 all(endpoints$terminal_status == "fixed_400_complete") &&
                 !any(endpoints$truth_used_for_fit_stopping_or_selection) &&
                 !any(endpoints$continuation_used),
               "Stage 6D endpoint eligibility or truth-isolation failed.")
  fit_paths <- file.path(stage6d_root, "fits", endpoints$fit_id, "fit.rds")
  record_paths <- file.path(stage6d_root, "fits", endpoints$fit_id,
                            "terminal_record.rds")
  holdout_paths <- file.path(stage6d_root, "fits", endpoints$fit_id,
                             "HOLDOUT_SUBJECT_TIME_ALL.csv")
  s6ea_assert(all(file.exists(fit_paths)) && all(file.exists(record_paths)) &&
                 all(file.exists(holdout_paths)),
               "Stage 6D fit, terminal or holdout files are missing.")
  list(endpoints = endpoints, winners = winners, science = science,
       dimensions = dimensions, fit_paths = fit_paths,
       record_paths = record_paths, holdout_paths = holdout_paths)
}

s6ea_read_holdout_losses <- function(stage6d_root, endpoints) {
  ordered <- endpoints[order(endpoints$data_id, endpoints$fit_config_id,
                             -endpoints$final_elbo, endpoints$fit_id),
                       , drop = FALSE]
  ordered$elbo_rank <- ave(seq_len(nrow(ordered)),
                           interaction(ordered$data_id,
                                       ordered$fit_config_id, drop = TRUE),
                           FUN = seq_along)
  rows <- lapply(seq_len(nrow(ordered)), function(index) {
    endpoint <- ordered[index, , drop = FALSE]
    path <- file.path(stage6d_root, "fits", endpoint$fit_id,
                      "HOLDOUT_SUBJECT_TIME_ALL.csv")
    value <- s6ea_read_csv(path)
    required <- c("study", "subject", "holdout_time", "mse", "mae",
                  "observed_second_moment")
    s6ea_assert(nrow(value) == 60L && all(required %in% names(value)),
                 paste0("Invalid holdout table: ", endpoint$fit_id))
    value$fit_id <- endpoint$fit_id
    value$data_id <- endpoint$data_id
    value$fit_config_id <- endpoint$fit_config_id
    value$seed_index <- endpoint$seed_index
    value$final_elbo <- endpoint$final_elbo
    value$elbo_rank <- as.integer(endpoint$elbo_rank)
    value
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  s6ea_assert(nrow(result) == 10800L &&
                 all(table(interaction(result$fit_id, drop = TRUE)) == 60L),
               "Holdout loss stack is incomplete.")
  key_checks <- split(result, result$data_id)
  for (data_rows in key_checks) {
    by_fit <- split(data_rows, data_rows$fit_id)
    keys <- lapply(by_fit, function(value) paste(value$study, value$subject,
                                                   sep = ":"))
    reference <- sort(keys[[1L]])
    s6ea_assert(all(vapply(keys, function(key) identical(sort(key), reference),
                             logical(1L))),
                 "Paired holdout keys differ within a data set.")
  }
  result
}

s6ea_make_risk_units <- function(losses) {
  methods <- data.frame(
    risk_method = c("max_elbo_start", "top3_elbo_median",
                    "top5_elbo_median"),
    top_k = c(1L, 3L, 5L), stringsAsFactors = FALSE
  )
  rows <- list()
  index <- 0L
  groups <- split(losses, interaction(losses$data_id,
                                      losses$fit_config_id, drop = TRUE))
  for (group in groups) {
    for (method_index in seq_len(nrow(methods))) {
      method <- methods[method_index, , drop = FALSE]
      chosen <- group[group$elbo_rank <= method$top_k, , drop = FALSE]
      unit_groups <- split(chosen, interaction(chosen$study, chosen$subject,
                                               drop = TRUE))
      unit_rows <- lapply(unit_groups, function(unit) {
        unit <- unit[order(unit$elbo_rank), , drop = FALSE]
        data.frame(
          risk_method = method$risk_method,
          top_k = method$top_k,
          data_id = unit$data_id[[1L]],
          fit_config_id = unit$fit_config_id[[1L]],
          study = unit$study[[1L]],
          subject = unit$subject[[1L]],
          holdout_time = unit$holdout_time[[1L]],
          source_fit_count = nrow(unit),
          mse = stats::median(unit$mse),
          mae = stats::median(unit$mae),
          observed_second_moment = unit$observed_second_moment[[1L]],
          stringsAsFactors = FALSE
        )
      })
      index <- index + 1L
      rows[[index]] <- do.call(rbind, unit_rows)
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  source_group_sizes <- table(interaction(losses$data_id,
                                          losses$fit_config_id,
                                          losses$study, losses$subject,
                                          drop = TRUE))
  expected_units <- length(unique(interaction(losses$data_id,
                                               losses$fit_config_id,
                                               losses$study,
                                               losses$subject,
                                               drop = TRUE))) * nrow(methods)
  s6ea_assert(all(source_group_sizes >= max(methods$top_k)) &&
                 nrow(result) == expected_units &&
                 all(table(interaction(result$risk_method, result$data_id,
                                       result$fit_config_id,
                                       drop = TRUE)) ==
                       length(unique(interaction(losses$study, losses$subject,
                                                 drop = TRUE)))),
               "Risk-unit summaries are incomplete.")
  result
}

s6ea_summarize_risk <- function(risk_units) {
  s6ea_group_apply(risk_units,
                    c("risk_method", "data_id", "fit_config_id"),
                    function(group) {
    mean_mse <- mean(group$mse)
    observed_rms <- sqrt(mean(group$observed_second_moment))
    data.frame(
      risk_method = group$risk_method[[1L]],
      data_id = group$data_id[[1L]],
      fit_config_id = group$fit_config_id[[1L]],
      heldout_units = nrow(group),
      source_fit_count = unique(group$source_fit_count)[[1L]],
      mean_subject_time_mse = mean_mse,
      stratified_se_mse = s6ea_stratified_mean_se(group$mse, group$study),
      heldout_rmse = sqrt(mean_mse),
      heldout_nrmse = sqrt(mean_mse) / observed_rms,
      heldout_mae = mean(group$mae),
      observed_rms = observed_rms,
      stringsAsFactors = FALSE
    )
  })
}

s6ea_compare_selection_rules <- function(risk_units) {
  candidate_rows <- list()
  selection_rows <- list()
  candidate_index <- 0L
  selection_index <- 0L
  for (method in unique(risk_units$risk_method)) {
    for (data_id in unique(risk_units$data_id)) {
      for (axis in c("factor", "fpca")) {
        metadata <- s6ea_axis_metadata(axis)
        rows <- risk_units[
          risk_units$risk_method == method & risk_units$data_id == data_id &
            risk_units$fit_config_id %in% metadata$fit_config_id,
          , drop = FALSE
        ]
        s6ea_assert(length(unique(rows$fit_config_id)) == 3L,
                     paste0("Missing ", axis, " configurations for ", data_id,
                            " and ", method, "."))
        summaries <- s6ea_summarize_risk(rows)
        summaries <- merge(summaries, metadata, by = "fit_config_id",
                           sort = FALSE)
        summaries <- summaries[order(summaries$complexity_value), , drop = FALSE]
        best <- summaries[order(summaries$mean_subject_time_mse,
                                summaries$complexity_value,
                                summaries$fit_config_id), , drop = FALSE][1L, ]
        best_units <- rows[rows$fit_config_id == best$fit_config_id,
                           , drop = FALSE]
        best_units <- best_units[order(best_units$study, best_units$subject), ]
        independent_threshold <- best$mean_subject_time_mse +
          best$stratified_se_mse
        summaries$mean_paired_difference_vs_min <- NA_real_
        summaries$paired_se_difference_vs_min <- NA_real_
        summaries$paired_ci95_low <- NA_real_
        summaries$paired_ci95_high <- NA_real_
        for (i in seq_len(nrow(summaries))) {
          config_units <- rows[rows$fit_config_id == summaries$fit_config_id[[i]],
                               , drop = FALSE]
          config_units <- config_units[order(config_units$study,
                                             config_units$subject), ]
          s6ea_assert(identical(paste(config_units$study, config_units$subject),
                                paste(best_units$study, best_units$subject)),
                       "Paired unit keys differ during rule comparison.")
          difference <- config_units$mse - best_units$mse
          difference_se <- s6ea_stratified_mean_se(difference,
                                                    config_units$study)
          summaries$mean_paired_difference_vs_min[[i]] <- mean(difference)
          summaries$paired_se_difference_vs_min[[i]] <- difference_se
          summaries$paired_ci95_low[[i]] <- mean(difference) -
            1.96 * difference_se
          summaries$paired_ci95_high[[i]] <- mean(difference) +
            1.96 * difference_se
        }
        summaries$minimum_mean_loss <- summaries$fit_config_id ==
          best$fit_config_id
        summaries$within_independent_one_se <-
          summaries$mean_subject_time_mse <= independent_threshold +
          10 * .Machine$double.eps
        summaries$within_paired_difference_one_se <-
          summaries$mean_paired_difference_vs_min <=
          summaries$paired_se_difference_vs_min + 10 * .Machine$double.eps
        summaries$selection_inputs_truth_free <- TRUE
        summaries$truth_used_only_for_posthoc_correctness <- TRUE
        candidate_index <- candidate_index + 1L
        candidate_rows[[candidate_index]] <- summaries

        for (rule in c("independent_one_se", "paired_difference_one_se")) {
          eligible <- if (rule == "independent_one_se") {
            summaries$within_independent_one_se
          } else summaries$within_paired_difference_one_se
          selected <- summaries[eligible, , drop = FALSE]
          selected <- selected[order(selected$complexity_value,
                                     selected$fit_config_id), , drop = FALSE][1L, ]
          selection_index <- selection_index + 1L
          selection_rows[[selection_index]] <- data.frame(
            risk_method = method,
            data_id = data_id,
            axis = axis,
            one_se_semantics = rule,
            minimum_loss_config_id = best$fit_config_id,
            selected_config_id = selected$fit_config_id,
            selected_complexity = selected$complexity_value,
            true_complexity_posthoc = selected$true_complexity_posthoc,
            correct_dimension_posthoc = selected$complexity_value ==
              selected$true_complexity_posthoc,
            selection_inputs_truth_free = TRUE,
            truth_used_only_for_posthoc_correctness = TRUE,
            rule_status = "exploratory_posthoc_not_validated",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  candidates <- do.call(rbind, candidate_rows)
  selections <- do.call(rbind, selection_rows)
  rownames(candidates) <- NULL
  rownames(selections) <- NULL
  list(candidates = candidates, selections = selections)
}

s6ea_top_set_diagnostics <- function(losses, science, risk_summary) {
  fit_risk <- s6ea_group_apply(
    losses, c("fit_id"), function(group) data.frame(
      fit_id = group$fit_id[[1L]],
      mean_holdout_mse = mean(group$mse),
      heldout_nrmse_recomputed = sqrt(mean(group$mse)) /
        sqrt(mean(group$observed_second_moment)),
      stringsAsFactors = FALSE
    )
  )
  starts <- merge(unique(losses[c("fit_id", "data_id", "fit_config_id",
                                  "seed_index", "final_elbo", "elbo_rank")]),
                  fit_risk, by = "fit_id", sort = FALSE)
  truth_columns <- c("fit_id", "core_truth_recovered",
                     "exact_factor_structure_recovered", "true_effective_M",
                     "selected_counts", "dense_signal_nrmse_ppi_selected",
                     "feature_total_ise_mean",
                     "covariance_operator_relative_error_mean")
  starts <- merge(starts, science[truth_columns], by = "fit_id", sort = FALSE)
  rows <- list()
  index <- 0L
  methods <- data.frame(
    risk_method = c("max_elbo_start", "top3_elbo_median",
                    "top5_elbo_median"), top_k = c(1L, 3L, 5L),
    stringsAsFactors = FALSE
  )
  for (group in split(starts, interaction(starts$data_id,
                                          starts$fit_config_id, drop = TRUE))) {
    for (i in seq_len(nrow(methods))) {
      method <- methods[i, ]
      selected <- group[group$elbo_rank <= method$top_k, , drop = FALSE]
      risk <- risk_summary[
        risk_summary$risk_method == method$risk_method &
          risk_summary$data_id == group$data_id[[1L]] &
          risk_summary$fit_config_id == group$fit_config_id[[1L]],
        , drop = FALSE
      ]
      index <- index + 1L
      rows[[index]] <- data.frame(
        data_id = group$data_id[[1L]],
        fit_config_id = group$fit_config_id[[1L]],
        risk_method = method$risk_method,
        top_k = method$top_k,
        start_count = nrow(selected),
        elbo_rank_range = paste(range(selected$elbo_rank), collapse = ":"),
        core_truth_hits_posthoc = sum(selected$core_truth_recovered),
        core_truth_rate_posthoc = mean(selected$core_truth_recovered),
        exact_factor_hits_posthoc =
          sum(selected$exact_factor_structure_recovered),
        exact_factor_rate_posthoc =
          mean(selected$exact_factor_structure_recovered),
        effective_M2_hits_posthoc = sum(selected$true_effective_M),
        effective_M2_rate_posthoc = mean(selected$true_effective_M),
        median_start_holdout_nrmse =
          stats::median(selected$heldout_nrmse_recomputed),
        robust_unit_heldout_nrmse = risk$heldout_nrmse,
        median_dense_nrmse_posthoc =
          stats::median(selected$dense_signal_nrmse_ppi_selected),
        median_feature_ise_posthoc =
          stats::median(selected$feature_total_ise_mean),
        median_covariance_error_posthoc = stats::median(
          selected$covariance_operator_relative_error_mean),
        selection_of_top_set_used_truth = FALSE,
        truth_metrics_joined_posthoc = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

s6ea_load_truth_free_fits <- function(stage6d_root, inputs) {
  candidate_rows <- list()
  component_rows <- list()
  match_rows <- list()
  loading_subspace_rows <- list()
  fpca_subspace_rows <- list()
  duplicate_rows <- list()
  leakage_rows <- list()
  candidate_index <- component_index <- match_index <- 0L
  loading_index <- fpca_index <- duplicate_index <- leakage_index <- 0L
  endpoint_groups <- split(inputs$endpoints,
                           interaction(inputs$endpoints$data_id,
                                       inputs$endpoints$fit_config_id,
                                       drop = TRUE))
  for (endpoint_group in endpoint_groups) {
    endpoint_group <- endpoint_group[order(endpoint_group$seed_index), ]
    extracted <- setNames(vector("list", nrow(endpoint_group)),
                          endpoint_group$fit_id)
    for (i in seq_len(nrow(endpoint_group))) {
      endpoint <- endpoint_group[i, ]
      fit <- readRDS(file.path(stage6d_root, "fits", endpoint$fit_id,
                               "fit.rds"))
      record <- readRDS(file.path(stage6d_root, "fits", endpoint$fit_id,
                                  "terminal_record.rds"))
      s6ea_assert(identical(record$fit_id, endpoint$fit_id) &&
                     isTRUE(record$objective_eligible) &&
                     isFALSE(record$truth_used_for_fit_stopping_or_selection),
                   "A terminal record failed truth-free eligibility.")
      value <- s6b_extract_fit(fit, record)
      value$candidates$data_id <- endpoint$data_id
      value$components$data_id <- endpoint$data_id
      candidate_index <- candidate_index + 1L
      component_index <- component_index + 1L
      candidate_rows[[candidate_index]] <- value$candidates
      component_rows[[component_index]] <- value$components
      extracted[[endpoint$fit_id]] <- value
      rm(fit)
      invisible(gc(FALSE))
    }
    winner_id <- inputs$winners$fit_id[
      inputs$winners$data_id == endpoint_group$data_id[[1L]] &
        inputs$winners$fit_config_id == endpoint_group$fit_config_id[[1L]]
    ]
    s6ea_assert(length(winner_id) == 1L && winner_id %in% names(extracted),
                 "The configuration winner is missing from extracted fits.")
    reference <- extracted[[winner_id]]
    roles <- sort(unique(reference$candidates$role_id))
    for (fit_id in names(extracted)) {
      candidate <- extracted[[fit_id]]
      for (role in roles) {
        reference_objects <- s6b_role_objects(reference, role)
        candidate_objects <- s6b_role_objects(candidate, role)
        s6ea_assert(length(reference_objects) == length(candidate_objects),
                     "Configured factor count changed across starts.")
        similarity <- outer(seq_along(reference_objects),
                            seq_along(candidate_objects),
                            Vectorize(function(left, right) {
          s6b_abs_cosine(reference_objects[[left]]$loading,
                         candidate_objects[[right]]$loading)
        }))
        assignment <- s6b_best_assignment(similarity)
        loading_similarity <- s6b_euclidean_subspace_similarity(
          do.call(cbind, lapply(reference_objects, `[[`, "loading")),
          do.call(cbind, lapply(candidate_objects, `[[`, "loading"))
        )
        loading_index <- loading_index + 1L
        loading_subspace_rows[[loading_index]] <- data.frame(
          data_id = endpoint_group$data_id[[1L]],
          fit_config_id = endpoint_group$fit_config_id[[1L]],
          reference_fit_id = winner_id,
          candidate_fit_id = fit_id,
          seed_index = candidate$record$seed_index,
          role_id = role,
          mean_squared_canonical_correlation =
            loading_similarity[["mean_squared_canonical_correlation"]],
          minimum_canonical_correlation =
            loading_similarity[["minimum_canonical_correlation"]],
          reference_rank = loading_similarity[["rank_left"]],
          candidate_rank = loading_similarity[["rank_right"]],
          stringsAsFactors = FALSE
        )
        for (left in seq_along(reference_objects)) {
          right <- assignment[[left]]
          object <- candidate_objects[[right]]
          candidate_row <- candidate$candidates[
            candidate$candidates$role_id == role &
              candidate$candidates$candidate_id == object$candidate_id,
            , drop = FALSE
          ]
          match_index <- match_index + 1L
          match_rows[[match_index]] <- data.frame(
            data_id = endpoint_group$data_id[[1L]],
            fit_config_id = endpoint_group$fit_config_id[[1L]],
            reference_fit_id = winner_id,
            candidate_fit_id = fit_id,
            seed_index = candidate$record$seed_index,
            role_id = role,
            reference_candidate_id = reference_objects[[left]]$candidate_id,
            candidate_id = object$candidate_id,
            loading_abs_cosine = similarity[left, right],
            candidate_contribution_share =
              candidate_row$weighted_contribution_share,
            candidate_factor_ppi = candidate_row$factor_ppi,
            stringsAsFactors = FALSE
          )
          fpca_similarity <- s6b_subspace_similarity(
            reference_objects[[left]]$functions, object$functions,
            reference$time_g
          )
          fpca_index <- fpca_index + 1L
          fpca_subspace_rows[[fpca_index]] <- data.frame(
            data_id = endpoint_group$data_id[[1L]],
            fit_config_id = endpoint_group$fit_config_id[[1L]],
            reference_fit_id = winner_id,
            candidate_fit_id = fit_id,
            seed_index = candidate$record$seed_index,
            role_id = role,
            reference_candidate_id = reference_objects[[left]]$candidate_id,
            candidate_id = object$candidate_id,
            mean_squared_canonical_correlation =
              fpca_similarity[["mean_squared_canonical_correlation"]],
            minimum_canonical_correlation =
              fpca_similarity[["minimum_canonical_correlation"]],
            reference_rank = fpca_similarity[["rank_left"]],
            candidate_rank = fpca_similarity[["rank_right"]],
            stringsAsFactors = FALSE
          )
        }
        if (length(candidate_objects) >= 2L) {
          pairs <- utils::combn(seq_along(candidate_objects), 2L)
          cosines <- apply(pairs, 2L, function(pair) {
            s6b_abs_cosine(candidate_objects[[pair[[1L]]]]$loading,
                           candidate_objects[[pair[[2L]]]]$loading)
          })
        } else cosines <- numeric()
        duplicate_index <- duplicate_index + 1L
        duplicate_rows[[duplicate_index]] <- data.frame(
          data_id = endpoint_group$data_id[[1L]],
          fit_config_id = endpoint_group$fit_config_id[[1L]],
          fit_id = fit_id,
          seed_index = candidate$record$seed_index,
          role_id = role,
          pair_count = length(cosines),
          maximum_within_role_abs_cosine =
            if (length(cosines)) max(cosines) else 0,
          pairs_ge_0_9 = sum(cosines >= 0.9),
          stringsAsFactors = FALSE
        )
      }
      shared <- s6b_role_objects(candidate, "A")
      for (study in 1:2) {
        specific <- s6b_role_objects(candidate, paste0("B", study))
        cosines <- outer(seq_along(shared), seq_along(specific),
                         Vectorize(function(left, right) {
          s6b_abs_cosine(shared[[left]]$loading,
                         specific[[right]]$loading)
        }))
        leakage_index <- leakage_index + 1L
        leakage_rows[[leakage_index]] <- data.frame(
          data_id = endpoint_group$data_id[[1L]],
          fit_config_id = endpoint_group$fit_config_id[[1L]],
          fit_id = fit_id,
          seed_index = candidate$record$seed_index,
          study = study,
          pair_count = length(cosines),
          maximum_cross_role_abs_cosine = max(cosines),
          pairs_ge_0_9 = sum(cosines >= 0.9),
          stringsAsFactors = FALSE
        )
      }
    }
    rm(extracted, reference)
    invisible(gc(FALSE))
  }
  list(
    candidates = do.call(rbind, candidate_rows),
    components = do.call(rbind, component_rows),
    matches = do.call(rbind, match_rows),
    loading_subspace = do.call(rbind, loading_subspace_rows),
    fpca_subspace = do.call(rbind, fpca_subspace_rows),
    duplicates = do.call(rbind, duplicate_rows),
    leakage = do.call(rbind, leakage_rows)
  )
}

s6ea_summarize_truth_free_complexity <- function(extracted) {
  candidates <- extracted$candidates
  components <- extracted$components
  fit_rows <- lapply(split(candidates, candidates$fit_id), function(group) {
    shares <- group$weighted_contribution_share
    role_tail <- sum(vapply(split(shares, group$role_id), function(value) {
      sum(value) - max(value)
    }, numeric(1L)))
    component_group <- components[components$fit_id == group$fit_id[[1L]], ]
    spectral_groups <- split(component_group,
                             interaction(component_group$role_id,
                                         component_group$candidate_id,
                                         drop = TRUE))
    spectral <- lapply(spectral_groups, function(value) {
      eigen <- as.numeric(value$normalized_eigenvalue)
      eigen[!is.finite(eigen)] <- 0
      eigen <- pmax(0, eigen)
      total <- sum(eigen, na.rm = TRUE)
      data.frame(
        last_share = if (total > 0) eigen[[nrow(value)]] / total else 0,
        beyond_two_share = if (total > 0) {
          sum(eigen[value$component >= 3L], na.rm = TRUE) / total
        } else 0
      )
    })
    spectral <- do.call(rbind, spectral)
    data.frame(
      fit_id = group$fit_id[[1L]],
      data_id = group$data_id[[1L]],
      fit_config_id = group$fit_config_id[[1L]],
      seed_index = group$seed_index[[1L]],
      final_elbo = group$final_elbo[[1L]],
      practical_converged = group$practical_converged[[1L]],
      candidate_count = nrow(group),
      ppi_retained_count = sum(group$ppi_retained_0_5),
      strongest_contribution_share = max(shares),
      within_role_tail_contribution_share = role_tail,
      inverse_simpson_effective_candidates =
        if (sum(shares^2) > 0) 1 / sum(shares^2) else 0,
      effective_M_mean = mean(group$effective_M_99pct),
      effective_M_at_cap_rate = mean(group$effective_M_99pct == group$rank_cap),
      last_component_spectral_share_mean = mean(spectral$last_share),
      beyond_M2_spectral_share_mean = mean(spectral$beyond_two_share),
      truth_fields_used = FALSE,
      stringsAsFactors = FALSE
    )
  })
  fit_summary <- do.call(rbind, fit_rows)
  rownames(fit_summary) <- NULL
  loading <- extracted$loading_subspace[
    extracted$loading_subspace$candidate_fit_id !=
      extracted$loading_subspace$reference_fit_id, , drop = FALSE]
  matches <- extracted$matches[
    extracted$matches$candidate_fit_id != extracted$matches$reference_fit_id,
    , drop = FALSE]
  fpca <- extracted$fpca_subspace[
    extracted$fpca_subspace$candidate_fit_id !=
      extracted$fpca_subspace$reference_fit_id, , drop = FALSE]
  config_rows <- lapply(split(fit_summary,
                              interaction(fit_summary$data_id,
                                          fit_summary$fit_config_id,
                                          drop = TRUE)), function(group) {
    data_id <- group$data_id[[1L]]
    config <- group$fit_config_id[[1L]]
    load_group <- loading[loading$data_id == data_id &
                            loading$fit_config_id == config, ]
    match_group <- matches[matches$data_id == data_id &
                             matches$fit_config_id == config, ]
    fpca_group <- fpca[fpca$data_id == data_id &
                         fpca$fit_config_id == config, ]
    duplicate_group <- extracted$duplicates[
      extracted$duplicates$data_id == data_id &
        extracted$duplicates$fit_config_id == config, ]
    leakage_group <- extracted$leakage[
      extracted$leakage$data_id == data_id &
        extracted$leakage$fit_config_id == config, ]
    direction_recurrence <- s6ea_group_apply(
      match_group, c("reference_candidate_id", "role_id"), function(value) {
        data.frame(recurrence = mean(value$loading_abs_cosine >= 0.9))
      })$recurrence
    data.frame(
      data_id = data_id,
      fit_config_id = config,
      fits = nrow(group),
      loading_subspace_mscc_median = stats::median(
        load_group$mean_squared_canonical_correlation),
      loading_subspace_min_cc_median = stats::median(
        load_group$minimum_canonical_correlation),
      direction_recurrence_ge_0_9_median = stats::median(direction_recurrence),
      direction_recurrence_ge_0_9_minimum = min(direction_recurrence),
      fit_role_duplicate_rate_ge_0_9 = mean(duplicate_group$pairs_ge_0_9 > 0),
      fit_study_cross_role_leakage_rate_ge_0_9 =
        mean(leakage_group$pairs_ge_0_9 > 0),
      fpca_subspace_mscc_median = stats::median(
        fpca_group$mean_squared_canonical_correlation),
      ppi_retained_count_median = stats::median(group$ppi_retained_count),
      strongest_contribution_share_median = stats::median(
        group$strongest_contribution_share),
      within_role_tail_contribution_share_median = stats::median(
        group$within_role_tail_contribution_share),
      inverse_simpson_effective_candidates_median = stats::median(
        group$inverse_simpson_effective_candidates),
      effective_M_at_cap_rate = mean(group$effective_M_at_cap_rate),
      last_component_spectral_share_median = stats::median(
        group$last_component_spectral_share_mean),
      beyond_M2_spectral_share_median = stats::median(
        group$beyond_M2_spectral_share_mean),
      truth_fields_used = FALSE,
      stringsAsFactors = FALSE
    )
  })
  config_summary <- do.call(rbind, config_rows)
  rownames(config_summary) <- NULL
  list(fit = fit_summary, config = config_summary)
}

s6ea_posthoc_complexity_association <- function(fit_summary, science) {
  columns <- c("fit_id", "core_truth_recovered",
               "exact_factor_structure_recovered", "true_effective_M",
               "heldout_nrmse_all", "dense_signal_nrmse_ppi_selected",
               "feature_total_ise_mean")
  joined <- merge(fit_summary, science[columns], by = "fit_id", sort = FALSE)
  metrics <- c("ppi_retained_count", "within_role_tail_contribution_share",
               "inverse_simpson_effective_candidates",
               "effective_M_at_cap_rate",
               "last_component_spectral_share_mean",
               "beyond_M2_spectral_share_mean")
  rows <- list()
  index <- 0L
  for (config in sort(unique(joined$fit_config_id))) {
    config_rows <- joined[joined$fit_config_id == config, ]
    for (metric in metrics) {
      hit <- config_rows[[metric]][config_rows$core_truth_recovered]
      miss <- config_rows[[metric]][!config_rows$core_truth_recovered]
      index <- index + 1L
      rows[[index]] <- data.frame(
        fit_config_id = config,
        diagnostic = metric,
        core_hit_n = length(hit),
        core_miss_n = length(miss),
        core_hit_mean_posthoc = if (length(hit)) mean(hit) else NA_real_,
        core_miss_mean_posthoc = if (length(miss)) mean(miss) else NA_real_,
        hit_minus_miss_posthoc = if (length(hit) && length(miss)) {
          mean(hit) - mean(miss)
        } else NA_real_,
        truth_join_status = "posthoc_exploratory",
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
