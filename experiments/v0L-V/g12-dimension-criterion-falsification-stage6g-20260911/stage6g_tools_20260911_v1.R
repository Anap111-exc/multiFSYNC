s6g_abort <- function(...) stop(paste0(...), call. = FALSE)
s6g_assert <- function(value, message) {
  if (!isTRUE(value)) s6g_abort(message)
  invisible(TRUE)
}

s6g_read_csv <- function(path) {
  s6g_assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

s6g_validate_registration <- function(scenarios, configs, splits, starts) {
  scenario_columns <- c(
    "scenario_index", "scenario_id", "data_seed", "truth_basis_seed",
    "expected_truth_fit_config_id", "truth_basis_method",
    "mean_truth_method", "truth_L_f", "truth_L_s_1", "truth_L_s_2",
    "truth_M_common", "S", "n_1", "n_2", "p", "d", "K",
    "n_time_min", "n_time_max", "sigma_eps", "mean_amp",
    "active_loading_fraction", "identifiability", "score_variance_rule",
    "formal_paper_mc", "execution_authorized"
  )
  config_columns <- c(
    "config_index", "fit_config_id", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_common", "fit_M_f", "fit_M_s_1",
    "fit_M_s_2", "total_candidate_factors", "complexity_proxy",
    "included_in_all_scenarios"
  )
  split_columns <- c(
    "split_index", "split_id", "holdout_rule", "paired_within_subject",
    "minimum_training_times"
  )
  start_columns <- c(
    "seed_index", "fit_seed",
    "paired_across_scenarios_splits_and_dimensions"
  )
  s6g_assert(identical(names(scenarios), scenario_columns),
             "Stage 6G scenario columns differ from registration.")
  s6g_assert(identical(names(configs), config_columns),
             "Stage 6G configuration columns differ from registration.")
  s6g_assert(identical(names(splits), split_columns),
             "Stage 6G split columns differ from registration.")
  s6g_assert(identical(names(starts), start_columns),
             "Stage 6G start columns differ from registration.")
  s6g_assert(
    nrow(scenarios) == 3L &&
      identical(as.integer(scenarios$scenario_index), 1:3) &&
      !anyDuplicated(scenarios$scenario_id) &&
      !anyDuplicated(scenarios$data_seed) &&
      !anyDuplicated(scenarios$truth_basis_seed) &&
      !any(scenarios$formal_paper_mc) &&
      !any(scenarios$execution_authorized),
    "Stage 6G must retain three unexecuted development scenarios."
  )
  s6g_assert(
    nrow(configs) == 5L &&
      identical(as.integer(configs$config_index), 1:5) &&
      !anyDuplicated(configs$fit_config_id) &&
      all(configs$included_in_all_scenarios) &&
      setequal(scenarios$expected_truth_fit_config_id,
               c("cand_L111_M2", "cand_L211_M2", "cand_L111_M3")),
    "Stage 6G must retain the registered five-candidate common grid."
  )
  expected_q <- as.integer(configs$fit_L_f) +
    as.integer(configs$fit_L_s_1) + as.integer(configs$fit_L_s_2)
  s6g_assert(identical(as.integer(configs$total_candidate_factors), expected_q),
             "Stage 6G candidate factor totals are inconsistent.")
  expected_complexity <- 500L * expected_q +
    as.integer(configs$fit_M_common) * (
      60L * as.integer(configs$fit_L_f) +
        30L * as.integer(configs$fit_L_s_1) +
        30L * as.integer(configs$fit_L_s_2) + 7L * expected_q
    )
  s6g_assert(identical(as.integer(configs$complexity_proxy),
                       expected_complexity),
             "Stage 6G complexity proxy differs from registration.")
  s6g_assert(
    nrow(splits) == 2L &&
      identical(as.integer(splits$split_index), 1:2) &&
      identical(as.character(splits$split_id),
                c("split1_middle", "split2_upper")) &&
      all(splits$paired_within_subject) &&
      all(as.integer(splits$minimum_training_times) == 5L),
    "Stage 6G holdout registration is invalid."
  )
  s6g_assert(
    nrow(starts) == 4L &&
      identical(as.integer(starts$seed_index), 1:4) &&
      !anyDuplicated(starts$fit_seed) &&
      all(starts$paired_across_scenarios_splits_and_dimensions),
    "Stage 6G must retain four paired starts."
  )
  s6g_assert(
    all(as.integer(scenarios$S) == 2L) &&
      all(as.integer(scenarios$n_1) == 30L) &&
      all(as.integer(scenarios$n_2) == 30L) &&
      all(as.integer(scenarios$p) == 500L) &&
      all(as.integer(scenarios$d) == 0L) &&
      all(as.integer(scenarios$K) == 5L) &&
      all(as.integer(scenarios$n_time_min) == 6L) &&
      all(as.integer(scenarios$n_time_max) == 9L) &&
      all(scenarios$truth_basis_method == "jaoua_splinets") &&
      all(scenarios$identifiability == "strong") &&
      all(abs(scenarios$active_loading_fraction - 0.1) < 1e-14),
    "Stage 6G common data design differs from the protocol."
  )
  invisible(TRUE)
}

s6g_build_fit_manifest <- function(scenarios, configs, splits, starts) {
  s6g_validate_registration(scenarios, configs, splits, starts)
  grid <- expand.grid(
    scenario_index = as.integer(scenarios$scenario_index),
    split_index = as.integer(splits$split_index),
    config_index = as.integer(configs$config_index),
    seed_index = as.integer(starts$seed_index),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$scenario_index, grid$split_index,
                     grid$config_index, grid$seed_index), , drop = FALSE]
  scenario_lookup <- scenarios[match(grid$scenario_index,
                                     scenarios$scenario_index), , drop = FALSE]
  split_lookup <- splits[match(grid$split_index,
                               splits$split_index), , drop = FALSE]
  config_lookup <- configs[match(grid$config_index,
                                 configs$config_index), , drop = FALSE]
  start_lookup <- starts[match(grid$seed_index,
                               starts$seed_index), , drop = FALSE]
  data_id <- as.character(scenario_lookup$scenario_id)
  split_id <- as.character(split_lookup$split_id)
  config_id <- as.character(config_lookup$fit_config_id)
  fit_id <- paste(
    data_id, split_id, config_id,
    sprintf("G12_%02d", as.integer(start_lookup$seed_index)), sep = "__"
  )
  result <- data.frame(
    experiment_id = "G12_STAGE6G_DIMENSION_CRITERION_FALSIFICATION_V1_20260911",
    fit_id = fit_id,
    data_id = data_id,
    data_index = as.integer(scenario_lookup$scenario_index),
    data_seed = as.integer(scenario_lookup$data_seed),
    scenario_id = data_id,
    split_id = split_id,
    split_index = as.integer(split_lookup$split_index),
    method_id = "G12",
    fit_config_id = config_id,
    selection_stratum_id = paste(data_id, split_id, config_id, sep = "__"),
    route_id = "random__pre1__jaoua_anneal__G12__fixed400",
    seed_index = as.integer(start_lookup$seed_index),
    fit_seed = as.integer(start_lookup$fit_seed),
    fit_L_f = as.integer(config_lookup$fit_L_f),
    fit_L_s_1 = as.integer(config_lookup$fit_L_s_1),
    fit_L_s_2 = as.integer(config_lookup$fit_L_s_2),
    fit_M_f = as.character(config_lookup$fit_M_f),
    fit_M_s_1 = as.character(config_lookup$fit_M_s_1),
    fit_M_s_2 = as.character(config_lookup$fit_M_s_2),
    fit_M_common = as.integer(config_lookup$fit_M_common),
    complexity_proxy = as.integer(config_lookup$complexity_proxy),
    n_cpus = 1L,
    planned_annealing_sweeps = 99L,
    planned_ordinary_t1_sweeps = 400L,
    objective_eligibility_required = TRUE,
    truth_free_winner_candidate = TRUE,
    cross_dimension_elbo_comparison_forbidden = TRUE,
    actual_racing = FALSE,
    continuation = FALSE,
    automatic_800 = FALSE,
    truth_available_to_fit = FALSE,
    truth_available_to_stopping = FALSE,
    truth_available_to_selection = FALSE,
    engineering_pilot = TRUE,
    formal_v0lv_result = FALSE,
    formal_paper_mc_result = FALSE,
    stringsAsFactors = FALSE
  )
  truth_columns <- c(
    "truth_basis_seed", "expected_truth_fit_config_id", "truth_L_f",
    "truth_L_s_1", "truth_L_s_2", "truth_M_common"
  )
  s6g_assert(
    nrow(result) == 120L && !anyDuplicated(result$fit_id) &&
      length(unique(result$selection_stratum_id)) == 30L &&
      all(table(result$selection_stratum_id) == 4L) &&
      all(table(result$data_id) == 40L) &&
      all(table(interaction(result$data_id, result$split_id,
                            result$fit_config_id, drop = TRUE)) == 4L) &&
      all(result$n_cpus == 1L) &&
      !any(result$truth_available_to_fit) &&
      !any(result$truth_available_to_stopping) &&
      !any(result$truth_available_to_selection) &&
      !length(intersect(names(result), truth_columns)),
    "Stage 6G fit-manifest accounting or truth isolation failed."
  )
  rownames(result) <- NULL
  result
}

s6g_stratified_mean <- function(loss, study) {
  groups <- split(as.numeric(loss), as.integer(study))
  s6g_assert(length(groups) == 2L && all(lengths(groups) >= 1L) &&
               all(vapply(groups, function(value) all(is.finite(value)),
                          logical(1L))),
             "Stage 6G stratified mean input is invalid.")
  mean(vapply(groups, mean, numeric(1L)))
}

s6g_stratified_se <- function(loss, study) {
  groups <- split(as.numeric(loss), as.integer(study))
  sizes <- lengths(groups)
  s6g_assert(length(groups) == 2L && all(sizes >= 2L) &&
               all(vapply(groups, function(value) all(is.finite(value)),
                          logical(1L))),
             "Stage 6G stratified SE input is invalid.")
  sqrt(sum((0.5^2) * vapply(groups, stats::var, numeric(1L)) / sizes))
}

s6g_select_one_se <- function(subject_losses, config_metadata, loss_column,
                              expected_subjects = 60L) {
  required <- c("fit_config_id", "study", "subject", loss_column)
  s6g_assert(all(required %in% names(subject_losses)) &&
               all(c("fit_config_id", "complexity_proxy") %in%
                     names(config_metadata)),
             "Stage 6G one-SE inputs are incomplete.")
  config_ids <- as.character(config_metadata$fit_config_id)
  s6g_assert(length(config_ids) == 5L && !anyDuplicated(config_ids),
             "Stage 6G one-SE selection requires five configurations.")
  reference_key <- NULL
  summaries <- lapply(config_ids, function(config_id) {
    rows <- subject_losses[
      subject_losses$fit_config_id == config_id, , drop = FALSE
    ]
    rows <- rows[order(rows$study, rows$subject), , drop = FALSE]
    key <- paste(rows$study, rows$subject, sep = ":")
    s6g_assert(nrow(rows) == as.integer(expected_subjects) &&
                 !anyDuplicated(key),
               paste0("Invalid Stage 6G subject losses for ", config_id, "."))
    if (is.null(reference_key)) {
      reference_key <<- key
    } else {
      s6g_assert(identical(key, reference_key),
                 "Subject keys differ across Stage 6G candidates.")
    }
    loss <- as.numeric(rows[[loss_column]])
    metadata <- config_metadata[
      config_metadata$fit_config_id == config_id, , drop = FALSE
    ]
    data.frame(
      fit_config_id = config_id,
      complexity_proxy = as.integer(metadata$complexity_proxy[[1L]]),
      subjects = nrow(rows), paired_splits_per_subject = 2L,
      mean_loss = s6g_stratified_mean(loss, rows$study),
      stratified_subject_se = s6g_stratified_se(loss, rows$study),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  best_order <- order(summary$mean_loss, summary$complexity_proxy,
                      summary$fit_config_id)
  best_index <- best_order[[1L]]
  threshold <- summary$mean_loss[[best_index]] +
    summary$stratified_subject_se[[best_index]]
  summary$minimum_loss <- seq_len(nrow(summary)) == best_index
  summary$one_se_threshold <- threshold
  summary$within_one_se <- summary$mean_loss <=
    threshold + 10 * .Machine$double.eps
  eligible <- which(summary$within_one_se)
  selected_index <- eligible[order(
    summary$complexity_proxy[eligible], summary$fit_config_id[eligible]
  )][[1L]]
  summary$selected_by_one_se <- seq_len(nrow(summary)) == selected_index
  summary <- summary[order(summary$complexity_proxy,
                           summary$fit_config_id), , drop = FALSE]
  rownames(summary) <- NULL
  list(
    summary = summary,
    selected_config_id = as.character(
      summary$fit_config_id[summary$selected_by_one_se]
    ),
    minimum_loss_config_id = as.character(
      summary$fit_config_id[summary$minimum_loss]
    ),
    threshold = threshold
  )
}

s6g_split_ranking_stability <- function(score_rows, loss_column) {
  required <- c("fit_config_id", "split_id", loss_column)
  s6g_assert(all(required %in% names(score_rows)),
             "Stage 6G split-ranking inputs are incomplete.")
  split_ids <- sort(unique(as.character(score_rows$split_id)))
  config_ids <- sort(unique(as.character(score_rows$fit_config_id)))
  s6g_assert(length(split_ids) == 2L && length(config_ids) == 5L,
             "Stage 6G ranking requires two splits and five candidates.")
  risks <- do.call(rbind, lapply(split_ids, function(split_id) {
    rows <- score_rows[score_rows$split_id == split_id, , drop = FALSE]
    values <- vapply(config_ids, function(config_id) {
      mean(as.numeric(rows[rows$fit_config_id == config_id, loss_column]))
    }, numeric(1L))
    data.frame(
      split_id = split_id, fit_config_id = config_ids,
      mean_loss = values, rank = rank(values, ties.method = "average"),
      stringsAsFactors = FALSE
    )
  }))
  wide <- merge(
    risks[risks$split_id == split_ids[[1L]], c("fit_config_id", "rank")],
    risks[risks$split_id == split_ids[[2L]], c("fit_config_id", "rank")],
    by = "fit_config_id", suffixes = c("_split1", "_split2"), sort = TRUE
  )
  spearman <- suppressWarnings(stats::cor(
    wide$rank_split1, wide$rank_split2, method = "spearman"
  ))
  minimum_by_split <- vapply(split_ids, function(split_id) {
    rows <- risks[risks$split_id == split_id, , drop = FALSE]
    rows$fit_config_id[order(rows$mean_loss, rows$fit_config_id)][[1L]]
  }, character(1L))
  pairs <- utils::combn(config_ids, 2L, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    differences <- vapply(split_ids, function(split_id) {
      rows <- risks[risks$split_id == split_id, , drop = FALSE]
      rows$mean_loss[rows$fit_config_id == pair[[1L]]] -
        rows$mean_loss[rows$fit_config_id == pair[[2L]]]
    }, numeric(1L))
    data.frame(
      config_a = pair[[1L]], config_b = pair[[2L]],
      difference_split1 = differences[[1L]],
      difference_split2 = differences[[2L]],
      any_tie = any(differences == 0),
      ordering_concordant = identical(sign(differences[[1L]]),
                                      sign(differences[[2L]])),
      stringsAsFactors = FALSE
    )
  }))
  summary <- data.frame(
    split_1 = split_ids[[1L]], split_2 = split_ids[[2L]],
    spearman_rank_correlation = spearman,
    minimum_loss_config_split_1 = minimum_by_split[[1L]],
    minimum_loss_config_split_2 = minimum_by_split[[2L]],
    minimum_loss_winner_agreement =
      identical(minimum_by_split[[1L]], minimum_by_split[[2L]]),
    pairwise_comparisons = nrow(pairwise),
    pairwise_ordering_concordance = mean(pairwise$ordering_concordant),
    stringsAsFactors = FALSE
  )
  list(risks = risks, rank_table = wide, pairwise = pairwise,
       summary = summary)
}

s6g_paired_contrasts <- function(subject_losses, selection) {
  reference_id <- selection$minimum_loss_config_id
  reference <- subject_losses[
    subject_losses$fit_config_id == reference_id, , drop = FALSE
  ]
  reference <- reference[order(reference$study, reference$subject), ]
  reference_key <- paste(reference$study, reference$subject, sep = ":")
  config_ids <- sort(unique(as.character(subject_losses$fit_config_id)))
  do.call(rbind, lapply(config_ids, function(config_id) {
    rows <- subject_losses[
      subject_losses$fit_config_id == config_id, , drop = FALSE
    ]
    rows <- rows[order(rows$study, rows$subject), ]
    key <- paste(rows$study, rows$subject, sep = ":")
    s6g_assert(identical(key, reference_key),
               "Stage 6G paired-contrast subject keys differ.")
    difference <- as.numeric(rows$mse) - as.numeric(reference$mse)
    full_mean <- s6g_stratified_mean(difference, rows$study)
    leave_one <- vapply(seq_along(difference), function(index) {
      s6g_stratified_mean(difference[-index], rows$study[-index])
    }, numeric(1L))
    influence <- abs(leave_one - full_mean)
    influential <- which.max(influence)
    data.frame(
      fit_config_id = config_id,
      reference_minimum_loss_config_id = reference_id,
      mean_paired_mse_difference = full_mean,
      paired_difference_se = s6g_stratified_se(difference, rows$study),
      maximum_leave_one_subject_influence = influence[[influential]],
      most_influential_study = rows$study[[influential]],
      most_influential_subject = rows$subject[[influential]],
      stringsAsFactors = FALSE
    )
  }))
}

s6g_leave_one_subject_out <- function(subject_losses, configs) {
  keys <- unique(subject_losses[c("study", "subject")])
  keys <- keys[order(keys$study, keys$subject), , drop = FALSE]
  decisions <- lapply(seq_len(nrow(keys)), function(index) {
    keep <- !(subject_losses$study == keys$study[[index]] &
                subject_losses$subject == keys$subject[[index]])
    selection <- s6g_select_one_se(
      subject_losses[keep, , drop = FALSE], configs, "mse",
      expected_subjects = nrow(keys) - 1L
    )
    data.frame(
      omitted_study = keys$study[[index]],
      omitted_subject = keys$subject[[index]],
      selected_config_id = selection$selected_config_id,
      minimum_loss_config_id = selection$minimum_loss_config_id,
      one_se_threshold = selection$threshold,
      stringsAsFactors = FALSE
    )
  })
  decisions <- do.call(rbind, decisions)
  count <- as.data.frame(table(
    factor(decisions$selected_config_id,
           levels = as.character(configs$fit_config_id))
  ), stringsAsFactors = FALSE)
  names(count) <- c("fit_config_id", "leave_one_subject_out_selections")
  count$leave_one_subject_out_selections <-
    as.integer(count$leave_one_subject_out_selections)
  list(decisions = decisions, counts = count)
}
