s6f_abort <- function(...) stop(paste0(...), call. = FALSE)
s6f_assert <- function(value, message) if (!isTRUE(value)) s6f_abort(message)

s6f_read_csv <- function(path) {
  s6f_assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

s6f_validate_registration <- function(configs, seeds, splits) {
  expected_configs <- c(
    "config_index", "fit_config_id", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2", "L_common",
    "M_common", "complexity_proxy", "stage6b_source_config_id"
  )
  expected_seeds <- c(
    "seed_index", "fit_seed", "paired_across_splits_and_dimensions"
  )
  expected_splits <- c(
    "split_index", "split_id", "holdout_rule", "reuses_stage6b_split"
  )
  s6f_assert(identical(names(configs), expected_configs),
             "Stage 6F configuration columns differ from registration.")
  s6f_assert(identical(names(seeds), expected_seeds),
             "Stage 6F seed columns differ from registration.")
  s6f_assert(identical(names(splits), expected_splits),
             "Stage 6F split columns differ from registration.")
  s6f_assert(
    nrow(configs) == 4L && !anyDuplicated(configs$fit_config_id) &&
      identical(as.integer(configs$config_index), 1:4) &&
      identical(as.character(configs$fit_config_id), c(
        "joint_L1_M2", "joint_L1_M4", "joint_L2_M2", "joint_L2_M4"
      )),
    "Stage 6F must contain the registered joint 2 by 2 dimension grid."
  )
  s6f_assert(
    all(as.integer(configs$fit_L_f) == as.integer(configs$L_common)) &&
      all(as.integer(configs$fit_L_s_1) == as.integer(configs$L_common)) &&
      all(as.integer(configs$fit_L_s_2) == as.integer(configs$L_common)),
    "Stage 6F L values are inconsistent."
  )
  expected_complexity <- 1500L * as.integer(configs$L_common) +
    141L * as.integer(configs$L_common) * as.integer(configs$M_common)
  s6f_assert(identical(as.integer(configs$complexity_proxy),
                       expected_complexity),
             "Stage 6F complexity proxy differs from registration.")
  for (index in seq_len(nrow(configs))) {
    L <- as.integer(configs$L_common[[index]])
    M <- as.integer(configs$M_common[[index]])
    expected_M <- paste(rep(M, L), collapse = ";")
    registered_M <- unlist(configs[index, c(
      "fit_M_f", "fit_M_s_1", "fit_M_s_2"
    )], use.names = FALSE)
    s6f_assert(all(as.character(registered_M) == expected_M), paste0(
      "Stage 6F M vector differs from registration for ",
      configs$fit_config_id[[index]], "."
    ))
  }
  s6f_assert(
    nrow(seeds) == 4L && !anyDuplicated(seeds$fit_seed) &&
      identical(as.integer(seeds$seed_index), 1:4) &&
      identical(as.integer(seeds$fit_seed), 83206101:83206104) &&
      all(seeds$paired_across_splits_and_dimensions),
    "Stage 6F must use the four paired Stage 6B start seeds."
  )
  s6f_assert(
    nrow(splits) == 2L && !anyDuplicated(splits$split_id) &&
      identical(as.integer(splits$split_index), 1:2) &&
      identical(as.character(splits$split_id),
                c("split1_middle", "split2_upper")) &&
      identical(as.logical(splits$reuses_stage6b_split), c(TRUE, FALSE)),
    "Stage 6F split registration is invalid."
  )
  invisible(TRUE)
}

s6f_build_manifest <- function(configs, seeds, splits,
                                data_id = "g12s6a_01") {
  s6f_validate_registration(configs, seeds, splits)
  grid <- merge(splits, configs, by = NULL, sort = FALSE)
  grid <- merge(grid, seeds, by = NULL, sort = FALSE)
  grid <- grid[order(grid$split_index, grid$config_index, grid$seed_index),
               , drop = FALSE]
  grid$experiment_id <- "G12_STAGE6F_JOINT_DIMENSION_SELECTION_V1_20260911"
  grid$data_id <- data_id
  grid$scenario_id <- "baseline_strong_stage6a_development_reuse"
  grid$method_id <- "G12"
  grid$selection_stratum_id <- paste(
    data_id, grid$split_id, grid$fit_config_id, sep = "__"
  )
  grid$fit_id <- paste0(
    data_id, "__s6f__", grid$split_id, "__", grid$fit_config_id, "__",
    sprintf("%02d", as.integer(grid$seed_index))
  )
  grid$route_id <- paste0(
    "random__gram_unit_energy__pre1__jaoua_default__multistart4__",
    "fixed400__", grid$split_id, "__", grid$fit_config_id
  )
  grid$reuse_stage6b <- grid$split_id == "split1_middle" &
    nzchar(as.character(grid$stage6b_source_config_id))
  grid$execute_new_fit <- !grid$reuse_stage6b
  grid$source_stage6b_fit_id <- ifelse(
    grid$reuse_stage6b,
    paste0(
      data_id, "__s6b__", grid$stage6b_source_config_id, "__",
      sprintf("%02d", as.integer(grid$seed_index))
    ),
    ""
  )
  grid$initialization <- "random"
  grid$function_initialization <- "gram_unit_energy"
  grid$pre_score_sweeps <- 1L
  grid$anneal <- "c(1,1.9,100)"
  grid$planned_ordinary_t1_sweeps <- 400L
  grid$total_maxit <- 499L
  grid$n_cpus <- 1L
  grid$objective_eligibility_required <- TRUE
  grid$truth_available_to_fit <- FALSE
  grid$truth_available_to_stopping <- FALSE
  grid$truth_available_to_start_selection <- FALSE
  grid$truth_available_to_dimension_selection <- FALSE
  grid$holdout_available_to_start_selection <- FALSE
  grid$continuation <- FALSE
  grid$development_only <- TRUE
  grid$formal_paper_mc_result <- FALSE
  preferred <- c(
    "experiment_id", "fit_id", "data_id", "scenario_id", "method_id",
    "split_index", "split_id", "holdout_rule", "fit_config_id",
    "config_index", "selection_stratum_id", "route_id", "seed_index",
    "fit_seed", "fit_L_f", "fit_L_s_1", "fit_L_s_2", "fit_M_f",
    "fit_M_s_1", "fit_M_s_2", "L_common", "M_common",
    "complexity_proxy", "reuse_stage6b", "execute_new_fit",
    "stage6b_source_config_id", "source_stage6b_fit_id", "initialization",
    "function_initialization", "pre_score_sweeps", "anneal",
    "planned_ordinary_t1_sweeps", "total_maxit", "n_cpus",
    "objective_eligibility_required", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_start_selection",
    "truth_available_to_dimension_selection",
    "holdout_available_to_start_selection", "continuation",
    "development_only", "formal_paper_mc_result"
  )
  missing_columns <- setdiff(preferred, names(grid))
  s6f_assert(!length(missing_columns), paste0(
    "Stage 6F manifest construction lacks columns: ",
    paste(missing_columns, collapse = ";")
  ))
  result <- grid[preferred]
  rownames(result) <- NULL
  s6f_assert(
    nrow(result) == 32L && !anyDuplicated(result$fit_id) &&
      all(table(result$selection_stratum_id) == 4L) &&
      sum(result$reuse_stage6b) == 12L &&
      sum(result$execute_new_fit) == 20L &&
      all(result$n_cpus == 1L) &&
      !any(result$truth_available_to_fit) &&
      !any(result$truth_available_to_start_selection) &&
      !any(result$truth_available_to_dimension_selection),
    "Stage 6F manifest must contain 12 reused and 20 new endpoints."
  )
  result
}

s6f_make_time_holdout <- function(observation,
                                   rule = c("middle", "upper"),
                                   min_train_times = 5L) {
  rule <- match.arg(rule)
  required <- c("Y", "time_obs", "Z", "dimensions")
  s6f_assert(all(required %in% names(observation)),
             "Observation bundle is incomplete.")
  if (exists("s6b_forbidden_observation_names", mode = "function")) {
    s6f_assert(!length(s6b_forbidden_observation_names(observation)),
               "Observation bundle contains a forbidden truth field.")
  }
  s6f_assert(is.null(observation$Z) &&
               as.integer(observation$dimensions$d) == 0L,
             "Stage 6F v1 supports the registered d=0 route only.")
  train <- observation
  records <- list()
  record_index <- 0L
  for (study in seq_along(observation$Y)) {
    s6f_assert(length(observation$Y[[study]]) ==
                 length(observation$time_obs[[study]]),
               "Subject counts differ between Y and time_obs.")
    for (subject in seq_along(observation$Y[[study]])) {
      time <- as.numeric(observation$time_obs[[study]][[subject]])
      n_time <- length(time)
      s6f_assert(n_time - 1L >= as.integer(min_train_times),
                 "Holdout would leave too few training times.")
      middle <- as.integer(floor((n_time + 1L) / 2L))
      holdout_index <- if (identical(rule, "middle")) middle else middle + 1L
      s6f_assert(holdout_index > 1L && holdout_index < n_time,
                 "Registered holdout is not an interior time.")
      responses <- observation$Y[[study]][[subject]]
      p <- length(responses)
      s6f_assert(p >= 1L &&
                   all(vapply(responses, length, integer(1L)) == n_time),
                 "Variables do not share the registered subject time grid.")
      heldout_y <- vapply(responses, function(value) {
        as.numeric(value[[holdout_index]])
      }, numeric(1L))
      train$Y[[study]][[subject]] <- lapply(responses, function(value) {
        as.numeric(value[-holdout_index])
      })
      train$time_obs[[study]][[subject]] <- time[-holdout_index]
      record_index <- record_index + 1L
      records[[record_index]] <- list(
        study = as.integer(study), subject = as.integer(subject),
        original_index = holdout_index, time = time[[holdout_index]],
        y = heldout_y, p = p, train_time_count = n_time - 1L
      )
    }
  }
  split_id <- if (identical(rule, "middle")) {
    "middle_interior_whole_time_v1"
  } else {
    "upper_adjacent_interior_whole_time_v1"
  }
  list(
    train_observation = train,
    holdout_plan = list(
      split_id = split_id, records = records,
      subject_time_count = length(records),
      complete_variable_vector_held_out = TRUE, truth_used = FALSE
    )
  )
}

s6f_subject_cluster_losses <- function(score_rows, loss_columns,
                                       expected_splits = 2L) {
  required <- c("fit_config_id", "split_id", "study", "subject",
                loss_columns)
  s6f_assert(all(required %in% names(score_rows)),
             "Stage 6F score rows are incomplete.")
  key <- interaction(
    score_rows$fit_config_id, score_rows$study, score_rows$subject,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(score_rows, key)
  result <- lapply(groups, function(rows) {
    s6f_assert(nrow(rows) == expected_splits &&
                 !anyDuplicated(rows$split_id),
               "A subject/configuration lacks the registered paired splits.")
    base <- data.frame(
      fit_config_id = as.character(rows$fit_config_id[[1L]]),
      study = as.integer(rows$study[[1L]]),
      subject = as.integer(rows$subject[[1L]]),
      paired_splits = nrow(rows), stringsAsFactors = FALSE
    )
    for (column in loss_columns) {
      value <- as.numeric(rows[[column]])
      s6f_assert(all(is.finite(value)),
                 paste0("Non-finite Stage 6F loss: ", column))
      base[[column]] <- mean(value)
    }
    base
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result[order(result$fit_config_id, result$study, result$subject),
         , drop = FALSE]
}

s6f_stratified_se <- function(loss, study) {
  loss <- as.numeric(loss)
  study <- as.integer(study)
  s6f_assert(length(loss) == length(study) && length(loss) >= 4L &&
               all(is.finite(loss)) && all(is.finite(study)),
             "Invalid Stage 6F subject-clustered losses.")
  strata <- split(loss, study)
  sizes <- vapply(strata, length, integer(1L))
  s6f_assert(all(sizes >= 2L),
             "Each study needs at least two subject-clustered losses.")
  weights <- sizes / sum(sizes)
  variances <- vapply(strata, stats::var, numeric(1L))
  sqrt(sum(weights^2 * variances / sizes))
}

s6f_select_one_se <- function(subject_losses, config_metadata,
                              loss_column, criterion_id) {
  required <- c("fit_config_id", "study", "subject", loss_column)
  s6f_assert(all(required %in% names(subject_losses)) &&
               all(c("fit_config_id", "complexity_proxy") %in%
                     names(config_metadata)),
             "Stage 6F one-SE inputs are incomplete.")
  config_ids <- as.character(config_metadata$fit_config_id)
  s6f_assert(length(config_ids) == 4L && !anyDuplicated(config_ids),
             "Stage 6F one-SE selection requires four configurations.")
  reference_key <- NULL
  summaries <- lapply(config_ids, function(config_id) {
    rows <- subject_losses[
      subject_losses$fit_config_id == config_id, , drop = FALSE
    ]
    rows <- rows[order(rows$study, rows$subject), , drop = FALSE]
    key <- paste(rows$study, rows$subject, sep = ":")
    s6f_assert(nrow(rows) == 60L && !anyDuplicated(key),
               paste0("Invalid subject losses for ", config_id, "."))
    if (is.null(reference_key)) reference_key <<- key else
      s6f_assert(identical(key, reference_key),
                 "Subject keys differ across dimension configurations.")
    loss <- as.numeric(rows[[loss_column]])
    metadata <- config_metadata[
      config_metadata$fit_config_id == config_id, , drop = FALSE
    ]
    data.frame(
      criterion_id = criterion_id,
      fit_config_id = config_id,
      complexity_proxy = as.integer(metadata$complexity_proxy[[1L]]),
      subjects = nrow(rows), paired_splits_per_subject = 2L,
      mean_loss = mean(loss),
      stratified_subject_se = s6f_stratified_se(loss, rows$study),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  best_order <- order(summary$mean_loss, summary$complexity_proxy,
                      summary$fit_config_id)
  best_index <- best_order[[1L]]
  threshold <- summary$mean_loss[[best_index]] +
    summary$stratified_subject_se[[best_index]]
  summary$minimum_loss <- FALSE
  summary$minimum_loss[[best_index]] <- TRUE
  summary$one_se_threshold <- threshold
  summary$within_one_se <- summary$mean_loss <=
    threshold + 10 * .Machine$double.eps
  eligible <- which(summary$within_one_se)
  selected_index <- eligible[order(
    summary$complexity_proxy[eligible], summary$fit_config_id[eligible]
  )][[1L]]
  summary$selected_by_one_se <- FALSE
  summary$selected_by_one_se[[selected_index]] <- TRUE
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

s6f_split_ranking_stability <- function(score_rows, loss_column,
                                        criterion_id) {
  required <- c("fit_config_id", "split_id", loss_column)
  s6f_assert(all(required %in% names(score_rows)),
             "Stage 6F split-ranking inputs are incomplete.")
  split_ids <- sort(unique(as.character(score_rows$split_id)))
  config_ids <- sort(unique(as.character(score_rows$fit_config_id)))
  s6f_assert(length(split_ids) == 2L && length(config_ids) == 4L,
             "Stage 6F ranking requires two splits and four configurations.")
  risks <- do.call(rbind, lapply(split_ids, function(split_id) {
    rows <- score_rows[score_rows$split_id == split_id, , drop = FALSE]
    values <- vapply(config_ids, function(config_id) {
      mean(as.numeric(rows[rows$fit_config_id == config_id, loss_column]))
    }, numeric(1L))
    data.frame(
      criterion_id = criterion_id, split_id = split_id,
      fit_config_id = config_ids, mean_loss = values,
      rank = rank(values, ties.method = "average"),
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
      criterion_id = criterion_id,
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
    criterion_id = criterion_id,
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
