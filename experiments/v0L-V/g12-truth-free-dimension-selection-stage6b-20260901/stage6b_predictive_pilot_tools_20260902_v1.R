s6bp_read_csv <- function(path) {
  s6b_assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

s6bp_set_warning_phase <- function(warnings, phase) {
  if (!is.data.frame(warnings) || !"phase" %in% names(warnings)) {
    return(warnings)
  }
  warnings$phase <- rep(phase, nrow(warnings))
  warnings
}

s6bp_validate_registration <- function(configs, seeds) {
  config_columns <- c(
    "config_index", "phase", "fit_config_id", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2",
    "complexity_value", "run_condition"
  )
  seed_columns <- c("seed_index", "fit_seed",
                    "paired_across_dimension_configs")
  s6b_assert(identical(names(configs), config_columns),
             "Stage 6B fit-config columns differ from registration.")
  s6b_assert(identical(names(seeds), seed_columns),
             "Stage 6B start-seed columns differ from registration.")
  s6b_assert(nrow(configs) == 9L && !anyDuplicated(configs$fit_config_id) &&
               identical(as.integer(configs$config_index), seq_len(9L)),
             "Stage 6B must register exactly nine unique configurations.")
  s6b_assert(nrow(seeds) == 4L && !anyDuplicated(seeds$fit_seed) &&
               identical(as.integer(seeds$seed_index), seq_len(4L)) &&
               all(seeds$paired_across_dimension_configs),
             "Stage 6B must register four unique paired start seeds.")

  factor <- configs[configs$phase == "factor_screen", , drop = FALSE]
  fpca <- configs[configs$phase == "fpca_screen", , drop = FALSE]
  s6b_assert(nrow(factor) == 3L && nrow(fpca) == 6L,
             "Stage 6B phase configuration counts must be 3 and 6.")
  s6b_assert(identical(as.integer(factor$fit_L_f), 1:3) &&
               identical(as.integer(factor$fit_L_s_1), 1:3) &&
               identical(as.integer(factor$fit_L_s_2), 1:3) &&
               all(vapply(strsplit(as.character(factor$fit_M_f), ";", TRUE),
                          function(value) all(as.integer(value) == 2L),
                          logical(1L))),
             "The factor screen must be common L=1:3 at M=2.")
  s6b_assert(all(as.integer(fpca$fit_L_f) %in% 1:3) &&
               all(as.integer(fpca$complexity_value) %in% 3:4) &&
               all(table(fpca$fit_L_f) == 2L),
             "The conditional FPCA screen must register M=3/4 for each L.")
  invisible(TRUE)
}

s6bp_build_fit_manifest <- function(configs, seeds,
                                     data_id = "g12s6a_01") {
  s6bp_validate_registration(configs, seeds)
  grid <- merge(configs, seeds, by = NULL, sort = FALSE)
  grid <- grid[order(grid$config_index, grid$seed_index), , drop = FALSE]
  grid$experiment_id <- "G12_STAGE6B_PREDICTIVE_PILOT_V1_20260902"
  grid$data_id <- data_id
  grid$scenario_id <- "baseline_strong_stage6a_development_reuse"
  grid$method_id <- "G12"
  grid$selection_stratum_id <- paste(data_id, grid$fit_config_id, sep = "__")
  grid$route_id <- paste0(
    "random__gram_unit_energy__pre1__jaoua_default__multistart4__",
    "fixed400__middle_time_holdout__", grid$fit_config_id
  )
  grid$fit_id <- paste0(
    data_id, "__s6b__", grid$fit_config_id, "__",
    sprintf("%02d", as.integer(grid$seed_index))
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
  grid$continuation <- FALSE
  grid$development_only <- TRUE
  grid$formal_paper_mc_result <- FALSE
  preferred <- c(
    "experiment_id", "fit_id", "data_id", "scenario_id", "method_id",
    "phase", "fit_config_id", "config_index", "selection_stratum_id",
    "route_id", "seed_index", "fit_seed", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2",
    "complexity_value", "run_condition", "initialization",
    "function_initialization", "pre_score_sweeps", "anneal",
    "planned_ordinary_t1_sweeps", "total_maxit", "n_cpus",
    "objective_eligibility_required", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_start_selection",
    "truth_available_to_dimension_selection", "continuation",
    "development_only", "formal_paper_mc_result"
  )
  result <- grid[preferred]
  rownames(result) <- NULL
  s6b_assert(nrow(result) == 36L && !anyDuplicated(result$fit_id) &&
               all(table(result$fit_config_id) == 4L),
             "Stage 6B expanded fit manifest must contain 36 conditional rows.")
  result
}

s6bp_phase_manifest <- function(manifest, phase,
                                selected_factor_count = NULL) {
  if (identical(phase, "factor_screen")) {
    result <- manifest[manifest$phase == phase, , drop = FALSE]
    s6b_assert(nrow(result) == 12L,
               "The Stage 6B factor phase must contain 12 fits.")
    return(result)
  }
  s6b_assert(identical(phase, "fpca_screen") &&
               length(selected_factor_count) == 1L &&
               as.integer(selected_factor_count) %in% 1:3,
             "The FPCA phase requires a selected factor count in 1:3.")
  result <- manifest[
    manifest$phase == phase &
      as.integer(manifest$fit_L_f) == as.integer(selected_factor_count),
    , drop = FALSE
  ]
  s6b_assert(nrow(result) == 8L &&
               length(unique(result$fit_config_id)) == 2L,
             "The conditional Stage 6B FPCA phase must contain eight fits.")
  result
}

s6bp_stratified_mean_se <- function(loss, study) {
  loss <- as.numeric(loss)
  study <- as.integer(study)
  s6b_assert(length(loss) == length(study) && length(loss) >= 4L &&
               all(is.finite(loss)) && all(is.finite(study)),
             "Invalid subject-time losses for stratified standard error.")
  split_loss <- split(loss, study)
  sizes <- vapply(split_loss, length, integer(1L))
  s6b_assert(all(sizes >= 2L),
             "Every study stratum needs at least two held-out subjects.")
  weights <- sizes / sum(sizes)
  variances <- vapply(split_loss, stats::var, numeric(1L))
  sqrt(sum(weights^2 * variances / sizes))
}

s6bp_select_one_se <- function(subject_losses, config_metadata,
                               phase = c("factor_screen", "fpca_screen")) {
  phase <- match.arg(phase)
  required_losses <- c("fit_config_id", "study", "subject", "mse", "mae",
                       "observed_second_moment")
  required_configs <- c("fit_config_id", "complexity_value")
  s6b_assert(all(required_losses %in% names(subject_losses)) &&
               all(required_configs %in% names(config_metadata)),
             "Stage 6B one-SE inputs are incomplete.")
  config_ids <- unique(as.character(config_metadata$fit_config_id))
  s6b_assert(length(config_ids) == 3L && !anyDuplicated(config_ids) &&
               all(config_ids %in% subject_losses$fit_config_id),
             "Stage 6B one-SE selection requires exactly three configurations.")

  key_reference <- NULL
  summaries <- lapply(config_ids, function(config_id) {
    rows <- subject_losses[subject_losses$fit_config_id == config_id,
                           , drop = FALSE]
    rows <- rows[order(rows$study, rows$subject), , drop = FALSE]
    key <- paste(rows$study, rows$subject, sep = ":")
    s6b_assert(nrow(rows) > 0L && !anyDuplicated(key),
               paste0("Invalid held-out units for ", config_id, "."))
    if (is.null(key_reference)) key_reference <<- key else
      s6b_assert(identical(key, key_reference),
                 "Held-out subject-time keys differ across configurations.")
    mean_mse <- mean(rows$mse)
    denominator <- sqrt(mean(rows$observed_second_moment))
    metadata <- config_metadata[
      config_metadata$fit_config_id == config_id, , drop = FALSE
    ]
    s6b_assert(nrow(metadata) == 1L,
               paste0("Configuration metadata is not unique: ", config_id))
    data.frame(
      phase = phase,
      fit_config_id = config_id,
      complexity_value = as.integer(metadata$complexity_value[[1L]]),
      heldout_units = nrow(rows),
      mean_subject_time_mse = mean_mse,
      stratified_se_mse = s6bp_stratified_mean_se(rows$mse, rows$study),
      heldout_rmse = sqrt(mean_mse),
      heldout_nrmse = if (denominator > 0) sqrt(mean_mse) / denominator else NA_real_,
      heldout_mae = mean(rows$mae),
      observed_rms = denominator,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  s6b_assert(max(summary$observed_rms) - min(summary$observed_rms) < 1e-12,
             "Held-out response scale differs across dimension configurations.")
  best_order <- order(summary$mean_subject_time_mse,
                      summary$complexity_value, summary$fit_config_id)
  best_index <- best_order[[1L]]
  threshold <- summary$mean_subject_time_mse[[best_index]] +
    summary$stratified_se_mse[[best_index]]
  summary$minimum_loss <- FALSE
  summary$minimum_loss[[best_index]] <- TRUE
  summary$one_se_threshold_mse <- threshold
  summary$within_one_se <- summary$mean_subject_time_mse <=
    threshold + 10 * .Machine$double.eps
  eligible <- which(summary$within_one_se)
  selected_order <- eligible[order(summary$complexity_value[eligible],
                                    summary$fit_config_id[eligible])]
  selected_index <- selected_order[[1L]]
  summary$selected_by_one_se <- FALSE
  summary$selected_by_one_se[[selected_index]] <- TRUE
  summary <- summary[order(summary$complexity_value, summary$fit_config_id),
                     , drop = FALSE]
  rownames(summary) <- NULL
  list(
    summary = summary,
    selected_config_id = as.character(
      summary$fit_config_id[summary$selected_by_one_se]
    ),
    minimum_loss_config_id = as.character(
      summary$fit_config_id[summary$minimum_loss]
    ),
    threshold_mse = threshold
  )
}
