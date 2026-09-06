s6d_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing Stage 6D CSV: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

s6d_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

s6d_set_warning_phase <- function(warnings, phase) {
  if (is.data.frame(warnings) && "phase" %in% names(warnings)) {
    warnings$phase <- rep(phase, nrow(warnings))
  }
  warnings
}

s6d_validate_registration <- function(data, configs, seeds) {
  data_columns <- c(
    "data_index", "data_id", "data_seed", "scenario_id", "S", "n_s_1",
    "n_s_2", "p", "d", "true_L_f", "true_L_s_1", "true_L_s_2",
    "true_M", "K", "n_obs_min", "n_obs_max", "sigma_eps", "mean_amp",
    "truth_unseal_authorized"
  )
  config_columns <- c(
    "config_index", "fit_config_id", "fit_L_f", "fit_L_s_1",
    "fit_L_s_2", "fit_M_f", "fit_M_s_1", "fit_M_s_2",
    "factor_screen", "fpca_screen"
  )
  seed_columns <- c(
    "seed_index", "fit_seed", "paired_across_data_and_dimension_configs"
  )
  s6d_assert(identical(names(data), data_columns),
             "Stage 6D data columns differ from the frozen registration.")
  s6d_assert(identical(names(configs), config_columns),
             "Stage 6D config columns differ from the frozen registration.")
  s6d_assert(identical(names(seeds), seed_columns),
             "Stage 6D seed columns differ from the frozen registration.")

  s6d_assert(nrow(data) == 3L &&
               identical(as.integer(data$data_index), 1:3) &&
               identical(as.character(data$data_id),
                         sprintf("g12s6d_%02d", 1:3)) &&
               identical(as.integer(data$data_seed), 83606001:83606003) &&
               !anyDuplicated(data$data_id) && !anyDuplicated(data$data_seed),
             "Stage 6D must contain exactly three pre-registered new data seeds.")
  s6d_assert(all(data$S == 2L & data$n_s_1 == 30L & data$n_s_2 == 30L &
                   data$p == 500L & data$d == 0L & data$true_L_f == 1L &
                   data$true_L_s_1 == 1L & data$true_L_s_2 == 1L &
                   data$true_M == 2L & data$K == 5L &
                   data$n_obs_min == 6L & data$n_obs_max == 9L &
                   abs(data$sigma_eps - 0.3) < 1e-15 &
                   abs(data$mean_amp - 0.6) < 1e-15) &&
               !any(data$truth_unseal_authorized),
             "Stage 6D data-generating settings or truth lock changed.")

  expected_ids <- c(
    "factor_L1_M2", "factor_L2_M2", "factor_L3_M2",
    "fpca_L1_M3", "fpca_L1_M4"
  )
  s6d_assert(nrow(configs) == 5L &&
               identical(as.integer(configs$config_index), 1:5) &&
               identical(as.character(configs$fit_config_id), expected_ids) &&
               !anyDuplicated(configs$fit_config_id),
             "Stage 6D must contain exactly the five registered configurations.")
  factor <- configs[configs$factor_screen, , drop = FALSE]
  fpca <- configs[configs$fpca_screen, , drop = FALSE]
  s6d_assert(nrow(factor) == 3L &&
               identical(as.integer(factor$fit_L_f), 1:3) &&
               all(factor$fit_L_f == factor$fit_L_s_1) &&
               all(factor$fit_L_f == factor$fit_L_s_2) &&
               all(vapply(strsplit(factor$fit_M_f, ";", fixed = TRUE),
                          function(value) all(as.integer(value) == 2L),
                          logical(1L))),
             "The factor axis must be L=1:3 at M=2.")
  s6d_assert(nrow(fpca) == 3L &&
               identical(as.character(fpca$fit_config_id),
                         c("factor_L1_M2", "fpca_L1_M3", "fpca_L1_M4")) &&
               all(fpca$fit_L_f == 1L & fpca$fit_L_s_1 == 1L &
                     fpca$fit_L_s_2 == 1L) &&
               identical(as.integer(fpca$fit_M_f), 2:4),
             "The FPCA axis must be M=2:4 at L=1.")

  s6d_assert(nrow(seeds) == 12L &&
               identical(as.integer(seeds$seed_index), 1:12) &&
               identical(as.integer(seeds$fit_seed), 83606101:83606112) &&
               !anyDuplicated(seeds$fit_seed) &&
               all(seeds$paired_across_data_and_dimension_configs),
             "Stage 6D must contain twelve unique paired fit seeds.")
  invisible(TRUE)
}

s6d_build_fit_manifest <- function(data, configs, seeds) {
  s6d_validate_registration(data, configs, seeds)
  grid <- merge(merge(data, configs, by = NULL, sort = FALSE), seeds,
                by = NULL, sort = FALSE)
  grid <- grid[order(grid$data_index, grid$config_index, grid$seed_index),
               , drop = FALSE]
  grid$experiment_id <- "G12_STAGE6D_INDEPENDENT_DIMENSION_SELECTION_V1_20260906"
  grid$method_id <- "G12"
  grid$selection_stratum_id <- paste(grid$data_id, grid$fit_config_id,
                                     sep = "__")
  grid$route_id <- paste0(
    "random__gram_unit_energy__pre1__jaoua_default__multistart12__",
    "fixed400__middle_time_holdout__", grid$fit_config_id
  )
  grid$fit_id <- paste0(
    grid$data_id, "__s6d__", grid$fit_config_id, "__",
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
    "experiment_id", "fit_id", "data_id", "data_index", "data_seed",
    "scenario_id", "method_id", "fit_config_id", "config_index",
    "selection_stratum_id", "route_id", "seed_index", "fit_seed",
    "fit_L_f", "fit_L_s_1", "fit_L_s_2", "fit_M_f", "fit_M_s_1",
    "fit_M_s_2", "factor_screen", "fpca_screen", "initialization",
    "function_initialization", "pre_score_sweeps", "anneal",
    "planned_ordinary_t1_sweeps", "total_maxit", "n_cpus",
    "objective_eligibility_required", "truth_available_to_fit",
    "truth_available_to_stopping", "truth_available_to_start_selection",
    "truth_available_to_dimension_selection", "continuation",
    "development_only", "formal_paper_mc_result"
  )
  result <- grid[preferred]
  rownames(result) <- NULL
  s6d_assert(nrow(result) == 180L && !anyDuplicated(result$fit_id) &&
               all(table(result$data_id) == 60L) &&
               all(table(result$selection_stratum_id) == 12L) &&
               all(result$n_cpus == 1L) &&
               !any(result$truth_available_to_fit) &&
               !any(result$truth_available_to_stopping) &&
               !any(result$truth_available_to_start_selection) &&
               !any(result$truth_available_to_dimension_selection),
             "Stage 6D expanded manifest must contain 180 truth-isolated fits.")
  result
}

s6d_axis_metadata <- function(configs, axis = c("factor", "fpca")) {
  axis <- match.arg(axis)
  if (axis == "factor") {
    result <- configs[configs$factor_screen,
                      c("fit_config_id", "fit_L_f"), drop = FALSE]
    names(result)[[2L]] <- "complexity_value"
  } else {
    result <- configs[configs$fpca_screen,
                      c("fit_config_id", "fit_M_f"), drop = FALSE]
    result$complexity_value <- as.integer(result$fit_M_f)
    result$fit_M_f <- NULL
  }
  result$complexity_value <- as.integer(result$complexity_value)
  s6d_assert(nrow(result) == 3L &&
               identical(result$complexity_value, 1:3 + (axis == "fpca")),
             paste0("Invalid Stage 6D ", axis, "-axis metadata."))
  result
}

s6d_select_axes <- function(winner_losses, configs, data_id) {
  rows <- winner_losses[winner_losses$data_id == data_id, , drop = FALSE]
  s6d_assert(length(unique(rows$fit_config_id)) == 5L,
             paste0("Winner losses are incomplete for ", data_id, "."))
  factor <- s6bp_select_one_se(
    rows[rows$fit_config_id %in%
           s6d_axis_metadata(configs, "factor")$fit_config_id, , drop = FALSE],
    s6d_axis_metadata(configs, "factor"), "factor_screen"
  )
  fpca <- s6bp_select_one_se(
    rows[rows$fit_config_id %in%
           s6d_axis_metadata(configs, "fpca")$fit_config_id, , drop = FALSE],
    s6d_axis_metadata(configs, "fpca"), "fpca_screen"
  )
  list(factor = factor, fpca = fpca)
}
