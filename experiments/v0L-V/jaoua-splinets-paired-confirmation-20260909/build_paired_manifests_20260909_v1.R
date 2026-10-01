#!/usr/bin/env Rscript

arguments <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(arguments) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(sub("^--file=", "", arguments[[1L]]))
if (!dir.exists(root)) {
  root <- dirname(normalizePath(sub("^--file=", "", arguments[[1L]]),
                                winslash = "/", mustWork = TRUE))
}

experiment_id <- "JAOUA_SPLINETS_PAIRED_CONFIRMATION_V1_20260909"
logical <- data.frame(
  logical_data_id = sprintf("jspc_%02d", 1:3),
  logical_data_index = 1:3,
  data_seed = 9091001:9091003,
  truth_basis_seed = 9091101:9091103,
  stringsAsFactors = FALSE
)
arms <- data.frame(
  basis_arm = c("legacy", "jaoua"),
  truth_basis_method = c("legacy_dense_gram", "jaoua_splinets"),
  stringsAsFactors = FALSE
)
data_manifest <- merge(logical, arms, by = NULL, sort = FALSE)
data_manifest <- data_manifest[order(
  data_manifest$logical_data_index,
  match(data_manifest$basis_arm, c("legacy", "jaoua"))
), ]
data_manifest$data_id <- paste(data_manifest$logical_data_id,
                               data_manifest$basis_arm, sep = "__")
data_manifest$experiment_id <- experiment_id
data_manifest$scenario_id <- "baseline_strong"
data_manifest$S <- 2L
data_manifest$n_s <- "30;30"
data_manifest$p <- 500L
data_manifest$d <- 0L
data_manifest$truth_L_f <- 1L
data_manifest$truth_L_s <- "1;1"
data_manifest$truth_M_f <- "2"
data_manifest$truth_M_s <- "2;2"
data_manifest$K <- 5L
data_manifest$n_obs_min <- 6L
data_manifest$n_obs_max <- 9L
data_manifest$sigma_eps <- 0.3
data_manifest$mean_amp <- 0.6
data_manifest$paired_nonbasis_rng_required <- TRUE
data_manifest$new_unseen_at_registration <- TRUE
data_manifest$truth_available_to_fit <- FALSE
data_manifest$truth_available_to_stopping <- FALSE
data_manifest$truth_available_to_selection <- FALSE
data_manifest$truth_unseal_authorized <- FALSE
data_manifest$development_only <- TRUE
data_manifest$formal_v0lv_result <- FALSE
data_manifest$formal_paper_mc_result <- FALSE
data_manifest <- data_manifest[c(
  "experiment_id", "data_id", "logical_data_id", "logical_data_index",
  "basis_arm", "truth_basis_method", "data_seed", "truth_basis_seed",
  "scenario_id", "S", "n_s", "p", "d", "truth_L_f", "truth_L_s",
  "truth_M_f", "truth_M_s", "K", "n_obs_min", "n_obs_max", "sigma_eps",
  "mean_amp", "paired_nonbasis_rng_required", "new_unseen_at_registration",
  "truth_available_to_fit", "truth_available_to_stopping",
  "truth_available_to_selection", "truth_unseal_authorized",
  "development_only", "formal_v0lv_result", "formal_paper_mc_result"
)]

fit_rows <- vector("list", nrow(data_manifest) * 12L)
cursor <- 0L
for (logical_index in 1:3) {
  paired_seeds <- 9092000L + (logical_index - 1L) * 12L + 1:12
  for (arm in c("legacy", "jaoua")) {
    data_row <- data_manifest[
      data_manifest$logical_data_index == logical_index &
        data_manifest$basis_arm == arm, , drop = FALSE]
    for (start_index in 1:12) {
      cursor <- cursor + 1L
      fit_rows[[cursor]] <- data.frame(
        experiment_id = experiment_id,
        fit_id = sprintf("%s__%s__g12_start_%02d",
                         data_row$logical_data_id, arm, start_index),
        data_id = data_row$data_id,
        logical_data_id = data_row$logical_data_id,
        logical_data_index = logical_index,
        basis_arm = arm,
        truth_basis_method = data_row$truth_basis_method,
        scenario_id = "baseline_strong", method_id = "G12",
        fit_config_id = "correct_L1_M2",
        selection_stratum_id = data_row$data_id,
        paired_start_id = sprintf("%s__start_%02d",
                                  data_row$logical_data_id, start_index),
        seed_index = start_index, fit_seed = paired_seeds[[start_index]],
        fit_L_f = 1L, fit_L_s_1 = 1L, fit_L_s_2 = 1L,
        fit_M_f = "2", fit_M_s_1 = "2", fit_M_s_2 = "2",
        route_id = "random_g12_prescore1_jaoua_fixed400",
        initialization = "random",
        function_initialization = "gram_unit_energy",
        pre_score_sweeps = 1L, anneal = "1;1.9;100",
        planned_annealing_sweeps = 99L,
        planned_ordinary_t1_sweeps = 400L, n_cpus = 1L,
        technical_gate = logical_index == 1L && start_index <= 4L,
        objective_eligibility_required = TRUE,
        truth_free_winner_candidate = TRUE,
        cross_arm_elbo_comparison_forbidden = TRUE,
        actual_racing = FALSE, continuation = FALSE, automatic_800 = FALSE,
        truth_available_to_fit = FALSE,
        truth_available_to_stopping = FALSE,
        truth_available_to_selection = FALSE,
        development_only = TRUE, formal_v0lv_result = FALSE,
        formal_paper_mc_result = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
}
fit_manifest <- do.call(rbind, fit_rows)

stopifnot(
  nrow(data_manifest) == 6L,
  nrow(fit_manifest) == 72L,
  !anyDuplicated(data_manifest$data_id),
  !anyDuplicated(fit_manifest$fit_id),
  all(table(data_manifest$logical_data_id) == 2L),
  all(table(fit_manifest$paired_start_id) == 2L),
  all(table(fit_manifest$selection_stratum_id) == 12L),
  sum(fit_manifest$technical_gate) == 8L,
  all(fit_manifest$n_cpus == 1L),
  !any(fit_manifest$truth_available_to_fit |
         fit_manifest$truth_available_to_stopping |
         fit_manifest$truth_available_to_selection)
)

write.csv(data_manifest, file.path(root, "DATA_MANIFEST.csv"),
          row.names = FALSE, quote = TRUE, na = "")
write.csv(fit_manifest, file.path(root, "FIT_MANIFEST.csv"),
          row.names = FALSE, quote = TRUE, na = "")
qc <- data.frame(
  check_id = c(
    "six_data_arms", "three_logical_pairs", "seventy_two_fits",
    "twelve_starts_per_arm", "thirty_six_paired_starts",
    "eight_gate_fits", "unique_fit_ids", "paired_fit_seeds",
    "single_core", "fixed_route", "truth_isolated", "no_continuation",
    "development_only", "not_formal"
  ),
  passed = c(
    nrow(data_manifest) == 6L,
    all(table(data_manifest$logical_data_id) == 2L),
    nrow(fit_manifest) == 72L,
    all(table(fit_manifest$selection_stratum_id) == 12L),
    length(unique(fit_manifest$paired_start_id)) == 36L &&
      all(table(fit_manifest$paired_start_id) == 2L),
    sum(fit_manifest$technical_gate) == 8L,
    !anyDuplicated(fit_manifest$fit_id),
    all(vapply(split(fit_manifest$fit_seed, fit_manifest$paired_start_id),
               function(value) length(unique(value)) == 1L, logical(1L))),
    all(fit_manifest$n_cpus == 1L),
    all(fit_manifest$route_id ==
          "random_g12_prescore1_jaoua_fixed400"),
    !any(fit_manifest$truth_available_to_fit |
           fit_manifest$truth_available_to_stopping |
           fit_manifest$truth_available_to_selection),
    !any(fit_manifest$continuation | fit_manifest$automatic_800 |
           fit_manifest$actual_racing),
    all(fit_manifest$development_only),
    !any(fit_manifest$formal_v0lv_result |
           fit_manifest$formal_paper_mc_result)
  ),
  detail = c("6", "3 x 2", "72", "12", "36 x 2", "8", "unique",
             "same within each arm pair", "1", "G12 fixed 400",
             "fit/stopping/selection FALSE", "all FALSE", "TRUE",
             "v0L-V and paper MC FALSE"),
  stringsAsFactors = FALSE
)
stopifnot(all(qc$passed))
write.csv(qc, file.path(root, "MANIFEST_QC.csv"),
          row.names = FALSE, quote = TRUE, na = "")
cat("MANIFEST_BUILD_PASS data=6 fits=72 gate=8 paired_starts=36\n")
