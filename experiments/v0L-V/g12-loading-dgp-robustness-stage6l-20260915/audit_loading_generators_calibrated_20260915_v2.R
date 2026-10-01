#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
output_root <- get_arg("output-root")
if (is.null(output_root) || !nzchar(output_root)) {
  stop("Supply --output-root=<new directory>.", call. = FALSE)
}
if (dir.exists(output_root) || file.exists(output_root)) {
  stop("Refusing to overwrite output-root: ", output_root, call. = FALSE)
}
if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
  stop("Install the current multiFSYNC source before running the audit.",
       call. = FALSE)
}
if (!requireNamespace("Splinets", quietly = TRUE)) {
  stop("Splinets is required for the registered truth basis.", call. = FALSE)
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name) {
  utils::write.csv(x, file.path(output_root, name), row.names = FALSE,
                   na = "")
}
rms <- function(x) sqrt(mean(as.numeric(x)^2))
flatten <- function(x) unlist(x, recursive = TRUE, use.names = FALSE)
safe_jaccard <- function(x, y) {
  union_n <- sum(x | y)
  if (union_n == 0L) return(NA_real_)
  sum(x & y) / union_n
}
safe_cosine <- function(x, y) {
  denom <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  abs(sum(x * y) / denom)
}

scenario_table <- data.frame(
  scenario_id = c("controlled_identifiable_fixed90_calibrated",
                  "jaoua_beta10_blockwise_calibrated"),
  identified_loadings = c(TRUE, FALSE),
  sparsity_mode = c("fixed", "jaoua_beta_iid_normal"),
  stringsAsFactors = FALSE
)
seed_table <- data.frame(
  audit_index = seq_len(200L),
  data_seed = 96150000L + seq_len(200L),
  truth_basis_seed = 97150000L + seq_len(200L),
  truth_loading_seed = 98150000L + seq_len(200L),
  stringsAsFactors = FALSE
)
write_csv(seed_table, "GENERATION_AUDIT_SEEDS.csv")

base_config <- list(
  S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
  L_f = 1L, L_s = 1L, M_f = 2L, M_s = list(2L, 2L), K = 5L,
  n_obs = 7L, common_grid = FALSE, n_obs_range = c(6L, 9L),
  sigma_eps = 0.3, bool_sparse_loadings = TRUE, prop_sparse = 0.9,
  score_var_decay = TRUE, n_dense = 201L, mean_amp = 0.6,
  truth_basis_method = "jaoua_splinets", omega_beta = c(1, 10)
)

calibrate_constant_per_active <- function(data, anchor_active = 10L,
                                          anchor_norm = 3) {
  truth <- data$true_params
  truth$generator_omega_a_true <- truth$omega_a_true
  truth$generator_omega_b_true <- truth$omega_b_true
  calibrate_column <- function(value) {
    active_count <- sum(value != 0)
    current_norm <- sqrt(sum(value^2))
    if (active_count < 1L || !is.finite(current_norm) || current_norm <= 1e-14) {
      stop("Cannot calibrate an empty or non-finite loading column.")
    }
    target_norm <- anchor_norm * sqrt(active_count / anchor_active)
    target_norm * value / current_norm
  }
  truth$a_true[, 1L] <- calibrate_column(truth$a_true[, 1L])
  for (study in seq_len(truth$S)) {
    truth$b_true[[study]][, 1L] <-
      calibrate_column(truth$b_true[[study]][, 1L])
  }
  truth$gamma_a_true <- 1L * (truth$a_true != 0)
  truth$gamma_b_true <- lapply(truth$b_true, function(value) 1L * (value != 0))
  truth$omega_a_true <- colMeans(truth$gamma_a_true)
  truth$omega_b_true <- lapply(truth$gamma_b_true, colMeans)
  truth$loading_calibration <- "v0f_constant_per_active_loading"
  truth$v0f_anchor_active <- as.integer(anchor_active)
  truth$v0f_anchor_loading_norm <- anchor_norm

  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      signal <- truth$mu_true_values[[study]][[subject]] +
        truth$beta_true_values[[study]][[subject]] +
        truth$f_true_values[[study]][[subject]] %*% t(truth$a_true) +
        truth$g_true_values[[study]][[subject]] %*% t(truth$b_true[[study]])
      truth$signal_true_values[[study]][[subject]] <- signal
      noise <- truth$noise_true_values[[study]][[subject]]
      data$Y[[study]][[subject]] <- lapply(
        seq_len(truth$p),
        function(variable) as.numeric(signal[, variable] + noise[, variable])
      )
    }
  }
  data$true_params <- truth
  data
}

paired_top_fields <- c("time_obs", "C", "Z")
paired_truth_fields <- c(
  "phase_mu_mat", "phase_beta", "zeta_true", "xi_true", "eta_true",
  "chi_true", "zeta_fpca_true", "xi_fpca_true", "phi_dense_list",
  "psi_dense_list", "theta_dense_list", "kappa_dense_list",
  "lambda_phi_true", "lambda_psi_true", "mu_true_values",
  "f_true_values", "g_true_values", "beta_true_values",
  "noise_true_values", "n_obs_by_subject", "int_knots",
  "sigma2_eps_true")

block_rows <- list()
dataset_rows <- list()
pair_rows <- list()
row_block <- 0L
row_dataset <- 0L
row_pair <- 0L

summarize_arm <- function(dat, scenario_id, audit_index) {
  tp <- dat$true_params
  blocks <- list(
    shared_1 = list(value = tp$a_true[, 1L],
                    gamma = tp$gamma_a_true[, 1L],
                    omega = tp$generator_omega_a_true[[1L]]),
    specific_s1_1 = list(value = tp$b_true[[1L]][, 1L],
                         gamma = tp$gamma_b_true[[1L]][, 1L],
                         omega = tp$generator_omega_b_true[[1L]][[1L]]),
    specific_s2_1 = list(value = tp$b_true[[2L]][, 1L],
                         gamma = tp$gamma_b_true[[2L]][, 1L],
                         omega = tp$generator_omega_b_true[[2L]][[1L]])
  )
  block_table <- do.call(rbind, lapply(names(blocks), function(block_id) {
    block <- blocks[[block_id]]
    active <- block$gamma != 0
    data.frame(
      audit_index = audit_index, scenario_id = scenario_id,
      block_id = block_id, omega = block$omega,
      active_count = sum(active), active_fraction = mean(active),
      loading_l2 = sqrt(sum(block$value^2)),
      loading_l2_per_sqrt_active =
        sqrt(sum(block$value^2)) / sqrt(sum(active)),
      nonzero_abs_mean = mean(abs(block$value[active])),
      nonzero_mean = mean(block$value[active]),
      finite = all(is.finite(c(block$value, block$gamma, block$omega))),
      stringsAsFactors = FALSE
    )
  }))

  study_geometry <- lapply(1:2, function(s) {
    loading <- cbind(tp$a_true[, 1L], tp$b_true[[s]][, 1L])
    gram <- crossprod(loading)
    eigenvalues <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
    data.frame(
      study = s,
      shared_specific_support_jaccard = safe_jaccard(
        tp$gamma_a_true[, 1L] != 0, tp$gamma_b_true[[s]][, 1L] != 0),
      shared_specific_abs_cosine = safe_cosine(loading[, 1L], loading[, 2L]),
      loading_gram_min_eigenvalue = min(eigenvalues),
      loading_gram_max_eigenvalue = max(eigenvalues),
      loading_gram_condition = if (min(eigenvalues) > 0) {
        max(eigenvalues) / min(eigenvalues)
      } else Inf,
      stringsAsFactors = FALSE
    )
  })
  geometry <- do.call(rbind, study_geometry)
  signal <- flatten(tp$signal_true_values)
  mean_part <- flatten(tp$mu_true_values)
  beta_part <- flatten(tp$beta_true_values)
  noise <- flatten(tp$noise_true_values)
  factor_part <- signal - mean_part - beta_part
  response_error <- max(abs(
    flatten(dat$Y) - signal - noise
  ))
  dataset <- data.frame(
    audit_index = audit_index, scenario_id = scenario_id,
    shared_specific_jaccard_mean = mean(geometry$shared_specific_support_jaccard),
    shared_specific_abs_cosine_mean = mean(geometry$shared_specific_abs_cosine),
    loading_gram_min_eigenvalue = min(geometry$loading_gram_min_eigenvalue),
    loading_gram_condition_max = max(geometry$loading_gram_condition),
    signal_rms = rms(signal), factor_signal_rms = rms(factor_part),
    mean_rms = rms(mean_part), noise_rms = rms(noise),
    factor_to_noise_rms_ratio = rms(factor_part) / rms(noise),
    response_identity_max_abs_error = response_error,
    all_blocks_nonempty = all(block_table$active_count >= 1L),
    all_metrics_finite = all(is.finite(unlist(
      c(block_table[c("omega", "active_count", "active_fraction",
                      "loading_l2", "loading_l2_per_sqrt_active",
                      "nonzero_abs_mean", "nonzero_mean")],
        geometry[c("shared_specific_support_jaccard",
                   "shared_specific_abs_cosine",
                   "loading_gram_min_eigenvalue",
                   "loading_gram_max_eigenvalue",
                   "loading_gram_condition")])
    ))),
    stringsAsFactors = FALSE
  )
  list(block = block_table, dataset = dataset)
}

for (index in seed_table$audit_index) {
  seed_row <- seed_table[index, ]
  generated <- vector("list", 2L)
  names(generated) <- scenario_table$scenario_id
  for (arm in seq_len(nrow(scenario_table))) {
    scenario <- scenario_table[arm, ]
    config <- c(base_config, list(
      seed = seed_row$data_seed,
      truth_basis_seed = seed_row$truth_basis_seed,
      truth_loading_seed = seed_row$truth_loading_seed,
      identified_loadings = scenario$identified_loadings,
      sparsity_mode = scenario$sparsity_mode
    ))
    raw_generated <- do.call(
      multiFSYNC::simulate_multi_study_structured, config)
    generated[[scenario$scenario_id]] <-
      calibrate_constant_per_active(raw_generated)
    one <- summarize_arm(generated[[scenario$scenario_id]],
                         scenario$scenario_id, index)
    row_block <- row_block + 1L
    block_rows[[row_block]] <- one$block
    row_dataset <- row_dataset + 1L
    dataset_rows[[row_dataset]] <- one$dataset
  }

  controlled <- generated[["controlled_identifiable_fixed90_calibrated"]]
  jaoua <- generated[["jaoua_beta10_blockwise_calibrated"]]
  checks <- c(
    setNames(vapply(paired_top_fields, function(field) {
      identical(controlled[[field]], jaoua[[field]])
    }, logical(1L)), paste0("identical_", paired_top_fields)),
    setNames(vapply(paired_truth_fields, function(field) {
      identical(controlled$true_params[[field]],
                jaoua$true_params[[field]])
    }, logical(1L)), paste0("identical_", paired_truth_fields)),
    loading_differs = !identical(controlled$true_params$a_true,
                                 jaoua$true_params$a_true),
    response_differs = !identical(controlled$Y, jaoua$Y)
  )
  row_pair <- row_pair + 1L
  pair_rows[[row_pair]] <- data.frame(
    audit_index = index, check_id = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE)
  if (index %% 25L == 0L) message("Completed paired generator index ", index)
}

block_table <- do.call(rbind, block_rows)
dataset_table <- do.call(rbind, dataset_rows)
pair_table <- do.call(rbind, pair_rows)
write_csv(block_table, "BLOCK_LEVEL_AUDIT.csv")
write_csv(dataset_table, "DATASET_LEVEL_AUDIT.csv")
write_csv(pair_table, "PAIRING_CHECKS.csv")

quantile_value <- function(x, probability) {
  as.numeric(stats::quantile(x, probability, names = FALSE, na.rm = TRUE))
}
numeric_fields <- c(
  "shared_specific_jaccard_mean", "shared_specific_abs_cosine_mean",
  "loading_gram_min_eigenvalue", "loading_gram_condition_max",
  "signal_rms", "factor_signal_rms", "mean_rms", "noise_rms",
  "factor_to_noise_rms_ratio", "response_identity_max_abs_error")
summary_rows <- lapply(split(dataset_table, dataset_table$scenario_id),
                       function(x) {
  do.call(rbind, lapply(numeric_fields, function(field) {
    value <- x[[field]]
    data.frame(
      scenario_id = x$scenario_id[[1L]], block_id = "all_blocks",
      metric = field,
      n = sum(is.finite(value)), mean = mean(value), sd = stats::sd(value),
      q05 = quantile_value(value, 0.05), median = stats::median(value),
      q95 = quantile_value(value, 0.95), min = min(value), max = max(value),
      stringsAsFactors = FALSE)
  }))
})
summary_table <- do.call(rbind, summary_rows)
block_metrics <- c(
  "active_count", "omega", "loading_l2", "loading_l2_per_sqrt_active",
  "nonzero_abs_mean", "nonzero_mean")
block_summary <- do.call(rbind, lapply(
  split(block_table, interaction(block_table$scenario_id,
                                 block_table$block_id, drop = TRUE)),
  function(x) do.call(rbind, lapply(block_metrics, function(metric) {
    value <- x[[metric]]
    data.frame(
      scenario_id = x$scenario_id[[1L]], block_id = x$block_id[[1L]],
      metric = metric, n = nrow(x), mean = mean(value), sd = stats::sd(value),
      q05 = quantile_value(value, 0.05), median = stats::median(value),
      q95 = quantile_value(value, 0.95), min = min(value), max = max(value),
      stringsAsFactors = FALSE)
  })))
)
write_csv(rbind(summary_table, block_summary), "GENERATOR_AUDIT_SUMMARY.csv")

validity <- data.frame(
  check_id = c(
    "paired_dataset_count_200", "scenario_dataset_count_400",
    "block_row_count_1200", "all_loading_blocks_nonempty",
    "all_dataset_metrics_finite", "all_nonloading_pair_checks_pass",
    "response_identity_below_1e_12",
    "constant_per_active_calibration_below_1e_12"),
  passed = c(
    length(unique(dataset_table$audit_index)) == 200L,
    nrow(dataset_table) == 400L,
    nrow(block_table) == 1200L,
    all(dataset_table$all_blocks_nonempty),
    all(dataset_table$all_metrics_finite),
    all(pair_table$passed),
    max(dataset_table$response_identity_max_abs_error) < 1e-12,
    max(abs(block_table$loading_l2_per_sqrt_active - 3 / sqrt(10))) < 1e-12),
  stringsAsFactors = FALSE)
write_csv(validity, "GENERATOR_AUDIT_VALIDITY.csv")
writeLines(c(
  paste0("status=", if (all(validity$passed)) "COMPLETE" else "FAILED"),
  "fit_count=0", "truth_used_for_generator_audit=TRUE",
  "formal_v0lv_result=FALSE", "formal_paper_mc_result=FALSE",
  paste0("completed_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(output_root, "GENERATION_AUDIT_COMPLETE.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_root, "SESSION_INFO.txt"))
if (!all(validity$passed)) stop("Generator audit validity failed.", call. = FALSE)
