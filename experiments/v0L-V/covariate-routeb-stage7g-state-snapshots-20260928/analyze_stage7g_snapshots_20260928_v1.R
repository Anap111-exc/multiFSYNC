#!/usr/bin/env Rscript

# Read-only Stage 7G mechanism evaluation. Truth is opened only after the
# independent fit and endpoint-parity gate has passed.
options(warn = 1, stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analyze_stage7g_snapshots_20260928_v1.R STAGE7C_ROOT STAGE7G_ROOT")
}
stage7c_root <- normalizePath(args[[1L]], mustWork = TRUE)
stage7g_root <- normalizePath(args[[2L]], mustWork = TRUE)
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
sha256 <- function(path) {
  assert(file.exists(path), paste0("Missing bound file: ", path))
  line <- system2("sha256sum", shQuote(path), stdout = TRUE)
  assert(length(line) == 1L, paste0("sha256sum failed: ", path))
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
parity_root <- file.path(stage7g_root, "endpoint_parity_v1")
assert(file.exists(file.path(parity_root, "PARITY_PASS.txt")) &&
         !file.exists(file.path(parity_root, "PARITY_FAILED.txt")) &&
         file.exists(file.path(stage7g_root,
                               "FITS_AND_ENDPOINT_PARITY_COMPLETE.txt")),
       "Truth remains sealed until Stage 7G endpoint parity passes.")
parity_tables <- c("PARAMETER_PARITY.csv", "CONFIGURATION_PARITY.csv",
                   "ORIGINAL_TRACE_PARITY.csv",
                   "STAGE7F_DENSE_TRACE_PARITY.csv",
                   "SNAPSHOT_FINAL_STATE_PARITY.csv")
for (name in parity_tables) {
  table <- read.csv(file.path(parity_root, name))
  pass_column <- if ("pass" %in% names(table)) "pass" else
    "configuration_and_sweeps_equal"
  assert(nrow(table) > 0L && all(table[[pass_column]]),
         paste0("Parity table failed: ", name))
}
manifest <- read.csv(file.path(stage7g_root, "FIT_MANIFEST.csv"))
assert(nrow(manifest) == 3L &&
         identical(manifest$fit_id,
                   c("covrb7c_01__C1__start06",
                     "covrb7c_02__C1__start11",
                     "covrb7c_03__C1__start06")),
       "Three-fit manifest changed.")
truth_sha <- c(
  "2101ae33f95597b8756437e673a93e690a9d5e7872199307fe6735ed38b46934",
  "b6e843ee9d0fc6b346d618670ef16db513dd7eb10b9ef272bacf4d743c4aa86d",
  "1a0c472a6986768e571796cac4d28beb67f29a6ab20c9b05212dca77e60d88cf"
)
out <- file.path(stage7g_root, "analysis_v1")
assert(!file.exists(out) && !dir.exists(out),
       "Refusing to overwrite Stage 7G analysis.")

trap_weights <- function(grid) {
  delta <- diff(as.numeric(grid))
  assert(length(delta) >= 1L && all(is.finite(delta)) && all(delta > 0),
         "Invalid integration grid.")
  c(delta[[1L]] / 2,
    (head(delta, -1L) + tail(delta, -1L)) / 2,
    tail(delta, 1L) / 2)
}
abs_cosine <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  abs(sum(x * y) / denominator)
}
sum_trace_product <- function(x, y) sum(diag(as.matrix(x) %*% as.matrix(y)))
true_function_on_grid <- function(dense, dense_grid, fit_grid) {
  dense <- as.matrix(dense)
  values <- vapply(seq_len(ncol(dense)), function(m)
    approx(dense_grid, dense[, m], xout = fit_grid,
           rule = 2, ties = "ordered")$y,
    numeric(length(fit_grid)))
  matrix(values, nrow = length(fit_grid), ncol = ncol(dense))
}
time_score_l2_rms <- function(scores, functions, weights) {
  gram <- crossprod(functions, sweep(functions, 1L, weights, "*"))
  sqrt(mean(pmax(rowSums((scores %*% gram) * scores), 0)))
}
block_state <- function(snapshot, block, study) {
  if (block == "shared") {
    list(
      loading = snapshot$mu_q_a[, 1L],
      slab = snapshot$mu_q_normal_a[, 1L],
      loading_second = snapshot$term_a[, 1L],
      ppi = snapshot$mu_q_gamma_a[, 1L],
      coefficients = snapshot$mu_q_nu_phi[[1L]],
      coefficient_covariances = snapshot$Sigma_q_nu_phi[[1L]],
      scores = snapshot$mu_q_zeta[[study]][[1L]],
      score_covariances = snapshot$Sigma_q_zeta[[study]][[1L]]
    )
  } else {
    list(
      loading = snapshot$mu_q_b_specific[[study]][, 1L],
      slab = snapshot$mu_q_normal_b[[study]][, 1L],
      loading_second = snapshot$term_b_specific[[study]][, 1L],
      ppi = snapshot$mu_q_gamma_b[[study]][, 1L],
      coefficients = snapshot$mu_q_nu_psi[[study]][[1L]],
      coefficient_covariances = snapshot$Sigma_q_nu_psi[[study]][[1L]],
      scores = snapshot$mu_q_xi[[study]][[1L]],
      score_covariances = snapshot$Sigma_q_xi[[study]][[1L]]
    )
  }
}
truth_block <- function(truth, block, study, fit_grid) {
  dense <- if (block == "shared") truth$theta_dense_list[[1L]] else
    truth$kappa_dense_list[[study]][[1L]]
  list(
    loading = if (block == "shared") truth$a_true[, 1L] else
      truth$b_true[[study]][, 1L],
    scores = if (block == "shared") truth$eta_true[[study]][[1L]] else
      truth$chi_true[[study]][[1L]],
    functions = true_function_on_grid(dense, truth$t_grid_dense, fit_grid)
  )
}

state_rows <- list()
truth_rows <- list()
binding_rows <- list()
row_index <- 0L
for (fit_index in seq_len(nrow(manifest))) {
  item <- manifest[fit_index, ]
  fit_dir <- file.path(stage7g_root, "fits", item$fit_id)
  assert(file.exists(file.path(fit_dir, "FIT_COMPLETE.txt")) &&
           !file.exists(file.path(fit_dir, "FIT_ERROR.txt")),
         paste0("Incomplete fit: ", item$fit_id))
  fit_path <- file.path(fit_dir, "fit.rds")
  snapshots_path <- file.path(fit_dir, "RAW_STATE_SNAPSHOTS.rds")
  truth_path <- file.path(stage7c_root, "data", item$data_id,
                          "sealed_truth", "truth_bundle.rds")
  assert(identical(sha256(truth_path), truth_sha[[fit_index]]),
         paste0("Sealed truth hash changed: ", item$data_id))
  fit <- readRDS(fit_path)
  snapshots <- readRDS(snapshots_path)
  truth_bundle <- readRDS(truth_path)
  assert(identical(truth_bundle$bundle_class,
                   "covariate_routeb_pilot_sealed_truth") &&
           identical(truth_bundle$data_id, item$data_id),
         "Sealed truth identity failed.")
  truth <- truth_bundle$data$true_params
  assert(identical(vapply(snapshots, `[[`, integer(1L), "iteration"),
                   as.integer(0:499)) &&
           is.matrix(fit$C_g) &&
           identical(dim(fit$C_g), c(51L, 11L)) &&
           length(fit$time_g) == 51L &&
           length(fit$response_scale) == 500L &&
           all(fit$response_scale == 1),
         "Snapshot, basis or working-scale contract failed.")
  weights <- trap_weights(fit$time_g)
  basis_gram <- crossprod(fit$C_g,
                          sweep(fit$C_g, 1L, weights, "*"))
  target_loading <- list(shared = as.numeric(truth$a_true[, 1L]),
                         specific_1 = as.numeric(truth$b_true[[1L]][, 1L]),
                         specific_2 = as.numeric(truth$b_true[[2L]][, 1L]))
  target <- list()
  for (study in 1:2) for (block in c("shared", "specific")) {
    key <- paste(block, study, sep = "_")
    target[[key]] <- truth_block(truth, block, study, fit$time_g)
  }
  binding_rows[[fit_index]] <- data.frame(
    fit_id = item$fit_id, data_id = item$data_id,
    fit_sha256 = sha256(fit_path),
    snapshot_sha256 = sha256(snapshots_path),
    truth_sha256 = truth_sha[[fit_index]])
  for (snapshot in snapshots) {
    for (study in 1:2) for (block in c("shared", "specific")) {
      row_index <- row_index + 1L
      state <- block_state(snapshot, block, study)
      target_state <- target[[paste(block, study, sep = "_")]]
      loading <- as.numeric(state$loading)
      slab <- as.numeric(state$slab)
      loading_second <- as.numeric(state$loading_second)
      ppi <- as.numeric(state$ppi)
      coefficients <- as.matrix(state$coefficients)
      scores <- as.matrix(state$scores)
      function_mean <- fit$C_g %*% coefficients
      function_gram <- crossprod(
        function_mean, sweep(function_mean, 1L, weights, "*"))
      assert(length(loading) == 500L && length(ppi) == 500L &&
               length(loading_second) == 500L &&
               all(is.finite(loading)) && all(is.finite(ppi)) &&
               all(is.finite(loading_second)) &&
               nrow(scores) == 30L && ncol(scores) == 2L &&
               nrow(coefficients) == ncol(fit$C_g) &&
               ncol(coefficients) == 2L &&
               length(state$score_covariances) == 30L &&
               length(state$coefficient_covariances) == 2L,
             "Raw state shape or finite-value contract failed.")
      score_mean_sq <- rowSums(scores^2)
      score_cov_trace <- vapply(state$score_covariances,
                                function(x) sum(diag(as.matrix(x))),
                                numeric(1L))
      time_score_mean_sq <- rowSums((scores %*% function_gram) * scores)
      time_score_cov_sq <- vapply(state$score_covariances,
                                 sum_trace_product, numeric(1L),
                                 y = function_gram)
      function_mean_sq <- sum(sweep(function_mean^2, 1L, weights, "*"))
      function_cov_sq <- sum(vapply(state$coefficient_covariances,
                                   sum_trace_product, numeric(1L),
                                   y = basis_gram))
      loading_norm <- sqrt(sum(loading^2))
      loading_second_norm <- sqrt(sum(pmax(loading_second, 0)))
      score_mean_rms <- sqrt(mean(pmax(score_mean_sq, 0)))
      score_second_rms <- sqrt(mean(pmax(score_mean_sq +
                                           score_cov_trace, 0)))
      function_mean_l2 <- sqrt(max(function_mean_sq, 0))
      function_second_l2 <- sqrt(max(function_mean_sq +
                                       function_cov_sq, 0))
      time_score_mean_l2 <- sqrt(mean(pmax(time_score_mean_sq, 0)))
      time_score_second_l2 <- sqrt(mean(pmax(time_score_mean_sq +
                                               time_score_cov_sq, 0)))
      mean_signal_rms <- sqrt(mean(loading^2)) * time_score_mean_l2
      conditional_second_signal_rms <- sqrt(mean(pmax(loading_second, 0))) *
        time_score_second_l2
      phase <- if (snapshot$iteration == 0L) "initial" else
        if (isTRUE(snapshot$annealing)) "annealing" else "ordinary_T1"
      state_rows[[row_index]] <- data.frame(
        fit_id = item$fit_id, data_id = item$data_id,
        iteration = snapshot$iteration, temperature = snapshot$temperature,
        phase = phase, block = block, study = study,
        ppi_mean = mean(ppi), ppi_max = max(ppi),
        support_count_at_0_5 = sum(ppi >= 0.5),
        loading_mean_norm = loading_norm,
        loading_slab_norm = sqrt(sum(slab^2)),
        loading_second_norm = loading_second_norm,
        loading_mean_energy_fraction = if (loading_second_norm > 0)
          loading_norm^2 / loading_second_norm^2 else NA_real_,
        score_mean_rms = score_mean_rms,
        score_second_rms = score_second_rms,
        score_mean_energy_fraction = if (score_second_rms > 0)
          score_mean_rms^2 / score_second_rms^2 else NA_real_,
        function_mean_l2 = function_mean_l2,
        function_second_l2 = function_second_l2,
        time_score_mean_l2_rms = time_score_mean_l2,
        time_score_second_l2_rms = time_score_second_l2,
        mean_signal_l2_rms_over_variables = mean_signal_rms,
        conditional_second_signal_l2_rms_over_variables =
          conditional_second_signal_rms)

      true_loading <- as.numeric(target_state$loading)
      true_scores <- as.matrix(target_state$scores)
      true_functions <- as.matrix(target_state$functions)
      assert(length(true_loading) == length(loading) &&
               identical(dim(true_scores), dim(scores)) &&
               identical(dim(true_functions), dim(function_mean)),
             "Truth factor dimensions differ from fitted factor.")
      true_loading_norm <- sqrt(sum(true_loading^2))
      true_score_rms <- sqrt(mean(rowSums(true_scores^2)))
      true_function_l2 <- sqrt(sum(sweep(true_functions^2, 1L,
                                          weights, "*")))
      true_signal_rms <- sqrt(mean(true_loading^2)) *
        time_score_l2_rms(true_scores, true_functions, weights)
      true_support <- abs(true_loading) > 0
      selected_support <- ppi >= 0.5
      direction <- true_loading / true_loading_norm
      loading_variance <- pmax(loading_second - loading^2, 0)
      projection_fraction <- if (sum(loading_second) > 0) {
        (sum(loading * direction)^2 +
           sum(loading_variance * direction^2)) / sum(loading_second)
      } else NA_real_
      truth_rows[[row_index]] <- data.frame(
        fit_id = item$fit_id, data_id = item$data_id,
        iteration = snapshot$iteration, block = block, study = study,
        raw_loading_abs_cosine_to_target =
          abs_cosine(loading, true_loading),
        slab_loading_abs_cosine_to_target =
          abs_cosine(slab, true_loading),
        raw_loading_abs_cosine_to_true_shared =
          abs_cosine(loading, target_loading$shared),
        raw_loading_abs_cosine_to_true_specific_1 =
          abs_cosine(loading, target_loading$specific_1),
        raw_loading_abs_cosine_to_true_specific_2 =
          abs_cosine(loading, target_loading$specific_2),
        posterior_loading_direction_energy_fraction =
          projection_fraction,
        true_support_loading_second_energy_fraction =
          if (sum(loading_second) > 0)
            sum(loading_second[true_support]) / sum(loading_second) else
              NA_real_,
        true_support_count = sum(true_support),
        support_true_positive = sum(true_support & selected_support),
        support_false_positive = sum(!true_support & selected_support),
        loading_norm_ratio_to_truth = loading_norm / true_loading_norm,
        score_mean_rms_ratio_to_truth = score_mean_rms / true_score_rms,
        score_second_rms_ratio_to_truth = score_second_rms /
          true_score_rms,
        function_mean_l2_ratio_to_truth = function_mean_l2 /
          true_function_l2,
        mean_signal_ratio_to_truth = mean_signal_rms / true_signal_rms,
        true_signal_l2_rms_over_variables = true_signal_rms,
        loading_direction_interpretable_at_1pct_truth_norm =
          loading_norm >= 0.01 * true_loading_norm)
    }
  }
  rm(fit, snapshots, truth_bundle, truth)
  gc(FALSE)
}
state_table <- do.call(rbind, state_rows)
truth_table <- do.call(rbind, truth_rows)
binding <- do.call(rbind, binding_rows)
assert(nrow(state_table) == 6000L && nrow(truth_table) == 6000L &&
         all(is.finite(state_table$temperature)) &&
         all(is.finite(state_table$mean_signal_l2_rms_over_variables)),
       "Stage 7G analysis table contract failed.")
key_sweeps <- c(0:10, 15, 20, 25, 30, 35, 40, 45, 50, 55,
                60:75, 80, 90, 99, 100, 110, 150, 200, 300, 400, 499)
key <- merge(state_table[state_table$block == "specific" &
                           state_table$study == 2L &
                           state_table$iteration %in% key_sweeps, ],
             truth_table[truth_table$block == "specific" &
                           truth_table$study == 2L &
                           truth_table$iteration %in% key_sweeps, ],
             by = c("fit_id", "data_id", "iteration", "block", "study"),
             sort = FALSE)
key <- key[order(key$fit_id, key$iteration), ]
summary_rows <- lapply(manifest$fit_id, function(fit_id) {
  x <- key[key$fit_id == fit_id, ]
  early <- x[x$iteration <= 60L, ]
  interpretable <- early[early$loading_direction_interpretable_at_1pct_truth_norm, ]
  data.frame(
    fit_id = fit_id,
    early_max_raw_loading_abs_cosine_to_true_B2 =
      max(early$raw_loading_abs_cosine_to_target, na.rm = TRUE),
    early_max_interpretable_cosine_to_true_B2 =
      if (nrow(interpretable))
        max(interpretable$raw_loading_abs_cosine_to_target,
            na.rm = TRUE) else NA_real_,
    first_sweep_cosine_at_least_0_5 = {
      j <- which(early$raw_loading_abs_cosine_to_target >= 0.5 &
                   early$loading_direction_interpretable_at_1pct_truth_norm)
      if (length(j)) early$iteration[[j[[1L]]]] else NA_integer_
    },
    first_sweep_cosine_at_least_0_9 = {
      j <- which(early$raw_loading_abs_cosine_to_target >= 0.9 &
                   early$loading_direction_interpretable_at_1pct_truth_norm)
      if (length(j)) early$iteration[[j[[1L]]]] else NA_integer_
    },
    first_sweep_mean_signal_ratio_at_least_0_1 = {
      j <- which(early$mean_signal_ratio_to_truth >= 0.1)
      if (length(j)) early$iteration[[j[[1L]]]] else NA_integer_
    },
    endpoint_raw_loading_abs_cosine_to_true_B2 =
      x$raw_loading_abs_cosine_to_target[x$iteration == 499L],
    endpoint_mean_signal_ratio_to_truth =
      x$mean_signal_ratio_to_truth[x$iteration == 499L])
})
summary_table <- do.call(rbind, summary_rows)

dir.create(out)
write_new <- function(value, name) {
  path <- file.path(out, name)
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  write.csv(value, path, row.names = FALSE, na = "")
}
write_new(binding, "INPUT_BINDING.csv")
write_new(state_table, "STATE_METRICS_ALL_SWEEPS.csv")
write_new(truth_table, "TRUTH_ALIGNMENT_ALL_SWEEPS.csv")
write_new(key, "B2_KEY_SWEEPS.csv")
write_new(summary_table, "B2_DIRECTION_SUMMARY.csv")
writeLines(c("STAGE7G_ANALYSIS_COMPLETE", "endpoint_parity_pass_required=TRUE",
             "truth_used_only_after_parity_for_scientific_evaluation=TRUE",
             "truth_used_for_fit_or_selection=FALSE", "fits=3",
             "states_per_fit=500", "blocks_per_state=4",
             "signal_energy_is_posterior_mean_l2_over_subjects_variables_and_grid=TRUE",
             "conditional_second_signal_omits_time_function_covariance=TRUE"),
           file.path(out, "ANALYSIS_COMPLETE.txt"))
cat("STAGE7G_ANALYSIS_COMPLETE rows=", nrow(state_table), "\n", sep = "")
