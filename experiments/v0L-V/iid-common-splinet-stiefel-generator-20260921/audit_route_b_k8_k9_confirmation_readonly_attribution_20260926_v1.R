#!/usr/bin/env Rscript

# Read-only post-unseal attribution for the Route-B K=8/9 independent confirmation.
# This script never fits, continues, or rewrites an existing endpoint.
options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x
abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) abort(message)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  assert(length(hit) == 1L, paste0("Duplicate --", name, "."))
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
input_root <- arg_value("input-root")
audit_root <- arg_value("audit-root")
assert(nzchar(input_root) && dir.exists(input_root),
       "A completed --input-root is required.")
assert(nzchar(audit_root) && !dir.exists(audit_root),
       "--audit-root must be a new directory.")
assert(requireNamespace("digest", quietly = TRUE), "digest is required.")

read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing CSV: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
}
write_lines <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  writeLines(x, path, useBytes = TRUE)
}
sha256 <- function(path) tolower(digest::digest(
  file = path, algo = "sha256", serialize = FALSE))

trap_weights <- function(grid) {
  grid <- as.numeric(grid)
  assert(length(grid) >= 2L && all(diff(grid) > 0), "Invalid grid.")
  delta <- diff(grid)
  c(delta[[1L]] / 2,
    (delta[-length(delta)] + delta[-1L]) / 2,
    delta[[length(delta)]] / 2)
}
weighted_norm <- function(x, w) sqrt(sum(w * as.numeric(x)^2))
normalize_function <- function(x, w) {
  value <- as.numeric(x)
  norm <- weighted_norm(value, w)
  assert(is.finite(norm) && norm > .Machine$double.eps,
         "Cannot normalize a zero function.")
  value / norm
}
weighted_ise <- function(x, y, w) sum(w * (as.numeric(x) - as.numeric(y))^2)
abs_cosine <- function(x, y) {
  nx <- sqrt(sum(as.numeric(x)^2)); ny <- sqrt(sum(as.numeric(y)^2))
  if (!is.finite(nx) || !is.finite(ny) || nx <= 0 || ny <= 0) return(0)
  abs(sum(as.numeric(x) * as.numeric(y)) / (nx * ny))
}
same_numeric <- function(x, y, tolerance = 0) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}
bind_rows <- function(values) {
  values <- Filter(Negate(is.null), values)
  if (!length(values)) return(data.frame())
  all_names <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(x) {
    missing <- setdiff(all_names, names(x))
    for (name in missing) x[[name]] <- NA
    x[all_names]
  })
  do.call(rbind, values)
}
finite_summary <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(c(median = NA, min = NA, max = NA))
  c(median = stats::median(x), min = min(x), max = max(x))
}

project_function <- function(truth_value, truth_grid, fit_grid, fit_design) {
  target_fit <- as.numeric(stats::approx(
    truth_grid, truth_value, xout = fit_grid, rule = 2, ties = mean)$y)
  fit_weights <- trap_weights(fit_grid)
  weighted_design <- sqrt(fit_weights) * fit_design
  coefficient <- qr.solve(weighted_design, sqrt(fit_weights) * target_fit)
  projected_fit <- as.numeric(fit_design %*% coefficient)
  list(
    fit_grid_value = projected_fit,
    truth_grid_value = as.numeric(stats::approx(
      fit_grid, projected_fit, xout = truth_grid,
      rule = 2, ties = mean)$y))
}

true_functions <- function(truth, role, study) {
  tp <- truth$data$true_params
  if (role == "shared") tp$phi_dense_list[[1L]] else
    tp$psi_dense_list[[as.integer(study)]][[1L]]
}
true_loading <- function(truth, role, study) {
  tp <- truth$data$true_params
  if (role == "shared") tp$a_true[, 1L] else
    tp$b_true[[as.integer(study)]][, 1L]
}
estimated_functions <- function(fit, role, study, factor) {
  if (role == "shared") fit$list_Phi_hat[[factor]] else
    fit$list_Phi_hat_spec[[as.integer(study)]][[factor]]
}
estimated_loading <- function(fit, role, study, factor) {
  if (role == "shared") fit$mu_q_a_hat[, factor] else
    fit$mu_q_b_specific_hat[[as.integer(study)]][, factor]
}

canonical_correlation <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  qx <- qr.Q(qr(x)); qy <- qr.Q(qr(y))
  max(svd(crossprod(qx, qy), nu = 0L, nv = 0L)$d)
}

continuous_orthonormal_coordinates <- function(g_cont, g_obs) {
  eig <- eigen((g_cont + t(g_cont)) / 2, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-10
  assert(any(keep), "Continuous fitted-basis Gram has no positive rank.")
  transform <- eig$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(eig$values[keep]), sum(keep))
  observed <- crossprod(transform, g_obs %*% transform)
  values <- eigen((observed + t(observed)) / 2,
                  symmetric = TRUE, only.values = TRUE)$values
  values <- pmax(values, .Machine$double.eps)
  list(
    rank = sum(keep), min_obs_eigen = min(values),
    max_obs_eigen = max(values),
    observed_condition = max(values) / min(values),
    max_continuous_to_observed_leverage = 1 / min(values))
}

manifest <- read_csv(file.path(input_root, "FIT_MANIFEST.csv"))
endpoints <- read_csv(file.path(
  input_root, "truth_free_selection", "ALL_TRUTH_FREE_ENDPOINTS.csv"))
winners <- read_csv(file.path(
  input_root, "truth_free_selection", "TRUTH_FREE_WINNERS.csv"))
science <- read_csv(file.path(
  input_root, "scientific_evaluation_20260926_v1",
  "ALL_ROUTE_B_SCIENTIFIC_RESULTS.csv"))
eval_status <- read_csv(file.path(
  input_root, "scientific_evaluation_20260926_v1", "EVALUATION_STATUS.csv"))
frozen_marker <- readLines(file.path(
  input_root, "truth_free_selection", "TRUTH_FREE_SELECTION_FROZEN.txt"),
  warn = FALSE)

assert(nrow(manifest) == 48L && nrow(endpoints) == 48L &&
         nrow(winners) == 6L && nrow(science) == 48L &&
         setequal(manifest$fit_id, science$fit_id),
       "Completed Route-B confirmation accounting changed.")
assert(all(as.logical(eval_status$ok)), "An existing evaluation is not OK.")

data_ids <- sort(unique(manifest$data_id))
truths <- setNames(lapply(data_ids, function(data_id) readRDS(file.path(
  input_root, "data", data_id, "sealed_truth", "truth_bundle.rds"))), data_ids)
observations <- setNames(lapply(data_ids, function(data_id) readRDS(file.path(
  input_root, "data", data_id, "observation_bundle.rds"))), data_ids)

binding_paths <- c(
  file.path(input_root, "FIT_MANIFEST.csv"),
  file.path(input_root, "truth_free_selection", "ALL_TRUTH_FREE_ENDPOINTS.csv"),
  file.path(input_root, "truth_free_selection", "TRUTH_FREE_WINNERS.csv"),
  unlist(lapply(data_ids, function(data_id) c(
    file.path(input_root, "data", data_id, "observation_bundle.rds"),
    file.path(input_root, "data", data_id, "sealed_truth", "truth_bundle.rds"))),
    use.names = FALSE),
  file.path(input_root, "fits", winners$fit_id, "fit.rds"))
binding_before <- data.frame(
  path = binding_paths, bytes = file.info(binding_paths)$size,
  sha256_before = vapply(binding_paths, sha256, character(1L)),
  stringsAsFactors = FALSE)

implementation_checks <- list()
add_check <- function(check_id, passed, detail) {
  implementation_checks[[length(implementation_checks) + 1L]] <<- data.frame(
    check_id = check_id, passed = isTRUE(passed), detail = as.character(detail),
    stringsAsFactors = FALSE)
}

add_check("selection_frozen_before_unseal",
          all(c("status=TRUTH_FREE_SELECTION_FROZEN",
                "truth_used_for_selection=FALSE",
                "truth_unseal_authorized=FALSE") %in% frozen_marker),
          "frozen marker retained and truth-free")
add_check("evaluation_accounting", all(as.logical(eval_status$ok)),
          paste0(sum(as.logical(eval_status$ok)), "/", nrow(eval_status)))

observability_rows <- list(); overlap_rows <- list(); obs_position <- 0L
for (data_id in data_ids) {
  truth <- truths[[data_id]]; observation <- observations[[data_id]]
  data <- truth$data; tp <- data$true_params
  parity <- identical(observation$Y, data$Y) &&
    identical(observation$Z, data$Z) &&
    identical(observation$time_obs, data$time_obs)
  add_check(paste0(data_id, "_sealed_observation_parity"), parity,
            "Y/Z/time_obs equal between observation-only and sealed truth")

  max_response_error <- 0
  for (s in seq_len(tp$S)) for (i in seq_len(tp$n_s[[s]])) {
    observed <- do.call(cbind, data$Y[[s]][[i]])
    max_response_error <- max(max_response_error, abs(
      observed - tp$signal_true_values[[s]][[i]] -
        tp$noise_true_values[[s]][[i]]))
  }
  add_check(paste0(data_id, "_response_decomposition"),
            max_response_error < 1e-10,
            format(max_response_error, scientific = TRUE))

  grid <- as.numeric(tp$t_grid_dense); w <- trap_weights(grid)
  blocks <- list(
    list(id = "A0_1", role = "shared", study = 0L,
         studies = seq_len(tp$S)),
    list(id = "A0_1", role = "shared", study = 1L, studies = 1L),
    list(id = "A0_1", role = "shared", study = 2L, studies = 2L),
    list(id = "B10_1", role = "specific", study = 1L, studies = 1L),
    list(id = "B20_1", role = "specific", study = 2L, studies = 2L))
  grid_gram_errors <- numeric(); loading_cosines <- numeric()
  basis_meta <- c(tp$phi_basis_meta,
                  unlist(tp$psi_basis_meta, recursive = FALSE))
  continuous_gram_errors <- vapply(basis_meta, function(meta) {
    max(abs(meta$gram - diag(nrow(meta$gram))))
  }, numeric(1L))
  for (block in blocks) {
    phi <- true_functions(truth, block$role,
                          if (block$study == 0L) 1L else block$study)
    grid_gram_errors <- c(grid_gram_errors,
      max(abs(crossprod(phi, w * phi) - diag(ncol(phi)))))
    times <- unlist(lapply(block$studies, function(s) {
      unlist(data$time_obs[[s]], use.names = FALSE)
    }), use.names = FALSE)
    observed_phi <- vapply(seq_len(ncol(phi)), function(component) {
      stats::approx(grid, phi[, component], xout = times,
                    rule = 2, ties = mean)$y
    }, numeric(length(times)))
    if (!is.matrix(observed_phi)) observed_phi <- matrix(observed_phi, ncol = 1L)
    obs_gram <- crossprod(observed_phi) / nrow(observed_phi)
    eig <- eigen((obs_gram + t(obs_gram)) / 2,
                 symmetric = TRUE, only.values = TRUE)$values
    edge <- grid <= 0.1 | grid >= 0.9
    edge_energy <- colSums(w[edge] * phi[edge, , drop = FALSE]^2)
    obs_position <- obs_position + 1L
    observability_rows[[obs_position]] <- data.frame(
      data_id = data_id, true_id = block$id, true_role = block$role,
      study_scope = block$study,
      observation_count = length(times),
      observed_boundary_fraction = mean(times <= 0.1 | times >= 0.9),
      observed_gram_min_eigen = min(eig),
      observed_gram_max_eigen = max(eig),
      observed_gram_condition = max(eig) / max(min(eig), .Machine$double.eps),
      max_continuous_to_observed_leverage =
        1 / max(min(eig), .Machine$double.eps),
      component_1_observed_energy = obs_gram[1L, 1L],
      component_2_observed_energy = obs_gram[2L, 2L],
      component_1_boundary_energy = edge_energy[[1L]],
      component_2_boundary_energy = edge_energy[[2L]],
      dense_201_trapezoid_gram_error = max(abs(
        crossprod(phi, w * phi) - diag(ncol(phi)))),
      stringsAsFactors = FALSE)
  }
  for (s in seq_len(tp$S)) {
    times <- unlist(data$time_obs[[s]], use.names = FALSE)
    shared <- vapply(seq_len(ncol(tp$phi_dense_list[[1L]])), function(m) {
      stats::approx(grid, tp$phi_dense_list[[1L]][, m], xout = times,
                    rule = 2, ties = mean)$y
    }, numeric(length(times)))
    specific_phi <- tp$psi_dense_list[[s]][[1L]]
    specific <- vapply(seq_len(ncol(specific_phi)), function(m) {
      stats::approx(grid, specific_phi[, m], xout = times,
                    rule = 2, ties = mean)$y
    }, numeric(length(times)))
    continuous_corr <- canonical_correlation(
      sqrt(w) * tp$phi_dense_list[[1L]], sqrt(w) * specific_phi)
    overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
      data_id = data_id, study = s,
      continuous_principal_correlation = continuous_corr,
      observed_principal_correlation = canonical_correlation(shared, specific),
      stringsAsFactors = FALSE)
    loading_cosines <- c(loading_cosines, abs_cosine(
      tp$a_true[, 1L], tp$b_true[[s]][, 1L]))
  }
  add_check(paste0(data_id, "_continuous_truth_gram"),
            max(continuous_gram_errors) < 1e-8,
            paste0("saved_continuous=",
                   format(max(continuous_gram_errors), scientific = TRUE),
                   ";dense_201_trapezoid=",
                   format(max(grid_gram_errors), scientific = TRUE)))
  add_check(paste0(data_id, "_strong_loading_identification"),
            max(loading_cosines) < 1e-10,
            format(max(loading_cosines), scientific = TRUE))
}
observability <- bind_rows(observability_rows)
overlap <- bind_rows(overlap_rows)

# K-specific fitted-space observability, computed in coordinates where the
# dense continuous Gram is the identity. This avoids dependence on the raw
# O'Sullivan coefficient scaling.
fit_space_rows <- list()
for (index in seq_len(nrow(winners))) {
  row <- winners[index, , drop = FALSE]
  fit <- readRDS(file.path(input_root, "fits", row$fit_id, "fit.rds"))
  w_fit <- trap_weights(fit$time_g)
  g_cont <- crossprod(fit$C_g, w_fit * fit$C_g)
  for (scope in 0:2) {
    if (scope == 0L) {
      g_obs <- Reduce(`+`, fit$sum_list_cp_C) /
        sum(vapply(fit$list_cp_C, length, integer(1L)))
      label <- "shared_pooled"
    } else {
      g_obs <- fit$sum_list_cp_C[[scope]] / length(fit$list_cp_C[[scope]])
      label <- paste0("study_", scope)
    }
    visibility <- continuous_orthonormal_coordinates(g_cont, g_obs)
    fit_space_rows[[length(fit_space_rows) + 1L]] <- data.frame(
      data_id = row$data_id, fit_K = row$fit_K, fit_id = row$fit_id,
      scope = label, coefficient_rank = visibility$rank,
      min_observed_eigen = visibility$min_obs_eigen,
      max_observed_eigen = visibility$max_obs_eigen,
      observed_condition = visibility$observed_condition,
      max_continuous_to_observed_leverage =
        visibility$max_continuous_to_observed_leverage,
      stringsAsFactors = FALSE)
  }
}
fit_space_visibility <- bind_rows(fit_space_rows)

# Independently reproduce the frozen winner component metrics, then partition
# both total and estimate-to-projection ISE into boundary and interior pieces.
component_rows <- list(); metric_differences <- numeric()
loading_differences <- numeric()
for (index in seq_len(nrow(winners))) {
  row <- winners[index, , drop = FALSE]
  data_id <- row$data_id[[1L]]; fit_id <- row$fit_id[[1L]]
  truth <- truths[[data_id]]; fit <- readRDS(file.path(
    input_root, "fits", fit_id, "fit.rds"))
  evaluation <- readRDS(file.path(
    input_root, "scientific_evaluation_20260926_v1", fit_id,
    "evaluation.rds"))
  feature <- evaluation$feature_component_errors
  feature <- feature[feature$component_matched, , drop = FALSE]
  grid <- truth$data$true_params$t_grid_dense; w <- trap_weights(grid)
  edge <- grid <= 0.1 | grid >= 0.9
  for (j in seq_len(nrow(feature))) {
    item <- feature[j, , drop = FALSE]
    role <- as.character(item$true_role[[1L]])
    study <- as.integer(item$study[[1L]])
    truth_matrix <- true_functions(truth, role, study)
    truth_value <- normalize_function(
      truth_matrix[, as.integer(item$true_component[[1L]])], w)
    projected <- project_function(truth_value, grid, fit$time_g, fit$C_g)
    projected_value <- normalize_function(projected$truth_grid_value, w)
    projection_sign <- sign(sum(w * truth_value * projected_value))
    if (!is.finite(projection_sign) || projection_sign == 0) projection_sign <- 1
    projected_value <- projection_sign * projected_value
    estimated_matrix <- estimated_functions(
      fit, as.character(item$estimated_block[[1L]]), study,
      as.integer(item$estimated_factor[[1L]]))
    estimate <- as.numeric(stats::approx(
      fit$time_g,
      estimated_matrix[, as.integer(item$estimated_component[[1L]])],
      xout = grid, rule = 2, ties = mean)$y)
    estimate <- normalize_function(estimate, w)
    estimate_sign <- sign(sum(w * truth_value * estimate))
    if (!is.finite(estimate_sign) || estimate_sign == 0) estimate_sign <- 1
    estimate <- estimate_sign * estimate
    total <- weighted_ise(estimate, truth_value, w)
    floor <- weighted_ise(projected_value, truth_value, w)
    excess <- weighted_ise(estimate, projected_value, w)
    metric_differences <- c(metric_differences,
      abs(total - item$total_sign_aligned_ise[[1L]]),
      abs(floor - item$projection_floor_ise[[1L]]),
      abs(excess - item$estimate_to_projection_ise[[1L]]))
    component_rows[[length(component_rows) + 1L]] <- data.frame(
      data_id = data_id, fit_K = row$fit_K[[1L]], fit_id = fit_id,
      truth_role = role, study = study, true_id = item$true_id[[1L]],
      true_component = item$true_component[[1L]],
      estimated_component = item$estimated_component[[1L]],
      total_ise = total, projection_floor_ise = floor,
      estimate_to_projection_ise = excess,
      total_ise_boundary = sum(w[edge] *
        (estimate[edge] - truth_value[edge])^2),
      total_ise_interior = sum(w[!edge] *
        (estimate[!edge] - truth_value[!edge])^2),
      estimate_to_projection_boundary = sum(w[edge] *
        (estimate[edge] - projected_value[edge])^2),
      estimate_to_projection_interior = sum(w[!edge] *
        (estimate[!edge] - projected_value[!edge])^2),
      stringsAsFactors = FALSE)
  }
  scale <- evaluation$loading_scale_by_block
  for (j in seq_len(nrow(scale))) {
    item <- scale[j, , drop = FALSE]
    if (!isTRUE(item$matched_by_loading[[1L]])) next
    role <- as.character(item$true_role[[1L]])
    study <- as.integer(item$study[[1L]])
    manual <- abs_cosine(
      estimated_loading(fit, as.character(item$estimated_block[[1L]]),
                        study, as.integer(item$estimated_factor[[1L]])),
      true_loading(truth, role, study))
    loading_differences <- c(loading_differences,
      abs(manual - item$canonical_loading_abs_cosine[[1L]]))
  }
}
component_decomposition <- bind_rows(component_rows)
add_check("independent_component_metric_reproduction",
          max(metric_differences) < 1e-8,
          format(max(metric_differences), scientific = TRUE))
add_check("independent_loading_match_reproduction",
          max(loading_differences) < 1e-10,
          format(max(loading_differences), scientific = TRUE))

# All correct endpoints, not merely winners: scale/path/contribution evidence.
correct <- science[as.logical(science$core_truth_recovered), , drop = FALSE]
scale_rows <- list()
for (index in seq_len(nrow(correct))) {
  row <- correct[index, , drop = FALSE]
  evaluation <- readRDS(file.path(
    input_root, "scientific_evaluation_20260926_v1", row$fit_id,
    "evaluation.rds"))
  scale <- evaluation$loading_scale_by_block
  process <- evaluation$process_contribution_by_role
  key <- c("study", "true_role", "true_factor", "estimated_block",
           "estimated_factor")
  joined <- merge(scale, process, by = key, all = FALSE, sort = FALSE)
  feature <- evaluation$feature_component_errors
  feature <- feature[feature$component_matched, , drop = FALSE]
  for (j in seq_len(nrow(joined))) {
    item <- joined[j, , drop = FALSE]
    f <- feature[feature$study == item$study &
                   feature$true_role == item$true_role &
                   feature$true_factor == item$true_factor, , drop = FALSE]
    scale_rows[[length(scale_rows) + 1L]] <- data.frame(
      data_id = row$data_id, fit_K = row$fit_K, fit_id = row$fit_id,
      start_index = row$start_index,
      study = item$study, true_role = item$true_role,
      canonical_loading_norm_ratio = item$canonical_loading_norm_ratio,
      factor_scale = item$factor_scale,
      true_path_rms = item$true_path_rms,
      estimated_path_rms = item$estimated_path_rms,
      path_rms_ratio = item$estimated_path_rms / item$true_path_rms,
      amplitude_compensation_product =
        item$canonical_loading_norm_ratio *
        item$estimated_path_rms / item$true_path_rms,
      factor_process_nrmse = item$factor_process_nrmse_matched,
      complete_contribution_nrmse = item$complete_contribution_nrmse_matched,
      feature_projection_floor_mean = mean(f$projection_floor_ise),
      feature_estimate_to_projection_mean =
        mean(f$estimate_to_projection_ise),
      feature_estimate_to_projection_max =
        max(f$estimate_to_projection_ise),
      stringsAsFactors = FALSE)
  }
}
scale_audit <- bind_rows(scale_rows)
scale_groups <- split(scale_audit, interaction(
  scale_audit$data_id, scale_audit$fit_K, scale_audit$study,
  scale_audit$true_role, drop = TRUE))
scale_summary <- bind_rows(lapply(scale_groups, function(x) {
  data.frame(
    data_id = x$data_id[[1L]], fit_K = x$fit_K[[1L]],
    study = x$study[[1L]], true_role = x$true_role[[1L]],
    endpoints = nrow(x),
    loading_norm_ratio_median = stats::median(x$canonical_loading_norm_ratio),
    loading_norm_ratio_range = paste(range(x$canonical_loading_norm_ratio),
                                     collapse = ";"),
    path_rms_ratio_median = stats::median(x$path_rms_ratio),
    path_rms_ratio_range = paste(range(x$path_rms_ratio), collapse = ";"),
    compensation_product_median =
      stats::median(x$amplitude_compensation_product),
    contribution_nrmse_median =
      stats::median(x$complete_contribution_nrmse),
    estimate_to_projection_mean_median =
      stats::median(x$feature_estimate_to_projection_mean),
    estimate_to_projection_max_max =
      max(x$feature_estimate_to_projection_max),
    stringsAsFactors = FALSE)
}))
scale_summary <- scale_summary[order(
  scale_summary$data_id, scale_summary$fit_K,
  scale_summary$study, scale_summary$true_role), , drop = FALSE]

dir.create(audit_root, recursive = TRUE, mode = "0700")
write_csv(observability, file.path(audit_root,
  "TRUTH_BLOCK_OBSERVABILITY.csv"))
write_csv(overlap, file.path(audit_root,
  "OBSERVED_AND_CONTINUOUS_BLOCK_OVERLAP.csv"))
write_csv(fit_space_visibility, file.path(audit_root,
  "FITTED_SPACE_OBSERVABILITY.csv"))
write_csv(component_decomposition, file.path(audit_root,
  "WINNER_COMPONENT_ERROR_DECOMPOSITION.csv"))
write_csv(scale_audit, file.path(audit_root,
  "CORRECT_ENDPOINT_SCALE_AUDIT.csv"))
write_csv(scale_summary, file.path(audit_root,
  "CORRECT_ENDPOINT_SCALE_SUMMARY.csv"))

checks <- bind_rows(implementation_checks)
write_csv(checks, file.path(audit_root, "IMPLEMENTATION_CHECKS.csv"))

binding_after <- binding_before
binding_after$sha256_after <- vapply(binding_paths, sha256, character(1L))
binding_after$unchanged <- binding_after$sha256_before ==
  binding_after$sha256_after
write_csv(binding_after, file.path(audit_root, "READONLY_INPUT_BINDING.csv"))
add_readonly <- data.frame(
  check_id = "input_hashes_unchanged",
  passed = all(binding_after$unchanged),
  detail = paste0(sum(binding_after$unchanged), "/", nrow(binding_after)),
  stringsAsFactors = FALSE)
checks <- rbind(checks, add_readonly)
# Rewrite is deliberately avoided: the final table is a distinct terminal copy.
write_csv(checks, file.path(audit_root, "IMPLEMENTATION_CHECKS_FINAL.csv"))

gate <- all(checks$passed)
write_csv(data.frame(
  gate_id = "NO_IMPLEMENTATION_ERROR_FOUND_IN_READONLY_AUDIT",
  passed = gate,
  existing_fits_read = 48L,
  new_fits_started = FALSE,
  continuations_started = FALSE,
  truth_used_post_selection_only = TRUE,
  independent_confirmation_may_proceed = gate,
  stringsAsFactors = FALSE),
  file.path(audit_root, "INDEPENDENT_CONFIRMATION_GATE.csv"))
write_lines(c(
  "status=ROUTE_B_K8_K9_CONFIRMATION_READONLY_ATTRIBUTION_COMPLETE",
  paste0("implementation_gate_passed=", gate),
  paste0("checks_passed=", sum(checks$passed), "/", nrow(checks)),
  "existing_fits_read=48",
  "new_fit_or_continuation_started=FALSE",
  "truth_used_post_selection_only=TRUE"),
  file.path(audit_root, "READONLY_ATTRIBUTION_COMPLETE.txt"))

cat("ROUTE_B_K8_K9_CONFIRMATION_READONLY_ATTRIBUTION_COMPLETE gate=", gate,
    " checks=", sum(checks$passed), "/", nrow(checks), "\n", sep = "")


