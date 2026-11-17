# Isolated evaluation adapters. Frozen legacy sources remain byte-identical.
nc_one_to_one_match_pruned <- function(similarity, threshold = .8) {
  # Same exhaustive traversal/objective/tolerance, with mathematically safe upper-bound pruning.
  # This avoids evaluating unmatched branches after a maximal-cardinality match exists.
  similarity <- as.matrix(similarity); nt <- nrow(similarity); nc <- ncol(similarity)
  if (!nt) return(list(assignment = integer(), match_count = 0L, total_similarity = 0))
  best <- rep.int(nc + 1L, nt); best_count <- -1L; best_score <- -Inf
  current <- rep.int(nc + 1L, nt); used <- rep(FALSE, nc)
  eligible_matrix <- is.finite(similarity) & similarity >= threshold
  lex_less <- function(x, y) {k <- which(x != y); length(k) > 0L && x[k[1L]] < y[k[1L]]}
  visit <- function(row, count, score) {
    if (row > nt) {
      better <- count > best_count || (count == best_count && score > best_score + 1e-12) ||
        (count == best_count && abs(score - best_score) <= 1e-12 && lex_less(current, best))
      if (better) {best <<- current; best_count <<- count; best_score <<- score}
      return(invisible(NULL))
    }
    left <- row:nt; available <- which(!used)
    valid <- eligible_matrix[left, available, drop = FALSE]
    upper_count <- count + min(length(left), length(available), sum(rowSums(valid) > 0))
    if (upper_count < best_count) return(invisible(NULL))
    if (upper_count == best_count && is.finite(best_score)) {
      # Ignoring competition yields an upper bound; include a generous roundoff margin.
      gains <- if (length(available)) vapply(left, function(k) {
        okay <- available[eligible_matrix[k, available]]
        if (length(okay)) max(similarity[k, okay]) else 0
      }, numeric(1)) else numeric(length(left))
      needed <- best_count - count
      upper_score <- score + sum(head(sort(gains, decreasing = TRUE), needed))
      margin <- 1e-12 + 100 * .Machine$double.eps * (1 + abs(upper_score) + abs(best_score))
      if (upper_score < best_score - margin) return(invisible(NULL))
    }
    eligible <- which(!used & eligible_matrix[row, ])
    for (candidate in eligible) {
      used[candidate] <<- TRUE; current[row] <<- candidate
      visit(row + 1L, count + 1L, score + similarity[row, candidate])
      used[candidate] <<- FALSE
    }
    current[row] <<- nc + 1L; visit(row + 1L, count, score)
    invisible(NULL)
  }
  visit(1L, 0L, 0); best[best > nc] <- 0L
  list(assignment = as.integer(best), match_count = as.integer(best_count), total_similarity = best_score)
}

nc_v0g_shared_candidate_metrics <- function(fit, paths, cross_study_threshold = .10) {
  rows <- lapply(seq_len(fit$L_f), function(factor) {
    loading <- fit$mu_q_a_hat[, factor]
    contribution <- vapply(seq_len(fit$S), function(study)
      sqrt(mean(paths$shared[[study]][, factor]^2)) * sqrt(mean(loading^2)), numeric(1))
    pooled <- sqrt(mean(unlist(lapply(paths$shared, function(x) x[, factor]))^2)) * sqrt(mean(loading^2))
    ans <- data.frame(candidate_factor = factor, factor_ppi = as.numeric(fit$factor_ppi_shared[factor]),
      retained = fit$factor_ppi_shared[factor] >= .5, loading_norm = sqrt(sum(loading^2)),
      contribution_rms_pooled = pooled, study_min_max_ratio = if (max(contribution) > 0) min(contribution)/max(contribution) else 0)
    for (s in seq_len(fit$S)) ans[[paste0('contribution_rms_study_', s)]] <- contribution[s]
    ans$nc_min_contribution_rms_all_studies <- min(contribution)
    ans
  })
  result <- do.call(rbind, rows)
  retained_pool_max <- max(result$contribution_rms_pooled[result$retained], 0)
  result$standardized_min_contribution <- if (retained_pool_max > 0)
    result$nc_min_contribution_rms_all_studies / retained_pool_max else 0
  result$cross_study_supported_0_1 <- result$retained &
    result$standardized_min_contribution >= cross_study_threshold
  result
}

nc_v0i_complete_operator_error <- function(fit, data) {
  tr <- data$true_params; w <- stage7m_weights(fit$time_g)
  interpolate <- function(x) vapply(seq_len(ncol(x)), function(j)
    approx(tr$t_grid_dense, x[, j], xout = fit$time_g, rule = 2)$y, numeric(length(fit$time_g)))
  truth_shared <- lapply(tr$theta_dense_list, interpolate)
  truth_specific <- lapply(tr$kappa_dense_list, function(x) lapply(x, interpolate))
  estimated_shared <- lapply(seq_along(fit$list_Phi_hat), function(l)
    fit$list_Phi_hat[[l]] %*% diag(sqrt(fit$list_eigenvalues[[l]]), nrow = length(fit$list_eigenvalues[[l]])))
  estimated_specific <- lapply(seq_len(fit$S), function(s)
    lapply(seq_along(fit$list_Phi_hat_spec[[s]]), function(l)
      fit$list_Phi_hat_spec[[s]][[l]] %*%
        diag(sqrt(fit$list_eigenvalues_spec[[s]][[l]]), nrow = length(fit$list_eigenvalues_spec[[s]][[l]]))))
  atom_inner <- function(a, f, b, g) sum(a * b)^2 * sum(crossprod(f, w * g)^2)
  all_inner <- function(load1, fun1, load2, fun2) {
    sum(vapply(seq_along(load1), function(i) sum(vapply(seq_along(load2), function(j)
      atom_inner(load1[[i]], fun1[[i]], load2[[j]], fun2[[j]]), numeric(1))), numeric(1)))
  }
  do.call(rbind, lapply(seq_len(fit$S), function(s) {
    tl <- c(lapply(seq_len(tr$L_f), function(l) tr$a_true[, l]),
      lapply(seq_len(ncol(tr$b_true[[s]])), function(l) tr$b_true[[s]][, l]))
    tf <- c(truth_shared, truth_specific[[s]])
    el <- c(lapply(seq_len(fit$L_f), function(l) fit$mu_q_a_hat[, l]),
      lapply(seq_len(ncol(fit$mu_q_b_specific_hat[[s]])), function(l) fit$mu_q_b_specific_hat[[s]][, l]))
    ef <- c(estimated_shared, estimated_specific[[s]])
    tn <- all_inner(tl, tf, tl, tf); en <- all_inner(el, ef, el, ef)
    cross <- all_inner(tl, tf, el, ef); sqerr <- max(0, tn + en - 2 * cross)
    data.frame(study = s, truth_operator_hs_norm = sqrt(tn), estimated_operator_hs_norm = sqrt(en),
      covariance_operator_relative_error = sqrt(sqerr / max(tn, .Machine$double.eps)),
      nc_operator_estimand = 'finite-rank training FPCA summary plug-in operator; not new-subject posterior covariance')
  }))
}

nc_evaluation_setup <- function(exec_dir) {
  dep <- file.path(exec_dir, 'dependencies'); legacy <- file.path(dep, 'legacy')
  expected <- c(
    v0g_evaluation.R = '6b765191657f88db46a027a95acd3917c602b35cb6cee8947d917652f9c6b271',
    v0h_evaluation.R = '45639c2d9f14e4df0187d519a62837ef3b3be090ed315123cd97c725f8576dd4',
    v0i_truth_diagnostic.R = '127615551ef1e4dbe14b8c486e3e3ad27ba45eea77f764c7ed7c2fba8318b194',
    v0j_evaluation.R = '3c91af1e07b0ba1aac4539a2c1e28d8abcfe2f5ee239a2a2daed7fcb2714a833',
    candidate2_evaluation.R = '5ef50d51572e946e158a6b67ee28601c4b8d8847f9830fa11d881e682df24d1e',
    stage6a_common.R = 'bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6',
    stage6h_helpers.R = '5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221',
    stage7c_evaluator.R = '1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8')
  for (nm in names(expected)) check_hash(file.path(legacy, nm), expected[[nm]])
  extract_functions(file.path(legacy, 'stage6a_common.R'),
    c('pp_abort', 'pp_assert', 'pp_nrmse', 'pp_vector_l2', 'pp_vector_rms', 'pp_abs_cosine', 'pp_dimension_activity_metrics'))
  source(file.path(legacy, 'candidate2_evaluation.R'), local = .GlobalEnv)
  source(file.path(legacy, 'stage6h_helpers.R'), local = .GlobalEnv)
  assign('PP_SNAPSHOT', legacy, .GlobalEnv)
  env <- new.env(parent = .GlobalEnv)
  for (nm in c('v0g_evaluation.R', 'v0h_evaluation.R', 'v0i_truth_diagnostic.R', 'v0j_evaluation.R'))
    sys.source(file.path(legacy, nm), envir = env, keep.source = TRUE)
  evaluation_env <- env
  env$exists <- function(x, where = -1, envir = NULL, mode = 'any', inherits = TRUE)
    base::exists(x, envir = evaluation_env, mode = mode, inherits = inherits)
  # Carry forward the already verified independent-mean truth evaluator.
  extract_functions(file.path(legacy, 'stage7c_evaluator.R'),
    c('trapezoid_weights', 'interp', 'mean_design', 'true_mean_dense',
      'estimated_mean_dense', 'install_dense_residual_override'))
  install_dense_residual_override(env)
  # Generalize finite-rank atom summation and min-over-study loops only.
  f <- nc_v0g_shared_candidate_metrics; environment(f) <- env
  env$v0g_shared_candidate_metrics <- f
  f <- nc_v0i_complete_operator_error; environment(f) <- env
  env$v0i_complete_operator_error <- f
  f <- nc_one_to_one_match_pruned; environment(f) <- env
  env$v0h_one_to_one_match <- f
  attr(env, 'nc_source_hashes') <- expected
  attr(env, 'nc_override_contract') <- data.frame(
    binding = c('v0j_dense_signal_metrics', 'v0g_shared_candidate_metrics', 'v0i_complete_operator_error', 'v0h_one_to_one_match'),
    reason = c('Reuse frozen Stage7C independent Splinets mean assembly',
      'Min contribution over S studies; same .10 threshold and RMS formula',
      'Same Hilbert-Schmidt atom inner product summed over all L_f and L_s columns',
      'Same exhaustive traversal/count-score-lex objective and 1e-12 tie tolerance; only upper-bound pruning'),
    legacy_source_unchanged = TRUE, model_updates_changed = FALSE)
  env
}

nc_general_structure <- function(fit, truth, historical, env) {
  tr <- truth$data$true_params; z <- historical$summary
  ls <- as.integer(tr$L_s_by_study)
  selected <- c(sum(fit$factor_ppi_shared >= .5),
    vapply(fit$factor_ppi_specific, function(x) sum(x >= .5), integer(1)))
  required <- c(tr$L_f, ls)
  core <- z$global_correct == z$truth_global_direction_count && z$global_misplaced == 0L &&
    z$global_missing == 0L && z$study_correct == z$truth_study_role_count &&
    z$study_misplaced == 0L && z$study_missing == 0L
  exact_counts <- identical(as.integer(selected), as.integer(required))
  exact <- core && exact_counts && z$global_extra == 0L && z$global_duplicate == 0L &&
    z$study_candidate_extra == 0L && z$study_candidate_duplicate == 0L
  components <- historical$feature_component_errors
  assert(nrow(components) == 2L * (tr$S * tr$L_f + sum(ls)), 'Legacy study-role component count omitted roles')
  full_components <- all(components$component_matched)
  block <- historical$primary$block_subspace
  current <- block[block$scope == 'ppi_selected', , drop = FALSE]
  correct_blocks <- current[current$metric_class == 'correct_block', , drop = FALSE]
  leak <- current[current$metric_class == 'cross_block_leakage', , drop = FALSE]
  assert(nrow(correct_blocks) == fit$S + 1L && nrow(leak) == 2L * fit$S,
         'R/P/L tables omitted a study')
  shared <- components[components$true_role == 'shared', , drop = FALSE]
  shared_consistent <- TRUE
  for (l in seq_len(tr$L_f)) for (m in 1:2) {
    rows <- shared[shared$true_factor == l & shared$true_component == m, , drop = FALSE]
    matched <- rows$component_matched
    if (any(matched)) shared_consistent <- shared_consistent &&
      all(matched) && all(rows$estimated_block == 'shared') &&
      length(unique(rows$estimated_factor)) == 1L && length(unique(rows$estimated_component)) == 1L
  }
  unique_keys <- paste(components$true_role, ifelse(components$true_role == 'shared', 0L, components$study),
                       components$true_factor, components$true_component, sep = ':')
  unique_coverage <- vapply(split(seq_len(nrow(components)), unique_keys),
    function(index) all(components$component_matched[index]), logical(1))
  summary <- data.frame(nc_truth_L_f = tr$L_f, nc_truth_S = tr$S,
    nc_truth_L_s_by_study = paste(ls, collapse = ';'), nc_selected_counts = paste(selected, collapse = ';'),
    nc_ppi_factor_counts_exact = exact_counts, nc_core_truth_recovered = core,
    nc_exact_thresholded_factor_structure = exact,
    nc_full_true_feature_component_coverage = full_components,
    nc_full_structural_recovery = exact && full_components,
    nc_study_role_features_total = nrow(components),
    nc_global_unique_features_total = length(unique_coverage),
    nc_global_unique_features_matched = sum(unique_coverage),
    nc_global_unique_feature_coverage = all(unique_coverage),
    nc_shared_function_component_map_consistent_all_studies = shared_consistent,
    nc_min_correct_block_R_all_studies = min(correct_blocks$subspace_recall),
    nc_min_correct_block_P_all_studies = min(correct_blocks$subspace_purity),
    nc_max_cross_block_L_all_studies = max(leak$subspace_recall),
    nc_legacy_B1_B2_summary_scope = if (fit$S == 2L) 'all studies' else 'first two studies only; use nc all-study summaries',
    nc_beta_metrics_applicability = 'NA: generated d=0 and fitted d=0')
  by_block <- do.call(rbind, lapply(0:tr$S, function(s) {
    role <- if (s == 0L) 'shared' else 'specific'
    table <- historical$primary$study_role_matches
    study_rows <- if (s == 0L) table[table$scope == 'ppi_selected' & table$true_block == 'shared', , drop = FALSE] else
      table[table$scope == 'ppi_selected' & table$unit_study == s & table$true_block == 'specific', , drop = FALSE]
    global <- historical$primary$global_direction_matches
    rows <- if (s == 0L) global[global$scope == 'ppi_selected' & global$true_block == 'shared', , drop = FALSE] else study_rows
    data.frame(role = role, study = s, expected_factors = if (s == 0L) tr$L_f else ls[s],
      selected_factors = selected[s + 1L], classification_unit = if (s == 0L) 'global_direction' else 'study_role',
      study_role_rows = nrow(study_rows), study_role_correct = sum(study_rows$status == 'correct'),
      study_role_missing = sum(study_rows$status == 'missing'), study_role_misplaced = sum(study_rows$status == 'misplaced'),
      correct = sum(rows$status == 'correct'), missing = sum(rows$status == 'missing'),
      misplaced = sum(rows$status == 'misplaced'))
  }))
  list(summary = summary, by_block = by_block, block_subspace = block,
    feature_component_coverage = data.frame(global_unique_key = names(unique_coverage), matched_all_required_studies = unique_coverage))
}
