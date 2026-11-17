v0h_default_thresholds <- function() {
  list(
    factor_ppi = 0.5,
    loading_abs_cosine = 0.8,
    duplicate_abs_cosine = 0.8,
    trajectory_abs_cor = 0.8,
    svd_relative_tolerance = 1e-8
  )
}

v0h_safe_cosine <- function(x, y, absolute = TRUE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) return(0)
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) return(0)
  value <- sum(x * y) / denominator
  value <- max(-1, min(1, value))
  if (absolute) abs(value) else value
}

v0h_safe_cor <- function(x, y, absolute = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L ||
      !all(is.finite(x)) || !all(is.finite(y)) ||
      stats::sd(x) <= 0 || stats::sd(y) <= 0) return(0)
  value <- suppressWarnings(stats::cor(x, y))
  if (!is.finite(value)) return(0)
  if (absolute) abs(value) else value
}

v0h_nrmse <- function(estimate, truth) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
  if (length(estimate) != length(truth) || !length(truth) ||
      !all(is.finite(estimate)) || !all(is.finite(truth))) return(NA_real_)
  denominator <- sqrt(mean(truth^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sqrt(mean((estimate - truth)^2)) / denominator
}

v0h_loading_sign <- function(estimate, truth) {
  value <- sum(as.numeric(estimate) * as.numeric(truth))
  if (!is.finite(value) || value == 0) 1 else sign(value)
}

v0h_empty_matrix <- function(nrow) {
  matrix(numeric(), nrow = as.integer(nrow), ncol = 0L)
}

v0h_subset_columns <- function(value, keep) {
  value <- as.matrix(value)
  keep <- as.integer(keep)
  if (!length(keep)) return(v0h_empty_matrix(nrow(value)))
  value[, keep, drop = FALSE]
}

v0h_trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) == 1L) return(1)
  if (length(grid) < 1L || any(!is.finite(grid)) ||
      any(diff(grid) <= 0)) stop("grid must be finite and strictly increasing")
  delta <- diff(grid)
  c(delta[[1L]] / 2,
    (delta[-1L] + delta[-length(delta)]) / 2,
    delta[[length(delta)]] / 2)
}

v0h_orthonormal_basis <- function(value, weights = NULL,
                                  relative_tolerance = 1e-8) {
  value <- as.matrix(value)
  if (!nrow(value)) stop("value must have at least one row")
  if (!ncol(value)) {
    return(list(
      basis = v0h_empty_matrix(nrow(value)), rank = 0L,
      singular_values = numeric()
    ))
  }
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(value) || any(weights < 0) ||
        any(!is.finite(weights))) stop("invalid weights")
    value <- sqrt(weights) * value
  }
  decomposition <- svd(value, nu = min(dim(value)), nv = 0L)
  maximum <- max(decomposition$d, 0)
  rank <- if (maximum <= 0) 0L else {
    sum(decomposition$d > maximum * relative_tolerance)
  }
  basis <- if (rank) {
    decomposition$u[, seq_len(rank), drop = FALSE]
  } else {
    v0h_empty_matrix(nrow(value))
  }
  list(
    basis = basis,
    rank = as.integer(rank),
    singular_values = decomposition$d
  )
}

v0h_subspace_metrics <- function(truth, estimate, weights = NULL,
                                 relative_tolerance = 1e-8) {
  truth <- as.matrix(truth)
  estimate <- as.matrix(estimate)
  if (nrow(truth) != nrow(estimate)) {
    stop("truth and estimate must have the same number of rows")
  }
  truth_working <- truth
  estimate_working <- estimate
  if (!is.null(weights)) {
    truth_working <- sqrt(as.numeric(weights)) * truth_working
    estimate_working <- sqrt(as.numeric(weights)) * estimate_working
  }
  truth_basis <- v0h_orthonormal_basis(
    truth, weights, relative_tolerance
  )
  estimate_basis <- v0h_orthonormal_basis(
    estimate, weights, relative_tolerance
  )
  cross <- crossprod(truth_basis$basis, estimate_basis$basis)
  overlap <- if (length(cross)) sum(cross^2) else 0
  recall <- if (truth_basis$rank) overlap / truth_basis$rank else 0
  purity <- if (estimate_basis$rank) overlap / estimate_basis$rank else 0
  projection <- if (ncol(truth_working)) {
    vapply(seq_len(ncol(truth_working)), function(column) {
      vector <- truth_working[, column]
      denominator <- sum(vector^2)
      if (!is.finite(denominator) || denominator <= 0 ||
          !estimate_basis$rank) return(0)
      sum(crossprod(estimate_basis$basis, vector)^2) / denominator
    }, numeric(1))
  } else numeric()
  energy_denominator <- sum(truth_working^2)
  energy_coverage <- if (energy_denominator > 0 && estimate_basis$rank) {
    sum(crossprod(estimate_basis$basis, truth_working)^2) /
      energy_denominator
  } else 0
  principal_cosine_sq <- if (length(cross)) {
    pmin(1, svd(cross, nu = 0L, nv = 0L)$d^2)
  } else numeric()
  list(
    truth_rank = truth_basis$rank,
    estimate_rank = estimate_basis$rank,
    overlap_trace = overlap,
    subspace_recall = recall,
    subspace_purity = purity,
    truth_column_coverage_mean = if (length(projection)) mean(projection) else 0,
    truth_column_coverage_min = if (length(projection)) min(projection) else 0,
    truth_column_coverage_max = if (length(projection)) max(projection) else 0,
    truth_energy_trace_coverage = energy_coverage,
    principal_cosine_sq_mean = if (length(principal_cosine_sq)) {
      mean(principal_cosine_sq)
    } else 0,
    principal_cosine_sq_min = if (length(principal_cosine_sq)) {
      min(principal_cosine_sq)
    } else 0
  )
}

v0h_lexicographically_less <- function(x, y) {
  difference <- which(x != y)
  if (!length(difference)) return(FALSE)
  x[[difference[[1L]]]] < y[[difference[[1L]]]]
}

v0h_one_to_one_match <- function(similarity, threshold = 0.8) {
  similarity <- as.matrix(similarity)
  n_truth <- nrow(similarity)
  n_candidate <- ncol(similarity)
  if (!n_truth) {
    return(list(
      assignment = integer(), match_count = 0L, total_similarity = 0
    ))
  }
  best_assignment <- rep.int(n_candidate + 1L, n_truth)
  best_count <- -1L
  best_score <- -Inf
  current <- rep.int(n_candidate + 1L, n_truth)
  used <- rep(FALSE, n_candidate)
  visit <- function(row, count, score) {
    if (row > n_truth) {
      better <- count > best_count ||
        (count == best_count && score > best_score + 1e-12) ||
        (count == best_count && abs(score - best_score) <= 1e-12 &&
           v0h_lexicographically_less(current, best_assignment))
      if (better) {
        best_assignment <<- current
        best_count <<- count
        best_score <<- score
      }
      return(invisible(NULL))
    }
    eligible <- which(!used & is.finite(similarity[row, ]) &
                        similarity[row, ] >= threshold)
    if (length(eligible)) {
      for (candidate in eligible) {
        used[[candidate]] <<- TRUE
        current[[row]] <<- candidate
        visit(
          row + 1L, count + 1L,
          score + similarity[row, candidate]
        )
        used[[candidate]] <<- FALSE
      }
    }
    current[[row]] <<- n_candidate + 1L
    visit(row + 1L, count, score)
    invisible(NULL)
  }
  visit(1L, 0L, 0)
  best_assignment[best_assignment > n_candidate] <- 0L
  list(
    assignment = as.integer(best_assignment),
    match_count = as.integer(best_count),
    total_similarity = best_score
  )
}

v0h_catalog <- function(table, vectors) {
  if (nrow(table) != length(vectors)) stop("catalog rows/vectors mismatch")
  list(table = table, vectors = vectors)
}

v0h_truth_global_catalog <- function(data) {
  true <- data$true_params
  rows <- list()
  vectors <- list()
  index <- 0L
  for (factor in seq_len(ncol(true$a_true))) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      id = sprintf("A0_%d", factor), block = "shared", study = 0L,
      factor = factor, ppi = 1, stringsAsFactors = FALSE
    )
    vectors[[index]] <- true$a_true[, factor]
  }
  for (study in seq_along(true$b_true)) {
    for (factor in seq_len(ncol(true$b_true[[study]]))) {
      index <- index + 1L
      rows[[index]] <- data.frame(
        id = sprintf("B%d0_%d", study, factor), block = "specific",
        study = study, factor = factor, ppi = 1,
        stringsAsFactors = FALSE
      )
      vectors[[index]] <- true$b_true[[study]][, factor]
    }
  }
  v0h_catalog(do.call(rbind, rows), vectors)
}

v0h_candidate_global_catalog <- function(fit, scope, factor_ppi = 0.5) {
  rows <- list()
  vectors <- list()
  index <- 0L
  add_block <- function(loading, ppi, block, study) {
    loading <- as.matrix(loading)
    keep <- if (scope == "ppi_selected") {
      which(as.numeric(ppi) >= factor_ppi)
    } else seq_len(ncol(loading))
    if (!length(keep)) return(invisible(NULL))
    for (factor in keep) {
      index <<- index + 1L
      prefix <- if (block == "shared") "A" else sprintf("B%d", study)
      rows[[index]] <<- data.frame(
        id = sprintf("%s_%d", prefix, factor), block = block,
        study = as.integer(study), factor = as.integer(factor),
        ppi = as.numeric(ppi[[factor]]), stringsAsFactors = FALSE
      )
      vectors[[index]] <<- loading[, factor]
    }
    invisible(NULL)
  }
  add_block(fit$mu_q_a_hat, fit$factor_ppi_shared, "shared", 0L)
  for (study in seq_len(fit$S)) {
    add_block(
      fit$mu_q_b_specific_hat[[study]],
      fit$factor_ppi_specific[[study]], "specific", study
    )
  }
  if (!length(rows)) {
    table <- data.frame(
      id = character(), block = character(), study = integer(),
      factor = integer(), ppi = numeric(), stringsAsFactors = FALSE
    )
  } else table <- do.call(rbind, rows)
  v0h_catalog(table, vectors)
}

v0h_subset_catalog <- function(catalog, keep) {
  keep <- as.integer(keep)
  v0h_catalog(
    catalog$table[keep, , drop = FALSE],
    catalog$vectors[keep]
  )
}

v0h_catalog_similarity <- function(truth, candidate) {
  result <- matrix(0, nrow = nrow(truth$table), ncol = nrow(candidate$table))
  if (!length(result)) return(result)
  for (row in seq_len(nrow(result))) {
    for (column in seq_len(ncol(result))) {
      result[row, column] <- v0h_safe_cosine(
        truth$vectors[[row]], candidate$vectors[[column]]
      )
    }
  }
  result
}

v0h_correct_block <- function(truth_row, candidate_row, global = TRUE) {
  if (truth_row$block != candidate_row$block) return(FALSE)
  if (truth_row$block == "specific") {
    return(as.integer(truth_row$study) == as.integer(candidate_row$study))
  }
  TRUE
}

v0h_classify_catalog_match <- function(
    truth, candidate, scope, matching_unit, unit_study = 0L,
    loading_threshold = 0.8, duplicate_threshold = 0.8) {
  similarity <- v0h_catalog_similarity(truth, candidate)
  optimum <- v0h_one_to_one_match(similarity, loading_threshold)
  truth_rows <- vector("list", nrow(truth$table))
  for (row in seq_len(nrow(truth$table))) {
    candidate_index <- optimum$assignment[[row]]
    matched <- candidate_index > 0L
    correct <- matched && v0h_correct_block(
      truth$table[row, , drop = FALSE],
      candidate$table[candidate_index, , drop = FALSE]
    )
    status <- if (!matched) "missing" else if (correct) "correct" else "misplaced"
    truth_rows[[row]] <- data.frame(
      scope = scope, matching_unit = matching_unit,
      unit_study = as.integer(unit_study),
      true_id = truth$table$id[[row]],
      true_block = truth$table$block[[row]],
      true_study = truth$table$study[[row]],
      true_factor = truth$table$factor[[row]],
      candidate_id = if (matched) candidate$table$id[[candidate_index]] else NA_character_,
      candidate_block = if (matched) candidate$table$block[[candidate_index]] else NA_character_,
      candidate_study = if (matched) candidate$table$study[[candidate_index]] else NA_integer_,
      candidate_factor = if (matched) candidate$table$factor[[candidate_index]] else NA_integer_,
      candidate_ppi = if (matched) candidate$table$ppi[[candidate_index]] else NA_real_,
      loading_abs_cosine = if (matched) similarity[row, candidate_index] else 0,
      matched = matched, status = status,
      stringsAsFactors = FALSE
    )
  }
  matched_candidate <- optimum$assignment[optimum$assignment > 0L]
  matched_truth <- which(optimum$assignment > 0L)
  candidate_rows <- vector("list", nrow(candidate$table))
  for (column in seq_len(nrow(candidate$table))) {
    truth_index <- which(optimum$assignment == column)
    matched <- length(truth_index) == 1L
    if (matched) {
      correct <- v0h_correct_block(
        truth$table[truth_index, , drop = FALSE],
        candidate$table[column, , drop = FALSE]
      )
      status <- if (correct) "correct" else "misplaced"
      maximum_matched_similarity <- similarity[truth_index, column]
      matched_truth_id <- truth$table$id[[truth_index]]
    } else {
      maximum_matched_similarity <- if (length(matched_truth)) {
        max(similarity[matched_truth, column], 0)
      } else 0
      status <- if (maximum_matched_similarity >= duplicate_threshold) {
        "duplicate"
      } else "extra"
      matched_truth_id <- NA_character_
    }
    candidate_rows[[column]] <- data.frame(
      scope = scope, matching_unit = matching_unit,
      unit_study = as.integer(unit_study),
      candidate_id = candidate$table$id[[column]],
      candidate_block = candidate$table$block[[column]],
      candidate_study = candidate$table$study[[column]],
      candidate_factor = candidate$table$factor[[column]],
      candidate_ppi = candidate$table$ppi[[column]],
      matched_truth_id = matched_truth_id,
      maximum_matched_truth_abs_cosine = maximum_matched_similarity,
      matched = matched, status = status,
      stringsAsFactors = FALSE
    )
  }
  empty_candidate <- data.frame(
    scope = character(), matching_unit = character(), unit_study = integer(),
    candidate_id = character(), candidate_block = character(),
    candidate_study = integer(), candidate_factor = integer(),
    candidate_ppi = numeric(), matched_truth_id = character(),
    maximum_matched_truth_abs_cosine = numeric(), matched = logical(),
    status = character(), stringsAsFactors = FALSE
  )
  list(
    similarity = similarity,
    truth = do.call(rbind, truth_rows),
    candidate = if (length(candidate_rows)) {
      do.call(rbind, candidate_rows)
    } else empty_candidate,
    objective = data.frame(
      scope = scope, matching_unit = matching_unit,
      unit_study = as.integer(unit_study),
      match_count = optimum$match_count,
      total_similarity = optimum$total_similarity,
      stringsAsFactors = FALSE
    )
  )
}

v0h_study_catalogs <- function(truth_global, candidate_global, study) {
  truth_keep <- which(
    truth_global$table$block == "shared" |
      (truth_global$table$block == "specific" &
         truth_global$table$study == study)
  )
  candidate_keep <- which(
    candidate_global$table$block == "shared" |
      (candidate_global$table$block == "specific" &
         candidate_global$table$study == study)
  )
  list(
    truth = v0h_subset_catalog(truth_global, truth_keep),
    candidate = v0h_subset_catalog(candidate_global, candidate_keep)
  )
}

v0h_selected_columns <- function(value, ppi, scope, threshold) {
  keep <- if (scope == "ppi_selected") {
    which(as.numeric(ppi) >= threshold)
  } else seq_len(ncol(as.matrix(value)))
  v0h_subset_columns(value, keep)
}

v0h_block_subspace_metrics <- function(fit, data, scope, thresholds) {
  true <- data$true_params
  shared_estimate <- v0h_selected_columns(
    fit$mu_q_a_hat, fit$factor_ppi_shared,
    scope, thresholds$factor_ppi
  )
  definitions <- list(list(
    metric = "R_A", metric_class = "correct_block",
    truth_block = "shared", truth_study = 0L,
    estimated_block = "shared", estimated_study = 0L,
    truth = true$a_true, estimate = shared_estimate
  ))
  for (study in seq_len(fit$S)) {
    specific_estimate <- v0h_selected_columns(
      fit$mu_q_b_specific_hat[[study]],
      fit$factor_ppi_specific[[study]],
      scope, thresholds$factor_ppi
    )
    definitions <- c(definitions, list(
      list(
        metric = sprintf("R_B%d", study), metric_class = "correct_block",
        truth_block = "specific", truth_study = study,
        estimated_block = "specific", estimated_study = study,
        truth = true$b_true[[study]], estimate = specific_estimate
      ),
      list(
        metric = sprintf("L_B%d_to_A", study),
        metric_class = "cross_block_leakage",
        truth_block = "specific", truth_study = study,
        estimated_block = "shared", estimated_study = 0L,
        truth = true$b_true[[study]], estimate = shared_estimate
      ),
      list(
        metric = sprintf("L_A_to_B%d", study),
        metric_class = "cross_block_leakage",
        truth_block = "shared", truth_study = 0L,
        estimated_block = "specific", estimated_study = study,
        truth = true$a_true, estimate = specific_estimate
      )
    ))
  }
  do.call(rbind, lapply(definitions, function(definition) {
    metric <- v0h_subspace_metrics(
      definition$truth, definition$estimate,
      relative_tolerance = thresholds$svd_relative_tolerance
    )
    data.frame(
      scope = scope,
      metric = definition$metric,
      metric_class = definition$metric_class,
      truth_block = definition$truth_block,
      truth_study = definition$truth_study,
      estimated_block = definition$estimated_block,
      estimated_study = definition$estimated_study,
      as.data.frame(metric, stringsAsFactors = FALSE),
      stringsAsFactors = FALSE
    )
  }))
}

v0h_truth_paths_by_study <- function(data) {
  true <- data$true_params
  bind_paths <- function(values, study) {
    observation_counts <- lengths(data$time_obs[[study]])
    if (!length(values) || all(vapply(values, is.null, logical(1)))) {
      return(v0h_empty_matrix(sum(observation_counts)))
    }
    if (length(values) != length(observation_counts)) {
      stop("truth paths and observed subjects have incompatible lengths")
    }
    do.call(rbind, lapply(seq_along(values), function(subject) {
      value <- values[[subject]]
      if (is.null(value)) {
        v0h_empty_matrix(observation_counts[[subject]])
      } else {
        as.matrix(value)
      }
    }))
  }
  list(
    shared = lapply(seq_len(true$S), function(study) {
      bind_paths(true$f_true_values[[study]], study)
    }),
    specific = lapply(seq_len(true$S), function(study) {
      bind_paths(true$g_true_values[[study]], study)
    })
  )
}

v0h_interpolate_matrix <- function(grid, value, target) {
  value <- as.matrix(value)
  if (!ncol(value)) return(v0h_empty_matrix(length(target)))
  vapply(seq_len(ncol(value)), function(column) {
    as.numeric(stats::approx(
      grid, value[, column], xout = target,
      rule = 2, ties = mean
    )$y)
  }, numeric(length(target)))
}

v0h_true_function_and_score <- function(data, study, block, factor) {
  if (block == "shared") {
    list(
      functions = data$true_params$phi_dense_list[[factor]],
      scores = data$true_params$zeta_fpca_true[[study]][[factor]]
    )
  } else {
    list(
      functions = data$true_params$psi_dense_list[[study]][[factor]],
      scores = data$true_params$xi_fpca_true[[study]][[factor]]
    )
  }
}

v0h_estimated_function_and_score <- function(fit, study, block, factor) {
  if (block == "shared") {
    offset <- sum(fit$n_s[seq_len(study - 1L)])
    rows <- offset + seq_len(fit$n_s[[study]])
    effective <- as.integer(fit$list_effective_M[[factor]])
    functions <- fit$list_Phi_hat[[factor]]
    scores <- fit$list_Zeta_hat[[factor]][rows, , drop = FALSE]
  } else {
    effective <- as.integer(fit$list_effective_M_spec[[study]][[factor]])
    functions <- fit$list_Phi_hat_spec[[study]][[factor]]
    scores <- fit$list_Zeta_hat_spec[[study]][[factor]]
  }
  effective <- max(0L, min(
    effective, ncol(functions), ncol(scores)
  ))
  keep <- if (effective) seq_len(effective) else integer()
  list(
    functions = v0h_subset_columns(functions, keep),
    scores = v0h_subset_columns(scores, keep),
    effective_M = effective
  )
}

v0h_weighted_cosine <- function(x, y, weights) {
  v0h_safe_cosine(sqrt(weights) * x, sqrt(weights) * y)
}

v0h_component_recovery <- function(
    true_object, estimated_object, loading_sign, grid_true, grid_estimate,
    metadata) {
  truth_function <- as.matrix(true_object$functions)
  truth_score <- as.matrix(true_object$scores)
  estimated_function <- v0h_interpolate_matrix(
    grid_estimate, estimated_object$functions, grid_true
  )
  estimated_score <- as.matrix(estimated_object$scores)
  weights <- v0h_trapezoid_weights(grid_true)
  similarity <- matrix(
    0, nrow = ncol(truth_function), ncol = ncol(estimated_function)
  )
  if (length(similarity)) {
    for (row in seq_len(nrow(similarity))) {
      for (column in seq_len(ncol(similarity))) {
        similarity[row, column] <- v0h_weighted_cosine(
          truth_function[, row], estimated_function[, column], weights
        )
      }
    }
  }
  optimum <- v0h_one_to_one_match(similarity, threshold = 0)
  rows <- list()
  index <- 0L
  score_abs_by_truth <- rep(0, ncol(truth_function))
  function_abs_by_truth <- rep(0, ncol(truth_function))
  for (component in seq_len(ncol(truth_function))) {
    estimated_component <- optimum$assignment[[component]]
    matched <- estimated_component > 0L
    if (matched) {
      raw_inner <- sum(
        weights * truth_function[, component] *
          estimated_function[, estimated_component]
      )
      function_sign <- if (!is.finite(raw_inner) || raw_inner == 0) {
        1
      } else sign(raw_inner)
      aligned_score <- loading_sign * function_sign *
        estimated_score[, estimated_component]
      signed_score_cor <- v0h_safe_cor(
        aligned_score, truth_score[, component]
      )
      score_abs_by_truth[[component]] <- abs(signed_score_cor)
      function_abs_by_truth[[component]] <-
        similarity[component, estimated_component]
    } else {
      function_sign <- NA_real_
      signed_score_cor <- NA_real_
    }
    index <- index + 1L
    rows[[index]] <- cbind(
      metadata,
      data.frame(
        component_record = "true_component",
        true_component = component,
        estimated_component = if (matched) estimated_component else NA_integer_,
        component_abs_inner_product = if (matched) {
          similarity[component, estimated_component]
        } else 0,
        loading_sign = loading_sign,
        function_sign = function_sign,
        score_signed_cor = signed_score_cor,
        score_abs_cor = if (matched) abs(signed_score_cor) else 0,
        component_status = if (matched) "matched" else "missing",
        stringsAsFactors = FALSE
      )
    )
  }
  extra <- setdiff(
    seq_len(ncol(estimated_function)),
    optimum$assignment[optimum$assignment > 0L]
  )
  if (length(extra)) {
    for (estimated_component in extra) {
      index <- index + 1L
      rows[[index]] <- cbind(
        metadata,
        data.frame(
          component_record = "estimated_component",
          true_component = NA_integer_,
          estimated_component = estimated_component,
          component_abs_inner_product = 0,
          loading_sign = loading_sign,
          function_sign = NA_real_,
          score_signed_cor = NA_real_, score_abs_cor = NA_real_,
          component_status = "extra", stringsAsFactors = FALSE
        )
      )
    }
  }
  list(
    rows = if (length(rows)) do.call(rbind, rows) else data.frame(),
    function_abs_by_truth = function_abs_by_truth,
    score_abs_by_truth = score_abs_by_truth,
    matched_components = optimum$match_count,
    interpolated_function = estimated_function,
    weights = weights
  )
}

v0h_functional_recovery <- function(
    fit, data, study_matches, reported_paths, thresholds) {
  truth_paths <- v0h_truth_paths_by_study(data)
  factor_rows <- list()
  component_rows <- list()
  factor_index <- 0L
  component_index <- 0L
  for (row in seq_len(nrow(study_matches))) {
    match_row <- study_matches[row, , drop = FALSE]
    factor_index <- factor_index + 1L
    metadata <- data.frame(
      scope = match_row$scope,
      study = match_row$unit_study,
      true_id = match_row$true_id,
      true_role = match_row$true_block,
      true_factor = match_row$true_factor,
      candidate_id = match_row$candidate_id,
      estimated_block = match_row$candidate_block,
      estimated_factor = match_row$candidate_factor,
      loading_structure_status = match_row$status,
      stringsAsFactors = FALSE
    )
    if (!isTRUE(match_row$matched[[1L]])) {
      factor_rows[[factor_index]] <- cbind(
        metadata,
        data.frame(
          matched_by_loading = FALSE,
          loading_abs_cosine = match_row$loading_abs_cosine,
          loading_sign = NA_real_, trajectory_signed_cor = NA_real_,
          trajectory_abs_cor = NA_real_, trajectory_nrmse = NA_real_,
          true_M = if (match_row$true_block == "shared") {
            ncol(data$true_params$phi_dense_list[[match_row$true_factor]])
          } else {
            ncol(data$true_params$psi_dense_list[[match_row$unit_study]][[
              match_row$true_factor
            ]])
          },
          estimated_effective_M = 0L,
          function_subspace_recall = 0,
          function_subspace_purity = 0,
          function_min_principal_cosine_sq = 0,
          function_component_abs_inner_mean = 0,
          function_component_abs_inner_min = 0,
          score_abs_cor_mean = 0, score_abs_cor_min = 0,
          strict_joint_recovery = FALSE,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    study <- match_row$unit_study[[1L]]
    true_block <- match_row$true_block[[1L]]
    true_factor <- match_row$true_factor[[1L]]
    estimated_block <- match_row$candidate_block[[1L]]
    estimated_factor <- match_row$candidate_factor[[1L]]
    true_loading <- if (true_block == "shared") {
      data$true_params$a_true[, true_factor]
    } else data$true_params$b_true[[study]][, true_factor]
    estimated_loading <- if (estimated_block == "shared") {
      fit$mu_q_a_hat[, estimated_factor]
    } else fit$mu_q_b_specific_hat[[study]][, estimated_factor]
    loading_sign <- v0h_loading_sign(estimated_loading, true_loading)
    estimated_path <- reported_paths[[estimated_block]][[study]][,
      estimated_factor]
    true_path <- truth_paths[[true_block]][[study]][, true_factor]
    aligned_path <- loading_sign * estimated_path
    trajectory_cor <- v0h_safe_cor(aligned_path, true_path)
    true_object <- v0h_true_function_and_score(
      data, study, true_block, true_factor
    )
    estimated_object <- v0h_estimated_function_and_score(
      fit, study, estimated_block, estimated_factor
    )
    component <- v0h_component_recovery(
      true_object, estimated_object, loading_sign,
      data$true_params$t_grid_dense, fit$time_g, metadata
    )
    if (nrow(component$rows)) {
      component_index <- component_index + 1L
      component_rows[[component_index]] <- component$rows
    }
    function_metrics <- v0h_subspace_metrics(
      true_object$functions, component$interpolated_function,
      component$weights, thresholds$svd_relative_tolerance
    )
    factor_rows[[factor_index]] <- cbind(
      metadata,
      data.frame(
        matched_by_loading = TRUE,
        loading_abs_cosine = match_row$loading_abs_cosine,
        loading_sign = loading_sign,
        trajectory_signed_cor = trajectory_cor,
        trajectory_abs_cor = abs(trajectory_cor),
        trajectory_nrmse = v0h_nrmse(aligned_path, true_path),
        true_M = ncol(true_object$functions),
        estimated_effective_M = estimated_object$effective_M,
        function_subspace_recall = function_metrics$subspace_recall,
        function_subspace_purity = function_metrics$subspace_purity,
        function_min_principal_cosine_sq =
          function_metrics$principal_cosine_sq_min,
        function_component_abs_inner_mean =
          mean(component$function_abs_by_truth),
        function_component_abs_inner_min =
          min(component$function_abs_by_truth),
        score_abs_cor_mean = mean(component$score_abs_by_truth),
        score_abs_cor_min = min(component$score_abs_by_truth),
        strict_joint_recovery =
          match_row$loading_abs_cosine >= thresholds$loading_abs_cosine &&
          abs(trajectory_cor) >= thresholds$trajectory_abs_cor,
        stringsAsFactors = FALSE
      )
    )
  }
  empty_components <- data.frame(
    scope = character(), study = integer(), true_id = character(),
    true_role = character(), true_factor = integer(),
    candidate_id = character(), estimated_block = character(),
    estimated_factor = integer(), loading_structure_status = character(),
    component_record = character(),
    true_component = integer(), estimated_component = integer(),
    component_abs_inner_product = numeric(), loading_sign = numeric(),
    function_sign = numeric(), score_signed_cor = numeric(),
    score_abs_cor = numeric(), component_status = character(),
    stringsAsFactors = FALSE
  )
  list(
    factors = do.call(rbind, factor_rows),
    components = if (length(component_rows)) {
      do.call(rbind, component_rows)
    } else empty_components
  )
}

v0h_extract_metric <- function(block_table, scope, metric, field) {
  row <- block_table$scope == scope & block_table$metric == metric
  if (sum(row) != 1L) return(NA_real_)
  as.numeric(block_table[[field]][row])
}

v0h_count_status <- function(table, scope, status) {
  sum(table$scope == scope & table$status == status)
}

v0h_last_finite <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value)) tail(value, 1L) else NA_real_
}

v0h_evaluate_fit <- function(fit, data, thresholds = NULL) {
  if (is.null(thresholds)) thresholds <- v0h_default_thresholds()
  required_v0g <- c("v0g_reported_paths_by_study", "v0g_prediction_diagnostics")
  if (!all(vapply(required_v0g, exists, logical(1), mode = "function"))) {
    stop("source v0g_evaluation.R before calling v0h_evaluate_fit")
  }
  scopes <- c("ppi_selected", "all_candidates")
  truth_global <- v0h_truth_global_catalog(data)
  reported_paths <- v0g_reported_paths_by_study(fit, data)
  prediction <- v0g_prediction_diagnostics(fit, data, reported_paths)
  block_rows <- list()
  global_truth_rows <- list()
  global_candidate_rows <- list()
  study_truth_rows <- list()
  study_candidate_rows <- list()
  objective_rows <- list()
  counter <- 0L
  for (scope in scopes) {
    candidate_global <- v0h_candidate_global_catalog(
      fit, scope, thresholds$factor_ppi
    )
    counter <- counter + 1L
    block_rows[[counter]] <- v0h_block_subspace_metrics(
      fit, data, scope, thresholds
    )
    global_match <- v0h_classify_catalog_match(
      truth_global, candidate_global, scope, "global_direction", 0L,
      thresholds$loading_abs_cosine,
      thresholds$duplicate_abs_cosine
    )
    global_truth_rows[[length(global_truth_rows) + 1L]] <- global_match$truth
    global_candidate_rows[[length(global_candidate_rows) + 1L]] <-
      global_match$candidate
    objective_rows[[length(objective_rows) + 1L]] <- global_match$objective
    for (study in seq_len(fit$S)) {
      catalogs <- v0h_study_catalogs(truth_global, candidate_global, study)
      study_match <- v0h_classify_catalog_match(
        catalogs$truth, catalogs$candidate, scope, "study_role", study,
        thresholds$loading_abs_cosine,
        thresholds$duplicate_abs_cosine
      )
      study_truth_rows[[length(study_truth_rows) + 1L]] <- study_match$truth
      study_candidate_rows[[length(study_candidate_rows) + 1L]] <-
        study_match$candidate
      objective_rows[[length(objective_rows) + 1L]] <- study_match$objective
    }
  }
  block <- do.call(rbind, block_rows)
  global_truth <- do.call(rbind, global_truth_rows)
  global_candidate <- do.call(rbind, global_candidate_rows)
  study_truth <- do.call(rbind, study_truth_rows)
  study_candidate <- do.call(rbind, study_candidate_rows)
  functional <- v0h_functional_recovery(
    fit, data, study_truth, reported_paths, thresholds
  )
  selected_shared <- sum(fit$factor_ppi_shared >= thresholds$factor_ppi)
  selected_specific <- vapply(
    fit$factor_ppi_specific,
    function(value) sum(value >= thresholds$factor_ppi), integer(1)
  )
  truth_specific_counts <- vapply(
    data$true_params$b_true, ncol, integer(1)
  )
  candidate_specific_counts <- vapply(
    fit$mu_q_b_specific_hat, ncol, integer(1)
  )
  if (length(truth_specific_counts) != data$true_params$S ||
      length(candidate_specific_counts) != fit$S) {
    stop("study-specific loading blocks do not match the number of studies")
  }
  primary <- "ppi_selected"
  primary_functional <- functional$factors[
    functional$factors$scope == primary, , drop = FALSE
  ]
  summary <- data.frame(
    truth_global_direction_count = nrow(truth_global$table),
    truth_study_role_count = data$true_params$S *
      data$true_params$L_f + sum(truth_specific_counts),
    candidate_shared = fit$L_f,
    candidate_specific_study_1 = candidate_specific_counts[[1L]],
    candidate_specific_study_2 = candidate_specific_counts[[2L]],
    selected_shared = selected_shared,
    selected_specific_study_1 = selected_specific[[1L]],
    selected_specific_study_2 = selected_specific[[2L]],
    selected_counts = paste(c(selected_shared, selected_specific), collapse = ","),
    selected_loading_rank_shared = v0h_extract_metric(
      block, primary, "R_A", "estimate_rank"
    ),
    selected_loading_rank_specific_study_1 = v0h_extract_metric(
      block, primary, "R_B1", "estimate_rank"
    ),
    selected_loading_rank_specific_study_2 = v0h_extract_metric(
      block, primary, "R_B2", "estimate_rank"
    ),
    R_A = v0h_extract_metric(block, primary, "R_A", "subspace_recall"),
    R_B1 = v0h_extract_metric(block, primary, "R_B1", "subspace_recall"),
    R_B2 = v0h_extract_metric(block, primary, "R_B2", "subspace_recall"),
    L_B1_to_A = v0h_extract_metric(
      block, primary, "L_B1_to_A", "subspace_recall"
    ),
    L_B2_to_A = v0h_extract_metric(
      block, primary, "L_B2_to_A", "subspace_recall"
    ),
    L_A_to_B1 = v0h_extract_metric(
      block, primary, "L_A_to_B1", "subspace_recall"
    ),
    L_A_to_B2 = v0h_extract_metric(
      block, primary, "L_A_to_B2", "subspace_recall"
    ),
    global_correct = v0h_count_status(global_truth, primary, "correct"),
    global_misplaced = v0h_count_status(global_truth, primary, "misplaced"),
    global_missing = v0h_count_status(global_truth, primary, "missing"),
    global_extra = v0h_count_status(global_candidate, primary, "extra"),
    global_duplicate = v0h_count_status(global_candidate, primary, "duplicate"),
    study_correct = v0h_count_status(study_truth, primary, "correct"),
    study_misplaced = v0h_count_status(study_truth, primary, "misplaced"),
    study_missing = v0h_count_status(study_truth, primary, "missing"),
    study_candidate_extra = v0h_count_status(
      study_candidate, primary, "extra"
    ),
    study_candidate_duplicate = v0h_count_status(
      study_candidate, primary, "duplicate"
    ),
    strict_joint_recovery = sum(primary_functional$strict_joint_recovery),
    strict_joint_correct_type = sum(
      primary_functional$strict_joint_recovery &
        primary_functional$loading_structure_status == "correct"
    ),
    strict_joint_misplaced = sum(
      primary_functional$strict_joint_recovery &
        primary_functional$loading_structure_status == "misplaced"
    ),
    matched_function_recall_mean = if (any(primary_functional$matched_by_loading)) {
      mean(primary_functional$function_subspace_recall[
        primary_functional$matched_by_loading
      ])
    } else 0,
    matched_trajectory_abs_cor_mean = if (any(primary_functional$matched_by_loading)) {
      mean(primary_functional$trajectory_abs_cor[
        primary_functional$matched_by_loading
      ])
    } else 0,
    matched_score_abs_cor_mean = if (any(primary_functional$matched_by_loading)) {
      mean(primary_functional$score_abs_cor_mean[
        primary_functional$matched_by_loading
      ])
    } else 0,
    correctly_placed_function_recall_mean = if (any(
      primary_functional$loading_structure_status == "correct"
    )) {
      mean(primary_functional$function_subspace_recall[
        primary_functional$loading_structure_status == "correct"
      ])
    } else 0,
    correctly_placed_trajectory_abs_cor_mean = if (any(
      primary_functional$loading_structure_status == "correct"
    )) {
      mean(primary_functional$trajectory_abs_cor[
        primary_functional$loading_structure_status == "correct"
      ])
    } else 0,
    correctly_placed_score_abs_cor_mean = if (any(
      primary_functional$loading_structure_status == "correct"
    )) {
      mean(primary_functional$score_abs_cor_mean[
        primary_functional$loading_structure_status == "correct"
      ])
    } else 0,
    shared_effective_M = paste(unlist(fit$list_effective_M), collapse = ";"),
    specific_effective_M_study_1 = paste(
      unlist(fit$list_effective_M_spec[[1L]]), collapse = ";"
    ),
    specific_effective_M_study_2 = paste(
      unlist(fit$list_effective_M_spec[[2L]]), collapse = ";"
    ),
    signal_nrmse_all_candidates = prediction[["signal_nrmse_all_candidates"]],
    signal_nrmse_ppi_selected = prediction[["signal_nrmse_ppi_selected"]],
    signal_nrmse_mean_only = prediction[["signal_nrmse_mean_only"]],
    final_elbo = v0h_last_finite(fit$ELBO),
    iterations = as.integer(fit$i_iter),
    t1_sweeps = as.integer(fit$t1_sweeps),
    continuation_used = isTRUE(fit$continuation_diagnostics$used),
    stringsAsFactors = FALSE
  )
  list(
    thresholds = thresholds,
    block_subspace = block,
    global_direction_matches = global_truth,
    global_candidate_classification = global_candidate,
    study_role_matches = study_truth,
    study_candidate_classification = study_candidate,
    matching_objectives = do.call(rbind, objective_rows),
    functional_recovery = functional$factors,
    fpca_component_recovery = functional$components,
    summary = summary
  )
}
