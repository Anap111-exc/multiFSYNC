s6c_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

s6c_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

s6c_write_lines <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(as.character(value)), path, useBytes = TRUE)
  invisible(path)
}

s6c_sha256 <- function(path) {
  s6c_assert(requireNamespace("digest", quietly = TRUE),
             "Package digest is required for the truth-free freeze hash.")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

s6c_no_truth_fields <- function(value) {
  !any(grepl("(^|_)(truth|true)(_|$)", names(value), ignore.case = TRUE))
}

s6c_matrix_basis <- function(value, tolerance = 1e-10) {
  value <- as.matrix(value)
  if (!nrow(value) || !ncol(value) || max(abs(value)) <= tolerance) {
    return(matrix(numeric(), nrow = nrow(value), ncol = 0L))
  }
  decomposition <- svd(value, nu = min(dim(value)), nv = 0L)
  rank <- sum(decomposition$d > tolerance * decomposition$d[[1L]])
  decomposition$u[, seq_len(rank), drop = FALSE]
}

s6c_subspace_similarity <- function(left, right) {
  left_basis <- s6c_matrix_basis(left)
  right_basis <- s6c_matrix_basis(right)
  left_rank <- ncol(left_basis)
  right_rank <- ncol(right_basis)
  if (!left_rank && !right_rank) {
    return(c(mean_squared_canonical_correlation = 1,
             minimum_canonical_correlation = 1,
             rank_left = 0, rank_right = 0))
  }
  if (!left_rank || !right_rank) {
    return(c(mean_squared_canonical_correlation = 0,
             minimum_canonical_correlation = 0,
             rank_left = left_rank, rank_right = right_rank))
  }
  singular_values <- svd(crossprod(left_basis, right_basis),
                         nu = 0L, nv = 0L)$d
  singular_values <- pmin(1, pmax(0, singular_values))
  target <- max(left_rank, right_rank)
  padded <- c(singular_values, rep(0, target - length(singular_values)))
  c(mean_squared_canonical_correlation = mean(padded^2),
    minimum_canonical_correlation = min(padded),
    rank_left = left_rank, rank_right = right_rank)
}

s6c_factor_path <- function(scores, functions) {
  scores <- as.matrix(scores)
  functions <- as.matrix(functions)
  s6c_assert(ncol(scores) == ncol(functions),
             "Score and function component counts differ.")
  scores %*% t(functions)
}

s6c_role_descriptor <- function(loadings, score_list, function_list,
                                 factor_ppi, effective_m) {
  loadings <- as.matrix(loadings)
  factors <- ncol(loadings)
  s6c_assert(factors >= 1L && length(score_list) == factors &&
               length(function_list) == factors,
             "Role dimensions are inconsistent.")
  paths <- lapply(seq_len(factors), function(index) {
    as.numeric(s6c_factor_path(score_list[[index]], function_list[[index]]))
  })
  path_matrix <- do.call(rbind, paths)
  candidate_energy <- vapply(seq_len(factors), function(index) {
    sum(loadings[, index]^2) * sum(path_matrix[index, ]^2)
  }, numeric(1L))
  total_candidate_energy <- sum(candidate_energy)
  list(
    loadings = loadings,
    paths = path_matrix,
    factor_ppi = as.numeric(factor_ppi),
    effective_m = as.integer(effective_m),
    candidate_energy_share = if (total_candidate_energy > 0) {
      candidate_energy / total_candidate_energy
    } else rep(0, factors)
  )
}

s6c_fit_descriptor <- function(fit, record, cohort_id, fit_path,
                                terminal_path) {
  required <- c(
    "S", "n_s", "p", "L_f", "L_s_by_study", "time_g",
    "mu_q_a_hat", "factor_ppi_shared", "list_Zeta_hat", "list_Phi_hat",
    "list_effective_M", "mu_q_b_specific_hat", "factor_ppi_specific",
    "list_Zeta_hat_spec", "list_Phi_hat_spec", "list_effective_M_spec"
  )
  s6c_assert(all(required %in% names(fit)),
             "Fit lacks fields required by Stage 6C.")
  s6c_assert(s6c_no_truth_fields(fit),
             "A fit contains a forbidden truth-named top-level field.")
  s6c_assert(isTRUE(record$objective_eligible),
             "A terminal record is not objective eligible.")
  s6c_assert(isFALSE(record$truth_used_for_fit_stopping_or_selection),
             "A terminal record indicates truth use.")
  s6c_assert(isFALSE(record$continuation_used),
             "Stage 6C excludes continuation endpoints.")
  s6c_assert(identical(as.character(record$fit_id), basename(dirname(fit_path))),
             "Fit directory and terminal fit_id differ.")

  shared <- s6c_role_descriptor(
    fit$mu_q_a_hat,
    fit$list_Zeta_hat,
    fit$list_Phi_hat,
    fit$factor_ppi_shared,
    fit$list_effective_M
  )
  specific <- lapply(seq_len(as.integer(fit$S)), function(study) {
    s6c_role_descriptor(
      fit$mu_q_b_specific_hat[[study]],
      fit$list_Zeta_hat_spec[[study]],
      fit$list_Phi_hat_spec[[study]],
      fit$factor_ppi_specific[[study]],
      fit$list_effective_M_spec[[study]]
    )
  })
  roles <- c(list(A = shared), setNames(specific, paste0("B", seq_along(specific))))
  active_counts <- vapply(roles, function(role) {
    sum(role$factor_ppi >= 0.5)
  }, integer(1L))
  effective_m <- vapply(roles, function(role) {
    if (length(role$effective_m)) max(role$effective_m) else 0L
  }, integer(1L))
  endpoint <- data.frame(
    cohort_id = cohort_id,
    config_group_id = paste(cohort_id, record$fit_config_id, sep = "::"),
    fit_id = as.character(record$fit_id),
    fit_config_id = as.character(record$fit_config_id),
    seed_index = as.integer(record$seed_index),
    fit_seed = as.integer(record$fit_seed),
    final_elbo = as.numeric(record$final_elbo),
    practical_converged = isTRUE(record$practical_converged),
    fit_L_f = as.integer(fit$L_f),
    fit_L_s_1 = as.integer(fit$L_s_by_study[[1L]]),
    fit_L_s_2 = as.integer(fit$L_s_by_study[[2L]]),
    ppi_active_shared = active_counts[["A"]],
    ppi_active_specific_1 = active_counts[["B1"]],
    ppi_active_specific_2 = active_counts[["B2"]],
    maximum_effective_M_shared = effective_m[["A"]],
    maximum_effective_M_specific_1 = effective_m[["B1"]],
    maximum_effective_M_specific_2 = effective_m[["B2"]],
    truth_used = FALSE,
    stringsAsFactors = FALSE
  )
  list(endpoint = endpoint, roles = roles, time_g = as.numeric(fit$time_g),
       input = data.frame(
         cohort_id = cohort_id,
         fit_id = as.character(record$fit_id),
         input_type = c("fit", "terminal"),
         path = c(fit_path, terminal_path),
         bytes = c(file.info(fit_path)$size, file.info(terminal_path)$size),
         stringsAsFactors = FALSE
       ))
}

s6c_contribution_pair <- function(left, right) {
  s6c_assert(nrow(left$loadings) == nrow(right$loadings) &&
               ncol(left$paths) == ncol(right$paths),
             "Role contribution dimensions differ.")
  cross_loading <- crossprod(left$loadings, right$loadings)
  cross_path <- tcrossprod(left$paths, right$paths)
  inner <- sum(cross_loading * cross_path)
  left_sq <- sum(crossprod(left$loadings) * tcrossprod(left$paths))
  right_sq <- sum(crossprod(right$loadings) * tcrossprod(right$paths))
  left_sq <- max(0, left_sq)
  right_sq <- max(0, right_sq)
  if (left_sq <= 1e-20 && right_sq <= 1e-20) {
    cosine <- 1
    amplitude_ratio <- 1
    relative_difference <- 0
  } else if (left_sq <= 1e-20 || right_sq <= 1e-20) {
    cosine <- 0
    amplitude_ratio <- 0
    relative_difference <- 1
  } else {
    cosine <- inner / sqrt(left_sq * right_sq)
    cosine <- min(1, max(-1, cosine))
    amplitude_ratio <- min(sqrt(left_sq), sqrt(right_sq)) /
      max(sqrt(left_sq), sqrt(right_sq))
    denominator <- left_sq + right_sq
    relative_difference <- sqrt(max(0, (denominator - 2 * inner) / denominator))
  }
  c(contribution_cosine = cosine,
    contribution_similarity = max(0, cosine),
    contribution_amplitude_ratio = amplitude_ratio,
    contribution_relative_difference = relative_difference)
}

s6c_pair_row <- function(left, right) {
  roles <- c("A", "B1", "B2")
  values <- list()
  for (role in roles) {
    contribution <- s6c_contribution_pair(left$roles[[role]],
                                          right$roles[[role]])
    loading <- s6c_subspace_similarity(left$roles[[role]]$loadings,
                                       right$roles[[role]]$loadings)
    values[[paste0(role, "_contribution_cosine")]] <-
      contribution[["contribution_cosine"]]
    values[[paste0(role, "_contribution_similarity")]] <-
      contribution[["contribution_similarity"]]
    values[[paste0(role, "_contribution_amplitude_ratio")]] <-
      contribution[["contribution_amplitude_ratio"]]
    values[[paste0(role, "_contribution_relative_difference")]] <-
      contribution[["contribution_relative_difference"]]
    values[[paste0(role, "_loading_subspace_mscc")]] <-
      loading[["mean_squared_canonical_correlation"]]
    values[[paste0(role, "_loading_subspace_min_cc")]] <-
      loading[["minimum_canonical_correlation"]]
  }
  contribution_similarity <- unlist(values[paste0(
    roles, "_contribution_similarity"
  )], use.names = FALSE)
  loading_similarity <- unlist(values[paste0(
    roles, "_loading_subspace_mscc"
  )], use.names = FALSE)
  amplitude_similarity <- unlist(values[paste0(
    roles, "_contribution_amplitude_ratio"
  )], use.names = FALSE)
  data.frame(
    cohort_id = left$endpoint$cohort_id,
    config_group_id = left$endpoint$config_group_id,
    fit_config_id = left$endpoint$fit_config_id,
    fit_id_left = left$endpoint$fit_id,
    fit_id_right = right$endpoint$fit_id,
    seed_index_left = left$endpoint$seed_index,
    seed_index_right = right$endpoint$seed_index,
    as.data.frame(values, check.names = FALSE),
    contribution_similarity_min = min(contribution_similarity),
    contribution_similarity_mean = mean(contribution_similarity),
    loading_subspace_similarity_min = min(loading_similarity),
    loading_subspace_similarity_mean = mean(loading_similarity),
    amplitude_similarity_min = min(amplitude_similarity),
    structural_similarity_min = min(c(contribution_similarity,
                                      loading_similarity)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

s6c_read_cohort <- function(root, cohort_id, expected_fits) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  fit_root <- file.path(root, "fits")
  directories <- sort(list.dirs(fit_root, recursive = FALSE,
                                full.names = FALSE))
  directories <- directories[nzchar(directories)]
  s6c_assert(length(directories) == expected_fits,
             paste0(cohort_id, " expected ", expected_fits,
                    " fit directories, found ", length(directories), "."))
  descriptors <- setNames(vector("list", length(directories)), directories)
  inputs <- vector("list", length(directories))
  for (index in seq_along(directories)) {
    fit_id <- directories[[index]]
    fit_path <- normalizePath(file.path(fit_root, fit_id, "fit.rds"),
                              winslash = "/", mustWork = TRUE)
    terminal_path <- normalizePath(
      file.path(fit_root, fit_id, "terminal_record.rds"),
      winslash = "/", mustWork = TRUE
    )
    fit <- readRDS(fit_path)
    record <- readRDS(terminal_path)
    descriptors[[fit_id]] <- s6c_fit_descriptor(
      fit, record, cohort_id, fit_path, terminal_path
    )
    inputs[[index]] <- descriptors[[fit_id]]$input
    rm(fit, record)
    invisible(gc(FALSE))
  }
  list(descriptors = descriptors, inputs = do.call(rbind, inputs))
}

s6c_cluster_group <- function(endpoints, pairs, mode, threshold) {
  fit_ids <- endpoints$fit_id
  count <- length(fit_ids)
  similarity_column <- if (mode == "contribution") {
    "contribution_similarity_min"
  } else "structural_similarity_min"
  similarity <- diag(1, count)
  dimnames(similarity) <- list(fit_ids, fit_ids)
  if (count >= 2L) {
    for (row in seq_len(nrow(pairs))) {
      left <- match(pairs$fit_id_left[[row]], fit_ids)
      right <- match(pairs$fit_id_right[[row]], fit_ids)
      similarity[left, right] <- pairs[[similarity_column]][[row]]
      similarity[right, left] <- pairs[[similarity_column]][[row]]
    }
    clustering <- stats::hclust(stats::as.dist(1 - similarity),
                                method = "complete")
    raw_cluster <- stats::cutree(clustering, h = 1 - threshold + 1e-12)
  } else raw_cluster <- setNames(1L, fit_ids)
  raw_cluster <- raw_cluster[fit_ids]
  cluster_order <- unique(raw_cluster[order(fit_ids)])
  cluster <- match(raw_cluster, cluster_order)
  cluster_size <- as.integer(table(cluster)[as.character(cluster)])
  cluster_support <- cluster_size / count
  max_size <- max(cluster_size)
  is_modal <- cluster_size == max_size
  elbo_winner <- endpoints$fit_id[[which.max(endpoints$final_elbo)]]
  representative <- vapply(seq_along(cluster), function(index) {
    rows <- which(cluster == cluster[[index]])
    endpoints$fit_id[rows[[which.max(endpoints$final_elbo[rows])]]]
  }, character(1L))
  data.frame(
    endpoints,
    basin_view = mode,
    similarity_threshold = threshold,
    basin_id = paste0("B", cluster),
    basin_size = cluster_size,
    basin_support = cluster_support,
    is_modal_basin = is_modal,
    basin_representative_fit_id = representative,
    is_basin_representative = endpoints$fit_id == representative,
    raw_elbo_winner_fit_id = elbo_winner,
    is_raw_elbo_winner = endpoints$fit_id == elbo_winner,
    raw_elbo_winner_in_modal_basin = is_modal[
      match(elbo_winner, endpoints$fit_id)
    ],
    stringsAsFactors = FALSE
  )
}

s6c_audit <- function(stage6a_root, stage6b_root, output_root,
                      thresholds = c(0.8, 0.9, 0.95)) {
  s6c_assert(!dir.exists(output_root) ||
               !length(list.files(output_root, all.files = TRUE,
                                  no.. = TRUE)),
             "Stage 6C output directory already contains files.")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  cohorts <- list(
    s6a_full_data = s6c_read_cohort(stage6a_root, "s6a_full_data", 36L),
    s6b_holdout_training = s6c_read_cohort(
      stage6b_root, "s6b_holdout_training", 20L
    )
  )
  descriptors <- unlist(lapply(cohorts, `[[`, "descriptors"),
                        recursive = FALSE)
  endpoints <- do.call(rbind, lapply(descriptors, `[[`, "endpoint"))
  rownames(endpoints) <- NULL
  inputs <- do.call(rbind, lapply(cohorts, `[[`, "inputs"))
  rownames(inputs) <- NULL
  s6c_assert(nrow(endpoints) == 56L && nrow(inputs) == 112L,
             "Unexpected Stage 6C endpoint or input count.")

  pair_rows <- list()
  groups <- split(descriptors, vapply(descriptors, function(value) {
    value$endpoint$config_group_id[[1L]]
  }, character(1L)))
  for (group in groups) {
    time_lengths <- vapply(group, function(value) length(value$time_g),
                           integer(1L))
    s6c_assert(length(unique(time_lengths)) == 1L,
               "Dense-grid lengths differ within an audit stratum.")
    if (length(group) >= 2L) {
      indices <- utils::combn(seq_along(group), 2L)
      for (column in seq_len(ncol(indices))) {
        pair_rows[[length(pair_rows) + 1L]] <- s6c_pair_row(
          group[[indices[1L, column]]], group[[indices[2L, column]]]
        )
      }
    }
  }
  pairs <- do.call(rbind, pair_rows)
  rownames(pairs) <- NULL

  assignments <- list()
  endpoint_groups <- split(endpoints, endpoints$config_group_id)
  for (group_id in names(endpoint_groups)) {
    group_endpoints <- endpoint_groups[[group_id]]
    group_pairs <- pairs[pairs$config_group_id == group_id, , drop = FALSE]
    for (mode in c("contribution", "structural")) {
      for (threshold in thresholds) {
        assignments[[length(assignments) + 1L]] <- s6c_cluster_group(
          group_endpoints, group_pairs, mode, threshold
        )
      }
    }
  }
  assignments <- do.call(rbind, assignments)
  rownames(assignments) <- NULL

  assignment_groups <- split(
    assignments,
    interaction(assignments$config_group_id, assignments$basin_view,
                assignments$similarity_threshold, drop = TRUE)
  )
  basin_summary <- do.call(rbind, lapply(assignment_groups, function(rows) {
    representatives <- rows[rows$is_basin_representative, , drop = FALSE]
    modal_representatives <- representatives[
      representatives$is_modal_basin, , drop = FALSE
    ]
    modal_choice <- modal_representatives$fit_id[[
      which.max(modal_representatives$final_elbo)
    ]]
    data.frame(
      cohort_id = rows$cohort_id[[1L]],
      config_group_id = rows$config_group_id[[1L]],
      fit_config_id = rows$fit_config_id[[1L]],
      starts = nrow(rows),
      basin_view = rows$basin_view[[1L]],
      similarity_threshold = rows$similarity_threshold[[1L]],
      basin_count = length(unique(rows$basin_id)),
      singleton_basin_count = length(unique(
        rows$basin_id[rows$basin_size == 1L]
      )),
      supported_basin_count_ge_2 = length(unique(
        rows$basin_id[rows$basin_size >= 2L]
      )),
      modal_basin_size = max(rows$basin_size),
      modal_basin_support = max(rows$basin_support),
      raw_elbo_winner_fit_id = rows$raw_elbo_winner_fit_id[[1L]],
      raw_elbo_winner_basin_size = rows$basin_size[
        rows$is_raw_elbo_winner][[1L]],
      raw_elbo_winner_basin_support = rows$basin_support[
        rows$is_raw_elbo_winner][[1L]],
      raw_elbo_winner_in_modal_basin = rows$raw_elbo_winner_in_modal_basin[[1L]],
      modal_basin_then_elbo_fit_id = modal_choice,
      selector_changes_endpoint = modal_choice !=
        rows$raw_elbo_winner_fit_id[[1L]],
      truth_used = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(basin_summary) <- NULL
  basin_summary <- basin_summary[order(
    basin_summary$cohort_id, basin_summary$fit_config_id,
    basin_summary$basin_view, basin_summary$similarity_threshold
  ), , drop = FALSE]

  integrity <- data.frame(
    check_id = c(
      "endpoint_count_56", "input_count_112", "stage6a_count_36",
      "stage6b_count_20", "eight_audit_strata", "all_objective_eligible",
      "all_truth_flags_false", "all_continuation_flags_false",
      "pair_similarities_finite_bounded", "truth_or_evaluation_not_read"
    ),
    passed = c(
      nrow(endpoints) == 56L,
      nrow(inputs) == 112L,
      sum(endpoints$cohort_id == "s6a_full_data") == 36L,
      sum(endpoints$cohort_id == "s6b_holdout_training") == 20L,
      length(unique(endpoints$config_group_id)) == 8L,
      TRUE,
      all(!endpoints$truth_used),
      TRUE,
      all(is.finite(pairs$contribution_similarity_min) &
            pairs$contribution_similarity_min >= 0 &
            pairs$contribution_similarity_min <= 1 &
            is.finite(pairs$structural_similarity_min) &
            pairs$structural_similarity_min >= 0 &
            pairs$structural_similarity_min <= 1),
      all(inputs$input_type %in% c("fit", "terminal")) &&
        !any(grepl("/evaluation/|sealed_truth|truth_bundle",
                   inputs$path, ignore.case = TRUE))
    ),
    detail = c(
      nrow(endpoints), nrow(inputs),
      sum(endpoints$cohort_id == "s6a_full_data"),
      sum(endpoints$cohort_id == "s6b_holdout_training"),
      length(unique(endpoints$config_group_id)),
      "enforced while reading each terminal record",
      "all endpoint truth flags FALSE",
      "enforced while reading each terminal record",
      paste0("contribution range=", paste(
        range(pairs$contribution_similarity_min), collapse = ";"
      ), "; structural range=", paste(
        range(pairs$structural_similarity_min), collapse = ";"
      )),
      "allowlist contains fit.rds and terminal_record.rds only"
    ),
    stringsAsFactors = FALSE
  )
  s6c_assert(all(integrity$passed), paste0(
    "Stage 6C truth-free integrity failure: ",
    paste(integrity$check_id[!integrity$passed], collapse = ";")
  ))

  output_files <- c(
    "INPUT_INVENTORY.csv", "FIT_ENDPOINTS.csv",
    "PAIRWISE_BASIN_SIMILARITY.csv", "BASIN_ASSIGNMENTS_TRUTH_FREE.csv",
    "BASIN_SUMMARY_TRUTH_FREE.csv", "TRUTH_FREE_INTEGRITY.csv"
  )
  s6c_write_csv(inputs, file.path(output_root, output_files[[1L]]))
  s6c_write_csv(endpoints, file.path(output_root, output_files[[2L]]))
  s6c_write_csv(pairs, file.path(output_root, output_files[[3L]]))
  s6c_write_csv(assignments, file.path(output_root, output_files[[4L]]))
  s6c_write_csv(basin_summary, file.path(output_root, output_files[[5L]]))
  s6c_write_csv(integrity, file.path(output_root, output_files[[6L]]))
  hashes <- data.frame(
    relative_path = output_files,
    sha256 = vapply(file.path(output_root, output_files), s6c_sha256,
                    character(1L)),
    stringsAsFactors = FALSE
  )
  s6c_write_csv(hashes, file.path(output_root,
                                  "TRUTH_FREE_OUTPUT_HASHES.csv"))
  s6c_write_lines(c(
    "status=STAGE6C_TRUTH_FREE_READONLY_AUDIT_COMPLETE",
    "fit_endpoints=56",
    "audit_strata=8",
    paste0("pairwise_rows=", nrow(pairs)),
    "truth_or_evaluation_read=FALSE",
    "fit_or_continuation_started=FALSE",
    "primary_similarity_threshold=0.90",
    "thresholds_reported=0.80;0.90;0.95"
  ), file.path(output_root, "TRUTH_FREE_AUDIT_COMPLETE.txt"))
  invisible(list(endpoints = endpoints, pairs = pairs,
                 assignments = assignments, summary = basin_summary,
                 integrity = integrity, hashes = hashes))
}

s6c_verify_truth_free_freeze <- function(output_root) {
  marker <- file.path(output_root, "TRUTH_FREE_AUDIT_COMPLETE.txt")
  hash_path <- file.path(output_root, "TRUTH_FREE_OUTPUT_HASHES.csv")
  s6c_assert(file.exists(marker) && file.exists(hash_path),
             "Truth-free Stage 6C output is not frozen.")
  hashes <- utils::read.csv(hash_path, stringsAsFactors = FALSE,
                            check.names = FALSE)
  paths <- file.path(output_root, hashes$relative_path)
  s6c_assert(all(file.exists(paths)), "A frozen truth-free file is missing.")
  observed <- vapply(paths, s6c_sha256, character(1L))
  s6c_assert(identical(unname(tolower(observed)),
                       tolower(as.character(hashes$sha256))),
             "A frozen truth-free hash does not match.")
  invisible(TRUE)
}
