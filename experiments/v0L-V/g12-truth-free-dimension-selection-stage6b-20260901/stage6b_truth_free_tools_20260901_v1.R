s6b_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

s6b_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(value, path, row.names = FALSE, na = "")
  invisible(path)
}

s6b_write_lines <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(as.character(value)), path, useBytes = TRUE)
  invisible(path)
}

s6b_abs_cosine <- function(left, right) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  s6b_assert(length(left) == length(right), "Cosine vectors differ in length.")
  denominator <- sqrt(sum(left^2) * sum(right^2))
  if (!is.finite(denominator) || denominator <= 1e-14) return(0)
  value <- abs(sum(left * right) / denominator)
  min(1, max(0, value))
}

s6b_rms <- function(value) {
  value <- as.numeric(value)
  if (!length(value)) return(0)
  sqrt(mean(value^2))
}

s6b_permutations <- function(n) {
  n <- as.integer(n)
  if (n == 1L) return(matrix(1L, nrow = 1L))
  previous <- s6b_permutations(n - 1L)
  rows <- vector("list", n * nrow(previous))
  index <- 0L
  for (i in seq_len(nrow(previous))) {
    for (position in seq_len(n)) {
      index <- index + 1L
      old <- previous[i, ]
      rows[[index]] <- append(old, n, after = position - 1L)
    }
  }
  do.call(rbind, rows)
}

s6b_best_assignment <- function(similarity) {
  similarity <- as.matrix(similarity)
  s6b_assert(nrow(similarity) == ncol(similarity) && nrow(similarity) >= 1L,
             "Matching requires a non-empty square similarity matrix.")
  permutations <- s6b_permutations(nrow(similarity))
  scores <- apply(permutations, 1L, function(permutation) {
    sum(similarity[cbind(seq_len(nrow(similarity)), permutation)])
  })
  permutations[which.max(scores), ]
}

s6b_trapezoid_weights <- function(time) {
  time <- as.numeric(time)
  s6b_assert(length(time) >= 2L && all(is.finite(time)) &&
               all(diff(time) > 0), "Invalid dense time grid.")
  difference <- diff(time)
  c(difference[[1L]] / 2,
    (difference[-length(difference)] + difference[-1L]) / 2,
    difference[[length(difference)]] / 2)
}

s6b_weighted_basis <- function(functions, time, tolerance = 1e-10) {
  functions <- as.matrix(functions)
  s6b_assert(nrow(functions) == length(time),
             "Function grid and time grid differ in length.")
  weighted <- sqrt(s6b_trapezoid_weights(time)) * functions
  decomposition <- svd(weighted, nu = min(dim(weighted)), nv = 0L)
  if (!length(decomposition$d) || decomposition$d[[1L]] <= 0) {
    return(matrix(numeric(), nrow = nrow(functions), ncol = 0L))
  }
  rank <- sum(decomposition$d > tolerance * decomposition$d[[1L]])
  decomposition$u[, seq_len(rank), drop = FALSE]
}

s6b_subspace_similarity <- function(left, right, time) {
  left_basis <- s6b_weighted_basis(left, time)
  right_basis <- s6b_weighted_basis(right, time)
  rank_left <- ncol(left_basis)
  rank_right <- ncol(right_basis)
  if (!rank_left || !rank_right) {
    return(c(mean_squared_canonical_correlation = 0,
             minimum_canonical_correlation = 0,
             rank_left = rank_left, rank_right = rank_right))
  }
  singular_values <- svd(crossprod(left_basis, right_basis),
                         nu = 0L, nv = 0L)$d
  k <- max(rank_left, rank_right)
  padded <- c(singular_values, rep(0, k - length(singular_values)))
  c(mean_squared_canonical_correlation = mean(padded^2),
    minimum_canonical_correlation = min(padded),
    rank_left = rank_left, rank_right = rank_right)
}

s6b_euclidean_basis <- function(value, tolerance = 1e-10) {
  value <- as.matrix(value)
  decomposition <- svd(value, nu = min(dim(value)), nv = 0L)
  if (!length(decomposition$d) || decomposition$d[[1L]] <= 0) {
    return(matrix(numeric(), nrow = nrow(value), ncol = 0L))
  }
  rank <- sum(decomposition$d > tolerance * decomposition$d[[1L]])
  decomposition$u[, seq_len(rank), drop = FALSE]
}

s6b_euclidean_subspace_similarity <- function(left, right) {
  left_basis <- s6b_euclidean_basis(left)
  right_basis <- s6b_euclidean_basis(right)
  rank_left <- ncol(left_basis)
  rank_right <- ncol(right_basis)
  if (!rank_left || !rank_right) {
    return(c(mean_squared_canonical_correlation = 0,
             minimum_canonical_correlation = 0,
             rank_left = rank_left, rank_right = rank_right))
  }
  singular_values <- svd(crossprod(left_basis, right_basis),
                         nu = 0L, nv = 0L)$d
  k <- max(rank_left, rank_right)
  padded <- c(singular_values, rep(0, k - length(singular_values)))
  c(mean_squared_canonical_correlation = mean(padded^2),
    minimum_canonical_correlation = min(padded),
    rank_left = rank_left, rank_right = rank_right)
}

s6b_fit_has_forbidden_truth_fields <- function(fit) {
  any(grepl("(^|_)(truth|true)(_|$)", names(fit), ignore.case = TRUE))
}

s6b_allowed_fit_file <- function(stage6a_root, fit_id, filename) {
  s6b_assert(filename %in% c("fit.rds", "terminal_record.rds"),
             "Only fit and terminal RDS files are allowed.")
  root <- normalizePath(stage6a_root, winslash = "/", mustWork = TRUE)
  path <- file.path(root, "fits", fit_id, filename)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  allowed_prefix <- paste0(root, "/fits/", fit_id, "/")
  s6b_assert(startsWith(resolved, allowed_prefix),
             "Resolved fit input escaped the Stage 6A fits directory.")
  resolved
}

s6b_factor_paths <- function(scores, functions) {
  scores <- as.matrix(scores)
  functions <- as.matrix(functions)
  s6b_assert(ncol(scores) == ncol(functions),
             "Score and eigenfunction component counts differ.")
  scores %*% t(functions)
}

s6b_extract_fit <- function(fit, record) {
  required <- c(
    "S", "n_s", "p", "L_f", "L_s_by_study", "time_g",
    "mu_q_a_hat", "mu_q_gamma_a_hat", "factor_ppi_shared",
    "list_Zeta_hat", "list_Phi_hat", "list_eigenvalues",
    "list_cumulated_pve", "list_effective_M", "list_rank_cap",
    "mu_q_b_specific_hat", "mu_q_gamma_b_hat", "factor_ppi_specific",
    "list_Zeta_hat_spec", "list_Phi_hat_spec", "list_eigenvalues_spec",
    "list_pve_spec", "list_effective_M_spec", "list_rank_cap_spec"
  )
  s6b_assert(all(required %in% names(fit)),
             "Fit is missing fields required for the truth-free audit.")
  s6b_assert(!s6b_fit_has_forbidden_truth_fields(fit),
             "Fit contains a forbidden truth-named top-level field.")
  s6b_assert(isTRUE(record$objective_eligible) &&
               isFALSE(record$truth_used_for_fit_stopping_or_selection),
             "Terminal record fails objective or truth-isolation eligibility.")

  candidate_rows <- list()
  candidate_objects <- list()
  component_rows <- list()
  index <- 0L

  add_candidate <- function(role_id, block, study, factor, loading,
                            loading_ppi, factor_ppi, scores, functions,
                            eigenvalues, cumulative_pve, effective_m,
                            rank_cap, subject_weight) {
    index <<- index + 1L
    candidate_id <- if (block == "shared") {
      paste0("A_", factor)
    } else {
      paste0("B", study, "_", factor)
    }
    paths <- s6b_factor_paths(scores, functions)
    loading_l2 <- sqrt(sum(as.numeric(loading)^2))
    loading_rms <- s6b_rms(loading)
    process_rms <- s6b_rms(paths)
    contribution_rms <- loading_rms * process_rms
    weighted_energy <- subject_weight * contribution_rms^2
    candidate_rows[[index]] <<- data.frame(
      fit_id = record$fit_id,
      fit_config_id = record$fit_config_id,
      seed_index = record$seed_index,
      final_elbo = record$final_elbo,
      practical_converged = record$practical_converged,
      role_id = role_id,
      block = block,
      study = study,
      candidate_id = candidate_id,
      candidate_factor = factor,
      factor_ppi = as.numeric(factor_ppi),
      ppi_retained_0_5 = as.numeric(factor_ppi) >= 0.5,
      expected_active_loading_count = sum(as.numeric(loading_ppi)),
      loading_l2 = loading_l2,
      zero_loading_direction = loading_l2 <= 1e-14,
      loading_rms = loading_rms,
      factor_process_rms = process_rms,
      complete_contribution_rms = contribution_rms,
      subject_weight = subject_weight,
      weighted_contribution_energy = weighted_energy,
      rank_cap = as.integer(rank_cap),
      effective_M_99pct = as.integer(effective_m),
      stringsAsFactors = FALSE
    )
    candidate_objects[[paste(role_id, candidate_id, sep = "::")]] <<- list(
      role_id = role_id,
      candidate_id = candidate_id,
      loading = as.numeric(loading),
      functions = as.matrix(functions)
    )
    eigenvalues <- as.numeric(eigenvalues)
    cumulative_pve <- as.numeric(cumulative_pve)
    component_rows[[index]] <<- data.frame(
      fit_id = record$fit_id,
      fit_config_id = record$fit_config_id,
      seed_index = record$seed_index,
      role_id = role_id,
      candidate_id = candidate_id,
      component = seq_len(as.integer(rank_cap)),
      normalized_eigenvalue = c(eigenvalues, rep(NA_real_,
        max(0L, as.integer(rank_cap) - length(eigenvalues))))[
          seq_len(as.integer(rank_cap))],
      cumulative_pve = c(cumulative_pve, rep(NA_real_,
        max(0L, as.integer(rank_cap) - length(cumulative_pve))))[
          seq_len(as.integer(rank_cap))],
      retained_after_99pct = seq_len(as.integer(rank_cap)) <=
        as.integer(effective_m),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  total_subjects <- sum(as.integer(fit$n_s))
  for (factor in seq_len(as.integer(fit$L_f))) {
    add_candidate(
      "A", "shared", 0L, factor,
      fit$mu_q_a_hat[, factor], fit$mu_q_gamma_a_hat[, factor],
      fit$factor_ppi_shared[[factor]], fit$list_Zeta_hat[[factor]],
      fit$list_Phi_hat[[factor]], fit$list_eigenvalues[[factor]],
      fit$list_cumulated_pve[[factor]], fit$list_effective_M[[factor]],
      fit$list_rank_cap[[factor]], 1
    )
  }
  for (study in seq_len(as.integer(fit$S))) {
    for (factor in seq_len(as.integer(fit$L_s_by_study[[study]]))) {
      add_candidate(
        paste0("B", study), "specific", study, factor,
        fit$mu_q_b_specific_hat[[study]][, factor],
        fit$mu_q_gamma_b_hat[[study]][, factor],
        fit$factor_ppi_specific[[study]][[factor]],
        fit$list_Zeta_hat_spec[[study]][[factor]],
        fit$list_Phi_hat_spec[[study]][[factor]],
        fit$list_eigenvalues_spec[[study]][[factor]],
        fit$list_pve_spec[[study]][[factor]],
        fit$list_effective_M_spec[[study]][[factor]],
        fit$list_rank_cap_spec[[study]][[factor]],
        as.integer(fit$n_s[[study]]) / total_subjects
      )
    }
  }
  candidates <- do.call(rbind, candidate_rows)
  total_energy <- sum(candidates$weighted_contribution_energy)
  candidates$weighted_contribution_share <- if (total_energy > 0) {
    candidates$weighted_contribution_energy / total_energy
  } else 0
  list(
    record = record,
    time_g = as.numeric(fit$time_g),
    candidates = candidates,
    components = do.call(rbind, component_rows),
    objects = candidate_objects
  )
}

s6b_role_objects <- function(extracted, role_id) {
  objects <- extracted$objects
  objects[vapply(objects, function(value) identical(value$role_id, role_id),
                 logical(1L))]
}

s6b_read_stage6a_fits <- function(stage6a_root) {
  winner_path <- file.path(stage6a_root, "truth_free_selection",
                           "G12_TRUTH_FREE_WINNERS.csv")
  s6b_assert(file.exists(winner_path), "Missing truth-free Stage 6A winner table.")
  winners <- utils::read.csv(winner_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  s6b_assert(nrow(winners) == 3L && length(unique(winners$fit_config_id)) == 3L,
             "Expected exactly three Stage 6A truth-free winners.")
  fit_directories <- list.dirs(file.path(stage6a_root, "fits"),
                               recursive = FALSE, full.names = FALSE)
  fit_directories <- sort(fit_directories[nzchar(fit_directories)])
  s6b_assert(length(fit_directories) == 36L,
             "Expected exactly 36 Stage 6A fit directories.")
  extracted <- setNames(vector("list", length(fit_directories)), fit_directories)
  input_log <- list()
  for (i in seq_along(fit_directories)) {
    fit_id <- fit_directories[[i]]
    fit_path <- s6b_allowed_fit_file(stage6a_root, fit_id, "fit.rds")
    record_path <- s6b_allowed_fit_file(stage6a_root, fit_id,
                                        "terminal_record.rds")
    fit <- readRDS(fit_path)
    record <- readRDS(record_path)
    s6b_assert(identical(record$fit_id, fit_id),
               "Fit directory and terminal fit_id differ.")
    extracted[[fit_id]] <- s6b_extract_fit(fit, record)
    input_log[[length(input_log) + 1L]] <- data.frame(
      fit_id = fit_id, input_type = "fit", relative_path =
        file.path("fits", fit_id, "fit.rds"), stringsAsFactors = FALSE
    )
    input_log[[length(input_log) + 1L]] <- data.frame(
      fit_id = fit_id, input_type = "terminal", relative_path =
        file.path("fits", fit_id, "terminal_record.rds"),
      stringsAsFactors = FALSE
    )
    rm(fit)
    invisible(gc(FALSE))
  }
  input_log[[length(input_log) + 1L]] <- data.frame(
    fit_id = "", input_type = "truth_free_winner_table",
    relative_path = file.path("truth_free_selection",
                              "G12_TRUTH_FREE_WINNERS.csv"),
    stringsAsFactors = FALSE
  )
  list(extracted = extracted, winners = winners,
       input_log = do.call(rbind, input_log))
}

s6b_audit_stage6a <- function(stage6a_root, output_root,
                              thresholds = c(0.8, 0.9, 0.95)) {
  loaded <- s6b_read_stage6a_fits(stage6a_root)
  fits <- loaded$extracted
  winners <- loaded$winners
  thresholds <- as.numeric(thresholds)
  s6b_assert(all(thresholds > 0 & thresholds <= 1),
             "Similarity thresholds must lie in (0,1].")

  candidate_activity <- do.call(rbind, lapply(fits, `[[`, "candidates"))
  fpca_components <- do.call(rbind, lapply(fits, `[[`, "components"))
  match_rows <- list()
  direction_rows <- list()
  duplicate_rows <- list()
  leakage_rows <- list()
  loading_subspace_rows <- list()
  fpca_subspace_rows <- list()

  for (config in sort(unique(candidate_activity$fit_config_id))) {
    config_fits <- fits[vapply(fits, function(value) {
      identical(value$record$fit_config_id, config)
    }, logical(1L))]
    config_fits <- config_fits[order(vapply(config_fits, function(value) {
      as.integer(value$record$seed_index)
    }, integer(1L)))]
    winner_id <- winners$fit_id[winners$fit_config_id == config]
    s6b_assert(length(winner_id) == 1L && winner_id %in% names(config_fits),
               paste0("Missing truth-free winner fit for ", config, "."))
    reference <- config_fits[[winner_id]]
    roles <- sort(unique(reference$candidates$role_id))

    for (role in roles) {
      reference_objects <- s6b_role_objects(reference, role)
      reference_ids <- vapply(reference_objects, `[[`, character(1L),
                              "candidate_id")
      for (fit_id in names(config_fits)) {
        candidate <- config_fits[[fit_id]]
        candidate_objects <- s6b_role_objects(candidate, role)
        candidate_ids <- vapply(candidate_objects, `[[`, character(1L),
                                "candidate_id")
        s6b_assert(length(reference_objects) == length(candidate_objects),
                   "Candidate count changed within a fit configuration and role.")
        similarity <- outer(seq_along(reference_objects),
                            seq_along(candidate_objects), Vectorize(function(i, j) {
          s6b_abs_cosine(reference_objects[[i]]$loading,
                         candidate_objects[[j]]$loading)
        }))
        assignment <- s6b_best_assignment(similarity)
        for (i in seq_along(reference_objects)) {
          j <- assignment[[i]]
          reference_row <- reference$candidates[
            reference$candidates$role_id == role &
              reference$candidates$candidate_id == reference_ids[[i]],
            , drop = FALSE
          ]
          candidate_row <- candidate$candidates[
            candidate$candidates$role_id == role &
              candidate$candidates$candidate_id == candidate_ids[[j]],
            , drop = FALSE
          ]
          cosine <- similarity[i, j]
          match_rows[[length(match_rows) + 1L]] <- data.frame(
            fit_config_id = config,
            reference_fit_id = winner_id,
            reference_candidate_id = reference_ids[[i]],
            role_id = role,
            candidate_fit_id = fit_id,
            candidate_seed_index = candidate$record$seed_index,
            candidate_id = candidate_ids[[j]],
            loading_abs_cosine = cosine,
            reference_contribution_share =
              reference_row$weighted_contribution_share,
            candidate_contribution_share =
              candidate_row$weighted_contribution_share,
            candidate_complete_contribution_rms =
              candidate_row$complete_contribution_rms,
            candidate_factor_ppi = candidate_row$factor_ppi,
            stringsAsFactors = FALSE
          )
          fpca_similarity <- s6b_subspace_similarity(
            reference_objects[[i]]$functions,
            candidate_objects[[j]]$functions,
            reference$time_g
          )
          fpca_subspace_rows[[length(fpca_subspace_rows) + 1L]] <- data.frame(
            fit_config_id = config,
            reference_fit_id = winner_id,
            reference_candidate_id = reference_ids[[i]],
            role_id = role,
            candidate_fit_id = fit_id,
            candidate_seed_index = candidate$record$seed_index,
            candidate_id = candidate_ids[[j]],
            mean_squared_canonical_correlation =
              fpca_similarity[["mean_squared_canonical_correlation"]],
            minimum_canonical_correlation =
              fpca_similarity[["minimum_canonical_correlation"]],
            reference_rank = fpca_similarity[["rank_left"]],
            candidate_rank = fpca_similarity[["rank_right"]],
            stringsAsFactors = FALSE
          )
        }
        loading_similarity <- s6b_euclidean_subspace_similarity(
          do.call(cbind, lapply(reference_objects, `[[`, "loading")),
          do.call(cbind, lapply(candidate_objects, `[[`, "loading"))
        )
        loading_subspace_rows[[length(loading_subspace_rows) + 1L]] <- data.frame(
          fit_config_id = config,
          reference_fit_id = winner_id,
          role_id = role,
          candidate_fit_id = fit_id,
          candidate_seed_index = candidate$record$seed_index,
          mean_squared_canonical_correlation =
            loading_similarity[["mean_squared_canonical_correlation"]],
          minimum_canonical_correlation =
            loading_similarity[["minimum_canonical_correlation"]],
          reference_rank = loading_similarity[["rank_left"]],
          candidate_rank = loading_similarity[["rank_right"]],
          stringsAsFactors = FALSE
        )

        if (length(candidate_objects) >= 2L) {
          pairs <- utils::combn(seq_along(candidate_objects), 2L)
          pair_cosines <- apply(pairs, 2L, function(pair) {
            s6b_abs_cosine(candidate_objects[[pair[[1L]]]]$loading,
                           candidate_objects[[pair[[2L]]]]$loading)
          })
        } else {
          pair_cosines <- numeric()
        }
        duplicate_row <- data.frame(
          fit_config_id = config,
          fit_id = fit_id,
          seed_index = candidate$record$seed_index,
          role_id = role,
          candidate_count = length(candidate_objects),
          pair_count = length(pair_cosines),
          maximum_within_role_abs_cosine = if (length(pair_cosines)) {
            max(pair_cosines)
          } else 0,
          stringsAsFactors = FALSE
        )
        for (threshold in thresholds) {
          duplicate_row[[paste0("pairs_ge_", gsub("[.]", "_", threshold))]] <-
            sum(pair_cosines >= threshold)
        }
        duplicate_rows[[length(duplicate_rows) + 1L]] <- duplicate_row
      }
    }

    for (fit_id in names(config_fits)) {
      candidate <- config_fits[[fit_id]]
      shared <- s6b_role_objects(candidate, "A")
      for (study in seq_len(2L)) {
        specific <- s6b_role_objects(candidate, paste0("B", study))
        cosines <- outer(seq_along(shared), seq_along(specific),
                         Vectorize(function(i, j) {
          s6b_abs_cosine(shared[[i]]$loading, specific[[j]]$loading)
        }))
        leakage_row <- data.frame(
          fit_config_id = config,
          fit_id = fit_id,
          seed_index = candidate$record$seed_index,
          study = study,
          pair_count = length(cosines),
          maximum_cross_role_abs_cosine = max(cosines),
          stringsAsFactors = FALSE
        )
        for (threshold in thresholds) {
          leakage_row[[paste0("pairs_ge_", gsub("[.]", "_", threshold))]] <-
            sum(cosines >= threshold)
        }
        leakage_rows[[length(leakage_rows) + 1L]] <- leakage_row
      }
    }
  }

  matches <- do.call(rbind, match_rows)
  duplicates <- do.call(rbind, duplicate_rows)
  leakage <- do.call(rbind, leakage_rows)
  loading_subspace <- do.call(rbind, loading_subspace_rows)
  fpca_subspace <- do.call(rbind, fpca_subspace_rows)

  match_groups <- split(matches, interaction(matches$fit_config_id,
                                              matches$reference_candidate_id,
                                              drop = TRUE))
  direction_summary <- do.call(rbind, lapply(match_groups, function(group) {
    row <- data.frame(
      fit_config_id = group$fit_config_id[[1L]],
      reference_fit_id = group$reference_fit_id[[1L]],
      role_id = group$role_id[[1L]],
      reference_candidate_id = group$reference_candidate_id[[1L]],
      starts = nrow(group),
      loading_abs_cosine_median = stats::median(group$loading_abs_cosine),
      loading_abs_cosine_minimum = min(group$loading_abs_cosine),
      contribution_share_median =
        stats::median(group$candidate_contribution_share),
      contribution_share_minimum = min(group$candidate_contribution_share),
      factor_ppi_median = stats::median(group$candidate_factor_ppi),
      stringsAsFactors = FALSE
    )
    for (threshold in thresholds) {
      row[[paste0("recurrence_ge_", gsub("[.]", "_", threshold))]] <-
        mean(group$loading_abs_cosine >= threshold)
    }
    row
  }))
  rownames(direction_summary) <- NULL

  fit_groups <- split(candidate_activity, candidate_activity$fit_id)
  fit_summary <- do.call(rbind, lapply(fit_groups, function(group) {
    data.frame(
      fit_id = group$fit_id[[1L]],
      fit_config_id = group$fit_config_id[[1L]],
      seed_index = group$seed_index[[1L]],
      final_elbo = group$final_elbo[[1L]],
      practical_converged = group$practical_converged[[1L]],
      candidate_count = nrow(group),
      ppi_retained_count = sum(group$ppi_retained_0_5),
      contribution_share_max = max(group$weighted_contribution_share),
      contribution_share_min = min(group$weighted_contribution_share),
      effective_M_mean = mean(group$effective_M_99pct),
      effective_M_at_cap = sum(group$effective_M_99pct == group$rank_cap),
      stringsAsFactors = FALSE
    )
  }))
  rownames(fit_summary) <- NULL

  config_summary <- do.call(rbind, lapply(sort(unique(fit_summary$fit_config_id)),
                                         function(config) {
    fit_group <- fit_summary[fit_summary$fit_config_id == config, , drop = FALSE]
    direction_group <- direction_summary[
      direction_summary$fit_config_id == config, , drop = FALSE
    ]
    duplicate_group <- duplicates[duplicates$fit_config_id == config,
                                  , drop = FALSE]
    leakage_group <- leakage[leakage$fit_config_id == config, , drop = FALSE]
    data.frame(
      fit_config_id = config,
      fits = nrow(fit_group),
      candidates_per_fit = unique(fit_group$candidate_count)[[1L]],
      ppi_retained_per_fit_median = stats::median(fit_group$ppi_retained_count),
      minimum_reference_direction_recurrence_ge_0_9 =
        min(direction_group$recurrence_ge_0_9),
      median_reference_direction_recurrence_ge_0_9 =
        stats::median(direction_group$recurrence_ge_0_9),
      fit_role_duplicate_rate_ge_0_9 =
        mean(duplicate_group$pairs_ge_0_9 > 0),
      fit_study_cross_role_leakage_rate_ge_0_9 =
        mean(leakage_group$pairs_ge_0_9 > 0),
      median_maximum_cross_role_abs_cosine =
        stats::median(leakage_group$maximum_cross_role_abs_cosine),
      median_effective_M = stats::median(
        candidate_activity$effective_M_99pct[
          candidate_activity$fit_config_id == config]
      ),
      stringsAsFactors = FALSE
    )
  }))

  integrity <- data.frame(
    check_id = c(
      "fit_count_36", "winner_count_3", "candidate_rows_complete",
      "input_allowlist_only", "truth_use_flags_false",
      "match_similarities_bounded", "contribution_shares_sum_one"
    ),
    passed = c(
      length(fits) == 36L,
      nrow(winners) == 3L,
      nrow(candidate_activity) == 180L,
      all(loaded$input_log$input_type %in%
            c("fit", "terminal", "truth_free_winner_table")),
      all(vapply(fits, function(value) {
        isFALSE(value$record$truth_used_for_fit_stopping_or_selection)
      }, logical(1L))),
      all(is.finite(matches$loading_abs_cosine) &
            matches$loading_abs_cosine >= 0 &
            matches$loading_abs_cosine <= 1),
      all(vapply(split(candidate_activity$weighted_contribution_share,
                       candidate_activity$fit_id), function(value) {
        abs(sum(value) - 1) <= 1e-10
      }, logical(1L)))
    ),
    detail = c(
      as.character(length(fits)), as.character(nrow(winners)),
      as.character(nrow(candidate_activity)),
      paste(sort(unique(loaded$input_log$input_type)), collapse = ";"),
      "all terminal truth-use flags FALSE",
      paste0("range=", paste(range(matches$loading_abs_cosine), collapse = ";"),
             ";zero_loading_directions=",
             sum(candidate_activity$zero_loading_direction)),
      "per-fit weighted shares"
    ),
    stringsAsFactors = FALSE
  )
  s6b_assert(all(integrity$passed), paste0(
    "Stage 6B-A integrity failure: ",
    paste(integrity$check_id[!integrity$passed], collapse = ";")
  ))

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  s6b_write_csv(loaded$input_log, file.path(output_root, "READ_INPUT_LOG.csv"))
  s6b_write_csv(candidate_activity,
                file.path(output_root, "FIT_CANDIDATE_ACTIVITY.csv"))
  s6b_write_csv(fpca_components,
                file.path(output_root, "FPCA_COMPONENT_ACTIVITY.csv"))
  s6b_write_csv(matches,
                file.path(output_root, "REFERENCE_MATCH_STABILITY.csv"))
  s6b_write_csv(direction_summary,
                file.path(output_root, "REFERENCE_DIRECTION_SUMMARY.csv"))
  s6b_write_csv(duplicates,
                file.path(output_root, "WITHIN_ROLE_DUPLICATION.csv"))
  s6b_write_csv(leakage,
                file.path(output_root, "CROSS_ROLE_LEAKAGE.csv"))
  s6b_write_csv(loading_subspace,
                file.path(output_root, "LOADING_SUBSPACE_STABILITY.csv"))
  s6b_write_csv(fpca_subspace,
                file.path(output_root, "FPCA_SUBSPACE_STABILITY.csv"))
  s6b_write_csv(fit_summary,
                file.path(output_root, "FIT_ACTIVITY_SUMMARY.csv"))
  s6b_write_csv(config_summary,
                file.path(output_root, "CONFIG_AUDIT_SUMMARY.csv"))
  s6b_write_csv(integrity,
                file.path(output_root, "AUDIT_INTEGRITY.csv"))
  s6b_write_lines(c(
    "status=STAGE6B_A_READONLY_AUDIT_COMPLETE",
    paste0("input_fit_count=", length(fits)),
    paste0("candidate_rows=", nrow(candidate_activity)),
    "truth_fields_read=FALSE",
    "fit_or_continuation_started=FALSE"
  ), file.path(output_root, "AUDIT_COMPLETE.txt"))
  invisible(list(
    candidate_activity = candidate_activity,
    fpca_components = fpca_components,
    matches = matches,
    direction_summary = direction_summary,
    duplicates = duplicates,
    leakage = leakage,
    loading_subspace = loading_subspace,
    fpca_subspace = fpca_subspace,
    fit_summary = fit_summary,
    config_summary = config_summary,
    integrity = integrity
  ))
}

s6b_forbidden_observation_names <- function(value) {
  found <- character()
  walk <- function(object) {
    object_names <- names(object)
    if (!is.null(object_names)) {
      found <<- c(found, object_names[
        grepl("(^|_)(truth|true)(_|$)|signal_true|noise_true",
              object_names, ignore.case = TRUE)
      ])
    }
    if (is.list(object)) for (item in object) walk(item)
    invisible(NULL)
  }
  walk(value)
  unique(found)
}

s6b_make_middle_time_holdout <- function(observation, min_train_times = 5L) {
  required <- c("Y", "time_obs", "Z", "dimensions")
  s6b_assert(all(required %in% names(observation)),
             "Observation bundle is incomplete.")
  s6b_assert(!length(s6b_forbidden_observation_names(observation)),
             "Observation bundle contains a forbidden truth field.")
  s6b_assert(is.null(observation$Z) &&
               as.integer(observation$dimensions$d) == 0L,
             "Stage 6B-B v1 supports the registered d=0 route only.")
  train <- observation
  records <- list()
  index <- 0L
  for (study in seq_along(observation$Y)) {
    s6b_assert(length(observation$Y[[study]]) ==
                 length(observation$time_obs[[study]]),
               "Subject counts differ between Y and time_obs.")
    for (subject in seq_along(observation$Y[[study]])) {
      time <- as.numeric(observation$time_obs[[study]][[subject]])
      n_time <- length(time)
      s6b_assert(n_time - 1L >= as.integer(min_train_times),
                 "Holdout would leave too few training times.")
      holdout_index <- as.integer(floor((n_time + 1L) / 2L))
      s6b_assert(holdout_index > 1L && holdout_index < n_time,
                 "Middle holdout is not an interior time.")
      responses <- observation$Y[[study]][[subject]]
      p <- length(responses)
      s6b_assert(p >= 1L && all(vapply(responses, length, integer(1L)) == n_time),
                 "Variables do not share the registered subject time grid.")
      heldout_y <- vapply(responses, function(value) {
        as.numeric(value[[holdout_index]])
      }, numeric(1L))
      train$Y[[study]][[subject]] <- lapply(responses, function(value) {
        as.numeric(value[-holdout_index])
      })
      train$time_obs[[study]][[subject]] <- time[-holdout_index]
      index <- index + 1L
      records[[index]] <- list(
        study = as.integer(study),
        subject = as.integer(subject),
        original_index = holdout_index,
        time = time[[holdout_index]],
        y = heldout_y,
        p = p,
        train_time_count = n_time - 1L
      )
    }
  }
  plan <- list(
    split_id = "middle_interior_whole_time_v1",
    records = records,
    subject_time_count = length(records),
    complete_variable_vector_held_out = TRUE,
    truth_used = FALSE
  )
  list(train_observation = train, holdout_plan = plan)
}

s6b_interpolate <- function(time_grid, values, time) {
  stats::approx(time_grid, as.numeric(values), xout = time,
                method = "linear", rule = 2, ties = "ordered")$y
}

s6b_predict_subject_time <- function(fit, study, subject, time,
                                     factor_mode = c("all", "ppi_0.5")) {
  factor_mode <- match.arg(factor_mode)
  s6b_assert(as.integer(fit$d) == 0L,
             "Stage 6B-B v1 prediction supports d=0 fits only.")
  p <- as.integer(fit$p)
  prediction <- vapply(seq_len(p), function(variable) {
    s6b_interpolate(fit$time_g, fit$list_mu_hat[[study]][[variable]], time)
  }, numeric(1L))
  include <- function(ppi) factor_mode == "all" || as.numeric(ppi) >= 0.5
  global_subject <- sum(as.integer(fit$n_s)[seq_len(study - 1L)]) + subject
  for (factor in seq_len(as.integer(fit$L_f))) {
    if (include(fit$factor_ppi_shared[[factor]])) {
      phi <- vapply(seq_len(ncol(fit$list_Phi_hat[[factor]])), function(m) {
        s6b_interpolate(fit$time_g, fit$list_Phi_hat[[factor]][, m], time)
      }, numeric(1L))
      process <- sum(fit$list_Zeta_hat[[factor]][global_subject, ] * phi)
      prediction <- prediction + fit$mu_q_a_hat[, factor] * process
    }
  }
  for (factor in seq_len(as.integer(fit$L_s_by_study[[study]]))) {
    if (include(fit$factor_ppi_specific[[study]][[factor]])) {
      functions <- fit$list_Phi_hat_spec[[study]][[factor]]
      phi <- vapply(seq_len(ncol(functions)), function(m) {
        s6b_interpolate(fit$time_g, functions[, m], time)
      }, numeric(1L))
      process <- sum(fit$list_Zeta_hat_spec[[study]][[factor]][subject, ] * phi)
      prediction <- prediction +
        fit$mu_q_b_specific_hat[[study]][, factor] * process
    }
  }
  as.numeric(prediction)
}

s6b_score_holdout <- function(fit, holdout_plan,
                              factor_mode = c("all", "ppi_0.5"),
                              variance_floor = 1e-8) {
  factor_mode <- match.arg(factor_mode)
  s6b_assert(is.list(holdout_plan) &&
               identical(holdout_plan$truth_used, FALSE),
             "Invalid or truth-contaminated holdout plan.")
  s6b_assert(is.matrix(fit$sigsq_eps_hat) &&
               nrow(fit$sigsq_eps_hat) == fit$S &&
               ncol(fit$sigsq_eps_hat) == fit$p,
             "Fit lacks residual variance estimates.")
  rows <- lapply(holdout_plan$records, function(record) {
    prediction <- s6b_predict_subject_time(
      fit, record$study, record$subject, record$time, factor_mode
    )
    observed <- as.numeric(record$y)
    s6b_assert(length(prediction) == length(observed) &&
                 all(is.finite(prediction)) && all(is.finite(observed)),
               "Held-out observation or prediction is invalid.")
    residual <- observed - prediction
    variance <- pmax(as.numeric(fit$sigsq_eps_hat[record$study, ]),
                     variance_floor)
    residual_only_nlpd <- mean(
      0.5 * log(2 * pi * variance) + 0.5 * residual^2 / variance
    )
    data.frame(
      study = record$study,
      subject = record$subject,
      holdout_time = record$time,
      original_index = record$original_index,
      train_time_count = record$train_time_count,
      variables = length(observed),
      mse = mean(residual^2),
      mae = mean(abs(residual)),
      observed_second_moment = mean(observed^2),
      residual_only_nlpd = residual_only_nlpd,
      stringsAsFactors = FALSE
    )
  })
  subject_time <- do.call(rbind, rows)
  mse <- mean(subject_time$mse)
  denominator <- sqrt(mean(subject_time$observed_second_moment))
  aggregate <- data.frame(
    factor_mode = factor_mode,
    subject_time_count = nrow(subject_time),
    variables_per_time_min = min(subject_time$variables),
    variables_per_time_max = max(subject_time$variables),
    heldout_rmse = sqrt(mse),
    heldout_nrmse = if (denominator > 0) sqrt(mse) / denominator else NA_real_,
    heldout_mae = mean(subject_time$mae),
    residual_only_nlpd = mean(subject_time$residual_only_nlpd),
    predictive_variance_scope = "residual_noise_only_not_calibrated_full_posterior",
    primary_dimension_score_ready = FALSE,
    truth_used = FALSE,
    stringsAsFactors = FALSE
  )
  list(subject_time = subject_time, aggregate = aggregate)
}
