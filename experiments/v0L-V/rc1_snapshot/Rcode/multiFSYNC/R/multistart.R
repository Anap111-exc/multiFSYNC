# Deterministic endpoint selection: among numerically valid ordinary-objective
# endpoints, choose the largest ELBO, then the earliest explicit start order
# for exact ties. Practical convergence is reported separately; a valid slow
# endpoint may be returned but can never be labelled production-ready.
.select_multistart_summary <- function(summary_table) {
  candidate_valid <- if ("candidate_valid" %in% names(summary_table)) {
    summary_table$candidate_valid
  } else {
    summary_table$eligible
  }
  candidate_ids <- which(candidate_valid)
  if (!length(candidate_ids)) return(NA_integer_)
  candidate_elbo <- summary_table$final_elbo[candidate_ids]
  candidate_ids[which.max(candidate_elbo)]
}

.count_independent_multistart_support <- function(summary_table, member) {
  if (!length(member) || !any(member)) return(0L)
  if (!"initialization_independence_id" %in% names(summary_table)) {
    return(as.integer(sum(member)))
  }
  ids <- as.character(summary_table$initialization_independence_id)
  valid <- member & !is.na(ids) & nzchar(ids)
  as.integer(length(unique(ids[valid])))
}

.evaluate_multistart_state <- function(
    summary_table, stability_snapshots,
    output_control, basin_rel_tol, basin_abs_tol,
    min_exploration_starts = 1L) {
  if (!"initialization_independence_id" %in% names(summary_table)) {
    summary_table$initialization_independence_id <-
      paste0("legacy_start_", seq_len(nrow(summary_table)))
  }
  if (!"candidate_valid" %in% names(summary_table)) {
    summary_table$candidate_valid <- summary_table$eligible
  }
  selected_start <- .select_multistart_summary(summary_table)
  summary_table$selected <- FALSE
  summary_table$basin_absolute_gap <- NA_real_
  summary_table$basin_relative_gap <- NA_real_
  summary_table$objective_basin_member <- FALSE
  summary_table$fitted_nrmse_to_selected <- NA_real_
  summary_table$rss_rel_to_selected <- NA_real_
  summary_table$ppi_max_abs_to_selected <- NA_real_
  summary_table$output_aligned <- FALSE
  summary_table$higher_than_selected_unfinished <- FALSE
  summary_table$higher_than_selected_unreliable <- FALSE

  best_basin_support <- 0L
  best_basin_output_support <- 0L
  higher_ineligible_competitor <- FALSE
  higher_ineligible_start_ids <- integer()
  higher_ineligible_seeds <- integer()
  higher_unreliable_competitor <- FALSE
  higher_unreliable_start_ids <- integer()
  higher_unreliable_seeds <- integer()
  within_basin_output_disagreement <- FALSE
  selected_snapshot_missing <- FALSE
  selected_endpoint_unfinished <- FALSE

  if (!is.na(selected_start)) {
    summary_table$selected[selected_start] <- TRUE
    selected_endpoint_unfinished <-
      !isTRUE(summary_table$eligible[selected_start])
    selected_elbo <- summary_table$final_elbo[selected_start]
    selected_scale <- max(1, abs(selected_elbo))
    basin_tolerance <- max(
      basin_abs_tol, basin_rel_tol * selected_scale)
    summary_table$basin_absolute_gap <-
      selected_elbo - summary_table$final_elbo
    summary_table$basin_relative_gap <-
      summary_table$basin_absolute_gap / selected_scale
    summary_table$objective_basin_member <-
      summary_table$candidate_valid &
      is.finite(summary_table$basin_absolute_gap) &
      abs(summary_table$basin_absolute_gap) <= basin_tolerance
    best_basin_support <- .count_independent_multistart_support(
      summary_table, summary_table$objective_basin_member)

    # Kept as backward-compatible columns. Valid unfinished endpoints now
    # participate directly in selection, so none can be higher than the
    # selected valid endpoint.
    summary_table$higher_than_selected_unfinished <- FALSE
    unreliable_ineligible <- !summary_table$candidate_valid &
      summary_table$success &
      is.finite(summary_table$final_elbo) &
      TRUE
    summary_table$higher_than_selected_unreliable <-
      unreliable_ineligible &
      summary_table$final_elbo > selected_elbo + basin_tolerance
    higher_unreliable_start_ids <- which(
      summary_table$higher_than_selected_unreliable)
    higher_unreliable_seeds <-
      summary_table$fit_seed[higher_unreliable_start_ids]
    higher_unreliable_competitor <-
      length(higher_unreliable_start_ids) > 0L

    selected_snapshot <- stability_snapshots[[selected_start]]
    selected_snapshot_missing <- is.null(selected_snapshot)
    if (!selected_snapshot_missing) {
      for (start_id in seq_len(nrow(summary_table))) {
        snapshot_now <- stability_snapshots[[start_id]]
        if (is.null(snapshot_now)) next
        comparison <- .practical_compare(
          selected_snapshot, snapshot_now, output_control)
        summary_table$fitted_nrmse_to_selected[start_id] <-
          comparison$fitted_nrmse
        summary_table$rss_rel_to_selected[start_id] <-
          comparison$rss_rel
        summary_table$ppi_max_abs_to_selected[start_id] <-
          comparison$ppi_max_abs
        summary_table$output_aligned[start_id] <- comparison$pass
      }
    }
    aligned_basin <- summary_table$objective_basin_member &
      summary_table$output_aligned
    best_basin_output_support <- .count_independent_multistart_support(
      summary_table, aligned_basin)
    within_basin_output_disagreement <- any(
      summary_table$objective_basin_member &
        !summary_table$output_aligned
    )
  }

  objective_competition_unresolved <- !is.na(selected_start) &&
    (selected_endpoint_unfinished ||
       best_basin_support < 2L ||
       higher_ineligible_competitor ||
       higher_unreliable_competitor)
  output_competition_unresolved <- !is.na(selected_start) &&
    (selected_snapshot_missing ||
       best_basin_output_support < 2L ||
       within_basin_output_disagreement)
  conditional_initialization_stable <- !is.na(selected_start) &&
    !objective_competition_unresolved &&
    !output_competition_unresolved
  exploration_sufficient <-
    nrow(summary_table) >= min_exploration_starts
  initialization_stable <-
    conditional_initialization_stable && exploration_sufficient

  unresolved_reasons <- character()
  if (is.na(selected_start)) {
    unresolved_reasons <- "no_valid_endpoint"
  } else {
    if (selected_endpoint_unfinished) {
      unresolved_reasons <- c(
        unresolved_reasons, "selected_endpoint_unfinished")
    }
    if (best_basin_support < 2L) {
      unresolved_reasons <- c(
        unresolved_reasons, "insufficient_objective_support")
    }
    if (selected_snapshot_missing ||
        best_basin_output_support < 2L) {
      unresolved_reasons <- c(
        unresolved_reasons, "insufficient_output_support")
    }
    if (within_basin_output_disagreement) {
      unresolved_reasons <- c(
        unresolved_reasons, "within_basin_output_disagreement")
    }
    if (higher_ineligible_competitor) {
      unresolved_reasons <- c(
        unresolved_reasons, "higher_unfinished_competitor")
    }
    if (higher_unreliable_competitor) {
      unresolved_reasons <- c(
        unresolved_reasons, "higher_unreliable_competitor")
    }
  }
  if (!exploration_sufficient) {
    unresolved_reasons <- c(
      unresolved_reasons, "insufficient_exploration_starts")
  }
  initialization_status <- if (initialization_stable) {
    "stable"
  } else if (is.na(selected_start)) {
    "no_valid_endpoint"
  } else if (selected_endpoint_unfinished) {
    "best_valid_unfinished"
  } else {
    "best_converged_unresolved"
  }

  diagnostic_finite <- which(is.finite(summary_table$final_elbo))
  diagnostic_best_start <- if (length(diagnostic_finite)) {
    diagnostic_finite[
      which.max(summary_table$final_elbo[diagnostic_finite])]
  } else {
    NA_integer_
  }

  list(
    summary_table = summary_table,
    selected_start = selected_start,
    selected_endpoint_unfinished = selected_endpoint_unfinished,
    best_basin_support = best_basin_support,
    best_basin_output_support = best_basin_output_support,
    higher_ineligible_competitor = higher_ineligible_competitor,
    higher_ineligible_start_ids = higher_ineligible_start_ids,
    higher_ineligible_seeds = higher_ineligible_seeds,
    higher_unreliable_competitor = higher_unreliable_competitor,
    higher_unreliable_start_ids = higher_unreliable_start_ids,
    higher_unreliable_seeds = higher_unreliable_seeds,
    within_basin_output_disagreement =
      within_basin_output_disagreement,
    objective_competition_unresolved =
      objective_competition_unresolved,
    output_competition_unresolved =
      output_competition_unresolved,
    conditional_initialization_stable =
      conditional_initialization_stable,
    min_exploration_starts = min_exploration_starts,
    exploration_sufficient = exploration_sufficient,
    initialization_stable = initialization_stable,
    initialization_status = initialization_status,
    unresolved_reasons = unresolved_reasons,
    diagnostic_best_start = diagnostic_best_start
  )
}

.multistart_stage_row <- function(
    state, stage_id, attempted_starts, final_stage,
    decision_override = NULL) {
  selected_start <- state$selected_start
  summary_table <- state$summary_table
  decision <- if (!is.null(decision_override)) {
    decision_override
  } else if (state$initialization_stable) {
    "stop_stable"
  } else if (final_stage) {
    "stop_budget_exhausted"
  } else {
    "expand"
  }
  data.frame(
    stage_id = as.integer(stage_id),
    attempted_starts = as.integer(attempted_starts),
    candidate_valid_starts = sum(summary_table$candidate_valid),
    eligible_starts = sum(summary_table$eligible),
    selected_start = selected_start,
    selected_seed = if (is.na(selected_start)) {
      NA_integer_
    } else {
      summary_table$fit_seed[selected_start]
    },
    selected_elbo = if (is.na(selected_start)) {
      NA_real_
    } else {
      summary_table$final_elbo[selected_start]
    },
    selected_endpoint_unfinished =
      state$selected_endpoint_unfinished,
    best_basin_support = state$best_basin_support,
    best_basin_output_support = state$best_basin_output_support,
    higher_ineligible_competitor =
      state$higher_ineligible_competitor,
    higher_ineligible_start_ids =
      paste(state$higher_ineligible_start_ids, collapse = ","),
    higher_ineligible_seeds =
      paste(state$higher_ineligible_seeds, collapse = ","),
    higher_unreliable_competitor =
      state$higher_unreliable_competitor,
    higher_unreliable_start_ids =
      paste(state$higher_unreliable_start_ids, collapse = ","),
    higher_unreliable_seeds =
      paste(state$higher_unreliable_seeds, collapse = ","),
    within_basin_output_disagreement =
      state$within_basin_output_disagreement,
    objective_competition_unresolved =
      state$objective_competition_unresolved,
    output_competition_unresolved =
      state$output_competition_unresolved,
    conditional_initialization_stable =
      state$conditional_initialization_stable,
    min_exploration_starts =
      state$min_exploration_starts,
    exploration_sufficient =
      state$exploration_sufficient,
    initialization_stable = state$initialization_stable,
    unresolved_reasons =
      paste(state$unresolved_reasons, collapse = " | "),
    decision = decision,
    stringsAsFactors = FALSE
  )
}

#' Fit multiFSYNC with reproducible multiple starts
#'
#' By default every pre-specified seed is fitted and the largest ordinary ELBO
#' among numerically valid endpoints is selected without simulated truth. A
#' valid endpoint that reaches the hard T=1 budget before practical stopping
#' remains selectable, but is returned with
#' \code{selected_endpoint_unfinished = TRUE} and can never be
#' production-ready. With \code{retry_only = TRUE}, an eligible first start
#' may be returned immediately; that fast path is explicitly marked as
#' unresolved because its selected basin has no independent support.
#'
#' @param ... Named arguments passed to \code{\link{bayesSYNC_multi}}. The
#'   convergence rule must be \code{"practical"} (it is supplied automatically
#'   when omitted).
#' @param start_seeds Explicit, distinct non-negative integer seeds.
#' @param retry_only If \code{TRUE}, return an eligible first start immediately;
#'   after an ineligible first start, run every supplied retry and select the
#'   eligible fit with the largest final ELBO. If \code{FALSE}, always run all
#'   supplied starts (the default and recommended stability audit).
#' @param keep_fits Retain all full fitted objects in the returned audit object.
#' @param basin_rel_tol Relative final-ELBO tolerance used to count independent
#'   support for the selected objective basin.
#' @param basin_abs_tol Absolute final-ELBO tolerance floor used together with
#'   `basin_rel_tol`. This prevents near-cancellation in the total ELBO from
#'   splitting scientifically indistinguishable solutions into artificial
#'   basins.
#' @param stage_sizes Optional strictly increasing cumulative start counts.
#'   For example, \code{c(3L, 5L, 10L)} evaluates stability after 3, 5, and
#'   10 starts, but the default minimum-exploration gate prevents stopping at
#'   3 or 5. The final count must equal \code{length(start_seeds)}. The default
#'   \code{NULL} runs every supplied start and preserves the original
#'   behaviour. Early stability is conditional on the attempted seeds and
#'   cannot exclude a better basin among untried seeds.
#' @param min_exploration_starts Minimum number of attempted starts required
#'   before a staged fit may be declared stable. The staged default is
#'   \code{min(8L, length(start_seeds))}; set this explicitly to a smaller value
#'   only for a diagnostic reproduction of the earlier conditional rule.
#'   Non-staged fits still run every requested seed.
#' @param start_initializations Optional character vector, with one entry per
#'   seed, choosing \code{"random"} or \code{"residual_fpca"} for each start.
#'   A length-one value is recycled. If omitted, a single
#'   \code{initialization} value supplied through \code{...} is recycled, and
#'   otherwise every start uses \code{"random"}. Deterministic residual-FPCA
#'   replicas share one independence identifier and therefore count only once
#'   toward basin support.
#'
#' @return A list containing \code{fit}, \code{multistart_summary}, selection
#'   metadata, objective- and output-competition diagnostics, and optionally
#'   every fitted object. \code{selected_endpoint_unfinished} distinguishes a
#'   valid hard-budget endpoint from a practically converged endpoint.
#'   \code{production_ready} is the stability decision
#'   conditional on the attempted seeds; \code{stage_history} records each
#'   boundary decision and \code{unresolved_reasons} explains a failed
#'   decision. \code{higher_ineligible_start_ids} contains summary-row/start
#'   indices, while \code{higher_ineligible_seeds} contains the corresponding
#'   seeds. The analogous \code{higher_unreliable_*} fields report finite but
#'   invalid anomalously high fits. \code{conditional_initialization_stable}
#'   reports the former basin/output decision before the minimum-exploration
#'   production gate. Output alignment uses the fitted-value, expected-RSS,
#'   and PPI tolerances in the resolved practical control.
#'
#' @export
bayesSYNC_multi_multistart <- function(
    ..., start_seeds,
    retry_only = FALSE,
    keep_fits = FALSE,
    basin_rel_tol = 1e-4,
    stage_sizes = NULL,
    start_initializations = NULL,
    basin_abs_tol = 0.05,
    min_exploration_starts = NULL) {
  dots <- list(...)
  if ("seed" %in% names(dots)) {
    stop("Supply initialization seeds only through start_seeds.")
  }
  if (!is.numeric(start_seeds) || !length(start_seeds) ||
      any(!is.finite(start_seeds)) ||
      any(!vapply(start_seeds, is_int, logical(1))) ||
      any(start_seeds < 0) ||
      any(start_seeds > .Machine$integer.max)) {
    stop("start_seeds must contain non-negative integers.")
  }
  start_seeds <- as.integer(start_seeds)
  if (anyDuplicated(start_seeds)) {
    stop("start_seeds must be distinct.")
  }
  initialization_from_dots <- NULL
  if ("initialization" %in% names(dots)) {
    if (!is.null(start_initializations)) {
      stop(
        "Supply initialization methods through either ",
        "start_initializations or initialization in ..., not both."
      )
    }
    initialization_from_dots <- dots$initialization
    dots$initialization <- NULL
  }
  if (is.null(start_initializations)) {
    start_initializations <- if (is.null(initialization_from_dots)) {
      "random"
    } else {
      initialization_from_dots
    }
  }
  if (!is.character(start_initializations) ||
      !length(start_initializations) ||
      anyNA(start_initializations) ||
      any(!start_initializations %in% c("random", "residual_fpca"))) {
    stop(
      "start_initializations must contain only 'random' or ",
      "'residual_fpca'."
    )
  }
  if (length(start_initializations) == 1L) {
    start_initializations <- rep(
      start_initializations, length(start_seeds))
  }
  if (length(start_initializations) != length(start_seeds)) {
    stop(
      "start_initializations must have length one or ",
      "length(start_seeds)."
    )
  }
  if (length(retry_only) != 1L || !is.logical(retry_only) ||
      is.na(retry_only)) {
    stop("retry_only must be TRUE or FALSE.")
  }
  if (length(keep_fits) != 1L || !is.logical(keep_fits) ||
      is.na(keep_fits)) {
    stop("keep_fits must be TRUE or FALSE.")
  }
  if (length(basin_rel_tol) != 1L || !is.finite(basin_rel_tol) ||
      basin_rel_tol < 0) {
    stop("basin_rel_tol must be one finite non-negative number.")
  }
  if (length(basin_abs_tol) != 1L || !is.finite(basin_abs_tol) ||
      basin_abs_tol < 0) {
    stop("basin_abs_tol must be one finite non-negative number.")
  }
  staged <- !is.null(stage_sizes)
  if (staged) {
    if (!is.numeric(stage_sizes) || !length(stage_sizes) ||
        any(!is.finite(stage_sizes)) ||
        any(!vapply(stage_sizes, is_int, logical(1))) ||
        any(stage_sizes < 1L)) {
      stop("stage_sizes must contain positive integers.")
    }
    stage_sizes <- as.integer(stage_sizes)
    if (stage_sizes[1L] < 3L ||
        any(diff(stage_sizes) <= 0L)) {
      stop("stage_sizes must be strictly increasing and start at 3 or later.")
    }
    if (utils::tail(stage_sizes, 1L) != length(start_seeds)) {
      stop("The final stage_sizes value must equal length(start_seeds).")
    }
    if (retry_only) {
      stop("stage_sizes cannot be combined with retry_only = TRUE.")
    }
  }
  if (is.null(min_exploration_starts)) {
    min_exploration_starts <- if (staged) {
      min(8L, length(start_seeds))
    } else {
      1L
    }
  } else {
    if (!is.numeric(min_exploration_starts) ||
        length(min_exploration_starts) != 1L ||
        !is.finite(min_exploration_starts) ||
        !is_int(min_exploration_starts) ||
        min_exploration_starts < 1L ||
        min_exploration_starts > length(start_seeds)) {
      stop(
        "min_exploration_starts must be one positive integer no larger ",
        "than length(start_seeds)."
      )
    }
    min_exploration_starts <- as.integer(min_exploration_starts)
  }
  if (is.null(dots$convergence_rule)) {
    dots$convergence_rule <- "practical"
  }
  if (!identical(dots$convergence_rule, "practical")) {
    stop("bayesSYNC_multi_multistart requires convergence_rule='practical'.")
  }
  if (!is.null(dots$lambda_orth) &&
      (length(dots$lambda_orth) != 1L ||
       !is.finite(dots$lambda_orth) ||
       dots$lambda_orth != 0)) {
    stop("Multi-start ELBO selection requires lambda_orth = 0.")
  }
  output_control <- .validate_practical_control(dots$practical_control)

  fits <- if (keep_fits) vector("list", length(start_seeds)) else NULL
  stability_snapshots <- vector("list", length(start_seeds))
  component_snapshots <- vector("list", length(start_seeds))
  summaries <- vector("list", length(start_seeds))
  stage_history_rows <- list()
  attempted <- 0L
  best_candidate_fit <- NULL
  best_candidate_start <- NA_integer_
  best_candidate_elbo <- -Inf

  for (start_id in seq_along(start_seeds)) {
    attempted <- attempted + 1L
    warnings_now <- character()
    error_now <- ""
    fit_now <- NULL
    elapsed <- system.time({
      fit_now <- tryCatch(
        withCallingHandlers(
          do.call(
            bayesSYNC_multi,
            c(
              dots,
              list(
                seed = start_seeds[start_id],
                initialization = start_initializations[start_id]
              )
            )
          ),
          warning = function(condition) {
            warnings_now <<- c(
              warnings_now, conditionMessage(condition))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(condition) {
          error_now <<- conditionMessage(condition)
          NULL
        }
      )
    })[["elapsed"]]
    # Single-bracket list assignment preserves an explicit NULL placeholder
    # when a start fails; [[<- NULL would delete the list element and shift
    # all subsequent start-to-fit indices.
    if (keep_fits) fits[start_id] <- list(fit_now)
    if (!is.null(fit_now) && !is.null(fit_now$stability_snapshot)) {
      stability_snapshots[[start_id]] <- fit_now$stability_snapshot
    }
    if (!is.null(fit_now) &&
        !is.null(fit_now$elbo_result$components)) {
      component_snapshots[[start_id]] <-
        fit_now$elbo_result$components
    }

    final_elbo <- if (!is.null(fit_now) &&
        length(fit_now$ELBO) &&
        all(is.finite(fit_now$ELBO))) {
      utils::tail(fit_now$ELBO, 1L)
    } else {
      NA_real_
    }
    objective_valid <- !is.null(fit_now) &&
      isTRUE(fit_now$elbo_objective_valid)
    no_decrease <- !is.null(fit_now) &&
      identical(as.integer(fit_now$ELBO_decrease_count), 0L)
    no_t1_jitter <- !is.null(fit_now) &&
      identical(as.integer(fit_now$elbo_t1_jitter_count), 0L)
    practical_ok <- !is.null(fit_now) &&
      isTRUE(fit_now$practical_converged) &&
      isTRUE(fit_now$converged)
    candidate_valid <- is.finite(final_elbo) &&
      objective_valid && no_decrease && no_t1_jitter
    eligible <- candidate_valid && practical_ok
    initialization_used <- if (
        !is.null(fit_now) &&
        length(fit_now$initialization) == 1L &&
        !is.na(fit_now$initialization)) {
      as.character(fit_now$initialization)
    } else {
      start_initializations[start_id]
    }
    initialization_independence_id <- if (
        !is.null(fit_now) &&
        length(fit_now$initialization_independence_id) == 1L &&
        !is.na(fit_now$initialization_independence_id) &&
        nzchar(fit_now$initialization_independence_id)) {
      as.character(fit_now$initialization_independence_id)
    } else {
      paste0(initialization_used, "_seed_", start_seeds[start_id])
    }
    initialization_diagnostics <- if (is.null(fit_now)) {
      NULL
    } else {
      fit_now$initialization_diagnostics
    }
    if (candidate_valid && final_elbo > best_candidate_elbo) {
      best_candidate_fit <- fit_now
      best_candidate_start <- start_id
      best_candidate_elbo <- final_elbo
    }

    stage_id_now <- if (staged) {
      which(start_id <= stage_sizes)[1L]
    } else {
      1L
    }
    summaries[[start_id]] <- data.frame(
      start_id = start_id,
      fit_seed = start_seeds[start_id],
      initialization = initialization_used,
      initialization_independence_id =
        initialization_independence_id,
      initialization_complete = if (is.null(fit_now)) {
        NA
      } else {
        isTRUE(initialization_diagnostics$complete)
      },
      initialization_fallback_count =
        if (is.null(initialization_diagnostics) ||
            is.null(initialization_diagnostics$fallback_count)) {
          NA_integer_
        } else {
          as.integer(initialization_diagnostics$fallback_count)
        },
      initialization_explained_fraction =
        if (is.null(initialization_diagnostics) ||
            is.null(initialization_diagnostics$explained_fraction)) {
          NA_real_
        } else {
          as.numeric(initialization_diagnostics$explained_fraction)
        },
      initialization_elapsed_seconds =
        if (is.null(initialization_diagnostics) ||
            is.null(initialization_diagnostics$elapsed_seconds)) {
          NA_real_
        } else {
          as.numeric(initialization_diagnostics$elapsed_seconds)
        },
      stage_id = as.integer(stage_id_now),
      success = !is.null(fit_now),
      candidate_valid = candidate_valid,
      eligible = eligible,
      converged = if (is.null(fit_now)) FALSE else isTRUE(fit_now$converged),
      practical_converged =
        if (is.null(fit_now)) FALSE else
          isTRUE(fit_now$practical_converged),
      slow_case =
        if (is.null(fit_now)) FALSE else isTRUE(fit_now$slow_case),
      objective_valid = objective_valid,
      elbo_decrease_count =
        if (is.null(fit_now)) NA_integer_ else
          as.integer(fit_now$ELBO_decrease_count),
      t1_jitter_count =
        if (is.null(fit_now)) NA_integer_ else
          as.integer(fit_now$elbo_t1_jitter_count),
      iterations_total =
        if (is.null(fit_now)) NA_integer_ else as.integer(fit_now$i_iter),
      annealing_sweeps =
        if (is.null(fit_now)) NA_integer_ else
          as.integer(fit_now$annealing_sweeps),
      t1_sweeps =
        if (is.null(fit_now)) NA_integer_ else
          as.integer(fit_now$t1_sweeps),
      final_elbo = final_elbo,
      expected_rss_total =
        if (is.null(fit_now)) NA_real_ else
          sum(fit_now$expected_rss_sum),
      elapsed_seconds = as.numeric(elapsed),
      warning_text = paste(unique(warnings_now), collapse = " | "),
      error_text = error_now,
      stringsAsFactors = FALSE
    )

    if (retry_only && start_id == 1L && eligible) break
    if (staged && start_id %in% stage_sizes) {
      stage_state <- .evaluate_multistart_state(
        do.call(rbind, summaries[seq_len(attempted)]),
        stability_snapshots[seq_len(attempted)],
        output_control,
        basin_rel_tol,
        basin_abs_tol,
        min_exploration_starts
      )
      stage_number <- match(start_id, stage_sizes)
      stage_history_rows[[length(stage_history_rows) + 1L]] <-
        .multistart_stage_row(
          stage_state,
          stage_id = stage_number,
          attempted_starts = attempted,
          final_stage = stage_number == length(stage_sizes)
        )
      if (stage_state$initialization_stable) break
    }
  }

  summary_table <- do.call(rbind, summaries[seq_len(attempted)])
  state <- .evaluate_multistart_state(
    summary_table,
    stability_snapshots[seq_len(attempted)],
    output_control,
    basin_rel_tol,
    basin_abs_tol,
    min_exploration_starts
  )
  summary_table <- state$summary_table
  legacy_summary_columns <- c(
    "start_id", "fit_seed", "success", "eligible", "converged",
    "practical_converged", "slow_case", "objective_valid",
    "elbo_decrease_count", "t1_jitter_count", "iterations_total",
    "annealing_sweeps", "t1_sweeps", "final_elbo",
    "expected_rss_total", "elapsed_seconds", "warning_text",
    "error_text", "selected", "basin_relative_gap",
    "objective_basin_member", "fitted_nrmse_to_selected",
    "rss_rel_to_selected", "ppi_max_abs_to_selected",
    "output_aligned"
  )
  new_summary_columns <- c(
    "candidate_valid",
    "initialization", "initialization_independence_id",
    "initialization_complete", "initialization_fallback_count",
    "initialization_explained_fraction",
    "initialization_elapsed_seconds",
    "stage_id", "basin_absolute_gap",
    "higher_than_selected_unfinished",
    "higher_than_selected_unreliable"
  )
  summary_table <- summary_table[
    c(legacy_summary_columns, new_summary_columns)]
  selected_start <- state$selected_start
  selected_fit <- NULL
  selection_reason <- "no_valid_endpoint"
  if (!is.na(selected_start)) {
    if (!identical(selected_start, best_candidate_start)) {
      stop("Internal multi-start selection state is inconsistent.")
    }
    selected_fit <- best_candidate_fit
    selection_reason <- if (attempted == 1L &&
                            !state$selected_endpoint_unfinished) {
      "first_start_practically_converged"
    } else if (state$selected_endpoint_unfinished) {
      "highest_elbo_valid_unfinished_endpoint"
    } else {
      "highest_elbo_among_converged_valid_starts"
    }
  }

  if (state$objective_competition_unresolved ||
       state$output_competition_unresolved) {
    warning(
      "The selected multi-start solution is not independently resolved ",
      "(selected endpoint unfinished=",
      state$selected_endpoint_unfinished,
      ", objective support=", state$best_basin_support,
      ", output-aligned support=", state$best_basin_output_support,
      ", higher unfinished competitor=",
      state$higher_ineligible_competitor,
      ", higher unreliable competitor=",
      state$higher_unreliable_competitor,
      ", output competition=", state$output_competition_unresolved,
      "). The highest-ELBO valid endpoint is returned with explicit ",
      "diagnostic flags."
    )
  } else if (staged && is.na(selected_start)) {
    warning(
      "The staged multi-start budget ended without a valid ordinary-ELBO ",
      "endpoint. No fit is returned; inspect stage_history and ",
      "multistart_summary."
    )
  }

  if (!length(stage_history_rows)) {
    nonstaged_decision <- if (
        retry_only && attempted == 1L && !is.na(selected_start)) {
      "first_eligible_fast_path"
    } else {
      "all_requested_starts_completed"
    }
    stage_history_rows[[1L]] <- .multistart_stage_row(
      state,
      stage_id = 1L,
      attempted_starts = attempted,
      final_stage = TRUE,
      decision_override = nonstaged_decision
    )
  }
  stage_history <- do.call(rbind, stage_history_rows)
  completed_stage <- if (staged) {
    match(attempted, stage_sizes)
  } else {
    1L
  }
  multistart_stop_reason <- if (staged) {
    if (state$initialization_stable) {
      paste0("stable_after_", attempted, "_starts")
    } else if (is.na(selected_start)) {
      "budget_exhausted_no_valid_endpoint"
    } else {
      "budget_exhausted_unresolved"
    }
  } else if (retry_only && attempted == 1L &&
             !is.na(selected_start)) {
    "first_eligible_fast_path"
  } else {
    "all_requested_starts_completed"
  }

  component_rows <- lapply(seq_len(attempted), function(start_id) {
    components <- component_snapshots[[start_id]]
    if (is.null(components) || !length(components)) return(NULL)
    data.frame(
      start_id = start_id,
      fit_seed = start_seeds[start_id],
      initialization =
        summary_table$initialization[start_id],
      initialization_independence_id =
        summary_table$initialization_independence_id[start_id],
      component = names(components),
      value = as.numeric(components),
      stringsAsFactors = FALSE
    )
  })
  component_rows <- component_rows[
    !vapply(component_rows, is.null, logical(1))]
  multistart_elbo_components <- if (length(component_rows)) {
    do.call(rbind, component_rows)
  } else {
    data.frame(
      start_id = integer(),
      fit_seed = integer(),
      initialization = character(),
      initialization_independence_id = character(),
      component = character(),
      value = double(),
      stringsAsFactors = FALSE
    )
  }

  result <- list(
    fit = selected_fit,
    multistart_summary = summary_table,
    attempted_starts = attempted,
    selected_start = selected_start,
    selected_seed = if (is.na(selected_start)) {
      NA_integer_
    } else {
      summary_table$fit_seed[selected_start]
    },
    selection_reason = selection_reason,
    initialization_unresolved = is.null(selected_fit),
    selected_endpoint_unfinished =
      state$selected_endpoint_unfinished,
    initialization_stable = state$initialization_stable,
    conditional_initialization_stable =
      state$conditional_initialization_stable,
    min_exploration_starts = min_exploration_starts,
    exploration_sufficient = state$exploration_sufficient,
    objective_competition_unresolved =
      state$objective_competition_unresolved,
    output_competition_unresolved =
      state$output_competition_unresolved,
    best_basin_support = state$best_basin_support,
    best_basin_output_support = state$best_basin_output_support,
    higher_ineligible_competitor =
      state$higher_ineligible_competitor,
    higher_unreliable_competitor =
      state$higher_unreliable_competitor,
    basin_rel_tol = basin_rel_tol,
    basin_abs_tol = basin_abs_tol,
    diagnostic_best_start = state$diagnostic_best_start,
    start_seeds = start_seeds,
    start_initializations = start_initializations,
    initialization_status = state$initialization_status,
    production_ready = state$initialization_stable,
    unresolved_reasons = state$unresolved_reasons,
    within_basin_output_disagreement =
      state$within_basin_output_disagreement,
    higher_ineligible_start_ids =
      state$higher_ineligible_start_ids,
    higher_ineligible_seeds =
      state$higher_ineligible_seeds,
    higher_unreliable_start_ids =
      state$higher_unreliable_start_ids,
    higher_unreliable_seeds =
      state$higher_unreliable_seeds,
    attempted_start_seeds = start_seeds[seq_len(attempted)],
    attempted_start_initializations =
      start_initializations[seq_len(attempted)],
    stage_sizes = stage_sizes,
    stage_history = stage_history,
    completed_stage = completed_stage,
    multistart_stop_reason = multistart_stop_reason,
    multistart_elbo_components = multistart_elbo_components
  )
  if (keep_fits) {
    result$fits <- fits[seq_len(attempted)]
  }
  class(result) <- c("multiFSYNC_multistart", "list")
  result
}
