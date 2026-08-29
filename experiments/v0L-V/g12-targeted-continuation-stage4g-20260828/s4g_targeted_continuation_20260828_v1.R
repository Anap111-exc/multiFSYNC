#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
))
stage4e_root <- Sys.getenv(
  "STAGE4G_STAGE4E_ROOT",
  unset = "/root/v0lv-g12-focused-reachability-stage4e-20260827-v1"
)
snapshot_root <- Sys.getenv(
  "STAGE4G_SNAPSHOT_ROOT",
  unset = "/root/multiFSYNC-git-rc1/experiments/v0L-V/rc1_snapshot"
)
development_repo <- Sys.getenv(
  "STAGE4G_DEVELOPMENT_REPO",
  unset = "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5"
)
development_library <- Sys.getenv(
  "STAGE4G_DEVELOPMENT_LIBRARY",
  unset = "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5-lib"
)
expected_commit <- "72d9a53f0ef5e9d3e5f49d9cf837958207b3c980"
source_t1 <- 400L
additional_t1 <- 400L
local_checkpoints <- seq.int(20L, 400L, 20L)
cumulative_checkpoints <- c(400L, 400L + local_checkpoints)

abort <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(value, message) if (!isTRUE(value)) abort(message)
`%||%` <- function(left, right) if (is.null(left)) right else left
iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}
read_csv <- function(path) {
  assert(file.exists(path), paste0("Missing input: ", path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
sha256 <- function(path) {
  assert(file.exists(path), paste0("Missing source: ", path))
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ))
  }
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  assert(is.null(status) || status == 0L, "sha256sum failed.")
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
find_one <- function(path, pattern) {
  result <- list.files(
    path, pattern = pattern, recursive = TRUE, full.names = TRUE
  )
  assert(length(result) == 1L,
         paste0("Expected one ", pattern, "; found ", length(result), "."))
  normalizePath(result[[1L]], winslash = "/", mustWork = TRUE)
}
parse_cli <- function(arguments) {
  result <- list()
  for (item in arguments) {
    assert(startsWith(item, "--") && grepl("=", item, fixed = TRUE),
           "Arguments must use --key=value syntax.")
    pieces <- strsplit(sub("^--", "", item), "=", fixed = TRUE)[[1L]]
    result[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
  }
  result
}
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
check_only <- identical(tolower(cli$`check-only` %||% "false"), "true")

assert(dir.exists(stage4e_root) && dir.exists(snapshot_root) &&
         dir.exists(development_repo) && dir.exists(development_library),
       "A bound Stage-4E source, snapshot, repo or library is missing.")
assert(file.exists(file.path(stage4e_root, "STAGE4E_COMPLETE.txt")),
       "Stage 4E is incomplete.")
manifest <- read_csv(file.path(root, "CONTINUATION_MANIFEST.csv"))
required_manifest <- c(
  "experiment_id", "continuation_id", "source_fit_id", "data_id",
  "scenario_id", "seed_index", "fit_seed", "selection_role",
  "adaptive_post_truth_target", "source_endpoint_role",
  "source_fit_relative_path", "source_fit_sha256",
  "observation_relative_path", "observation_sha256",
  "source_final_elbo", "source_objective_eligible",
  "source_structure_exact", "source_slow_case", "source_t1_sweeps",
  "additional_t1_sweeps", "cumulative_endpoint",
  "cumulative_checkpoints", "n_cpus", "outer_parallel_fit_limit",
  "truth_available_to_continuation", "winner_selection_reopened",
  "automatic_further_extension", "development_only",
  "formal_v0lv_result"
)
assert(
  nrow(manifest) == 2L && all(required_manifest %in% names(manifest)) &&
    !anyDuplicated(manifest$continuation_id) &&
    !anyDuplicated(manifest$source_fit_id) &&
    identical(as.character(manifest$source_fit_id), c(
      "g12ss4c_base_01__G12__04",
      "g12ss4c_sparse_02__G12__07"
    )) &&
    all(manifest$source_objective_eligible) &&
    all(manifest$source_structure_exact) && all(manifest$source_slow_case) &&
    all(manifest$source_t1_sweeps == source_t1) &&
    all(manifest$additional_t1_sweeps == additional_t1) &&
    all(manifest$cumulative_endpoint == source_t1 + additional_t1) &&
    all(manifest$n_cpus == 1L) &&
    all(manifest$outer_parallel_fit_limit == 2L) &&
    all(manifest$adaptive_post_truth_target) &&
    !any(manifest$truth_available_to_continuation) &&
    !any(manifest$winner_selection_reopened) &&
    !any(manifest$automatic_further_extension) &&
    all(manifest$development_only) && !any(manifest$formal_v0lv_result),
  "Stage-4G manifest identity, scope or fixed horizon is invalid."
)

source_paths <- file.path(stage4e_root, manifest$source_fit_relative_path)
observation_paths <- file.path(
  stage4e_root, manifest$observation_relative_path
)
assert(all(file.exists(source_paths)) && all(file.exists(observation_paths)),
       "A Stage-4G source fit or observation bundle is missing.")
assert(identical(
  unname(vapply(source_paths, sha256, character(1L))),
  as.character(manifest$source_fit_sha256)
) && identical(
  unname(vapply(observation_paths, sha256, character(1L))),
  as.character(manifest$observation_sha256)
), "A registered Stage-4G source hash changed.")

source_binding <- read_csv(file.path(root, "SOURCE_BINDING.csv"))
assert(nrow(source_binding) == 1L &&
         identical(source_binding$dev_commit_recorded[[1L]],
                   expected_commit),
       "Stage-4G source binding commit changed.")
bound_paths <- c(
  continuation_state_sha256 = "R/continuation_state.R",
  scale_trace_sha256 = "R/diagnostic_scale_trace.R",
  multi_core_sha256 = "R/multi_core.R",
  practical_stopping_sha256 = "R/convergence_practical.R",
  g12_stopping_control_sha256 = "R/g12_stopping_control.R",
  description_sha256 = "DESCRIPTION",
  namespace_sha256 = "NAMESPACE"
)
for (field in names(bound_paths)) {
  assert(identical(
    sha256(file.path(development_repo, bound_paths[[field]])),
    source_binding[[field]][[1L]]
  ), paste0("Stage-4G bound source changed: ", bound_paths[[field]]))
}
if (dir.exists(file.path(development_repo, ".git"))) {
  git_head <- system2(
    "git", c("-C", development_repo, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = TRUE
  )
  assert(length(git_head) == 1L && identical(git_head[[1L]], expected_commit),
         "Development package checkout is not the Stage-4E commit.")
  git_dirty <- system2(
    "git", c("-C", development_repo, "status", "--porcelain"),
    stdout = TRUE, stderr = TRUE
  )
  assert(!length(git_dirty), "Development package checkout is dirty.")
}

.libPaths(unique(c(development_library, .libPaths())))
assert(requireNamespace("multiFSYNC", quietly = TRUE),
       "Bound development multiFSYNC package is unavailable.")
assert(identical(
  normalizePath(find.package("multiFSYNC"), winslash = "/"),
  normalizePath(file.path(development_library, "multiFSYNC"),
                winslash = "/")
), "multiFSYNC did not load from the bound development library.")
source(find_one(snapshot_root, "^v0lv_candidate2_runtime[.]R$"),
       local = FALSE)
source(find_one(snapshot_root, "^r_route_v3_candidate2_runner[.]R$"),
       local = FALSE)

copy_matrix <- function(value) {
  value <- as.matrix(value)
  matrix(as.numeric(value), nrow = nrow(value), ncol = ncol(value))
}
score_kernel <- function(means, covariances) {
  means <- copy_matrix(means)
  result <- means %*% t(means)
  if (length(covariances)) {
    assert(length(covariances) == nrow(means),
           "Score covariance length differs from subjects.")
    diag(result) <- diag(result) + vapply(
      covariances, function(value) sum(diag(as.matrix(value))), numeric(1L)
    )
  }
  result
}
make_direction_snapshot <- function(cumulative_sweep, C_g, state,
                                    contributions) {
  list(
    cumulative_sweep = as.integer(cumulative_sweep),
    shared_loading = copy_matrix(state$mu_q_a),
    specific_loading = lapply(state$mu_q_b_specific, copy_matrix),
    shared_feature = lapply(
      state$mu_q_nu_phi,
      function(value) copy_matrix(C_g %*% as.matrix(value))
    ),
    specific_feature = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_nu_psi[[study]], function(value) {
        copy_matrix(C_g %*% as.matrix(value))
      })
    }),
    shared_score_mean = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_zeta[[study]], copy_matrix)
    }),
    specific_score_mean = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_xi[[study]], copy_matrix)
    }),
    shared_score_kernel = lapply(seq_len(state$S), function(study) {
      Map(score_kernel, state$mu_q_zeta[[study]],
          state$Sigma_q_zeta[[study]])
    }),
    specific_score_kernel = lapply(seq_len(state$S), function(study) {
      Map(score_kernel, state$mu_q_xi[[study]],
          state$Sigma_q_xi[[study]])
    }),
    shared_contribution = as.numeric(contributions$shared),
    specific_contribution = as.numeric(contributions$specific)
  )
}
snapshot_values <- function(item) {
  as.numeric(c(
    item$shared_loading,
    unlist(item$specific_loading, use.names = FALSE),
    unlist(item$shared_feature, use.names = FALSE),
    unlist(item$specific_feature, use.names = FALSE),
    unlist(item$shared_score_mean, use.names = FALSE),
    unlist(item$specific_score_mean, use.names = FALSE),
    unlist(item$shared_score_kernel, use.names = FALSE),
    unlist(item$specific_score_kernel, use.names = FALSE),
    item$shared_contribution, item$specific_contribution
  ))
}
with_private_traces <- function(code) {
  code <- substitute(code)
  namespace <- asNamespace("multiFSYNC")
  option_name <- get(
    ".multiFSYNC_scale_trace_option", envir = namespace, inherits = FALSE
  )
  validator <- get(
    ".validate_scale_trace_control", envir = namespace, inherits = FALSE
  )
  original <- get(
    ".make_scale_trace_row", envir = namespace, inherits = FALSE
  )
  contribution_function <- get(
    ".scale_trace_mean_contributions", envir = namespace, inherits = FALSE
  )
  trace_control <- validator(list(
    include_initial = TRUE, include_annealing = FALSE,
    t1_sweeps = local_checkpoints, block_annealing_sweeps = integer()
  ))
  capture <- new.env(parent = emptyenv())
  capture$snapshots <- list()
  wrapper <- function(...) {
    arguments <- list(...)
    result <- do.call(original, arguments)
    stage <- as.character(arguments$stage)
    local_sweep <- as.integer(arguments$t1_sweep)
    should_capture <- identical(stage, "initial") ||
      (identical(stage, "t1") && local_sweep %in% local_checkpoints)
    if (should_capture) {
      cumulative <- if (identical(stage, "initial")) {
        source_t1
      } else source_t1 + local_sweep
      contributions <- contribution_function(
        Y = arguments$Y, C = arguments$C, Z = arguments$Z,
        state = arguments$state
      )
      capture$snapshots[[as.character(cumulative)]] <-
        make_direction_snapshot(
          cumulative, arguments$C_g, arguments$state, contributions
        )
    }
    result
  }
  binding_locked <- bindingIsLocked(".make_scale_trace_row", namespace)
  old_options <- options()
  option_existed <- option_name %in% names(old_options)
  old_option <- old_options[[option_name]]
  if (binding_locked) unlockBinding(".make_scale_trace_row", namespace)
  assign(".make_scale_trace_row", wrapper, envir = namespace)
  if (binding_locked) lockBinding(".make_scale_trace_row", namespace)
  options(stats::setNames(list(trace_control), option_name))
  on.exit({
    if (bindingIsLocked(".make_scale_trace_row", namespace)) {
      unlockBinding(".make_scale_trace_row", namespace)
    }
    assign(".make_scale_trace_row", original, envir = namespace)
    if (binding_locked) lockBinding(".make_scale_trace_row", namespace)
    options(stats::setNames(
      list(if (option_existed) old_option else NULL), option_name
    ))
  }, add = TRUE)
  value <- eval(code, envir = parent.frame())
  list(fit = value, direction_snapshots = capture$snapshots)
}

validate_observation <- function(observation, data_id) {
  forbidden <- if (exists("candidate2_forbidden_observation_names",
                          mode = "function")) {
    candidate2_forbidden_observation_names(observation)
  } else {
    names(observation)[grepl("truth|true_params|signal_true|noise_true",
                             names(observation), ignore.case = TRUE)]
  }
  assert(
    identical(observation$bundle_class, "v0lv_observation_only_candidate2") &&
      identical(observation$data_id, data_id) && !length(forbidden),
    paste0("Observation bundle identity or truth isolation failed: ", data_id)
  )
  invisible(TRUE)
}

make_call_args <- function(observation, source_fit, state) {
  spec <- r_route_v3_candidate2_spec("formal")
  practical <- multiFSYNC::g12_stopping_control()
  practical$consecutive <- additional_t1 + 1L
  practical$checkpoints <- local_checkpoints
  list(
    Y = observation$Y, Z = observation$Z,
    time_obs = observation$time_obs,
    L_f = state$L_f, L_s = state$L_s,
    M_f = state$M_f, M_s = state$M_s, K = state$K,
    anneal = NULL, list_hyper = source_fit$list_hyper,
    n_g = length(state$time_g), time_g = state$time_g,
    tol_abs = spec$tol_abs, tol_rel = spec$tol_rel,
    maxit = additional_t1, n_cpus = 1L, verbose = FALSE,
    seed = as.integer(source_fit$fit_seed), bool_scale = FALSE,
    bool_var_spec_prob = state$bool_var_spec_prob,
    d_0 = as.integer(source_fit$d_0), convergence_rule = "practical",
    lambda_orth = 0, practical_control = practical,
    initialization = "random",
    initialization_control = source_fit$initialization_control,
    continuation_state = state
  )
}

terminal_summary <- function(record) {
  data.frame(
    continuation_id = record$continuation_id,
    source_fit_id = record$source_fit_id,
    data_id = record$data_id, scenario_id = record$scenario_id,
    seed_index = record$seed_index, fit_seed = record$fit_seed,
    selection_role = record$selection_role,
    terminal_status = record$terminal_status,
    source_t1_sweeps = record$source_t1_sweeps,
    continuation_t1_sweeps = record$continuation_t1_sweeps,
    cumulative_t1_sweeps = record$cumulative_t1_sweeps,
    source_final_elbo = record$source_final_elbo,
    final_elbo = record$final_elbo,
    total_elbo_gain = record$total_elbo_gain,
    bridge_difference = record$bridge_difference,
    bridge_monotone = record$bridge_monotone,
    source_snapshot_max_abs_difference =
      record$source_snapshot_max_abs_difference,
    objective_eligible = record$objective_eligible,
    warning_count = if (is.data.frame(record$warnings))
      nrow(record$warnings) else 0L,
    elapsed_seconds = record$elapsed_seconds,
    peak_memory_bytes = record$peak_memory_bytes,
    direction_snapshot_count = record$direction_snapshot_count,
    truth_used = record$truth_used,
    winner_selection_reopened = record$winner_selection_reopened,
    automatic_further_extension = record$automatic_further_extension,
    formal_v0lv_result = record$formal_v0lv_result,
    stringsAsFactors = FALSE
  )
}

run_one <- function(index) {
  row <- manifest[index, , drop = FALSE]
  continuation_id <- row$continuation_id[[1L]]
  output_dir <- file.path(root, "continuations", continuation_id)
  if (dir.exists(output_dir)) {
    assert(file.exists(file.path(output_dir, "FIT_COMPLETE.txt")) &&
             file.exists(file.path(output_dir, "terminal_record.rds")),
           paste0("Incomplete pre-existing continuation: ", continuation_id))
    return(readRDS(file.path(output_dir, "terminal_record.rds")))
  }
  dir.create(output_dir, recursive = TRUE, mode = "0700")
  source_fit <- readRDS(source_paths[[index]])
  source_directions <- readRDS(file.path(
    dirname(source_paths[[index]]), "DIRECTION_SNAPSHOTS.rds"
  ))
  observation <- readRDS(observation_paths[[index]])
  validate_observation(observation, row$data_id[[1L]])
  assert(
    identical(as.integer(source_fit$t1_sweeps), source_t1) &&
      isTRUE(all.equal(
        as.numeric(tail(source_fit$ELBO, 1L)),
        as.numeric(row$source_final_elbo[[1L]]), tolerance = 1e-10
      )) && "400" %in% names(source_directions),
    paste0("Source endpoint mismatch: ", row$source_fit_id[[1L]])
  )
  namespace <- asNamespace("multiFSYNC")
  make_state <- get(
    ".make_reduced_continuation_state", envir = namespace,
    inherits = FALSE
  )
  all_finite <- get(
    ".continuation_all_finite", envir = namespace, inherits = FALSE
  )
  state <- make_state(
    source_fit, shared_hat_indices = seq_len(source_fit$L_f)
  )
  assert(inherits(state, "multiFSYNC_continuation_state") &&
           all(vapply(state$parameters, all_finite, logical(1L))),
         paste0("Continuation state is invalid: ", continuation_id))

  fit_function <- getExportedValue("multiFSYNC", "bayesSYNC_multi")
  captured <- c2_capture_conditions(
    with_private_traces(do.call(
      fit_function, make_call_args(observation, source_fit, state)
    )), phase = "stage4g_fixed_400_to_800_continuation"
  )
  result <- captured$value
  fit <- if (is.null(result)) NULL else result$fit
  directions <- if (is.null(result)) list() else result$direction_snapshots
  if (!is.null(fit) && exists("c2_annotate_fit_warnings", mode = "function")) {
    captured$warnings <- c2_annotate_fit_warnings(captured$warnings, fit)
  }

  validation_error <- NULL
  bridge_difference <- bridge_tolerance <- snapshot_difference <- NA_real_
  bridge_monotone <- FALSE
  objective <- list(
    eligible = FALSE, checks = logical(), invalid_reasons = "terminal_error"
  )
  if (is.null(captured$error)) {
    validation_error <- tryCatch({
      assert(is.list(fit) && is.list(directions),
             "Continuation return lacks fit or direction snapshots.")
      assert(
        isTRUE(fit$continuation_diagnostics$used) &&
          identical(fit$initialization_diagnostics$used_method,
                    "continuation_state") &&
          fit$annealing_sweeps == 0L && fit$t1_sweeps == additional_t1 &&
          fit$n_cpus_used == 1L &&
          identical(as.numeric(fit$lambda_orth), 0) &&
          identical(fit$convergence_rule, "practical") &&
          fit$practical_control$consecutive == additional_t1 + 1L &&
          identical(names(fit$practical_checkpoints),
                    as.character(local_checkpoints)) &&
          nrow(fit$practical_diagnostics) == additional_t1,
        "Continuation violates the fixed 400-to-800 route."
      )
      assert(identical(names(directions),
                       as.character(cumulative_checkpoints)),
             "Continuation checkpoint sequence is incomplete.")
      source_values <- snapshot_values(source_directions[["400"]])
      continuation_values <- snapshot_values(directions[["400"]])
      assert(length(source_values) == length(continuation_values),
             "Source snapshot bridge shape changed.")
      snapshot_difference <- max(abs(
        source_values - continuation_values
      ))
      assert(is.finite(snapshot_difference) && snapshot_difference <= 1e-10,
             "Source snapshot bridge is not exact.")
      objective <- r_route_v3_candidate2_objective_status(fit)
      assert(objective$eligible, paste0(
        "Continuation objective is ineligible: ",
        paste(objective$invalid_reasons, collapse = ";")
      ))
      source_elbo <- as.numeric(row$source_final_elbo[[1L]])
      first_elbo <- as.numeric(fit$ELBO[[1L]])
      bridge_tolerance <- sqrt(.Machine$double.eps) *
        (1 + max(abs(c(source_elbo, first_elbo))))
      bridge_difference <- first_elbo - source_elbo
      assert(is.finite(bridge_difference) &&
               bridge_difference >= -bridge_tolerance,
             "Source-to-continuation ELBO bridge decreased.")
      bridge_monotone <- TRUE
      NULL
    }, error = function(condition) condition)
  }
  if (!is.null(validation_error)) {
    captured$error <- list(
      class = paste(class(validation_error), collapse = ";"),
      message = conditionMessage(validation_error),
      call = "stage4g continuation return validation",
      trace_summary = "s4g_targeted_continuation_20260828_v1"
    )
    fit <- NULL
    directions <- list()
  }

  terminal_status <- if (is.null(captured$error)) {
    "fixed_800_complete"
  } else "error"
  record <- list(
    experiment_id = "G12_TARGETED_CONTINUATION_STAGE4G_V1_20260828",
    continuation_id = continuation_id,
    source_fit_id = row$source_fit_id[[1L]],
    data_id = row$data_id[[1L]],
    scenario_id = row$scenario_id[[1L]],
    seed_index = as.integer(row$seed_index[[1L]]),
    fit_seed = as.integer(row$fit_seed[[1L]]),
    selection_role = row$selection_role[[1L]],
    adaptive_post_truth_target = TRUE,
    source_fit_sha256 = row$source_fit_sha256[[1L]],
    package_commit = expected_commit,
    started_at_utc = captured$started_at_utc,
    ended_at_utc = captured$ended_at_utc,
    elapsed_seconds = captured$elapsed_seconds,
    peak_memory_bytes = captured$peak_memory_bytes,
    terminal_status = terminal_status,
    error = captured$error %||% list(
      class = "", message = "", call = "", trace_summary = ""
    ),
    warnings = captured$warnings,
    source_t1_sweeps = source_t1,
    continuation_t1_sweeps = if (is.null(fit)) NA_integer_ else
      as.integer(fit$t1_sweeps),
    cumulative_t1_sweeps = if (is.null(fit)) NA_integer_ else
      source_t1 + as.integer(fit$t1_sweeps),
    source_final_elbo = as.numeric(row$source_final_elbo[[1L]]),
    first_continuation_elbo = if (is.null(fit)) NA_real_ else
      as.numeric(fit$ELBO[[1L]]),
    bridge_difference = bridge_difference,
    bridge_tolerance = bridge_tolerance,
    bridge_monotone = bridge_monotone,
    source_snapshot_max_abs_difference = snapshot_difference,
    final_elbo = if (is.null(fit)) NA_real_ else
      as.numeric(tail(fit$ELBO, 1L)),
    total_elbo_gain = if (is.null(fit)) NA_real_ else
      as.numeric(tail(fit$ELBO, 1L) - row$source_final_elbo[[1L]]),
    objective_eligible = isTRUE(objective$eligible),
    objective_checks = objective$checks,
    objective_invalid_reasons = objective$invalid_reasons,
    direction_snapshot_count = length(directions),
    truth_used = FALSE, winner_selection_reopened = FALSE,
    automatic_further_extension = FALSE, formal_v0lv_result = FALSE
  )
  saveRDS(record, file.path(output_dir, "terminal_record.rds"), version = 3)
  if (!is.null(fit)) {
    saveRDS(fit, file.path(output_dir, "fit.rds"), version = 3)
    saveRDS(directions,
            file.path(output_dir, "DIRECTION_SNAPSHOTS.rds"), version = 3)
    utils::write.csv(rbind(
      data.frame(
        stage = "stage4e_source", local_t1_sweep = seq_along(source_fit$ELBO),
        cumulative_t1_sweep = seq_along(source_fit$ELBO),
        elbo = as.numeric(source_fit$ELBO)
      ),
      data.frame(
        stage = "stage4g_continuation", local_t1_sweep = seq_along(fit$ELBO),
        cumulative_t1_sweep = source_t1 + seq_along(fit$ELBO),
        elbo = as.numeric(fit$ELBO)
      )
    ), file.path(output_dir, "ELBO_1_TO_800.csv"), row.names = FALSE,
    quote = TRUE)
    diagnostics <- fit$practical_diagnostics
    diagnostics$cumulative_t1_sweep <- source_t1 + diagnostics$t1_sweep
    utils::write.csv(
      diagnostics,
      file.path(output_dir, "PRACTICAL_DIAGNOSTICS_401_TO_800.csv"),
      row.names = FALSE, quote = TRUE, na = ""
    )
    scale <- fit$scale_trace
    scale$cumulative_t1_sweep <- ifelse(
      scale$stage == "initial", source_t1,
      source_t1 + scale$t1_sweep
    )
    utils::write.csv(
      scale, file.path(output_dir, "SCALE_TRACE_400_TO_800.csv"),
      row.names = FALSE, quote = TRUE, na = ""
    )
  }
  if (is.data.frame(record$warnings)) {
    utils::write.csv(
      record$warnings, file.path(output_dir, "WARNINGS.csv"),
      row.names = FALSE, quote = TRUE, na = ""
    )
  }
  if (terminal_status == "error") {
    writeLines(c(
      "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN",
      paste0("continuation_id=", continuation_id),
      paste0("message=", record$error$message),
      paste0("ended_utc=", iso_time())
    ), file.path(output_dir, "FIT_ERROR.txt"), useBytes = TRUE)
  } else {
    writeLines(c(
      "status=COMPLETE", paste0("continuation_id=", continuation_id),
      paste0("source_fit_id=", row$source_fit_id[[1L]]),
      "source_T1_sweeps=400", "continuation_T1_sweeps=400",
      "cumulative_T1_sweeps=800", "direction_snapshots=21",
      "truth_used=FALSE", "winner_selection_reopened=FALSE",
      "automatic_further_extension=FALSE", "formal_v0lv_result=FALSE",
      paste0("objective_eligible=", record$objective_eligible),
      paste0("ended_utc=", iso_time())
    ), file.path(output_dir, "FIT_COMPLETE.txt"), useBytes = TRUE)
  }
  record
}

if (check_only) {
  first <- manifest[1L, , drop = FALSE]
  source_fit <- readRDS(source_paths[[1L]])
  observation <- readRDS(observation_paths[[1L]])
  validate_observation(observation, first$data_id[[1L]])
  state <- getFromNamespace(
    ".make_reduced_continuation_state", "multiFSYNC"
  )(source_fit, shared_hat_indices = seq_len(source_fit$L_f))
  assert(inherits(state, "multiFSYNC_continuation_state"),
         "Check-only continuation state construction failed.")
  cat("STAGE4G_CHECK_ONLY_PASS continuations=2 source=400 target=800 truth=0\n")
  quit(save = "no", status = 0L)
}

assert(!file.exists(file.path(root, "CONTINUATIONS_COMPLETE.txt")),
       "Stage-4G continuations are already complete.")
dir.create(file.path(root, "continuations"), recursive = TRUE,
           mode = "0700")
records <- parallel::mclapply(
  seq_len(nrow(manifest)), run_one,
  mc.cores = 2L, mc.preschedule = FALSE, mc.set.seed = FALSE
)
valid_records <- vapply(records, is.list, logical(1L))
assert(all(valid_records), "A continuation worker did not return a record.")
statuses <- vapply(records, `[[`, character(1L), "terminal_status")
terminals <- do.call(rbind, lapply(records, terminal_summary))
terminals <- terminals[match(manifest$continuation_id,
                             terminals$continuation_id), , drop = FALSE]
utils::write.csv(
  terminals, file.path(root, "ALL_2_TERMINALS.csv"),
  row.names = FALSE, quote = TRUE, na = ""
)
if (any(statuses != "fixed_800_complete")) {
  writeLines(c(
    "status=ERROR_RETAINED_NO_AUTOMATIC_RERUN",
    paste0("completed_utc=", iso_time()),
    paste0("failed=", paste(
      manifest$continuation_id[statuses != "fixed_800_complete"],
      collapse = ";"
    ))
  ), file.path(root, "CONTINUATION_ERROR.txt"), useBytes = TRUE)
  abort("Stage-4G continuation error retained without selective rerun.")
}
assert(all(terminals$objective_eligible) && all(terminals$bridge_monotone) &&
         all(terminals$source_snapshot_max_abs_difference <= 1e-10) &&
         all(terminals$cumulative_t1_sweeps == 800L) &&
         !any(terminals$truth_used) &&
         !any(terminals$winner_selection_reopened) &&
         !any(terminals$automatic_further_extension),
       "A completed continuation violates the registered contract.")
writeLines(c(
  "status=CONTINUATIONS_COMPLETE", paste0("completed_utc=", iso_time()),
  "continuations=2", "source_T1_sweeps=400",
  "additional_T1_sweeps=400", "cumulative_T1_sweeps=800",
  "outer_workers=2", "per_fit_n_cpus=1", "truth_read=FALSE",
  "winner_selection_reopened=FALSE",
  "automatic_further_extension=FALSE", "formal_v0lv_result=FALSE"
), file.path(root, "CONTINUATIONS_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE4G_CONTINUATIONS_PASS fits=2 cumulative=800 truth=0\n")
