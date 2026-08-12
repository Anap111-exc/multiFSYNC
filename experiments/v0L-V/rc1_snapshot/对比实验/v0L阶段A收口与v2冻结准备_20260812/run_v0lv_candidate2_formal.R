# Formal candidate2 execution entry point. The current candidate binding is
# deliberately unauthorized, so every mutating action refuses until a separate
# frozen binding is created and explicitly authorized by the user.

candidate2_formal_script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(argument)) stop("Run this entry point with Rscript.")
  normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/",
                mustWork = TRUE)
}

candidate2_parse_cli <- function(arguments) {
  result <- list()
  for (argument in arguments) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Arguments must use --key=value syntax.")
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = "=")
    if (!nzchar(key) || key %in% names(result)) stop("Duplicate/empty CLI key.")
    result[[key]] <- value
  }
  result
}

script_path <- candidate2_formal_script_path()
base_dir <- dirname(script_path)
project_root <- normalizePath(file.path(base_dir, "..", ".."),
                              winslash = "/", mustWork = TRUE)
source(file.path(base_dir, "v0lv_candidate2_runtime.R"), local = FALSE)
source(file.path(base_dir, "r_route_v3_candidate2_runner.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_controller.R"), local = FALSE)
source(file.path(base_dir, "v0lv_candidate2_evaluation.R"), local = FALSE)

cli <- candidate2_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "status"
binding_dir <- normalizePath(
  cli$`binding-dir` %||%
    file.path(base_dir, "protocol_manifests_20260812_v3_candidate2"),
  winslash = "/", mustWork = TRUE
)
formal_root <- normalizePath(
  cli$`formal-root` %||%
    file.path(base_dir, "v0lv_candidate2_formal_results_NOT_STARTED"),
  winslash = "/", mustWork = FALSE
)

candidate2_manifest_paths <- function() {
  list(
    data = file.path(binding_dir, "V0LV_DATA_SEEDS_CANDIDATE2.csv"),
    fit = file.path(binding_dir, "V0LV_FIT_MANIFEST_CANDIDATE2.csv"),
    audit = file.path(binding_dir, "V0LV_AUDIT_MANIFEST_CANDIDATE2.csv")
  )
}

candidate2_read_manifests <- function() {
  paths <- candidate2_manifest_paths()
  if (any(!file.exists(unlist(paths)))) stop("Formal manifest files are missing.")
  result <- lapply(paths, utils::read.csv, stringsAsFactors = FALSE,
                   check.names = FALSE)
  if (any(result$data$manifest_schema_id !=
          V0LV_CANDIDATE2_MANIFEST_SCHEMA_ID) ||
      any(result$fit$total_protocol_id !=
          V0LV_CANDIDATE2_TOTAL_PROTOCOL_ID) ||
      any(result$audit$r_only_protocol_id !=
          R_ROUTE_V3_CANDIDATE2_PROTOCOL_ID)) {
    stop("Formal manifest version identity mismatch.")
  }
  result
}

candidate2_formal_preflight <- function() {
  verification <- c2_verify_formal_binding(binding_dir, project_root)
  manifests <- candidate2_read_manifests()
  if (nrow(manifests$data) != 10L || nrow(manifests$fit) != 210L ||
      nrow(manifests$audit) != 72L) {
    stop("Formal manifest row-count contract failed.")
  }
  invisible(list(verification = verification, manifests = manifests))
}

candidate2_existing_terminal <- function(directory) {
  terminal <- file.path(directory, "terminal_record.rds")
  marker <- file.path(directory, "TERMINAL_COMPLETE.txt")
  if (!file.exists(terminal) && !file.exists(marker)) return(NULL)
  if (!file.exists(terminal) || !file.exists(marker)) {
    stop("Incomplete prior output requires manual quarantine: ", directory)
  }
  record <- readRDS(terminal); c2_validate_terminal_record(record)
  if (!paste0("terminal_record_sha256=", c2_sha256(terminal)) %in%
      readLines(marker, warn = FALSE)) stop("Existing terminal hash mismatch.")
  record
}

candidate2_formal_hashes <- function(observation_path) {
  list(
    input = c2_sha256(observation_path),
    config = c2_sha256(file.path(binding_dir, "V0LV_CONFIG_CANDIDATE2.csv")),
    source = c2_sha256(file.path(
      binding_dir, "FULL_RUNTIME_SOURCE_BINDINGS.csv"
    ))
  )
}

candidate2_formal_generate <- function(data_id) {
  preflight <- candidate2_formal_preflight()
  row <- preflight$manifests$data[preflight$manifests$data$data_id == data_id,
                                  , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown formal data_id.")
  directory <- file.path(formal_root, "data", data_id)
  marker <- file.path(directory, "DATA_SEAL_COMPLETE.txt")
  if (file.exists(marker)) {
    observation <- file.path(directory, "observation_bundle.rds")
    truth <- file.path(directory, "sealed_truth", "truth_bundle.rds")
    if (!file.exists(observation) || !file.exists(truth)) {
      stop("Existing data seal is incomplete.")
    }
    return(invisible(list(reused = TRUE, observation = observation, truth = truth)))
  }
  candidate2_generate_and_seal(
    data_id = data_id, data_seed = row$data_seed[[1L]], output_dir = directory,
    mode = "formal", binding_dir = binding_dir, project_root = project_root
  )
}

candidate2_formal_fit <- function(fit_id) {
  preflight <- candidate2_formal_preflight()
  row <- preflight$manifests$fit[preflight$manifests$fit$fit_id == fit_id,
                                 , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown formal fit_id.")
  output_dir <- file.path(formal_root, "fits", fit_id)
  existing <- candidate2_existing_terminal(output_dir)
  if (!is.null(existing)) return(invisible(list(reused = TRUE, record = existing)))
  observation_path <- file.path(formal_root, "data", row$data_id[[1L]],
                                "observation_bundle.rds")
  if (!file.exists(observation_path)) stop("Observation bundle is not generated.")
  observation <- readRDS(observation_path)
  hashes <- candidate2_formal_hashes(observation_path)
  result <- if (row$path_id[[1L]] == "R") {
    fit_R_route_v3_candidate2(
      Y = observation$Y, Z = observation$Z,
      time_obs = observation$time_obs,
      L_f = observation$dimensions$L_f, L_s = observation$dimensions$L_s,
      M_f = observation$dimensions$M_f, M_s = observation$dimensions$M_s,
      K = observation$dimensions$K, fit_id = row$fit_id[[1L]],
      data_id = row$data_id[[1L]], fit_seed = row$fit_seed[[1L]],
      seed_index = row$seed_index[[1L]],
      initialization_independence_id =
        row$initialization_independence_id[[1L]],
      mode = "formal", binding_dir = binding_dir,
      project_root = project_root, hashes = hashes
    )
  } else if (row$path_id[[1L]] %in% c("B0", "A")) {
    candidate2_run_multi_main(
      observation = observation, method = row$path_id[[1L]],
      fit_id = row$fit_id[[1L]], fit_seed = row$fit_seed[[1L]],
      seed_index = row$seed_index[[1L]],
      initialization_independence_id =
        row$initialization_independence_id[[1L]],
      mode = "formal", binding_dir = binding_dir,
      project_root = project_root, hashes = hashes
    )
  } else if (row$path_id[[1L]] == "pooled_bayesSYNC") {
    candidate2_run_pooled_main(
      observation = observation, fit_id = row$fit_id[[1L]],
      fit_seed = row$fit_seed[[1L]], seed_index = row$seed_index[[1L]],
      initialization_independence_id =
        row$initialization_independence_id[[1L]], mode = "formal",
      binding_dir = binding_dir, project_root = project_root, hashes = hashes
    )
  } else stop("Unsupported formal path_id.")
  c2_write_terminal_bundle(result$terminal_record, result$fit, output_dir)
}

candidate2_load_fit_records <- function(rows) {
  records <- lapply(rows$fit_id, function(fit_id) {
    directory <- file.path(formal_root, "fits", fit_id)
    record <- candidate2_existing_terminal(directory)
    if (is.null(record)) stop("Required terminal is missing: ", fit_id)
    record
  })
  names(records) <- rows$fit_id
  records
}

candidate2_formal_select <- function(data_id) {
  preflight <- candidate2_formal_preflight()
  rows <- preflight$manifests$fit[
    preflight$manifests$fit$data_id == data_id, , drop = FALSE
  ]
  if (nrow(rows) != 21L) stop("Formal data has the wrong fit manifest size.")
  root <- file.path(formal_root, "selection", data_id)
  if (file.exists(file.path(root, "all", "TRUTH_FREE_SELECTION_FROZEN.txt"))) {
    return(invisible(list(reused = TRUE)))
  }
  components <- list()
  for (method in c("B0", "A")) {
    method_rows <- rows[rows$path_id == method, , drop = FALSE]
    records <- candidate2_load_fit_records(method_rows)
    directory <- file.path(root, method)
    components[[method]] <- directory
    candidate2_freeze_generic_truth_free_selection(
      rows, records, method = method, output_dir = directory,
      expected_starts = 4L, formal_experiment = TRUE
    )
  }
  r_rows <- rows[rows$path_id == "R", , drop = FALSE]
  r_records <- candidate2_load_fit_records(r_rows)
  components$R <- file.path(root, "R")
  candidate2_freeze_truth_free_selection(
    r_rows, r_records, components$R,
    expected_starts = 12L, formal_experiment = TRUE
  )
  pooled_rows <- rows[rows$path_id == "pooled_bayesSYNC", , drop = FALSE]
  pooled_record <- candidate2_load_fit_records(pooled_rows)[[1L]]
  pooled_summary <- c2_terminal_summary(pooled_record)
  pooled_manifest_key <- paste(
    pooled_rows$fit_id, pooled_rows$data_id, pooled_rows$fit_seed,
    pooled_rows$seed_index, pooled_rows$initialization_independence_id,
    sep = "\r"
  )
  pooled_terminal_key <- paste(
    pooled_summary$fit_id, pooled_summary$data_id, pooled_summary$fit_seed,
    pooled_summary$seed_index, pooled_summary$initialization_independence_id,
    sep = "\r"
  )
  if (!identical(pooled_manifest_key, pooled_terminal_key)) {
    stop("Pooled terminal does not exactly join to the manifest.")
  }
  candidate2_freeze_combined_truth_free_selection(
    components, pooled_record, file.path(root, "all"),
    formal_experiment = TRUE
  )
}

candidate2_formal_offline_racing <- function(data_id) {
  preflight <- candidate2_formal_preflight()
  rows <- preflight$manifests$fit[
    preflight$manifests$fit$data_id == data_id &
      preflight$manifests$fit$path_id == "R", , drop = FALSE
  ]
  records <- candidate2_load_fit_records(rows)
  summaries <- do.call(rbind, lapply(records, c2_terminal_summary))
  selection_path <- file.path(formal_root, "selection", data_id, "R",
                              "TRUTH_FREE_SELECTION.csv")
  if (!file.exists(selection_path)) stop("R truth-free winner is not frozen.")
  winner <- utils::read.csv(selection_path, stringsAsFactors = FALSE)
  traces <- lapply(rows$fit_id, function(fit_id) {
    fit <- readRDS(file.path(formal_root, "fits", fit_id, "fit.rds"))
    data.frame(t1_sweep = seq_along(fit$ELBO), elbo = as.numeric(fit$ELBO))
  }); names(traces) <- rows$fit_id
  replay <- r_route_v3_candidate2_offline_racing_replay(
    summaries, traces, winner, screen_t1 = 100L,
    minimum_keep = 4L, elbo_margin = 500
  )
  directory <- file.path(formal_root, "selection", data_id, "offline_racing")
  if (dir.exists(directory)) stop("Offline racing output already exists.")
  dir.create(directory, recursive = TRUE)
  c2_atomic_write_csv(replay$per_fit, file.path(directory, "PER_FIT.csv"))
  c2_atomic_write_csv(replay$per_data, file.path(directory, "PER_DATA.csv"))
}

candidate2_formal_audit <- function(audit_id) {
  preflight <- candidate2_formal_preflight()
  row <- preflight$manifests$audit[
    preflight$manifests$audit$audit_id == audit_id, , drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Unknown formal audit_id.")
  output_dir <- file.path(formal_root, "audits", audit_id)
  existing <- candidate2_existing_terminal(output_dir)
  if (!is.null(existing)) return(invisible(list(reused = TRUE, record = existing)))
  base_fit_id <- sprintf("%s_R_%02d", row$data_id[[1L]], row$seed_index[[1L]])
  source_id <- if (row$audit_stage[[1L]] == "anchor_200") base_fit_id else
    if (row$audit_stage[[1L]] == "to_400") {
      sprintf("%s_R_%02d_anchor_200", row$data_id[[1L]], row$seed_index[[1L]])
    } else sprintf("%s_R_%02d_to_400", row$data_id[[1L]], row$seed_index[[1L]])
  source_dir <- if (row$audit_stage[[1L]] == "anchor_200") {
    file.path(formal_root, "fits", source_id)
  } else file.path(formal_root, "audits", source_id)
  source_record <- candidate2_existing_terminal(source_dir)
  if (is.null(source_record)) stop("Audit source terminal is missing.")
  source_fit_path <- file.path(source_dir, "fit.rds")
  source_fit <- if (file.exists(source_fit_path)) readRDS(source_fit_path) else NULL
  observation_path <- file.path(formal_root, "data", row$data_id[[1L]],
                                "observation_bundle.rds")
  observation <- readRDS(observation_path)
  candidate2_run_R_audit_horizon(
    observation, source_fit, source_record,
    target_cumulative_T1 = row$target_cumulative_T1[[1L]],
    audit_id = audit_id, output_dir = output_dir, mode = "formal",
    binding_dir = binding_dir, project_root = project_root,
    hashes = candidate2_formal_hashes(observation_path)
  )
}

candidate2_formal_audit_freeze <- function(data_id, target) {
  preflight <- candidate2_formal_preflight()
  rows <- preflight$manifests$audit[
    preflight$manifests$audit$data_id == data_id &
      preflight$manifests$audit$target_cumulative_T1 == target, , drop = FALSE
  ]
  records <- lapply(rows$audit_id, function(id) {
    value <- candidate2_existing_terminal(file.path(formal_root, "audits", id))
    if (is.null(value)) stop("Audit terminal missing: ", id)
    value
  }); names(records) <- rows$audit_id
  fits <- lapply(rows$audit_id, function(id) {
    path <- file.path(formal_root, "audits", id, "fit.rds")
    if (file.exists(path)) readRDS(path) else NULL
  }); names(fits) <- rows$audit_id
  candidate2_freeze_audit_horizon_truth_free(
    rows, records, fits, target,
    file.path(formal_root, "audit_selection", data_id, paste0("T1_", target)),
    expected_starts = 12L, formal_experiment = TRUE
  )
}

candidate2_formal_authorize_unseal <- function(data_id, operator_authorized) {
  candidate2_formal_preflight()
  candidate2_authorize_unseal(
    selection_dir = file.path(formal_root, "selection", data_id, "all"),
    sealed_truth_path = file.path(formal_root, "data", data_id, "sealed_truth",
                                  "truth_bundle.rds"),
    authorization_dir = file.path(formal_root, "unseal_authorization", data_id),
    formal_experiment = TRUE, operator_authorized = operator_authorized
  )
}

candidate2_formal_evaluate <- function(fit_id) {
  preflight <- candidate2_formal_preflight()
  row <- preflight$manifests$fit[preflight$manifests$fit$fit_id == fit_id,
                                 , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown evaluation fit_id.")
  data_id <- row$data_id[[1L]]
  truth <- candidate2_unseal_truth(
    file.path(formal_root, "data", data_id, "sealed_truth", "truth_bundle.rds"),
    file.path(formal_root, "selection", data_id, "all"),
    file.path(formal_root, "unseal_authorization", data_id)
  )
  fit <- readRDS(file.path(formal_root, "fits", fit_id, "fit.rds"))
  evaluation <- if (row$path_id[[1L]] == "pooled_bayesSYNC") {
    candidate2_evaluate_pooled(
      fit, truth, project_root, mode = "formal", binding_dir = binding_dir
    )
  } else candidate2_evaluate_multi(
    fit, truth, project_root, mode = "formal", binding_dir = binding_dir
  )
  directory <- file.path(formal_root, "evaluation", fit_id)
  if (dir.exists(directory)) stop("Evaluation output already exists.")
  dir.create(directory, recursive = TRUE)
  c2_atomic_save_rds(evaluation, file.path(directory, "evaluation.rds"))
  for (name in names(evaluation)) if (is.data.frame(evaluation[[name]])) {
    c2_atomic_write_csv(evaluation[[name]],
                        file.path(directory, paste0(name, ".csv")))
  }
}

candidate2_formal_summarize <- function() {
  candidate2_formal_preflight()
  path <- file.path(formal_root, "evaluation", "PAIRED_METRICS.csv")
  if (!file.exists(path)) stop("PAIRED_METRICS.csv has not been assembled.")
  summary <- candidate2_summarize_paired_differences(
    utils::read.csv(path, stringsAsFactors = FALSE)
  )
  c2_atomic_write_csv(
    summary, file.path(formal_root, "evaluation", "SUCCESS_GATE_SUMMARY.csv")
  )
}

if (action == "status") {
  status <- utils::read.csv(file.path(binding_dir, "BINDING_STATUS.csv"),
                            stringsAsFactors = FALSE)
  cat("binding_status=", status$status[[1L]], "\n",
      "formal_execution_authorized=", status$formal_execution_authorized[[1L]],
      "\nformal_root=", formal_root, "\n", sep = "")
} else if (action == "generate") {
  candidate2_formal_generate(c2_scalar_string(cli$id, "--id"))
} else if (action == "fit") {
  candidate2_formal_fit(c2_scalar_string(cli$id, "--id"))
} else if (action == "select") {
  candidate2_formal_select(c2_scalar_string(cli$id, "--id"))
} else if (action == "offline-racing") {
  candidate2_formal_offline_racing(c2_scalar_string(cli$id, "--id"))
} else if (action == "audit") {
  candidate2_formal_audit(c2_scalar_string(cli$id, "--id"))
} else if (action == "audit-freeze") {
  candidate2_formal_audit_freeze(
    c2_scalar_string(cli$id, "--id"),
    c2_integer(as.numeric(cli$target), "--target", 1L)
  )
} else if (action == "authorize-unseal") {
  candidate2_formal_authorize_unseal(
    c2_scalar_string(cli$id, "--id"),
    identical(toupper(cli$`operator-authorized` %||% "FALSE"), "TRUE")
  )
} else if (action == "evaluate") {
  candidate2_formal_evaluate(c2_scalar_string(cli$id, "--id"))
} else if (action == "summarize") {
  candidate2_formal_summarize()
} else stop("Unknown action: ", action)
