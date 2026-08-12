# Runtime, hashing, condition-capture and atomic-artifact helpers for the
# third-round v0L-V candidate. This file contains no model update equations.

V0LV_CANDIDATE2_TERMINAL_SCHEMA <- "v0lv_terminal_v3_candidate2_1.0.0"

c2_scalar_string <- function(value, name, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!allow_empty && !nzchar(value))) {
    stop(name, " must be one ", if (!allow_empty) "non-empty " else "",
         "string.")
  }
  value
}

c2_integer <- function(value, name, minimum = 0L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      abs(value - round(value)) > sqrt(.Machine$double.eps) ||
      value < minimum || value > .Machine$integer.max) {
    stop(name, " must be one integer of at least ", minimum, ".")
  }
  as.integer(round(value))
}

c2_iso_time <- function(value = Sys.time()) {
  format(value, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC", usetz = FALSE)
}

c2_sha256 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (file.info(path)$isdir) stop("SHA256 requires a file: ", path)
  if (file.info(path)$size == 0) {
    return("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
  }

  valid_hash <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) &&
      grepl("^[0-9a-f]{64}$", tolower(value), useBytes = TRUE)
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    candidate <- tryCatch(
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
      error = function(error) NA_character_
    )
    if (valid_hash(candidate)) return(tolower(candidate))
  }

  executable <- Sys.which("sha256sum")
  if (nzchar(executable)) {
    output <- tryCatch(
      system2(executable, shQuote(path), stdout = TRUE, stderr = TRUE),
      error = function(error) character()
    )
    candidate <- sub("[[:space:]].*$", "", trimws(output))
    candidate <- candidate[
      grepl("^[0-9A-Fa-f]{64}$", candidate, useBytes = TRUE)
    ]
    if (length(candidate)) return(tolower(candidate[[1L]]))
  }

  executable <- Sys.which("certutil")
  if (nzchar(executable)) {
    output <- tryCatch(
      system2(
        executable, c("-hashfile", shQuote(path), "SHA256"),
        stdout = TRUE, stderr = TRUE
      ),
      error = function(error) character()
    )
    candidate <- output[
      grepl("^[0-9A-Fa-f ]{64,}$", output, useBytes = TRUE)
    ]
    candidate <- tolower(gsub(" ", "", candidate, fixed = TRUE))
    candidate <- candidate[
      grepl("^[0-9a-f]{64}$", candidate, useBytes = TRUE)
    ]
    if (length(candidate)) return(candidate[[1L]])
  }

  stop(
    "No supported SHA256 backend succeeded for: ", path,
    ". Install the R package 'digest' or provide sha256sum/certutil."
  )
}

c2_relative_path <- function(path, project_root) {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(tolower(path), tolower(prefix))) {
    stop("Path is outside project root: ", path)
  }
  substring(path, nchar(prefix) + 1L)
}

c2_hash_table <- function(paths, project_root, role = "artifact") {
  normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  if (length(role) == 1L) {
    roles <- rep(as.character(role), length(normalized))
  } else if (length(role) == length(normalized)) {
    roles <- as.character(role)
  } else stop("role must be scalar or parallel to paths.")
  if (any(!nzchar(roles))) stop("Hash-table roles must be non-empty.")
  input <- data.frame(path = normalized, role = roles,
                      stringsAsFactors = FALSE)
  duplicated_paths <- unique(input$path[duplicated(input$path)])
  for (path in duplicated_paths) {
    if (length(unique(input$role[input$path == path])) != 1L) {
      stop("One source path has conflicting roles: ", path)
    }
  }
  input <- input[!duplicated(input$path), , drop = FALSE]
  input <- input[!file.info(input$path)$isdir, , drop = FALSE]
  input <- input[order(input$path), , drop = FALSE]
  paths <- input$path
  data.frame(
    role = input$role,
    relative_path = vapply(
      paths, c2_relative_path, character(1), project_root = project_root
    ),
    bytes = unname(file.info(paths)$size),
    sha256 = vapply(paths, c2_sha256, character(1)),
    stringsAsFactors = FALSE
  )
}

c2_verify_hash_index <- function(index, project_root) {
  required <- c("relative_path", "sha256")
  if (!is.data.frame(index) || !all(required %in% names(index)) ||
      !nrow(index) || anyDuplicated(index$relative_path)) {
    stop("Hash index is missing, empty, duplicated or malformed.")
  }
  paths <- file.path(project_root, index$relative_path)
  exists <- file.exists(paths) & !dir.exists(paths)
  actual <- rep(NA_character_, length(paths))
  actual[exists] <- vapply(paths[exists], c2_sha256, character(1))
  result <- data.frame(
    relative_path = as.character(index$relative_path),
    exists = exists,
    expected_sha256 = tolower(as.character(index$sha256)),
    actual_sha256 = actual,
    matches = exists & !is.na(actual) &
      actual == tolower(as.character(index$sha256)),
    stringsAsFactors = FALSE
  )
  result
}

c2_atomic_target <- function(path) {
  directory <- dirname(path)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  file.path(
    directory,
    paste0(".", basename(path), ".tmp_", Sys.getpid(), "_",
           format(Sys.time(), "%Y%m%d%H%M%OS6"))
  )
}

c2_atomic_commit <- function(temp, target) {
  if (file.exists(target)) stop("Refusing to overwrite artifact: ", target)
  if (!file.rename(temp, target)) {
    stop("Atomic rename failed from ", temp, " to ", target)
  }
  invisible(normalizePath(target, winslash = "/", mustWork = TRUE))
}

c2_atomic_save_rds <- function(value, path, compress = "xz") {
  if (file.exists(path)) stop("Refusing to overwrite artifact: ", path)
  temp <- c2_atomic_target(path)
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  saveRDS(value, temp, compress = compress)
  check <- readRDS(temp)
  if (!identical(value, check)) stop("Temporary RDS verification failed: ", path)
  hash <- c2_sha256(temp)
  c2_atomic_commit(temp, path)
  list(path = normalizePath(path, winslash = "/"), sha256 = hash,
       bytes = unname(file.info(path)$size))
}

c2_atomic_write_csv <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite artifact: ", path)
  temp <- c2_atomic_target(path)
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  utils::write.csv(value, temp, row.names = FALSE, na = "")
  check <- utils::read.csv(temp, check.names = FALSE,
                           stringsAsFactors = FALSE)
  if (nrow(check) != nrow(value) || ncol(check) != ncol(value) ||
      !identical(names(check), names(value))) {
    stop("Temporary CSV shape/name verification failed: ", path)
  }
  hash <- c2_sha256(temp)
  c2_atomic_commit(temp, path)
  list(path = normalizePath(path, winslash = "/"), sha256 = hash,
       bytes = unname(file.info(path)$size))
}

c2_atomic_write_lines <- function(value, path) {
  if (file.exists(path)) stop("Refusing to overwrite artifact: ", path)
  temp <- c2_atomic_target(path)
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  writeLines(as.character(value), temp, useBytes = TRUE)
  check <- readLines(temp, warn = FALSE)
  if (!identical(check, as.character(value))) {
    stop("Temporary text verification failed: ", path)
  }
  hash <- c2_sha256(temp)
  c2_atomic_commit(temp, path)
  list(path = normalizePath(path, winslash = "/"), sha256 = hash,
       bytes = unname(file.info(path)$size))
}

c2_peak_memory_bytes <- function(gc_result = gc()) {
  # R reports Ncells and Vcells. Their conventional byte widths are 56 and 8.
  column <- if ("max used" %in% colnames(gc_result)) "max used" else "used"
  cells <- as.numeric(gc_result[, column])
  as.numeric(cells[[1L]] * 56 + cells[[2L]] * 8)
}

c2_empty_warnings <- function() {
  data.frame(
    warning_index = integer(), condition_class = character(),
    message = character(), phase = character(), sweep = integer(),
    iteration = integer(), observed_at_utc = character(),
    elapsed_seconds = numeric(), call = character(),
    stringsAsFactors = FALSE
  )
}

c2_condition_call <- function(condition) {
  call <- conditionCall(condition)
  if (is.null(call)) "" else paste(deparse(call, width.cutoff = 200L),
                                    collapse = " ")
}

c2_capture_conditions <- function(expression, phase = "runtime") {
  phase <- c2_scalar_string(phase, "phase")
  warnings <- list()
  started_time <- Sys.time()
  started_elapsed <- proc.time()[["elapsed"]]
  gc(reset = TRUE)
  error <- NULL
  value <- tryCatch(
    withCallingHandlers(
      force(expression),
      warning = function(condition) {
        index <- length(warnings) + 1L
        warnings[[index]] <<- data.frame(
          warning_index = index,
          condition_class = paste(class(condition), collapse = ";"),
          message = conditionMessage(condition),
          phase = phase,
          sweep = NA_integer_, iteration = NA_integer_,
          observed_at_utc = c2_iso_time(),
          elapsed_seconds = proc.time()[["elapsed"]] - started_elapsed,
          call = c2_condition_call(condition), stringsAsFactors = FALSE
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      calls <- sys.calls()
      error <<- list(
        class = paste(class(condition), collapse = ";"),
        message = conditionMessage(condition),
        call = c2_condition_call(condition),
        trace_summary = paste(
          vapply(tail(calls, 12L), function(value) {
            paste(deparse(value, width.cutoff = 160L), collapse = " ")
          }, character(1)), collapse = " <- "
        )
      )
      NULL
    }
  )
  ended_time <- Sys.time()
  warning_table <- if (length(warnings)) do.call(rbind, warnings) else
    c2_empty_warnings()
  list(
    value = value, error = error, warnings = warning_table,
    started_at_utc = c2_iso_time(started_time),
    ended_at_utc = c2_iso_time(ended_time),
    elapsed_seconds = as.numeric(difftime(
      ended_time, started_time, units = "secs"
    )),
    peak_memory_bytes = c2_peak_memory_bytes(gc())
  )
}

c2_annotate_fit_warnings <- function(warnings, fit) {
  if (!is.data.frame(warnings) || !nrow(warnings) || !is.list(fit)) {
    return(warnings)
  }
  for (index in seq_len(nrow(warnings))) {
    message <- warnings$message[[index]]
    if (grepl("practical fit", message, fixed = TRUE)) {
      warnings$phase[[index]] <- "ordinary_T1_stopping"
      warnings$sweep[[index]] <- as.integer(fit$t1_sweeps)
      warnings$iteration[[index]] <- as.integer(fit$i_iter)
    } else if (grepl("Max iterations", message, fixed = TRUE)) {
      warnings$phase[[index]] <- "maximum_iteration_boundary"
      warnings$sweep[[index]] <- as.integer(fit$t1_sweeps)
      warnings$iteration[[index]] <- as.integer(fit$i_iter)
    } else if (grepl("ELBO", message, fixed = TRUE)) {
      warnings$phase[[index]] <- "objective_diagnostic"
      warnings$sweep[[index]] <- as.integer(fit$t1_sweeps)
      warnings$iteration[[index]] <- as.integer(fit$i_iter)
    }
  }
  warnings
}

c2_terminal_summary <- function(record) {
  warning_count <- if (is.data.frame(record$warnings)) nrow(record$warnings) else 0L
  data.frame(
    terminal_schema = as.character(record$terminal_schema),
    protocol_id = as.character(record$protocol_id %||% ""),
    runner_version = as.character(record$runner_version %||% ""),
    fit_id = as.character(record$fit_id),
    data_id = as.character(record$data_id),
    method = as.character(record$method),
    route = as.character(record$route),
    fit_seed = as.integer(record$fit_seed),
    seed_index = as.integer(record$seed_index),
    initialization_independence_id =
      as.character(record$initialization_independence_id),
    terminal_status = as.character(record$terminal_status),
    terminal_reason = as.character(record$terminal_reason),
    actual_annealing_sweeps = as.integer(record$actual_annealing_sweeps),
    actual_T1_sweeps = as.integer(record$actual_T1_sweeps),
    first_strict_convergence_sweep =
      as.integer(record$first_strict_convergence_sweep),
    strict_practical_converged = isTRUE(record$strict_practical_converged),
    objective_eligible = isTRUE(record$objective_eligible),
    objective_invalid_reasons = paste(
      as.character(record$objective_invalid_reasons), collapse = ";"
    ),
    final_elbo = as.numeric(record$final_elbo),
    selected_endpoint_unfinished =
      !isTRUE(record$strict_practical_converged),
    warning_count = warning_count,
    error_class = as.character(record$error$class %||% ""),
    error_message = as.character(record$error$message %||% ""),
    started_at_utc = as.character(record$started_at_utc),
    ended_at_utc = as.character(record$ended_at_utc),
    elapsed_seconds = as.numeric(record$elapsed_seconds),
    peak_memory_bytes = as.numeric(record$peak_memory_bytes),
    formal_experiment = isTRUE(record$formal_experiment),
    truth_used_for_fit_or_selection =
      isTRUE(record$truth_used_for_fit_or_selection),
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

c2_validate_terminal_record <- function(record) {
  required <- c(
    "terminal_schema", "protocol_id", "runner_version", "fit_id",
    "data_id", "method", "route",
    "fit_seed", "seed_index", "initialization_independence_id",
    "started_at_utc", "ended_at_utc", "elapsed_seconds",
    "peak_memory_bytes", "terminal_status", "terminal_reason", "error",
    "warnings", "actual_annealing_sweeps", "actual_T1_sweeps",
    "first_strict_convergence_sweep", "strict_practical_converged",
    "convergence_diagnostics", "objective_eligible",
    "objective_invalid_reasons", "final_elbo", "hashes",
    "formal_experiment", "truth_used_for_fit_or_selection"
  )
  if (!is.list(record) || !all(required %in% names(record))) {
    stop("Terminal record is incomplete; missing: ",
         paste(setdiff(required, names(record)), collapse = ", "))
  }
  c2_scalar_string(record$fit_id, "record$fit_id")
  c2_scalar_string(record$data_id, "record$data_id")
  c2_scalar_string(record$initialization_independence_id,
                   "record$initialization_independence_id")
  if (!record$terminal_status %in% c(
      "strict_practical_converged", "max_budget_reached", "error",
      "audit_horizon_reached", "audit_unavailable_due_to_main_error")) {
    stop("Unknown terminal status: ", record$terminal_status)
  }
  if (!is.data.frame(record$warnings) ||
      !all(names(c2_empty_warnings()) %in% names(record$warnings))) {
    stop("Terminal warning table is malformed.")
  }
  if (!is.list(record$hashes) ||
      !all(c("input", "config", "source") %in% names(record$hashes))) {
    stop("Terminal record lacks input/config/source hashes.")
  }
  invisible(TRUE)
}

c2_write_terminal_bundle <- function(record, fit, directory) {
  if (dir.exists(directory)) {
    marker <- file.path(directory, "TERMINAL_COMPLETE.txt")
    terminal_path <- file.path(directory, "terminal_record.rds")
    if (file.exists(marker) && file.exists(terminal_path)) {
      existing <- readRDS(terminal_path)
      c2_validate_terminal_record(existing)
      marker_text <- readLines(marker, warn = FALSE)
      expected <- paste0("terminal_record_sha256=", c2_sha256(terminal_path))
      if (!expected %in% marker_text) stop("Existing terminal marker hash mismatch.")
      fit_path <- file.path(directory, "fit.rds")
      fit_marker <- grep("^fit_sha256=", marker_text, value = TRUE)
      if (length(fit_marker) != 1L) {
        stop("Existing terminal marker lacks one fit hash.")
      }
      expected_fit <- sub("^fit_sha256=", "", fit_marker[[1L]])
      if (identical(expected_fit, "NA")) {
        if (file.exists(fit_path)) stop("Unexpected fit exists beside NA marker.")
      } else if (!file.exists(fit_path) || c2_sha256(fit_path) != expected_fit) {
        stop("Existing fit hash does not match terminal marker.")
      }
      return(list(reused = TRUE, record = existing,
                  terminal_path = normalizePath(terminal_path, winslash = "/")))
    }
    stop("Incomplete pre-existing fit directory requires manual quarantine: ",
         directory)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  c2_validate_terminal_record(record)
  fit_hash <- NULL
  if (!is.null(fit)) {
    fit_write <- c2_atomic_save_rds(fit, file.path(directory, "fit.rds"))
    fit_hash <- fit_write$sha256
    record$hashes$output_fit <- fit_hash
  } else {
    record$hashes$output_fit <- NA_character_
  }
  terminal_write <- c2_atomic_save_rds(
    record, file.path(directory, "terminal_record.rds")
  )
  c2_atomic_write_csv(
    c2_terminal_summary(record), file.path(directory, "terminal_summary.csv")
  )
  c2_atomic_write_lines(c(
    "V0LV_TERMINAL_RECORD_COMPLETE",
    paste0("terminal_status=", record$terminal_status),
    paste0("terminal_record_sha256=", terminal_write$sha256),
    paste0("fit_sha256=", if (is.null(fit_hash)) "NA" else fit_hash),
    "success_inference_requires_verified_terminal_record=TRUE"
  ), file.path(directory, "TERMINAL_COMPLETE.txt"))
  list(reused = FALSE, record = record,
       terminal_path = normalizePath(
         file.path(directory, "terminal_record.rds"), winslash = "/"
       ))
}

c2_verify_formal_binding <- function(binding_dir, project_root) {
  status_path <- file.path(binding_dir, "BINDING_STATUS.csv")
  source_path <- file.path(binding_dir, "FULL_RUNTIME_SOURCE_BINDINGS.csv")
  if (!file.exists(status_path) || !file.exists(source_path)) {
    stop("Formal source binding is missing status or source index.")
  }
  status <- utils::read.csv(status_path, stringsAsFactors = FALSE)
  if (nrow(status) != 1L || status$status[[1L]] != "frozen" ||
      !isTRUE(status$formal_execution_authorized[[1L]])) {
    stop("Formal execution requires a separately authorized frozen binding.")
  }
  index <- utils::read.csv(source_path, stringsAsFactors = FALSE)
  verification <- c2_verify_hash_index(index, project_root)
  if (!all(verification$matches)) {
    failed <- verification$relative_path[!verification$matches]
    stop("Formal source hash mismatch: ", paste(failed, collapse = ", "))
  }
  environment_path <- file.path(binding_dir, "ENVIRONMENT_BINDING.csv")
  if (!file.exists(environment_path)) {
    stop("Formal environment binding is missing.")
  }
  expected_environment <- utils::read.csv(
    environment_path, stringsAsFactors = FALSE
  )
  actual_environment <- c2_current_environment_binding(
    unique(expected_environment$package[
      expected_environment$binding_type == "package"
    ])
  )
  expected_key <- paste(expected_environment$binding_type,
                        expected_environment$package,
                        expected_environment$key, sep = "\r")
  actual_key <- paste(actual_environment$binding_type,
                      actual_environment$package,
                      actual_environment$key, sep = "\r")
  if (anyDuplicated(expected_key) || anyDuplicated(actual_key) ||
      !setequal(expected_key, actual_key)) {
    stop("Formal environment binding keys do not match.")
  }
  actual_environment <- actual_environment[match(expected_key, actual_key), ]
  if (any(as.character(expected_environment$value) !=
          as.character(actual_environment$value))) {
    failed <- expected_key[
      as.character(expected_environment$value) !=
        as.character(actual_environment$value)
    ]
    stop("Formal R/platform/BLAS/dependency mismatch: ",
         paste(failed, collapse = ", "))
  }
  namespace_path <- file.path(
    binding_dir, "MULTIFSYNC_NAMESPACE_FUNCTION_BINDINGS.csv"
  )
  if (!file.exists(namespace_path)) {
    stop("Formal multiFSYNC namespace function binding is missing.")
  }
  if (!requireNamespace("multiFSYNC", quietly = TRUE)) {
    stop("Formal execution requires the bound multiFSYNC namespace.")
  }
  namespace_index <- utils::read.csv(namespace_path, stringsAsFactors = FALSE)
  namespace <- asNamespace("multiFSYNC")
  current_names <- sort(ls(namespace, all.names = TRUE)[vapply(
    ls(namespace, all.names = TRUE),
    function(name) is.function(get(name, envir = namespace, inherits = FALSE)),
    logical(1)
  )])
  if (!identical(sort(as.character(namespace_index$function_name)),
                 current_names)) {
    stop("Formal multiFSYNC namespace function inventory mismatch.")
  }
  current_hash <- vapply(current_names, function(name) {
    c2_function_signature_sha256(
      get(name, envir = namespace, inherits = FALSE)
    )
  }, character(1))
  expected_hash <- setNames(
    as.character(namespace_index$sha256), namespace_index$function_name
  )[current_names]
  if (any(current_hash != expected_hash)) {
    stop("Formal multiFSYNC namespace function hash mismatch: ",
         paste(current_names[current_hash != expected_hash], collapse = ", "))
  }
  invisible(verification)
}

c2_current_environment_binding <- function(packages = character()) {
  packages <- sort(unique(as.character(packages[nzchar(packages)])))
  info <- sessionInfo()
  base <- data.frame(
    binding_type = "runtime",
    package = "",
    key = c("R_version", "platform", "BLAS", "LAPACK"),
    value = c(
      as.character(getRversion()), R.version$platform,
      as.character(info$BLAS %||% ""), as.character(info$LAPACK %||% "")
    ),
    stringsAsFactors = FALSE
  )
  package_rows <- lapply(packages, function(package) {
    version <- if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else "NOT_INSTALLED"
    data.frame(binding_type = "package", package = package,
               key = "version", value = version,
               stringsAsFactors = FALSE)
  })
  if (length(package_rows)) rbind(base, do.call(rbind, package_rows)) else base
}

c2_function_signature_sha256 <- function(fun) {
  if (!is.function(fun)) stop("fun must be a function.")
  value <- list(
    formals = formals(fun),
    body = paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
  )
  temp <- tempfile(fileext = ".rds")
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  saveRDS(value, temp, version = 3L)
  c2_sha256(temp)
}
