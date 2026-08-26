#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4c_common_20260826_v1.R"),
       local = FALSE)

cli <- uc_parse_cli(commandArgs(trailingOnly = TRUE))
check_only <- identical(tolower(cli$`check-only` %||% "false"), "true")
workers <- UC_OUTER_WORKERS

s4c_verify_environment(load_runtime = TRUE, package_role = "development")
manifest <- uc_read_csv(file.path(UC_ROOT, "FIT_MANIFEST.csv"))
s4c_validate_fit_manifest(manifest)
uc_assert(file.exists(file.path(UC_ROOT, "DATA_GENERATION_COMPLETE.txt")),
          "Data generation/sealing is incomplete.")
for (data_id in UC_DATA_IDS) uc_load_observation(data_id)

worker_path <- file.path(UC_ROOT, "s4c_worker_20260826_v1.R")
uc_assert(file.exists(worker_path), "Stage-4C worker script is missing.")
if (check_only) {
  status <- system2("Rscript", c(
    worker_path, paste0("--fit-id=", manifest$fit_id[[1L]]),
    "--check-only=true"
  ))
  uc_assert(as.integer(status) == 0L, "Worker check-only failed.")
  cat("SUPERVISOR_CHECK_ONLY_PASS rows=12 workers=4 fixed_t1=400\n")
  quit(save = "no", status = 0L)
}

log_dir <- file.path(UC_ROOT, "supervisor_logs")
fit_root <- file.path(UC_ROOT, "fits")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
dir.create(fit_root, recursive = TRUE, showWarnings = FALSE, mode = "0700")

write_state <- function(state, detail = "") {
  complete_count <- length(list.files(
    fit_root, pattern = "^FIT_COMPLETE[.]txt$", recursive = TRUE
  ))
  error_count <- length(list.files(
    fit_root, pattern = "^FIT_ERROR[.]txt$", recursive = TRUE
  ))
  value <- data.frame(
    timestamp_utc = uc_iso_time(), state = state, detail = detail,
    total_fits = UC_EXPECTED_FITS, complete_fits = complete_count,
    error_fits = error_count, outer_workers = workers,
    per_fit_n_cpus = 1L, fixed_T1_sweeps = 400L,
    public_g12_stopping_profile = TRUE,
    truth_read_by_worker_or_supervisor = FALSE,
    continuation_started = FALSE, automatic_800_started = FALSE,
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )
  utils::write.csv(value, file.path(UC_ROOT, "CURRENT_STATUS.csv"),
                   row.names = FALSE, quote = TRUE, na = "")
  cat(sprintf("[%s] state=%s complete=%d error=%d %s\n",
              value$timestamp_utc, state, complete_count, error_count, detail),
      file = file.path(UC_ROOT, "SUPERVISOR.log"), append = TRUE)
}

run_one <- function(fit_id) {
  key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
  exit_path <- file.path(log_dir, paste0(key, ".exit.csv"))
  if (file.exists(exit_path)) {
    prior <- uc_read_csv(exit_path)
    uc_assert(prior$exit_status[[1L]] == 0L &&
                file.exists(file.path(fit_root, fit_id, "FIT_COMPLETE.txt")),
              paste0("Prior non-success action record exists: ", fit_id))
    return(list(fit_id = fit_id, status = 0L,
                elapsed_seconds = prior$elapsed_seconds[[1L]], reused = TRUE))
  }
  started <- Sys.time()
  status <- system2(
    "Rscript", c(worker_path, paste0("--fit-id=", fit_id)),
    stdout = file.path(log_dir, paste0(key, ".stdout.log")),
    stderr = file.path(log_dir, paste0(key, ".stderr.log"))
  )
  ended <- Sys.time()
  record <- data.frame(
    fit_id = fit_id, started_utc = uc_iso_time(started),
    ended_utc = uc_iso_time(ended),
    elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
    exit_status = as.integer(status), stringsAsFactors = FALSE
  )
  utils::write.csv(record, exit_path, row.names = FALSE, quote = TRUE, na = "")
  list(fit_id = fit_id, status = as.integer(status),
       elapsed_seconds = record$elapsed_seconds[[1L]], reused = FALSE)
}

write_state("RUNNING", "12 registered independent fixed-400 G12 fits")
results <- parallel::mclapply(
  manifest$fit_id, run_one, mc.cores = workers,
  mc.preschedule = FALSE, mc.set.seed = FALSE
)
statuses <- vapply(results, function(value) value$status, integer(1L))
if (any(statuses != 0L)) {
  failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                   character(1L))
  write_state("ERROR", paste0(
    "retained_without_automatic_rerun=", paste(failed, collapse = ";")
  ))
  uc_abort("Fit errors retained without selective rerun: ",
           paste(failed, collapse = ", "))
}

terminal_paths <- file.path(fit_root, manifest$fit_id, "terminal_record.rds")
complete_paths <- file.path(fit_root, manifest$fit_id, "FIT_COMPLETE.txt")
uc_assert(all(file.exists(terminal_paths)) && all(file.exists(complete_paths)),
          "Successful worker exits lack complete fit bundles.")
records <- lapply(terminal_paths, readRDS)
uc_assert(all(vapply(records, function(record) {
  identical(record$terminal_status, "fixed_400_complete") &&
    isTRUE(record$objective_eligible) && record$actual_T1_sweeps == 400L &&
    isTRUE(record$fixed_400_endpoint) &&
    isTRUE(record$public_g12_stopping_profile_applied) &&
    !isTRUE(record$truth_used) && !isTRUE(record$continuation_used) &&
    !isTRUE(record$automatic_800_used) && !isTRUE(record$formal_v0lv_result)
}, logical(1L))),
"A completed endpoint violates fixed-horizon, eligibility, or isolation.")

terminals <- do.call(rbind, lapply(records, s4c_terminal_summary))
terminals <- terminals[match(manifest$fit_id, terminals$fit_id), , drop = FALSE]
uc_assert(nrow(terminals) == UC_EXPECTED_FITS &&
            !anyDuplicated(terminals$fit_id) &&
            identical(terminals$fit_id, manifest$fit_id),
          "Terminal records do not exactly join to the manifest.")
utils::write.csv(terminals, file.path(UC_ROOT, "ALL_12_TERMINALS.csv"),
                 row.names = FALSE, quote = TRUE, na = "")
writeLines(c(
  "status=FITS_COMPLETE_FIXED_400",
  paste0("completed_utc=", uc_iso_time()),
  "fits=12", "objective_eligible=12", "outer_workers=4",
  "per_fit_n_cpus=1", "ordinary_T1_sweeps=400",
  "public_g12_stopping_profile=TRUE",
  "truth_read_by_fit_stopping_or_supervisor=FALSE",
  "truth_unseal_authorized=FALSE", "continuation_started=FALSE",
  "automatic_800_started=FALSE", "formal_v0lv_result=FALSE"
), file.path(UC_ROOT, "FITS_COMPLETE.txt"), useBytes = TRUE)
write_state("FITS_COMPLETE_WAITING_FOR_TRUTH_FREE_ANALYSIS",
            "12/12 objective-eligible fixed-400 fits")
cat("FITS_COMPLETE fixed_400=12/12 objective_eligible=12 truth=0\n")
