#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4e_common_20260827_v1.R"),
       local = FALSE)

inputs <- s4e_verify_inputs(load_runtime = TRUE)
observations <- lapply(UC_DATA_IDS, uc_load_observation)
worker <- file.path(UC_ROOT, "s4e_worker_20260827_v1.R")
status <- system2("Rscript", c(
  worker, paste0("--fit-id=", inputs$fit_manifest$fit_id[[1L]]),
  "--check-only=true"
))

toy <- data.frame(
  fit_id = sprintf("toy_%02d", 1:8), final_elbo = c(1:7, 9),
  objective_eligible = TRUE, stringsAsFactors = FALSE
)
toy <- toy[toy$objective_eligible & is.finite(toy$final_elbo), ]
toy <- toy[order(-toy$final_elbo, toy$fit_id), ]
checks <- c(
  three_bound_observations = length(observations) == 3L,
  observation_only_class = all(vapply(observations, function(value) {
    identical(value$bundle_class, "v0lv_observation_only_candidate2")
  }, logical(1L))),
  eighteen_registered_new_fits = nrow(inputs$fit_manifest) == 18L,
  six_registered_inherited_endpoints = nrow(inputs$old_manifest) == 6L,
  eight_total_starts_per_data = all(table(c(
    inputs$fit_manifest$data_id, inputs$old_manifest$data_id
  )) == 8L),
  new_seed_indices_3_to_8 = setequal(
    inputs$fit_manifest$seed_index, 3:8
  ),
  single_cpu_per_fit = all(inputs$fit_manifest$n_cpus == 1L),
  fixed400_no_continuation = all(inputs$fit_manifest$fixed_400_endpoint) &&
    !any(inputs$fit_manifest$continuation) &&
    !any(inputs$fit_manifest$automatic_800),
  truth_free_fit_stop_selection =
    !any(inputs$fit_manifest$truth_available_to_fit) &&
    !any(inputs$fit_manifest$truth_available_to_stopping) &&
    !any(inputs$fit_manifest$truth_available_to_selection),
  worker_check_only = identical(as.integer(status), 0L),
  truth_free_elbo_rule = identical(toy$fit_id[[1L]], "toy_08"),
  development_not_formal = all(inputs$fit_manifest$development_only) &&
    !any(inputs$fit_manifest$formal_v0lv_result)
)
result <- data.frame(check_id = names(checks), passed = unname(checks),
                     stringsAsFactors = FALSE)
utils::write.csv(result, file.path(UC_ROOT, "NONFORMAL_SMOKE_STATUS.csv"),
                 row.names = FALSE, quote = TRUE)
uc_assert(all(checks), "Stage-4E nonformal orchestration smoke failed.")
writeLines(c(
  "status=PASS", paste0("completed_utc=", uc_iso_time()),
  "checks=12/12", "full_size_fit_started=FALSE",
  "truth_evaluated=FALSE", "formal_v0lv_result=FALSE"
), file.path(UC_ROOT, "NONFORMAL_SMOKE_COMPLETE.txt"), useBytes = TRUE)
cat("STAGE4E_NONFORMAL_SMOKE_PASS checks=12/12 fits_started=0 truth=0\n")
