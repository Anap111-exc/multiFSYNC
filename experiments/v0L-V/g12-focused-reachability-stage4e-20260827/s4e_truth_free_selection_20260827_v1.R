#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
source(file.path(dirname(script_path), "s4e_common_20260827_v1.R"),
       local = FALSE)

output_inventory <- file.path(UC_ROOT, "ALL_24_TRUTH_FREE_ENDPOINTS.csv")
output_winners <- file.path(UC_ROOT, "TRUTH_FREE_ELBO_WINNERS_8_STARTS.csv")
output_marker <- file.path(UC_ROOT, "TRUTH_FREE_SELECTION_COMPLETE.txt")
uc_assert(!any(file.exists(c(output_inventory, output_winners, output_marker))),
          "Refusing to overwrite a frozen Stage-4E selection output.")

inputs <- s4e_verify_inputs(load_runtime = TRUE)
uc_assert(file.exists(file.path(UC_ROOT, "FITS_COMPLETE.txt")),
          "The 18 new Stage-4E fits are incomplete.")

new_paths <- file.path(
  UC_ROOT, "fits", inputs$fit_manifest$fit_id, "terminal_record.rds"
)
uc_assert(all(file.exists(new_paths)), "A new terminal record is missing.")
new_records <- lapply(new_paths, readRDS)
old_records <- lapply(inputs$old_terminal_paths, readRDS)

old <- do.call(rbind, lapply(old_records, s4e_terminal_summary,
                            source_stage = "stage4c_inherited"))
new <- do.call(rbind, lapply(new_records, s4e_terminal_summary,
                            source_stage = "stage4e_new"))
inventory <- rbind(old, new)
inventory <- inventory[order(
  match(inventory$data_id, UC_DATA_IDS), inventory$seed_index
), , drop = FALSE]
rownames(inventory) <- NULL

uc_assert(
  nrow(inventory) == 24L && !anyDuplicated(inventory$fit_id) &&
    all(table(inventory$data_id) == 8L) &&
    all(inventory$terminal_status == "fixed_400_complete") &&
    all(inventory$objective_eligible) &&
    all(is.finite(inventory$final_elbo)) &&
    !any(inventory$truth_used) && !any(inventory$continuation_used) &&
    !any(inventory$automatic_800_used) &&
    !any(inventory$formal_v0lv_result),
  "The combined 24-endpoint truth-free inventory is invalid."
)

winner_rows <- lapply(UC_DATA_IDS, function(data_id) {
  group <- inventory[inventory$data_id == data_id, , drop = FALSE]
  group <- group[order(-group$final_elbo, group$fit_id), , drop = FALSE]
  uc_assert(nrow(group) == 8L, "Every target data set must have eight starts.")
  data.frame(
    experiment_id = UC_EXPERIMENT_ID,
    data_id = data_id, scenario_id = group$scenario_id[[1L]],
    eligible_endpoint_count = nrow(group),
    selected_fit_id = group$fit_id[[1L]],
    selected_seed_index = group$seed_index[[1L]],
    selected_fit_seed = group$fit_seed[[1L]],
    selected_source_stage = group$source_stage[[1L]],
    selected_final_elbo = group$final_elbo[[1L]],
    runner_up_fit_id = group$fit_id[[2L]],
    runner_up_final_elbo = group$final_elbo[[2L]],
    elbo_gap = group$final_elbo[[1L]] - group$final_elbo[[2L]],
    selection_rule = "max_finite_objective_eligible_ordinary_T1_ELBO",
    deterministic_tie_break = "fit_id_ascending",
    selection_used_truth = FALSE, selection_used_structure = FALSE,
    selection_used_NRMSE_ISE_RPL = FALSE,
    adaptive_post_truth_design = TRUE, development_only = TRUE,
    formal_v0lv_result = FALSE, stringsAsFactors = FALSE
  )
})
winners <- do.call(rbind, winner_rows)
rownames(winners) <- NULL

utils::write.csv(inventory, output_inventory, row.names = FALSE,
                 quote = TRUE, na = "")
utils::write.csv(winners, output_winners, row.names = FALSE,
                 quote = TRUE, na = "")
writeLines(c(
  "status=TRUTH_FREE_SELECTION_FROZEN",
  paste0("completed_utc=", uc_iso_time()),
  "target_data_sets=3", "eligible_endpoints=24",
  "starts_per_data=8", "selection=max_valid_ordinary_T1_ELBO",
  "selection_used_truth=FALSE", "selection_used_structure=FALSE",
  "selection_used_NRMSE_ISE_RPL=FALSE",
  "truth_evaluation_started=FALSE", "formal_v0lv_result=FALSE"
), output_marker, useBytes = TRUE)
cat("STAGE4E_TRUTH_FREE_SELECTION_PASS data=3 endpoints=24 starts=8\n")
