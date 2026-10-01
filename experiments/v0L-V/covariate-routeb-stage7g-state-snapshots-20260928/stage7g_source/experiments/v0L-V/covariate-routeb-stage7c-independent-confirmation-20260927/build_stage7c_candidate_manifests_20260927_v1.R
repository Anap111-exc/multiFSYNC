#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
protocol_dir <- dirname(file_arg)
if (!dir.exists(protocol_dir)) stop("Protocol directory does not exist.")

data_path <- file.path(protocol_dir, "DATA_MANIFEST_CANDIDATE_V1.csv")
fit_path <- file.path(protocol_dir, "FIT_MANIFEST_CANDIDATE_V1.csv")
status_path <- file.path(protocol_dir, "DRAFT_STATUS_V1_20260927.txt")
if (any(file.exists(c(data_path, fit_path, status_path)))) {
  stop("Versioned Stage 7C candidate outputs already exist; refusing overwrite.")
}

protocol_id <- "COVARIATE_ROUTEB_STAGE7C_V1_20260927"
logical_data_id <- sprintf("covrb7c_%02d", 1:3)
seed_table <- data.frame(
  logical_data_id = logical_data_id,
  data_index = 1:3,
  data_seed = 99711001:99711003,
  truth_basis_seed = 99712001:99712003,
  truth_loading_seed = 99713001:99713003,
  mean_truth_seed = 99714001:99714003,
  beta_truth_seed = 99715001:99715003,
  stringsAsFactors = FALSE)

arms <- data.frame(
  beta_arm = c("C0", "C1", "C2"),
  beta_truth = c("zero", "nonzero_bspline", "reuse_C1_nonzero"),
  d_generated = c(2L, 2L, 2L),
  d_fit = c(2L, 2L, 0L),
  stringsAsFactors = FALSE)

data_rows <- vector("list", nrow(seed_table) * nrow(arms))
k <- 0L
for (i in seq_len(nrow(seed_table))) {
  for (a in seq_len(nrow(arms))) {
    k <- k + 1L
    logical_id <- seed_table$logical_data_id[[i]]
    arm <- arms$beta_arm[[a]]
    data_id <- paste(logical_id, arm, sep = "__")
    source_id <- if (arm == "C2") {
      paste(logical_id, "C1", sep = "__")
    } else data_id
    data_rows[[k]] <- data.frame(
      protocol_id = protocol_id,
      status = "draft_not_authorized",
      logical_data_id = logical_id,
      data_index = seed_table$data_index[[i]],
      beta_arm = arm,
      data_id = data_id,
      observation_source_data_id = source_id,
      data_seed = seed_table$data_seed[[i]],
      truth_basis_seed = seed_table$truth_basis_seed[[i]],
      truth_loading_seed = seed_table$truth_loading_seed[[i]],
      mean_truth_seed = seed_table$mean_truth_seed[[i]],
      beta_truth_seed = seed_table$beta_truth_seed[[i]],
      beta_truth = arms$beta_truth[[a]],
      d_generated = arms$d_generated[[a]],
      d_fit = arms$d_fit[[a]],
      new_unseen_at_registration = TRUE,
      truth_available_to_fit = FALSE,
      truth_available_to_stopping = FALSE,
      truth_available_to_selection = FALSE,
      protocol_confirmed = FALSE,
      seed_schedule_confirmed = FALSE,
      execution_authorized = FALSE,
      data_generated = FALSE,
      formal_fit_started = FALSE,
      stringsAsFactors = FALSE)
  }
}
data_manifest <- do.call(rbind, data_rows)

fit_rows <- vector("list", nrow(data_manifest) * 12L)
k <- 0L
for (i in seq_len(nrow(data_manifest))) {
  row <- data_manifest[i, , drop = FALSE]
  for (start_index in 1:12) {
    k <- k + 1L
    fit_seed <- 99716000L + row$data_index * 100L + start_index
    fit_rows[[k]] <- data.frame(
      protocol_id = protocol_id,
      status = "draft_not_authorized",
      fit_id = sprintf("%s__%s__start%02d",
                       row$logical_data_id, row$beta_arm, start_index),
      logical_data_id = row$logical_data_id,
      data_index = row$data_index,
      beta_arm = row$beta_arm,
      data_id = row$data_id,
      observation_source_data_id = row$observation_source_data_id,
      data_seed = row$data_seed,
      truth_basis_seed = row$truth_basis_seed,
      truth_loading_seed = row$truth_loading_seed,
      mean_truth_seed = row$mean_truth_seed,
      beta_truth_seed = row$beta_truth_seed,
      d_generated = row$d_generated,
      d_fit = row$d_fit,
      fit_K = 9L,
      n_g = 51L,
      lambda_orth = 0,
      start_index = start_index,
      fit_seed = fit_seed,
      initialization = "random",
      function_initialization = "gram_unit_energy",
      pre_score_sweeps = 1L,
      anneal = "1;1.9;100",
      planned_annealing_sweeps = 99L,
      planned_t1_sweeps = 400L,
      maxit = 499L,
      n_cpus = 1L,
      objective_eligibility_required = TRUE,
      truth_available_to_fit = FALSE,
      truth_available_to_stopping = FALSE,
      truth_available_to_selection = FALSE,
      protocol_confirmed = FALSE,
      seed_schedule_confirmed = FALSE,
      execution_authorized = FALSE,
      data_generated = FALSE,
      fit_started = FALSE,
      stringsAsFactors = FALSE)
  }
}
fit_manifest <- do.call(rbind, fit_rows)

stopifnot(
  nrow(data_manifest) == 9L,
  nrow(fit_manifest) == 108L,
  length(unique(data_manifest$data_id)) == 9L,
  length(unique(fit_manifest$fit_id)) == 108L,
  all(table(data_manifest$logical_data_id) == 3L),
  all(table(fit_manifest$logical_data_id, fit_manifest$beta_arm) == 12L),
  all(!data_manifest$execution_authorized),
  all(!data_manifest$data_generated),
  all(!fit_manifest$execution_authorized),
  all(!fit_manifest$fit_started),
  all(!fit_manifest$truth_available_to_fit),
  all(!fit_manifest$truth_available_to_stopping),
  all(!fit_manifest$truth_available_to_selection))

for (logical_id in logical_data_id) {
  group <- fit_manifest[fit_manifest$logical_data_id == logical_id, ]
  seed_by_arm <- split(group$fit_seed, group$beta_arm)
  stopifnot(identical(seed_by_arm$C0, seed_by_arm$C1),
            identical(seed_by_arm$C1, seed_by_arm$C2))
  dgroup <- data_manifest[data_manifest$logical_data_id == logical_id, ]
  stopifnot(
    dgroup$observation_source_data_id[dgroup$beta_arm == "C2"] ==
      dgroup$data_id[dgroup$beta_arm == "C1"])
}

write.csv(data_manifest, data_path, row.names = FALSE, na = "")
write.csv(fit_manifest, fit_path, row.names = FALSE, na = "")
writeLines(c(
  "status=DRAFT_NOT_AUTHORIZED",
  paste0("protocol_id=", protocol_id),
  "logical_data_seeds=3",
  "data_arm_rows=9",
  "planned_fits=108",
  "starts_per_data_arm=12",
  "seed_collision_search=NO_EXACT_COLLISIONS_FOUND_BEFORE_DRAFT",
  "protocol_confirmed=FALSE",
  "seed_schedule_confirmed=FALSE",
  "execution_authorized=FALSE",
  "stage7c_data_generated=FALSE",
  "stage7c_fits_started=FALSE",
  "truth_unseal_authorized=FALSE",
  "scientific_evaluation_started=FALSE",
  "ecs_started_by_this_draft=FALSE"), status_path, useBytes = TRUE)

cat("Stage 7C candidate data rows:", nrow(data_manifest), "\n")
cat("Stage 7C candidate fit rows:", nrow(fit_manifest), "\n")
cat("No data or fit was generated.\n")
