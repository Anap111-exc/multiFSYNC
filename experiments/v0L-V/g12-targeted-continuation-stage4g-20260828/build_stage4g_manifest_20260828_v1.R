#!/usr/bin/env Rscript

options(warn = 1, stringsAsFactors = FALSE)
argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
root <- dirname(normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
))
stage4e_root <- Sys.getenv(
  "STAGE4G_STAGE4E_ROOT",
  unset = "/root/v0lv-g12-focused-reachability-stage4e-20260827-v1"
)
output <- file.path(root, "CONTINUATION_MANIFEST.csv")
if (file.exists(output)) {
  stop("Refusing to overwrite CONTINUATION_MANIFEST.csv.", call. = FALSE)
}
assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
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

inventory <- read_csv(file.path(
  stage4e_root, "ALL_24_TRUTH_FREE_ENDPOINTS.csv"
))
reachability <- read_csv(file.path(
  stage4e_root, "analysis_20260827_v1", "PER_ENDPOINT_REACHABILITY.csv"
))
winners <- read_csv(file.path(
  stage4e_root, "TRUTH_FREE_ELBO_WINNERS_8_STARTS.csv"
))
target_ids <- c(
  "g12ss4c_base_01__G12__04",
  "g12ss4c_sparse_02__G12__07"
)
selection_role <- c(
  "exact_baseline_marginal_short_rss_control",
  "exact_sparse_winner_late_scale_movement"
)
target <- inventory[match(target_ids, inventory$fit_id), , drop = FALSE]
truth_rows <- reachability[match(target_ids, reachability$fit_id),
                           , drop = FALSE]
assert(!anyNA(target$fit_id) && !anyNA(truth_rows$fit_id) &&
         identical(as.character(target$fit_id), target_ids) &&
         identical(as.character(truth_rows$fit_id), target_ids),
       "Target identities do not join exactly.")
assert(all(target$terminal_status == "fixed_400_complete") &&
         all(target$actual_T1_sweeps == 400L) &&
         all(target$objective_eligible) && all(target$slow_case) &&
         !any(target$truth_used) && !any(target$continuation_used) &&
         all(truth_rows$structure_exact),
       "Targets are not the registered exact slow fixed-400 endpoints.")
assert(
  identical(winners$selected_fit_id[
    winners$data_id == "g12ss4c_sparse_02"
  ], target_ids[[2L]]) &&
    identical(winners$runner_up_fit_id[
      winners$data_id == "g12ss4c_base_01"
    ], target_ids[[1L]]),
  "Target winner/control roles differ from frozen Stage-4E selection."
)

source_relative <- file.path("fits", target$fit_id, "fit.rds")
observation_relative <- file.path(
  "data", target$data_id, "observation_bundle.rds"
)
source_paths <- file.path(stage4e_root, source_relative)
observation_paths <- file.path(stage4e_root, observation_relative)
assert(all(file.exists(source_paths)) && all(file.exists(observation_paths)),
       "A registered source fit or observation bundle is missing.")

manifest <- data.frame(
  experiment_id = "G12_TARGETED_CONTINUATION_STAGE4G_V1_20260828",
  continuation_id = c(
    "g12s4g_base_control__G12__04__to800",
    "g12s4g_sparse02_winner__G12__07__to800"
  ),
  source_fit_id = target$fit_id,
  data_id = target$data_id,
  scenario_id = target$scenario_id,
  seed_index = as.integer(target$seed_index),
  fit_seed = as.integer(target$fit_seed),
  selection_role = selection_role,
  adaptive_post_truth_target = TRUE,
  source_endpoint_role = c("runner_up", "winner"),
  source_fit_relative_path = source_relative,
  source_fit_sha256 = vapply(source_paths, sha256, character(1L)),
  observation_relative_path = observation_relative,
  observation_sha256 = vapply(
    observation_paths, sha256, character(1L)
  ),
  source_final_elbo = as.numeric(target$final_elbo),
  source_objective_eligible = TRUE,
  source_structure_exact = TRUE,
  source_slow_case = TRUE,
  source_t1_sweeps = 400L,
  additional_t1_sweeps = 400L,
  cumulative_endpoint = 800L,
  cumulative_checkpoints = paste(
    c(400L, seq.int(420L, 800L, 20L)), collapse = ";"
  ),
  n_cpus = 1L, outer_parallel_fit_limit = 2L,
  truth_available_to_continuation = FALSE,
  winner_selection_reopened = FALSE,
  automatic_further_extension = FALSE,
  development_only = TRUE, formal_v0lv_result = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest, output, row.names = FALSE, quote = TRUE, na = ""
)
writeLines(c(
  "status=REGISTERED", "continuations=2",
  "source_cumulative_T1=400", "target_cumulative_T1=800",
  "adaptive_post_truth_target=TRUE",
  "truth_available_to_continuation=FALSE",
  "winner_selection_reopened=FALSE", "formal_v0lv_result=FALSE"
), file.path(root, "MANIFEST_REGISTERED.txt"), useBytes = TRUE)
cat("STAGE4G_MANIFEST_REGISTERED continuations=2 target=800\n")
