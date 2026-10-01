#!/usr/bin/env Rscript

# Developmental role-rotation audit. Reuses the frozen Stage 7A fit functions
# without reimplementing their optimization or objective-eligibility rules.
options(warn = 1, stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"),
      silent = TRUE)
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")

`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
assert(length(script_arg) == 1L, "Run this file through Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
script_root <- dirname(script_path)
stage7a_path <- normalizePath(file.path(
  script_root, "..", "covariate-effect-stage7a-20260919",
  "run_covariate_stage7a_20260919_v1.R"),
  winslash = "/", mustWork = TRUE)
stage7a_sha256 <- "7b99b9030433d7033eec74cd3fa12736a4952b22a315d7b358ce881c9d5e3732"
generator_sha256 <- "2c1fb6ee1281e243b7a67752287f1f057dbb1d2e4239e6ac26369823371ee3c5"
assert(requireNamespace("digest", quietly = TRUE), "digest is required.")
sha256 <- function(path) tolower(digest::digest(
  file = path, algo = "sha256", serialize = FALSE))
assert(identical(sha256(stage7a_path), stage7a_sha256),
       "Stage 7A runner source SHA256 differs from the registered source.")
generator_path <- normalizePath(file.path(
  dirname(stage7a_path), "..", "..", "..", "R",
  "generate_data_structured.R"), winslash = "/", mustWork = TRUE)
assert(identical(sha256(generator_path), generator_sha256),
       "Stage 7A generator source SHA256 differs from the registered source.")

# Evaluate declarations only; never run the old runner's command dispatch.
stage7a <- new.env(parent = globalenv())
expressions <- parse(stage7a_path)
cli_index <- which(vapply(expressions, function(x) {
  is.call(x) && length(x) >= 3L && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name("cli"))
}, logical(1L)))
assert(length(cli_index) == 1L && cli_index > 1L,
       "Cannot locate the Stage 7A dispatch boundary.")
for (i in seq_len(cli_index - 1L)) eval(expressions[[i]], envir = stage7a)
stage7a$script_root <- dirname(stage7a_path)
stage7a$experiment_id <- "MULTIFSYNC_STAGE7B_ROLE_ROTATION_V1_20260921"
experiment_id <- stage7a$experiment_id

role_maps <- list(rot231 = c(2L, 3L, 1L),
                  rot312 = c(3L, 1L, 2L))
source_data_row <- data.frame(
  data_seed = 97201003L, truth_basis_seed = 97202003L,
  truth_loading_seed = 97203003L, beta_bspline_seed = 97204003L)

data_manifest <- function() {
  data.frame(
    experiment_id = experiment_id,
    data_id = paste0("s7b_03__", names(role_maps)),
    logical_data_id = "s7b_03", beta_arm = "sine", data_index = 3L,
    data_seed = source_data_row$data_seed[[1L]],
    truth_basis_seed = source_data_row$truth_basis_seed[[1L]],
    truth_loading_seed = source_data_row$truth_loading_seed[[1L]],
    shared_block = vapply(role_maps, `[[`, integer(1L), 1L),
    study1_block = vapply(role_maps, `[[`, integer(1L), 2L),
    study2_block = vapply(role_maps, `[[`, integer(1L), 3L),
    truth_available_to_fit = FALSE,
    truth_available_to_selection = FALSE,
    truth_unseal_authorized = FALSE, row.names = NULL)
}

fit_manifest <- function(data_rows, smoke = FALSE) {
  starts <- if (smoke) 1L else seq_len(12L)
  rows <- lapply(seq_len(nrow(data_rows)), function(i) {
    row <- data_rows[i, , drop = FALSE]
    data.frame(
      fit_id = sprintf("%s_start%02d", row$data_id[[1L]], starts),
      data_id = row$data_id[[1L]],
      logical_data_id = row$logical_data_id[[1L]],
      beta_arm = row$beta_arm[[1L]], start_index = starts,
      fit_seed = 97210300L + starts, n_cpus = 1L,
      truth_available_to_fit = FALSE,
      truth_available_to_selection = FALSE)
  })
  do.call(rbind, rows)
}

flatten_numeric <- function(x) unlist(x, recursive = TRUE, use.names = FALSE)
max_abs_difference <- function(x, y) {
  a <- flatten_numeric(x); b <- flatten_numeric(y)
  assert(length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)),
         "Paired numeric fields differ in shape or finiteness.")
  max(abs(a - b))
}

rotate_truth <- function(data, mapping, blocks) {
  assert(identical(sort(as.integer(mapping)), 1:3),
         "The role map must permute the three truth blocks.")
  truth <- data$true_params
  assert(truth$S == 2L && truth$L_f == 1L &&
           identical(as.integer(truth$L_s_by_study), c(1L, 1L)) &&
           identical(as.integer(truth$M_f), 2L),
         "Unexpected Stage 7A dimensions in rotation input.")
  metadata_keys <- c("block_id", "components", "gram", "grid_gram",
                     "method", "degree", "knots", "periodic",
                     "selected_indices", "common_selection_pool_size",
                     "source_binding")
  block_meta <- function(block) block[intersect(metadata_keys, names(block))]
  eval_block <- getFromNamespace("eval_jaoua_splinet_block", "multiFSYNC")

  shared <- blocks[[mapping[[1L]]]]
  truth$phi_dense_list[[1L]] <- shared$phi
  truth$theta_dense_list[[1L]] <- sweep(
    shared$phi, 2L, sqrt(truth$lambda_phi_true[[1L]]), `*`)
  truth$phi_basis_meta[[1L]] <- block_meta(shared)
  for (s in 1:2) {
    specific <- blocks[[mapping[[s + 1L]]]]
    truth$psi_dense_list[[s]][[1L]] <- specific$phi
    truth$kappa_dense_list[[s]][[1L]] <- sweep(
      specific$phi, 2L, sqrt(truth$lambda_psi_true[[s]][[1L]]), `*`)
    truth$psi_basis_meta[[s]][[1L]] <- block_meta(specific)
    for (i in seq_len(truth$n_s[[s]])) {
      t_i <- data$time_obs[[s]][[i]]
      theta <- sweep(eval_block(t_i, shared$evaluation), 2L,
                     sqrt(truth$lambda_phi_true[[1L]]), `*`)
      kappa <- sweep(eval_block(t_i, specific$evaluation), 2L,
                     sqrt(truth$lambda_psi_true[[s]][[1L]]), `*`)
      f <- matrix(as.vector(theta %*% truth$zeta_true[[s]][[1L]][i, ]),
                  ncol = 1L)
      g <- matrix(as.vector(kappa %*% truth$xi_true[[s]][[1L]][i, ]),
                  ncol = 1L)
      truth$f_true_values[[s]][[i]] <- f
      truth$g_true_values[[s]][[i]] <- g
      signal <- truth$mu_true_values[[s]][[i]] +
        truth$beta_true_values[[s]][[i]] +
        f %*% t(truth$a_true) + g %*% t(truth$b_true[[s]])
      truth$signal_true_values[[s]][[i]] <- signal
      noise <- truth$noise_true_values[[s]][[i]]
      data$Y[[s]][[i]] <- lapply(seq_len(truth$p), function(j) {
        signal[, j] + noise[, j]
      })
    }
  }
  truth$truth_role_to_block <- as.integer(mapping)
  data$true_params <- truth
  data
}

pairing_checks <- function(original, rotated, mapping, label) {
  a <- original$true_params; b <- rotated$true_params
  equal_fields <- c(
    time_obs = identical(original$time_obs, rotated$time_obs),
    Z = identical(original$Z, rotated$Z),
    mu = identical(a$mu_true_values, b$mu_true_values),
    beta = identical(a$beta_true_values, b$beta_true_values),
    noise = identical(a$noise_true_values, b$noise_true_values),
    shared_scores = identical(a$zeta_true, b$zeta_true),
    specific_scores = identical(a$xi_true, b$xi_true),
    shared_loadings = identical(a$a_true, b$a_true),
    specific_loadings = identical(a$b_true, b$b_true),
    beta_grid = identical(a$beta_true_grid, b$beta_true_grid),
    block_map = identical(b$truth_role_to_block, as.integer(mapping)))
  assert(all(equal_fields), paste0(label, " breaks pairing: ",
                                  paste(names(equal_fields)[!equal_fields],
                                        collapse = ";")))
  reconstruction_error <- 0
  for (s in 1:2) for (i in seq_len(b$n_s[[s]])) {
    y <- do.call(cbind, rotated$Y[[s]][[i]])
    reconstruction_error <- max(reconstruction_error, max(abs(
      y - b$signal_true_values[[s]][[i]] -
        b$noise_true_values[[s]][[i]])))
  }
  assert(reconstruction_error < 1e-10,
         paste0(label, " response decomposition failed."))
  data.frame(map = label, preserved_fields = length(equal_fields),
             max_response_reconstruction_error = reconstruction_error,
             max_Y_difference_from_baseline =
               max_abs_difference(original$Y, rotated$Y))
}

save_dataset <- function(data, row, output_root, smoke) {
  data_root <- file.path(output_root, "data", row$data_id[[1L]])
  assert(!dir.exists(data_root), "Refusing to overwrite a data directory.")
  dir.create(file.path(data_root, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  observation <- list(
    bundle_class = "stage7a_observation_only", experiment_id = experiment_id,
    data_id = row$data_id[[1L]], Y = data$Y, Z = data$Z,
    time_obs = data$time_obs,
    dimensions = stage7a$generator_config(smoke)[c(
      "S", "n_s", "p", "d", "L_f", "L_s", "M_f", "M_s", "K")],
    truth_available_to_fit = FALSE)
  sealed <- list(
    bundle_class = "stage7b_sealed_truth", experiment_id = experiment_id,
    data_id = row$data_id[[1L]], role_to_block = as.integer(c(
      row$shared_block, row$study1_block, row$study2_block)),
    truth_unseal_authorized = FALSE, data = data)
  stage7a$save_new(observation, file.path(data_root, "observation_bundle.rds"))
  stage7a$save_new(sealed, file.path(data_root, "sealed_truth", "truth_bundle.rds"))
  stage7a$lines_new(c("status=DATA_SEAL_COMPLETE",
                      "truth_unseal_authorized=FALSE",
                      "scientific_evaluation_started=FALSE"),
                    file.path(data_root, "DATA_SEAL_COMPLETE.txt"))
}

generate <- function(output_root, smoke, baseline_root) {
  stage7a$load_package()
  assert(!dir.exists(output_root), "Refusing to overwrite output root.")
  rows <- data_manifest()
  fits <- fit_manifest(rows, smoke)
  assert(nrow(rows) == 2L && nrow(fits) == if (smoke) 2L else 24L,
         "Unexpected Stage 7B manifest size.")
  original <- stage7a$paired_beta_data(source_data_row, smoke)$sine
  n_grid <- stage7a$generator_config(smoke)$n_dense
  set.seed(source_data_row$truth_basis_seed[[1L]])
  block_maker <- getFromNamespace("construct_jaoua_splinet_blocks", "multiFSYNC")
  blocks <- block_maker(3L, 2L, n_g = n_grid)
  identity <- rotate_truth(original, 1:3, blocks)
  identity_diff <- max_abs_difference(original$Y, identity$Y)
  assert(identity_diff < 1e-10,
         "Identity role mapping does not reproduce generated Stage 7A Y.")
  if (!smoke) {
    baseline_path <- file.path(baseline_root, "data", "s7a_03__sine",
                               "sealed_truth", "truth_bundle.rds")
    assert(file.exists(baseline_path), "Stage 7A baseline truth is required for QC.")
    prior <- readRDS(baseline_path)$data
    baseline_diff <- max_abs_difference(prior$Y, original$Y)
    assert(baseline_diff < 1e-10 &&
             identical(prior$time_obs, original$time_obs) &&
             identical(prior$Z, original$Z) &&
             identical(prior$true_params$noise_true_values,
                       original$true_params$noise_true_values),
           "The registered Stage 7A seed-03 baseline was not reproduced.")
  } else baseline_diff <- NA_real_
  rotated <- lapply(seq_along(role_maps), function(i) {
    rotate_truth(original, role_maps[[i]], blocks)
  })
  checks <- do.call(rbind, lapply(seq_along(rotated), function(i) {
    pairing_checks(original, rotated[[i]], role_maps[[i]],
                   names(role_maps)[[i]])
  }))
  assert(all(checks$max_Y_difference_from_baseline > 0.01),
         "A role rotation did not change the response.")
  dir.create(output_root, recursive = TRUE, mode = "0700")
  stage7a$csv_new(rows, file.path(output_root, "DATA_MANIFEST.csv"))
  stage7a$csv_new(fits, file.path(output_root, "FIT_MANIFEST.csv"))
  stage7a$csv_new(checks, file.path(output_root, "PAIRING_QC.csv"))
  stage7a$csv_new(data.frame(
    original_stage7a_Y_max_difference = baseline_diff,
    identity_role_map_Y_max_difference = identity_diff,
    stage7a_runner_sha256 = stage7a_sha256,
    generator_sha256 = generator_sha256),
    file.path(output_root, "BASELINE_REPRODUCTION_QC.csv"))
  for (i in seq_len(nrow(rows))) {
    save_dataset(rotated[[i]], rows[i, , drop = FALSE], output_root, smoke)
  }
  stage7a$lines_new(c("status=DATA_GENERATION_AND_SEAL_COMPLETE",
                      paste0("planned_fits=", nrow(fits)),
                      "baseline_reproduced=TRUE",
                      "truth_unseal_authorized=FALSE",
                      "formal_execution=FALSE"),
                    file.path(output_root, "DATA_GENERATION_COMPLETE.txt"))
  cat("STAGE7B_GENERATION_COMPLETE datasets=2 fits=", nrow(fits), "\n", sep = "")
}

batch <- function(output_root, workers) {
  assert(.Platform$OS.type == "unix", "Batch requires Linux.")
  assert(is.finite(workers) && workers >= 1L && workers <= 8L,
         "workers must be 1 through 8.")
  assert(file.exists(file.path(output_root, "DATA_GENERATION_COMPLETE.txt")),
         "Generate and seal before fitting.")
  manifest <- stage7a$read_manifest(file.path(output_root, "FIT_MANIFEST.csv"))
  assert(nrow(manifest) == 24L && !anyDuplicated(manifest$fit_id),
         "The registered batch requires 24 unique fits.")
  log_root <- file.path(output_root, "supervisor_logs")
  assert(!dir.exists(log_root), "Refusing to rerun the batch.")
  dir.create(log_root, recursive = TRUE, mode = "0700")
  run_one <- function(fit_id) {
    code <- system2("Rscript", c(script_path, "--action=fit",
                                 paste0("--fit-id=", fit_id),
                                 paste0("--output-root=", output_root)),
                    stdout = file.path(log_root, paste0(fit_id, ".stdout.log")),
                    stderr = file.path(log_root, paste0(fit_id, ".stderr.log")))
    data.frame(fit_id = fit_id, exit_status = as.integer(code))
  }
  results <- parallel::mclapply(manifest$fit_id, run_one,
                                mc.cores = workers, mc.preschedule = FALSE,
                                mc.set.seed = FALSE)
  stage7a$csv_new(do.call(rbind, results),
                  file.path(output_root, "BATCH_PROCESS_STATUS.csv"))
  complete <- length(list.files(file.path(output_root, "fits"),
                                "^FIT_COMPLETE[.]txt$", recursive = TRUE))
  errors <- length(list.files(file.path(output_root, "fits"),
                              "^FIT_ERROR[.]txt$", recursive = TRUE))
  assert(complete + errors == 24L,
         "Some registered fits have no terminal record; selection withheld.")
  stage7a$select_winners(output_root, FALSE)
  stage7a$lines_new(c("status=TRUTH_FREE_ENDPOINTS_AND_SELECTION_COMPLETE",
                      paste0("complete=", complete),
                      paste0("errors=", errors),
                      paste0("workers=", workers),
                      "per_fit_n_cpus=1", "truth_unseal_authorized=FALSE",
                      "scientific_evaluation_started=FALSE"),
                    file.path(output_root, "FITS_AND_SELECTION_COMPLETE.txt"))
  cat("STAGE7B_BATCH_COMPLETE complete=", complete,
      " errors=", errors, "\n", sep = "")
}

cli <- stage7a$parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% ""
if (action == "check") {
  stage7a$load_package()
  rows <- data_manifest(); fits <- fit_manifest(rows)
  assert(nrow(rows) == 2L && nrow(fits) == 24L &&
           all(table(fits$data_id) == 12L) &&
           identical(sort(unique(fits$fit_seed)), 97210301:97210312),
         "Stage 7B manifest or seed pairing is invalid.")
  cat("STAGE7B_CHECK_PASS datasets=2 fits=24\n")
} else if (action == "generate") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, FALSE, cli$`baseline-root` %||% "")
} else if (action == "fit") {
  assert(nzchar(output_root) && nzchar(cli$`fit-id` %||% ""),
         "--output-root and --fit-id are required.")
  stage7a$fit_one(output_root, cli$`fit-id`, FALSE)
} else if (action == "batch") {
  assert(nzchar(output_root), "--output-root is required.")
  batch(output_root, as.integer(cli$workers %||% "8"))
} else if (action == "smoke") {
  assert(nzchar(output_root), "--output-root is required.")
  generate(output_root, TRUE, "")
  for (fit_id in fit_manifest(data_manifest(), TRUE)$fit_id) {
    stage7a$fit_one(output_root, fit_id, TRUE)
  }
  stage7a$select_winners(output_root, TRUE)
  stage7a$lines_new(c("status=MICRO_SMOKE_PASS", "datasets=2", "fits=2",
                      "truth_unseal_authorized=FALSE"),
                    file.path(output_root, "MICRO_SMOKE_COMPLETE.txt"))
  cat("STAGE7B_MICRO_SMOKE_PASS\n")
} else if (action == "status") {
  assert(nzchar(output_root), "--output-root is required.")
  cat("complete=", length(list.files(file.path(output_root, "fits"),
                                        "^FIT_COMPLETE[.]txt$", recursive = TRUE)),
      " errors=", length(list.files(file.path(output_root, "fits"),
                                      "^FIT_ERROR[.]txt$", recursive = TRUE)),
      " selected=", file.exists(file.path(output_root,
                                          "truth_free_selection",
                                          "TRUTH_FREE_SELECTION_FROZEN.txt")),
      "\n", sep = "")
} else stop("Unknown action: ", action, call. = FALSE)
