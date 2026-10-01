#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", script_argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(
    script_file, winslash = "/", mustWork = TRUE
  ))
}

experiments_root <- dirname(script_root)
stage6a_root <- file.path(
  experiments_root, "g12-dimension-overspec-stage6a-20260831"
)
stage6b_root <- file.path(
  experiments_root, "g12-truth-free-dimension-selection-stage6b-20260901"
)
stage6f_root <- file.path(
  experiments_root, "g12-joint-dimension-selection-stage6f-20260911"
)
source(file.path(stage6a_root, "stage6a_common_20260831_v1.R"),
       local = FALSE)
source(file.path(stage6b_root, "stage6b_truth_free_tools_20260901_v1.R"),
       local = FALSE)
source(file.path(stage6f_root,
                 "stage6f_joint_selection_tools_20260911_v1.R"),
       local = FALSE)
source(file.path(script_root, "stage6g_tools_20260911_v1.R"),
       local = FALSE)

PP_EXPERIMENT_ID <-
  "G12_STAGE6G_DIMENSION_CRITERION_FALSIFICATION_V1_20260911"
stage6g_library <- pp_env("STAGE6G_R_LIB")
if (nzchar(stage6g_library)) Sys.setenv(STAGE6A_R_LIB = stage6g_library)
STAGE6G_WORKERS <- as.integer(pp_env("STAGE6G_WORKERS", "7"))
pp_assert(length(STAGE6G_WORKERS) == 1L && is.finite(STAGE6G_WORKERS) &&
            STAGE6G_WORKERS %in% 1:7,
          "STAGE6G_WORKERS must be an integer from 1 to 7.")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

cli <- pp_parse_cli(commandArgs(trailingOnly = TRUE))
action <- cli$action %||% "check"
output_root <- cli$`output-root` %||% pp_env("STAGE6G_OUTPUT_ROOT")

scenario_path <- file.path(script_root, "SCENARIOS_STAGE6G_V1.csv")
config_path <- file.path(script_root, "FIT_CONFIGS_STAGE6G_V1.csv")
split_path <- file.path(script_root, "HOLDOUT_SPLITS_STAGE6G_V1.csv")
start_path <- file.path(script_root, "START_SEEDS_STAGE6G_V1.csv")
package_binding_path <- file.path(
  script_root, "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv"
)
authorization_path <- file.path(
  script_root, "STAGE6G_EXECUTION_AUTHORIZATION_20260911.txt"
)

assert_execution_authorized <- function() {
  pp_assert(file.exists(authorization_path),
            "Stage 6G truth-free execution authorization is missing.")
  lines <- readLines(authorization_path, warn = FALSE)
  pp_assert(
    any(lines == "status=AUTHORIZED_TRUTH_FREE_STAGE6G_SCREEN") &&
      any(lines == "truth_unseal_authorized=FALSE") &&
      any(lines == "formal_paper_mc=FALSE"),
    "Stage 6G authorization file has the wrong scope."
  )
  invisible(TRUE)
}

require_output_root <- function(create = FALSE) {
  pp_assert(nzchar(output_root), "--output-root is required.")
  if (create && !dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
  }
  output_root
}

load_protocol_registration <- function() {
  scenarios <- s6g_read_csv(scenario_path)
  configs <- s6g_read_csv(config_path)
  splits <- s6g_read_csv(split_path)
  starts <- s6g_read_csv(start_path)
  s6g_validate_registration(scenarios, configs, splits, starts)
  list(
    scenarios = scenarios, configs = configs, splits = splits,
    starts = starts,
    fits = s6g_build_fit_manifest(scenarios, configs, splits, starts)
  )
}

runtime_paths <- function() list(
  fit_manifest = file.path(output_root, "FIT_MANIFEST.csv"),
  configs = file.path(output_root, "FIT_CONFIGS.csv"),
  splits = file.path(output_root, "HOLDOUT_SPLITS.csv"),
  starts = file.path(output_root, "START_SEEDS.csv"),
  data_ids = file.path(output_root, "DATA_IDS.csv"),
  prepared = file.path(output_root, "STAGE6G_PREPARED.txt")
)

load_runtime_registration <- function() {
  require_output_root(FALSE)
  paths <- runtime_paths()
  pp_assert(all(file.exists(unlist(paths[c(
    "fit_manifest", "configs", "splits", "starts", "data_ids"
  )], use.names = FALSE))),
  "Prepared Stage 6G truth-free manifests are incomplete.")
  fits <- pp_read_csv(paths$fit_manifest)
  configs <- pp_read_csv(paths$configs)
  splits <- pp_read_csv(paths$splits)
  starts <- pp_read_csv(paths$starts)
  data_ids <- pp_read_csv(paths$data_ids)
  forbidden <- c(
    "truth_basis_seed", "expected_truth_fit_config_id", "truth_L_f",
    "truth_L_s_1", "truth_L_s_2", "truth_M_common"
  )
  pp_assert(
    nrow(fits) == 120L && !anyDuplicated(fits$fit_id) &&
      length(unique(fits$selection_stratum_id)) == 30L &&
      all(table(fits$selection_stratum_id) == 4L) &&
      !length(intersect(names(fits), forbidden)) &&
      nrow(configs) == 5L && nrow(splits) == 2L &&
      nrow(starts) == 4L && nrow(data_ids) == 3L,
    "Prepared Stage 6G runtime registration is invalid or truth-bearing."
  )
  list(fits = fits, configs = configs, splits = splits,
       starts = starts, data_ids = data_ids)
}

verify_package_binding <- function() {
  pp_assert(file.exists(package_binding_path),
            "Stage 6G package-source binding is missing.")
  source_root <- pp_env("STAGE6G_PACKAGE_SOURCE_ROOT", PP_REPO_ROOT)
  pp_assert(dir.exists(source_root),
            "STAGE6G_PACKAGE_SOURCE_ROOT does not exist.")
  binding <- pp_read_csv(package_binding_path)
  pp_assert(nrow(binding) == 29L &&
              identical(names(binding), c("relative_path", "sha256", "bytes")) &&
              !anyDuplicated(binding$relative_path),
            "Stage 6G package-source binding is malformed.")
  paths <- file.path(source_root, binding$relative_path)
  info <- file.info(paths)
  matches <- file.exists(paths) & is.finite(info$size) &
    as.numeric(info$size) == as.numeric(binding$bytes) &
    vapply(paths, pp_sha256, character(1L)) == binding$sha256
  pp_assert(all(matches), paste0(
    "Stage 6G package-source binding mismatch: ",
    paste(binding$relative_path[!matches], collapse = ";")
  ))
  list(binding = binding, source_root = source_root)
}

verify_environment <- function(read_truth_registration = FALSE,
                               require_prepared = FALSE) {
  assert_execution_authorized()
  pp_assert(requireNamespace("digest", quietly = TRUE),
            "The digest package is required.")
  pp_load_package()
  pp_source_runtime(FALSE)
  package_binding <- verify_package_binding()
  pp_assert(requireNamespace("Splinets", quietly = TRUE) &&
              as.character(utils::packageVersion("Splinets")) == "1.5.1",
            "Stage 6G requires the bound Splinets 1.5.1 environment.")
  registration <- if (read_truth_registration) {
    load_protocol_registration()
  } else if (require_prepared) {
    load_runtime_registration()
  } else NULL
  list(registration = registration,
       package_binding = package_binding)
}

data_paths <- function(data_id) list(
  observation = file.path(output_root, "data", data_id,
                          "observation_bundle.rds"),
  truth = file.path(output_root, "data", data_id, "sealed_truth",
                    "truth_bundle.rds"),
  assertions = file.path(output_root, "data", data_id, "sealed_truth",
                         "truth_assertions.csv"),
  seal = file.path(output_root, "data", data_id,
                   "DATA_SEAL_COMPLETE.txt")
)

split_paths <- function(data_id, split_id) list(
  train = file.path(output_root, "data", data_id, "splits", split_id,
                    "train_observation.rds"),
  holdout = file.path(output_root, "data", data_id, "splits", split_id,
                      "holdout_plan.rds")
)

assert_truth_locked <- function() {
  pp_assert(!file.exists(file.path(output_root,
                                   "TRUTH_UNSEAL_AUTHORIZATION.txt")),
            "Refusing Stage 6G truth-free work after truth authorization.")
  invisible(TRUE)
}

s6g_generator_config <- function(row) {
  L_f <- as.integer(row$truth_L_f[[1L]])
  L_s <- c(as.integer(row$truth_L_s_1[[1L]]),
           as.integer(row$truth_L_s_2[[1L]]))
  M <- as.integer(row$truth_M_common[[1L]])
  list(
    S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
    L_f = L_f, L_s = L_s, M_f = rep(M, L_f),
    M_s = list(rep(M, L_s[[1L]]), rep(M, L_s[[2L]])),
    K = 5L, n_obs = 8L, common_grid = FALSE,
    n_obs_range = c(6L, 9L), sigma_eps = 0.3,
    bool_sparse_loadings = TRUE, prop_sparse = 0.9,
    score_var_decay = TRUE, n_dense = 201L,
    mean_amp = 0.6, beta_amp = 0, bs_degree = c(2L, 3L),
    identified_loadings = TRUE, sparsity_mode = "fixed",
    seed = as.integer(row$data_seed[[1L]]),
    truth_basis_method = "jaoua_splinets",
    truth_basis_seed = as.integer(row$truth_basis_seed[[1L]])
  )
}

s6g_truth_assertions <- function(data, row, config) {
  truth <- data$true_params
  L_f <- as.integer(row$truth_L_f[[1L]])
  L_s <- c(as.integer(row$truth_L_s_1[[1L]]),
           as.integer(row$truth_L_s_2[[1L]]))
  M <- as.integer(row$truth_M_common[[1L]])
  counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  metadata <- c(
    truth$phi_basis_meta,
    unlist(truth$psi_basis_meta, recursive = FALSE)
  )
  continuous_error <- max(vapply(metadata, function(value) {
    max(abs(value$gram - diag(nrow(value$gram))))
  }, numeric(1L)))
  loadings <- c(
    lapply(seq_len(ncol(truth$a_true)), function(index) {
      truth$a_true[, index]
    }),
    unlist(lapply(truth$b_true, function(matrix_value) {
      lapply(seq_len(ncol(matrix_value)), function(index) {
        matrix_value[, index]
      })
    }), recursive = FALSE)
  )
  active_counts <- vapply(loadings, function(value) sum(value != 0),
                          integer(1L))
  loading_gram_error <- max(vapply(seq_len(truth$S), function(study) {
    combined <- cbind(truth$a_true, truth$b_true[[study]])
    gram <- crossprod(combined)
    max(abs(gram - diag(diag(gram))))
  }, numeric(1L)))
  response_error <- 0
  signal_error <- 0
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      signal <- truth$mu_true_values[[study]][[subject]] +
        truth$beta_true_values[[study]][[subject]] +
        truth$f_true_values[[study]][[subject]] %*% t(truth$a_true) +
        truth$g_true_values[[study]][[subject]] %*%
          t(truth$b_true[[study]])
      observed <- do.call(cbind, data$Y[[study]][[subject]])
      signal_error <- max(
        signal_error,
        abs(signal - truth$signal_true_values[[study]][[subject]])
      )
      response_error <- max(
        response_error,
        abs(observed - signal - truth$noise_true_values[[study]][[subject]])
      )
    }
  }
  expected_lambda <- 1 / seq_len(M)^2
  expected_lambda <- expected_lambda / sum(expected_lambda)
  all_lambda <- c(
    truth$lambda_phi_true,
    unlist(truth$lambda_psi_true, recursive = FALSE)
  )
  lambda_error <- max(vapply(all_lambda, function(value) {
    max(abs(as.numeric(value) - expected_lambda))
  }, numeric(1L)))
  checks <- data.frame(
    assertion_id = c(
      "dimensions", "truth_basis_and_seed", "factor_fpca_counts",
      "continuous_gram", "observation_counts", "active_loading_counts",
      "loading_gram_identification", "score_variance_rule",
      "signal_identity", "response_identity"
    ),
    passed = c(
      truth$S == 2L && all(truth$n_s == c(30L, 30L)) &&
        truth$p == 500L && truth$d_use == 0L,
      identical(truth$truth_basis_method, "jaoua_splinets") &&
        as.integer(truth$truth_basis_seed) ==
          as.integer(row$truth_basis_seed[[1L]]),
      truth$L_f == L_f &&
        identical(as.integer(truth$L_s_by_study), L_s) &&
        identical(as.integer(truth$M_f), rep(M, L_f)) &&
        identical(lapply(truth$M_s, as.integer),
                  list(rep(M, L_s[[1L]]), rep(M, L_s[[2L]]))) &&
        length(metadata) == L_f + sum(L_s),
      continuous_error < 1e-12,
      all(counts >= 6L & counts <= 9L),
      all(active_counts == 50L),
      isTRUE(truth$identified_loadings) && loading_gram_error < 1e-12,
      lambda_error < 1e-14,
      signal_error < 1e-12,
      response_error < 1e-12
    ),
    detail = c(
      "S=2;n=(30,30);p=500;d=0",
      paste0("jaoua_splinets;seed=", row$truth_basis_seed[[1L]]),
      paste0("L=", L_f, "/", paste(L_s, collapse = "/"), ";M=", M),
      format(continuous_error, digits = 17),
      paste(range(counts), collapse = "-"),
      paste(active_counts, collapse = ";"),
      format(loading_gram_error, digits = 17),
      format(lambda_error, digits = 17),
      format(signal_error, digits = 17),
      format(response_error, digits = 17)
    ),
    stringsAsFactors = FALSE
  )
  checks
}

write_runtime_bindings <- function(package_binding) {
  relative <- c(
    "SCENARIOS_STAGE6G_V1.csv", "FIT_CONFIGS_STAGE6G_V1.csv",
    "HOLDOUT_SPLITS_STAGE6G_V1.csv", "START_SEEDS_STAGE6G_V1.csv",
    "BUDGET_STAGE6G_V1.csv",
    "PROTOCOL_G12_STAGE6G_DIMENSION_CRITERION_FALSIFICATION_V1_20260911.md",
    "STAGE6G_EXECUTION_AUTHORIZATION_20260911.txt",
    "PACKAGE_SOURCE_BINDING_STAGE6G_V1.csv",
    "stage6g_tools_20260911_v1.R", "run_stage6g_20260911_v1.R",
    file.path("..", "g12-dimension-overspec-stage6a-20260831",
              "stage6a_common_20260831_v1.R"),
    file.path("..", "g12-truth-free-dimension-selection-stage6b-20260901",
              "stage6b_truth_free_tools_20260901_v1.R"),
    file.path("..", "g12-joint-dimension-selection-stage6f-20260911",
              "stage6f_joint_selection_tools_20260911_v1.R")
  )
  paths <- file.path(script_root, relative)
  pp_assert(all(file.exists(paths)),
            "A Stage 6G runtime source file is missing.")
  info <- file.info(paths)
  pp_write_csv(data.frame(
    relative_to_stage6g_script_root = gsub("\\\\", "/", relative),
    bytes = as.numeric(info$size),
    sha256 = vapply(paths, pp_sha256, character(1L)),
    stringsAsFactors = FALSE
  ), file.path(output_root, "RUNTIME_SOURCE_BINDING.csv"))
  copied <- package_binding$binding
  copied$package_source_root_at_runtime <- package_binding$source_root
  pp_write_csv(copied,
               file.path(output_root, "PACKAGE_SOURCE_BINDING_VERIFIED.csv"))
  session <- capture.output(utils::sessionInfo())
  pp_write_lines(session, file.path(output_root, "SESSION_INFO.txt"))
  invisible(TRUE)
}

s6g_generate_and_seal <- function(row) {
  data_id <- as.character(row$scenario_id[[1L]])
  paths <- data_paths(data_id)
  output_dir <- dirname(paths$observation)
  if (dir.exists(output_dir)) {
    pp_assert(all(file.exists(unlist(paths, use.names = FALSE))),
              paste0("Incomplete existing Stage 6G data directory: ", data_id))
    return(invisible(FALSE))
  }
  config <- s6g_generator_config(row)
  generator <- getExportedValue("multiFSYNC", "simulate_multi_study_structured")
  captured <- c2_capture_conditions(
    do.call(generator, config), phase = "stage6g_data_generation"
  )
  pp_assert(is.null(captured$error), paste0(
    "Stage 6G data generation failed for ", data_id, ": ",
    captured$error$message %||% "unknown"
  ))
  data <- captured$value
  assertions <- s6g_truth_assertions(data, row, config)
  pp_assert(all(assertions$passed), paste0(
    "Stage 6G truth assertions failed for ", data_id, ": ",
    paste(assertions$assertion_id[!assertions$passed], collapse = ";")
  ))
  observation <- list(
    bundle_class = "stage6g_observation_only_v1",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_id,
    data_index = as.integer(row$scenario_index[[1L]]),
    scenario_label_available_to_model = FALSE,
    engineering_pilot = TRUE,
    formal_experiment = FALSE,
    Y = data$Y, time_obs = data$time_obs, Z = data$Z,
    dimensions = list(
      S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L, K = 5L
    )
  )
  pp_assert(!length(candidate2_forbidden_observation_names(observation)) &&
              !length(s6b_forbidden_observation_names(observation)) &&
              !length(intersect(
                names(observation$dimensions),
                c("L_f", "L_s", "M_f", "M_s", "truth_L", "truth_M")
              )),
            "Stage 6G observation bundle leaks structural truth.")
  sealed <- list(
    bundle_class = "stage6g_sealed_truth_v1",
    experiment_id = PP_EXPERIMENT_ID,
    data_id = data_id,
    data_index = as.integer(row$scenario_index[[1L]]),
    engineering_pilot = TRUE,
    formal_experiment = FALSE,
    data = data,
    sealed_generation_provenance = list(
      generator_config = config,
      generation_warnings = captured$warnings,
      truth_assertions = assertions,
      generated_at_utc = captured$ended_at_utc
    )
  )
  dir.create(file.path(output_dir, "sealed_truth"), recursive = TRUE,
             showWarnings = FALSE, mode = "0700")
  pp_save_rds(observation, paths$observation)
  pp_save_rds(sealed, paths$truth)
  pp_write_csv(assertions, paths$assertions)
  warnings <- if (is.data.frame(captured$warnings)) captured$warnings else
    data.frame(class = character(), message = character(),
               call = character(), phase = character(),
               stringsAsFactors = FALSE)
  pp_write_csv(warnings, file.path(output_dir, "GENERATION_WARNINGS.csv"))
  pp_write_lines(c(
    "status=DATA_SEAL_COMPLETE", paste0("data_id=", data_id),
    paste0("observation_sha256=", pp_sha256(paths$observation)),
    paste0("sealed_truth_sha256=", pp_sha256(paths$truth)),
    "observation_contains_L_or_M=FALSE",
    "truth_unseal_authorized=FALSE", "engineering_pilot=TRUE",
    "formal_paper_mc_result=FALSE"
  ), paths$seal)
  invisible(TRUE)
}

s6g_prepare_splits <- function(data_id) {
  observation <- readRDS(data_paths(data_id)$observation)
  pp_assert(identical(observation$bundle_class,
                      "stage6g_observation_only_v1") &&
              !length(s6b_forbidden_observation_names(observation)),
            paste0("Invalid Stage 6G observation bundle: ", data_id))
  middle <- s6f_make_time_holdout(observation, "middle")
  upper <- s6f_make_time_holdout(observation, "upper")
  middle_index <- vapply(middle$holdout_plan$records, `[[`, integer(1L),
                         "original_index")
  upper_index <- vapply(upper$holdout_plan$records, `[[`, integer(1L),
                        "original_index")
  pp_assert(
    length(middle_index) == 60L && length(upper_index) == 60L &&
      all(upper_index == middle_index + 1L) &&
      all(vapply(middle$holdout_plan$records, `[[`, integer(1L),
                 "train_time_count") >= 5L) &&
      all(vapply(upper$holdout_plan$records, `[[`, integer(1L),
                 "train_time_count") >= 5L),
    paste0("Stage 6G holdout contract failed for ", data_id)
  )
  values <- list(split1_middle = middle, split2_upper = upper)
  summary_rows <- list()
  for (split_id in names(values)) {
    paths <- split_paths(data_id, split_id)
    pp_save_rds(values[[split_id]]$train_observation, paths$train)
    pp_save_rds(values[[split_id]]$holdout_plan, paths$holdout)
    records <- values[[split_id]]$holdout_plan$records
    summary_rows[[split_id]] <- data.frame(
      data_id = data_id, split_id = split_id,
      study = vapply(records, `[[`, integer(1L), "study"),
      subject = vapply(records, `[[`, integer(1L), "subject"),
      original_index = vapply(records, `[[`, integer(1L), "original_index"),
      holdout_time = vapply(records, `[[`, numeric(1L), "time"),
      train_time_count = vapply(records, `[[`, integer(1L),
                                "train_time_count"),
      variables = vapply(records, `[[`, integer(1L), "p"),
      truth_used = FALSE, stringsAsFactors = FALSE
    )
  }
  pp_write_csv(
    do.call(rbind, summary_rows),
    file.path(output_root, "data", data_id, "splits",
              "HOLDOUT_SPLIT_SUMMARY.csv")
  )
  invisible(TRUE)
}

prepare_stage6g <- function() {
  require_output_root(TRUE)
  assert_truth_locked()
  verification <- verify_environment(read_truth_registration = TRUE)
  registration <- verification$registration
  paths <- runtime_paths()
  if (file.exists(paths$prepared)) {
    expected <- c(
      unlist(paths[c("fit_manifest", "configs", "splits", "starts",
                     "data_ids")], use.names = FALSE),
      unlist(lapply(registration$scenarios$scenario_id, function(data_id) {
        c(unlist(data_paths(data_id), use.names = FALSE),
          unlist(split_paths(data_id, "split1_middle"), use.names = FALSE),
          unlist(split_paths(data_id, "split2_upper"), use.names = FALSE))
      }), use.names = FALSE)
    )
    pp_assert(all(file.exists(expected)),
              "Existing Stage 6G preparation is incomplete.")
    cat("STAGE6G_ALREADY_PREPARED\n")
    return(invisible(TRUE))
  }
  pp_assert(!dir.exists(file.path(output_root, "fits")) &&
              !dir.exists(file.path(output_root, "truth_free_selection")),
            "Refusing Stage 6G preparation over fit/selection artifacts.")
  pp_write_csv(registration$fits, paths$fit_manifest)
  pp_write_csv(registration$configs, paths$configs)
  pp_write_csv(registration$splits, paths$splits)
  pp_write_csv(registration$starts, paths$starts)
  pp_write_csv(data.frame(
    data_index = as.integer(registration$scenarios$scenario_index),
    data_id = as.character(registration$scenarios$scenario_id),
    stringsAsFactors = FALSE
  ), paths$data_ids)
  sealed_registration <- file.path(output_root, "sealed_registration")
  dir.create(sealed_registration, recursive = TRUE, showWarnings = FALSE,
             mode = "0700")
  pp_write_csv(registration$scenarios, file.path(
    sealed_registration, "SCENARIOS_WITH_GENERATING_TRUTH.csv"
  ))
  write_runtime_bindings(verification$package_binding)
  for (index in seq_len(nrow(registration$scenarios))) {
    row <- registration$scenarios[index, , drop = FALSE]
    s6g_generate_and_seal(row)
    s6g_prepare_splits(row$scenario_id[[1L]])
  }
  pp_write_lines(c(
    "status=STAGE6G_PREPARED", paste0("prepared_utc=", pp_iso_time()),
    "data_sets=3", "registered_endpoints=120",
    "configuration_split_strata=30", "starts_per_stratum=4",
    "observation_bundles_exclude_generating_L_and_M=TRUE",
    "truth_read_by_fit_or_selection=FALSE",
    "truth_unseal_authorized=FALSE", "continuation_authorized=FALSE",
    "formal_paper_mc_result=FALSE"
  ), paths$prepared)
  cat("STAGE6G_PREPARED data=3 endpoints=120 truth_sealed=TRUE\n")
  invisible(TRUE)
}

load_split <- function(data_id, split_id) {
  pp_assert(split_id %in% c("split1_middle", "split2_upper"),
            "Unknown Stage 6G split_id.")
  paths <- split_paths(data_id, split_id)
  pp_assert(file.exists(paths$train) && file.exists(paths$holdout),
            paste0("Missing Stage 6G split: ", data_id, "/", split_id))
  train <- readRDS(paths$train)
  holdout <- readRDS(paths$holdout)
  pp_assert(!length(s6b_forbidden_observation_names(train)) &&
              identical(holdout$truth_used, FALSE) &&
              holdout$subject_time_count == 60L,
            "Stage 6G split is invalid or truth-bearing.")
  list(train = train, holdout = holdout)
}

s6g_set_warning_phase <- function(warnings, phase) {
  if (!is.data.frame(warnings) || !nrow(warnings)) {
    return(data.frame(
      class = character(), message = character(), call = character(),
      phase = character(), stringsAsFactors = FALSE
    ))
  }
  warnings$phase <- phase
  warnings
}

write_fit_error <- function(fit_dir, fit_id, message) {
  pp_write_lines(c(
    "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
    paste0("message=", message), paste0("ended_utc=", pp_iso_time()),
    "truth_read=FALSE"
  ), file.path(fit_dir, "FIT_ERROR.txt"))
}

run_fit <- function(fit_id) {
  require_output_root(FALSE)
  assert_truth_locked()
  verify_environment(require_prepared = TRUE)
  registration <- load_runtime_registration()
  pp_assert(file.exists(runtime_paths()$prepared),
            "Stage 6G must be prepared before fitting.")
  row <- registration$fits[registration$fits$fit_id == fit_id,
                           , drop = FALSE]
  pp_assert(nrow(row) == 1L, "Stage 6G fit_id is not registered.")
  fit_dir <- file.path(output_root, "fits", fit_id)
  if (dir.exists(fit_dir)) {
    complete <- file.exists(file.path(fit_dir, "FIT_COMPLETE.txt"))
    error <- file.exists(file.path(fit_dir, "FIT_ERROR.txt"))
    pp_assert(xor(complete, error) &&
                file.exists(file.path(fit_dir, "terminal_record.rds")),
              paste0("Incomplete existing Stage 6G fit directory: ", fit_id))
    cat(if (complete) "FIT_ALREADY_COMPLETE fit_id=" else
          "FIT_ERROR_RETAINED fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  split <- load_split(row$data_id[[1L]], row$split_id[[1L]])
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  result <- pp_fit_g12(split$train, row, smoke = FALSE)
  if (!is.null(result$fit)) {
    pp_save_rds(result$fit, file.path(fit_dir, "fit.rds"))
  }
  result$captured$warnings <- s6g_set_warning_phase(
    result$captured$warnings, "stage6g_g12_fixed400"
  )
  record <- pp_make_terminal_record(result, row, smoke = FALSE)
  record$terminal_schema <- "G12_STAGE6G_TERMINAL_V1"
  record$experiment_id <- PP_EXPERIMENT_ID
  record$split_id <- as.character(row$split_id[[1L]])
  record$holdout_scoring_complete <- FALSE
  record$holdout_available_to_start_selection <- FALSE
  record$truth_used_for_holdout_scoring <- FALSE
  pp_save_rds(record, file.path(fit_dir, "terminal_record.rds"))
  pp_write_csv(result$captured$warnings, file.path(fit_dir, "WARNINGS.csv"))
  if (is.null(result$fit)) {
    write_fit_error(fit_dir, fit_id,
                    record$error$message %||% "unknown endpoint error")
    cat("FIT_ERROR_RETAINED fit_id=", fit_id, "\n", sep = "")
    return(invisible(TRUE))
  }
  pp_write_csv(data.frame(
    t1_sweep = seq_along(result$fit$ELBO),
    ordinary_elbo = as.numeric(result$fit$ELBO)
  ), file.path(fit_dir, "ORDINARY_T1_ELBO.csv"))
  if (is.data.frame(result$fit$practical_diagnostics)) {
    pp_write_csv(result$fit$practical_diagnostics,
                 file.path(fit_dir, "PRACTICAL_DIAGNOSTICS.csv"))
  }
  if (is.data.frame(result$fit$random_scale_calibration_diagnostics)) {
    pp_write_csv(result$fit$random_scale_calibration_diagnostics,
                 file.path(fit_dir, "INITIALIZATION_DIAGNOSTICS.csv"))
  }
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    "terminal_status=fixed_400_complete",
    paste0("objective_eligible=", record$objective_eligible),
    paste0("elapsed_seconds=", format(record$elapsed_seconds, digits = 12)),
    "ordinary_T1_sweeps=400", "n_cpus=1", "holdout_scored=FALSE",
    "truth_read=FALSE", "continuation_used=FALSE",
    paste0("ended_utc=", pp_iso_time())
  ), file.path(fit_dir, "FIT_COMPLETE.txt"))
  cat("FIT_COMPLETE fit_id=", fit_id, " elapsed_seconds=",
      format(record$elapsed_seconds, digits = 8), "\n", sep = "")
  invisible(TRUE)
}

collect_endpoints <- function(manifest) {
  record_paths <- file.path(output_root, "fits", manifest$fit_id,
                            "terminal_record.rds")
  complete_paths <- file.path(output_root, "fits", manifest$fit_id,
                              "FIT_COMPLETE.txt")
  error_paths <- file.path(output_root, "fits", manifest$fit_id,
                           "FIT_ERROR.txt")
  pp_assert(all(file.exists(record_paths)) &&
              all(file.exists(complete_paths) | file.exists(error_paths)),
            "Stage 6G lacks one or more terminal endpoint records.")
  records <- lapply(record_paths, readRDS)
  pp_assert(all(vapply(records, function(record) {
    identical(record$terminal_schema, "G12_STAGE6G_TERMINAL_V1") &&
      !isTRUE(record$truth_used_for_fit_stopping_or_selection) &&
      !isTRUE(record$truth_used_for_holdout_scoring) &&
      !isTRUE(record$holdout_scoring_complete)
  }, logical(1L))), "Stage 6G terminal truth/score contract failed.")
  rows <- do.call(rbind, lapply(records, pp_terminal_row))
  rows <- rows[match(manifest$fit_id, rows$fit_id), , drop = FALSE]
  rows$split_id <- manifest$split_id
  rows$complexity_proxy <- manifest$complexity_proxy
  rows$endpoint_disposition <- ifelse(file.exists(complete_paths),
                                      "complete", "error_retained")
  pp_assert(identical(as.character(rows$fit_id),
                      as.character(manifest$fit_id)),
            "Stage 6G terminal-to-manifest join failed.")
  list(records = records, rows = rows)
}

freeze_start_winners <- function(registration) {
  selection_root <- file.path(output_root, "truth_free_selection")
  marker <- file.path(selection_root, "START_SELECTION_FROZEN.txt")
  winner_path <- file.path(selection_root, "CONFIG_SPLIT_WINNERS.csv")
  if (file.exists(marker)) return(pp_read_csv(winner_path))
  pp_assert(!dir.exists(selection_root),
            "Incomplete Stage 6G truth-free selection directory exists.")
  endpoints <- collect_endpoints(registration$fits)$rows
  groups <- split(endpoints, endpoints$selection_stratum_id)
  eligible_counts <- vapply(groups, function(rows) {
    sum(rows$objective_eligible & is.finite(rows$final_elbo) &
          rows$terminal_status == "fixed_400_complete")
  }, integer(1L))
  if (any(eligible_counts == 0L)) {
    dir.create(selection_root, recursive = TRUE, showWarnings = FALSE,
               mode = "0700")
    pp_write_csv(data.frame(
      selection_stratum_id = names(eligible_counts),
      eligible_endpoints = unname(eligible_counts),
      stringsAsFactors = FALSE
    ), file.path(selection_root, "ELIGIBLE_ENDPOINT_COUNTS.csv"))
    pp_write_csv(endpoints, file.path(selection_root, "ALL_ENDPOINTS.csv"))
    pp_abort("At least one Stage 6G stratum has no eligible endpoint.")
  }
  winners <- pp_select_g12_winners(endpoints)
  pp_assert(nrow(winners) == 30L &&
              length(unique(winners$selection_stratum_id)) == 30L &&
              all(table(winners$data_id) == 10L),
            "Stage 6G start selection did not yield 30 winners.")
  dir.create(selection_root, recursive = TRUE, showWarnings = FALSE,
             mode = "0700")
  pp_write_csv(endpoints, file.path(selection_root, "ALL_ENDPOINTS.csv"))
  pp_write_csv(data.frame(
    selection_stratum_id = names(eligible_counts),
    eligible_endpoints = unname(eligible_counts),
    stringsAsFactors = FALSE
  ), file.path(selection_root, "ELIGIBLE_ENDPOINT_COUNTS.csv"))
  pp_write_csv(winners, winner_path)
  pp_write_lines(c(
    "status=START_SELECTION_FROZEN", paste0("completed_utc=", pp_iso_time()),
    "endpoints=120", "split_configuration_winners=30",
    paste0("winner_table_sha256=", pp_sha256(winner_path)),
    "winner_rule=maximum_eligible_ordinary_T1_ELBO_within_data_split_configuration",
    "heldout_response_read_for_start_selection=FALSE",
    "cross_dimension_elbo_comparison=FALSE", "selection_used_truth=FALSE"
  ), marker)
  winners
}

score_winner <- function(winner) {
  fit_id <- winner$fit_id[[1L]]
  score_dir <- file.path(output_root, "truth_free_scores", fit_id)
  complete <- file.path(score_dir, "SCORE_COMPLETE.txt")
  if (dir.exists(score_dir)) {
    pp_assert(file.exists(complete),
              paste0("Incomplete Stage 6G score directory: ", fit_id))
    return(invisible(TRUE))
  }
  split <- load_split(winner$data_id[[1L]], winner$split_id[[1L]])
  fit_path <- file.path(output_root, "fits", fit_id, "fit.rds")
  pp_assert(file.exists(fit_path), paste0("Winner fit is missing: ", fit_id))
  fit <- readRDS(fit_path)
  score <- tryCatch(
    s6b_score_holdout(fit, split$holdout, "all"),
    error = function(condition) condition
  )
  dir.create(score_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  if (inherits(score, "condition")) {
    pp_write_lines(c(
      "status=ERROR_RETAINED_NO_SELECTIVE_RERUN", paste0("fit_id=", fit_id),
      paste0("message=", conditionMessage(score)), "truth_read=FALSE"
    ), file.path(score_dir, "SCORE_ERROR.txt"))
    pp_abort("Stage 6G winner scoring failed: ", fit_id)
  }
  decorate <- function(value) {
    value$data_id <- winner$data_id[[1L]]
    value$fit_id <- fit_id
    value$fit_config_id <- winner$fit_config_id[[1L]]
    value$split_id <- winner$split_id[[1L]]
    value
  }
  pp_write_csv(decorate(score$subject_time),
               file.path(score_dir, "HOLDOUT_SUBJECT_TIME_MSE.csv"))
  pp_write_csv(decorate(score$aggregate),
               file.path(score_dir, "HOLDOUT_SUMMARY.csv"))
  pp_write_lines(c(
    "status=COMPLETE", paste0("fit_id=", fit_id),
    paste0("fit_config_id=", winner$fit_config_id[[1L]]),
    paste0("split_id=", winner$split_id[[1L]]),
    "fit_selected_before_holdout_score=TRUE", "truth_read=FALSE",
    paste0("completed_utc=", pp_iso_time())
  ), complete)
  invisible(TRUE)
}

score_winners <- function(winners) {
  for (index in seq_len(nrow(winners))) {
    score_winner(winners[index, , drop = FALSE])
  }
  pp_write_lines(c(
    "status=WINNER_SCORING_COMPLETE", "winners_scored=30",
    "score_errors=0", "truth_read=FALSE",
    paste0("completed_utc=", pp_iso_time())
  ), file.path(output_root, "truth_free_scores",
               "WINNER_SCORING_COMPLETE.txt"))
  invisible(TRUE)
}

read_winner_scores <- function(winners) {
  values <- lapply(seq_len(nrow(winners)), function(index) {
    path <- file.path(
      output_root, "truth_free_scores", winners$fit_id[[index]],
      "HOLDOUT_SUBJECT_TIME_MSE.csv"
    )
    pp_assert(file.exists(path), paste0("Missing Stage 6G winner score: ", path))
    pp_read_csv(path)
  })
  do.call(rbind, values)
}

freeze_dimensions <- function(registration, winners) {
  selection_root <- file.path(output_root, "truth_free_selection")
  final_marker <- file.path(selection_root,
                            "DIMENSION_SELECTION_FROZEN.txt")
  if (file.exists(final_marker)) return(invisible(TRUE))
  score_marker <- file.path(output_root, "truth_free_scores",
                            "WINNER_SCORING_COMPLETE.txt")
  pp_assert(file.exists(score_marker) &&
              file.exists(file.path(selection_root,
                                    "START_SELECTION_FROZEN.txt")),
            "Stage 6G start selection or scoring is incomplete.")
  scores <- read_winner_scores(winners)
  pp_assert(nrow(scores) == 1800L && !any(scores$truth_used),
            "Stage 6G must contain 30 x 60 truth-free score rows.")
  selections <- list()
  for (data_index in seq_len(nrow(registration$data_ids))) {
    data_id <- registration$data_ids$data_id[[data_index]]
    rows <- scores[scores$data_id == data_id, , drop = FALSE]
    pp_assert(nrow(rows) == 600L &&
                length(unique(rows$fit_config_id)) == 5L &&
                length(unique(rows$split_id)) == 2L,
              paste0("Incomplete Stage 6G scores for ", data_id))
    subject_losses <- s6f_subject_cluster_losses(
      rows, c("mse", "mae", "observed_second_moment"),
      expected_splits = 2L
    )
    selection <- s6g_select_one_se(
      subject_losses, registration$configs, "mse", expected_subjects = 60L
    )
    reference_observed <- subject_losses[
      subject_losses$fit_config_id ==
        registration$configs$fit_config_id[[1L]], , drop = FALSE
    ]
    observed_rms <- sqrt(s6g_stratified_mean(
      reference_observed$observed_second_moment, reference_observed$study
    ))
    selection$summary$heldout_rmse <- sqrt(selection$summary$mean_loss)
    selection$summary$heldout_nrmse <-
      selection$summary$heldout_rmse / observed_rms
    ranking <- s6g_split_ranking_stability(rows, "mse")
    contrasts <- s6g_paired_contrasts(subject_losses, selection)
    loo <- s6g_leave_one_subject_out(subject_losses, registration$configs)
    data_dir <- file.path(selection_root, data_id)
    pp_write_csv(subject_losses,
                 file.path(data_dir, "SUBJECT_CLUSTERED_MSE.csv"))
    pp_write_csv(selection$summary,
                 file.path(data_dir, "MSE_ONE_SE_CANDIDATES.csv"))
    pp_write_csv(contrasts,
                 file.path(data_dir, "PAIRED_MSE_DIFFERENCES.csv"))
    pp_write_csv(ranking$risks,
                 file.path(data_dir, "SPLIT_SPECIFIC_RISKS.csv"))
    pp_write_csv(ranking$summary,
                 file.path(data_dir, "CROSS_SPLIT_RANK_STABILITY.csv"))
    pp_write_csv(ranking$pairwise,
                 file.path(data_dir, "CROSS_SPLIT_PAIRWISE_ORDERING.csv"))
    pp_write_csv(loo$decisions,
                 file.path(data_dir, "LEAVE_ONE_SUBJECT_OUT_SELECTIONS.csv"))
    pp_write_csv(loo$counts,
                 file.path(data_dir, "LEAVE_ONE_SUBJECT_OUT_COUNTS.csv"))
    final <- data.frame(
      data_id = data_id,
      selected_config_id = selection$selected_config_id,
      minimum_loss_config_id = selection$minimum_loss_config_id,
      one_se_threshold = selection$threshold,
      observed_rms = observed_rms,
      split_rank_spearman = ranking$summary$spearman_rank_correlation,
      split_minimum_winner_agreement =
        ranking$summary$minimum_loss_winner_agreement,
      within_configuration_start_rule =
        "maximum_eligible_ordinary_T1_ELBO",
      cross_dimension_rule =
        "paired_two_split_subject_clustered_MSE_one_SE_simplest",
      selection_used_truth = FALSE,
      continuation_started = FALSE,
      stringsAsFactors = FALSE
    )
    pp_write_csv(final,
                 file.path(data_dir, "FINAL_TRUTH_FREE_SELECTION.csv"))
    selections[[data_index]] <- final
  }
  selections <- do.call(rbind, selections)
  final_path <- file.path(selection_root,
                          "FINAL_TRUTH_FREE_DIMENSION_SELECTIONS.csv")
  pp_write_csv(selections, final_path)
  artifact_paths <- c(
    list.files(selection_root, pattern = "[.](csv|txt)$", recursive = TRUE,
               full.names = TRUE, ignore.case = TRUE),
    list.files(file.path(output_root, "truth_free_scores"),
               pattern = "[.](csv|txt)$", recursive = TRUE,
               full.names = TRUE, ignore.case = TRUE)
  )
  artifact_paths <- sort(unique(artifact_paths))
  hash_path <- file.path(selection_root, "TRUTH_FREE_ARTIFACT_HASHES.csv")
  artifact_paths <- setdiff(artifact_paths, hash_path)
  pp_write_csv(data.frame(
    relative_to_output_root = gsub(
      "\\\\", "/", substring(
        artifact_paths, nchar(normalizePath(output_root)) + 2L
      )
    ),
    bytes = as.numeric(file.info(artifact_paths)$size),
    sha256 = vapply(artifact_paths, pp_sha256, character(1L)),
    stringsAsFactors = FALSE
  ), hash_path)
  pp_write_lines(c(
    "status=DIMENSION_SELECTION_FROZEN",
    paste0("completed_utc=", pp_iso_time()),
    paste0("start_winner_table_sha256=", pp_sha256(file.path(
      selection_root, "CONFIG_SPLIT_WINNERS.csv"
    ))),
    paste0("dimension_selection_sha256=", pp_sha256(final_path)),
    paste0("artifact_hash_manifest_sha256=", pp_sha256(hash_path)),
    "data_sets=3", "candidate_dimensions_per_data=5",
    "dimension_selection_used_truth=FALSE", "truth_read=FALSE",
    "new_fits=120", "continuations=0", "formal_paper_mc_result=FALSE"
  ), final_marker)
  pp_write_lines(c(
    "status=STAGE6G_TRUTH_FREE_COMPLETE", "registered_endpoints=120",
    paste0("endpoint_errors=", sum(collect_endpoints(
      registration$fits
    )$rows$endpoint_disposition == "error_retained")),
    "split_configuration_winners=30", "dimension_selections=3",
    "truth_read=FALSE", "truth_unseal_authorized=FALSE",
    "continuations=0", paste0("completed_utc=", pp_iso_time())
  ), file.path(output_root, "STAGE6G_TRUTH_FREE_COMPLETE.txt"))
  invisible(selections)
}

write_status <- function(state, detail = "") {
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  pp_write_csv(data.frame(
    timestamp_utc = pp_iso_time(), state = state, detail = detail,
    registered_endpoints = 120L, complete_endpoints = complete,
    error_endpoints = errors, outer_workers = STAGE6G_WORKERS,
    per_fit_n_cpus = 1L, truth_read_by_fit_or_selection = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "CURRENT_STATUS.csv"))
}

summarize_resources <- function(registration) {
  endpoints <- collect_endpoints(registration$fits)$rows
  groups <- split(endpoints, endpoints$fit_config_id)
  runtime <- do.call(rbind, lapply(groups, function(rows) data.frame(
    fit_config_id = rows$fit_config_id[[1L]], endpoints = nrow(rows),
    complete = sum(rows$endpoint_disposition == "complete"),
    errors = sum(rows$endpoint_disposition == "error_retained"),
    elapsed_median_seconds = stats::median(rows$elapsed_seconds, na.rm = TRUE),
    elapsed_max_seconds = max(rows$elapsed_seconds, na.rm = TRUE),
    peak_memory_max_bytes = max(rows$peak_memory_bytes, na.rm = TRUE),
    warning_total = sum(rows$warning_count, na.rm = TRUE),
    objective_eligible = sum(rows$objective_eligible, na.rm = TRUE),
    practical_converged = sum(rows$practical_converged, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
  pp_write_csv(runtime,
               file.path(output_root, "RUNTIME_AND_MEMORY_SUMMARY.csv"))
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  sizes <- file.info(files)$size
  pp_write_csv(data.frame(
    generated_at_utc = pp_iso_time(), files = length(files),
    disk_bytes = sum(sizes[is.finite(sizes)], na.rm = TRUE),
    data_sets = 3L, endpoints = nrow(endpoints),
    objective_eligible = sum(endpoints$objective_eligible, na.rm = TRUE),
    worker_count = STAGE6G_WORKERS, n_cpus_per_fit = 1L,
    truth_read = FALSE, continuation_started = FALSE,
    formal_paper_mc_result = FALSE, stringsAsFactors = FALSE
  ), file.path(output_root, "RESOURCE_SUMMARY.csv"))
  invisible(TRUE)
}

run_supervisor <- function() {
  require_output_root(FALSE)
  assert_truth_locked()
  verify_environment(require_prepared = TRUE)
  registration <- load_runtime_registration()
  pp_assert(file.exists(runtime_paths()$prepared),
            "Prepare Stage 6G before starting the supervisor.")
  pp_assert(.Platform$OS.type == "unix",
            "The Stage 6G full supervisor requires Unix.")
  log_dir <- file.path(output_root, "supervisor_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  write_status("RUNNING", "120 registered endpoints; no selective rerun")
  run_one <- function(fit_id) {
    key <- gsub("[^A-Za-z0-9_.-]", "_", fit_id)
    stdout <- file.path(log_dir, paste0(key, ".stdout.log"))
    stderr <- file.path(log_dir, paste0(key, ".stderr.log"))
    started <- Sys.time()
    status <- tryCatch(system2(
      "Rscript", c(PP_SCRIPT_FILE, "--action=run-fit",
                   paste0("--fit-id=", fit_id),
                   paste0("--output-root=", output_root)),
      stdout = stdout, stderr = stderr
    ), error = function(condition) 999L)
    ended <- Sys.time()
    pp_write_csv(data.frame(
      fit_id = fit_id, started_utc = pp_iso_time(started),
      ended_utc = pp_iso_time(ended),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      exit_status = as.integer(status), stringsAsFactors = FALSE
    ), file.path(log_dir, paste0(key, ".exit.csv")))
    list(fit_id = fit_id, status = as.integer(status))
  }
  results <- parallel::mclapply(
    registration$fits$fit_id, run_one, mc.cores = STAGE6G_WORKERS,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
  statuses <- vapply(results, function(value) value$status, integer(1L))
  if (any(statuses != 0L)) {
    failed <- vapply(results[statuses != 0L], function(value) value$fit_id,
                     character(1L))
    write_status("INCOMPLETE_PROCESS_FAILURE",
                 paste(failed, collapse = ";"))
    pp_write_lines(c(
      "status=STAGE6G_INCOMPLETE_PROCESS_FAILURE",
      paste0("failed_fit_processes=", paste(failed, collapse = ";")),
      "selective_rerun_started=FALSE", "truth_unseal_authorized=FALSE",
      paste0("ended_utc=", pp_iso_time())
    ), file.path(output_root, "STAGE6G_INCOMPLETE.txt"))
    pp_abort("Stage 6G child-process failures retained: ",
             paste(failed, collapse = ", "))
  }
  tryCatch({
    winners <- freeze_start_winners(registration)
    score_winners(winners)
    freeze_dimensions(registration, winners)
    summarize_resources(registration)
    write_status("COMPLETE_TRUTH_REMAINS_SEALED",
                 "120 endpoints; 30 winners; three selections frozen")
  }, error = function(condition) {
    write_status("INCOMPLETE_FINALIZATION", conditionMessage(condition))
    pp_write_lines(c(
      "status=STAGE6G_INCOMPLETE_FINALIZATION",
      paste0("message=", conditionMessage(condition)),
      "selective_rerun_started=FALSE", "truth_unseal_authorized=FALSE",
      paste0("ended_utc=", pp_iso_time())
    ), file.path(output_root, "STAGE6G_INCOMPLETE.txt"))
    stop(condition)
  })
  cat("STAGE6G_TRUTH_FREE_COMPLETE endpoints=120 winners=30 ",
      "selections=3 truth_remains_sealed=TRUE\n", sep = "")
  invisible(TRUE)
}

status_stage6g <- function() {
  require_output_root(FALSE)
  fit_root <- file.path(output_root, "fits")
  complete <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_COMPLETE[.]txt$", recursive = TRUE
  )) else 0L
  errors <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^FIT_ERROR[.]txt$", recursive = TRUE
  )) else 0L
  terminal <- if (dir.exists(fit_root)) length(list.files(
    fit_root, "^terminal_record[.]rds$", recursive = TRUE
  )) else 0L
  cat("registered_endpoints=120\n")
  cat("terminal_records=", terminal, "\n", sep = "")
  cat("fit_complete=", complete, "\n", sep = "")
  cat("fit_errors=", errors, "\n", sep = "")
  cat("start_selection_frozen=", file.exists(file.path(
    output_root, "truth_free_selection", "START_SELECTION_FROZEN.txt"
  )), "\n", sep = "")
  cat("dimension_selection_frozen=", file.exists(file.path(
    output_root, "truth_free_selection", "DIMENSION_SELECTION_FROZEN.txt"
  )), "\n", sep = "")
  cat("truth_unseal_authorized=", file.exists(file.path(
    output_root, "TRUTH_UNSEAL_AUTHORIZATION.txt"
  )), "\n", sep = "")
  invisible(TRUE)
}

if (identical(action, "check")) {
  verification <- verify_environment(read_truth_registration = TRUE)
  pp_assert(nrow(verification$registration$fits) == 120L,
            "Stage 6G check found the wrong endpoint count.")
  cat("STAGE6G_CHECK_PASS endpoints=120 strata=30 workers=",
      STAGE6G_WORKERS, "\n", sep = "")
} else if (identical(action, "prepare")) {
  prepare_stage6g()
} else if (identical(action, "run-fit")) {
  fit_id <- cli$`fit-id` %||% ""
  pp_assert(nzchar(fit_id), "--fit-id is required for run-fit.")
  run_fit(fit_id)
} else if (identical(action, "run")) {
  run_supervisor()
} else if (identical(action, "status")) {
  status_stage6g()
} else {
  pp_abort("Unknown --action: ", action)
}
