argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
s4c_script_path <- normalizePath(
  sub("^--file=", "", argument[[1L]]), winslash = "/", mustWork = TRUE
)

S4C_UPSTREAM_ROOT <- dirname(s4c_script_path)
source(file.path(
  S4C_UPSTREAM_ROOT, "upstream_v0lv_iguc_common_20260819_v1.R"
), local = FALSE)

UC_EXPERIMENT_ID <- "G12_STOPPING_SEMANTICS_STAGE4C_V1_20260826"
UC_DEV_BRANCH <- "dev/v0lv-init-geometry-20260818"
UC_DEV_COMMIT <- "72d9a53f0ef5e9d3e5f49d9cf837958207b3c980"
UC_DEV_REPO <- uc_env(
  "G12SS4C_DEV_REPO", "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5"
)
UC_DEV_LIB <- uc_env(
  "G12SS4C_DEV_LIB", "/root/multiFSYNC-g12-fixed400-72d9a53f0ef5-lib"
)
UC_DATA_IDS <- c(
  "g12ss4c_base_01", "g12ss4c_base_02",
  "g12ss4c_weak_01", "g12ss4c_weak_02",
  "g12ss4c_sparse_01", "g12ss4c_sparse_02"
)
UC_SCENARIO_IDS <- c(
  "baseline_strong", "weak_loading_separation", "sparse_observation"
)
UC_METHOD_IDS <- "G12"
UC_EXPECTED_FITS <- 12L
UC_OUTER_WORKERS <- 4L
S4C_FIXED_T1 <- 400L
S4C_MAXIT <- 499L
S4C_CHECKPOINTS <- seq.int(20L, 400L, by = 20L)
S4C_STOPPING_MIN_T1 <- 396L
S4C_STOPPING_CONSECUTIVE <- 5L

s4c_verify_environment <- function(
    load_runtime = TRUE,
    package_role = c("development", "frozen_generation")) {
  package_role <- match.arg(package_role)
  uc_assert(dir.exists(UC_DEV_REPO) && dir.exists(UC_FROZEN_REPO),
            "Development source archive or frozen Git checkout is missing.")
  uc_assert(dir.exists(UC_PARENT_FORMAL_ROOT) && dir.exists(UC_BINDING_DIR) &&
              dir.exists(UC_SNAPSHOT),
            "Frozen parent result, binding, or snapshot is missing.")
  uc_assert(dir.exists(UC_DEV_LIB) && dir.exists(UC_FROZEN_LIB),
            "Required R library is missing.")
  binding <- uc_read_csv(file.path(UC_ROOT, "SOURCE_BINDING.csv"))
  required <- c(
    "package_code_commit", "package_source_sha256",
    "public_interface_sha256",
    "generator_source_sha256",
    "calibration_source_sha256", "scale_trace_sha256",
    "multi_core_sha256", "practical_stopping_sha256",
    "initialization_sha256", "g12_stopping_control_sha256",
    "description_sha256", "namespace_sha256",
    "upstream_common_sha256", "package_version"
  )
  uc_assert(nrow(binding) == 1L && all(required %in% names(binding)) &&
              identical(binding$package_code_commit[[1L]], UC_DEV_COMMIT),
            "Stage-4C source binding is incomplete.")
  package_paths <- c(
    package_source_sha256 = "R/diagnostic_driver_trace.R",
    public_interface_sha256 = "R/pre_score_interface.R",
    generator_source_sha256 = "R/generate_data_structured.R",
    scale_trace_sha256 = "R/diagnostic_scale_trace.R",
    multi_core_sha256 = "R/multi_core.R",
    practical_stopping_sha256 = "R/convergence_practical.R",
    initialization_sha256 = "R/initialization.R",
    g12_stopping_control_sha256 = "R/g12_stopping_control.R",
    description_sha256 = "DESCRIPTION",
    namespace_sha256 = "NAMESPACE"
  )
  for (field in names(package_paths)) {
    uc_assert(identical(
      uc_sha256(file.path(UC_DEV_REPO, package_paths[[field]])),
      binding[[field]][[1L]]
    ), paste0("Stage-4C source hash mismatch: ", package_paths[[field]]))
  }
  calibration_path <- uc_find_one(UC_SNAPSHOT, "^v0f_data_calibration[.]R$")
  uc_assert(identical(
    uc_sha256(calibration_path), binding$calibration_source_sha256[[1L]]
  ), "Loading calibration source hash mismatch.")
  uc_assert(identical(
    uc_sha256(file.path(
      UC_ROOT, "upstream_v0lv_iguc_common_20260819_v1.R"
    )), binding$upstream_common_sha256[[1L]]
  ), "Copied upstream common source hash mismatch.")

  selected_library <- if (package_role == "development") UC_DEV_LIB else
    UC_FROZEN_LIB
  .libPaths(unique(c(selected_library, .libPaths())))
  uc_assert(requireNamespace("multiFSYNC", quietly = TRUE),
            paste0(package_role, " multiFSYNC package is unavailable."))
  package_path <- normalizePath(find.package("multiFSYNC"), winslash = "/")
  expected_library <- normalizePath(selected_library, winslash = "/")
  uc_assert(startsWith(package_path, paste0(expected_library, "/")),
            paste0("multiFSYNC was not loaded from ", package_role,
                   " library."))

  if (isTRUE(load_runtime)) {
    source(uc_find_one(UC_SNAPSHOT, "^v0lv_candidate2_runtime[.]R$"),
           local = FALSE)
    source(uc_find_one(
      UC_SNAPSHOT, "^r_route_v3_candidate2_runner[.]R$"
    ), local = FALSE)
    source(uc_find_one(
      UC_SNAPSHOT, "^v0lv_candidate2_controller[.]R$"
    ), local = FALSE)
    if (package_role == "frozen_generation") {
      c2_verify_formal_binding(UC_BINDING_DIR, UC_SNAPSHOT)
    } else {
      status <- uc_read_csv(file.path(UC_BINDING_DIR, "BINDING_STATUS.csv"))
      index <- uc_read_csv(file.path(
        UC_BINDING_DIR, "FULL_RUNTIME_SOURCE_BINDINGS.csv"
      ))
      verification <- c2_verify_hash_index(index, UC_SNAPSHOT)
      uc_assert(nrow(status) == 1L && status$status[[1L]] == "frozen" &&
                  isTRUE(status$formal_execution_authorized[[1L]]) &&
                  all(verification$matches),
                "Frozen parent file-level binding failed.")
    }
  }
  if (package_role == "development") {
    uc_assert(identical(
      as.character(utils::packageVersion("multiFSYNC")),
      binding$package_version[[1L]]
    ), "Installed development package version differs from binding.")
    interface <- getExportedValue("multiFSYNC", "bayesSYNC_multi_pre_score")
    uc_assert(is.function(interface), "Public G12 pre-score interface missing.")
    stopping <- getExportedValue("multiFSYNC", "g12_stopping_control")
    uc_assert(is.function(stopping), "Public G12 stopping helper missing.")
    namespace <- asNamespace("multiFSYNC")
    private <- c(
      ".make_scale_trace_row", ".scale_trace_mean_contributions",
      ".multiFSYNC_scale_trace_option", ".validate_scale_trace_control"
    )
    uc_assert(all(vapply(
      private, exists, logical(1L), envir = namespace, inherits = FALSE
    )), "Private read-only trajectory capture helpers are unavailable.")
  }
  invisible(binding)
}

s4c_validate_data_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "data_id", "data_index", "scenario_id",
    "scenario_replicate", "data_seed", "n_obs_min", "n_obs_max",
    "target_shared_specific_abs_cosine", "new_unseen_at_registration",
    "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_analysis", "truth_unseal_authorized",
    "development_only", "formal_v0lv_result"
  )
  uc_assert(
    is.data.frame(manifest) && nrow(manifest) == 6L &&
      all(required %in% names(manifest)) &&
      all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
      identical(as.character(manifest$data_id), UC_DATA_IDS) &&
      identical(as.integer(manifest$data_index), 1:6) &&
      identical(as.integer(manifest$data_seed), 82726001:82726006) &&
      identical(as.character(manifest$scenario_id),
                rep(UC_SCENARIO_IDS, each = 2L)) &&
      identical(as.integer(manifest$scenario_replicate),
                rep(1:2, times = 3L)) &&
      !anyDuplicated(manifest$data_id) && !anyDuplicated(manifest$data_seed) &&
      all(manifest$new_unseen_at_registration) &&
      !any(manifest$truth_available_to_fit) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_analysis) &&
      !any(manifest$truth_unseal_authorized) &&
      all(manifest$development_only) && !any(manifest$formal_v0lv_result),
    "Stage-4C data manifest identity or truth isolation failed."
  )
  baseline <- manifest$scenario_id == "baseline_strong"
  weak <- manifest$scenario_id == "weak_loading_separation"
  sparse <- manifest$scenario_id == "sparse_observation"
  uc_assert(
    all(manifest$n_obs_min[baseline | weak] == 6L) &&
      all(manifest$n_obs_max[baseline | weak] == 9L) &&
      all(manifest$n_obs_min[sparse] == 3L) &&
      all(manifest$n_obs_max[sparse] == 5L) &&
      all(manifest$target_shared_specific_abs_cosine[baseline | sparse] == 0) &&
      all(abs(manifest$target_shared_specific_abs_cosine[weak] - 0.6) < 1e-15),
    "Stage-4C scenario parameter registration failed."
  )
  invisible(TRUE)
}

s4c_validate_fit_manifest <- function(manifest) {
  required <- c(
    "experiment_id", "fit_id", "data_id", "data_index", "data_seed",
    "scenario_id", "method_id", "function_initialization",
    "calibration_mode", "seed_index", "fit_seed", "n_cpus",
    "initialization", "pre_score_sweeps", "anneal",
    "planned_annealing_sweeps", "fixed_ordinary_T1_sweeps", "total_maxit",
    "stopping_profile", "stopping_min_t1", "stopping_max_t1",
    "stopping_gate", "stopping_consecutive", "checkpoint_sweeps",
    "fixed_400_endpoint",
    "objective_eligibility_required", "outer_parallel_fit_limit",
    "actual_racing", "continuation", "automatic_800",
    "truth_available_to_fit", "truth_available_to_stopping",
    "truth_available_to_analysis", "development_only", "formal_v0lv_result"
  )
  expected_ids <- unlist(lapply(
    UC_DATA_IDS, function(data_id) sprintf("%s__G12__%02d", data_id, 1:2)
  ), use.names = FALSE)
  expected_seed <- 87266000L +
    10L * as.integer(manifest$data_index) + as.integer(manifest$seed_index)
  uc_assert(
    is.data.frame(manifest) && nrow(manifest) == UC_EXPECTED_FITS &&
      all(required %in% names(manifest)) &&
      all(manifest$experiment_id == UC_EXPERIMENT_ID) &&
      identical(as.character(manifest$fit_id), expected_ids) &&
      !anyDuplicated(manifest$fit_id) && !anyDuplicated(manifest$fit_seed) &&
      all(table(manifest$data_id) == 2L) &&
      identical(as.integer(manifest$fit_seed), expected_seed) &&
      all(manifest$method_id == "G12") &&
      all(manifest$function_initialization == "gram_unit_energy") &&
      all(manifest$calibration_mode == "function") &&
      all(manifest$n_cpus == 1L) &&
      all(manifest$initialization == "random") &&
      all(manifest$pre_score_sweeps == 1L) &&
      all(manifest$anneal == "c(1,1.9,100)") &&
      all(manifest$planned_annealing_sweeps == 99L) &&
      all(manifest$fixed_ordinary_T1_sweeps == S4C_FIXED_T1) &&
      all(manifest$total_maxit == S4C_MAXIT) &&
      all(manifest$stopping_profile == "g12_stopping_control") &&
      all(manifest$stopping_min_t1 == S4C_STOPPING_MIN_T1) &&
      all(manifest$stopping_max_t1 == S4C_FIXED_T1) &&
      all(manifest$stopping_gate == "quantile_factor") &&
      all(manifest$stopping_consecutive == S4C_STOPPING_CONSECUTIVE) &&
      all(manifest$checkpoint_sweeps ==
            paste(S4C_CHECKPOINTS, collapse = ";")) &&
      all(manifest$fixed_400_endpoint) &&
      all(manifest$objective_eligibility_required) &&
      all(manifest$outer_parallel_fit_limit == UC_OUTER_WORKERS) &&
      !any(manifest$actual_racing) && !any(manifest$continuation) &&
      !any(manifest$automatic_800) &&
      !any(manifest$truth_available_to_fit) &&
      !any(manifest$truth_available_to_stopping) &&
      !any(manifest$truth_available_to_analysis) &&
      all(manifest$development_only) && !any(manifest$formal_v0lv_result),
    "Stage-4C fit route, budget, seed, or isolation contract failed."
  )
  invisible(TRUE)
}

s4c_base_config <- function(data_row, smoke = FALSE) {
  scenario <- as.character(data_row$scenario_id[[1L]])
  if (smoke) {
    return(list(
      S = 2L, n_s = c(4L, 4L), p = 20L, d = 0L,
      L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
      M_s = list(2L, 2L), K = 3L, n_obs = 5L,
      common_grid = FALSE,
      n_obs_range = if (scenario == "sparse_observation") c(3L, 4L)
      else c(4L, 6L),
      sigma_eps = 0.3, bool_sparse_loadings = TRUE,
      prop_sparse = 0.8, score_var_decay = TRUE,
      n_dense = 51L, mean_amp = 0.6, beta_amp = 0,
      bs_degree = c(2L, 3L), identified_loadings = TRUE,
      sparsity_mode = "fixed"
    ))
  }
  list(
    S = 2L, n_s = c(30L, 30L), p = 500L, d = 0L,
    L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
    M_s = list(2L, 2L), K = 5L, n_obs = 8L,
    common_grid = FALSE,
    n_obs_range = if (scenario == "sparse_observation") c(3L, 5L)
    else c(6L, 9L),
    sigma_eps = 0.3, bool_sparse_loadings = TRUE,
    prop_sparse = 0.9, score_var_decay = TRUE,
    n_dense = 201L, mean_amp = 0.6, beta_amp = 0,
    bs_degree = c(2L, 3L), identified_loadings = TRUE,
    sparsity_mode = "fixed"
  )
}

s4c_abs_cosine <- function(left, right) {
  denominator <- sqrt(sum(left^2) * sum(right^2))
  if (!is.finite(denominator) || denominator <= 1e-14) return(NA_real_)
  abs(sum(left * right) / denominator)
}

s4c_apply_weak_loading_transform <- function(data, rho = 0.6) {
  truth <- data$true_params
  uc_assert(truth$L_f == 1L && all(truth$L_s_by_study == 1L) &&
              isTRUE(truth$identified_loadings),
            "Weak-loading transform requires identified 1/1/1 truth.")
  a <- as.numeric(truth$a_true[, 1L])
  a_support <- which(a != 0)
  a_unit <- a / sqrt(sum(a^2))
  observed <- numeric(truth$S)
  support_overlap <- numeric(truth$S)
  for (study in seq_len(truth$S)) {
    old_b <- as.numeric(truth$b_true[[study]][, 1L])
    old_values <- old_b[old_b != 0]
    uc_assert(length(old_values) == length(a_support),
              "Weak-loading transform requires equal active counts.")
    raw <- numeric(truth$p)
    raw[a_support] <- old_values
    orthogonal <- raw - a_unit * sum(a_unit * raw)
    orthogonal_norm <- sqrt(sum(orthogonal^2))
    if (!is.finite(orthogonal_norm) || orthogonal_norm <= 1e-12) {
      raw[a_support] <- c(tail(old_values, -1L), old_values[[1L]])
      orthogonal <- raw - a_unit * sum(a_unit * raw)
      orthogonal_norm <- sqrt(sum(orthogonal^2))
    }
    uc_assert(is.finite(orthogonal_norm) && orthogonal_norm > 1e-12,
              "Could not construct a nondegenerate weak-loading direction.")
    new_b <- sqrt(sum(old_b^2)) *
      (rho * a_unit + sqrt(1 - rho^2) * orthogonal / orthogonal_norm)
    new_b[-a_support] <- 0
    truth$b_true[[study]][, 1L] <- new_b
    observed[[study]] <- s4c_abs_cosine(a, new_b)
    support_overlap[[study]] <- length(intersect(
      a_support, which(new_b != 0)
    )) / length(a_support)
  }
  truth$gamma_b_true <- lapply(truth$b_true, function(value) 1L * (value != 0))
  truth$omega_b_true <- lapply(truth$gamma_b_true, colMeans)
  truth$identified_loadings <- FALSE
  truth$s4c_loading_scenario <- "weak_loading_separation"
  truth$s4c_target_shared_specific_abs_cosine <- rho
  truth$s4c_observed_shared_specific_abs_cosine <- observed
  truth$s4c_shared_specific_support_overlap_fraction <- support_overlap
  data$true_params <- truth
  v0f_rebuild_observations(data)
}

s4c_truth_assertions <- function(data, data_row, config) {
  truth <- data$true_params
  scenario <- as.character(data_row$scenario_id[[1L]])
  checks <- list()
  add <- function(id, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      assertion_id = id, passed = isTRUE(passed), detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  add("dimensions", identical(
    as.integer(c(truth$S, truth$n_s, truth$p, truth$d_use)),
    as.integer(c(config$S, config$n_s, config$p, config$d))
  ), "S,n_s,p,d match scenario config")
  add("factor_and_fpca_counts",
      truth$L_f == 1L && all(truth$L_s_by_study == 1L) &&
        identical(as.integer(truth$M_f), 2L) &&
        all(vapply(truth$M_s, function(value) identical(
          as.integer(value), 2L
        ), logical(1L))),
      "truth is 1/1/1 with M=2")
  observed_counts <- unlist(lapply(data$time_obs, lengths), use.names = FALSE)
  add("observation_counts",
      all(observed_counts >= config$n_obs_range[[1L]] &
            observed_counts <= config$n_obs_range[[2L]]),
      paste0("range=", paste(range(observed_counts), collapse = "-")))
  a <- as.numeric(truth$a_true[, 1L])
  b <- lapply(truth$b_true, function(value) as.numeric(value[, 1L]))
  active <- c(sum(a != 0), vapply(b, function(value) sum(value != 0), integer(1L)))
  target_active <- as.integer(round(0.1 * config$p))
  if (config$p == 20L) target_active <- as.integer(round(0.2 * config$p))
  add("active_loading_counts", all(active == target_active),
      paste0("active=", paste(active, collapse = ",")))
  norms <- c(sqrt(sum(a^2)),
             vapply(b, function(value) sqrt(sum(value^2)), numeric(1L)))
  add("loading_norms", max(abs(norms - norms[[1L]])) <= 1e-10,
      paste0("norms=", paste(format(norms, digits = 10), collapse = ",")))
  cosines <- vapply(b, function(value) s4c_abs_cosine(a, value), numeric(1L))
  if (scenario == "weak_loading_separation") {
    overlap <- vapply(b, function(value) {
      length(intersect(which(a != 0), which(value != 0))) / sum(a != 0)
    }, numeric(1L))
    add("weak_loading_cosine", max(abs(cosines - 0.6)) <= 1e-10,
        paste0("cosines=", paste(format(cosines, digits = 10), collapse = ",")))
    add("weak_loading_support_overlap", min(overlap) == 1,
        paste0("overlap=", paste(overlap, collapse = ",")))
  } else {
    disjoint <- vapply(b, function(value) {
      length(intersect(which(a != 0), which(value != 0))) == 0L
    }, logical(1L))
    add("strong_loading_separation", all(disjoint) && max(cosines) <= 1e-12,
        paste0("cosines=", paste(format(cosines, digits = 5), collapse = ",")))
  }
  reconstruction_error <- 0
  observation_error <- 0
  for (study in seq_len(truth$S)) {
    for (subject in seq_len(truth$n_s[[study]])) {
      shared <- truth$f_true_values[[study]][[subject]] %*% t(truth$a_true)
      specific <- truth$g_true_values[[study]][[subject]] %*%
        t(truth$b_true[[study]])
      signal <- truth$mu_true_values[[study]][[subject]] +
        truth$beta_true_values[[study]][[subject]] + shared + specific
      reconstruction_error <- max(
        reconstruction_error,
        abs(signal - truth$signal_true_values[[study]][[subject]])
      )
      observed <- do.call(cbind, data$Y[[study]][[subject]])
      observation_error <- max(
        observation_error,
        abs(observed - signal - truth$noise_true_values[[study]][[subject]])
      )
    }
  }
  add("signal_and_noise_reconstruction",
      max(reconstruction_error, observation_error) <= 1e-12,
      paste0("signal=", reconstruction_error, ";Y=", observation_error))
  result <- do.call(rbind, checks)
  rownames(result) <- NULL
  result
}

s4c_copy_matrix <- function(value) {
  value <- as.matrix(value)
  matrix(as.numeric(value), nrow = nrow(value), ncol = ncol(value))
}

s4c_score_kernel <- function(means, covariances) {
  means <- s4c_copy_matrix(means)
  result <- means %*% t(means)
  if (length(covariances)) {
    uc_assert(length(covariances) == nrow(means),
              "Score covariance list length differs from subjects.")
    diag(result) <- diag(result) + vapply(
      covariances, function(value) sum(diag(as.matrix(value))), numeric(1L)
    )
  }
  result
}

s4c_make_direction_snapshot <- function(cumulative_sweep, C_g, state,
                                         contributions) {
  list(
    cumulative_sweep = as.integer(cumulative_sweep),
    shared_loading = s4c_copy_matrix(state$mu_q_a),
    specific_loading = lapply(state$mu_q_b_specific, s4c_copy_matrix),
    shared_feature = lapply(
      state$mu_q_nu_phi,
      function(value) s4c_copy_matrix(C_g %*% as.matrix(value))
    ),
    specific_feature = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_nu_psi[[study]], function(value) {
        s4c_copy_matrix(C_g %*% as.matrix(value))
      })
    }),
    shared_score_mean = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_zeta[[study]], s4c_copy_matrix)
    }),
    specific_score_mean = lapply(seq_len(state$S), function(study) {
      lapply(state$mu_q_xi[[study]], s4c_copy_matrix)
    }),
    shared_score_kernel = lapply(seq_len(state$S), function(study) {
      Map(s4c_score_kernel, state$mu_q_zeta[[study]],
          state$Sigma_q_zeta[[study]])
    }),
    specific_score_kernel = lapply(seq_len(state$S), function(study) {
      Map(s4c_score_kernel, state$mu_q_xi[[study]],
          state$Sigma_q_xi[[study]])
    }),
    shared_contribution = as.numeric(contributions$shared),
    specific_contribution = as.numeric(contributions$specific)
  )
}

s4c_trace_control <- function() {
  list(
    include_initial = TRUE, include_annealing = FALSE,
    t1_sweeps = S4C_CHECKPOINTS, block_annealing_sweeps = integer()
  )
}

s4c_with_private_traces <- function(trace_control, code) {
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
  trace_control <- validator(trace_control)
  capture <- new.env(parent = emptyenv())
  capture$snapshots <- list()
  wrapper <- function(...) {
    arguments <- list(...)
    result <- do.call(original, arguments)
    stage <- as.character(arguments$stage)
    sweep <- as.integer(arguments$t1_sweep)
    should_capture <- identical(stage, "initial") ||
      (identical(stage, "t1") && sweep %in% trace_control$t1_sweeps)
    if (should_capture) {
      cumulative <- if (identical(stage, "initial")) 0L else sweep
      contributions <- contribution_function(
        Y = arguments$Y, C = arguments$C, Z = arguments$Z,
        state = arguments$state
      )
      capture$snapshots[[as.character(cumulative)]] <-
        s4c_make_direction_snapshot(
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

s4c_fit_control <- function() {
  spec <- r_route_v3_candidate2_spec("formal")
  spec$t1_budget <- S4C_FIXED_T1
  spec$maxit <- S4C_MAXIT
  spec$practical_control <- multiFSYNC::g12_stopping_control()
  spec$practical_control$checkpoints <- S4C_CHECKPOINTS
  spec
}

s4c_terminal_summary <- function(record) {
  data.frame(
    fit_id = record$fit_id, data_id = record$data_id,
    scenario_id = record$scenario_id, seed_index = record$seed_index,
    fit_seed = record$fit_seed, terminal_status = record$terminal_status,
    actual_annealing_sweeps = record$actual_annealing_sweeps,
    actual_T1_sweeps = record$actual_T1_sweeps,
    objective_eligible = record$objective_eligible,
    final_elbo = record$final_elbo,
    practical_converged = record$practical_converged,
    slow_case = record$slow_case,
    convergence_status = record$convergence_status,
    convergence_reason = record$convergence_reason,
    warning_count = if (is.data.frame(record$warnings))
      nrow(record$warnings) else 0L,
    elapsed_seconds = record$elapsed_seconds,
    peak_memory_bytes = record$peak_memory_bytes,
    direction_snapshot_count = record$direction_snapshot_count,
    truth_used = record$truth_used,
    continuation_used = record$continuation_used,
    automatic_800_used = record$automatic_800_used,
    formal_v0lv_result = record$formal_v0lv_result,
    stringsAsFactors = FALSE
  )
}
