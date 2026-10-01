#!/usr/bin/env Rscript

stage6g_root <- Sys.getenv(
  "STAGE6I_STAGE6G_OUTPUT_ROOT",
  "/root/multifsync-results/stage6g_20260911_v1"
)
stage6h_root <- Sys.getenv(
  "STAGE6I_STAGE6H_OUTPUT_ROOT",
  "/root/multifsync-results/stage6h_20260913_v1"
)
script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]),
  winslash = "/", mustWork = TRUE
)
source(file.path(dirname(script_path), "stage6i_helpers_20260913_v1.R"))

fit_paths <- list.files(
  file.path(stage6h_root, "fits"), pattern = "fit[.]rds$",
  recursive = TRUE, full.names = TRUE
)
paths <- list(
  fit = fit_paths[[1L]],
  observation = file.path(
    stage6g_root, "data", "s6g_minimal_L111_M2", "observation_bundle.rds"
  ),
  truth = file.path(
    stage6g_root, "data", "s6g_minimal_L111_M2", "sealed_truth",
    "truth_bundle.rds"
  )
)
stopifnot(length(fit_paths) >= 1L, all(file.exists(unlist(paths))))

for (name in names(paths)) {
  value <- readRDS(paths[[name]])
  cat("===", name, "===\n")
  cat("path=", paths[[name]], "\n", sep = "")
  cat("class=", paste(class(value), collapse = ";"), "\n", sep = "")
  cat("names=", paste(names(value), collapse = ";"), "\n", sep = "")
  if (name == "fit") {
    cat("fit_C_g_dim=", paste(dim(value$C_g), collapse = "x"), "\n",
        sep = "")
    cat("fit_time_g_length=", length(value$time_g), "\n", sep = "")
    cat("nu_phi_dim=",
        paste(dim(value$mu_q_nu_phi[[1L]]), collapse = "x"),
        "\n", sep = "")
    cat("Sigma_phi_component_dim=",
        paste(dim(value$Sigma_q_nu_phi[[1L]][[1L]]), collapse = "x"),
        "\n", sep = "")
    cat("list_cp_C_class=", paste(class(value$list_cp_C), collapse = ";"),
        " length=", length(value$list_cp_C), "\n", sep = "")
    cat("list_cp_C_study1_length=", length(value$list_cp_C[[1L]]),
        " first_dim=", paste(dim(value$list_cp_C[[1L]][[1L]]), collapse = "x"),
        "\n", sep = "")
    cat("sum_list_cp_C_class=",
        paste(class(value$sum_list_cp_C), collapse = ";"),
        " length=", length(value$sum_list_cp_C), "\n", sep = "")
    cat("sum_list_cp_C_dims=",
        paste(vapply(value$sum_list_cp_C, function(x) {
          paste(dim(x), collapse = "x")
        }, character(1L)), collapse = ";"), "\n", sep = "")
    cat("driver_trace_names=",
        paste(names(value$driver_trace), collapse = ";"), "\n", sep = "")
    cat("driver_trace_rows=", NROW(value$driver_trace), "\n", sep = "")
    if (is.data.frame(value$driver_trace)) {
      print(value$driver_trace[, intersect(
        c("iteration", "temperature", "phase", "block", "study", "factor",
          "component", "posterior_mean_norm", "ppi_mean", "ppi_max",
          "gate_active"), names(value$driver_trace)
      )], row.names = FALSE)
    }
    cat("scale_trace_names=",
        paste(names(value$scale_trace), collapse = ";"), "\n", sep = "")
    cat("scale_trace_rows=", NROW(value$scale_trace), "\n", sep = "")
  } else if (name == "observation") {
    cat("dimension_names=", paste(names(value$dimensions), collapse = ";"),
        "\n", sep = "")
    cat("dimension_values=", paste(unlist(value$dimensions), collapse = ";"),
        "\n", sep = "")
    cat("study1_subjects=", length(value$Y[[1L]]), "\n", sep = "")
    cat("study2_subjects=", length(value$Y[[2L]]), "\n", sep = "")
    cat("Y11_dim=", paste(dim(value$Y[[1L]][[1L]]), collapse = "x"),
        "\n", sep = "")
    cat("time_count_head=",
        paste(vapply(value$time_obs[[1L]][1:5], length, integer(1L)),
              collapse = ";"), "\n", sep = "")
  } else {
    cat("truth_bundle_names=", paste(names(value$data), collapse = ";"),
        "\n", sep = "")
    cat("true_parameter_names=",
        paste(names(value$data$true_params), collapse = ";"), "\n", sep = "")
    cat("true_theta_dim=",
        paste(dim(value$data$true_params$theta_dense_list[[1L]]),
              collapse = "x"), "\n", sep = "")
    for (field in c("b_true", "omega_b_true", "zeta_true", "xi_true",
                    "eta_true", "chi_true", "psi_basis_meta", "M_s",
                    "n_obs_by_subject", "phase_mu_mat")) {
      item <- value$data$true_params[[field]]
      cat("truth_field_", field, "_class=",
          paste(class(item), collapse = ";"), " length=", length(item),
          " dim=", paste(dim(item), collapse = "x"), "\n", sep = "")
    }
    swapped <- s6i_swap_truth_bundle(value)
    restored <- s6i_swap_truth_bundle(swapped)
    cat("truth_double_swap_key_equal=",
        identical(restored$data$true_params$b_true,
                  value$data$true_params$b_true) &&
          identical(restored$data$Y, value$data$Y), "\n", sep = "")
  }
}
