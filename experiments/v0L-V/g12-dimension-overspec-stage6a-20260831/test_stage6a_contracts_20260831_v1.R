#!/usr/bin/env Rscript

argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(argument) != 1L) stop("Run with Rscript.", call. = FALSE)
script_file <- sub("^--file=", "", argument[[1L]])
script_root <- dirname(script_file)
if (!dir.exists(script_root)) {
  script_root <- dirname(normalizePath(script_file, winslash = "/",
                                      mustWork = TRUE))
}
source(file.path(script_root, "stage6a_common_20260831_v1.R"), local = FALSE)

manifests <- pp_load_manifests()
stopifnot(
  nrow(manifests$data) == 1L,
  nrow(manifests$fits) == 36L,
  sum(manifests$fits$capacity_probe) == 3L,
  all(table(manifests$fits$data_id, manifests$fits$fit_config_id) == 12L)
)

dimensions <- lapply(PP_FIT_CONFIG_IDS, function(config_id) {
  row <- manifests$fits[manifests$fits$fit_config_id == config_id,
                        , drop = FALSE][1L, , drop = FALSE]
  pp_fit_dimensions(row)
})
names(dimensions) <- PP_FIT_CONFIG_IDS
stopifnot(
  identical(dimensions$truth_L1_M2,
            list(L_f = 1L, L_s = c(1L, 1L), M_f = 2L,
                 M_s = list(2L, 2L))),
  identical(dimensions$factor_L3_M2,
            list(L_f = 3L, L_s = c(3L, 3L), M_f = c(2L, 2L, 2L),
                 M_s = list(c(2L, 2L, 2L), c(2L, 2L, 2L)))),
  identical(dimensions$fpca_L1_M4,
            list(L_f = 1L, L_s = c(1L, 1L), M_f = 4L,
                 M_s = list(4L, 4L)))
)

synthetic <- manifests$fits[
  manifests$fits$data_id == PP_DATA_IDS[[1L]] &
    manifests$fits$seed_index <= 2L, , drop = FALSE
]
synthetic$terminal_status <- "fixed_400_complete"
synthetic$objective_eligible <- TRUE
synthetic$truth_used_for_fit_stopping_or_selection <- FALSE
synthetic$final_elbo <- ifelse(synthetic$seed_index == 2L, 10, 1)
winners <- pp_select_g12_winners(synthetic)
stopifnot(
  nrow(winners) == 3L,
  all(winners$seed_index == 2L),
  all(table(winners$selection_stratum_id) == 1L),
  all(winners$selection_rule ==
        "maximum_eligible_ordinary_T1_ELBO_within_data_and_fit_config")
)

cat("STAGE6A_CONTRACT_TEST_PASS data=1 configs=3 fits=36 winners_per_data=3\n")
