#!/usr/bin/env Rscript

tools_path <- "stage6c_basin_tools_20260904_v1.R"
if (!file.exists(tools_path)) {
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                          value = TRUE)[[1L]])
  tools_path <- file.path(dirname(script_path), basename(tools_path))
}
source(tools_path, local = FALSE)

set.seed(61004)
loading <- matrix(rnorm(30), 10, 3)
paths <- matrix(rnorm(60), 3, 20)
rotation <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
left <- list(loadings = loading, paths = paths)
right <- list(loadings = loading %*% rotation,
              paths = t(rotation) %*% paths)
pair <- s6c_contribution_pair(left, right)
stopifnot(abs(pair[["contribution_cosine"]] - 1) < 1e-10,
          pair[["contribution_relative_difference"]] < 1e-8)
subspace <- s6c_subspace_similarity(left$loadings, right$loadings)
stopifnot(abs(subspace[["mean_squared_canonical_correlation"]] - 1) < 1e-10)

zero <- matrix(0, 10, 2)
stopifnot(s6c_subspace_similarity(zero, zero)[[1L]] == 1,
          s6c_subspace_similarity(zero, loading)[[1L]] == 0)

endpoints <- data.frame(
  cohort_id = "test", config_group_id = "test::config", fit_id = c("a", "b", "c"),
  fit_config_id = "config", seed_index = 1:3, fit_seed = 11:13,
  final_elbo = c(1, 3, 2), practical_converged = TRUE,
  fit_L_f = 1L, fit_L_s_1 = 1L, fit_L_s_2 = 1L,
  ppi_active_shared = 1L, ppi_active_specific_1 = 1L,
  ppi_active_specific_2 = 1L, maximum_effective_M_shared = 2L,
  maximum_effective_M_specific_1 = 2L,
  maximum_effective_M_specific_2 = 2L, truth_used = FALSE,
  stringsAsFactors = FALSE
)
pairs <- data.frame(
  fit_id_left = c("a", "a", "b"), fit_id_right = c("b", "c", "c"),
  contribution_similarity_min = c(0.95, 0.2, 0.2),
  structural_similarity_min = c(0.94, 0.1, 0.1),
  stringsAsFactors = FALSE
)
assignment <- s6c_cluster_group(endpoints, pairs, "contribution", 0.9)
stopifnot(assignment$basin_size[assignment$fit_id == "a"] == 2L,
          assignment$basin_size[assignment$fit_id == "b"] == 2L,
          assignment$basin_size[assignment$fit_id == "c"] == 1L,
          unique(assignment$basin_representative_fit_id[
            assignment$basin_size == 2L]) == "b")

freeze_root <- tempfile("stage6c-freeze-")
dir.create(freeze_root)
s6c_write_lines("frozen", file.path(freeze_root, "result.txt"))
s6c_write_csv(data.frame(
  relative_path = "result.txt",
  sha256 = s6c_sha256(file.path(freeze_root, "result.txt")),
  stringsAsFactors = FALSE
), file.path(freeze_root, "TRUTH_FREE_OUTPUT_HASHES.csv"))
s6c_write_lines("status=complete",
                file.path(freeze_root, "TRUTH_FREE_AUDIT_COMPLETE.txt"))
stopifnot(isTRUE(s6c_verify_truth_free_freeze(freeze_root)))

cat("STAGE6C_TARGETED_TESTS_PASS tests=8\n")
