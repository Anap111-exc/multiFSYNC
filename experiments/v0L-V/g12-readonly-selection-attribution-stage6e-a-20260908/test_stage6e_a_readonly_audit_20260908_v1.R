script_filename <- "test_stage6e_a_readonly_audit_20260908_v1.R"
if (file.exists(file.path(getwd(), script_filename))) {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                               winslash = "/", mustWork = TRUE)
  script_dir <- dirname(script_path)
}
source(file.path(script_dir, "stage6e_a_tools_20260908_v1.R"),
       local = FALSE)

results <- list()
run_test <- function(name, expression) {
  error <- NULL
  tryCatch(force(expression), error = function(condition) {
    error <<- conditionMessage(condition)
  })
  results[[length(results) + 1L]] <<- data.frame(
    test = name, passed = is.null(error), detail = if (is.null(error)) "PASS" else error,
    stringsAsFactors = FALSE
  )
}

run_test("stratified_se_matches_manual", {
  loss <- c(1, 2, 3, 4, 2, 4, 6, 8)
  study <- rep(1:2, each = 4)
  expected <- sqrt(sum((c(0.5, 0.5)^2) *
                         c(stats::var(loss[1:4]), stats::var(loss[5:8])) / 4))
  stopifnot(abs(s6ea_stratified_mean_se(loss, study) - expected) < 1e-12)
})

run_test("top_k_unit_median", {
  rows <- do.call(rbind, lapply(1:5, function(rank) data.frame(
    fit_id = paste0("f", rank), data_id = "d", fit_config_id = "c",
    seed_index = rank, final_elbo = -rank, elbo_rank = rank,
    study = rep(1:2, each = 2), subject = rep(1:2, 2),
    holdout_time = 0.5, mse = rank, mae = rank,
    observed_second_moment = 4, stringsAsFactors = FALSE
  )))
  result <- s6ea_make_risk_units(rows)
  stopifnot(nrow(result) == 12L,
            all(result$mse[result$risk_method == "max_elbo_start"] == 1),
            all(result$mse[result$risk_method == "top3_elbo_median"] == 2),
            all(result$mse[result$risk_method == "top5_elbo_median"] == 3))
})

run_test("paired_rule_uses_difference_se", {
  make_rows <- function(config, complexity, base, delta) data.frame(
    risk_method = "max_elbo_start", top_k = 1L, data_id = "d",
    fit_config_id = config, study = rep(1:2, each = 4),
    subject = rep(1:4, 2), holdout_time = 0.5, source_fit_count = 1L,
    mse = base + delta, mae = sqrt(base + delta),
    observed_second_moment = 4, stringsAsFactors = FALSE
  )
  base <- c(1, 2, 4, 8, 1, 2, 4, 8)
  rows <- rbind(
    make_rows("factor_L1_M2", 1, base, 0.20),
    make_rows("factor_L2_M2", 2, base, 0.02),
    make_rows("factor_L3_M2", 3, base, 0),
    make_rows("fpca_L1_M3", 3, base, 0.30),
    make_rows("fpca_L1_M4", 4, base, 0.40)
  )
  compared <- s6ea_compare_selection_rules(rows)
  paired <- compared$selections[
    compared$selections$one_se_semantics == "paired_difference_one_se" &
      compared$selections$axis == "factor", ]
  stopifnot(nrow(paired) == 1L,
            paired$minimum_loss_config_id == "factor_L3_M2",
            paired$selected_config_id != "factor_L1_M2")
})

result <- do.call(rbind, results)
print(result)
utils::write.csv(result, file.path(script_dir, "STAGE6E_A_TEST_RESULTS.csv"),
                 row.names = FALSE, na = "")
if (!all(result$passed)) quit(status = 1L)
