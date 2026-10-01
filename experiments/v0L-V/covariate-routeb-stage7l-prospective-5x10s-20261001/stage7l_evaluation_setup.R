evaluation_setup <- function(exec_dir,lib,support,load=TRUE) {
  parent <- file.path(support,"experiments","v0L-V")
  stage6a <- file.path(parent,"g12-dimension-overspec-stage6a-20260831","stage6a_common_20260831_v1.R")
  stage6h <- file.path(parent,"g12-reachability-causal-audit-stage6h-20260913","stage6h_evaluation_helpers_20260913_v1.R")
  check_hash(stage6a,"bc8e692899c792d440f1e470648f6be7dc0feafb6aab6b4fbe89b630e1841ae6")
  check_hash(stage6h,"5d6327289579a0624b058aafe6d6db5164b377e89d4efb18fe26abc1bb57f221")
  snapshot <- file.path(parent,"rc1_snapshot")
  expected <- c(
    v0lv_candidate2_runtime.R="6f0eb7d9c1c01ec3d83074eee43b167cec96c53e7d681ed4ec852f4f791f76bb",
    r_route_v3_candidate2_runner.R="e7d5808bfdd6af756c840eb84b02d87fc3883823742bd93675d1d4d5451003b6",
    v0lv_candidate2_controller.R="ae88153e32bf15455e03e98ade5f64a91a06a47c2947adc38f3da2472a22c5ab",
    v0lv_candidate2_evaluation.R="5ef50d51572e946e158a6b67ee28601c4b8d8847f9830fa11d881e682df24d1e",
    v0g_evaluation.R="6b765191657f88db46a027a95acd3917c602b35cb6cee8947d917652f9c6b271",
    v0h_evaluation.R="45639c2d9f14e4df0187d519a62837ef3b3be090ed315123cd97c725f8576dd4",
    v0i_truth_diagnostic.R="127615551ef1e4dbe14b8c486e3e3ad27ba45eea77f764c7ed7c2fba8318b194",
    v0j_evaluation.R="3c91af1e07b0ba1aac4539a2c1e28d8abcfe2f5ee239a2a2daed7fcb2714a833")
  paths <- vapply(names(expected),function(nm) {
    hits <- list.files(snapshot,pattern=paste0("^",nm,"$"),recursive=TRUE,full.names=TRUE)
    assert(length(hits)==1L,paste("Ambiguous frozen evaluator",nm))
    check_hash(hits,expected[[nm]]);hits
  },character(1))
  check_hash(file.path(exec_dir,"dependencies","stage7c_evaluator.R"),
             "1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8")
  if (!load) {
    paths <- c(stage6a=stage6a,stage6h=stage6h,paths)
    new_csv(data.frame(role=names(paths),path=unname(paths),sha256=vapply(paths,sha,character(1))),
            file.path(exec_dir,"EVALUATOR_SOURCE_BINDING.csv"))
    return(invisible(TRUE))
  }
  Sys.setenv(STAGE6A_R_LIB=lib,STAGE6A_USE_LOAD_ALL="false")
  source(stage6a,local=.GlobalEnv);source(stage6h,local=.GlobalEnv)
  assign("PP_SNAPSHOT",snapshot,.GlobalEnv)
  pp_load_package();pp_source_runtime(TRUE)
  env <- candidate2_load_frozen_multi_evaluator(snapshot)
  extract_functions(file.path(exec_dir,"dependencies","stage7c_evaluator.R"),
    c("trapezoid_weights","interp","mean_design","true_mean_dense",
      "estimated_mean_dense","mean_recovery","install_dense_residual_override",
      "new_error_accumulator","accumulate_error","finish_rmse","finish_nrmse",
      "new_projection_accumulator","accumulate_projection","finish_projection",
      "covariate_and_signal_metrics","add_recovery_flags"))
  install_dense_residual_override(env)
  env
}
