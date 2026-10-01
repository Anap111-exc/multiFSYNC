#!/usr/bin/env Rscript
args <- commandArgs(TRUE);stopifnot(length(args)==4L)
exec_dir<-args[1L];output<-args[2L];lib<-args[3L];support<-args[4L]
source(file.path(exec_dir,"stage7l_common.R"));load_stage7l(exec_dir,lib)
freeze<-file.path(output,"blind_freeze")
marker<-jsonlite::fromJSON(file.path(freeze,"BLIND_ENDPOINTS_FROZEN.txt"))
assert(identical(marker$status,"STAGE7L_BLIND_ENDPOINTS_FROZEN") &&
       isFALSE(marker$truth_opened) && marker$endpoints==55L,"Blind freeze absent.")
for(nm in names(marker$hashes)) check_hash(file.path(freeze,nm),marker$hashes[[nm]])
hashes<-read_csv(file.path(freeze,"BLIND_FILE_HASHES.csv"))
for(i in seq_len(nrow(hashes))) check_hash(hashes$path[i],hashes$sha256[i])
source(file.path(exec_dir,"stage7l_evaluation_setup.R"))
evaluator_environment<-evaluation_setup(exec_dir,lib,support)
endpoints<-read_csv(file.path(freeze,"ALL_ENDPOINTS_BLIND.csv"))
selection<-read_csv(file.path(freeze,"PAIRED_SELECTION_BLIND.csv"))
root<-file.path(output,"scientific_evaluation")
assert(!dir.exists(root),"Refusing evaluation overwrite.");dir.create(root,mode="0700")
new_lines(c("POST_BLIND_FREEZE_TRUTH_UNSEAL_AUTHORIZED","authorization_date=2026-10-01",
 "truth_used_for_fit_stopping_or_selection=FALSE",paste0("freeze_sha256=",sha(
 file.path(freeze,"BLIND_ENDPOINTS_FROZEN.txt")))),file.path(root,"TRUTH_UNSEAL_AUTHORIZATION.txt"))
evaluate_one<-function(i) {
  row<-endpoints[i,,drop=FALSE]; target<-file.path(root,row$fit_id)
  dir.create(target,mode="0700");warnings<-character();begin<-proc.time()[["elapsed"]]
  if(row$status!="COMPLETE") return(list(status=data.frame(fit_id=row$fit_id,
    evaluated=FALSE,evaluation_error=paste("No complete endpoint:",row$status)),compact=NULL,role=NULL))
  set.seed(as.integer(row$fit_seed))
  result<-tryCatch(withCallingHandlers({
    fit<-readRDS(row$fit_path)
    assert(fit$i_iter==499L && length(fit$ELBO)==400L &&
      abs(tail(fit$ELBO,1L)-row$final_ordinary_elbo)<1e-8,"Frozen endpoint differs.")
    sealed<-readRDS(file.path(output,"data",row$data_id,"sealed_truth","truth_bundle.rds"))
    assert(identical(sealed$experiment_id,experiment_id) && identical(sealed$data_id,row$data_id) &&
      identical(sealed$data$true_params$truth_basis_method,"iid_common_splinet_stiefel") &&
      identical(sealed$data$true_params$beta_truth_method,
                "independent_cubic_bspline_d2_matched_mean_centered_energy"),"Truth identity changed.")
    historical<-s6h_evaluate_fit(fit,sealed,evaluator_environment)
    mean<-mean_recovery(fit,sealed$data)
    covariate<-covariate_and_signal_metrics(fit,sealed,evaluator_environment,"C1")
    for(nm in names(historical)) if(is.data.frame(historical[[nm]]))
      new_csv(historical[[nm]],file.path(target,paste0(nm,".csv")))
    new_csv(covariate$beta_errors,file.path(target,"BETA_COEFFICIENT_ERRORS.csv"))
    new_csv(mean$by_study,file.path(target,"MEAN_BY_STUDY.csv"))
    new_rds(list(historical=historical,mean=mean,covariate=covariate),file.path(target,"evaluation.rds"))
    legacy<-historical$summary
    legacy<-legacy[,!names(legacy)%in%c("signal_nrmse_all_candidates",
      "signal_nrmse_ppi_selected","signal_nrmse_mean_only"),drop=FALSE]
    names(legacy)[grepl("^dense_signal_",names(legacy))]<-paste0("legacy_residual_",
      names(legacy)[grepl("^dense_signal_",names(legacy))])
    activity<-historical$candidate_activity_summary
    names(activity)<-paste0("activity_",names(activity))
    base<-data.frame(fit_id=row$fit_id,data_id=row$data_id,branch=row$branch,
      logical_data_id=sub("__C1$","",row$data_id),beta_arm="C1",
      truth_L_f=1L,truth_L_s_1=1L,truth_L_s_2=1L,truth_M_common=2L)
    science<-add_recovery_flags(cbind(base,covariate$summary,legacy,
      historical$process_contribution_summary,historical$loading_scale_summary,mean$summary,activity))
    science$full_structural_recovery<-with(science,core_truth_recovered &
      exact_thresholded_factor_structure & full_true_feature_component_coverage)
    roles<-historical$process_contribution_by_role
    roles$fit_id<-row$fit_id;roles$data_id<-row$data_id;roles$branch<-row$branch
    bad<-roles[roles$loading_structure_status!="correct",,drop=FALSE]
    science$only_B2_loading_failure<-nrow(bad)>0 &&
      all(bad$study==2L & bad$true_role=="specific")
    science$other_block_loading_failure<-nrow(bad)>0 &&
      any(bad$study!=2L | bad$true_role!="specific")
    new_csv(science,file.path(target,"science.csv"))
    list(compact=science,role=roles)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
   error=function(e)e)
  new_lines(warnings,file.path(target,"EVALUATION_WARNINGS.txt"))
  okay<-!inherits(result,"error")
  list(status=data.frame(fit_id=row$fit_id,evaluated=okay,
    evaluation_error=if(okay) "" else conditionMessage(result),
    evaluation_seconds=proc.time()[["elapsed"]]-begin),
    compact=if(okay)result$compact else NULL,role=if(okay)result$role else NULL)
}
results<-parallel::mclapply(seq_len(nrow(endpoints)),evaluate_one,
                          mc.cores=8L,mc.preschedule=FALSE,mc.set.seed=FALSE)
assert(!any(vapply(results,inherits,logical(1),"try-error")),"Evaluation worker crashed; evidence retained.")
status<-pp_rbind_fill(lapply(results,`[[`,"status"))
science<-pp_rbind_fill(lapply(results,`[[`,"compact"))
roles<-pp_rbind_fill(lapply(results,`[[`,"role"))
new_csv(status,file.path(root,"EVALUATION_STATUS.csv"))
all<-merge(endpoints,science,by=c("fit_id","data_id","branch"),all.x=TRUE,sort=FALSE)
all<-merge(all,status,by="fit_id",all.x=TRUE,sort=FALSE)
all<-all[match(endpoints$fit_id,all$fit_id),,drop=FALSE]
new_csv(all,file.path(root,"ALL_ENDPOINT_SCIENCE.csv"))
new_csv(roles,file.path(root,"ALL_ROLE_SCIENCE.csv"))
paired<-do.call(rbind,lapply(seq_len(nrow(selection)),function(i) {
  sel<-selection[i,,drop=FALSE];pool<-all[all$data_id==sel$data_id & all$branch=="U",,drop=FALSE]
  result<-sel
  result$U_pool_core_count<-sum(pool$core_truth_recovered%in%TRUE)
  result$U_pool_exact_count<-sum(pool$exact_thresholded_factor_structure%in%TRUE)
  result$U_pool_full_recovery_count<-sum(pool$full_structural_recovery%in%TRUE)
  result$U_pool_ever_correct<-any(pool$full_structural_recovery%in%TRUE)
  result$U_cpu_seconds<-sum(pool$elapsed_seconds,na.rm=TRUE)
  indicators<-c("core_truth_recovered","exact_thresholded_factor_structure",
    "full_true_feature_component_coverage","full_structural_recovery","objective_eligible",
    "strict_practical_converged","final_ordinary_elbo","elapsed_seconds","status","evaluated",
    "only_B2_loading_failure","other_block_loading_failure","R_A","R_B1","R_B2",
    "global_missing","global_misplaced","global_extra","global_duplicate",
    "dense_full_signal_nrmse_selected","complete_contribution_nrmse_mean_missing_as_zero",
    "feature_total_ise_mean","raw_loading_norm_ratio_mean","canonical_loading_norm_ratio_mean")
  for(role in c("U","S","enhanced")) {
    id<-sel[[switch(role,U="U_winner_fit_id",S="S_fit_id",enhanced="enhanced_winner_fit_id")]]
    row<-all[all$fit_id==id,,drop=FALSE]
    for(nm in indicators) result[[paste0(role,"_",nm)]]<-if(nrow(row)) row[[nm]][1L] else NA
  }
  result$S_correct_but_lower_ELBO_unselected<-isTRUE(result$S_full_structural_recovery) &&
    !isTRUE(sel$selected_S) && isTRUE(result$S_final_ordinary_elbo<result$U_final_ordinary_elbo)
  # selected_S is a logical column when read by R.
  result$S_selected_but_full_structure_worse<-isTRUE(sel$selected_S) &&
    isTRUE(result$U_full_structural_recovery) && !isTRUE(result$S_full_structural_recovery)
  result$S_selected_but_structure_worse<-isTRUE(sel$selected_S) && any(vapply(
    c("core_truth_recovered","exact_thresholded_factor_structure","full_true_feature_component_coverage"),
    function(nm) isTRUE(result[[paste0("U_",nm)]]) &&
      identical(result[[paste0("S_",nm)]],FALSE),logical(1)))
  for(nm in c("final_ordinary_elbo","R_A","R_B1","R_B2","dense_full_signal_nrmse_selected",
    "complete_contribution_nrmse_mean_missing_as_zero","feature_total_ise_mean"))
    result[[paste0("S_minus_U_",nm)]]<-result[[paste0("S_",nm)]]-result[[paste0("U_",nm)]]
  result
}))
new_csv(paired,file.path(root,"PAIRED_DATASET_SCIENCE.csv"))
fixed_winners<-all[all$fit_id%in%c(selection$U_winner_fit_id,selection$S_fit_id,
                                selection$enhanced_winner_fit_id),,drop=FALSE]
new_csv(fixed_winners,file.path(root,"WINNER_AND_S_SCIENCE.csv"))
role_delta<-do.call(rbind,lapply(seq_len(nrow(selection)),function(i) {
  sel<-selection[i,,drop=FALSE]
  u<-roles[roles$fit_id==sel$U_winner_fit_id,,drop=FALSE]
  s<-roles[roles$fit_id==sel$S_fit_id,,drop=FALSE]
  if(!nrow(u)||!nrow(s)) return(NULL)
  m<-merge(u,s,by=c("data_id","study","true_role","true_factor","true_id"),suffixes=c("_U","_S"))
  for(nm in c("loading_norm_ratio","loading_abs_cosine","factor_process_nrmse_missing_as_zero",
              "complete_contribution_nrmse_missing_as_zero"))
    m[[paste0("S_minus_U_",nm)]]<-m[[paste0(nm,"_S")]]-m[[paste0(nm,"_U")]]
  m$correct_U_damaged_in_S<-m$loading_structure_status_U=="correct" &
    m$loading_structure_status_S!="correct"
  m
}))
new_csv(if(is.null(role_delta)) data.frame(data_id=character()) else role_delta,
        file.path(root,"PAIRED_ROLE_DAMAGE.csv"))
nflag<-function(nm)sum(paired[[nm]]%in%TRUE)
lines<-c("# Stage 7L 前瞻 5×10U + S 科学报告", "",
 "五个预登记新数据，仅 C1。科学评价在全部既定端点及选择冻结后统一解封。",
 "模型、先验、CAVI、ELBO、生成器和 S 数学规则保持绑定版本。", "",
 sprintf("原池曾有完整正确端点：%d/5；U winner：%d/5；S：%d/5；增强 winner：%d/5。",
  nflag("U_pool_ever_correct"),nflag("U_full_structural_recovery"),
  nflag("S_full_structural_recovery"),nflag("enhanced_full_structural_recovery")),
 sprintf("选择S：%d/5；S正确但ELBO较低未选：%d；选S却损伤完整结构：%d。",
  nflag("selected_S"),nflag("S_correct_but_lower_ELBO_unselected"),
  nflag("S_selected_but_structure_worse")),"",
 "| 数据 | U池正确条数 | U winner完整正确 | S完整正确 | 选S | 最终完整正确 | S−U ELBO |",
 "| --- | ---: | --- | --- | --- | --- | ---: |")
for(i in seq_len(nrow(paired))) {
 r<-paired[i,,drop=FALSE]
 lines<-c(lines,sprintf("| %s | %d | %s | %s | %s | %s | %.6f |",r$data_id,
  r$U_pool_full_recovery_count,r$U_full_structural_recovery,r$S_full_structural_recovery,
  r$selected_S,r$enhanced_full_structural_recovery,r$S_minus_U_final_ordinary_elbo))
}
lines<-c(lines,"",sprintf("技术失败端点：%d/55；完整端点评价失败：%d。",
 sum(endpoints$status!="COMPLETE"),sum(endpoints$status=="COMPLETE" & !status$evaluated)),
 sprintf("U累计fit墙钟（单核）：%.2f小时；额外S累计：%.2f小时。",
  sum(all$elapsed_seconds[all$branch=="U"],na.rm=TRUE)/3600,
  sum(all$elapsed_seconds[all$branch=="S"],na.rm=TRUE)/3600),"",
 "## 解释边界", "",
 "完整正确同时要求core、exact及完整真特征覆盖；未完成499端点与评价失败不计为科学成功，空指标不代表结构错误。",
 "增强池包含原池，最高ELBO不下降不是改进证据。表中结构结果和shared/B1损伤决定是否有科学收益。",
 "PAIRED_ROLE_DAMAGE.csv逐真角色报告匹配状态和完整贡献误差；缺失块按0计入，不能只比较已匹配块。",
 "ALL_ENDPOINT_SCIENCE.csv保留全部55个计划端点、R/P/L、missing/misplaced/extra/duplicate、函数及尺度指标。",
 "独立重复数为5，不能将55候选或历史诊断当作55个独立重复。固定winner的sweep50重启效果不等于从0初始化有效、全局最优、一般识别或跨维度选择成功。",
 "S失败时不更换起点、时点或幅值；本批结束不自动补实验。","")
if(nflag("U_full_structural_recovery")==5L) lines<-c(lines,"本批五个U winner均正确，缺少winner失败救援信息。")
if(any(status$evaluated==FALSE & endpoints$status=="COMPLETE"))
 lines<-c(lines,"存在科学评价错误：本报告为带评价缺失的结果，须查EVALUATION_STATUS.csv，不能以缺失评价认定方法失败。")
new_lines(lines,file.path(output,"REPORT_STAGE7L_PROSPECTIVE_20261001_v1.md"))
new_lines(c("STAGE7L_EVALUATION_COMPLETE","fit_or_selection_changed=FALSE"),
          file.path(root,"EVALUATION_COMPLETE.txt"))
