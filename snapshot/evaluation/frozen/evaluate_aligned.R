# Local evaluation only: both 72-endpoint stages must be frozen and verified.
# Sourcing this file defines pure helpers and does not read fit/truth values.
# Direct execution requires FSYNC_EXEC_ROOT and FSYNC_OUTPUT_ROOT. No fitting.
ae_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
ae_close <- function(x,y,tolerance=1e-12) isTRUE(all.equal(x,y,tolerance=tolerance,check.attributes=FALSE))
ae_sha <- function(path) digest::digest(file=path,algo='sha256',serialize=FALSE)
ae_canonical <- function(path) normalizePath(path,winslash='/',mustWork=TRUE)
ae_inside <- function(path,root) identical(path,root)||startsWith(path,paste0(root,'/'))
ae_relative_path <- function(path) {
  ae_assert(length(path)==1L&&!is.na(path)&&nzchar(path)&&!grepl('(^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))',path),'Unsafe relative manifest path')
  path
}
ae_manifest_lines <- function(lines) {
  ae_assert(length(lines)>0L&&all(grepl('^[0-9a-fA-F]{64} [ *].+$',lines)),'Malformed SHA256 manifest')
  paths<-substring(lines,67L);invisible(lapply(paths,ae_relative_path))
  ae_assert(!anyDuplicated(paths),'Duplicate manifest paths')
  data.frame(relative_path=paths,sha256=tolower(substr(lines,1L,64L)),stringsAsFactors=FALSE)
}
ae_verify_record <- function(path,expected,bytes=NULL) {
  ae_assert(file.exists(path)&&!isTRUE(file.info(path)$isdir),paste('Missing file:',path))
  ae_assert(identical(ae_sha(path),tolower(as.character(expected))),paste('SHA256 mismatch:',path))
  if(!is.null(bytes))ae_assert(file.info(path)$size==as.numeric(bytes),paste('Byte-count mismatch:',path))
  invisible(TRUE)
}
ae_verify_manifest <- function(path,base) {
  ae_assert(file.exists(path),paste('Missing SHA256 manifest:',path))
  base<-ae_canonical(base);tab<-ae_manifest_lines(readLines(path,warn=FALSE))
  tab$path<-vapply(tab$relative_path,function(x)ae_canonical(file.path(base,x)),character(1))
  ae_assert(!anyDuplicated(tab$path)&&all(vapply(tab$path,ae_inside,logical(1),root=base)),'Manifest canonical path escapes/aliases authorized base')
  for(i in seq_len(nrow(tab)))ae_verify_record(tab$path[i],tab$sha256[i])
  tab
}
ae_monotonic_pass <- function(elbo) {
  is.numeric(elbo)&&length(elbo)>=2L&&all(is.finite(elbo))&&all(diff(elbo)>=-1e-7*pmax(1,abs(head(elbo,-1L))))
}
ae_valid_endpoints <- function(x) {
  ok<-x$status=='complete' & x$objective_valid & x$finite & x$metadata_pass & x$monotonic_pass & is.finite(x$ELBO)
  ok[is.na(ok)]<-FALSE;ok
}
ae_verify_selection <- function(plan,endpoints,winners,expected_n=144L,starts=12L) {
  keys<-c('stage','dataset_id','arm','start','seed','fit_id','L_f','L_s','M','reference_task_id')
  needed<-c(keys,'status','objective_valid','finite','metadata_pass','monotonic_pass','ELBO','fit_path','fit_sha256')
  ae_assert(all(keys%in%names(plan))&&all(needed%in%names(endpoints))&&all(needed%in%names(winners)),'Incomplete selection columns')
  ae_assert(nrow(plan)==expected_n&&nrow(endpoints)==expected_n&&!anyDuplicated(plan$fit_id)&&!anyDuplicated(endpoints$fit_id)&&!anyDuplicated(winners$fit_id),'Invalid endpoint cardinality')
  ae_assert(setequal(plan$fit_id,endpoints$fit_id)&&all(winners$fit_id%in%endpoints$fit_id),'Endpoint identity mismatch')
  for(key in keys)ae_assert(identical(as.character(plan[[key]]),as.character(endpoints[[key]][match(plan$fit_id,endpoints$fit_id)])),paste('Plan column mismatch:',key))
  ae_assert(all(endpoints$status%in%c('complete','failed','timeout','terminated','missing')),'Every planned endpoint must have a frozen terminal status')
  valid<-ae_valid_endpoints(endpoints)
  ae_assert(all(valid[endpoints$status=='complete']),'Complete endpoint lacks required validity gates')
  group<-function(x)paste(x$dataset_id,x$arm,sep='::')
  ae_assert(!anyDuplicated(group(winners))&&setequal(group(winners),unique(group(endpoints[valid,,drop=FALSE]))),'Winner group mismatch')
  for(key in needed)ae_assert(identical(as.character(winners[[key]]),as.character(endpoints[[key]][match(winners$fit_id,endpoints$fit_id)])),paste('Winner row differs from endpoint:',key))
  for(i in seq_len(nrow(winners))) {
    pool<-endpoints[group(endpoints)==group(winners[i,,drop=FALSE]),,drop=FALSE]
    ae_assert(nrow(pool)==starts&&setequal(pool$start,seq_len(starts)),'Incomplete within-arm start pool')
    pool<-pool[ae_valid_endpoints(pool),,drop=FALSE]
    winner<-pool$fit_id[order(-pool$ELBO,pool$fit_id)[1L]]
    ae_assert(identical(as.character(winner),as.character(winners$fit_id[i])),'Winner violates frozen maximum-ELBO/fit-id tie rule')
  }
  invisible(TRUE)
}
ae_verify_freezes <- function(exec_root,output_root,code_manifest) {
  # All prerequisites are tested before loading any .rds file.
  ae_assert(all(file.exists(file.path(output_root,paste0('stage',1:2,'_frozen'),'SELECTION_FROZEN.json'))),'Both stage1 and stage2 selection freezes are required; no truth has been read')
  code<-ae_verify_manifest(code_manifest,exec_root)
  plan<-read.csv(file.path(output_root,'PLAN.csv'),stringsAsFactors=FALSE)
  plan_local<-read.csv(file.path(exec_root,'config','PLAN.csv'),stringsAsFactors=FALSE)
  ae_assert(identical(plan,plan_local),'Returned plan differs from preregistered local plan')
  stages<-lapply(1:2,function(stage) {
    d<-file.path(output_root,paste0('stage',stage,'_frozen'))
    marker<-jsonlite::fromJSON(file.path(d,'SELECTION_FROZEN.json'))
    ae_assert(identical(as.integer(marker$stage),as.integer(stage))&&isFALSE(marker$truth_read)&&isFALSE(marker$test_response_read)&&marker$planned==72L&&marker$complete>=0L&&marker$complete<=72L&&marker$winners>=0L&&marker$winners<=6L,'Invalid/incomplete stage freeze marker')
    ae_verify_record(file.path(d,'FROZEN_SHA256SUMS.txt'),marker$freeze_manifest_sha256)
    ae_verify_record(code_manifest,marker$source_manifest_sha256)
    tab<-ae_verify_manifest(file.path(d,'FROZEN_SHA256SUMS.txt'),output_root)
    e<-read.csv(file.path(d,'ENDPOINTS.csv'),stringsAsFactors=FALSE)
    w<-read.csv(file.path(d,'WINNERS.csv'),stringsAsFactors=FALSE)
    ae_assert(marker$complete==sum(e$status=='complete')&&marker$winners==nrow(w),'Freeze marker/table count mismatch')
    available<-!is.na(e$fit_sha256)&nzchar(e$fit_sha256)
    need<-c('PLAN.csv',file.path(paste0('stage',stage,'_frozen'),c('ENDPOINTS.csv','WINNERS.csv')),e$fit_path[available],file.path('fits',e$fit_id[e$status!='missing'],'STATUS.csv'))
    ae_assert(all(gsub('\\\\','/',need)%in%gsub('\\\\','/',tab$relative_path)),'Stage SHA manifest lacks plan/table/fit/status coverage')
    for(i in which(available)) {
      ae_assert(identical(gsub('\\\\','/',e$fit_path[i]),paste0('fits/',e$fit_id[i],'/fit.rds')),'Unexpected fit path')
      ix<-match(ae_canonical(file.path(output_root,e$fit_path[i])),tab$path)
      ae_assert(!is.na(ix)&&identical(tolower(e$fit_sha256[i]),tab$sha256[ix]),'Fit hash column differs from frozen manifest')
    }
    ae_assert(all(e$stage==stage)&&all(w$stage==stage),'Stage table contains another stage')
    list(endpoints=e,winners=w,hashes=tab,marker=marker)
  })
  endpoints<-do.call(rbind,lapply(stages,`[[`,'endpoints'));winners<-do.call(rbind,lapply(stages,`[[`,'winners'))
  ae_verify_selection(plan,endpoints,winners)
  ae_assert(setequal(plan$dataset_id,paste0('NC01_BASE__rep0',1:3))&&setequal(plan$arm,c('original','mgp','original_hansen','mgp_hansen'))&&all(plan$L_f==2L)&&all(plan$L_s=='2;2')&&all(plan$M==3L),'Unexpected aligned experiment dimensions')
  for(id in unique(plan$dataset_id))for(k in 1:12)ae_assert(length(unique(plan$seed[plan$dataset_id==id&plan$start==k]))==1L,'Unpaired start seeds')
  list(plan=plan,endpoints=endpoints,winners=winners,code=code,stages=stages)
}
ae_verify_evaluation_lock <- function(exec_root,lock_path=file.path(exec_root,'audit','EVALUATION_FROZEN.json')) {
  ae_assert(file.exists(lock_path),'Missing pre-truth local EVALUATION_FROZEN.json; evaluator definition is not locked')
  lock<-jsonlite::fromJSON(lock_path)
  ae_assert(lock$schema_version==1L&&isFALSE(lock$truth_values_read)&&is.data.frame(lock$files),'Invalid evaluation definition freeze')
  needed<-c('scripts/evaluate_aligned.R','tests/test_evaluation_contract.R','audit/EVALUATION_DESIGN.md',
    'config/target_baseline.json','config/INPUT_BINDINGS.json','audit/LOCAL_EVALUATION_BINDINGS.json')
  files<-lock$files
  ae_assert(all(c('path','sha256','bytes')%in%names(files))&&!anyDuplicated(files$path)&&all(needed%in%files$path),'Evaluation lock lacks required source/input binding coverage')
  for(i in seq_len(nrow(files))) {
    rel<-ae_relative_path(files$path[i]);path<-ae_canonical(file.path(exec_root,rel))
    ae_assert(ae_inside(path,ae_canonical(exec_root)),'Evaluation lock path leaves research root')
    ae_verify_record(path,files$sha256[i],files$bytes[i])
  }
  invisible(list(path=lock_path,sha256=ae_sha(lock_path),files=files))
}
ae_legacy_environment <- function(formal_root) {
  dep<-file.path(formal_root,'scripts','dependencies');legacy<-file.path(dep,'legacy')
  paths<-c(file.path(dep,'stage7m_reconstruct.R'),file.path(legacy,'v0g_evaluation.R'),file.path(legacy,'v0h_evaluation.R'),file.path(formal_root,'scripts','nc_evaluation_setup.R'))
  expected<-c('51d83c1864dbc346bb4fee0e25feee7e099698080e2cfc38199970f7fd92ea27','6b765191657f88db46a027a95acd3917c602b35cb6cee8947d917652f9c6b271','45639c2d9f14e4df0187d519a62837ef3b3be090ed315123cd97c725f8576dd4',NA_character_)
  # nc_evaluation_setup is bound in target_baseline.json and independently checked by caller.
  for(i in 1:3)ae_verify_record(paths[i],expected[i])
  en<-new.env(parent=globalenv());en$assert<-ae_assert;en$`%||%`<-function(a,b)if(is.null(a))b else a
  for(path in paths)sys.source(path,envir=en,keep.source=FALSE)
  en$v0h_one_to_one_match<-en$nc_one_to_one_match_pruned
  attr(en,'source_bindings')<-data.frame(path=paths,sha256=vapply(paths,ae_sha,character(1)),stringsAsFactors=FALSE)
  en
}
ae_cosine <- function(x,y) {
  den<-sqrt(sum(x^2)*sum(y^2));if(!is.finite(den)||den<=0)return(NA_real_)
  min(1,abs(sum(x*y))/den)
}
ae_signed_correlation <- function(x,y) {
  if(length(x)<2L||length(x)!=length(y)||!all(is.finite(c(x,y)))||sd(x)==0||sd(y)==0)return(NA_real_)
  cor(x,y)
}
ae_correlation <- function(x,y) abs(ae_signed_correlation(x,y))
ae_factor_ppi <- function(gamma) {
  ae_assert(is.numeric(gamma)&&length(gamma)>0L&&all(is.finite(gamma))&&all(gamma>=0&gamma<=1),'Invalid loading inclusion probabilities')
  -expm1(sum(log1p(-gamma)))
}
ae_rank_to_raw <- function(order,ranks,L) {
  ranks<-as.integer(unlist(ranks));order<-as.integer(order)
  ae_assert(length(order)==L&&!anyDuplicated(order)&&setequal(order,seq_len(L))&&length(ranks)==L&&!anyNA(ranks)&&all(ranks>=0),'Invalid reported/raw factor permutation or rank')
  raw<-integer(L);raw[order]<-ranks;raw
}
ae_tau_moments <- function(fit) {
  st<-fit$factor_shrinkage_state
  if(is.null(st))return(list(shared=list(mean=rep(1,fit$L_f),elog=rep(0,fit$L_f)),specific=lapply(fit$L_s_by_study,function(L)list(mean=rep(1,L),elog=rep(0,L)))))
  ae_assert(st$version==1L&&st$control$mode=='mgp'&&st$p==fit$p&&st$L_f==fit$L_f&&identical(as.integer(st$L_s),as.integer(fit$L_s_by_study)),'Unexpected shrinkage state')
  calc<-function(z,L){ae_assert(length(z$shape)==L&&length(z$rate)==L&&all(is.finite(c(z$shape,z$rate)))&&all(z$shape>0)&all(z$rate>0),'Invalid Gamma shape/rate');list(mean=exp(cumsum(log(z$shape/z$rate))),elog=cumsum(digamma(z$shape)-log(z$rate)))}
  list(shared=calc(st$shared,fit$L_f),specific=lapply(seq_len(fit$S),function(s)calc(st$specific[[s]],fit$L_s_by_study[s])))
}
ae_basis_check <- function(fit,obs,fixed,legacy,meta=NULL,formal=TRUE) {
  fields<-c('S','n_s','p','d','K','K_total','C_g','time_g')
  ae_assert(all(vapply(fields,function(nm)ae_close(fit[[nm]],fixed$fit_metadata[[nm]]),logical(1))),'Fit/frozen metadata invariants differ')
  ae_assert(is.null(obs$Z)&&obs$dimensions$d==0L&&isFALSE(obs$truth_available_to_fit)&&isFALSE(obs$test_response_available_to_fit),'Observation truth/covariate contract failed')
  err<-max(abs(legacy$stage7m_C(fit$time_g,fixed$basis)-fit$C_g));ae_assert(err<=1e-12,'Frozen Linux basis differs from fitted C_g')
  ce<-max(vapply(seq_len(fit$S),function(s)max(vapply(seq_len(fit$n_s[s]),function(i)max(abs(legacy$stage7m_C(obs$time_obs[[s]][[i]],fixed$basis)-fixed$C[[s]][[i]])),numeric(1))),numeric(1)))
  ae_assert(ce<=1e-12,'Frozen C/time/design mismatch')
  ae_assert(fit$L_f==2L&&all(fit$L_s_by_study==2L)&&all(fit$M_f==3L)&&all(unlist(fit$M_s)==3L),'Unexpected fitted capacity')
  ae_rank_to_raw(fit$factor_order_shared,fit$list_effective_M,fit$L_f)
  for(s in seq_len(fit$S))ae_rank_to_raw(fit$factor_order_specific[[s]],fit$list_effective_M_spec[[s]],fit$L_s_by_study[s])
  pp<-apply(fit$mu_q_gamma_a,2,ae_factor_ppi)
  ae_assert(ae_close(pp[fit$factor_order_shared],fit$factor_ppi_shared),'Reported shared PPI permutation differs from raw PPI')
  for(s in seq_len(fit$S)) {
    pp<-apply(fit$mu_q_gamma_b[[s]],2,ae_factor_ppi)
    ae_assert(ae_close(pp[fit$factor_order_specific[[s]]],fit$factor_ppi_specific[[s]]),'Reported specific PPI permutation differs from raw PPI')
  }
  ae_assert(all(is.finite(unlist(fit[c('mu_q_nu_mu','mu_q_nu_phi','mu_q_nu_psi','mu_q_zeta','mu_q_xi','mu_q_a','mu_q_b_specific','mu_q_gamma_a','mu_q_gamma_b','ELBO')]))),'Nonfinite endpoint state')
  if(formal) {
    ae_assert(ae_monotonic_pass(fit$ELBO),'Ordinary ELBO violates the frozen cloud 1e-7 relative monotonicity gate')
    ae_assert(fit$S==2L&&fit$p==500L&&all(fit$n_s==30L)&&fit$K==9L&&fit$K_total==11L&&fit$d==0L&&fit$annealing_sweeps==99L&&fit$t1_sweeps==400L&&length(fit$ELBO)==400L&&isTRUE(fit$elbo_objective_valid),'Formal fit control mismatch')
    if(!is.null(meta)) {
      ae_assert(ae_close(tail(fit$ELBO,1L),meta$ELBO),'Endpoint ELBO differs from frozen table')
      ae_assert(identical(as.integer(fit$fit_seed),as.integer(meta$seed)),'Fit seed differs from frozen paired plan')
      ae_assert(grepl('^mgp',meta$arm)==!is.null(fit$factor_shrinkage_state),'Arm/state mismatch')
      fc<-fit$fit_control;pc<-fc$practical_control;pi<-fit$pre_score_interface
      ae_assert(ae_close(fc$anneal,c(1,1.9,100))&&fc$maxit==499L&&fc$tol_abs==1e-3&&fc$tol_rel==1e-5&&fc$lambda_orth==0&&fc$convergence_rule_requested=='practical'&&isFALSE(fc$continuation_used),'Frozen ordinary/anneal control mismatch')
      ae_assert(pc$min_t1==396L&&pc$consecutive==5L&&pc$max_t1==400L&&ae_close(pc$checkpoints,c(380L,400L)),'Frozen G12 control mismatch')
      ae_assert(pi$pre_score_sweeps==1L&&pi$function_initialization=='gram_unit_energy'&&pi$random_scale_calibration=='function'&&ae_close(pi$trace_sweeps,c(1L,380L,400L)),'Frozen pre-score/Gram calibration mismatch')
      ae_assert(fc$initialization==if(grepl('_hansen$',meta$arm))'residual_fpca'else'random','Arm initialization differs from plan')
      ae_assert(fit$d_0==500&&fit$n_cpus_used==1L&&all(fit$response_center==0)&&all(fit$response_scale==1),'Frozen no-scaling/prior/thread contract mismatch')
      if(grepl('^mgp',meta$arm))ae_assert(ae_close(fit$factor_shrinkage_state$control,list(mode='mgp',a1=2.1,a2=3.1,rate=1)),'MGP hyperparameters differ from registered arm')
    }
  }
  ae_tau_moments(fit)
  c(C_g_max_absolute_error=err,C_max_absolute_error=ce)
}
ae_weighted_loss <- function(pred,target,studies=seq_along(pred)) {
  mse<-energy<-mae<-0
  for(s in studies) {
    ae_assert(length(pred[[s]])==length(target[[s]])&&length(pred[[s]])>0L,'Loss subject-count mismatch')
    for(i in seq_along(pred[[s]])) {
      a<-as.matrix(pred[[s]][[i]]);b<-as.matrix(target[[s]][[i]])
      ae_assert(identical(dim(a),dim(b))&&all(is.finite(c(a,b))),'Loss dimensions/nonfinite mismatch')
      weight<-1/length(studies)/length(pred[[s]])
      mse<-mse+weight*mean((a-b)^2);energy<-energy+weight*mean(b^2);mae<-mae+weight*mean(abs(a-b))
    }
  }
  c(mse=mse,truth_energy=energy,rmse=sqrt(mse),mae=mae,nrmse=if(energy>0)sqrt(mse/energy)else NA_real_)
}
ae_formal_structure <- function(fit,truth,obs,legacy) {
  data<-list(true_params=truth,time_obs=obs$time_obs)
  tg<-legacy$v0h_truth_global_catalog(data);cg<-legacy$v0h_candidate_global_catalog(fit,'ppi_selected',.5)
  gl<-legacy$v0h_classify_catalog_match(tg,cg,'ppi_selected','global_direction',0L,.8,.8)
  st<-lapply(seq_len(fit$S),function(s) {
    ca<-legacy$v0h_study_catalogs(tg,cg,s)
    legacy$v0h_classify_catalog_match(ca$truth,ca$candidate,'ppi_selected','study_role',s,.8,.8)
  })
  tt<-do.call(rbind,lapply(st,`[[`,'truth'));cc<-do.call(rbind,lapply(st,`[[`,'candidate'))
  core<-all(gl$truth$status=='correct')&&all(tt$status=='correct')
  counts<-c(sum(fit$factor_ppi_shared>=.5),vapply(fit$factor_ppi_specific,function(x)sum(x>=.5),integer(1)))
  exactcounts<-identical(as.integer(counts),c(1L,1L,1L))
  exact<-core&&exactcounts&&!any(gl$candidate$status%in%c('extra','duplicate'))&&!any(cc$status%in%c('extra','duplicate'))
  sumrow<-data.frame(formal_nc_loading_core_08=core,formal_nc_exact_loading_structure_08=exact,formal_nc_ppi_ge_half_counts_exact=exactcounts,
    formal_global_correct=sum(gl$truth$status=='correct'),formal_global_missing=sum(gl$truth$status=='missing'),formal_global_misplaced=sum(gl$truth$status=='misplaced'),formal_global_duplicate=sum(gl$candidate$status=='duplicate'),formal_global_extra=sum(gl$candidate$status=='extra'),
    formal_study_correct=sum(tt$status=='correct'),formal_study_missing=sum(tt$status=='missing'),formal_study_misplaced=sum(tt$status=='misplaced'),formal_study_duplicate=sum(cc$status=='duplicate'),formal_study_extra=sum(cc$status=='extra'))
  # These are exactly the formal catalog/one-to-one loading definitions. They
  # do not constitute a function/trajectory/score or full-basin recovery claim.
  paths<-legacy$v0g_reported_paths_by_study(fit,data)
  shared<-legacy$nc_v0g_shared_candidate_metrics(fit,paths,.10)
  shared$raw_factor<-fit$factor_order_shared[shared$candidate_factor]
  block<-legacy$v0h_block_subspace_metrics(fit,data,'ppi_selected',legacy$v0h_default_thresholds())
  list(summary=sumrow,global_truth=gl$truth,global_candidate=gl$candidate,study_truth=tt,study_candidate=cc,shared_support=shared,block_subspace=block)
}
ae_raw_metrics <- function(fit,truth,obs,fixed) {
  ranka<-ae_rank_to_raw(fit$factor_order_shared,fit$list_effective_M,fit$L_f)
  rankb<-lapply(seq_len(fit$S),function(s)ae_rank_to_raw(fit$factor_order_specific[[s]],fit$list_effective_M_spec[[s]],fit$L_s_by_study[s]))
  tau<-ae_tau_moments(fit);factors<-roles<-energies<-list()
  for(blockid in 0:fit$S) {
    sh<-blockid==0L;load<-if(sh)fit$mu_q_a_original else fit$mu_q_b_specific_original[[blockid]]
    gam<-if(sh)fit$mu_q_gamma_a else fit$mu_q_gamma_b[[blockid]]
    target<-if(sh)truth$a_true[,1L]else truth$b_true[[blockid]][,1L]
    ranks<-if(sh)ranka else rankb[[blockid]];tm<-if(sh)tau$shared else tau$specific[[blockid]]
    relevant<-if(sh)seq_len(fit$S)else blockid
    one<-list()
    for(l in seq_len(ncol(load))) {
      bystudy<-lapply(relevant,function(s) {
        coef<-if(sh)fit$mu_q_nu_phi[[l]]else fit$mu_q_nu_psi[[s]][[l]]
        scores<-if(sh)fit$mu_q_zeta[[s]][[l]]else fit$mu_q_xi[[s]][[l]]
        path<-lapply(seq_len(fit$n_s[s]),function(i)as.vector(fixed$C[[s]][[i]]%*%coef%*%scores[i,]))
        targetpath<-if(sh)truth$f_true_values[[s]]else truth$g_true_values[[s]]
        te<-unlist(lapply(targetpath,function(x)as.matrix(x)[,1L]),use.names=FALSE)
        pe<-mean(vapply(path,function(x)mean(x^2),numeric(1)))
        data.frame(block=if(sh)'shared'else'specific',block_study=blockid,raw_factor=l,evaluation_study=s,
          process_mean_energy_equal_subject=pe,contribution_mean_energy_equal_subject=pe*mean(load[,l]^2),
          own_role_trajectory_abs_cor_pooled_rows=ae_correlation(unlist(path,use.names=FALSE),te),
          own_role_trajectory_loading_aligned_cor_pooled_rows=sign(sum(load[,l]*target))*ae_signed_correlation(unlist(path,use.names=FALSE),te),
          function_mean_energy_equal_subject=mean(vapply(seq_len(fit$n_s[s]),function(i)mean(rowSums((fixed$C[[s]][[i]]%*%coef)^2)),numeric(1))),
          score_mean_squarednorm=mean(rowSums(scores^2)))
      })
      es<-do.call(rbind,bystudy);energies[[length(energies)+1L]]<-es
      ppi<-ae_factor_ppi(gam[,l]);co<-ae_cosine(load[,l],target)
      cosall<-vapply(c(list(truth$a_true[,1L]),lapply(truth$b_true,function(x)x[,1L])),function(t)ae_cosine(load[,l],t),numeric(1))
      validcos<-replace(cosall,!is.finite(cosall),-Inf);besttruth<-if(all(!is.finite(cosall)))NA_integer_ else which.max(validcos)-1L
      one[[l]]<-data.frame(block=if(sh)'shared'else'specific',study=blockid,raw_factor=l,reported_factor=match(l,if(sh)fit$factor_order_shared else fit$factor_order_specific[[blockid]]),
        factor_ppi=ppi,live_ppi_gt_half=ppi>.5,expected_support=sum(gam[,l]),loading_ppi_gt_half_count=sum(gam[,l]>.5),loading_mean_norm_original=sqrt(sum(load[,l]^2)),
        expected_tau=tm$mean[l],expected_log_tau=tm$elog[l],fpca_effective_M=ranks[l],own_role_loading_abs_cosine=co,
        own_role_loading_core_09=ppi>.5&&is.finite(co)&&co>=.9,own_role_loading_core_08=ppi>.5&&is.finite(co)&&co>=.8,
        min_own_role_trajectory_abs_cor=if(all(is.finite(es$own_role_trajectory_abs_cor_pooled_rows)))min(es$own_role_trajectory_abs_cor_pooled_rows)else NA_real_,
        min_own_role_loading_aligned_trajectory_cor=if(all(is.finite(es$own_role_trajectory_loading_aligned_cor_pooled_rows)))min(es$own_role_trajectory_loading_aligned_cor_pooled_rows)else NA_real_,
        mean_contribution_energy=mean(es$contribution_mean_energy_equal_subject),min_study_contribution_energy=min(es$contribution_mean_energy_equal_subject),
        cosine_true_shared=cosall[1],cosine_true_specific_1=cosall[2],cosine_true_specific_2=cosall[3],best_true_role_study=besttruth,best_true_role_is_offblock=!is.na(besttruth)&&besttruth!=blockid)
    }
    tab<-do.call(rbind,one);factors[[length(factors)+1L]]<-tab
    candidates<-which(tab$own_role_loading_core_09)
    matched<-if(length(candidates))candidates[order(-tab$own_role_loading_abs_cosine[candidates],tab$raw_factor[candidates])[1L]]else NA_integer_
    chosen<-if(is.na(matched))NA_integer_ else tab$raw_factor[matched]
    support<-target!=0;selected<-if(is.na(chosen))rep(FALSE,length(support))else gam[,chosen]>.5
    tp<-sum(selected&support);fp<-sum(selected&!support);fn<-sum(!selected&support)
    roles[[length(roles)+1L]]<-data.frame(block=if(sh)'shared'else'specific',study=blockid,truth_factor_count=1L,fitted_capacity=ncol(load),live_factor_count=sum(tab$live_ppi_gt_half),
      loading_core_present_09=any(tab$own_role_loading_core_09),loading_core_present_08=any(tab$own_role_loading_core_08),duplicate_own_role_core_count_09=max(0L,sum(tab$own_role_loading_core_09)-1L),
      best_own_role_loading_abs_cosine=if(any(is.finite(tab$own_role_loading_abs_cosine)))max(tab$own_role_loading_abs_cosine,na.rm=TRUE)else NA_real_,
      matched_raw_factor_09=chosen,matched_fpca_M_09=if(is.na(matched))NA_integer_ else tab$fpca_effective_M[matched],matched_fpca_exact_09=!is.na(matched)&&tab$fpca_effective_M[matched]==2L,
      loading_core_and_trajectory08_present=any(tab$own_role_loading_core_09&is.finite(tab$min_own_role_loading_aligned_trajectory_cor)&tab$min_own_role_loading_aligned_trajectory_cor>=.8),
      support_precision_ppi_gt_half=if(tp+fp>0)tp/(tp+fp)else 0,support_recall_ppi_gt_half=tp/sum(support),support_true_positive=tp,support_false_positive=fp,support_false_negative=fn)
  }
  fr<-do.call(rbind,factors);rr<-do.call(rbind,roles)
  summary<-data.frame(loading_role_core_09=all(rr$loading_core_present_09),loading_role_core_08=all(rr$loading_core_present_08),
    ppi_counts_gt_half_exact=all(rr$live_factor_count==1L),legacy_PPI_exact_09=all(rr$loading_core_present_09)&&all(rr$live_factor_count==1L),
    loading_core_and_trajectory08_all_roles=all(rr$loading_core_and_trajectory08_present),matched_fpca_all_roles_exact_09=all(rr$matched_fpca_exact_09),
    live_shared=rr$live_factor_count[1L],live_specific_1=rr$live_factor_count[2L],live_specific_2=rr$live_factor_count[3L],
    fpca_ranks_raw_shared=paste(ranka,collapse=';'),fpca_ranks_raw_specific_1=paste(rankb[[1]],collapse=';'),fpca_ranks_raw_specific_2=paste(rankb[[2]],collapse=';'))
  list(summary=summary,factors=fr,roles=rr,energies=do.call(rbind,energies))
}
ae_winner_losses <- function(fit,truth,obs,fixed,test_truth,legacy) {
  estimate<-target<-setNames(vector('list',4),c('signal','shared','specific','mean'))
  for(nm in names(estimate))estimate[[nm]]<-target[[nm]]<-vector('list',fit$S)
  for(s in seq_len(fit$S))for(i in seq_len(fit$n_s[s])) {
    a<-legacy$stage7m_estimate(fit,fixed$basis,s,i,obs$time_obs[[s]][[i]],'raw')
    estimate$shared[[s]][[i]]<-a$shared_contribution;estimate$specific[[s]][[i]]<-a$specific_contribution
    estimate$mean[[s]][[i]]<-a$mean;estimate$signal[[s]][[i]]<-a$mean+a$factor
    target$shared[[s]][[i]]<-truth$f_true_values[[s]][[i]]%*%t(truth$a_true)
    target$specific[[s]][[i]]<-truth$g_true_values[[s]][[i]]%*%t(truth$b_true[[s]])
    target$mean[[s]][[i]]<-truth$mu_true_values[[s]][[i]];target$signal[[s]][[i]]<-truth$signal_true_values[[s]][[i]]
    ae_assert(ae_close(target$mean[[s]][[i]]+target$shared[[s]][[i]]+target$specific[[s]][[i]],target$signal[[s]][[i]],1e-10),'Stored true signal decomposition fails')
  }
  rows<-list()
  for(nm in names(estimate))for(scope in c('all_studies',paste0('study_',seq_len(fit$S)))) {
    studies<-if(scope=='all_studies')seq_len(fit$S)else as.integer(sub('study_','',scope))
    rows[[length(rows)+1L]]<-cbind(data.frame(scope=scope,representation='raw_original',object=nm),as.data.frame(as.list(ae_weighted_loss(estimate[[nm]],target[[nm]],studies))))
  }
  pred<-legacy$stage7m_prediction(fit,obs,fixed$basis,list(time=test_truth$time))
  for(rep in c('raw','reported'))for(object in c('test_response','test_signal')) {
    yy<-if(object=='test_response')test_truth$Y else test_truth$signal
    rows[[length(rows)+1L]]<-cbind(data.frame(scope='heldout_all_studies',representation=rep,object=object),as.data.frame(as.list(ae_weighted_loss(pred[[rep]],yy))))
  }
  do.call(rbind,rows)
}
ae_main <- function() {
  options(warn=1);Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
  exec_root<-Sys.getenv('FSYNC_EXEC_ROOT',unset='');output_root<-Sys.getenv('FSYNC_OUTPUT_ROOT',unset='')
  ae_assert(nzchar(exec_root)&&nzchar(output_root),'Set FSYNC_EXEC_ROOT and FSYNC_OUTPUT_ROOT explicitly')
  exec_root<-ae_canonical(exec_root);output_root<-ae_canonical(output_root)
  result_dir<-Sys.getenv('FSYNC_EVALUATION_DIR',unset=file.path(exec_root,'results','aligned_evaluation_v1'))
  ae_assert(!file.exists(result_dir),'Refusing to overwrite an evaluation directory')
  # Existing parent must remain within the authorized E research root.
  ae_assert(ae_inside(ae_canonical(dirname(result_dir)),exec_root),'Evaluation parent leaves research root')
  manifest<-Sys.getenv('FSYNC_CODE_MANIFEST',unset=file.path(exec_root,'delivery','cloud_package_v1','CODE_SHA256SUMS.txt'))
  frozen<-ae_verify_freezes(exec_root,output_root,manifest)
  evaluation_lock<-ae_verify_evaluation_lock(exec_root)
  contract<-jsonlite::fromJSON(file.path(exec_root,'config','target_baseline.json'),simplifyVector=FALSE)
  formal_root<-contract$formal_root
  for(x in contract$evaluation_source_bindings)ae_verify_record(x$path,x$sha256,x$bytes)
  legacy<-ae_legacy_environment(formal_root)
  input<-jsonlite::fromJSON(file.path(exec_root,'config','INPUT_BINDINGS.json'))
  truth_bind<-jsonlite::fromJSON(file.path(exec_root,'audit','LOCAL_EVALUATION_BINDINGS.json'))
  expected_files<-unlist(lapply(contract$datasets,`[[`,'files'),recursive=FALSE)
  allbindings<-list()
  for(x in expected_files) {
    ae_verify_record(x$path,x$sha256,x$bytes)
    allbindings[[length(allbindings)+1L]]<-data.frame(kind=x$role,path=x$path,sha256=x$sha256,bytes=x$bytes)
  }
  for(i in seq_len(nrow(truth_bind))) {
    x<-Filter(function(z)identical(z$path,gsub('\\\\','/',truth_bind$local_path[i])),expected_files)
    ae_assert(length(x)==1L&&identical(x[[1L]]$sha256,truth_bind$sha256[i]),'Local truth binding differs from alignment contract')
  }
  ids<-sort(unique(frozen$plan$dataset_id));dataset<-list()
  for(id in ids) {
    row<-input[input$data_id==id,,drop=FALSE];ae_assert(nrow(row)==3L,'Incomplete observed input bindings')
    for(k in seq_len(nrow(row)))ae_verify_record(row$origin[k],row$sha256[k],row$bytes[k])
    obs<-readRDS(row$origin[row$role=='observation']);fixed<-readRDS(row$origin[row$role=='metadata']);blind<-readRDS(row$origin[row$role=='test_design'])
    ae_assert(identical(obs$data_id,id),'Observation dataset ID mismatch')
    dataset[[id]]<-list(obs=obs,fixed=fixed,blind=blind)
  }
  # PASS ONE: every valid terminal endpoint's state/basis checks precede truth.
  # Failed/missing endpoints remain in the planned denominator with recovery FALSE.
  valid_indices<-which(ae_valid_endpoints(frozen$endpoints));basis_rows<-list()
  for(i in valid_indices) {
    meta<-frozen$endpoints[i,,drop=FALSE];path<-file.path(output_root,meta$fit_path)
    ae_verify_record(path,meta$fit_sha256);fit<-readRDS(path);ds<-dataset[[meta$dataset_id]]
    checks<-ae_basis_check(fit,ds$obs,ds$fixed,legacy,meta)
    basis_rows[[length(basis_rows)+1L]]<-cbind(meta[c('fit_id','dataset_id','arm','start')],as.data.frame(as.list(checks)))
    rm(fit);if(i%%12L==0L){cat('Pre-truth basis gates',i,'/144\n');invisible(gc())}
  }
  ae_assert(dir.create(result_dir,recursive=FALSE),'Cannot create evaluation directory')
  write.csv(if(length(basis_rows))do.call(rbind,basis_rows)else data.frame(fit_id=character()),file.path(result_dir,'PRE_TRUTH_BASIS_CHECKS.csv'),row.names=FALSE)
  write.csv(do.call(rbind,allbindings),file.path(result_dir,'LOCAL_INPUT_HASH_VERIFICATION.csv'),row.names=FALSE)
  evaluator_sources<-rbind(attr(legacy,'source_bindings'),data.frame(path=file.path(exec_root,'scripts','evaluate_aligned.R'),sha256=ae_sha(file.path(exec_root,'scripts','evaluate_aligned.R'))))
  write.csv(evaluator_sources,file.path(result_dir,'EVALUATION_SOURCE_BINDINGS.csv'),row.names=FALSE)
  jsonlite::write_json(list(planned_endpoints=144L,valid_endpoints=length(valid_indices),frozen_winners=nrow(frozen$winners),all_frozen_endpoints_and_winners_verified=TRUE,both_stage_freezes_verified=TRUE,evaluation_definition_freeze_sha256=evaluation_lock$sha256,all_valid_bases_checked_before_truth=TRUE,truth_used_for_selection=FALSE,started_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(result_dir,'EVALUATION_STARTED.json'),auto_unbox=TRUE,pretty=TRUE)
  tables<-list();append<-function(name,x,meta) {
    if(!nrow(x))return(invisible(NULL))
    header<-meta[rep(1L,nrow(x)),c('fit_id','dataset_id','stage','arm','start','seed'),drop=FALSE]
    header$winner<-meta$fit_id%in%frozen$winners$fit_id
    tables[[name]][[length(tables[[name]])+1L]]<<-cbind(header,x)
  }
  for(id in ids) {
    if(!any(frozen$endpoints$dataset_id[valid_indices]==id))next
    ds<-dataset[[id]];bind<-truth_bind[truth_bind$data_id==id,,drop=FALSE]
    ae_assert(nrow(bind)==2L,'Expected truth and test_truth bindings')
    for(k in seq_len(nrow(bind)))ae_verify_record(bind$local_path[k],bind$sha256[k],bind$bytes[k])
    # FIRST TRUTH VALUE READ: all global hashes, selections, and basis gates above passed.
    sealed<-readRDS(bind$local_path[bind$role=='truth']);test<-readRDS(bind$local_path[bind$role=='test_truth']);tr<-sealed$data$true_params
    ae_assert(identical(sealed$data_id,id)&&tr$S==2L&&tr$p==500L&&tr$L_f==1L&&all(tr$L_s_by_study==1L)&&tr$d_use==0L&&all(tr$n_s==30L)&&all(lengths(tr$lambda_phi_true)==2L),'Frozen NC truth identity/dimensions mismatch')
    ae_assert(ae_close(test$time,ds$blind$time),'Heldout truth/design times mismatch')
    for(i in intersect(which(frozen$endpoints$dataset_id==id),valid_indices)) {
      meta<-frozen$endpoints[i,,drop=FALSE];path<-file.path(output_root,meta$fit_path);ae_verify_record(path,meta$fit_sha256);fit<-readRDS(path)
      raw<-ae_raw_metrics(fit,tr,ds$obs,ds$fixed);formal<-ae_formal_structure(fit,tr,ds$obs,legacy)
      append('ENDPOINT_SCIENCE',cbind(data.frame(evaluation_available=TRUE,terminal_status=meta$status),raw$summary,formal$summary),meta)
      append('RAW_FACTOR_DIAGNOSTICS',raw$factors,meta);append('ROLE_RECOVERY',raw$roles,meta);append('FACTOR_CONTRIBUTION_BY_STUDY',raw$energies,meta)
      for(nm in c('global_truth','global_candidate','study_truth','study_candidate','shared_support','block_subspace'))append(paste0('FORMAL_',toupper(nm)),formal[[nm]],meta)
      if(meta$fit_id%in%frozen$winners$fit_id)append('WINNER_LOSSES',ae_winner_losses(fit,tr,ds$obs,ds$fixed,test,legacy),meta)
      rm(fit);cat('Evaluated',meta$fit_id,'\n')
    }
    rm(sealed,test,tr);invisible(gc())
  }
  # Freeze failures contribute zero recovery to the planned denominator, but no
  # fabricated continuous errors or rank/count values are supplied.
  for(i in setdiff(seq_len(nrow(frozen$endpoints)),valid_indices)) {
    meta<-frozen$endpoints[i,,drop=FALSE]
    if(length(tables$ENDPOINT_SCIENCE)) {
      template<-tables$ENDPOINT_SCIENCE[[1L]][1L,,drop=FALSE]
      scientific<-template[,setdiff(names(template),c('fit_id','dataset_id','stage','arm','start','seed','winner')),drop=FALSE]
      for(nm in names(scientific))scientific[[nm]]<-if(is.logical(scientific[[nm]]))FALSE else if(is.character(scientific[[nm]]))NA_character_ else NA_real_
      scientific$evaluation_available<-FALSE;scientific$terminal_status<-meta$status
    } else scientific<-data.frame(evaluation_available=FALSE,terminal_status=meta$status,
      loading_role_core_09=FALSE,loading_role_core_08=FALSE,ppi_counts_gt_half_exact=FALSE,legacy_PPI_exact_09=FALSE,
      loading_core_and_trajectory08_all_roles=FALSE,matched_fpca_all_roles_exact_09=FALSE,
      formal_nc_loading_core_08=FALSE,formal_nc_exact_loading_structure_08=FALSE,formal_nc_ppi_ge_half_counts_exact=FALSE)
    append('ENDPOINT_SCIENCE',scientific,meta)
  }
  for(nm in names(tables)){tab<-do.call(rbind,tables[[nm]]);rownames(tab)<-NULL;write.csv(tab,file.path(result_dir,paste0(nm,'.csv')),row.names=FALSE)}
  endpoint<-do.call(rbind,tables$ENDPOINT_SCIENCE)
  ae_assert(nrow(endpoint)==144L&&sum(endpoint$winner)==nrow(frozen$winners),'Evaluation row count mismatch')
  write.csv(endpoint[endpoint$winner,,drop=FALSE],file.path(result_dir,'WINNER_SCIENCE.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(result_dir,'RUNTIME.txt'))
  writeLines(c('144 planned frozen endpoints; at most12 frozen within-dataset/arm winners. Failed/missing terminal endpoints count as no recovery in the planned denominator; their scientific continuous metrics are unavailable, not zero. Statistical replication unit: three datasets, not144 starts.',
    'Primary loading_role_core_09: same-role raw loading cosine >=0.9 and factor-PPI >0.5 in each true role; this alone is NOT full correct-basin or trajectory recovery.',
    'formal_nc_loading_core_08: exact historical global/study one-to-one loading catalog criterion, cosine >=0.8 and PPI >=0.5; it also does not require correct trajectories.',
    'formal_nc_exact_loading_structure_08 adds exact counts and no extra/duplicate; primary legacy_PPI_exact_09 is kept separately.',
    'FORMAL_SHARED_SUPPORT uses original reported paths interpolated from fitted grid and the original 0.10 minimum-study/maximum-retained-pooled contribution RMS ratio.',
    'FACTOR_CONTRIBUTION_BY_STUDY uses exact frozen C, raw means and equal-subject average energy. Trajectory correlations use pooled observed rows.',
    'loading_core_and_trajectory08_all_roles additionally requires one same-role core candidate with observed signed trajectory correlation times sign(dot(estimated loading,true loading)) >=0.8 in every relevant study; this is an added diagnostic, not a validated universal basin test.',
    'Matched FPCA: select highest-cosine retained >=0.9 candidate per true role, raw-index tie break, then map reported effective rank back through the frozen permutation. Missing role => NA rank and exact FALSE.',
    'Winner NRMSE uses all raw fitted factors without truth-selected pruning; equal studies, then equal subjects, then observed times and variables. Heldout RMSE is new-time prediction for existing subjects, not new-subject risk.',
    'Support precision/recall uses PPI>0.5 on the matched raw loading, missing roles treated as all-zero support.',
    'No cross-prior ELBO comparison. No true-data-selected new starts. No fitted basis regeneration on Windows.',
    'Unassessed: full eigenfunction/component score recovery, uncertainty coverage, population covariance error and prospective independent-data generalization.'),file.path(result_dir,'METRIC_DEFINITIONS.txt'))
  files<-list.files(result_dir,full.names=TRUE);writeLines(vapply(files,function(p)paste(ae_sha(p),basename(p),sep='  '),character(1)),file.path(result_dir,'RESULT_SHA256SUMS.txt'))
  writeLines('PASS: both frozen stages evaluated without fitting, reselection or basis regeneration.',file.path(result_dir,'PASS.txt'))
  cat('Evaluation completed:',result_dir,'\n')
  invisible(result_dir)
}
if(sys.nframe()==0L)ae_main()
