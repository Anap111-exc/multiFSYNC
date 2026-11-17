# Evaluation definitions only. Existing legacy metrics are left unchanged.
stage7m_auc<-function(ppi,support) {
 np<-sum(support);nn<-sum(!support)
 auc<-if(np&&nn)(sum(rank(ppi,ties.method='average')[support])-np*(np+1)/2)/(np*nn)else NA_real_
 ord<-order(-ppi);groups<-split(seq_along(ppi),factor(ppi[ord],levels=unique(ppi[ord])))
 tp<-fp<-0;ap<-0
 for(ix in groups){delta<-sum(support[ord[ix]]);tp<-tp+delta;fp<-fp+length(ix)-delta;if(np)ap<-ap+(delta/np)*tp/(tp+fp)}
 selected<-ppi>=.5;tp<-sum(selected&support);fp<-sum(selected&!support);fn<-sum(!selected&support)
 precision<-if(tp+fp)tp/(tp+fp)else 0;recall<-if(np)tp/np else NA_real_
 data.frame(roc_auc=auc,pr_average_precision=if(np)ap else NA_real_,precision_ppi05=precision,recall_ppi05=recall,f1_ppi05=if(is.finite(recall)&&precision+recall>0)2*precision*recall/(precision+recall)else 0,support_size=np,selected_size=sum(selected),true_positive=tp,false_positive=fp,false_negative=fn)
}
stage7m_curve_error<-function(est,true,w) {
 est<-as.matrix(est);true<-as.matrix(true);assert(identical(dim(est),dim(true)),'Curve error dimensions')
 mse<-sum(w*(est-true)^2)/ncol(true);energy<-sum(w*true^2)/ncol(true)
 c(mse=mse,truth_energy=energy,rmse=sqrt(mse/sum(w)),nrmse=if(energy>0)sqrt(mse/energy)else NA_real_)
}
stage7m_cos<-function(x,y,w=rep(1,length(x))) {
 den<-sqrt(sum(w*x^2)*sum(w*y^2));if(den>1e-20)sum(w*x*y)/den else NA_real_
}
stage7m_score_error<-function(est,true) {
 reason<-if(any(!is.finite(est))||any(!is.finite(true)))'nonfinite'else if(length(est)<2)'fewer_than_two_scores'else if(sd(est)==0)'zero_estimate_variance'else if(sd(true)==0)'zero_truth_variance'else ''
 data.frame(correlation=if(nzchar(reason))NA_real_ else cor(est,true),correlation_reason=reason,score_rmse=sqrt(mean((est-true)^2)),score_nrmse=if(sum(true^2)>0)sqrt(sum((est-true)^2)/sum(true^2))else NA_real_,calibration_slope=if(var(true)>0)cov(true,est)/var(true)else NA_real_,estimated_sd=sd(est),truth_sd=sd(true))
}
stage7m_component_map<-function(true,estimate,w) {
 mt<-ncol(true);me<-ncol(estimate);assert(mt==2L&&me%in%c(1L,2L),'Component map expects frozen M=2, PVE output<=2')
 sim<-matrix(0,mt,me)
 for(i in 1:mt)for(j in 1:me){v<-stage7m_cos(true[,i],estimate[,j],w);sim[i,j]<-if(is.finite(v))abs(v)else 0}
 if(me==1L){assignment<-c(0L,0L);assignment[which.max(sim[,1])]<-1L}else assignment<-if(sim[1,1]+sim[2,2]>=sim[1,2]+sim[2,1])c(1L,2L)else c(2L,1L)
 fs<-rep(1,mt)
 for(i in 1:mt)if(assignment[i]>0){val<-sum(w*true[,i]*estimate[,assignment[i]]);fs[i]<-if(val<0)-1 else 1}
 list(assignment=assignment,function_sign=fs,similarity=sim)
}
stage7m_grid_estimate<-function(fit,basis,t,representation,reported) {
 C<-stage7m_C(t,basis)
 mu<-lapply(fit$mu_q_nu_mu_original,function(x)C%*%do.call(cbind,x))
 beta<-lapply(seq_len(fit$d),function(r)C%*%do.call(cbind,lapply(fit$mu_q_nu_beta_original,`[[`,r)))
 if(representation=='raw') {
  shared<-lapply(fit$mu_q_nu_phi,function(x)C%*%x);specific<-lapply(fit$mu_q_nu_psi,function(x)lapply(x,function(y)C%*%y))
  scores_shared<-fit$mu_q_zeta;scores_specific<-fit$mu_q_xi;A<-fit$mu_q_a_original;B<-fit$mu_q_b_specific_original
 } else {
  shared<-lapply(reported$shared,function(x)C%*%x);specific<-lapply(reported$specific,function(x)lapply(x,function(y)C%*%y))
  scores_shared<-lapply(seq_len(fit$S),function(s)lapply(fit$list_Zeta_hat,function(x)x[sum(head(fit$n_s,s-1))+seq_len(fit$n_s[s]),,drop=FALSE]));scores_specific<-fit$list_Zeta_hat_spec;A<-fit$mu_q_a_hat;B<-fit$mu_q_b_specific_hat
 }
 list(mean=mu,beta=beta,shared=shared,specific=specific,scores_shared=scores_shared,scores_specific=scores_specific,A=A,B=B)
}
stage7m_grid_subject<-function(estimate,obs,s,i) {
 sh<-lapply(seq_along(estimate$shared),function(l)as.vector(estimate$shared[[l]]%*%estimate$scores_shared[[s]][[l]][i,]))
 sp<-lapply(seq_along(estimate$specific[[s]]),function(l)as.vector(estimate$specific[[s]][[l]]%*%estimate$scores_specific[[s]][[l]][i,]))
 nr<-nrow(estimate$mean[[s]]);p<-ncol(estimate$mean[[s]]);gc<-sc<-bc<-matrix(0,nr,p)
 for(l in seq_along(sh))gc<-gc+outer(sh[[l]],estimate$A[,l])
 for(l in seq_along(sp))sc<-sc+outer(sp[[l]],estimate$B[[s]][,l])
 for(r in seq_along(estimate$beta))bc<-bc+obs$Z[[s]][i,r]*estimate$beta[[r]]
 list(mean=estimate$mean[[s]],beta=bc,shared=sh,specific=sp,shared_contribution=gc,specific_contribution=sc,factor=gc+sc,signal=estimate$mean[[s]]+bc+gc+sc)
}
stage7m_endpoint_science<-function(fit,obs,fixed,truth,historical,pred,test_truth) {
 grid<-fixed$grid;w<-fixed$weights;tr<-truth$data$true_params;reported<-stage7m_reported_coefficients(fit)
 ftrue<-stage7m_truth_functions(tr,grid);tbas<-stage7m_truth_basis(tr)
 # The main grid is the prespecified full [0,1] interval, with no truth-driven trimming.
 raw<-stage7m_grid_estimate(fit,fixed$basis,grid,'raw',reported);canonical<-stage7m_grid_estimate(fit,fixed$basis,grid,'reported',reported)
 functions<-historical$primary$functional_recovery
 functions<-functions[functions$scope=='ppi_selected',,drop=FALSE]
 subject<-list();postdiff<-list();predrows<-list();roles<-list();components<-list();scorepoints<-list();loadingpoints<-list();loadrows<-list();maps<-list()
 append_metric<-function(study,person,scope,representation,object,est,target,weights,role_status='all_estimated_blocks') {
  er<-stage7m_curve_error(est,target,weights)
  subject[[length(subject)+1]]<<-data.frame(study=study,subject=person,scope=scope,representation=representation,object=object,role_status=role_status,mse=er['mse'],truth_energy=er['truth_energy'],rmse=er['rmse'],nrmse=er['nrmse'],relative_reason=if(er['truth_energy']==0)'zero_truth_energy'else'')
 }
 for(s in seq_len(fit$S))for(i in seq_len(fit$n_s[s])) {
  target<-stage7m_true_parts(tr,obs,s,i,grid,tbas,ftrue)
  a<-stage7m_grid_subject(raw,obs,s,i);b<-stage7m_grid_subject(canonical,obs,s,i)
  for(rep in c('raw','reported')) {
   est<-if(rep=='raw')a else b
   for(nm in c('signal','factor','mean','beta'))append_metric(s,i,'grid201',rep,nm,est[[nm]],target[[nm]],w)
  }
  postdiff[[length(postdiff)+1]]<-data.frame(study=s,subject=i,scope='grid201',raw_vs_reported_rmse=sqrt(sum(w*(a$signal-b$signal)^2)/fit$p),raw_vs_reported_max_abs=max(abs(a$signal-b$signal)))
  tobs<-obs$time_obs[[s]][[i]];observed_target<-list(signal=tr$signal_true_values[[s]][[i]],mean=tr$mu_true_values[[s]][[i]],beta=tr$beta_true_values[[s]][[i]],factor=tr$f_true_values[[s]][[i]]%*%t(tr$a_true)+tr$g_true_values[[s]][[i]]%*%t(tr$b_true[[s]]))
  rawobs<-reportobs<-NULL
  for(rep in c('raw','reported')) {
   e<-stage7m_estimate(fit,fixed$basis,s,i,tobs,rep,reported);bc<-matrix(0,length(tobs),fit$p)
   for(r in seq_len(fit$d))bc<-bc+obs$Z[[s]][i,r]*e$beta[[r]]
   e$beta<-bc;e$signal<-e$mean+bc+e$factor
   if(rep=='raw')rawobs<-e$signal else reportobs<-e$signal
   for(nm in c('signal','factor','mean','beta'))append_metric(s,i,'observed',rep,nm,e[[nm]],observed_target[[nm]],rep(1/length(tobs),length(tobs)))
  }
  postdiff[[length(postdiff)+1]]<-data.frame(study=s,subject=i,scope='observed',raw_vs_reported_rmse=sqrt(mean((rawobs-reportobs)^2)),raw_vs_reported_max_abs=max(abs(rawobs-reportobs)))
  for(rep in c('raw','reported')) {
   yp<-pred[[rep]][[s]][[i]];ys<-test_truth$signal[[s]][[i]];yn<-test_truth$Y[[s]][[i]]
   predrows[[length(predrows)+1]]<-data.frame(study=s,subject=i,method=paste0('multiFSYNC_',rep),noisy_mse=mean((yp-yn)^2),noisy_mae=mean(abs(yp-yn)),signal_mse=mean((yp-ys)^2),signal_energy=mean(ys^2),noisy_rmse=sqrt(mean((yp-yn)^2)),signal_rmse=sqrt(mean((yp-ys)^2)),signal_nrmse=sqrt(sum((yp-ys)^2)/sum(ys^2)))
  }
 }
 for(k in seq_len(nrow(functions))) {
  item<-functions[k,,drop=FALSE];s<-as.integer(item$study);role<-item$true_role;l<-as.integer(item$true_factor);matched<-isTRUE(item$matched_by_loading);block<-item$estimated_block;el<-as.integer(item$estimated_factor)
  truthload<-if(role=='shared')tr$a_true[,l]else tr$b_true[[s]][,l]
  truefun<-if(role=='shared')ftrue$shared[[l]]else ftrue$specific[[s]][[l]]
  truescores<-if(role=='shared')tr$zeta_fpca_true[[s]][[l]]else tr$xi_fpca_true[[s]][[l]]
  trueeig<-if(role=='shared')tr$lambda_phi_true[[l]]else tr$lambda_psi_true[[s]][[l]]
  estload<-rawload<-rep(0,fit$p);ppi<-rep(0,fit$p);loadsign<-rawsign<-1
  if(matched) {
   estload<-if(block=='shared')fit$mu_q_a_hat[,el]else fit$mu_q_b_specific_hat[[s]][,el]
   oldl<-if(block=='shared')fit$factor_order_shared[el]else fit$factor_order_specific[[s]][el]
   rawload<-if(block=='shared')fit$mu_q_a_original[,oldl]else fit$mu_q_b_specific_original[[s]][,oldl]
   ppi<-if(block=='shared')fit$mu_q_gamma_a_hat[,el]else fit$mu_q_gamma_b_hat[[s]][,el]
   loadsign<-if(sum(estload*truthload)<0)-1 else 1;rawsign<-if(sum(rawload*truthload)<0)-1 else 1
   estfun<-if(block=='shared')canonical$shared[[el]]else canonical$specific[[s]][[el]]
   estscores<-if(block=='shared')canonical$scores_shared[[s]][[el]]else canonical$scores_specific[[s]][[el]]
   esteig<-if(block=='shared')fit$list_eigenvalues[[el]]else fit$list_eigenvalues_spec[[s]][[el]]
   map<-stage7m_component_map(truefun,estfun,w)
  } else {oldl<-NA_integer_;estfun<-matrix(0,length(grid),2);estscores<-matrix(0,fit$n_s[s],2);esteig<-c(0,0);map<-list(assignment=c(0,0),function_sign=c(1,1))}
  maprow<-data.frame(study=s,true_role=role,true_factor=l,estimated_block=if(matched)block else NA_character_,estimated_factor=if(matched)el else NA_integer_,raw_factor_index=oldl,matched_by_loading=matched,loading_status=item$loading_structure_status,loading_sign=loadsign,raw_loading_sign=rawsign,shared_alignment_scope=if(role=='shared'&&matched&&block=='shared')'one_common_function_map_across_both_studies'else if(role=='shared')'shared_role_missing_or_misplaced'else'study_specific_function_map')
  maps[[length(maps)+1]]<-maprow
  auc<-stage7m_auc(ppi,truthload!=0)
  for(rep in c('raw_original','canonical_FPCA')) {
   est<-if(rep=='raw_original')rawload else estload;sig<-if(rep=='raw_original')rawsign else loadsign
   loadrows[[length(loadrows)+1]]<-cbind(maprow,representation=rep,data.frame(loading_cosine=abs(stage7m_cos(est,truthload)),loading_norm_ratio=sqrt(sum(est^2)/sum(truthload^2)),loading_relative_l2=sqrt(sum((sig*est-truthload)^2)/sum(truthload^2))),auc)
  }
  loadingpoints[[length(loadingpoints)+1]]<-data.frame(study=s,true_role=role,true_factor=l,variable=seq_len(fit$p),true_loading=truthload,raw_original_loading_aligned=rawsign*rawload,canonical_loading_aligned=loadsign*estload,ppi=ppi,matched_by_loading=matched,loading_status=item$loading_structure_status)
  for(i in seq_len(fit$n_s[s])) {
   path<-if(matched)loadsign*as.vector(estfun%*%estscores[i,])else rep(0,length(grid))
   truepath<-as.vector(truefun%*%truescores[i,]);contrib<-if(matched)outer(as.vector(estfun%*%estscores[i,]),estload)else matrix(0,length(grid),fit$p)
   for(nm in c('process','contribution')) {
    est<-if(nm=='process')matrix(path,ncol=1)else contrib;target<-if(nm=='process')matrix(truepath,ncol=1)else outer(truepath,truthload)
    er<-stage7m_curve_error(est,target,w)
    roles[[length(roles)+1]]<-cbind(maprow,subject=i,object=nm,data.frame(stage7m_mse=er['mse'],stage7m_truth_energy=er['truth_energy'],stage7m_rmse=er['rmse'],stage7m_nrmse_missing_as_zero=er['nrmse']))
   }
  }
  for(m in 1:2) {
   j<-map$assignment[m];has<-matched&&j>0;fs<-map$function_sign[m]
   ef<-if(has)fs*estfun[,j]else rep(0,length(grid));es<-if(has)loadsign*fs*estscores[,j]else rep(0,fit$n_s[s]);ee<-if(has)esteig[j]else 0
   fe<-stage7m_curve_error(matrix(ef,ncol=1),matrix(truefun[,m],ncol=1),w)
   scoremetrics<-stage7m_score_error(es,truescores[,m])
   components[[length(components)+1]]<-cbind(maprow,true_component=m,estimated_component=if(has)j else NA_integer_,component_sign=fs,component_matched=has,data.frame(stage7m_function_ise=fe['mse'],stage7m_function_relative_l2=fe['nrmse'],stage7m_function_norm_ratio=sqrt(sum(w*ef^2)/sum(w*truefun[,m]^2)),stage7m_function_abs_inner=abs(sum(w*ef*truefun[,m])),stage7m_function_cosine=abs(stage7m_cos(ef,truefun[,m],w)),stage7m_true_eigenvalue=trueeig[m],stage7m_estimated_eigenvalue=ee,stage7m_eigenvalue_error=ee-trueeig[m],stage7m_true_pve=trueeig[m]/sum(trueeig),stage7m_estimated_pve=if(sum(esteig)>0)ee/sum(esteig)else 0),scoremetrics)
   scorepoints[[length(scorepoints)+1]]<-data.frame(study=s,true_role=role,true_factor=l,true_component=m,subject=seq_len(fit$n_s[s]),true_score=truescores[,m],estimated_score_aligned=es,component_matched=has,loading_status=item$loading_structure_status)
  }
 }
 variable_metrics<-list()
 for(s in seq_len(fit$S)) {
  target<-ftrue$mean[[s]];estimate<-raw$mean[[s]]
  variable_metrics[[length(variable_metrics)+1]]<-data.frame(object='mean',study=s,covariate=NA_integer_,variable=seq_len(fit$p),stage7m_ise=colSums(w*(estimate-target)^2),stage7m_truth_energy=colSums(w*target^2),stage7m_estimated_energy=colSums(w*estimate^2),relative_reason='')
 }
 for(r in seq_len(fit$d)){target<-ftrue$beta[[r]];estimate<-raw$beta[[r]];variable_metrics[[length(variable_metrics)+1]]<-data.frame(object='beta',study=NA_integer_,covariate=r,variable=seq_len(fit$p),stage7m_ise=colSums(w*(estimate-target)^2),stage7m_truth_energy=colSums(w*target^2),stage7m_estimated_energy=colSums(w*estimate^2),relative_reason=if(identical(tr$beta_truth_method,'zero_effect_d2'))'zero_truth_energy'else'')}
 variable_metrics<-do.call(rbind,variable_metrics);variable_metrics$stage7m_relative_ise<-ifelse(variable_metrics$stage7m_truth_energy>0,variable_metrics$stage7m_ise/variable_metrics$stage7m_truth_energy,NA_real_)
 subjects<-do.call(rbind,subject);prediction<-do.call(rbind,predrows);summary<-data.frame()
 for(rep in c('raw','reported'))for(scope in c('grid201','observed'))for(object in c('signal','factor','mean','beta')) {
  z<-subjects[subjects$representation==rep & subjects$scope==scope & subjects$object==object,,drop=FALSE];stem<-paste('stage7m',object,scope,rep,sep='_')
  summary[1,paste0(stem,'_rmse')]<-sqrt(mean(z$mse));summary[1,paste0(stem,'_nrmse')]<-if(sum(z$truth_energy)>0)sqrt(sum(z$mse)/sum(z$truth_energy))else NA_real_
 }
 for(rep in c('raw','reported')){z<-prediction[prediction$method==paste0('multiFSYNC_',rep),];summary[1,paste0('stage7m_test_response_rmse_',rep)]<-sqrt(mean(z$noisy_mse));summary[1,paste0('stage7m_test_response_mae_',rep)]<-mean(z$noisy_mae);summary[1,paste0('stage7m_test_signal_rmse_',rep)]<-sqrt(mean(z$signal_mse));summary[1,paste0('stage7m_test_signal_nrmse_',rep)]<-sqrt(sum(z$signal_mse)/sum(z$signal_energy))}
 beta<-variable_metrics[variable_metrics$object=='beta',];summary$stage7m_beta_ise_mean<-mean(beta$stage7m_ise);summary$stage7m_beta_estimated_energy_mean<-mean(beta$stage7m_estimated_energy);summary$stage7m_beta_relative_ise<-if(sum(beta$stage7m_truth_energy)>0)sum(beta$stage7m_ise)/sum(beta$stage7m_truth_energy)else NA_real_
 pd<-do.call(rbind,postdiff);summary$stage7m_postprocess_grid_rmse<-sqrt(mean(pd$raw_vs_reported_rmse[pd$scope=='grid201']^2));summary$stage7m_postprocess_observed_max_abs<-max(pd$raw_vs_reported_max_abs[pd$scope=='observed'])
 list(summary=summary,subject_metrics=subjects,role_metrics=do.call(rbind,roles),component_metrics=do.call(rbind,components),score_points=do.call(rbind,scorepoints),loading_points=do.call(rbind,loadingpoints),loading_metrics=do.call(rbind,loadrows),alignment_maps=do.call(rbind,maps),variable_metrics=variable_metrics,prediction_metrics=prediction,postprocess_metrics=pd,reported_coefficients=reported)
}
stage7m_export_nested<-function(x,root,prefix='legacy',index=list()) {
 if(is.data.frame(x)) {
  path<-file.path(root,paste0(prefix,'.csv'));new_csv(x,path)
  return(list(data.frame(object=prefix,rows=nrow(x),columns=ncol(x),file=basename(path),sha256=sha(path))))
 }
 if(is.list(x)) {
  if(is.null(names(x)))names(x)<-sprintf('item%03d',seq_along(x))
  for(nm in names(x)) index<-c(index,stage7m_export_nested(x[[nm]],root,paste(prefix,nm,sep='__')))
 }
 index
}
