# Evaluation-only multi-factor adapter; frozen NC thresholds and estimands unchanged.
# All six global truth roles remain in the denominator, including absent roles.
science_match <- function(sim, live, threshold) {
  z <- as.matrix(sim)
  z[,!live] <- -Inf
  nc_one_to_one_match_pruned(z,threshold)$assignment
}

ae_raw_metrics <- function(fit,truth,obs,fixed) {
  ranka<-ae_rank_to_raw(fit$factor_order_shared,fit$list_effective_M,fit$L_f)
  rankb<-lapply(seq_len(fit$S),function(s)ae_rank_to_raw(fit$factor_order_specific[[s]],fit$list_effective_M_spec[[s]],fit$L_s_by_study[s]))
  tau<-ae_tau_moments(fit);factors<-roles<-energies<-list()
  targets<-c(list(as.matrix(truth$a_true)),lapply(truth$b_true,as.matrix))
  alltargets<-do.call(cbind,targets)
  global_block<-rep(0:fit$S,vapply(targets,ncol,integer(1)))
  global_factor<-unlist(lapply(targets,function(x)seq_len(ncol(x))))
  counts<-counts_gt<-integer(fit$S+1L)
  for(blockid in 0:fit$S) {
    sh<-blockid==0L
    load<-if(sh)fit$mu_q_a_original else fit$mu_q_b_specific_original[[blockid]]
    gam<-if(sh)fit$mu_q_gamma_a else fit$mu_q_gamma_b[[blockid]]
    target<-targets[[blockid+1L]];nt<-ncol(target);nc<-ncol(load)
    ranks<-if(sh)ranka else rankb[[blockid]]
    tm<-if(sh)tau$shared else tau$specific[[blockid]]
    relevant<-if(sh)seq_len(fit$S)else blockid
    ppi<-vapply(seq_len(nc),function(l)ae_factor_ppi(gam[,l]),numeric(1))
    live<-ppi>=.5;counts[blockid+1L]<-sum(live);counts_gt[blockid+1L]<-sum(ppi>.5)
    sim<-matrix(0,nt,nc)
    for(k in seq_len(nt))for(l in seq_len(nc))sim[k,l]<-ae_cosine(load[,l],target[,k])
    sim[!is.finite(sim)]<-0
    assignment09<-science_match(sim,live,.9);assignment08<-science_match(sim,live,.8)
    trajectory<-array(NA_real_,c(nt,nc,length(relevant)))
    for(l in seq_len(nc)) {
      cosall<-vapply(seq_len(ncol(alltargets)),function(k)ae_cosine(load[,l],alltargets[,k]),numeric(1))
      best<-if(any(is.finite(cosall)))which.max(replace(cosall,!is.finite(cosall),-Inf))else NA_integer_
      bestown<-which.max(sim[,l]);bes<-list()
      for(si in seq_along(relevant)) {
        s<-relevant[si]
        coef<-if(sh)fit$mu_q_nu_phi[[l]]else fit$mu_q_nu_psi[[s]][[l]]
        scores<-if(sh)fit$mu_q_zeta[[s]][[l]]else fit$mu_q_xi[[s]][[l]]
        path<-lapply(seq_len(fit$n_s[s]),function(i)as.vector(fixed$C[[s]][[i]]%*%coef%*%scores[i,]))
        targetpath<-if(sh)truth$f_true_values[[s]]else truth$g_true_values[[s]]
        for(k in seq_len(nt)) {
          te<-unlist(lapply(targetpath,function(x)as.matrix(x)[,k]),use.names=FALSE)
          trajectory[k,l,si]<-sign(sum(load[,l]*target[,k]))*ae_signed_correlation(unlist(path,use.names=FALSE),te)
        }
        pe<-mean(vapply(path,function(x)mean(x^2),numeric(1)))
        bes[[si]]<-data.frame(block=if(sh)'shared'else'specific',block_study=blockid,raw_factor=l,evaluation_study=s,
          process_mean_energy_equal_subject=pe,contribution_mean_energy_equal_subject=pe*mean(load[,l]^2),
          function_mean_energy_equal_subject=mean(vapply(seq_len(fit$n_s[s]),function(i)mean(rowSums((fixed$C[[s]][[i]]%*%coef)^2)),numeric(1))),
          score_mean_squarednorm=mean(rowSums(scores^2)))
      }
      es<-do.call(rbind,bes);energies[[length(energies)+1L]]<-es
      factors[[length(factors)+1L]]<-data.frame(block=if(sh)'shared'else'specific',study=blockid,raw_factor=l,
        reported_factor=match(l,if(sh)fit$factor_order_shared else fit$factor_order_specific[[blockid]]),
        factor_ppi=ppi[l],live_ppi_ge_half=live[l],live_ppi_gt_half=ppi[l]>.5,
        expected_support=sum(gam[,l]),loading_ppi_gt_half_count=sum(gam[,l]>.5),loading_mean_norm_original=sqrt(sum(load[,l]^2)),
        expected_tau=tm$mean[l],expected_log_tau=tm$elog[l],fpca_effective_M=ranks[l],
        best_own_true_factor=bestown,own_role_loading_abs_cosine=sim[bestown,l],
        own_role_loading_core_09=live[l]&&sim[bestown,l]>=.9,own_role_loading_core_08=live[l]&&sim[bestown,l]>=.8,
        mean_contribution_energy=mean(es$contribution_mean_energy_equal_subject),min_study_contribution_energy=min(es$contribution_mean_energy_equal_subject),
        best_true_role_study=if(is.na(best))NA_integer_ else global_block[best],
        best_true_factor=if(is.na(best))NA_integer_ else global_factor[best],
        best_true_role_is_offblock=!is.na(best)&&global_block[best]!=blockid)
    }
    for(k in seq_len(nt)) {
      chosen<-assignment09[k];present<-chosen>0L
      support<-target[,k]!=0;selected<-if(present)gam[,chosen]>.5 else rep(FALSE,nrow(target))
      tp<-sum(selected&support);fp<-sum(selected&!support);fn<-sum(!selected&support)
      trcor<-if(present)trajectory[k,chosen,]else NA_real_
      roles[[length(roles)+1L]]<-data.frame(block=if(sh)'shared'else'specific',study=blockid,true_factor=k,
        truth_factor_count=nt,fitted_capacity=nc,live_factor_count=sum(live),
        loading_core_present_09=present,loading_core_present_08=assignment08[k]>0L,
        duplicate_own_role_core_count_09=max(0L,sum(live&sim[k,]>=.9)-1L),
        best_own_role_loading_abs_cosine=max(sim[k,]),
        matched_raw_factor_09=if(present)chosen else NA_integer_,
        matched_fpca_M_09=if(present)ranks[chosen]else NA_integer_,
        matched_fpca_exact_09=present&&ranks[chosen]==2L,
        loading_core_and_trajectory08_present=present&&all(is.finite(trcor))&&all(trcor>=.8),
        support_precision_ppi_gt_half=if(tp+fp>0)tp/(tp+fp)else 0,
        support_recall_ppi_gt_half=if(any(support))tp/sum(support)else NA_real_,
        support_true_positive=tp,support_false_positive=fp,support_false_negative=fn)
    }
  }
  rr<-do.call(rbind,roles);truecounts<-as.integer(c(truth$L_f,truth$L_s_by_study))
  ae_assert(nrow(rr)==sum(truecounts),'Every true factor must have exactly one evaluation row')
  summary<-data.frame(true_global_roles=sum(truecounts),matched_global_roles_09=sum(rr$loading_core_present_09),
    loading_role_core_09=all(rr$loading_core_present_09),loading_role_core_08=all(rr$loading_core_present_08),
    ppi_counts_ge_half_exact=all(counts==truecounts),ppi_counts_gt_half_exact=all(counts_gt==truecounts),
    legacy_PPI_exact_09=all(rr$loading_core_present_09)&&all(counts==truecounts),
    loading_core_and_trajectory08_all_roles=all(rr$loading_core_and_trajectory08_present),
    matched_fpca_all_roles_exact_09=all(rr$matched_fpca_exact_09),
    live_shared=counts[1],live_specific_1=counts[2],live_specific_2=counts[3],
    fpca_ranks_raw_shared=paste(ranka,collapse=';'),fpca_ranks_raw_specific_1=paste(rankb[[1]],collapse=';'),fpca_ranks_raw_specific_2=paste(rankb[[2]],collapse=';'))
  list(summary=summary,factors=do.call(rbind,factors),roles=rr,energies=do.call(rbind,energies))
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
  exactcounts<-identical(as.integer(counts),as.integer(c(truth$L_f,truth$L_s_by_study)))
  exact<-core&&exactcounts&&!any(gl$candidate$status%in%c('extra','duplicate'))&&!any(cc$status%in%c('extra','duplicate'))
  sumrow<-data.frame(formal_nc_loading_core_08=core,formal_nc_exact_loading_structure_08=exact,formal_nc_ppi_ge_half_counts_exact=exactcounts,
    formal_true_global_roles=truth$L_f+sum(truth$L_s_by_study),formal_global_correct=sum(gl$truth$status=='correct'),formal_global_missing=sum(gl$truth$status=='missing'),formal_global_misplaced=sum(gl$truth$status=='misplaced'),formal_global_duplicate=sum(gl$candidate$status=='duplicate'),formal_global_extra=sum(gl$candidate$status=='extra'),
    formal_study_correct=sum(tt$status=='correct'),formal_study_missing=sum(tt$status=='missing'),formal_study_misplaced=sum(tt$status=='misplaced'),formal_study_duplicate=sum(cc$status=='duplicate'),formal_study_extra=sum(cc$status=='extra'))
  # These are exactly the formal catalog/one-to-one loading definitions. They
  # do not constitute a function/trajectory/score or full-basin recovery claim.
  paths<-legacy$v0g_reported_paths_by_study(fit,data)
  shared<-legacy$nc_v0g_shared_candidate_metrics(fit,paths,.10)
  shared$raw_factor<-fit$factor_order_shared[shared$candidate_factor]
  block<-legacy$v0h_block_subspace_metrics(fit,data,'ppi_selected',legacy$v0h_default_thresholds())
  list(summary=sumrow,global_truth=gl$truth,global_candidate=gl$candidate,study_truth=tt,study_candidate=cc,shared_support=shared,block_subspace=block)
}
