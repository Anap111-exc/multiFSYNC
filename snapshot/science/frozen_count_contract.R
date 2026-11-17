# Pure no-truth endpoint gates for the registered M=3 factor-count pilot.
# Existing M=3 helpers remain unchanged; this wrapper is specific to this study.
fc_valid_status <- function(z) {
  req <- c('status','objective_valid','finite','metadata_pass','monotonic_pass','ELBO')
  iv_assert(nrow(z)==1L&&all(req%in%names(z)),'Malformed terminal fit status')
  identical(as.character(z$status),'complete')&&isTRUE(z$objective_valid)&&
    isTRUE(z$finite)&&isTRUE(z$metadata_pass)&&isTRUE(z$monotonic_pass)&&is.finite(z$ELBO)
}
fc_check_endpoint <- function(fit,fixed,seed,elbo) {
  iv_assert(fit$S==2L&&fit$p==500L&&all(fit$n_s==30L)&&fit$d==0L&&
    fit$L_f==2L&&all(fit$L_s_by_study==2L)&&all(fit$M_f==3L)&&all(unlist(fit$M_s)==3L)&&
    fit$K==9L&&fit$K_total==11L&&fit$t1_sweeps==400L&&fit$annealing_sweeps==99L&&
    fit$i_iter==499L&&length(fit$ELBO)==400L&&isTRUE(fit$elbo_objective_valid)&&
    identical(as.integer(fit$fit_seed),as.integer(seed)),'Pilot endpoint design/counters differ')
  iv_assert(ae_close(tail(fit$ELBO,1L),elbo,1e-12)&&ae_monotonic_pass(fit$ELBO),
    'Complete ordinary ELBO/history differs')
  fields<-c('S','n_s','p','d','K','K_total','C_g','time_g')
  iv_assert(all(vapply(fields,function(nm)ae_close(fit[[nm]],fixed$fit_metadata[[nm]]),logical(1))),
    'Frozen observation metadata differs')
  fields<-c('mu_q_nu_mu','mu_q_nu_phi','mu_q_nu_psi','mu_q_zeta','mu_q_xi',
    'mu_q_a','mu_q_b_specific','mu_q_gamma_a','mu_q_gamma_b','mu_q_recip_sigsq_eps','ELBO')
  iv_assert(all(fields%in%names(fit))&&all(is.finite(unlist(fit[fields]))),'Nonfinite/missing endpoint fields')
  z<-fit$factor_shrinkage_state
  iv_assert(!is.null(z)&&z$version==1L&&z$control$mode=='mgp'&&z$control$a1==2.1&&
    z$control$a2==3.1&&z$control$rate==1&&z$p==500L&&z$L_f==2L&&all(z$L_s==2L),
    'Frozen MGP state differs')
  ae_tau_moments(fit)
  invisible(TRUE)
}
fc_basis_ppi_check <- function(fit,obs,fixed,legacy) {
  iv_assert(is.null(obs$Z)&&obs$dimensions$d==0L&&isFALSE(obs$truth_available_to_fit)&&
    isFALSE(obs$test_response_available_to_fit),'Observation truth contract differs')
  err<-max(abs(legacy$stage7m_C(fit$time_g,fixed$basis)-fit$C_g))
  ce<-max(vapply(seq_len(fit$S),function(s)max(vapply(seq_len(fit$n_s[s]),function(i)
    max(abs(legacy$stage7m_C(obs$time_obs[[s]][[i]],fixed$basis)-fixed$C[[s]][[i]])),numeric(1))),numeric(1)))
  iv_assert(err<=1e-12&&ce<=1e-12,'Frozen Linux basis/time design differs')
  ae_rank_to_raw(fit$factor_order_shared,fit$list_effective_M,fit$L_f)
  for(s in seq_len(fit$S))ae_rank_to_raw(fit$factor_order_specific[[s]],fit$list_effective_M_spec[[s]],fit$L_s_by_study[s])
  rows<-list()
  for(s in 0:fit$S) {
    gam<-if(s==0L)fit$mu_q_gamma_a else fit$mu_q_gamma_b[[s]]
    ord<-if(s==0L)fit$factor_order_shared else fit$factor_order_specific[[s]]
    saved<-if(s==0L)fit$factor_ppi_shared else fit$factor_ppi_specific[[s]]
    stable<-apply(gam,2,ae_factor_ppi)[ord];direct<-apply(gam,2,function(x)1-prod(1-x))[ord]
    iv_assert(length(saved)==length(stable)&&all(is.finite(c(saved,direct,stable)))&&
      max(abs(saved-direct))==0&&max(abs(saved-stable))<=1e-12,'PPI reporting roundoff')
    iv_assert(identical(saved>.5,stable>.5)&&identical(saved>=.5,stable>=.5),'PPI classification differs')
    rows[[s+1L]]<-data.frame(study=s,live_gt_half=sum(stable>.5),formal_ge_half=sum(stable>=.5),
      exact_half=sum(stable==.5),max_ppi_absolute_error=max(abs(saved-stable)))
  }
  list(basis=data.frame(C_g_max_absolute_error=err,C_max_absolute_error=ce),ppi=do.call(rbind,rows))
}
fc_splinets <- function(ex) {
  dep<-jsonlite::read_json(file.path(ex,'generation/DEPENDENCY_BINDING.json'),simplifyVector=FALSE)
  iv_assert(Sys.info()[['sysname']]=='Linux'&&as.character(getRversion())=='4.5.3','Registered Linux R 4.5.3 required')
  .libPaths(unique(c(dirname(dep$package_root),dep$R_LIBS_USER,.libPaths())))
  iv_assert(as.character(utils::packageVersion('Splinets'))=='1.5.1','Frozen Splinets version differs')
  iv_assert(identical(normalizePath(getNamespaceInfo(asNamespace('Splinets'),'path'),mustWork=TRUE),
    normalizePath(dep$package_root,mustWork=TRUE)),'Frozen Splinets path differs')
  for(entry in dep$package_files)iv_check(file.path(dep$package_root,entry$relative_path),entry$sha256)
  invisible(dep)
}
fc_choose <- function(attempts) {
  iv_assert(nrow(attempts)==12L&&!anyDuplicated(attempts$fit_id)&&setequal(attempts$start,1:12),
    'Selection requires all 12 registered terminal attempts')
  valid<-attempts[attempts$valid,,drop=FALSE]
  if(!nrow(valid))return(list(status='no_valid_candidate',selected_fit_id='',selected_ELBO=NULL,
    tie_fit_ids=character(),tie_size=0L,tie_rule='exact_equal_ELBO_then_smallest_start_then_fit_id'))
  best<-max(valid$ELBO);ties<-valid[valid$ELBO==best,,drop=FALSE]
  ties<-ties[order(ties$start,ties$fit_id),,drop=FALSE]
  list(status='selected',selected_fit_id=as.character(ties$fit_id[1L]),selected_ELBO=best,
    tie_fit_ids=as.character(ties$fit_id),tie_size=nrow(ties),
    tie_rule='exact_equal_ELBO_then_smallest_start_then_fit_id')
}
