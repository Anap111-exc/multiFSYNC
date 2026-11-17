# Evaluation adapters only. Frozen numerical helpers remain byte-identical.
# No inference update, threshold change or truth-based winner selection.
pair_load <- function(ex) {
  source(file.path(ex,'evaluation/io.R'),local=.GlobalEnv)
  source(file.path(ex,'evaluation/frozen/evaluate_aligned.R'),local=.GlobalEnv)
  source(file.path(ex,'evaluation/endpoint.R'),local=.GlobalEnv)
  source(file.path(ex,'science/frozen_count_contract.R'),local=.GlobalEnv)
  source(file.path(ex,'science/frozen_adaptation_science.R'),local=.GlobalEnv)
  env <- adapt_load_evaluator(ex)
  source(file.path(ex,'science/generalized_roles.R'),local=.GlobalEnv)
  # Extend the report-only two-component matching to zero through four estimates.
  assign('pair_legacy_component_map',stage7m_component_map,envir=.GlobalEnv)
  assign('stage7m_component_map',pair_component_map,envir=.GlobalEnv)
  env
}

pair_component_map <- function(true,estimate,w) {
  mt<-ncol(true);me<-ncol(estimate)
  stopifnot(mt==2L,me%in%0:4,length(w)==nrow(true),nrow(estimate)==nrow(true))
  if(me%in%c(1L,2L))return(pair_legacy_component_map(true,estimate,w))
  sim<-matrix(0,mt,me)
  for(i in seq_len(mt))for(j in seq_len(me)) {
    v<-stage7m_cos(true[,i],estimate[,j],w);sim[i,j]<-if(is.finite(v))abs(v)else 0
  }
  assignment<-c(0L,0L)
  if(me==1L)assignment[which.max(sim[,1])]<-1L
  if(me>=2L) {
    # Lexical exhaustive enumeration reproduces the original M2 exact-tie rule.
    best<--Inf
    for(a in seq_len(me))for(b in seq_len(me))if(a!=b) {
      score<-sim[1,a]+sim[2,b]
      if(score>best){best<-score;assignment<-as.integer(c(a,b))}
    }
  }
  fs<-rep(1,mt)
  for(i in seq_len(mt))if(assignment[i]>0) {
    val<-sum(w*true[,i]*estimate[,assignment[i]]);fs[i]<-if(val<0)-1 else 1
  }
  list(assignment=assignment,function_sign=fs,similarity=sim)
}

pair_check_endpoint <- function(fit,fixed,seed,elbo,arm) {
  iv_assert(identical(arm,'original'),'Unknown arm')
  iv_assert(fit$S==2L&&fit$p==500L&&all(fit$n_s==60L)&&fit$d==0L&&
    fit$L_f==4L&&all(fit$L_s_by_study==4L)&&all(fit$M_f==4L)&&all(unlist(fit$M_s)==4L)&&
    fit$K==9L&&fit$K_total==11L&&fit$t1_sweeps==400L&&fit$annealing_sweeps==99L&&
    fit$i_iter==499L&&length(fit$ELBO)==400L&&isTRUE(fit$elbo_objective_valid)&&
    identical(as.integer(fit$fit_seed),as.integer(seed)),'L444/M4 endpoint design/counters differ')
  if(arm=='original') {
    iv_assert(is.null(fit$factor_shrinkage_state)&&is.null(fit$factor_shrinkage),
      'Original CAVI contains an MGP state or control')
  }
  iv_assert(ae_close(tail(fit$ELBO,1),elbo,1e-12)&&ae_monotonic_pass(fit$ELBO),'Ordinary objective/history differs')
  fields<-c('S','n_s','p','d','K','K_total','C_g','time_g')
  iv_assert(all(vapply(fields,function(nm)ae_close(fit[[nm]],fixed$fit_metadata[[nm]]),logical(1))),
    'Frozen observation metadata differs')
  fields<-c('mu_q_nu_mu','mu_q_nu_phi','mu_q_nu_psi','mu_q_zeta','mu_q_xi',
    'mu_q_a','mu_q_b_specific','mu_q_gamma_a','mu_q_gamma_b','mu_q_recip_sigsq_eps','ELBO')
  iv_assert(all(fields%in%names(fit))&&all(is.finite(unlist(fit[fields]))),'Invalid endpoint fields')
  invisible(TRUE)
}

pair_historical <- function(fit,obs,fixed,truth,env) {
  assert(fit$S==2L&&fit$L_f==4L&&all(unlist(fit$L_s)==4L)&&all(fit$M_f==4L)&&
    all(unlist(fit$M_s)==4L),'This evaluation is registered for S2/L444/M4')
  assert(fit$d==0L&&truth$data$true_params$d_use==0L,'Unexpected covariate estimand')
  assert(max(abs(stage7m_C(fit$time_g,fixed$basis)-fit$C_g))<=1e-12,'Frozen basis differs')
  old<-env$estimated_mean_dense
  on.exit(assign('estimated_mean_dense',old,envir=env),add=TRUE)
  env$estimated_mean_dense<-function(fit,data,study) {
    stage7m_C(data$true_params$t_grid_dense,fixed$basis)%*%do.call(cbind,fit$mu_q_nu_mu_original[[study]])
  }
  s6h_evaluate_fit(fit,truth,env)
}

pair_supplement <- function(fit,obs,fixed,truth) {
  tr<-truth$data$true_params;t<-fixed$grid;w<-fixed$weights
  tb<-stage7m_truth_basis(tr);tf<-stage7m_truth_functions(tr,t,tb)
  a<-which(rowSums(abs(as.matrix(tr$a_true)))>0)
  b<-lapply(tr$b_true,function(x)which(rowSums(abs(as.matrix(x)))>0))
  support<-lapply(seq_len(fit$S),function(s)list(shared_contribution=a,specific_contribution=b[[s]],
    factor=sort(union(a,b[[s]])),signal=sort(union(a,b[[s]]))))
  pack<-list(estimate=fit,observation=obs,basis=fixed$basis);cache<-band_cache(pack,t)
  qi<-cr<-list()
  for(s in seq_len(fit$S)) {
    eg<-tg<-vector('list',fit$n_s[s])
    for(i in seq_len(fit$n_s[s])) {
      x<-band_subject(pack,cache,s,i);target<-stage7m_true_parts(tr,obs,s,i,t,tb,tf)
      eg[[i]]<-x$mean$signal;tg[[i]]<-target$signal
      for(object in names(support[[s]])) {
        cols<-support[[s]][[object]];assert(length(cols)>0,'Empty active support')
        met<-band_summary(x$mean[[object]],x$variance[[object]],target[[object]],w,qnorm(.975),cols)
        qi[[length(qi)+1L]]<-data.frame(study=s,subject=i,scope='grid201',object=object,
          mask='active',level=.95,n_variables=length(cols),weight=1/fit$S/fit$n_s[s],t(met),check.names=FALSE)
      }
    }
    em<-Reduce(`+`,eg)/length(eg);tm<-Reduce(`+`,tg)/length(tg)
    for(i in seq_along(eg)) {
      er<-stage7m_curve_error(eg[[i]]-em,tg[[i]]-tm,w)
      cr[[length(cr)+1L]]<-data.frame(study=s,subject=i,scope='grid201',object='within_study_centered_signal',
        mse=unname(er['mse']),truth_energy=unname(er['truth_energy']),weight=1/fit$S/fit$n_s[s])
    }
  }
  qi<-do.call(rbind,qi);cr<-do.call(rbind,cr)
  qsum<-do.call(rbind,lapply(unique(qi$object),function(object) {
    x<-qi[qi$object==object,];data.frame(scope='grid201',object=object,mask='active',level=.95,
      coverage=sum(x$weight*x$coverage),mean_width=sum(x$weight*x$width),
      all_grid_covered=sum(x$weight*x$all_grid_covered))
  }))
  mse<-sum(cr$weight*cr$mse);energy<-sum(cr$weight*cr$truth_energy)
  center<-data.frame(scope='grid201',object='within_study_centered_signal',mse=mse,
    truth_energy=energy,rmse=sqrt(mse),nrmse=sqrt(mse/energy))
  tv<-as.matrix(tr$sigma2_eps_true);ev<-as.matrix(fit$sigsq_eps_hat)
  assert(identical(dim(tv),dim(ev))&&all(ev>0)&&all(tv>0),'Noise dimension/positivity')
  assert(max(abs(ev-sweep(fit$sigsq_eps_hat_working,2,fit$response_scale^2,'*')))<1e-10,'Noise scale binding')
  noise<-list();variables<-list()
  for(s in seq_len(fit$S))variables[[s]]<-data.frame(study=s,variable=seq_len(fit$p),
    truth_variance=tv[s,],estimated_variance=ev[s,],shared_active=seq_len(fit$p)%in%a,
    specific_active=seq_len(fit$p)%in%b[[s]])
  for(mask in c('all','shared_active','specific_active','factor_union_active'))for(scale in c('variance','plugin_sd')) {
    parts<-lapply(seq_len(fit$S),function(s) {
      cols<-switch(mask,all=seq_len(fit$p),shared_active=a,specific_active=b[[s]],factor_union_active=union(a,b[[s]]))
      x<-ev[s,cols];y<-tv[s,cols];if(scale=='plugin_sd'){x<-sqrt(x);y<-sqrt(y)}
      c(truth_mean=mean(y),estimate_mean=mean(x),bias=mean(x-y),mse=mean((x-y)^2),mae=mean(abs(x-y)),
        mean_relative_bias=mean((x-y)/y),mean_absolute_relative_error=mean(abs(x-y)/y),relative_squared_error=mean(((x-y)/y)^2))
    })
    met<-Reduce(`+`,parts)/fit$S
    noise[[length(noise)+1L]]<-data.frame(mask=mask,scale=scale,t(met),rmse=sqrt(met['mse']),
      relative_rmse=sqrt(met['relative_squared_error']),check.names=FALSE)
  }
  list(active_intervals=qsum,active_interval_subjects=qi,centered_signal=center,centered_subjects=cr,
    noise=do.call(rbind,noise),noise_variables=do.call(rbind,variables))
}

pair_winner_science <- function(fit,obs,fixed,truth,test_truth,pred,env) {
  historical<-pair_historical(fit,obs,fixed,truth,env)
  result<-adapt_winner_science(fit,obs,fixed,truth,test_truth,pred,env,historical)
  result$supplement<-pair_supplement(fit,obs,fixed,truth)
  result
}

pair_assert_table <- function(new,old,keys,fields,tolerance=1e-10) {
  key<-function(x)do.call(paste,c(x[keys],sep='::'))
  a<-key(new);b<-key(old);assert(!anyDuplicated(a)&&!anyDuplicated(b)&&setequal(a,b),'Metric key mismatch')
  old<-old[match(a,b),,drop=FALSE]
  for(nm in fields) {
    x<-new[[nm]];y<-old[[nm]]
    if(is.numeric(x)&&is.numeric(y))assert(isTRUE(all.equal(x,y,tolerance=tolerance,check.attributes=FALSE)),paste('Metric mismatch',nm))
    else assert(identical(as.character(x),as.character(y)),paste('Metric mismatch',nm))
  }
  data.frame(metric=fields,pass=TRUE,tolerance=tolerance)
}
