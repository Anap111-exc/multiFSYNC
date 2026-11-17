# Pure raw endpoint acceptance and report-only postprocessing helpers.
iv_check_endpoint <- function(fit,metadata,seed,elbo,core) {
  iv_assert(fit$S==2L&&fit$p==500L&&all(fit$n_s==30L)&&fit$d==0L&&
    fit$L_f==2L&&all(fit$L_s_by_study==2L)&&all(fit$M_f==3L)&&all(unlist(fit$M_s)==3L)&&
    fit$K==9L&&fit$K_total==11L&&fit$t1_sweeps==400L&&fit$annealing_sweeps==99L&&
    fit$i_iter==499L&&length(fit$ELBO)==400L&&isTRUE(fit$elbo_objective_valid)&&
    isTRUE(all.equal(tail(fit$ELBO,1L),elbo,tolerance=1e-12,check.attributes=FALSE))&&
    identical(as.integer(fit$fit_seed),as.integer(seed)),'Frozen endpoint design/objective/counters mismatch')
  iv_assert(all(is.finite(fit$ELBO))&&all(diff(fit$ELBO)>=-1e-7*pmax(1,abs(head(fit$ELBO,-1L)))),
    'Frozen original monotonicity gate failed')
  ms<-fit$factor_shrinkage_state
  iv_assert(!is.null(ms)&&ms$version==1L&&ms$control$mode=='mgp'&&ms$control$a1==2.1&&
    ms$control$a2==3.1&&ms$control$rate==1&&ms$p==500L&&ms$L_f==2L&&all(ms$L_s==2L),'MGP state differs')
  q<-core$ns_make_global_q(fit,metadata);theta<-core$ns_theta(q,FALSE)
  iv_assert(!any(grepl('zeta|xi|score',names(q)))&&identical(theta$global_q_draw,FALSE),
    'Scoring adapter accesses fitted scores')
  invisible(TRUE)
}

iv_report_endpoint <- function(raw,en) {
  before<-iv_hash(raw);fit<-raw
  required<-setdiff(names(formals(en$orthonormalise_multi)),c('zeta_true','pve_cutoff'))
  iv_assert(all(required%in%names(fit)),'Endpoint lacks report-only postprocessor fields')
  args<-fit[required];args['zeta_true']<-list(NULL);args$pve_cutoff<-.99
  hats<-do.call(en$orthonormalise_multi,args);retained<-setdiff(names(raw),names(hats));fit[names(hats)]<-hats
  iv_assert(identical(fit[retained],raw[retained])&&identical(before,iv_hash(raw)),
    'Reporting changed raw endpoint')
  fit
}
