# Original-CAVI-only initialization validation. Zero CAVI sweeps; no truth.
# The imported capture helper copies functions without editing package namespaces.
original_initial_state_gate <- function(args,baseline,candidate,retain_state=FALSE) {
  stopifnot(is.null(args$factor_shrinkage),args$L_f==4L,all(args$L_s==4L),
            all(args$M_f==4L),all(unlist(args$M_s)==4L))
  caller_kind<-RNGkind()
  caller_had_seed<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  caller_seed<-if(caller_had_seed)get('.Random.seed',.GlobalEnv)else NULL
  on.exit({
    do.call(RNGkind,as.list(caller_kind))
    if(caller_had_seed)assign('.Random.seed',caller_seed,.GlobalEnv)
    else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  },add=TRUE)
  b<-adapt_capture_initial(baseline,args)
  c<-adapt_capture_initial(candidate,args)
  states<-list(baseline=b,candidate_original=c)
  rows<-do.call(rbind,lapply(names(states),function(nm){
    s<-states[[nm]]
    data.frame(arm=nm,common_q_identical=identical(b$common_q,s$common_q),
      fixed_metadata_identical=identical(b$fixed,s$fixed),
      rng_identical=identical(b$rng_seed,s$rng_seed)&&identical(b$rng_kind,s$rng_kind),
      common_q_sha256=adapt_sha(s$common_q),fixed_sha256=adapt_sha(s$fixed),rng_sha256=adapt_sha(s$rng_seed),
      extra_mgp_state=!is.null(s$factor_shrinkage_state),stringsAsFactors=FALSE)
  }))
  list(pass=all(rows$common_q_identical&rows$fixed_metadata_identical&rows$rng_identical)&!any(rows$extra_mgp_state),
    rows=rows,stage='after_random_and_G12_before_first_CAVI_loop',original_CAVI_sweeps=0L,
    pre_score_executed=FALSE,caller_rng_restore_guard=TRUE,mgp_state_initialized=FALSE,
    anchors=lapply(states,function(s)attr(s,'capture_anchor')),state=if(retain_state)states else NULL)
}
