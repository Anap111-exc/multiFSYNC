# Zero-CAVI argument transport and initial-state interception test.
root<-normalizePath(commandArgs(TRUE)[1L],winslash='/',mustWork=TRUE)
Sys.setenv(FSYNC_EXEC_ROOT=root)
source(file.path(root,'scripts/common.R'))
source(file.path(root,'scripts/validation/initialization_gate.R'))
en<-load_source('candidate')
obs<-tiny_observations();args<-fit_arguments(obs,200611171L,en,'original',smoke=TRUE)
args['anneal']<-list(NULL);args$maxit<-2L
args$convergence_rule<-'practical';args$practical_control<-en$g12_stopping_control()
args$practical_control$window<-1L;args$practical_control$long_window<-1L
args$practical_control$min_t1<-2L;args$practical_control$consecutive<-1L
args$practical_control$max_t1<-2L;args$practical_control$checkpoints<-2L
args$trace_sweeps<-c(1L,2L)
stopifnot('anneal'%in%names(args),is.null(args$anneal))
# Execute the actual public forwarding wrapper against a capture-only driver.
stub<-adapt_clone_functions(en)
stub$.bayesSYNC_multi_driver_diagnostic<-function(...,anneal=c(1,1.9,100),maxit=500L,control=NULL) {
  stopifnot(!missing(anneal),is.null(anneal),maxit==2L)
  list(explicit_null_forwarded=TRUE,CAVI_updates=0L)
}
transport<-do.call(stub$bayesSYNC_multi_pre_score,args)
stopifnot(isTRUE(transport$explicit_null_forwarded))
# Actual candidate validation and initialization must reach the interception
# immediately before the first for-loop. No CAVI update is executed.
options(multiFSYNC.__stage7g_snapshot_20260928__=NULL,multiFSYNC.__stage7i_resume_20260929__=NULL)
captured<-adapt_capture_initial(en,args)
stopifnot(is.list(captured$common_q),is.null(captured$factor_shrinkage_state))
out<-list(status='PASS',explicit_null_named_argument=TRUE,public_wrapper_forwarding=TRUE,
  actual_core_validation_and_initialization_passed=TRUE,CAVI_updates=0L,
  formal_fit_attempts=0L,new_scientific_starts=0L,source_modified=FALSE,
  test_data='existing deterministic tiny observation fixture, not production data',
  probe_seed_unchanged=200611171L)
jsonlite::write_json(out,file.path(root,'amendments/probe_call_v2/EXPLICIT_NULL_TEST_PASS.json'),
  auto_unbox=TRUE,pretty=TRUE)
cat(jsonlite::toJSON(out,auto_unbox=TRUE),'\n')
