source(file.path(Sys.getenv('FSYNC_EXEC_ROOT'),'scripts','common.R'))
stopifnot(Sys.info()[['sysname']]=='Linux');verify_gate()
id<-commandArgs(TRUE)[1L]
plan<-read.csv(file.path(output_root,'PLAN.csv'),stringsAsFactors=FALSE)
row<-plan[plan$fit_id==id,,drop=FALSE];stopifnot(nrow(row)==1L,row$arm=='original',row$stage==1L,file.exists(file.path(output_root,'DATA_READY.txt')))
config<-jsonlite::read_json(file.path(exec_root,'config/EXECUTION_CONFIG.json'),simplifyVector=FALSE)
stopifnot(identical(config$status,'FROZEN'))
out<-new_output(file.path(output_root,'fits',id))
writeLines(runtime_lines(),file.path(out,'RUNTIME.txt'))
write.csv(data.frame(fit_id=id,pid=Sys.getpid(),status='running',started_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'STARTED.csv'),row.names=FALSE)
warnings<-character();begin<-proc.time();fit<-NULL
status<-data.frame(fit_id=id,status='failed',seconds=NA_real_,cpu_seconds=NA_real_,ELBO=NA_real_,
  t1_sweeps=NA_integer_,annealing_sweeps=NA_integer_,converged=FALSE,convergence_status='',
  decreases=NA_integer_,objective_valid=FALSE,finite=FALSE,monotonic_pass=FALSE,metadata_pass=FALSE,
  initial_common_state_pass=FALSE,initial_rng_pass=FALSE,fitting_seconds=NA_real_,fit_sha256='',fit_bytes=NA_real_,error='')
err<-tryCatch(withCallingHandlers({
  folder<-file.path(output_root,'fit_inputs',row$dataset_id);verify_sha256_file(file.path(folder,'INPUT_SHA256SUMS.txt'),folder)
  obs<-readRDS(file.path(folder,'observation_bundle.rds'));fixed<-readRDS(file.path(folder,'model_metadata.rds'))
  stopifnot(obs$data_id==row$dataset_id,is.null(obs$Z),isFALSE(obs$truth_available_to_fit))
  stopifnot(isFALSE(obs$test_response_available_to_fit),isFALSE(fixed$truth_read),
    !any(c('L_f','L_s','M_f','M_s')%in%names(obs$dimensions)),
    !any(c('L_f','L_s','M_f','M_s')%in%names(fixed$fit_metadata)))
  source(file.path(exec_root,'scripts/validation/initialization_gate.R'))
  source(file.path(exec_root,'scripts/validation/original_initial_gate.R'))
  verify_baseline_package(config$baseline_library)
  baseline<-adapt_load_baseline(config$baseline_library)
  en<-load_source('candidate')
  args<-fit_arguments(obs,row$seed,en,row$arm)
  stopifnot(is.null(args$factor_shrinkage),obs$dimensions$S==2L,obs$dimensions$p==500L,all(obs$dimensions$n_s==60L))
  options(multiFSYNC.__stage7g_snapshot_20260928__=NULL,multiFSYNC.__stage7i_resume_20260929__=NULL,
    multiFSYNC.__private_scale_trace_8d72c46f__=list(include_initial=FALSE,include_annealing=TRUE,
      t1_sweeps=integer(),block_annealing_sweeps=integer()))
  gate<-original_initial_state_gate(args,baseline,en)
  write.csv(gate$rows,file.path(out,'INITIAL_STATE_GATE.csv'),row.names=FALSE)
  new_json(gate,file.path(out,'INITIAL_STATE_GATE.json'))
  status$initial_common_state_pass<-all(gate$rows$common_q_identical&gate$rows$fixed_metadata_identical)
  status$initial_rng_pass<-all(gate$rows$rng_identical)
  stopifnot(gate$pass,status$initial_common_state_pass,status$initial_rng_pass)
  patch_path <- file.path(exec_root,'runtime_patch','diagnostic_preallocation_patch.R')
  stopifnot(identical(sha256_file(patch_path),'e319d15c81ed7395a018eb55804d36b19a17c3cee46f7fab28a91c8f3c7dd51a'))
  patch_env <- new.env(parent=globalenv()); sys.source(patch_path,envir=patch_env,keep.source=FALSE)
  patch_receipt <- patch_env$install_diagnostic_preallocation_patch(en)
  patch_receipt$patch_sha256 <- sha256_file(patch_path)
  atomic_rds(patch_receipt,file.path(out,'DIAGNOSTIC_RUNTIME_PATCH.rds'))
  atomic_rds(args[setdiff(names(args),c('Y','time_obs'))],file.path(out,'resolved_fit_config.rds'))
  options(multiFSYNC.__stage7g_snapshot_20260928__=NULL,multiFSYNC.__stage7i_resume_20260929__=NULL,
    multiFSYNC.__private_scale_trace_8d72c46f__=list(include_initial=FALSE,include_annealing=TRUE,
      t1_sweeps=integer(),block_annealing_sweeps=integer()))
  new_json(list(fit_id=id,dataset_id=row$dataset_id,seed=row$seed,baseline_library=config$baseline_library,
    source_manifest_sha256=sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt')),
    model='original_CAVI',factor_shrinkage='none',true_values_read=FALSE,
    initial_state_gate_zero_CAVI_sweeps=TRUE,diagnostic_runtime_patch=patch_receipt),file.path(out,'SOURCE_BINDING.json'))
  fit_begin<-proc.time()[['elapsed']]
  fit<-do.call(en$bayesSYNC_multi_pre_score,args)
  status$fitting_seconds<-proc.time()[['elapsed']]-fit_begin
  storage<-save_full_fit(fit,file.path(out,'fit.rds'))
  new_json(storage,file.path(out,'STORAGE_ROUNDTRIP.json'))
  status$fit_sha256<-storage$sha256;status$fit_bytes<-storage$bytes
  invariant<-c('S','n_s','p','d','K','K_total','C_g','time_g')
  checks<-data.frame(field=invariant,pass=vapply(invariant,function(nm) isTRUE(all.equal(fit[[nm]],fixed$fit_metadata[[nm]],tolerance=1e-12,check.attributes=FALSE)),logical(1)))
  write.csv(checks,file.path(out,'METADATA_CHECKS.csv'),row.names=FALSE)
  status$metadata_pass<-all(checks$pass)
  stopifnot(status$metadata_pass,fit$L_f==4L,all(fit$L_s_by_study==4L),all(fit$M_f==4L),all(unlist(fit$M_s)==4L))
  stopifnot(is.null(fit$factor_shrinkage_state),is.null(fit$factor_shrinkage))
  status$ELBO<-tail(fit$ELBO,1L);status$t1_sweeps<-fit$t1_sweeps;status$annealing_sweeps<-fit$annealing_sweeps
  status$converged<-isTRUE(fit$converged);status$convergence_status<-fit$convergence_status
  status$decreases<-fit$ELBO_decrease_count;status$objective_valid<-isTRUE(fit$elbo_objective_valid)
  state_fields<-c('mu_q_nu_mu','mu_q_nu_phi','mu_q_nu_psi','mu_q_zeta','mu_q_xi','mu_q_a','mu_q_b_specific','mu_q_gamma_a','mu_q_gamma_b','mu_q_recip_sigsq_eps','ELBO')
  status$finite<-all(is.finite(unlist(fit[state_fields])))
  status$monotonic_pass<-all(diff(fit$ELBO)>=-1e-7*pmax(1,abs(head(fit$ELBO,-1L))))
  write.csv(data.frame(t1_sweep=seq_along(fit$ELBO),ELBO=fit$ELBO),file.path(out,'ORDINARY_ELBO.csv'),row.names=FALSE)
  write.csv(fit$practical_diagnostics,file.path(out,'PRACTICAL_DIAGNOSTICS.csv'),row.names=FALSE)
  write.csv(fit$scale_trace,file.path(out,'ANNEALING_OBJECTIVE.csv'),row.names=FALSE)
  compact_storage<-save_full_fit(list(initialization=fit$initialization_diagnostics,driver_trace=fit$driver_trace,
    practical_final=fit$practical_final,checkpoints=fit$practical_checkpoints),file.path(out,'COMPACT_DIAGNOSTICS.rds'))
  new_json(compact_storage,file.path(out,'COMPACT_STORAGE_ROUNDTRIP.json'))
  stopifnot(status$t1_sweeps==400L,status$annealing_sweeps==99L,length(fit$ELBO)==400L,
            status$finite,status$objective_valid,status$monotonic_pass)
  status$status<-'complete';NULL
},warning=function(w) {warnings<<-c(warnings,conditionMessage(w));cat(conditionMessage(w),'\n',file=file.path(out,'WARNINGS.txt'),append=TRUE);invokeRestart('muffleWarning')}),
  error=function(e) conditionMessage(e))
if(!is.null(err)) {status$error<-err;writeLines(err,file.path(out,'ERROR.txt'))}
elapsed<-proc.time()-begin;status$seconds<-unname(elapsed[['elapsed']]);status$cpu_seconds<-sum(elapsed[1:2])
write.csv(status,file.path(out,'STATUS.csv'),row.names=FALSE)
cat(status$status,id,status$seconds,'seconds\n')
if(status$status!='complete') quit(status=2L)
