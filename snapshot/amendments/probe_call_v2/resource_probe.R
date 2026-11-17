# Two ordinary updates for capacity/resource acceptance only. Not a formal candidate.
# Invocation: set FSYNC_EXEC_ROOT/FSYNC_OUTPUT_ROOT; Rscript --vanilla resource_probe.R
stopifnot(Sys.info()[['sysname']]=='Linux')
source(file.path(Sys.getenv('FSYNC_EXEC_ROOT'),'scripts/common.R'))
verify_gate()
acceptance<-jsonlite::read_json(file.path(output_root,'startup/RUNNER_ACCEPTANCE.json'),simplifyVector=FALSE)
stopifnot(identical(acceptance$status,'PASS'))
probe_seed<-200611171L
plan<-read_csv(file.path(output_root,'PLAN.csv'))
stopifnot(nrow(plan)==120L,!probe_seed%in%plan$seed,all(plan$arm=='original'))
id<-plan$dataset_id[1L]
folder<-file.path(output_root,'fit_inputs',id)
verify_sha256_file(file.path(folder,'INPUT_SHA256SUMS.txt'),folder)
obs<-readRDS(file.path(folder,'observation_bundle.rds'))
stopifnot(isFALSE(obs$truth_available_to_fit),isFALSE(obs$test_response_available_to_fit),
  is.null(obs$Z),all(obs$dimensions$n_s==60L),obs$dimensions$p==500L)
dest<-new_output(file.path(output_root,'resource_probe/original_L4M4_2sweep_v2'))
en<-load_source('candidate')
source(file.path(exec_root,'scripts/validation/initialization_gate.R'))
source(file.path(exec_root,'runtime_patch/diagnostic_preallocation_patch.R'))
patch<-install_diagnostic_preallocation_patch(en)
args<-fit_arguments(obs,probe_seed,en,'original')
args['anneal']<-list(NULL);args$maxit<-2L
stopifnot('anneal'%in%names(args),is.null(args$anneal))
args$practical_control$window<-1L;args$practical_control$long_window<-1L
args$practical_control$min_t1<-2L;args$practical_control$consecutive<-1L
args$practical_control$max_t1<-2L;args$practical_control$checkpoints<-2L
args$trace_sweeps<-c(1L,2L)
en$.validate_practical_control(args$practical_control)
options(multiFSYNC.__stage7g_snapshot_20260928__=NULL,multiFSYNC.__stage7i_resume_20260929__=NULL,
  multiFSYNC.__private_scale_trace_8d72c46f__=list(include_initial=FALSE,include_annealing=TRUE,
    t1_sweeps=integer(),block_annealing_sweeps=integer()))
new_json(list(status='FROZEN_BEFORE_EXECUTION',purpose='resource_acceptance_not_candidate',data_id=id,
  seed=probe_seed,formal_fit_budget_used=0L,extra_scientific_starts=0L,annealing_sweeps=0L,
  ordinary_sweeps=2L,capacity=c(4L,4L,4L),fpca_capacity=4L,
  formal_configuration_unchanged='99 annealing and400 ordinary; full practical20/60 windows',
  limitations='Short state size is a lower-bound measurement. Early iteration timing is not a full-fit ETA.',
  input_manifest_sha256=sha256_file(file.path(folder,'INPUT_SHA256SUMS.txt')),
  source_manifest_sha256=sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt'))),file.path(dest,'PROBE_FROZEN.json'))
writeLines(runtime_lines(),file.path(dest,'RUNTIME.txt'))
atomic_rds(args[setdiff(names(args),c('Y','time_obs'))],file.path(dest,'RESOLVED_CONFIG.rds'))
error<-NULL;report<-NULL
tryCatch({
  initial_begin<-proc.time()[['elapsed']]
  initial<-adapt_capture_initial(en,args)
  initial_seconds<-proc.time()[['elapsed']]-initial_begin
  initial_sizes<-list(snapshot_object_bytes=as.numeric(object.size(initial)),
    variational_state_object_bytes=as.numeric(object.size(initial$common_q)),
    fixed_metadata_object_bytes=as.numeric(object.size(initial$fixed)),
    elapsed_seconds=initial_seconds,CAVI_updates=0L)
  initial_fields<-rbind(
    data.frame(group='common_q',field=names(initial$common_q),object_bytes=vapply(initial$common_q,function(x)as.numeric(object.size(x)),numeric(1))),
    data.frame(group='fixed',field=names(initial$fixed),object_bytes=vapply(initial$fixed,function(x)as.numeric(object.size(x)),numeric(1))))
  write.csv(initial_fields,file.path(dest,'INITIAL_FIELD_OBJECT_SIZES.csv'),row.names=FALSE)
  rm(initial);invisible(gc(FALSE))
  # The actual fit resets its explicit independent seed, so the zero-update capture
  # does not alter the probe realization or any formal start's RNG state.
  begin<-proc.time()
  fit<-do.call(en$bayesSYNC_multi_pre_score,args)
  elapsed<-proc.time()-begin
  stopifnot(fit$L_f==4L,all(fit$L_s_by_study==4L),all(fit$M_f==4L),all(unlist(fit$M_s)==4L),
    fit$t1_sweeps==2L,fit$annealing_sweeps==0L,length(fit$ELBO)==2L,
    all(is.finite(fit$ELBO)),isTRUE(fit$elbo_objective_valid),is.null(fit$factor_shrinkage_state))
  storage<-save_full_fit(fit,file.path(dest,'PROBE_FIT.rds'))
  new_json(storage,file.path(dest,'PROBE_STORAGE.json'))
  compact<-list(initialization=fit$initialization_diagnostics,driver_trace=fit$driver_trace,
    practical_final=fit$practical_final,checkpoints=fit$practical_checkpoints)
  compact_storage<-save_full_fit(compact,file.path(dest,'COMPACT_DIAGNOSTICS.rds'))
  field_sizes<-data.frame(field=names(fit),object_bytes=vapply(fit,function(x)as.numeric(object.size(x)),numeric(1)))
  write.csv(field_sizes,file.path(dest,'FIELD_OBJECT_SIZES.csv'),row.names=FALSE)
  practical_fields<-c('stability_snapshot','practical_final','practical_checkpoints','practical_diagnostics','driver_trace')
  practical_sizes<-lapply(intersect(practical_fields,names(fit)),function(nm) {
    x<-fit[[nm]]
    list(field=nm,object_bytes=as.numeric(object.size(x)),
      component_object_bytes=if(is.list(x)) lapply(x,function(z)as.numeric(object.size(z))) else NULL)
  })
  new_json(practical_sizes,file.path(dest,'PRACTICAL_LAST_STATE_SIZES.json'))
  raw_response_count<-sum(vapply(obs$Y,function(study)sum(vapply(study,function(subject)sum(lengths(subject)),numeric(1))),numeric(1)))
  proc<-readLines('/proc/self/status',warn=FALSE)
  highwater_kib<-as.numeric(sub('^VmHWM:\\s+([0-9]+).*','\\1',grep('^VmHWM:',proc,value=TRUE)))
  report<-list(status='PASS',data_id=id,seed=probe_seed,formal_fit_budget_used=0L,
    fit_seconds=unname(elapsed[['elapsed']]),cpu_seconds=sum(elapsed[1:2]),
    initial_state=initial_sizes,fit_object_bytes=as.numeric(object.size(fit)),fit_xz_bytes=storage$bytes,
    compact_xz_bytes=compact_storage$bytes,peak_rss_mib=highwater_kib/1024,
    response_count=raw_response_count,full_history_fitted_vectors_uncompressed_estimate_bytes=8*raw_response_count*61,
    full_history_is_internal_not_serialized=TRUE,practical_last_state_fields=practical_sizes,
    probe_file_bytes_total=sum(file.info(list.files(dest,full.names=TRUE))$size,na.rm=TRUE),
    full_state_storage_not_upper_bounded=TRUE,ordinary_sweeps=2L,annealing_sweeps=0L,
    true_values_read=FALSE,endpoint_is_candidate=FALSE,patch=patch,
    limitation='No 499-step cost or full endpoint compression bound is established by two updates.')
  write.csv(fit$driver_trace,file.path(dest,'DRIVER_TRACE.csv'),row.names=FALSE)
},error=function(e){error<<-conditionMessage(e)})
if(!is.null(error))report<-list(status='FAIL',error=error,formal_fit_budget_used=0L,automatic_retry=FALSE)
new_json(report,file.path(dest,'RESOURCE_PROBE_COMPLETE.json'))
cat(jsonlite::toJSON(report,auto_unbox=TRUE),'\n')
if(!is.null(error))quit(status=2L)
