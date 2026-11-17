#!/usr/bin/env Rscript
# New registered Linux datasets only. No fitting, endpoint selection or model update.
options(warn=1,stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
           BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
stopifnot(identical(Sys.info()[['sysname']],'Linux'))
argv<-commandArgs(TRUE)
stopifnot(length(argv)%in%c(2L,3L))
ex<-normalizePath(argv[1L],winslash='/',mustWork=TRUE)
out<-normalizePath(argv[2L],winslash='/',mustWork=FALSE)
fit_lib<-if(length(argv)==3L)argv[3L] else '/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1'
gen<-file.path(ex,'generation')
stopifnot(startsWith(out,'/root/multifsync-results/'),
          basename(out)=='cavi_true222_fit444_m4_20261117_v1')
dep_libs<-c('/root/nocov-nc01-nc05-20261006-v1/generation/library',
            '/root/mgp-unseen-validation-dependencies-20261020-v1')
.libPaths(c(fit_lib,dep_libs[dir.exists(dep_libs)],.libPaths()))
stopifnot(requireNamespace('digest',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE),
          requireNamespace('Splinets',quietly=TRUE))
source(file.path(gen,'nc_common.R'),local=.GlobalEnv)
source(file.path(gen,'stage7m_reconstruct.R'),local=.GlobalEnv)
source(file.path(gen,'nc_fixed_metadata.R'),local=.GlobalEnv)
source(file.path(gen,'nc_generation_core.R'),local=.GlobalEnv)
source(file.path(gen,'project_fit_inputs.R'),local=.GlobalEnv)
new_json<-function(x,path){
  assert(!file.exists(path)&&!file.exists(paste0(path,'.part')),paste('Refusing overwrite',path))
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-paste0(path,'.part')
  jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
  assert(file.rename(tmp,path),'JSON atomic rename failed')
}
verify_source_manifest<-function(path,base){
  lines<-readLines(path,warn=FALSE)
  assert(length(lines)>0L&&all(grepl('^[0-9a-f]{64}  .+',lines)),'Malformed source manifest')
  for(line in lines){
    rel<-substring(line,67L)
    assert(!grepl('(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))',rel),'Unsafe source path')
    check_hash(file.path(base,rel),substr(line,1L,64L))
  }
}
contract<-jsonlite::read_json(file.path(gen,'GENERATION_CONTRACT.json'),simplifyVector=FALSE)
assert(contract$experiment_id==basename(out)&&contract$independent_datasets==10L&&
       contract$registered_fits==120L&&isFALSE(contract$automatic_retry),'Wrong generation contract')
frozen<-jsonlite::read_json(file.path(ex,'config/PROTOCOL_FROZEN.json'),simplifyVector=FALSE)
assert(identical(frozen$status,'FROZEN')&&is.character(frozen$protocol_sha256),'Protocol is not frozen')
check_hash(file.path(ex,'config/PROTOCOL.json'),frozen$protocol_sha256)
verify_source_manifest(file.path(ex,'CODE_SHA256SUMS.txt'),ex)
check_hash(file.path(ex,'config/DATASET_REGISTRATION.csv'),contract$registration_sha256)
check_hash(file.path(gen,'FIT_SEED_REFERENCE.csv'),contract$fit_seed_reference_sha256)
check_hash(file.path(gen,'SEED_REGISTRATION.csv'),contract$seed_registration_sha256)
check_hash(file.path(gen,'stage7m_reconstruct.R'),contract$reconstruct_sha256)
library(multiFSYNC,lib.loc=fit_lib)
loaded<-normalizePath(getNamespaceInfo('multiFSYNC','path'),winslash='/',mustWork=TRUE)
assert(identical(loaded,normalizePath(file.path(fit_lib,'multiFSYNC'),winslash='/',mustWork=TRUE)),
       'Unexpected generator package library')
generator_sha<-nc_canonical_function_sha256(multiFSYNC::simulate_multi_study_structured)
assert(identical(generator_sha,contract$generator_canonical_sha256),'Frozen NC generator binding differs')
rows<-read_csv(file.path(ex,'config/DATASET_REGISTRATION.csv'))
plan<-read_csv(file.path(ex,'config/PLAN.csv'))
starts<-read_csv(file.path(gen,'FIT_SEED_REFERENCE.csv'))
expected_ids<-sprintf('NC01_TRUE222_20261117__rep%02d',1:10)
assert(nrow(rows)==10L&&identical(as.character(rows$data_id),expected_ids),'Dataset registration changed')
assert(all(rows$S==2L)&all(rows$n_each==60L)&all(rows$p==500L)&all(rows$L_f==2L)&
       all(rows$M==2L)&all(rows$K==9L)&all(rows$d_generated==0L)&all(rows$d_fit==0L),'Registered dimensions changed')
assert(nrow(plan)==120L&&all(plan$arm=='original')&&!anyDuplicated(plan$fit_id)&&
       setequal(plan$dataset_id,expected_ids),'Run plan differs from 120 original-CAVI fits')
key<-paste(plan$dataset_id,plan$start);ix<-match(key,paste(starts$dataset_id,starts$start))
assert(!anyNA(ix)&&!anyDuplicated(key)&&identical(as.integer(plan$seed),as.integer(starts$seed[ix])),
       'Run plan start seeds differ from unchanged registered references')
assert(all(table(plan$dataset_id)==12L)&&all(vapply(split(plan$start,plan$dataset_id),
       function(x)setequal(x,1:12),logical(1))),'Expected twelve starts per dataset')
assert(!dir.exists(file.path(out,'data'))&&!dir.exists(file.path(out,'fit_inputs'))&&
       !file.exists(file.path(out,'GENERATION_STARTED.json'))&&!file.exists(file.path(out,'DATA_READY.txt')),
       'Generation destination exists; no automatic replacement or retry')
dir.create(out,recursive=TRUE,showWarnings=FALSE,mode='0700')
new_json(list(status='RUNNING',experiment_id=contract$experiment_id,new_datasets=0L,
              registered_fits=120L,automatic_retries=0L),file.path(out,'GENERATION_STARTED.json'))
began<-proc.time()[['elapsed']];done<-character();checks<-loadings<-bindings<-list();projections<-list()
error<-NULL
tryCatch({
  for(index in seq_len(nrow(rows))){
    row<-rows[index,,drop=FALSE];id<-as.character(row$data_id)
    pair<-nc_generate_dataset(row,contract$experiment_id);obs<-pair$observation;tr<-pair$data$true_params
    qc<-nc_generation_qc(pair,row)
    assert(all(qc$qc$pass),paste('Original NC generation QC failed for',id))
    assert(tr$L_f==2L&&identical(as.integer(tr$L_s_by_study),c(2L,2L))&&
           identical(as.integer(tr$n_s),c(60L,60L))&&tr$p==500L&&tr$truth_basis_dim==8L,
           'True dimensions differ from registration')
    disjoint<-all(vapply(seq_len(2L),function(s){
      support<-cbind(tr$a_true,tr$b_true[[s]])!=0
      all(colSums(support)==50L)&&all(rowSums(support)<=1L)
    },logical(1)))
    assert(disjoint,'Loading supports overlap or have unexpected active count')
    fixed<-nc_fixed_metadata(obs)
    assert(isFALSE(fixed$truth_read)&&all(fixed$fit_metadata$n_s==60L)&&fixed$fit_metadata$K_total==11L,
           'Invalid train-only fixed metadata')
    assert(length(unlist(pair$data$C))==length(unlist(fixed$C)),'Generator and fit C dimensions differ')
    basis_error<-max(abs(unlist(pair$data$C)-unlist(fixed$C)))
    assert(basis_error<1e-12,'Generated C differs from frozen fit metadata C')
    set.seed(as.integer(row$test_time_seed))
    test_time<-lapply(obs$time_obs,function(st)lapply(st,function(tt)sort(runif(3L,min(tt),max(tt)))))
    set.seed(as.integer(row$test_noise_seed))
    noise<-lapply(test_time,function(st)lapply(st,function(tt)matrix(rnorm(length(tt)*tr$p,sd=.3),length(tt),tr$p)))
    tb<-stage7m_truth_basis(tr)
    signal<-lapply(seq_len(2L),function(s)lapply(seq_len(60L),function(i)
      stage7m_true_parts(tr,obs,s,i,test_time[[s]][[i]],tb)$signal))
    test_truth<-list(time=test_time,noise=noise,signal=signal,Y=Map(function(a,b)Map(`+`,a,b),signal,noise))
    test_blind<-list(time=test_time,Z=NULL,subject_ids=lapply(tr$n_s,seq_len),
                     truth_read=FALSE,test_response_read=FALSE)
    assert(all_finite(test_truth)&&all(lengths(unlist(test_time,recursive=FALSE))==3L),
           'Invalid held-out same-subject times or values')
    sealed<-list(bundle_class='nc_sealed_truth_v1',experiment_id=contract$experiment_id,
      scenario_id=as.character(row$scenario_id),data_id=id,truth_unseal_authorized=FALSE,
      seed_record=as.list(row),data=list(true_params=tr))
    base<-file.path(out,'data',id)
    paths<-c(observation=file.path(base,'observation_bundle.rds'),metadata=file.path(base,'model_metadata.rds'),
      test_design=file.path(base,'test_design_blind.rds'),truth=file.path(base,'sealed_truth/truth_bundle.rds'),
      test_truth=file.path(base,'sealed_truth/test_truth.rds'))
    objects<-list(observation=obs,metadata=fixed,test_design=test_blind,truth=sealed,test_truth=test_truth)
    for(role in names(paths))new_rds(objects[[role]],paths[[role]])
    bind<-data.frame(data_id=id,scenario_id=as.character(row$scenario_id),role=names(paths),
      relative_path=substring(unname(paths),nchar(out)+2L),sha256=vapply(paths,sha,character(1)),bytes=file.info(paths)$size)
    bindings[[index]]<-bind
    projected<-project_unseen_fit_inputs(obs,fixed)
    fitdir<-file.path(out,'fit_inputs',id)
    new_rds(projected$observation,file.path(fitdir,'observation_bundle.rds'))
    new_rds(projected$metadata,file.path(fitdir,'model_metadata.rds'))
    fitnames<-c('observation_bundle.rds','model_metadata.rds')
    new_lines(vapply(fitnames,function(nm)paste(sha(file.path(fitdir,nm)),nm,sep='  '),character(1)),
              file.path(fitdir,'INPUT_SHA256SUMS.txt'))
    projections[[id]]<-c(list(data_id=id,source_observation_sha256=sha(paths[['observation']]),
      source_metadata_sha256=sha(paths[['metadata']]),projected_observation_sha256=sha(file.path(fitdir,fitnames[1L])),
      projected_metadata_sha256=sha(file.path(fitdir,fitnames[2L]))),projected$audit)
    allowed<-bind[bind$role%in%c('observation','metadata','test_design'),,drop=FALSE]
    new_lines(paste(allowed$sha256,basename(allowed$relative_path),sep='  '),file.path(base,'INPUT_SHA256SUMS.txt'))
    new_lines(c('NEW_REGISTERED_NC_DATA_SEALED','generated_d=0','fit_d=0','truth_and_targets_excluded_from_fit=TRUE'),
              file.path(base,'DATA_SEAL_COMPLETE.txt'))
    extra<-data.frame(data_id=id,check_id=c('registered_true_222','fixed_loading_support_disjoint',
      'linux_generated_C_matches_fixed_metadata','heldout_times_three_each','fit_projection_no_truth_dimensions'),
      pass=c(TRUE,disjoint,basis_error<1e-12,TRUE,TRUE))
    checks[[index]]<-rbind(qc$qc,extra);loadings[[index]]<-qc$loading_qc
    done<-c(done,id);cat('REGISTERED_DATA_COMPLETE',id,'elapsed_seconds=',proc.time()[['elapsed']]-began,'\n')
    rm(pair,obs,tr,qc,fixed,tb,signal,noise,test_truth,test_blind,sealed,objects,projected);gc(FALSE)
  }
  new_csv(do.call(rbind,checks),file.path(out,'DGP_QC.csv'))
  new_csv(do.call(rbind,loadings),file.path(out,'DGP_LOADING_QC.csv'))
  new_csv(do.call(rbind,bindings),file.path(out,'DATA_BINDINGS.csv'))
  assert(file.copy(file.path(ex,'config/PLAN.csv'),file.path(out,'PLAN.csv'),overwrite=FALSE),'Plan copy failed')
  new_json(list(status='PASS',datasets=10L,fit_input_projection=projections,
    true_dimensions_excluded_from_fit=TRUE,basis_and_observation_values_unchanged=TRUE),
    file.path(out,'FIT_INPUT_PROJECTION_AUDIT.json'))
  runtime<-list(actual_loaded_library=loaded,package_version=as.character(packageVersion('multiFSYNC')),
    generator_canonical_sha256=generator_sha,Splinets_version=as.character(packageVersion('Splinets')),
    R_version=R.version.string,RNGkind=as.list(RNGkind()),sys_info=as.list(Sys.info()),
    ext_software=as.list(extSoftVersion()),session_info=capture.output(sessionInfo()),libraries=.libPaths(),
    thread_environment=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))))
  new_json(runtime,file.path(out,'GENERATION_RUNTIME.json'))
  dir.create(file.path(out,'startup'),recursive=TRUE,showWarnings=FALSE)
  new_lines(readLines(file.path(ex,'CODE_SHA256SUMS.txt'),warn=FALSE),file.path(out,'startup/SOURCE_SHA256SUMS.txt'))
  new_lines('PASS: generation QC only; separate RUNNER_ACCEPTANCE required; ten NEW registered datasets',
            file.path(out,'startup/GATE_PASS.txt'))
  verify_source_manifest(file.path(ex,'CODE_SHA256SUMS.txt'),ex)
  new_json(list(status='PASS',experiment_id=contract$experiment_id,new_datasets=10L,reused_datasets=0L,
    registered_fits=120L,new_initializations=0L,completed_data_ids=done,truth_read_by_fit=FALSE,
    generating_truth_read_for_legality_qc_only=TRUE,test_response_read_by_fit=FALSE,
    actual_generator_library=loaded,canonical_formals_plus_body_sha256=list(generator=generator_sha),
    registration_sha256=contract$registration_sha256,generation_contract_sha256=sha(file.path(gen,'GENERATION_CONTRACT.json')),
    data_bindings_sha256=sha(file.path(out,'DATA_BINDINGS.csv')),elapsed_seconds=proc.time()[['elapsed']]-began,
    automatic_retries=0L,replacement_seeds=0L),file.path(out,'GENERATION_STATUS.json'))
  new_lines('TEN_NEW_REGISTERED_INDEPENDENT_DATASETS_READY',file.path(out,'DATA_READY.txt'))
},error=function(e)error<<-conditionMessage(e))
if(!is.null(error)){
  if(!file.exists(file.path(out,'GENERATION_STATUS.json')))new_json(list(status='FAIL',error=error,
    new_datasets=length(done),reused_datasets=0L,registered_fits=120L,completed_data_ids=done,
    automatic_retries=0L,replacement_seeds=0L),file.path(out,'GENERATION_STATUS.json'))
  new_lines(error,file.path(out,'GENERATION_FAILURE.txt'));stop(error,call.=FALSE)
}
cat('PASS_GENERATION',length(done),'new independent datasets, 120 registered original-CAVI fits, zero fits run\n')
