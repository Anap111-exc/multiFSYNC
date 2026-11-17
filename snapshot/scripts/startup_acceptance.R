# Linux startup gate: parse, source provenance, input legality and zero-sweep state equality.
# No CAVI updates, no truth deserialization and no posterior endpoint are produced.
stopifnot(Sys.info()[['sysname']]=='Linux')
source(file.path(Sys.getenv('FSYNC_EXEC_ROOT'),'scripts/common.R'))
stopifnot(file.exists(file.path(output_root,'DATA_READY.txt')))
config<-jsonlite::read_json(file.path(exec_root,'config/EXECUTION_CONFIG.json'),simplifyVector=FALSE)
stopifnot(identical(config$status,'FROZEN'),config$new_fit_count==120L,
          identical(as.integer(unlist(config$fitted_factor_capacities)),c(4L,4L,4L)),
          config$fitted_fpca_capacity==4L,config$annealing_sweeps==99L,config$ordinary_t1_sweeps==400L)
verify_sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt'),exec_root)
verify_sha256_file(file.path(exec_root,'scripts/validation/CANDIDATE_SOURCE_SHA256SUMS.txt'),exec_root)
paths<-c(list.files(file.path(exec_root,'candidate/R'),'[.]R$',full.names=TRUE),
         list.files(file.path(exec_root,'scripts'),'[.]R$',recursive=TRUE,full.names=TRUE),
         list.files(file.path(exec_root,'runtime_patch'),'[.]R$',full.names=TRUE))
invisible(lapply(paths,function(p)parse(p,keep.source=FALSE)))
plan<-read_csv(file.path(output_root,'PLAN.csv'))
stopifnot(nrow(plan)==120L,all(plan$arm=='original'),all(plan$stage==1L),
          length(unique(plan$dataset_id))==10L,!anyDuplicated(plan$fit_id),
          all(table(plan$dataset_id)==12L))
for(id in unique(plan$dataset_id)) {
  rows<-plan[plan$dataset_id==id,,drop=FALSE]
  stopifnot(identical(sort(as.integer(rows$start)),seq_len(12L)))
  folder<-file.path(output_root,'fit_inputs',id)
  verify_sha256_file(file.path(folder,'INPUT_SHA256SUMS.txt'),folder)
  obs<-readRDS(file.path(folder,'observation_bundle.rds'))
  fixed<-readRDS(file.path(folder,'model_metadata.rds'))
  stopifnot(identical(obs$data_id,id),is.null(obs$Z),isFALSE(obs$truth_available_to_fit),
    isFALSE(obs$test_response_available_to_fit),isFALSE(fixed$truth_read),
    obs$dimensions$S==2L,obs$dimensions$p==500L,all(obs$dimensions$n_s==60L),obs$dimensions$d==0L,
    !any(c('L_f','L_s','M_f','M_s')%in%names(obs$dimensions)),
    !any(c('L_f','L_s','M_f','M_s')%in%names(fixed$fit_metadata)),
    all(is.finite(unlist(obs$Y))),length(obs$Y)==2L,
    all(vapply(obs$Y,length,integer(1))==60L))
}
verify_baseline_package(config$baseline_library)
source(file.path(exec_root,'scripts/validation/initialization_gate.R'))
source(file.path(exec_root,'scripts/validation/original_initial_gate.R'))
baseline<-adapt_load_baseline(config$baseline_library);en<-load_source('candidate')
first<-plan[1L,,drop=FALSE]
obs<-readRDS(file.path(output_root,'fit_inputs',first$dataset_id,'observation_bundle.rds'))
args<-fit_arguments(obs,first$seed,en,'original')
options(multiFSYNC.__stage7g_snapshot_20260928__=NULL,multiFSYNC.__stage7i_resume_20260929__=NULL,
  multiFSYNC.__private_scale_trace_8d72c46f__=list(include_initial=FALSE,include_annealing=TRUE,
    t1_sweeps=integer(),block_annealing_sweeps=integer()))
begin<-proc.time()[['elapsed']]
gate<-original_initial_state_gate(args,baseline,en)
stopifnot(gate$pass)
patch_env<-new.env(parent=globalenv())
patch<-file.path(exec_root,'runtime_patch/diagnostic_preallocation_patch.R')
stopifnot(identical(sha256_file(patch),'e319d15c81ed7395a018eb55804d36b19a17c3cee46f7fab28a91c8f3c7dd51a'))
sys.source(patch,envir=patch_env)
receipt<-patch_env$install_diagnostic_preallocation_patch(en)
dir.create(file.path(output_root,'startup'),recursive=TRUE,showWarnings=FALSE)
new_json(list(status='PASS',parsed_R_files=length(paths),fit_count=120L,independent_datasets=10L,
  capacity=c(4L,4L,4L),FPCA_capacity=4L,initial_state_gate=gate,patch=receipt,
  truth_values_read=0L,CAVI_updates=0L,elapsed_seconds=proc.time()[['elapsed']]-begin,
  source_manifest_sha256=sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt'))),
  file.path(output_root,'startup/RUNNER_ACCEPTANCE.json'))
writeLines(runtime_lines(),file.path(output_root,'startup/RUNNER_RUNTIME.txt'))
manifest_out<-file.path(output_root,'startup/SOURCE_SHA256SUMS.txt')
if(!file.exists(manifest_out)) {
  new_lines(readLines(file.path(exec_root,'CODE_SHA256SUMS.txt'),warn=FALSE),manifest_out)
} else {
  stopifnot(identical(readLines(manifest_out,warn=FALSE),readLines(file.path(exec_root,'CODE_SHA256SUMS.txt'),warn=FALSE)))
}
gate_out<-file.path(output_root,'startup/GATE_PASS.txt')
if(!file.exists(gate_out))new_lines('PASS: generation, source, projected input and original initial-state gates passed',gate_out)
cat('RUNNER_ACCEPTANCE PASS; zero CAVI updates; original arm only\n')
