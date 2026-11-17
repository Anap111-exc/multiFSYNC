options(warn=1)
compiler::enableJIT(0)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
           BLIS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
exec_root <- normalizePath(Sys.getenv('FSYNC_EXEC_ROOT', unset='E:/multiFSYNC-results/cavi_true222_fit444_m4_20261117_v1'),winslash='/',mustWork=TRUE)
output_root <- Sys.getenv('FSYNC_OUTPUT_ROOT',unset=file.path(exec_root,'results','local_validation'))
load_source <- function(version='candidate') {
  stopifnot(version %in% c('baseline','candidate'))
  en <- new.env(parent=globalenv()); en$spline.des <- splines::spline.des
  for(path in sort(list.files(file.path(exec_root,version,'R'),'[.]R$',full.names=TRUE))) sys.source(path,envir=en,keep.source=FALSE)
  for(nm in ls(en,all.names=TRUE)) if(is.function(en[[nm]])) en[[nm]] <- compiler::cmpfun(en[[nm]])
  en
}
new_output <- function(path) {
  if(file.exists(path)) stop('Refusing overwrite: ',path)
  stopifnot(dir.create(path,recursive=TRUE));path
}
atomic_rds <- function(x,path) {
  stopifnot(!file.exists(path),!file.exists(paste0(path,'.partial')))
  saveRDS(x,paste0(path,'.partial'),version=3);stopifnot(file.rename(paste0(path,'.partial'),path))
}
sha256_file <- function(path) digest::digest(file=path,algo='sha256',serialize=FALSE)
verify_sha256_file <- function(manifest,directory) {
  lines <- readLines(manifest,warn=FALSE)
  stopifnot(length(lines)>0,all(grepl('^[0-9a-f]{64}  .+',lines)))
  for(line in lines) {
    rel <- substring(line,67)
    stopifnot(!grepl('(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))',rel))
    stopifnot(identical(sha256_file(file.path(directory,rel)),substr(line,1,64)))
  }
  invisible(TRUE)
}
verify_gate <- function() {
  stopifnot(file.exists(file.path(output_root,'startup','GATE_PASS.txt')))
  acceptance <- jsonlite::read_json(file.path(output_root,'startup','RUNNER_ACCEPTANCE.json'),simplifyVector=FALSE)
  stopifnot(identical(acceptance$status,'PASS'),acceptance$CAVI_updates==0L)
  a <- readLines(file.path(output_root,'startup','SOURCE_SHA256SUMS.txt'),warn=FALSE)
  stopifnot(identical(a,readLines(file.path(exec_root,'CODE_SHA256SUMS.txt'),warn=FALSE)))
  verify_sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt'),exec_root)
  stopifnot(identical(acceptance$source_manifest_sha256,sha256_file(file.path(exec_root,'CODE_SHA256SUMS.txt'))))
}
runtime_lines <- function() c(capture.output(sessionInfo()),capture.output(Sys.info()),capture.output(extSoftVersion()),
  capture.output(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'))))
assert <- function(value,message) if(!isTRUE(value)) stop(message,call.=FALSE)
check_hash <- function(path,expected) assert(identical(sha256_file(path),expected),paste('Hash mismatch:',path))
read_csv <- function(path) read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
new_lines <- function(value,path) {assert(!file.exists(path),paste('Refusing overwrite:',path));writeLines(value,path)}
new_json <- function(value,path) {
  assert(!file.exists(path),paste('Refusing overwrite:',path))
  jsonlite::write_json(value,path,pretty=TRUE,auto_unbox=TRUE,digits=17,na='null')
}
verify_baseline_package <- function(lib) {
  host <- if(Sys.info()[['sysname']]=='Linux') 'Linux_4.5.3' else 'Windows_4.5.1'
  hashes <- read_csv(file.path(exec_root,'scripts/validation/HOST_INSTALLED_PACKAGE_HASHES.csv'))
  hashes <- hashes[hashes$host_id==host,,drop=FALSE]
  assert(nrow(hashes)>4L,'Baseline host binding missing')
  for(i in seq_len(nrow(hashes))) check_hash(file.path(lib,'multiFSYNC',hashes$package_relative_path[i]),hashes$sha256[i])
  invisible(TRUE)
}
save_full_fit <- function(fit,path) {
  assert(!file.exists(path)&&!file.exists(paste0(path,'.partial')),'Endpoint storage path exists')
  before <- digest::digest(fit,algo='sha256',serialize=TRUE)
  saveRDS(fit,paste0(path,'.partial'),compress='xz',version=3L)
  recovered <- readRDS(paste0(path,'.partial'))
  assert(identical(fit,recovered),'Complete endpoint roundtrip differs')
  assert(identical(before,digest::digest(recovered,algo='sha256',serialize=TRUE)),'Recovered object hash differs')
  assert(file.rename(paste0(path,'.partial'),path),'Endpoint rename failed')
  list(path=path,sha256=sha256_file(path),bytes=unname(file.info(path)$size),object_sha256=before,
       compression='xz',serialization_version=3L,identical_roundtrip=TRUE)
}
tiny_observations <- function() {
  tt<-lapply(1:2,function(s) lapply(1:4,function(i) seq(0,1,length.out=9+i%%2)))
  Y<-lapply(1:2,function(s) lapply(1:4,function(i) lapply(1:6,function(j) {
    t<-tt[[s]][[i]]; .1*s+.03*j+sin((1+j%%3)*pi*t+.4*i)+.05*cos((s+i)*pi*t)
  })))
  list(Y=Y,Z=NULL,time_obs=tt)
}
fit_arguments <- function(obs,seed,en,arm,smoke=FALSE) {
  p<-length(obs$Y[[1]][[1]])
  control<-en$g12_stopping_control();control$checkpoints<-c(380L,400L)
  args<-list(Y=obs$Y,Z=NULL,time_obs=obs$time_obs,
    L_f=4L,L_s=c(4L,4L),M_f=rep(4L,4L),M_s=list(rep(4L,4L),rep(4L,4L)),
    K=if(smoke) 3L else as.integer(obs$dimensions$K),n_g=51L,anneal=c(1,1.9,100),maxit=499L,
    n_cpus=1L,verbose=FALSE,seed=as.integer(seed),tol_abs=1e-3,tol_rel=1e-5,
    bool_scale=FALSE,bool_var_spec_prob=FALSE,d_0=as.integer(p),convergence_rule='practical',
    lambda_orth=0,practical_control=control,initialization='random',
    initialization_control=list(grid_size=81L,perturb_sd=.05,rank_tol=1e-8),continuation_state=NULL,
    pre_score_sweeps=1L,trace_sweeps=c(1L,380L,400L),function_initialization='gram_unit_energy')
  stopifnot(identical(arm,'original'),is.null(args$factor_shrinkage))
  if(smoke) {
    args$anneal<-c(1,1.9,3);args$maxit<-10L;args$convergence_rule<-'parameters'
    args$practical_control<-NULL;args$tol_abs<-0;args$tol_rel<-0;args$trace_sweeps<-1L
  }
  args
}
