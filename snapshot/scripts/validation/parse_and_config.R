# Local read-only numerical-code inspection. Does not call a fitting function.
root<-normalizePath(commandArgs(TRUE)[1L],winslash='/',mustWork=TRUE)
Sys.setenv(FSYNC_EXEC_ROOT=root)
source(file.path(root,'scripts/common.R'))
verify_sha256_file(file.path(root,'scripts/validation/CANDIDATE_SOURCE_SHA256SUMS.txt'),root)
paths<-c(list.files(file.path(root,'candidate/R'),'[.]R$',full.names=TRUE),
  list.files(file.path(root,'scripts'),'[.]R$',recursive=TRUE,full.names=TRUE),
  list.files(file.path(root,'runtime_patch'),'[.]R$',full.names=TRUE))
for(p in paths)parse(p,keep.source=FALSE)
en<-load_source('candidate')
obs<-list(Y=list(list(vector('list',500L))),time_obs=NULL,dimensions=list(K=9L))
args<-fit_arguments(obs,20261117L,en,'original')
stopifnot(identical(args$L_f,4L),identical(args$L_s,c(4L,4L)),identical(args$M_f,rep(4L,4L)),
  identical(args$M_s,list(rep(4L,4L),rep(4L,4L))),is.null(args$factor_shrinkage),
  args$maxit==499L,identical(args$anneal,c(1,1.9,100)),args$practical_control$max_t1==400L,
  args$practical_control$min_t1==396L,args$practical_control$consecutive==5L,
  args$pre_score_sweeps==1L,args$function_initialization=='gram_unit_energy',args$initialization=='random')
source(file.path(root,'runtime_patch/diagnostic_preallocation_patch.R'))
patch<-install_diagnostic_preallocation_patch(en)
source(file.path(root,'scripts/validation/original_initial_gate.R'))
receipt<-list(status='PASS',parsed_files=length(paths),numerical_source_files=28L,
  capacity=c(4L,4L,4L),fpca_capacity=4L,annealing_sweeps=99L,ordinary_sweeps=400L,
  active_MGP=FALSE,random_numbers_drawn=FALSE,fit_calls=0L,CAVI_updates=0L,
  diagnostic_preallocation_patch=patch,R_version=R.version.string)
jsonlite::write_json(receipt,file.path(root,'scripts/validation/LOCAL_PARSE_ACCEPTANCE.json'),pretty=TRUE,auto_unbox=TRUE)
cat(jsonlite::toJSON(receipt,auto_unbox=TRUE),'\n')
