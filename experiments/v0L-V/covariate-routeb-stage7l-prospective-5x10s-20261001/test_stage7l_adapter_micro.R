#!/usr/bin/env Rscript
# Small adapter check using the existing p=3 toy; no prospective data.
args<-commandArgs(TRUE);stopifnot(length(args)==4L)
exec_dir<-args[1L];lib<-args[2L];toy<-args[3L];out<-args[4L]
source(file.path(exec_dir,"stage7l_common.R"));load_stage7l(exec_dir,lib)
assert(!dir.exists(out),"Micro output exists.");dir.create(out,recursive=TRUE)
obs<-readRDS(file.path(toy,"TOY_OBSERVATION.rds"))
obs$dimensions<-list(K=3L,p=as.integer(length(obs$Y[[1L]][[1L]])))
old<-readRDS(file.path(toy,"TOY_ORIGINAL_SNAPSHOTS.rds"))
find_state<-function(i) {
 hit<-Filter(function(x)x$iteration==i,old)
 assert(length(hit)==1L,"Toy snapshot missing.");hit[[1L]]
}
resume<-new.env(parent=emptyenv());resume$schema<-"stage7i_sweep50_resume_v1"
resume$state<-find_state(50L);resume$stop_after_iteration<-52L
sink<-make_snapshot_sink(c(51L,52L),out)
options(multiFSYNC.__stage7i_resume_20260929__=resume,
        multiFSYNC.__stage7g_snapshot_20260928__=sink,
        multiFSYNC.__private_scale_trace_8d72c46f__=NULL)
fa<-fit_arguments(obs,99102L)
fa$M_f<-1L;fa$M_s<-list(1L,1L);fa$n_g<-21L
fit<-do.call(multiFSYNC::bayesSYNC_multi_pre_score,fa)
assert(isTRUE(resume$applied) && fit$i_iter==52L && length(sink$records)==2L,
       "Micro resume integration failed.")
parity<-lapply(sink$records,function(x) {
 ref<-find_state(x$iteration)
 data.frame(iteration=x$iteration,field=names(x),pass=vapply(names(x),function(nm)
  isTRUE(all.equal(x[[nm]],ref[[nm]],tolerance=1e-10,check.attributes=TRUE)),logical(1)))
})
parity<-do.call(rbind,parity)
new_csv(parity,file.path(out,"TWO_SWEEP_TOY_PARITY.csv"))
assert(all(parity$pass),"Micro resumed raw states disagree with existing toy.")
assert(identical(readRDS(file.path(out,"STATE_SWEEP051.rds")),sink$records[[1L]]) &&
       identical(readRDS(file.path(out,"STATE_SWEEP052.rds")),sink$records[[2L]]),
       "Incremental snapshot persistence failed.")
fixed<-fixed_metadata(obs)
assert(all(vapply(fixed$C[[2L]],ncol,integer(1))==5L) &&
       isTRUE(all.equal(fixed$fit_metadata$list_cp_C,fit$list_cp_C,check.attributes=FALSE)),
       "Independent metadata designs disagree with fitted designs.")
new_lines(c("STAGE7L_ADAPTER_MICRO_PASS","new_prospective_data=FALSE",
 "historical_449_sweep_refit=FALSE","resumed_updates=2"),file.path(out,"MICRO_PASS.txt"))
cat("STAGE7L_ADAPTER_MICRO_PASS\n")
