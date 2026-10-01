#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
exec_dir <- args[1L]; lib <- args[2L]; support <- args[3L]
source(file.path(exec_dir,"stage7l_common.R"));load_stage7l(exec_dir,lib)
for(p in list.files(exec_dir,pattern="[.]R$",full.names=TRUE)) parse(p)
check_hash(file.path(exec_dir,"dependencies","stage7c_runner.R"),
 "56b13bb6ff3673c95449afd748976472ba9c3136feb435f78279a95e2d16c519")
check_hash(file.path(exec_dir,"dependencies","stage7c_evaluator.R"),
 "1ecb4901eac60e53f4e7969735847542a1d78c4b734cfb5611deadbacbbc43b8")
packages<-c("digest","jsonlite","Splinets")
present<-vapply(packages,requireNamespace,logical(1),quietly=TRUE)
assert(all(present),paste("Missing dependencies:",paste(packages[!present],collapse=",")))
splinets<-find.package("Splinets")
assert(packageVersion("Splinets")>=package_version("1.5.1"),"Unexpected Splinets version.")
paths<-file.path(splinets,c("DESCRIPTION","NAMESPACE","R/Splinets.rdb","R/Splinets.rdx"))
new_csv(data.frame(role=basename(paths),path=paths,sha256=vapply(paths,sha,character(1))),
 file.path(exec_dir,"DGP_DEPENDENCY_BINDING.csv"))
assert(identical(as.integer(formals(multiFSYNC::bayesSYNC_multi_pre_score)$pre_score_sweeps),1L),
       "Unexpected pre-score default interface.")
core <- paste(deparse(body(multiFSYNC:::bayesSYNC_multi_core)),collapse="\n")
assert(grepl("stage7i_sweep50_resume_v1",core,fixed=TRUE) &&
       grepl("stage7g_snapshot_v1",core,fixed=TRUE),"Required frozen hooks unavailable.")
source(file.path(exec_dir,"stage7l_evaluation_setup.R"))
evaluation_setup(exec_dir,lib,support,load=FALSE)
new_csv(as.data.frame(installed.packages()[,c("Package","Version","LibPath","Built")]),
 file.path(exec_dir,"DEPENDENCY_ENVIRONMENT.csv"))
new_lines(c(capture.output(sessionInfo()),capture.output(Sys.info()),
 paste0("RNGkind=",paste(RNGkind(),collapse=",")),
 paste(names(Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","LC_ALL"))),
       Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","LC_ALL")),sep="=")),
 file.path(exec_dir,"R_ENVIRONMENT.txt"))
cat("STAGE7L_R_PREFLIGHT_PASS fit_started=FALSE truth_opened=FALSE\n")
