#!/usr/bin/env Rscript
argv <- commandArgs(TRUE)
stopifnot(length(argv) >= 4L)
action <- argv[1L]; exec_dir <- argv[2L]; output <- argv[3L]; lib <- argv[4L]
source(file.path(exec_dir,"stage7l_common.R"))
load_stage7l(exec_dir,lib)
manifest <- read_csv(file.path(exec_dir,"FIT_MANIFEST.csv"))
data_rows <- read_csv(file.path(exec_dir,"DATA_MANIFEST.csv"))
assert(nrow(manifest)==50L && nrow(data_rows)==5L &&
         !anyDuplicated(manifest$fit_seed) && all(data_rows$beta_arm=="C1"),
       "Prospective manifest invalid.")
if (action=="generate") {
  assert(!dir.exists(output), "Refusing output-root reuse.")
  dir.create(output,recursive=TRUE,mode="0700")
  file.copy(file.path(exec_dir,c("DATA_MANIFEST.csv","FIT_MANIFEST.csv",
    "PROTOCOL_STAGE7L_PROSPECTIVE_20261001_v1.md","REGISTRATION_FROZEN.txt",
    "SOURCE_FILES.sha256")),output)
  qc <- list()
  for (i in seq_len(nrow(data_rows))) {
    row <- data_rows[i,,drop=FALSE]
    pair <- paired_beta_data(row,FALSE)
    checks <- pair_checks(pair,i)
    assert(all(checks$passed),paste("DGP check failed:",row$data_id))
    save_arm(pair$C1,row,output,FALSE)
    qc[[i]] <- checks
    rm(pair); gc()
  }
  new_csv(do.call(rbind,qc),file.path(output,"DGP_QC.csv"))
  paths <- unlist(lapply(data_rows$data_id,function(id) c(
    file.path(output,"data",id,"observation_bundle.rds"),
    file.path(output,"data",id,"sealed_truth","truth_bundle.rds"))))
  new_csv(data.frame(path=paths,sha256=vapply(paths,sha,character(1))),
          file.path(output,"DATA_INPUT_BINDING.csv"))
  new_lines(c("STAGE7L_DATA_GENERATED_AND_SEALED","datasets=5","arms=C1",
    "fit_reads_truth=FALSE"),file.path(output,"DATA_GENERATION_COMPLETE.txt"))
}
if (action=="fit") {
  assert(length(argv)==6L,"Expected branch and fit_id.")
  branch <- argv[5L]; fit_id <- argv[6L]
  assert(branch %in% c("U","S"),"Unsupported branch.")
  source_U <- ""
  if (branch=="U") row <- manifest[manifest$fit_id==fit_id,,drop=FALSE] else {
    winners <- read_csv(file.path(output,"U_WINNERS_FROZEN.csv"))
    row <- winners[winners$S_fit_id==fit_id,,drop=FALSE]
    assert(nrow(row)==1L && nzchar(row$fit_id),"S has no fixed U winner.")
    source_U <- row$fit_id
    row$fit_id <- fit_id
  }
  assert(nrow(row)==1L,"fit_id is not unique.")
  fit_dir <- file.path(output,"fits",fit_id)
  assert(!dir.exists(fit_dir),"Refusing task rerun.")
  dir.create(fit_dir,recursive=TRUE,mode="0700")
  term <- terminal_row(row,branch,source_U)
  warnings <- character(); start <- proc.time()[["elapsed"]]
  sink <- NULL; resume <- NULL
  outcome <- tryCatch(withCallingHandlers({
    obs_path <- file.path(output,"data",row$data_id,"observation_bundle.rds")
    inputs <- read_csv(file.path(output,"DATA_INPUT_BINDING.csv"))
    hit <- inputs[inputs$path==obs_path,,drop=FALSE]
    assert(nrow(hit)==1L,"Observation input binding missing.")
    check_hash(obs_path,hit$sha256)
    obs <- readRDS(obs_path)
    assert(identical(obs$experiment_id,experiment_id) && isFALSE(obs$truth_available_to_fit),
           "Input is not the bound observation-only bundle.")
    sweeps <- if (branch=="U") 50L else c(51L,55L,60L,70L,100L,499L)
    sink <- make_snapshot_sink(sweeps,fit_dir)
    resume <- NULL
    if (branch=="S") {
      state_path <- file.path(output,"S_construction",row$data_id,"STATE_S_SWEEP50.rds")
      resume <- new.env(parent=emptyenv())
      resume$schema <- "stage7i_sweep50_resume_v1"
      resume$state <- readRDS(state_path); resume$stop_after_iteration <- 499L
    }
    options(multiFSYNC.__stage7g_snapshot_20260928__=sink,
            multiFSYNC.__stage7i_resume_20260929__=resume,
            multiFSYNC.__private_scale_trace_8d72c46f__=NULL)
    do.call(multiFSYNC::bayesSYNC_multi_pre_score,fit_arguments(obs,row$fit_seed))
  },warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  }),error=function(e)e)
  term$elapsed_seconds <- proc.time()[["elapsed"]]-start
  new_lines(warnings,file.path(fit_dir,"WARNINGS.txt"))
  if (inherits(outcome,"error")) {
    term$error <- conditionMessage(outcome)
  } else {
    fit <- outcome
    new_rds(fit,file.path(fit_dir,"fit.rds"))
    term$fit_path <- file.path(fit_dir,"fit.rds")
    new_csv(data.frame(t1_sweep=seq_along(fit$ELBO),ordinary_elbo=fit$ELBO),
            file.path(fit_dir,"ORDINARY_ELBO.csv"))
    new_csv(fit$practical_diagnostics,file.path(fit_dir,"PRACTICAL_DIAGNOSTICS.csv"))
    term$global_iteration <- fit$i_iter
    term$ordinary_t1_sweeps <- fit$t1_sweeps
    term$actual_updates <- fit$i_iter-if (branch=="S") 50L else 0L
    term$final_ordinary_elbo <- if(length(fit$ELBO)) tail(fit$ELBO,1L) else NA_real_
    term$strict_practical_converged <- isTRUE(fit$practical_converged) &&
      isTRUE(fit$converged) && identical(fit$convergence_status,"converged_practical")
    objective <- objective_status(fit)
    complete <- fit$i_iter==499L && fit$t1_sweeps==400L &&
      fit$annealing_sweeps==99L && length(fit$ELBO)==400L &&
      all(is.finite(fit$ELBO)) && fit$n_cpus_used==1L &&
      length(sink$records)==length(sink$sweeps) &&
      identical(vapply(sink$records,`[[`,integer(1),"iteration"),sink$sweeps) &&
      (branch=="U" || isTRUE(resume$applied))
    term$status <- if(isTRUE(complete)) "COMPLETE" else "CONTRACT_FAILED"
    term$error <- if(isTRUE(complete)) "" else "Not a finite complete fixed-budget endpoint."
    term$objective_eligible <- isTRUE(complete) && objective$eligible
    term$invalid_reasons <- paste(objective$invalid_reasons,collapse=";")
    if (branch=="U") {
      term$snapshot50_path <- file.path(fit_dir,"STATE_SWEEP050.rds")
      fixed <- fixed_metadata(obs)
      compare <- c("S","L_f","L_s","M_s","K_total","K","p","n_s","d",
                   "C_g","time_g","list_cp_C")
      assert(all(vapply(compare,function(nm) isTRUE(all.equal(
        fixed$fit_metadata[[nm]],fit[[nm]],check.attributes=FALSE)),logical(1))) &&
        isTRUE(all.equal(fixed$fit_metadata$list_hyper$Sigma_beta,
                        fit$list_hyper$Sigma_beta,check.attributes=FALSE)),
        "Independent fixed metadata disagrees with fit static fields.")
      new_rds(fixed,file.path(fit_dir,"FIXED_MODEL_METADATA.rds"))
    }
  }
  if (branch=="U" && file.exists(file.path(fit_dir,"STATE_SWEEP050.rds")))
    term$snapshot50_path <- file.path(fit_dir,"STATE_SWEEP050.rds")
  new_csv(term[,terminal_columns],file.path(fit_dir,"terminal.csv"))
  new_lines(c(paste0("status=",term$status),"truth_read=FALSE",
              paste0("actual_updates=",term$actual_updates)),file.path(fit_dir,"TASK_FINISHED.txt"))
  if(term$status!="COMPLETE") quit(status=2L)
}
if (action=="select_u") {
  terminals <- do.call(rbind,lapply(manifest$fit_id,function(id)
    read_csv(file.path(output,"fits",id,"terminal.csv"))))
  assert(nrow(terminals)==50L,"U task accounting incomplete.")
  winners <- do.call(rbind,lapply(seq_len(nrow(data_rows)),function(i) {
    d <- data_rows$data_id[i]
    pool <- terminals[terminals$data_id==d & terminals$status=="COMPLETE" &
      terminals$objective_eligible & is.finite(terminals$final_ordinary_elbo),,drop=FALSE]
    if(nrow(pool)) row <- pool[order(-pool$final_ordinary_elbo,pool$fit_id)[1L],,drop=FALSE] else {
      row <- terminal_row(data.frame(fit_id="",data_id=d,start_index=0L,fit_seed=NA_integer_),"U")
      row$status <- "NO_ELIGIBLE_U"
    }
    row$S_fit_id <- paste0(d,"__S")
    row
  }))
  new_csv(terminals,file.path(output,"U_TERMINALS.csv"))
  new_csv(winners,file.path(output,"U_WINNERS_FROZEN.csv"))
  new_lines(c("STAGE7L_U_WINNERS_FROZEN", "selection_truth_read=FALSE",
    paste0("winner_sha256=",sha(file.path(output,"U_WINNERS_FROZEN.csv")))),
    file.path(output,"U_SELECTION_FROZEN.txt"))
}
if (action=="construct_s") {
  id <- argv[5L]
  row <- read_csv(file.path(output,"U_WINNERS_FROZEN.csv"))
  row <- row[row$data_id==id,,drop=FALSE]
  assert(nrow(row)==1L && nzchar(row$fit_id),"No registered U winner.")
  target <- file.path(output,"S_construction",id)
  assert(!dir.exists(target),"Refusing construction retry.")
  dir.create(target,recursive=TRUE,mode="0700")
  warnings <- character()
  outcome <- tryCatch(withCallingHandlers({
    paths <- c(observation=file.path(output,"data",id,"observation_bundle.rds"),
      snapshot50=row$snapshot50_path,
      metadata=file.path(output,"fits",row$fit_id,"FIXED_MODEL_METADATA.rds"),
      helper=file.path(exec_dir,"dependencies","uas_helper.R"),
      winner_selection=file.path(output,"U_WINNERS_FROZEN.csv"))
    inputs <- read_csv(file.path(output,"DATA_INPUT_BINDING.csv"))
    check_hash(paths[["observation"]],inputs$sha256[inputs$path==paths[["observation"]]])
    new_csv(data.frame(role=names(paths),path=unname(paths),sha256=vapply(paths,sha,character(1))),
            file.path(target,"S_INPUT_BINDING.csv"))
    obs <- readRDS(paths[["observation"]]); snap <- readRDS(paths[["snapshot50"]])
    fixed <- readRDS(paths[["metadata"]])
    uas <- stage7i_build_uas_states(obs,fixed$fit_metadata,snap,fixed$C[[2L]])
    assert(all_finite(uas[c("U","A","S")]),"Nonfinite constructed states.")
    new_rds(uas$S,file.path(target,"STATE_S_SWEEP50.rds"))
    new_rds(uas$A,file.path(target,"STATE_A_UNUSED_SWEEP50.rds"))
    new_rds(uas$blind_diagnostics,file.path(target,"S_BLIND_DIAGNOSTICS.rds"))
    flat <- unlist(uas$blind_diagnostics,recursive=TRUE)
    new_csv(data.frame(field=names(flat),value=as.character(flat)),
            file.path(target,"S_AMPLITUDE_AND_LEGALITY.csv"))
    TRUE
  },warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e)e)
  new_lines(warnings,file.path(target,"WARNINGS.txt"))
  if(inherits(outcome,"error")) {
    new_lines(conditionMessage(outcome),file.path(target,"CONSTRUCTION_FAILED.txt"))
    quit(status=2L)
  }
  new_lines(c("STAGE7L_S_CONSTRUCTION_COMPLETE","truth_read=FALSE","A_fitted=FALSE"),
            file.path(target,"CONSTRUCTION_COMPLETE.txt"))
}
