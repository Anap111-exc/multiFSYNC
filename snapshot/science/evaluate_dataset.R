#!/usr/bin/env Rscript
# Seal the original-CAVI twelve-start choice before reading either evaluation truth.
options(warn=1)
compiler::enableJIT(0)
arguments <- commandArgs(TRUE)
if(length(arguments)==6L) {
  stopifnot(setequal(arguments[c(1,3,5)],c('--exec-root','--output-root','--data-id')))
  values <- setNames(arguments[c(2,4,6)],arguments[c(1,3,5)])
  arguments <- unname(values[c('--exec-root','--output-root','--data-id')])
}
stopifnot(length(arguments)==3L)
ex <- normalizePath(arguments[1],winslash='/',mustWork=TRUE)
out <- normalizePath(arguments[2],winslash='/',mustWork=TRUE)
id <- arguments[3]
source(file.path(ex,'science/pair_science.R'))
env <- pair_load(ex)
iv_threads()
started <- proc.time()[['elapsed']]
truth_reads <- 0L
checks <- list()
dest <- file.path(out,'science',id)
selection_root <- file.path(out,'selection',id)
iv_assert(!dir.exists(dest)&&!dir.exists(selection_root),'Refuse overwrite or automatic science retry')
dir.create(dest,recursive=TRUE)
dir.create(selection_root,recursive=TRUE)

add <- function(name,pass) {
  checks[[length(checks)+1L]] <<- data.frame(check=name,pass=isTRUE(pass))
  iv_assert(pass,name)
}
record <- function(path) {
  list(path=normalizePath(path,winslash='/',mustWork=TRUE),sha256=iv_sha(path),bytes=unname(file.info(path)$size))
}
verify_record <- function(item) {
  iv_check(item$path,item$sha256)
  iv_assert(file.info(item$path)$size==item$bytes,'Bound artifact size changed')
}
bind_rows <- function(rows) {
  if(length(rows))do.call(rbind,rows)else data.frame()
}

status <- tryCatch({
  code_hash <- iv_verify_package(ex)
  protocol_path <- file.path(ex,'config/PROTOCOL.json')
  protocol_hash <- iv_sha(protocol_path)
  plan_path <- file.path(out,'PLAN.csv')
  plan <- read.csv(plan_path,stringsAsFactors=FALSE)
  add('plan_matches_registration',identical(plan,read.csv(file.path(ex,'config/PLAN.csv'),stringsAsFactors=FALSE)))
  registered <- plan[plan$dataset_id==id,,drop=FALSE]
  arms <- 'original'
  add('registered_original_12',nrow(registered)==12L&&!anyDuplicated(registered$fit_id)&&setequal(registered$arm,arms))
  for(arm in arms) {
    rows <- registered[registered$arm==arm,,drop=FALSE]
    add(paste0(arm,'_12_starts'),nrow(rows)==12L&&setequal(rows$start,1:12))
  }
  observation_binding <- iv_data_binding(out,id,'observation')
  metadata_binding <- iv_data_binding(out,id,'metadata')
  design_binding <- iv_data_binding(out,id,'test_design')
  obs <- iv_read_binding(observation_binding)
  fixed <- iv_read_binding(metadata_binding)
  test_design <- iv_read_binding(design_binding)
  add('observation_identity',identical(obs$data_id,id)&&is.null(obs$Z)&&isFALSE(obs$truth_available_to_fit))
  artifacts <- lapply(c(plan_path,file.path(out,'DATA_BINDINGS.csv'),observation_binding$path,
                       metadata_binding$path,design_binding$path),record)
  attempts <- lapply(seq_len(nrow(registered)),function(index) {
    row <- registered[index,,drop=FALSE]
    status_path <- file.path(out,'fits',row$fit_id,'STATUS.csv')
    iv_assert(file.exists(status_path),paste('Missing terminal status',row$fit_id))
    terminal <- read.csv(status_path,stringsAsFactors=FALSE)
    iv_assert(nrow(terminal)==1L&&terminal$fit_id==row$fit_id&&
      terminal$status%in%c('complete','failed','timeout','terminated','missing'),
      'All 12 registered attempts require terminal receipts before selection')
    artifacts[[length(artifacts)+1L]] <<- record(status_path)
    fit_path <- file.path(out,'fits',row$fit_id,'fit.rds')
    valid <- fc_valid_status(terminal)
    reason <- ''
    fit_hash <- ''
    if(valid) {
      error <- tryCatch({
        iv_assert(isTRUE(terminal$initial_common_state_pass)&&isTRUE(terminal$initial_rng_pass),
          'Initial state or RNG gate failed')
        iv_check(fit_path,terminal$fit_sha256)
        fit_hash <- iv_sha(fit_path)
        raw <- readRDS(fit_path)
        pair_check_endpoint(raw,fixed,row$seed,terminal$ELBO,row$arm)
        rm(raw)
        gc(FALSE)
        NULL
      },error=function(error)conditionMessage(error))
      if(!is.null(error)) {
        valid <- FALSE
        reason <- error
      }
    } else {
      reason <- if('error'%in%names(terminal)&&!is.na(terminal$error)&&nzchar(terminal$error))
        terminal$error else 'terminal_status_not_valid'
    }
    if(nzchar(fit_hash))artifacts[[length(artifacts)+1L]] <<- record(fit_path)
    data.frame(data_id=id,arm=row$arm,fit_id=row$fit_id,start=row$start,seed=row$seed,
      status=terminal$status,valid=valid,ELBO=if(is.finite(terminal$ELBO))terminal$ELBO else NA_real_,
      fit_path=fit_path,fit_sha256=fit_hash,reason=reason,
      fit_seconds=if('seconds'%in%names(terminal))terminal$seconds else NA_real_,stringsAsFactors=FALSE)
  })
  attempts <- bind_rows(attempts)
  choices <- winner_records <- list()
  for(arm in arms) {
    arm_attempts <- attempts[attempts$arm==arm,,drop=FALSE]
    choices[[arm]] <- fc_choose(arm_attempts)
    folder <- file.path(selection_root,arm)
    dir.create(folder)
    iv_csv(arm_attempts,file.path(folder,'CANDIDATE_VALIDITY.csv'))
    seal <- list(schema='true222_fit444_m4_selection_v1',data_id=id,arm=arm,status='FROZEN_BEFORE_TRUTH',
      registered=12L,valid=sum(arm_attempts$valid),all_12_attempts_terminal=TRUE,
      incomplete_valid_pool=sum(arm_attempts$valid)<12L,truth_reads=0L,fit_updates=0L,
      rule='maximum_complete_ordinary_ELBO',code_manifest_sha256=code_hash,
      protocol_sha256=protocol_hash,artifacts=artifacts,decision=choices[[arm]],cross_prior_ELBO_comparison=FALSE)
    path <- file.path(folder,'WINNER.json')
    iv_json(seal,path)
    winner_records[[arm]] <- record(path)
  }
  pair_seal_path <- file.path(selection_root,'WINNERS.json')
  iv_json(list(status='FROZEN_BEFORE_TRUTH',data_id=id,registered=12L,truth_reads=0L,
    fit_updates=0L,winners=winner_records,single_arm_original=TRUE),pair_seal_path)
  pair_seal_hash <- iv_sha(pair_seal_path)
  valid <- attempts[attempts$valid,,drop=FALSE]

  Sys.setenv(FSYNC_EXEC_ROOT=ex,FSYNC_OUTPUT_ROOT=out)
  source(file.path(ex,'scripts/common.R'))
  inference <- load_source('candidate')
  legacy <- ae_legacy_environment(file.path(ex,'evaluation/legacy'))
  fits <- predictions <- list()
  basis_checks <- ppi_checks <- list()
  for(index in seq_len(nrow(valid))) {
    row <- valid[index,,drop=FALSE]
    raw <- adapt_read_bound(row$fit_path,row$fit_sha256)
    fit <- iv_report_endpoint(raw,inference)
    gate <- fc_basis_ppi_check(fit,obs,fixed,legacy)
    basis_checks[[length(basis_checks)+1L]] <- cbind(fit_id=row$fit_id,arm=row$arm,gate$basis)
    ppi_checks[[length(ppi_checks)+1L]] <- cbind(fit_id=row$fit_id,arm=row$arm,gate$ppi)
    if(row$fit_id==choices[[row$arm]]$selected_fit_id)fits[[row$arm]] <- fit
    rm(raw,fit)
    gc(FALSE)
  }
  for(arm in names(fits))predictions[[arm]] <- stage7m_prediction(fits[[arm]],obs,fixed$basis,test_design)
  add('all_predictions_truth_blind',all(vapply(predictions,function(value)isFALSE(value$truth_read),logical(1))))
  iv_rds(predictions,file.path(selection_root,'WINNER_PREDICTIONS_BLIND.rds'))
  iv_csv(bind_rows(basis_checks),file.path(dest,'PRE_TRUTH_BASIS_CHECKS.csv'))
  iv_csv(bind_rows(ppi_checks),file.path(dest,'PRE_TRUTH_PPI_CHECKS.csv'))
  iv_json(list(status='PASS_BEFORE_TRUTH',truth_reads=0L,pair_winners_sha256=pair_seal_hash),
    file.path(dest,'PRE_TRUTH_GATE.json'))
  for(item in c(artifacts,winner_records))verify_record(item)
  iv_check(pair_seal_path,pair_seal_hash)

  truth_binding <- iv_data_binding(out,id,'truth')
  test_truth_binding <- iv_data_binding(out,id,'test_truth')
  sealed <- iv_read_binding(truth_binding)
  truth_reads <- truth_reads+1L
  test_truth <- iv_read_binding(test_truth_binding)
  truth_reads <- truth_reads+1L
  truth <- adapt_evaluation_truth(sealed,obs,fixed)
  true_parameters <- truth$data$true_params
  add('truth_design',true_parameters$S==2L&&true_parameters$p==500L&&all(true_parameters$n_s==60L)&&
    true_parameters$L_f==2L&&all(true_parameters$L_s_by_study==2L)&&true_parameters$d_use==0L&&
    all(lengths(true_parameters$lambda_phi_true)==2L)&&
    all(vapply(true_parameters$lambda_psi_true,function(x)all(lengths(x)==2L),logical(1))))
  add('reserved_times_binding',isTRUE(all.equal(test_design$time,test_truth$time,check.attributes=FALSE)))

  tables <- list()
  append_table <- function(name,value,row) {
    if(!is.data.frame(value)||!nrow(value))return(invisible(NULL))
    prefix <- data.frame(data_id=id,fit_id=row$fit_id,start=row$start,arm=row$arm,ELBO=row$ELBO)
    prefix <- prefix[rep(1L,nrow(value)),,drop=FALSE]
    rownames(prefix) <- NULL
    tables[[name]][[length(tables[[name]])+1L]] <<- cbind(prefix,value)
  }
  for(index in seq_len(nrow(valid))) {
    row <- valid[index,,drop=FALSE]
    raw <- adapt_read_bound(row$fit_path,row$fit_sha256)
    before <- iv_hash(raw)
    fit <- iv_report_endpoint(raw,inference)
    basic <- ae_raw_metrics(fit,true_parameters,obs,fixed)
    formal <- ae_formal_structure(fit,true_parameters,obs,legacy)
    prefix <- toupper(row$arm)
    append_table(paste0(prefix,'_ENDPOINT_SCIENCE'),cbind(basic$summary,formal$summary),row)
    append_table(paste0(prefix,'_RAW_FACTOR_DIAGNOSTICS'),basic$factors,row)
    append_table(paste0(prefix,'_ROLE_RECOVERY'),basic$roles,row)
    append_table(paste0(prefix,'_FACTOR_CONTRIBUTION_BY_STUDY'),basic$energies,row)
    for(nm in c('global_truth','global_candidate','study_truth','study_candidate','shared_support','block_subspace'))
      append_table(paste0(prefix,'_FORMAL_',toupper(nm)),formal[[nm]],row)
    add(paste(row$fit_id,'state_unchanged'),identical(before,iv_hash(raw)))
    rm(raw,fit,basic,formal)
    gc(FALSE)
  }

  summary_rows <- raw_rows <- interval_rows <- historical_rows <- list()
  for(arm in names(fits)) {
    fit <- fits[[arm]]
    before <- iv_hash(fit)
    row <- valid[valid$fit_id==choices[[arm]]$selected_fit_id,,drop=FALSE]
    set.seed(as.integer(row$seed))
    result <- pair_winner_science(fit,obs,fixed,truth,test_truth,predictions[[arm]],env)
    formal <- ae_formal_structure(fit,true_parameters,obs,legacy)
    basic <- ae_raw_metrics(fit,true_parameters,obs,fixed)
    summary <- cbind(data.frame(data_id=id,fit_id=row$fit_id,start=row$start,arm=arm,ELBO=row$ELBO),
      basic$summary,formal$summary)
    append_table('WINNER_LEGACY_LOSSES',ae_winner_losses(fit,true_parameters,obs,fixed,test_truth,legacy),row)
    append_table('WINNER_FPCA_ROLE_RANKS',basic$roles,row)
    folder <- file.path(dest,arm)
    dir.create(folder)
    iv_rds(result,file.path(folder,'WINNER_SCIENCE.rds'))
    for(name in c('subject_metrics','role_metrics','component_metrics','loading_metrics','variable_metrics',
      'prediction_metrics','alignment_maps','score_points','loading_points','postprocess_metrics'))
      iv_csv(result$extended[[name]],file.path(folder,paste0(toupper(name),'.csv')))
    iv_csv(result$raw$subject,file.path(folder,'RAW_COMPONENT_SUBJECT_METRICS.csv'))
    iv_csv(result$raw_roles,file.path(folder,'RAW_ALIGNED_ROLE_METRICS.csv'))
    iv_csv(result$intervals$subject,file.path(folder,'Q_INTERVAL_SUBJECT_METRICS.csv'))
    iv_csv(result$historical$primary$block_subspace,file.path(folder,'BLOCK_SUBSPACE.csv'))
    for(name in names(result$supplement))iv_csv(result$supplement[[name]],file.path(folder,paste0(toupper(name),'.csv')))
    summary_rows[[arm]] <- cbind(summary,result$structure$summary,result$extended$summary)
    raw_rows[[arm]] <- cbind(data.frame(data_id=id,arm=arm),result$raw$summary)
    interval_rows[[arm]] <- cbind(data.frame(data_id=id,arm=arm),result$intervals$summary)
    historical_rows[[arm]] <- cbind(data.frame(data_id=id,arm=arm),result$historical$summary)
    add(paste(arm,'winner_state_unchanged'),identical(before,iv_hash(fit)))
    cat('Winner full science',id,arm,'complete\n')
    flush.console()
    rm(result,fit)
    gc(FALSE)
  }
  for(name in names(tables))iv_csv(bind_rows(tables[[name]]),file.path(dest,paste0(name,'.csv')))
  # Empty files preserve the output contract for an arm with no valid winner.
  for(name in c('ORIGINAL_ENDPOINT_SCIENCE','WINNER_FPCA_ROLE_RANKS'))
    if(is.null(tables[[name]]))iv_csv(data.frame(),file.path(dest,paste0(name,'.csv')))
  iv_csv(bind_rows(summary_rows),file.path(dest,'WINNER_SCIENCE.csv'))
  iv_csv(bind_rows(raw_rows),file.path(dest,'WINNER_RAW_COMPONENT_METRICS.csv'))
  iv_csv(bind_rows(interval_rows),file.path(dest,'WINNER_Q_INTERVAL_METRICS.csv'))
  iv_csv(bind_rows(historical_rows),file.path(dest,'WINNER_HISTORICAL_SUMMARY.csv'))
  pool <- bind_rows(lapply(arms,function(arm) {
    endpoint <- bind_rows(tables[[paste0(toupper(arm),'_ENDPOINT_SCIENCE')]])
    winner <- if(nrow(endpoint))endpoint[endpoint$fit_id==choices[[arm]]$selected_fit_id,,drop=FALSE]else data.frame()
    data.frame(data_id=id,arm=arm,registered=12L,valid=nrow(endpoint),failed_or_invalid=12L-nrow(endpoint),
      correct_candidates=if(nrow(endpoint))sum(endpoint$formal_nc_exact_loading_structure_08)else 0L,
      pool_has_correct_structure=if(nrow(endpoint))any(endpoint$formal_nc_exact_loading_structure_08)else FALSE,
      candidate_fraction_of_registered=if(nrow(endpoint))sum(endpoint$formal_nc_exact_loading_structure_08)/12 else 0,
      winner_exists=nrow(winner)==1L,winner_correct=if(nrow(winner))winner$formal_nc_exact_loading_structure_08 else FALSE,
      winner_shared=if(nrow(winner))winner$live_shared else NA_integer_,
      winner_specific_1=if(nrow(winner))winner$live_specific_1 else NA_integer_,
      winner_specific_2=if(nrow(winner))winner$live_specific_2 else NA_integer_)
  }))
  iv_csv(pool,file.path(dest,'POOL_SCIENCE.csv'))
  for(item in c(artifacts,winner_records))verify_record(item)
  iv_check(pair_seal_path,pair_seal_hash)
  iv_check(truth_binding$path,truth_binding$sha256)
  iv_check(test_truth_binding$path,test_truth_binding$sha256)
  add('sources_unchanged',identical(code_hash,iv_verify_package(ex))&&identical(protocol_hash,iv_sha(protocol_path)))
  writeLines(capture.output(sessionInfo()),file.path(dest,'EVALUATION_SESSION_INFO.txt'))
  list(status='SCIENCE_COMPLETE',data_id=id,registered=12L,registered_per_arm=12L,
    valid_by_arm=as.list(setNames(vapply(arms,function(arm)sum(valid$arm==arm),integer(1)),arms)),
    selected_fit_ids=lapply(choices,function(choice)choice$selected_fit_id),
    selection_frozen_before_truth=TRUE,truth_reads=truth_reads,truth_used_for_selection=FALSE,
    fit_updates=0L,pair_winners_sha256=pair_seal_hash,automatic_retries=0L,seconds=proc.time()[['elapsed']]-started)
},error=function(error) {
  list(status='FAIL',data_id=id,error=conditionMessage(error),truth_reads=truth_reads,
    fit_updates=0L,automatic_retries=0L,seconds=proc.time()[['elapsed']]-started)
})
if(length(checks))iv_csv(bind_rows(checks),file.path(dest,'CHECKS.csv'))
iv_json(status,file.path(dest,'STATUS.json'))
if(status$status=='SCIENCE_COMPLETE')iv_json(status,file.path(dest,'SCIENCE_COMPLETE.json'))
iv_json(lapply(list.files(dest,recursive=TRUE,full.names=TRUE),record),file.path(dest,'OUTPUT_MANIFEST.json'))
cat(status$status,id,'elapsed',status$seconds,'seconds\n')
if(status$status=='FAIL')quit(status=1L)
