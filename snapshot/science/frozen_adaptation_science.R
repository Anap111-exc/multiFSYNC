# Evaluation adapters only. No model parameter update or candidate selection by truth.
`%||%` <- function(x, y) if (is.null(x)) y else x
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
sha <- function(path) digest::digest(file = path, algo = 'sha256', serialize = FALSE)
check_hash <- function(path, expected) {
  assert(file.exists(path) && identical(tolower(sha(path)), tolower(expected)), paste('Hash mismatch:', path))
}
extract_functions <- function(path, target, envir = .GlobalEnv) {
  found <- character()
  for (expr in parse(path, keep.source = FALSE)) {
    if (!is.call(expr) || !identical(expr[[1L]], as.name('<-')) ||
        !is.symbol(expr[[2L]]) || !is.call(expr[[3L]]) ||
        !identical(expr[[3L]][[1L]], as.name('function'))) next
    nm <- as.character(expr[[2L]])
    if (nm %in% target) {eval(expr, envir); found <- c(found, nm)}
  }
  assert(setequal(found, target) && !anyDuplicated(found), paste('Function extraction:', path))
}
adapt_new_csv <- function(x, path) {
  assert(!file.exists(path), paste('Refusing output overwrite:', path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = '')
}
adapt_new_rds <- function(x, path) {
  assert(!file.exists(path), paste('Refusing output overwrite:', path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = 'xz', version = 3L)
}
adapt_new_json <- function(x, path) {
  assert(!file.exists(path), paste('Refusing output overwrite:', path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = 'null', digits = 17)
}
adapt_read_bound <- function(path, expected) {check_hash(path, expected); readRDS(path)}
adapt_is_true <- function(x) tolower(as.character(x)) %in% c('true', '1')
adapt_load_evaluator <- function(exec_root) {
  sr <- file.path(exec_root, 'science'); root <- file.path(sr, 'frozen')
  manifest <- read.csv(file.path(sr, 'FROZEN_DEPENDENCIES.csv'), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) check_hash(file.path(sr, manifest$relative_path[i]), manifest$sha256[i])
  source(file.path(root, 'dependencies', 'stage7m_reconstruct.R'), local = .GlobalEnv)
  source(file.path(root, 'dependencies', 'stage7m_science.R'), local = .GlobalEnv)
  source(file.path(root, 'dependencies', 'auc_tie_groups.R'), local = .GlobalEnv)
  source(file.path(root, 'nc_evaluation_setup.R'), local = .GlobalEnv)
  source(file.path(root, 'nc_science.R'), local = .GlobalEnv)
  source(file.path(root, 'band_moments.R'), local = .GlobalEnv)
  nc_evaluation_setup(root)
}
adapt_evaluation_truth <- function(sealed, obs, fixed) {
  assert(identical(sealed$data_id, obs$data_id), 'Truth/observation identity differs')
  truth <- sealed
  truth$data$Y <- obs$Y; truth$data$Z <- obs$Z
  truth$data$time_obs <- obs$time_obs; truth$data$C <- fixed$C
  truth
}
adapt_historical <- function(fit, obs, fixed, truth, env) {
  assert(fit$S == 2L && fit$L_f == 1L && all(unlist(fit$L_s) == 1L),
         'This adapter is registered only for the correctly specified S2 L=(1,1,1) experiment')
  assert(fit$d == obs$dimensions$d && fit$d == truth$data$true_params$d_use, 'Covariate dimension mismatch')
  assert(max(abs(stage7m_C(fit$time_g, fixed$basis) - fit$C_g)) <= 1e-12,
         'Frozen Linux basis does not match fit: do not regenerate a platform basis')
  old <- env$estimated_mean_dense
  on.exit(assign('estimated_mean_dense', old, envir = env), add = TRUE)
  env$estimated_mean_dense <- function(fit, data, study) {
    stage7m_C(data$true_params$t_grid_dense, fixed$basis) %*%
      do.call(cbind, fit$mu_q_nu_mu_original[[study]])
  }
  s6h_evaluate_fit(fit, truth, env)
}
adapt_structure <- function(fit, truth, historical, env) {
  x <- nc_general_structure(fit, truth, historical, env)
  x$summary$nc_beta_metrics_applicability <- if (fit$d) 'nonzero covariate C1: beta evaluated separately' else 'NA: generated d=0 and fitted d=0'
  x$summary$effective_M_shared <- paste(unlist(fit$list_effective_M), collapse = ';')
  x$summary$effective_M_specific <- paste(vapply(fit$list_effective_M_spec, function(z) paste(unlist(z), collapse = ','), character(1)), collapse = ';')
  x$summary$factor_ppi_threshold <- 0.5
  x$summary$factor_ppi_comparison <- '>=' 
  x
}
adapt_raw_metrics <- function(fit, obs, fixed, truth) {
  tr <- truth$data$true_params; tb <- stage7m_truth_basis(tr)
  grid <- fixed$grid; w <- fixed$weights
  assert(length(grid) == 201L && max(abs(grid-seq(0,1,length.out=201L))) < 1e-12 && abs(sum(w)-1) < 1e-12, 'Grid201 contract')
  ft <- stage7m_truth_functions(tr, grid, tb)
  cache <- stage7m_grid_estimate(fit, fixed$basis, grid, 'raw', NULL)
  records <- list()
  for (s in seq_len(fit$S)) for (i in seq_len(fit$n_s[s])) for (scope in c('grid201', 'observed')) {
    tt <- if (scope == 'grid201') grid else obs$time_obs[[s]][[i]]
    ww <- if (scope == 'grid201') w else rep(1/length(tt), length(tt))
    target <- stage7m_true_parts(tr, obs, s, i, tt, tb, if (scope == 'grid201') ft else NULL)
    if (scope == 'grid201') estimate <- stage7m_grid_subject(cache, obs, s, i) else {
      estimate <- stage7m_estimate(fit, fixed$basis, s, i, tt, 'raw')
      bc <- matrix(0, length(tt), fit$p)
      for (r in seq_len(fit$d)) bc <- bc + obs$Z[[s]][i,r] * estimate$beta[[r]]
      estimate$beta <- bc; estimate$signal <- estimate$mean + estimate$factor + bc
    }
    objects <- c('shared_contribution','specific_contribution','factor','mean','signal')
    if (fit$d > 0L) objects <- c(objects, 'beta')
    for (object in objects) {
      er <- stage7m_curve_error(estimate[[object]], target[[object]], ww)
      records[[length(records)+1L]] <- data.frame(study=s,subject=i,scope=scope,object=object,
        representation='raw_all_estimated_blocks',mse=unname(er['mse']),truth_energy=unname(er['truth_energy']),
        rmse=unname(er['rmse']),nrmse=unname(er['nrmse']),weight=1/fit$S/fit$n_s[s])
    }
  }
  detail <- do.call(rbind, records)
  keys <- unique(detail[,c('scope','object')]); sums <- vector('list',nrow(keys))
  for (k in seq_len(nrow(keys))) {
    z <- detail[detail$scope==keys$scope[k] & detail$object==keys$object[k],]
    mse <- sum(z$weight*z$mse); energy <- sum(z$weight*z$truth_energy)
    sums[[k]] <- data.frame(scope=keys$scope[k],object=keys$object[k],mse=mse,truth_energy=energy,
      rmse=sqrt(mse),nrmse=if(energy>0)sqrt(mse/energy)else NA_real_,
      aggregation='equal_study_then_equal_subject_then_equal_variable')
  }
  list(subject=detail,summary=do.call(rbind,sums))
}
adapt_interval_metrics <- function(fit, obs, fixed, truth, test_truth) {
  pack <- list(estimate=fit,observation=obs,basis=fixed$basis)
  tr <- truth$data$true_params; tb <- stage7m_truth_basis(tr)
  t <- fixed$grid; w <- fixed$weights; z <- qnorm(.975)
  tf <- stage7m_truth_functions(tr,t,tb); cache <- band_cache(pack,t)
  records <- list()
  for(s in seq_len(fit$S))for(i in seq_len(fit$n_s[s]))for(scope in c('grid201','observed','test')) {
    tt <- if(scope=='grid201')t else if(scope=='observed')obs$time_obs[[s]][[i]] else test_truth$time[[s]][[i]]
    ww <- if(scope=='grid201')w else rep(1/length(tt),length(tt))
    cc <- if(scope=='grid201')cache else band_cache(pack,tt)
    aa <- band_subject(pack,cc,s,i)
    target <- stage7m_true_parts(tr,obs,s,i,tt,tb,if(scope=='grid201')tf else NULL)
    targets <- target[c('signal','factor','shared_contribution','specific_contribution','mean')]
    if(fit$d)targets$beta_contribution <- target$beta
    if(scope!='grid201')targets$response <- if(scope=='observed')do.call(cbind,obs$Y[[s]][[i]])else test_truth$Y[[s]][[i]]
    for(object in names(targets)) {
      met <- band_summary(aa$mean[[object]],aa$variance[[object]],targets[[object]],ww,z)
      records[[length(records)+1L]] <- data.frame(study=s,subject=i,scope=scope,object=object,
        interval_type=if(object=='response')'predictive_q_moment'else'latent_q_moment',level=.95,
        weight=1/fit$S/fit$n_s[s],t(met),check.names=FALSE)
    }
  }
  for(r in seq_len(fit$d)) {
    met <- band_summary(cache$beta[[r]],cache$beta_var[[r]],tf$beta[[r]],w,z)
    records[[length(records)+1L]] <- data.frame(study=0,subject=r,scope='grid201',object=paste0('beta_function_',r),
      interval_type='gaussian_coefficient_q',level=.95,weight=1,t(met),check.names=FALSE)
  }
  detail <- do.call(rbind,records)
  keys <- unique(detail[,c('scope','object','interval_type')]); sums <- vector('list',nrow(keys))
  for(k in seq_len(nrow(keys))) {
    a <- detail[detail$scope==keys$scope[k] & detail$object==keys$object[k],]
    sums[[k]] <- cbind(keys[k,,drop=FALSE],level=.95,
      coverage=sum(a$weight*a$coverage)/sum(a$weight),mean_width=sum(a$weight*a$width)/sum(a$weight),
      all_grid_covered=sum(a$weight*a$all_grid_covered)/sum(a$weight),
      interpretation='pointwise Gaussian moment interval under variational q; not a simultaneous band')
  }
  list(subject=detail,summary=do.call(rbind,sums))
}
adapt_winner_science <- function(fit, obs, fixed, truth, test_truth, prediction, env, historical=NULL) {
  if(is.null(historical))historical <- adapt_historical(fit,obs,fixed,truth,env)
  extended <- stage7m_endpoint_science(fit,obs,fixed,truth,historical,prediction,test_truth)
  if(fit$d==0L) {
    for(nm in grep('beta',names(extended$summary),value=TRUE))extended$summary[[nm]] <- NA_real_
    extended$subject_metrics <- extended$subject_metrics[extended$subject_metrics$object!='beta',,drop=FALSE]
    extended$variable_metrics <- extended$variable_metrics[extended$variable_metrics$object!='beta',,drop=FALSE]
  }
  list(structure=adapt_structure(fit,truth,historical,env),historical=historical,extended=extended,
       raw=adapt_raw_metrics(fit,obs,fixed,truth),
       raw_roles=nc_raw_role_metrics(fit,obs,fixed,truth,historical),
       intervals=adapt_interval_metrics(fit,obs,fixed,truth,test_truth))
}
