# New NC fields reuse named Stage7M definitions, with raw process/contribution additions.
nc_load_science <- function(exec_dir) {
  dep <- file.path(exec_dir, 'dependencies')
  check_hash(file.path(dep, 'stage7m_science.R'), '8fe1198843e7334aeaa3de6f3799d9e75d0f79436590be1eb10fc824e6c0d0d8')
  check_hash(file.path(dep, 'auc_tie_groups.R'), 'eb4740b78077ba4e5576eaf2aeb52a087d994aba784fbbab0eb3ead66bb914af')
  source(file.path(dep, 'stage7m_science.R'), local = .GlobalEnv)
  source(file.path(dep, 'auc_tie_groups.R'), local = .GlobalEnv)
  invisible(TRUE)
}

nc_prefix_fields <- function(x) {
  if (is.data.frame(x)) {
    names(x) <- sub('^stage7m_', 'nc_', names(x))
    return(x)
  }
  if (is.list(x)) return(lapply(x, nc_prefix_fields))
  x
}

nc_raw_role_metrics <- function(fit, obs, fixed, truth, historical) {
  tr <- truth$data$true_params; grid <- fixed$grid; w <- fixed$weights
  raw <- stage7m_grid_estimate(fit, fixed$basis, grid, 'raw', NULL)
  fun <- stage7m_truth_functions(tr, grid)
  rows <- historical$primary$functional_recovery
  rows <- rows[rows$scope == 'ppi_selected', , drop = FALSE]
  result <- list()
  for (k in seq_len(nrow(rows))) {
    item <- rows[k, , drop = FALSE]; s <- as.integer(item$study); l <- as.integer(item$true_factor)
    role <- as.character(item$true_role); matched <- isTRUE(item$matched_by_loading)
    truthload <- if (role == 'shared') tr$a_true[, l] else tr$b_true[[s]][, l]
    truefun <- if (role == 'shared') fun$shared[[l]] else fun$specific[[s]][[l]]
    truescores <- if (role == 'shared') tr$zeta_fpca_true[[s]][[l]] else tr$xi_fpca_true[[s]][[l]]
    if (matched) {
      block <- as.character(item$estimated_block); el <- as.integer(item$estimated_factor)
      oldl <- if (block == 'shared') fit$factor_order_shared[el] else fit$factor_order_specific[[s]][el]
      loading <- if (block == 'shared') raw$A[, oldl] else raw$B[[s]][, oldl]
      f <- if (block == 'shared') raw$shared[[oldl]] else raw$specific[[s]][[oldl]]
      scores <- if (block == 'shared') raw$scores_shared[[s]][[oldl]] else raw$scores_specific[[s]][[oldl]]
      loading_sign <- if (sum(loading * truthload) < 0) -1 else 1
    } else {
      block <- NA_character_; el <- oldl <- NA_integer_; loading <- rep(0, fit$p)
      f <- matrix(0, length(grid), 2L); scores <- matrix(0, fit$n_s[s], 2L); loading_sign <- 1
    }
    for (i in seq_len(fit$n_s[s])) {
      truepath <- as.vector(truefun %*% truescores[i, ])
      path <- as.vector(f %*% scores[i, ])
      estimates <- list(process = matrix(loading_sign * path, ncol = 1L), contribution = outer(path, loading))
      targets <- list(process = matrix(truepath, ncol = 1L), contribution = outer(truepath, truthload))
      for (object in names(estimates)) {
        er <- stage7m_curve_error(estimates[[object]], targets[[object]], w)
        result[[length(result) + 1L]] <- data.frame(study = s, subject = i, true_role = role,
          true_factor = l, estimated_block = block, estimated_factor = el, raw_factor_index = oldl,
          loading_status = item$loading_structure_status, matched_by_loading = matched,
          representation = 'raw_original', object = object, loading_sign = loading_sign,
          nc_mse = unname(er['mse']), nc_truth_energy = unname(er['truth_energy']),
          nc_rmse = unname(er['rmse']), nc_nrmse_missing_as_zero = unname(er['nrmse']),
          nc_nrmse_matched = if (matched) unname(er['nrmse']) else NA_real_,
          aggregation_weight = 1 / fit$S / fit$n_s[s],
          scale_convention = 'Sign and frozen factor permutation only; no truth regression rescaling')
      }
    }
  }
  do.call(rbind, result)
}

nc_endpoint_science <- function(fit, obs, fixed, truth, historical, pred, test_truth) {
  x <- stage7m_endpoint_science(fit, obs, fixed, truth, historical, pred, test_truth)
  x <- nc_prefix_fields(x)
  # d=0 has no beta estimand. Never present an identically-zero placeholder as recovery.
  beta_fields <- grep('beta', names(x$summary), value = TRUE)
  for (nm in beta_fields) x$summary[[nm]] <- NA_real_
  x$summary$nc_beta_metrics_reason <- 'not applicable: generated and fitted d=0'
  x$subject_metrics <- x$subject_metrics[x$subject_metrics$object != 'beta', , drop = FALSE]
  x$variable_metrics <- x$variable_metrics[x$variable_metrics$object != 'beta', , drop = FALSE]
  x$role_metrics$representation <- 'reported_original_FPCA'
  x$alignment_maps$shared_alignment_scope <- ifelse(x$alignment_maps$true_role == 'shared',
    'common shared function/component map across all studies when matched to shared block', 'study-specific map')
  raw_roles <- nc_raw_role_metrics(fit, obs, fixed, truth, historical)
  x$raw_role_metrics <- raw_roles
  for (object in c('process', 'contribution')) {
    z <- raw_roles[raw_roles$object == object, , drop = FALSE]
    mse <- sum(z$aggregation_weight * z$nc_mse)
    energy <- sum(z$aggregation_weight * z$nc_truth_energy)
    x$summary[[paste0('nc_', object, '_grid201_raw_missing_as_zero_nrmse')]] <- if (energy > 0) sqrt(mse / energy) else NA_real_
    good <- z$matched_by_loading
    x$summary[[paste0('nc_', object, '_grid201_raw_matched_nrmse')]] <-
      if (any(good) && sum(z$aggregation_weight[good] * z$nc_truth_energy[good]) > 0)
        sqrt(sum(z$aggregation_weight[good] * z$nc_mse[good]) /
          sum(z$aggregation_weight[good] * z$nc_truth_energy[good])) else NA_real_
  }
  x
}
