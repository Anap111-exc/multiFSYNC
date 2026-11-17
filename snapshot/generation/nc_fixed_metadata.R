nc_canonical_function_sha256 <- function(fn) {
  # Source references include platform paths/environments and are diagnostic metadata.
  # Preserve the complete formals/body syntax in UTF-8 LF text, as independently dumped.
  lines <- c(paste0('FORMALS\n', paste(suppressWarnings(deparse(formals(fn), width.cutoff = 500L, control = 'all')), collapse = '\n')),
    paste0('BODY\n', paste(suppressWarnings(deparse(body(fn), width.cutoff = 500L, control = 'all')), collapse = '\n')))
  bytes <- charToRaw(enc2utf8(paste0(paste(lines, collapse = '\n'), '\n')))
  digest::digest(bytes, algo = 'sha256', serialize = FALSE)
}

nc_fixed_metadata <- function(obs) {
  a <- obs$dimensions
  grid <- multiFSYNC:::get_grid_objects(unlist(obs$time_obs, recursive = FALSE),
    a$K, n_g = 51L, format_univ = TRUE)
  n_s <- as.integer(vapply(obs$Y, length, integer(1)))
  end <- cumsum(n_s); first <- c(1L, head(end, -1L) + 1L)
  C <- Map(function(x, y) grid$C[x:y], first, end)
  basis <- stage7m_basis(obs$time_obs, a$K)
  actual_columns <- as.integer(ncol(grid$C_g))
  assert(actual_columns == a$K + 2L, 'Actual spline columns differ from bound K+2')
  assert(max(abs(stage7m_C(grid$time_g, basis) - grid$C_g)) < 1e-12,
         'Standalone and installed fit basis differ')
  list(schema = 'nc_fixed_train_only_metadata_v1', basis = basis,
    grid = seq(0, 1, length.out = 201L), weights = stage7m_weights(seq(0, 1, length.out = 201L)),
    C = C, fit_metadata = list(S = as.integer(a$S), n_s = n_s, p = as.integer(a$p), d = 0L,
      L_f = as.integer(a$L_f), L_s = as.integer(a$L_s), M_f = as.integer(a$M_f), M_s = a$M_s,
      K = as.integer(a$K), K_total = actual_columns, C_g = grid$C_g, time_g = grid$time_g,
      list_cp_C = lapply(C, function(study) lapply(study, crossprod)),
      list_hyper = list(Sigma_beta = multiFSYNC:::set_hyper(d_0 = as.integer(a$p))$Sigma_beta)),
    aggregation = list(study_weight = 1 / a$S, subject_weight = 1 / n_s,
      variable_weight = 1 / a$p), truth_read = FALSE)
}
