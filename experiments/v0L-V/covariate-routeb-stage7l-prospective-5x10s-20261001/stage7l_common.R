options(warn = 1, stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
check_hash <- function(path, expected) {
  assert(file.exists(path) && identical(tolower(sha(path)), tolower(expected)),
         paste("Binding conflict:", path))
}
new_csv <- function(value, path) {
  assert(!file.exists(path), paste("Refusing overwrite:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(value, path, row.names = FALSE, na = "")
}
new_rds <- function(value, path) {
  assert(!file.exists(path), paste("Refusing overwrite:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path, version = 3L)
}
new_lines <- function(value, path) {
  assert(!file.exists(path), paste("Refusing overwrite:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(value, path, useBytes = TRUE)
}
read_csv <- function(path) read.csv(path, check.names = FALSE)
extract_functions <- function(path, target, envir = .GlobalEnv) {
  found <- character()
  for (expr in parse(path, keep.source = FALSE)) {
    if (!is.call(expr) || !identical(expr[[1L]], as.name("<-")) ||
        !is.symbol(expr[[2L]]) || !is.call(expr[[3L]]) ||
        !identical(expr[[3L]][[1L]], as.name("function"))) next
    nm <- as.character(expr[[2L]])
    if (nm %in% target) { eval(expr, envir); found <- c(found, nm) }
  }
  assert(setequal(found, target) && !anyDuplicated(found), "Function extraction failed.")
}
load_stage7l <- function(exec_dir, lib) {
  .libPaths(c(lib, Sys.getenv("STAGE7L_DGP_LIB",
    "/root/multiFSYNC-jaoua-paired-lib-20260909"), .libPaths()))
  library(multiFSYNC)
  assert(requireNamespace("digest", quietly = TRUE), "Missing digest.")
  check_hash(file.path(exec_dir, "dependencies", "pilot.R"),
             "5eaa571692b497f3abbf0c4a3083edaa11d9547d6d7f24d972c9a14e2458729d")
  check_hash(file.path(exec_dir, "dependencies", "uas_helper.R"),
             "03ae170103e2d70e63889b68553d9d13261233edfa1ec775ce185de39f840c99")
  extract_functions(file.path(exec_dir, "dependencies", "pilot.R"),
    c("generator_config", "calibrate_loadings", "bspline_design", "paired_beta_data",
      "pair_checks", "save_arm", "objective_status"))
  source(file.path(exec_dir, "dependencies", "uas_helper.R"), local = .GlobalEnv)
  assign("experiment_id", "COVARIATE_ROUTEB_STAGE7L_V1_20261001", .GlobalEnv)
  assign("save_new", new_rds, .GlobalEnv)
  assign("csv_new", new_csv, .GlobalEnv)
  assign("lines_new", new_lines, .GlobalEnv)
  assign("%||%", function(x,y) if (is.null(x)) y else x, .GlobalEnv)
}
fit_arguments <- function(obs, fit_seed) {
  control <- multiFSYNC::g12_stopping_control()
  control$checkpoints <- c(380L, 400L)
  assert(control$min_t1 == 396L && control$consecutive == 5L &&
           control$max_t1 == 400L, "G12 fixed-horizon binding changed.")
  list(Y = obs$Y, Z = obs$Z, time_obs = obs$time_obs,
       L_f = 1L, L_s = c(1L,1L), M_f = 2L, M_s = list(2L,2L),
       K = obs$dimensions$K, n_g = 51L, anneal = c(1,1.9,100), maxit = 499L,
       n_cpus = 1L, verbose = FALSE, seed = as.integer(fit_seed),
       tol_abs = 1e-3, tol_rel = 1e-5, bool_scale = FALSE,
       bool_var_spec_prob = FALSE, d_0 = as.integer(obs$dimensions$p),
       convergence_rule = "practical", lambda_orth = 0, practical_control = control,
       initialization = "random", initialization_control = list(
         grid_size = 81L, perturb_sd = 0.05, rank_tol = 1e-8),
       continuation_state = NULL, pre_score_sweeps = 1L,
       trace_sweeps = c(1L,380L,400L), function_initialization = "gram_unit_energy")
}
fixed_metadata <- function(obs) {
  grid <- multiFSYNC:::get_grid_objects(
    unlist(obs$time_obs, recursive = FALSE), obs$dimensions$K,
    n_g = 51L, format_univ = TRUE)
  n_s <- as.integer(vapply(obs$Y, length, integer(1)))
  end <- cumsum(n_s); first <- c(1L, head(end, -1L)+1L)
  C <- Map(function(a,b) grid$C[a:b], first, end)
  meta <- list(S = 2L, L_f = 1L, L_s = c(1L,1L), M_s = list(2L,2L),
       K_total = as.integer(obs$dimensions$K+2L), K = obs$dimensions$K,
       p = obs$dimensions$p, n_s = n_s, d = 2L,
       C_g = grid$C_g, time_g = grid$time_g,
       list_cp_C = lapply(C, function(study) lapply(study, crossprod)),
       list_hyper = list(Sigma_beta = multiFSYNC:::set_hyper(
         d_0 = as.integer(obs$dimensions$p))$Sigma_beta))
  list(fit_metadata = meta, C = C)
}
all_finite <- function(x) {
  if (is.numeric(x)) return(all(is.finite(x)))
  if (is.list(x)) return(all(vapply(x, all_finite, logical(1))))
  TRUE
}
make_snapshot_sink <- function(sweeps, output_dir) {
  sink <- new.env(parent = emptyenv())
  sink$schema <- "stage7g_snapshot_v1"; sink$sweeps <- as.integer(sweeps)
  saved <- list()
  makeActiveBinding("records", function(value) {
    if (missing(value)) return(saved)
    assert(length(value) == length(saved)+1L, "Unexpected snapshot mutation.")
    state <- tail(value, 1L)[[1L]]
    assert(all_finite(state), "Nonfinite raw snapshot.")
    new_rds(state, file.path(output_dir, sprintf("STATE_SWEEP%03d.rds",state$iteration)))
    cat("SNAPSHOT_SAVED sweep=", state$iteration, "\n", sep="")
    saved <<- value
  }, sink)
  sink
}
terminal_columns <- c("fit_id","data_id","branch","start_index","fit_seed",
  "source_U_fit_id","status","error","elapsed_seconds","objective_eligible",
  "invalid_reasons","final_ordinary_elbo","global_iteration","ordinary_t1_sweeps",
  "actual_updates","strict_practical_converged","fit_path","snapshot50_path")
terminal_row <- function(row, branch, source_U = "") {
  data.frame(fit_id=row$fit_id, data_id=row$data_id, branch=branch,
    start_index=row$start_index,fit_seed=row$fit_seed,source_U_fit_id=source_U,
    status="ERROR",error="",elapsed_seconds=0,objective_eligible=FALSE,
    invalid_reasons="fit_error",final_ordinary_elbo=NA_real_,
    global_iteration=NA_integer_,ordinary_t1_sweeps=NA_integer_,actual_updates=0L,
    strict_practical_converged=FALSE,fit_path="",snapshot50_path="")
}
