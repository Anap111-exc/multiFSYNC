# Audit-only adapters. Never modify a package namespace or frozen source.
# Pre-score executes WITHIN first CAVI sweep (after mu/beta), not at iteration 0.
# This gate captures G12-calibrated initial q before the first CAVI loop.
adapt_sha <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
adapt_function_text <- function(f) paste(c(deparse(formals(f), width.cutoff=500),
                                         deparse(body(f), width.cutoff=500)), collapse="\n")
adapt_function_sha <- function(f) digest::digest(adapt_function_text(f), algo="sha256", serialize=FALSE)

adapt_clone_functions <- function(original) {
  out <- new.env(parent=parent.env(original))
  for (nm in ls(original, all.names=TRUE)) {
    if (nm %in% c(".__NAMESPACE__.", ".__S3MethodsTable__.", ".packageName")) next
    x <- get(nm, original, inherits=FALSE)
    # Imported functions retain their own namespaces (notably splines::spline.des).
    if (is.function(x) && !is.primitive(x) && identical(environment(x),original)) environment(x) <- out
    if (is.environment(x)) {
      if (nm != ".multiFSYNC_spd_state") stop("Unexpected mutable source binding: ", nm)
      x <- list2env(as.list.environment(x, all.names=TRUE), parent=parent.env(x))
    }
    assign(nm, x, out)
  }
  out
}

adapt_load_baseline <- function(lib) {
  ns <- loadNamespace("multiFSYNC", lib.loc=lib)
  expected <- normalizePath(file.path(lib,"multiFSYNC"),winslash="/",mustWork=TRUE)
  actual <- normalizePath(getNamespaceInfo(ns,"path"),winslash="/",mustWork=TRUE)
  if (!identical(expected,actual)) stop("Already-loaded baseline namespace comes from a different library")
  # copy function values and rebind ONLY the copies, leaving locked namespace intact
  adapt_clone_functions(ns)
}

adapt_load_candidate <- function(root) {
  en <- new.env(parent=globalenv())
  en$spline.des <- splines::spline.des
  paths <- sort(list.files(file.path(root,"candidate","R"), "[.]R$", full.names=TRUE))
  stopifnot(length(paths)>0L)
  manifest <- file.path(root,"tests","CANDIDATE_SOURCE_SHA256SUMS.txt")
  lines <- readLines(manifest,warn=FALSE)
  expected_names <- substring(lines,67L)
  stopifnot(all(grepl("^[0-9a-f]{64}  candidate/R/[^/]+[.]R$",lines)),
            identical(sort(expected_names),sort(paste0("candidate/R/",basename(paths)))))
  for (i in seq_along(lines)) {
    observed <- digest::digest(file=file.path(root,expected_names[i]),algo="sha256",serialize=FALSE)
    if (!identical(observed,substr(lines[i],1L,64L))) stop("Candidate source hash mismatch: ",expected_names[i])
  }
  for (p in paths) sys.source(p, envir=en, keep.source=FALSE)
  en
}

adapt_is_loop <- function(expr) is.call(expr) && identical(expr[[1L]],as.name("for")) &&
  identical(expr[[2L]],as.name("i_iter"))

adapt_capture_snapshot <- function(frame) {
  available <- ls(frame, all.names=TRUE)
  qnames <- sort(unique(c(grep("^(mu_q_|Sigma_q_|inv_Sigma_q_|kappa_q_|lambda_q_|alpha_q_)", available, value=TRUE),
                         intersect(c("term_a","term_b_specific","c_1_omega_a","d_1_omega_a",
                                     "c_1_omega_b","d_1_omega_b"), available))))
  stopifnot(all(c("mu_q_nu_mu","mu_q_zeta","mu_q_xi","mu_q_a","mu_q_b_specific") %in% qnames))
  metadata <- c("S","n_s","p","d","L_f","L_s","M_f","M_s","K","C","C_g","Y","Z",
                "n_g","time_g","time_obs","list_hyper","list_cp_C","list_cp_C_Y","list_cp_Y",
                "sum_list_cp_C","total_obs_sj","inv_Sigma_beta","fit_seed","initialization",
                "initialization_control","initialization_diagnostics","c_val")
  metadata <- intersect(metadata, available)
  snap <- list(common_q=mget(qnames,frame,inherits=FALSE),
               fixed=mget(metadata,frame,inherits=FALSE),
               rng_kind=RNGkind(), rng_seed=get(".Random.seed",.GlobalEnv),
               factor_shrinkage_state=get0("factor_shrinkage_state",frame,inherits=FALSE))
  stop(structure(list(message="audit-only initial-state capture",call=NULL,snapshot=snap),
                 class=c("adapt_initial_capture","error","condition")))
}

adapt_capture_environment <- function(en) {
  out <- adapt_clone_functions(en)
  f <- out$bayesSYNC_multi_core
  expr <- as.list(body(f))
  anchor <- which(vapply(expr, adapt_is_loop, logical(1)))
  if (length(anchor)!=1L) stop("Expected one top-level CAVI for(i_iter ...) AST anchor")
  before <- paste(vapply(expr[seq_len(anchor-1L)],function(x) paste(deparse(x),collapse=" "),character(1)),collapse="\n")
  if (!grepl(".calibrate_random_factor_state",before,fixed=TRUE)) stop("Missing G12 calibration before anchor")
  # Inserting a call before the unique loop executes zero CAVI sweeps.
  body(f) <- as.call(append(expr,list(quote(.adapt_capture_initial(environment()))),after=anchor-1L))
  environment(f) <- out
  out$bayesSYNC_multi_core <- f
  out$.adapt_capture_initial <- adapt_capture_snapshot
  attr(out,"capture_anchor") <- list(position=anchor, loop=deparse(expr[[anchor]][1:3]),
                                   base_core_sha256=adapt_function_sha(en$bayesSYNC_multi_core))
  out
}

adapt_capture_initial <- function(en, args) {
  stopifnot(is.null(args$continuation_state), identical(args$initialization,"random"),
            identical(args$function_initialization,"gram_unit_energy"),
            identical(as.integer(args$pre_score_sweeps),1L))
  stopifnot(is.numeric(args$seed),length(args$seed)==1L,is.finite(args$seed),
            is.null(getOption("multiFSYNC.__stage7g_snapshot_20260928__")),
            is.null(getOption("multiFSYNC.__stage7i_resume_20260929__")))
  cap <- adapt_capture_environment(en)
  snapshot <- tryCatch({do.call(cap$bayesSYNC_multi_pre_score,args); stop("Capture did not intercept")},
                       adapt_initial_capture=function(e)e$snapshot)
  attr(snapshot,"capture_anchor") <- attr(cap,"capture_anchor")
  snapshot
}

adapt_initial_state_gate <- function(args, baseline, candidate, mgp=list(mode="mgp",a1=2.1,a2=3.1,rate=1),
                                     retain_state=FALSE, check_off=TRUE) {
  caller_kind <- RNGkind()
  caller_had_seed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  caller_seed <- if(caller_had_seed) get(".Random.seed",.GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind,as.list(caller_kind))
    if(caller_had_seed) assign(".Random.seed",caller_seed,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  args$factor_shrinkage <- NULL
  b <- adapt_capture_initial(baseline,args)
  states <- list(baseline=b)
  if (check_off) states$candidate_off <- adapt_capture_initial(candidate,args)
  args$factor_shrinkage <- mgp
  states$candidate_mgp <- adapt_capture_initial(candidate,args)
  rows <- lapply(names(states),function(nm) {
    s <- states[[nm]]
    data.frame(arm=nm, common_q_identical=identical(b$common_q,s$common_q),
      fixed_metadata_identical=identical(b$fixed,s$fixed),
      rng_identical=identical(b$rng_seed,s$rng_seed)&&identical(b$rng_kind,s$rng_kind),
      common_q_sha256=adapt_sha(s$common_q),fixed_sha256=adapt_sha(s$fixed),rng_sha256=adapt_sha(s$rng_seed),
      extra_mgp_state=!is.null(s$factor_shrinkage_state),stringsAsFactors=FALSE)
  })
  rows <- do.call(rbind,rows)
  pass <- all(rows$common_q_identical & rows$fixed_metadata_identical & rows$rng_identical) &&
    !rows$extra_mgp_state[rows$arm=="baseline"] && rows$extra_mgp_state[rows$arm=="candidate_mgp"]
  list(pass=pass, rows=rows, stage="after_random_and_G12_before_first_CAVI_loop",
       original_CAVI_sweeps=0L, pre_score_executed=FALSE,
       caller_rng_restore_guard=TRUE,
       anchors=lapply(states,function(s)attr(s,"capture_anchor")),
       state=if(retain_state)states else NULL)
}

adapt_sourcebench <- function(baseline,candidate) {
  common <- intersect(ls(baseline,all.names=TRUE),ls(candidate,all.names=TRUE))
  common <- sort(common[vapply(common,function(n)is.function(baseline[[n]])&&is.function(candidate[[n]]),logical(1))])
  rows <- lapply(common,function(n) data.frame(name=n,baseline_sha256=adapt_function_sha(baseline[[n]]),
    candidate_sha256=adapt_function_sha(candidate[[n]]),stringsAsFactors=FALSE))
  rows <- do.call(rbind,rows); rows$identical_canonical <- rows$baseline_sha256==rows$candidate_sha256
  list(functions=rows,changed=rows$name[!rows$identical_canonical],
       new_functions=setdiff(ls(candidate,all.names=TRUE),ls(baseline,all.names=TRUE)))
}
