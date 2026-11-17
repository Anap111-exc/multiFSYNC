# Exact storage-only replacement for the frozen scale diagnostic.
# The original loop, matrix products, additions and returned fields are unchanged.
# Pieces retain their original names; one final unlist reproduces default c()
# concatenation, including the original leading numeric() type seed.
diagnostic_function_sha256 <- function(fn) {
  stopifnot(is.function(fn))
  text <- paste(c(deparse(formals(fn), width.cutoff = 500L),
                  deparse(body(fn), width.cutoff = 500L)), collapse = "\n")
  digest::digest(text, algo = "sha256", serialize = FALSE)
}

install_diagnostic_preallocation_patch <- function(en) {
  stopifnot(is.environment(en),
            exists(".scale_trace_mean_contributions", envir = en, inherits = FALSE))
  original <- get(".scale_trace_mean_contributions", envir = en, inherits = FALSE)
  expected <- "8222f39307c7ea607bc886b9809921bdf618e3dc2e2d7e3a669b010290828000"
  stopifnot(identical(diagnostic_function_sha256(original), expected))
  text <- paste(deparse(body(original), width.cutoff = 500L), collapse = "\n")
  replace_once <- function(from, to) {
    positions <- gregexpr(from, text, fixed = TRUE)[[1L]]
    stopifnot(length(positions) == 1L, positions[[1L]] > 0L)
    text <<- sub(from, to, text, fixed = TRUE)
  }
  replace_once("observed <- numeric()", paste0(
    ".trace_chunk_count <- state$p * sum(state$n_s[seq_len(state$S)]); ",
    ".trace_chunk_index <- 1L; ",
    "observed <- vector(\"list\", .trace_chunk_count + 1L); ",
    "observed[[1L]] <- numeric()"))
  for (name in c("mean_part", "beta_part", "shared_part", "specific_part")) {
    replace_once(paste0(name, " <- numeric()"), paste0(
      name, " <- vector(\"list\", .trace_chunk_count + 1L); ",
      name, "[[1L]] <- numeric()"))
  }
  replace_once("observed_one <- as.numeric(Y[[s]][[i]][[j]])", paste0(
    ".trace_chunk_index <- .trace_chunk_index + 1L; ",
    "observed_one <- as.numeric(Y[[s]][[i]][[j]])"))
  for (name in c("observed", "mean_part", "beta_part", "shared_part", "specific_part")) {
    one <- if (name == "observed") "observed_one" else
      paste0(sub("_part$", "", name), "_one")
    replace_once(paste0(name, " <- c(", name, ", ", one, ")"),
                 paste0(name, "[.trace_chunk_index] <- list(", one, ")"))
  }
  target <- "fitted <- mean_part + beta_part + shared_part + specific_part"
  flatten <- paste(vapply(c("observed", "mean_part", "beta_part", "shared_part", "specific_part"),
    function(name) paste0(name, " <- unlist(", name,
      ", recursive = FALSE, use.names = TRUE)"), character(1L)), collapse = "; ")
  replace_once(target, paste0(flatten, "; ", target))
  replacement <- original
  body(replacement) <- parse(text = text, keep.source = FALSE)[[1L]]
  replacement <- compiler::cmpfun(replacement)
  stopifnot(identical(formals(replacement), formals(original)),
            identical(environment(replacement), environment(original)))
  assign(".scale_trace_mean_contributions", replacement, envir = en)
  invisible(list(patch_id = "diagnostic_piece_preallocation_v1",
    function_name = ".scale_trace_mean_contributions",
    original_function_sha256 = expected,
    replacement_function_sha256 = diagnostic_function_sha256(replacement),
    changed_estimation_functions = 0L, diagnostics_disabled = FALSE))
}
