# Experimental, explicit factor-level Gamma shrinkage for active loading slabs.
# NULL is the original fixed precision-one model. All state is local to a fit.

.validate_factor_shrinkage_control <- function(control) {
  if (is.null(control)) return(NULL)
  if (!is.list(control) || is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!nzchar(names(control)))) {
    stop("factor_shrinkage must be NULL or a named control list.")
  }
  unknown <- setdiff(names(control), c("mode", "a1", "a2", "rate"))
  if (length(unknown)) stop("Unknown factor_shrinkage controls: ",
                            paste(unknown, collapse = ", "), ".")
  out <- list(mode = "mgp", a1 = 2.1, a2 = 3.1, rate = 1)
  out[names(control)] <- control
  if (!is.character(out$mode) || length(out$mode) != 1L ||
      is.na(out$mode) || !out$mode %in% c("mgp", "ard")) {
    stop("factor_shrinkage mode must be 'mgp' or 'ard'.")
  }
  for (key in c("a1", "a2", "rate")) {
    z <- out[[key]]
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0) {
      stop("factor_shrinkage ", key, " must be finite and positive.")
    }
    out[[key]] <- as.numeric(z)
  }
  out
}

.factor_shrinkage_prior_shapes <- function(L, control) {
  if (L == 0L) return(numeric())
  if (control$mode == "ard") return(rep(control$a1, L))
  c(control$a1, rep(control$a2, L - 1L))
}

.factor_shrinkage_initialize <- function(control, L_f, L_s, p) {
  control <- .validate_factor_shrinkage_control(control)
  if (is.null(control)) return(NULL)
  if (any(!is.finite(c(L_f, L_s, p))) || length(L_f) != 1L ||
      length(p) != 1L || p < 1 || p != as.integer(p) ||
      any(c(L_f, L_s) < 0 | c(L_f, L_s) != as.integer(c(L_f, L_s)))) {
    stop("Invalid factor_shrinkage model dimensions.")
  }
  chain <- function(L) list(
    shape = .factor_shrinkage_prior_shapes(L, control),
    rate = rep(control$rate, L))
  state <- list(version = 1L, control = control, p = as.integer(p),
                L_f = as.integer(L_f), L_s = as.integer(L_s),
                shared = chain(L_f), specific = lapply(L_s, chain))
  .factor_shrinkage_validate_state(state)
  state
}

.factor_shrinkage_validate_state <- function(state, L_f = NULL, L_s = NULL,
                                             p = NULL) {
  required <- c("version", "control", "p", "L_f", "L_s", "shared", "specific")
  if (!is.list(state) || !identical(names(state), required) ||
      !identical(state$version, 1L) || is.null(state$control) ||
      !identical(state$control,
                 .validate_factor_shrinkage_control(state$control))) {
    stop("Incomplete or unsupported factor_shrinkage state.")
  }
  dims <- c(state$L_f, state$L_s, state$p)
  if (length(state$L_f) != 1L || length(state$p) != 1L ||
      !is.numeric(dims) || any(!is.finite(dims)) ||
      any(dims != as.integer(dims)) || any(dims < 0) || state$p < 1L) {
    stop("Invalid factor_shrinkage state dimensions.")
  }
  if ((!is.null(L_f) && !identical(as.integer(L_f), state$L_f)) ||
      (!is.null(L_s) && !identical(as.integer(L_s), state$L_s)) ||
      (!is.null(p) && !identical(as.integer(p), state$p))) {
    stop("factor_shrinkage state does not match model dimensions.")
  }
  check_chain <- function(chain, L) {
    if (!is.list(chain) || !identical(names(chain), c("shape", "rate")) ||
        !is.numeric(chain$shape) || !is.numeric(chain$rate) ||
        length(chain$shape) != L || length(chain$rate) != L ||
        any(!is.finite(c(chain$shape, chain$rate))) ||
        any(c(chain$shape, chain$rate) <= 0)) {
      stop("Invalid Gamma density in factor_shrinkage state.")
    }
  }
  check_chain(state$shared, state$L_f)
  if (!is.list(state$specific) || length(state$specific) != length(state$L_s)) {
    stop("factor_shrinkage state needs one independent chain per study.")
  }
  for (s in seq_along(state$L_s)) check_chain(state$specific[[s]], state$L_s[s])
  invisible(TRUE)
}

.factor_shrinkage_chain_moments <- function(chain, mode) {
  means <- chain$shape / chain$rate
  logmeans <- digamma(chain$shape) - log(chain$rate)
  out <- if (mode == "mgp") {
    list(mean = exp(cumsum(log(means))), elog = cumsum(logmeans))
  } else list(mean = means, elog = logmeans)
  if (any(!is.finite(c(out$mean, out$elog))) || any(out$mean <= 0)) {
    stop("Non-finite factor precisions; no automatic clamping is permitted.")
  }
  out
}

.factor_shrinkage_moments <- function(state) {
  if (is.null(state)) return(NULL)
  .factor_shrinkage_validate_state(state)
  list(shared = .factor_shrinkage_chain_moments(state$shared, state$control$mode),
       specific = lapply(state$specific, .factor_shrinkage_chain_moments,
                         mode = state$control$mode))
}

# Validate explicit active-slab sufficient statistics, including empty blocks.
.factor_shrinkage_loading_statistics <- function(gamma, second, p, L) {
  if (L == 0L && is.null(gamma) && is.null(second)) {
    return(list(count = numeric(), energy = numeric()))
  }
  if (!is.matrix(gamma) || !is.matrix(second) ||
      !identical(dim(gamma), c(as.integer(p), as.integer(L))) ||
      !identical(dim(second), dim(gamma)) ||
      any(!is.finite(c(gamma, second))) || any(gamma < 0 | gamma > 1) ||
      any(second < 0) || any(second[gamma == 0] != 0)) {
    stop("Invalid inclusion probabilities or loading second moments.")
  }
  list(count = colSums(gamma), energy = colSums(second))
}

.factor_shrinkage_all_statistics <- function(state, gamma_a, term_a,
                                              gamma_b, term_b) {
  shared <- .factor_shrinkage_loading_statistics(
    gamma_a, term_a, state$p, state$L_f)
  if (sum(state$L_s) == 0L && is.null(gamma_b) && is.null(term_b)) {
    specific <- lapply(state$L_s, function(L) list(count = numeric(), energy = numeric()))
  } else {
    if (!is.list(gamma_b) || !is.list(term_b) ||
        length(gamma_b) != length(state$L_s) ||
        length(term_b) != length(state$L_s)) {
      stop("Missing study-specific loading states for factor_shrinkage.")
    }
    specific <- lapply(seq_along(state$L_s), function(s) {
      .factor_shrinkage_loading_statistics(
        gamma_b[[s]], term_b[[s]], state$p, state$L_s[s])
    })
  }
  list(shared = shared, specific = specific)
}

.factor_shrinkage_update_chain <- function(chain, statistics, control, c_val) {
  L <- length(chain$shape)
  prior_shape <- .factor_shrinkage_prior_shapes(L, control)
  for (h in seq_len(L)) {
    affected <- if (control$mode == "mgp") seq.int(h, L) else h
    weights <- vapply(affected, function(l) {
      if (control$mode == "ard") return(1)
      other <- setdiff(seq_len(l), h)
      exp(sum(log(chain$shape[other] / chain$rate[other])))
    }, numeric(1))
    shape <- 1 + c_val * (prior_shape[h] - 1 +
                           0.5 * sum(statistics$count[affected]))
    rate <- c_val * (control$rate +
                      0.5 * sum(weights * statistics$energy[affected]))
    if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
      stop("Gamma coordinate is not a finite proper density at this temperature.")
    }
    # Gauss-Seidel in chain order: later coordinates use this latest density.
    chain$shape[h] <- shape
    chain$rate[h] <- rate
  }
  chain
}

.factor_shrinkage_update <- function(state, gamma_a, term_a, gamma_b, term_b,
                                      c_val = 1) {
  if (is.null(state)) return(NULL)
  .factor_shrinkage_validate_state(state)
  if (!is.numeric(c_val) || length(c_val) != 1L || !is.finite(c_val) || c_val <= 0) {
    stop("factor_shrinkage inverse temperature must be finite and positive.")
  }
  statistics <- .factor_shrinkage_all_statistics(
    state, gamma_a, term_a, gamma_b, term_b)
  state$shared <- .factor_shrinkage_update_chain(
    state$shared, statistics$shared, state$control, c_val)
  for (s in seq_along(state$L_s)) {
    state$specific[[s]] <- .factor_shrinkage_update_chain(
      state$specific[[s]], statistics$specific[[s]], state$control, c_val)
  }
  .factor_shrinkage_validate_state(state)
  .factor_shrinkage_moments(state) # Fail on precision overflow before returning.
  state
}

.factor_shrinkage_elbo <- function(state, gamma_a, term_a, gamma_b, term_b,
                                    temperature = 1) {
  if (is.null(state)) return(numeric())
  .factor_shrinkage_validate_state(state)
  if (length(temperature) != 1L || !is.finite(temperature) || temperature <= 0) {
    stop("factor_shrinkage temperature must be finite and positive.")
  }
  statistics <- .factor_shrinkage_all_statistics(
    state, gamma_a, term_a, gamma_b, term_b)
  contribution <- function(chain, statistics) {
    moments <- .factor_shrinkage_chain_moments(chain, state$control$mode)
    slab_delta <- 0.5 * sum(statistics$count * moments$elog -
                             statistics$energy * (moments$mean - 1))
    a <- .factor_shrinkage_prior_shapes(length(chain$shape), state$control)
    b <- state$control$rate
    mean_delta <- chain$shape / chain$rate
    elog_delta <- digamma(chain$shape) - log(chain$rate)
    expected_log_prior <- a * log(b) - lgamma(a) +
      (a - 1) * elog_delta - b * mean_delta
    entropy <- chain$shape - log(chain$rate) + lgamma(chain$shape) +
      (1 - chain$shape) * digamma(chain$shape)
    c(slab_delta = slab_delta,
      gamma_prior_entropy = sum(expected_log_prior + temperature * entropy))
  }
  shared <- contribution(state$shared, statistics$shared)
  specific <- c(slab_delta = 0, gamma_prior_entropy = 0)
  for (s in seq_along(state$L_s)) {
    specific <- specific + contribution(state$specific[[s]], statistics$specific[[s]])
  }
  out <- c(factor_slab_precision_shared = unname(shared[1]),
           factor_slab_precision_specific = unname(specific[1]),
           factor_gamma_shared = unname(shared[2]),
           factor_gamma_specific = unname(specific[2]))
  if (any(!is.finite(out))) stop("Non-finite factor_shrinkage objective.")
  out
}

.factor_shrinkage_numeric_state <- function(state) {
  list(shared = state$shared, specific = state$specific)
}

.factor_shrinkage_diagnostics <- function(state, gamma_a, term_a, gamma_b,
                                           term_b, iteration, temperature) {
  statistics <- .factor_shrinkage_all_statistics(state, gamma_a, term_a, gamma_b, term_b)
  rows <- list()
  append_chain <- function(chain, stat, block, study) {
    L <- length(chain$shape)
    if (!L) return(NULL)
    moments <- .factor_shrinkage_chain_moments(chain, state$control$mode)
    data.frame(iteration = iteration, temperature = temperature,
               mode = state$control$mode, block = block, study = study,
               factor = seq_len(L), gamma_shape = chain$shape,
               gamma_rate = chain$rate, expected_precision = moments$mean,
               expected_log_precision = moments$elog,
               expected_support = stat$count, loading_second_energy = stat$energy)
  }
  rows[[1L]] <- append_chain(state$shared, statistics$shared, "shared", 0L)
  for (s in seq_along(state$L_s)) {
    rows[[length(rows) + 1L]] <- append_chain(
      state$specific[[s]], statistics$specific[[s]], "specific", s)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.factor_shrinkage_loading_precision <- function(moments, L) {
  if (is.null(moments)) return(NULL)
  if (!is.list(moments) || !identical(names(moments), c("mean", "elog")) ||
      !is.numeric(moments$mean) || !is.numeric(moments$elog) ||
      length(moments$mean) != L || length(moments$elog) != L ||
      any(!is.finite(c(moments$mean, moments$elog))) || any(moments$mean <= 0) ||
      any(moments$elog > log(moments$mean) + 1e-10)) {
    stop("Invalid factor-level expected precision and expected log precision.")
  }
  moments
}
