# R-G12 next development route protocol v2 candidate

Date: 2026-08-29
Protocol ID: `R_G12_DEVELOPMENT_ROUTE_V2_CANDIDATE_20260829`
Status: approved development candidate for paper-pilot use
Formal v0L-V result: `FALSE`
Formal paper Monte Carlo authorized: `FALSE`

## 1. Purpose and inheritance

This candidate advances the post-v0L-V R route after the independent G12,
cross-scenario, stopping-semantics, Stage-4F/4G mechanism, and Stage-5A
correct-basin audits. It does not overwrite the 2026-08-24 v1 route, the
candidate2/v0L-V frozen protocols, or any historical result.

The statistical model, priors, CAVI updates, ordinary T=1 ELBO, bounded
pre-score update, Jaoua annealing schedule, loading-first identity evaluation,
and `lambda_orth=0` remain unchanged.

## 2. Registered R route

Route ID:

`random__gram_unit_energy__pre1__jaoua_default__multistart12__fixed400`

Every fit must explicitly use:

- `initialization = "random"`;
- `function_initialization = "gram_unit_energy"`;
- one bounded pre-score update;
- Jaoua annealing `c(1, 1.9, 100)`, i.e. 99 annealed sweeps;
- ordinary `T=1` CAVI after annealing;
- `lambda_orth = 0`;
- `n_cpus = 1` and one BLAS/OpenMP thread;
- the public `g12_stopping_control()` profile with `min_t1=396`,
  `consecutive=5`, and `max_t1=400`.

The public profile makes sweep 400 the first possible practical stop. Thus the
candidate records a common fixed-400 endpoint for scientific comparison. A
failure of the practical gate at sweep 400 is retained as a valid slow-case
status when the ordinary objective endpoint is otherwise eligible; it is not
relabeled as convergence and is not automatically continued.

## 3. Multistart and endpoint selection

- The main development and paper-pilot route uses 12 preregistered random fit
  seeds per data set.
- All starts are run independently to the registered fixed-400 endpoint. There
  is no live racing or truth-dependent cancellation.
- Objective eligibility and practical convergence are separate fields.
- A candidate is objective eligible only if the fit completed without error,
  the ordinary T=1 ELBO endpoint is finite, the required terminal object is
  complete, and no truth entered fitting, stopping, or selection.
- The unique data-level winner is the objective-eligible endpoint with maximum
  ordinary T=1 ELBO. A numerical tie is resolved by preregistered `fit_id`
  lexical order.
- Truth, structure labels, R/P/L, NRMSE, ISE, loading recovery, process
  recovery, or contribution recovery must not be used for stopping, endpoint
  eligibility, racing, or winner selection.

Fewer than 12 starts are allowed only in a separately versioned engineering
smoke or pilot protocol whose purpose is not to estimate scientific success.

## 4. Reporting contract

Every endpoint must retain its complete fit RDS, ordinary ELBO trace,
practical diagnostics, warnings/errors, elapsed time, peak memory when
available, objective eligibility, practical status, and COMPLETE/ERROR marker.

After truth-free winners are frozen and truth access is separately authorized,
the evaluation hierarchy is:

1. observed-time and dense-grid signal reconstruction;
2. loading-first one-to-one shared/specific identity, including missing,
   misplaced, duplicate, and extra components;
3. feature ISE, projection floor, estimate-to-projection excess ISE, and
   matched/total component coverage;
4. loading direction, support, and scale recovery;
5. factor-process and FPCA-score recovery;
6. complete-contribution and covariance-kernel/operator recovery;
7. ordinary ELBO, runtime, memory, warnings, objective status, practical
   status, and 380-to-400 scientific-output stability;
8. factor count and retained M only as auxiliary descriptions unless a
   separate dimension-selection experiment is registered.

Because recent audits found correct-direction scale sub-basins, reports must
include raw and canonical loading norms, factor scales, observed-to-dense
reconstruction gaps, process error, and complete-contribution error. These are
evaluation outputs, never winner-selection inputs.

## 5. Evidence-qualified interpretation

- G12 plus multistart and maximum ordinary ELBO is a supported mitigation for
  nonconvex reachability, not proof that optimization is solved.
- Practical convergence means local numerical change is small. It is not a
  structural or scientific-success criterion.
- Sweep 800 is not a default budget. It is allowed only in a separately
  preregistered, small mechanism audit and never silently replaces a frozen
  sweep-400 winner.
- Factor union-PPI or a correct factor count does not establish correct
  shared/specific identity. Loading-first matching and missing-aware recovery
  remain required.
- Evidence currently concerns the core `d=0`, mainly `S=2`, model. Covariate
  recovery, general dimension selection, larger `S/L/M`, uncertainty coverage,
  and broad computational scaling require separate validation.

## 6. Historical compatibility

The public function default remains `current_1_over_m`; this protocol opts in
to G12 explicitly. Historical scripts and results must not change behavior.
`current_1_over_m` remains available for reproduction and preregistered
ablation. It must not be merged with G12 into a truth-adaptive 24-start pool.

## 7. Authorization boundary

This file authorizes the route as a development candidate and its use in a
separately registered small paper-pilot. It does not authorize a full paper
Monte Carlo, a model/prior/ELBO change, a Git release tag, or retrospective
replacement of any v0L-V result.
