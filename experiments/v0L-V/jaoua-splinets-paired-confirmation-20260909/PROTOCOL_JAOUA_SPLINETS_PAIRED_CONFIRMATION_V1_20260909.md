# Jaoua/Splinets truth-generator paired confirmation protocol V1

Date: 2026-09-09

Status: registered development sensitivity experiment. This is neither a
formal v0L-V result nor a replacement for any frozen experiment.

## Question

Does changing only the simulation truth-function construction from the legacy
dense-grid Gram orthogonalisation to the upstream Jaoua/bayesSYNC Splinets
construction materially change multiFSYNC optimisation reachability or later
scientific recovery under the already selected G12 fitting route?

The experiment is a generator sensitivity check. It must not choose the
generator that makes multiFSYNC look better. Numerical/specification failure is
a valid reason not to adopt an arm; a recovery difference is evidence of model
sensitivity and must be reported as such.

## Paired design

- truth basis arms: `legacy_dense_gram` and `jaoua_splinets`;
- three new registered logical data seeds;
- twelve registered random fitting starts per logical data seed and arm;
- 2 x 3 x 12 = 72 G12 fits;
- a separate `truth_basis_seed` isolates basis construction and restores the
  main RNG stream, so paired arms have identical observation times, mean
  phases, covariates, scores, loadings, sparsity supports and measurement
  errors;
- the two arms may differ only in truth eigenfunctions, resulting factor
  processes, signals and observed responses;
- paired arms use the same fit seed for each logical data/start pair.

## Fixed data and model specification

`S=2`, `n_s=(30,30)`, `p=500`, `d=0`, true `L_f=1`, true
`L_s=(1,1)`, `M_f=2`, `M_s=(2,2)`, `K=5`, identified strong loading
geometry, fixed 10% active loadings, irregular 6--9 observation times,
`sigma_eps=0.3`, `mean_amp=0.6`, and score eigenvalue weights proportional to
`1/m^2`. The established v0F constant-per-active-loading calibration is
applied independently and deterministically to both arms.

## Fixed fitting route

Every fit uses random initialisation, Gram unit-function-energy calibration
(`G12`), one pre-score sweep, Jaoua annealing `c(1,1.9,100)`, and exactly 400
ordinary `T=1` sweeps. Each fit has `n_cpus=1`; only independent fits are
parallelised. There is no racing, continuation, automatic 800-sweep extension,
model/ELBO change, or alternative initialiser.

Within each data-arm stratum, the truth-free winner is the objective-eligible
endpoint with maximum finite ordinary `T=1` ELBO, with lexical `fit_id` as the
predeclared tie breaker. ELBO is never compared across basis arms or datasets.

## Execution gate

First run the eight registered fits formed by logical data seed 1, both arms,
and starts 1--4. They count toward the 72 fits. Continue to the remaining fits
only if all eight have complete terminal records, no fit error, an eligible
ordinary ELBO and verified single-core route provenance. A failed fit is
retained and the supervisor stops; it is not selectively rerun.

The target Linux server may use seven outer workers. Estimated compute is about
90 CPU-hours and 15--18 wall-clock hours, conditional on prior L1M2 timings.

## Truth isolation and endpoint

Generation writes observation-only bundles and separately sealed truth. Fitting,
stopping and winner selection cannot read sealed truth. After all 72 endpoints
and six winners are frozen, execution stops and waits for a separate explicit
truth-unseal/evaluation authorization. No scientific success threshold is
tested before unsealing.

Later paired evaluation will use the existing loading-first framework:
reconstruction, feature/kernel recovery, loading recovery, FPCA scores, factor
processes/contributions, structure errors, ELBO/runtime/warnings and convergence
descriptors. No new success criterion is introduced by this protocol.

