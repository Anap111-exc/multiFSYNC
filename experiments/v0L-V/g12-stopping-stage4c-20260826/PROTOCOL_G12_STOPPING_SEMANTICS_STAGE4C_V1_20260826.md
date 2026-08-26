# G12 stopping semantics stage 4C protocol v1

Registration time: 2026-08-26T08:22:45Z
Status: authorized development experiment; truth remains sealed

## Question

Stage 4C independently checks whether the public opt-in
`g12_stopping_control()` route can execute safely at full size and whether its
fixed-400 endpoints are locally stable from sweep 380 to 400. It is not a new
method comparison, not a truth-recovery evaluation, and not a formal v0L-V
result.

## Registered data and fits

- Six new data seeds: `82726001` through `82726006`.
- Three scenarios, two data sets per scenario: baseline strong identification,
  weak shared-specific loading separation (full support overlap and target
  absolute cosine 0.6), and sparse irregular observation.
- Fixed generator settings: `S=2`, `n_s=(30,30)`, `p=500`, `d=0`, true
  `L_f=1`, true `L_s=(1,1)`, `M=2`, `K=5`, `sigma_eps=0.3`, and
  `mean_amp=0.6`. Baseline and weak scenarios use 6-9 observations; sparse
  scenarios use 3-5.
- Two G12 starts per data set, for 12 fits. Fit seeds are
  `87266011/12, 87266021/22, ..., 87266061/62`.
- Each fit uses one CPU. Four independent fits may run concurrently.

Every fit uses random initialization, Gram-unit-energy function calibration,
one bounded pre-score sweep, Jaoua annealing `c(1,1.9,100)`, and the public G12
stopping profile. The profile has `min_t1=396`, `consecutive=5`,
`max_t1=400`, and the quantile-factor gate. Therefore the first possible
practical stop is sweep 400; failure of the gate at sweep 400 is retained as a
slow case rather than relabelled as convergence. Objective eligibility remains
separate from practical convergence.

## Truth isolation and selection

Data generation creates observation-only bundles and separately sealed truth
bundles. Truth is unavailable to fitting, stopping, analysis, and endpoint
selection. For each data set, the winner is selected from the two eligible
fixed-400 endpoints by maximum finite ordinary `T=1` ELBO. The six winner IDs
are frozen before any future truth authorization.

## Registered truth-free outputs

- completion, objective eligibility, convergence reason, slow-case status,
  warnings, elapsed time, and peak memory;
- the 380-to-400 change in fitted values, RSS, variable/factor PPI, loading
  subspaces, feature subspaces, score means, score posterior kernels, and full
  contributions;
- truth-free ELBO winner and winner-runner-up ELBO gap for each data set.

The previously registered development safety boundaries are retained:
fitted 0.003, RSS 0.006, loading 0.001, feature 0.02, score mean 0.02, score
kernel 0.05, and full contribution 0.01. These boundaries describe local
endpoint stability only; they do not establish truth recovery.

## Completion and decision rules

Execution integrity requires 12/12 complete, objective-eligible fixed-400
endpoints, exact public-profile provenance, six truth-free ELBO selections, and
zero truth/continuation/automatic-800 use. Any error is retained without a
selective automatic rerun.

If some endpoints are slow or fail a 380-to-400 stability boundary, Stage 4C
reports them without automatically extending to 800. A separate user decision
is required before truth unsealing, continuation, or a protocol change.

Expected ECS wall time is approximately 6.5-7.5 hours with four outer workers;
expected additional disk use is approximately 1-1.5 GiB.
