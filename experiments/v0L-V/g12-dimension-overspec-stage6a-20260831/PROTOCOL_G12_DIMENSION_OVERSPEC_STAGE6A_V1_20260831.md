# G12 dimension-overspecification Stage 6A protocol

Protocol ID: `G12_DIMENSION_OVERSPEC_STAGE6A_V1_20260831`

Registration date: 2026-08-31

Status: registered development experiment; not a formal v0L-V or paper Monte Carlo result

## Question and claim boundary

Stage 6A asks whether the established G12 route remains scientifically useful when
the fitted factor count or the fitted within-factor FPCA rank cap is deliberately
larger than its data-generating truth. It does not reselect the G12 optimization
route, modify the model or prior, or estimate a paper-level success probability.

One new data seed is used for a time-bounded screening experiment. It can expose
software failures, gross scientific degradation and dimension-specific failure
modes, but it cannot establish cross-data generality or a dimension-selection
success probability. Two additional candidate seeds remain unused and unseen;
they may be registered later only in a separate independent-confirmation stage.

## Data registration

The single registered data ID and seed are:

| data ID | seed | scenario |
|---|---:|---|
| `g12s6a_01` | 83131001 | baseline strong |

It uses the established standard-density generator: `S=2`,
`n_s=(30,30)`, `p=500`, `d=0`, true shared factor count 1, true
study-specific factor count 1 in each study, true `M=2` for every factor,
`K=5`, 6--9 irregular observations per subject, `sigma_eps=0.3`,
`mean_amp=0.6`, fixed 90% loading sparsity, and strong loading separation.
The data are marked unseen at registration. Observation-only bundles are emitted
separately from sealed truth bundles.

## Separate fitted-dimension arms

The same data set is fitted under three separate arms:

| fit config | fitted factor counts `(L_f,L_s1,L_s2)` | fitted FPCA caps | purpose |
|---|---|---|---|
| `truth_L1_M2` | `(1,1,1)` | `2` per factor | reference |
| `factor_L3_M2` | `(3,3,3)` | `2` per factor | factor-count overspecification only |
| `fpca_L1_M4` | `(1,1,1)` | `4` per factor | FPCA-rank overspecification only |

There is no joint `L=3, M=4` arm. This prevents a failure from being
unidentifiably attributed to two simultaneous changes and limits compute cost.
Ordinary ELBO values are comparable only among the 12 starts of the same data ID
and fit configuration. Cross-configuration ELBO comparison or selection is
forbidden.

## Optimization route and fit count

Every arm uses the unchanged development G12 route:

- random initialization with Gram unit-function-energy calibration;
- one pre-score sweep;
- Jaoua annealing `c(1,1.9,100)` (99 annealing sweeps);
- 12 independently seeded starts per data/configuration stratum;
- fixed 400 ordinary `T=1` sweeps, with checkpoints at 380 and 400;
- `lambda_orth=0`, `n_cpus=1`, no racing, no continuation, and no automatic 800;
- truth-free winner selection by maximum objective-eligible ordinary `T=1`
  ELBO within each data/configuration stratum, with lexicographic `fit_id` tie
  breaking.

This yields `1 data x 3 configurations x 12 starts = 36 fits` and three
truth-free winners. Fit seeds are bound in `FIT_MANIFEST.csv`; they are not
generated at run time.

## Capacity gate

The first start in each of the three arms is the registered three-fit capacity
probe. The probes run concurrently on the target Linux host, are part of the 36
fits, and are never selectively rerun. Full execution starts
only if all three processes exit successfully, all endpoints are objective
eligible, GNU time records RSS for all three, the conservative sum of their peak
RSS is below 6 GiB, and at least 10 GiB disk remains. The observed probe times
replace the preliminary 1--2 day wall-time estimate.

## Truth isolation and stopping point

Generation, fitting, stopping, capacity checks, and selection cannot read truth.
After all 36 fits and three selections are frozen, the runner writes
`WAITING_FOR_TRUTH_UNSEAL_AUTHORIZATION.txt` and stops. Evaluation requires a
separate authorization file bound to the SHA-256 of the frozen winner table.
Stage 6A does not authorize new fits or continuation after truth unsealing.

## Registered evaluation profile

After separate authorization, evaluation is loading-first and uses one-to-one
truth/candidate matching. It reports both every start and the three preselected
winners. The evidence profile is:

1. correct, missing, misplaced, duplicate and extra loading directions at global
   and study-role levels;
2. candidate factor PPI, expected active-loading count, loading energy, factor
   process energy and complete-contribution energy;
3. fitted rank cap, post-fit 99% PVE `effective_M`, retained eigenvalue shares,
   feature ISE, projection floor, covariance-kernel error and score recovery;
4. observed-point and dense-grid signal reconstruction;
5. loading, factor-process and complete-contribution recovery;
6. objective eligibility, within-stratum ELBO support, cross-start identity
   stability, warnings, runtime, memory and the 380-to-400 checkpoint behavior.

Factor PPI is auxiliary because its union probability can saturate at `p=500`.
The `effective_M` output is a post-fit 99% PVE truncation of a model fitted at the
registered rank cap; it is not a Bayesian discrete dimension-selection result.
No truth metric enters fitting or winner selection.

## Interpretation rule

The single data set supports a screening conclusion only. Evidence of tolerance
requires agreement across loading identity, candidate activity, reconstruction,
and functional/covariance recovery; a favorable reconstruction alone is not
enough. Conversely, an extra factor with high union PPI but negligible loading,
process and contribution energy is reported as a PPI-selection limitation rather
than automatically as a scientifically active factor. A favorable Stage 6A result
must be independently confirmed on unused data before any general claim.

## Prohibited actions

Stage 6A does not modify the model, priors, CAVI, ELBO, G12 initialization,
annealing, multistart count, or 400-sweep budget. It does not compare against
bayesSYNC, fit covariate effects, use racing, extend to 800, or claim formal v0L-V
or paper Monte Carlo validation.
