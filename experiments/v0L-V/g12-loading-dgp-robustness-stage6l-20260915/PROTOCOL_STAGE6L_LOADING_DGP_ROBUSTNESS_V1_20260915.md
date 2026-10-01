# Stage 6L two-loading-scenario audit and micro bridge protocol V1

Date: 2026-09-15

Status: registered development sensitivity study. This is not a formal v0L-V
result, not a paper-level Monte Carlo experiment, and not an authorization to
change the model, prior, CAVI, ELBO, G12 route or gate schedule.

## Scientific question and neutral interpretation

The primary question is whether the established multiFSYNC conclusions are
fragile to replacing the controlled strong-identifiability loading DGP with a
Jaoua-style random-sparsity loading DGP. The study does not assume that the new
DGP should improve performance. Better, similar and worse recovery are all
reportable evidence.

The secondary purpose is to provide a reproducible explanation to the
supervisor for retaining the controlled loading DGP as the main structure-
recovery scenario while optionally adding the Beta-random DGP as a robustness
scenario. A one-dataset micro bridge cannot establish broad performance.

## Frozen loading scenarios

### `controlled_identifiable_fixed90_calibrated`

- `identified_loadings=TRUE`;
- `bool_sparse_loadings=TRUE`, `prop_sparse=0.9`;
- shared and study-specific supports are disjoint within each study;
- each L1M2 block has exactly 50 active variables at `p=500`;
- loading columns are initially generated at the fixed norms supplied by
  `generate_identified_loadings`, then receive the already established v0F
  `constant_per_active_loading` calibration with `anchor_p=100`,
  `active_fraction=0.1`, and `anchor_norm=3`;
- with 50 active variables this gives every shared/specific loading column
  target norm `3*sqrt(50/10)`, rather than retaining the generator's raw
  1.6/0.8 norm difference;
- this is the controlled, strong-identifiability main DGP.

### `jaoua_beta10_blockwise_calibrated`

- `identified_loadings=FALSE` and `bool_sparse_loadings=TRUE`;
- independently for every shared or study-specific loading block,
  `omega_l ~ Beta(1,10)`;
- conditional on `omega_l`, each variable support indicator is independently
  `Bernoulli(omega_l)`;
- a zero-support block is forced to contain one randomly chosen active variable;
- all raw loading values are independently `N(0,1)` and inactive values are
  zero;
- the same v0F calibration is then applied separately to each nonempty loading
  column, with target norm `3*sqrt(active_count/10)`;
- shared and specific supports may overlap and loading norms are not forced to
  be ordered or orthogonal.

This is the model-intended loading law in Jaoua et al. adapted blockwise to a
multi-study shared/specific model. The bound upstream simulator is
`https://github.com/hruffieux/bayesSYNC`, file `R/simulate_data.R`, commit
`816db2fe27bb6f596012135caaea3c8852eb8da8`. That implementation calls
`rbinom(p*Q, prob=omega)` before reshaping to a `p x Q` matrix. R recycles the
probability vector elementwise, so for `Q>1` it does not necessarily apply one
fixed `omega_q` to an entire matrix column. Stage 6L follows the paper/model
meaning `gamma_jl | omega_l ~ Bernoulli(omega_l)` rather than duplicating that
likely vector-recycling artefact. Consequently the arm is called
"Jaoua-style/model-intended", not an exact seed-level reproduction.

The common calibration holds the intended per-active loading information
constant, but the scenarios still differ in active counts, support overlap,
direction geometry and total block information. Therefore any performance
difference is a loading-DGP sensitivity result; it is not isolated causal
evidence that the Beta distribution alone helps or hurts.

## Common data-generating configuration

Both arms use `S=2`, `n_s=(30,30)`, `p=500`, `d=0`, true `L_f=1`, true
`L_s=(1,1)`, `M_f=2`, `M_s=(2,2)`, `K=5`, irregular 6--9 observation times,
`sigma_eps=0.3`, `mean_amp=0.6`, score eigenvalue weights proportional to
`1/m^2`, `n_dense=201`, and `truth_basis_method="jaoua_splinets"`.

An optional `truth_loading_seed` isolates loading generation and restores the
main RNG stream. Paired arms therefore share times, fitted O'Sullivan design,
truth functions, scores, mean phases, covariates and measurement errors; only
loadings, loading-derived signals and responses differ. Historical calls with
`truth_loading_seed=NULL` retain their prior RNG sequence and output.

## Part A: generation-only audit

The audit uses 200 registered paired generator seeds, indices 1--200:

- `data_seed = 96150000 + index`;
- `truth_basis_seed = 97150000 + index`;
- `truth_loading_seed = 98150000 + index`.

No model is fitted. The audit records blockwise active counts, sampled omega,
calibrated loading norms, support overlaps, normalized loading cosines, studywise loading
Gram eigenvalues/condition numbers, signal and noise RMS, factor-only signal
RMS, empirical SNR and exact signal/noise decomposition. It also verifies exact
pairing of all non-loading random quantities.

Mechanical validity requires 200/200 pairs, no generator error, nonempty
loading blocks, finite metrics, exact non-loading pairing, response
decomposition error below `1e-12`, and target loading norm per square-root
active count equal to `3/sqrt(10)` within `1e-12`. No recovery threshold is
tested in Part A.

Before registering the common calibration, a 200-pair raw-scale preflight was
run and retained as `generation_audit_v1_20260915`. It showed that the raw
controlled generator (norms 1.6/0.8) and unscaled Jaoua-style standard-normal
loadings differed by several-fold in factor-signal RMS. That preflight is not
the fitting DGP and is retained to document why the calibrated V2 audit is
necessary rather than silently changing the comparison.

## Part B: extremely small paired fitting bridge

The bridge uses one new, separately registered logical dataset:

- `data_seed=96151001`;
- `truth_basis_seed=97151001`;
- `truth_loading_seed=98151001`;
- both loading scenarios above;
- four paired fit seeds `99151001:99151004`;
- 2 arms x 4 starts = 8 fits.

Every fit uses the already selected development route: random initialization,
Gram unit-function-energy initialization (`G12`), one pre-score sweep, Jaoua
annealing `c(1,1.9,100)`, exactly 400 ordinary `T=1` sweeps and `n_cpus=1`.
Independent fits may run in parallel. There is no racing, continuation, 800-
sweep extension, gate modification or alternative initializer.

Within each loading arm, the winner is the objective-eligible endpoint with
maximum finite ordinary `T=1` ELBO; lexical `fit_id` is the registered tie
breaker. ELBO is never compared across arms. Truth is unavailable to fitting,
stopping and selection. Only after both within-arm winners are frozen may the
sealed truth be read for the already established loading-first scientific
evaluation.

The bridge reports, without turning them into a paper-level success claim:

- how many of four starts reach the correct shared/specific core in each arm;
- whether each truth-free winner recovers that core;
- reconstruction, feature/kernel, loading, score, process/contribution and
  structural error metrics;
- ELBO, runtime, warnings and convergence descriptors;
- the loading-DGP scale/SNR context from Part A.

One paired dataset can falsify the claim that the Beta arm is mechanically
unusable, but it cannot establish general superiority, equivalence or
inferiority. An adverse result is not followed by gate tuning. The next main
experiment remains the supervisor-requested covariate-effect simulation.
