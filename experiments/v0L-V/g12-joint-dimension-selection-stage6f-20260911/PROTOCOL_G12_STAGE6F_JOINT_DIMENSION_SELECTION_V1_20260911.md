# G12 Stage 6F joint dimension-selection protocol v1

Date: 2026-09-11  
Status: development protocol; not a formal paper Monte Carlo experiment

## Question

Can a truth-free paired holdout criterion distinguish the correctly specified
`L1M2` fit from factor-count overspecification (`L2M2`), FPCA-cap
overspecification (`L1M4`), and joint overspecification (`L2M4`) after start
selection has been performed within every fixed dimension configuration?

The experiment separates two failure modes:

1. optimization reachability within a fixed `(L,M)` configuration, handled by
   four paired starts and ordinary-ELBO start selection;
2. selection across different `(L,M)` configurations, handled only by held-out
   predictive scores.

Ordinary ELBO values are never compared across dimension configurations.

## Frozen inputs and fitted grid

- data: existing Stage 6A development data `g12s6a_01`;
- generated truth: shared factor count 1, each study-specific factor count 1,
  and two FPCA components per factor;
- observed design: two studies, 30 subjects per study, 500 variables, no
  covariates, 6--9 irregular times, strong-identifiability setting;
- configurations: `L1M2`, `L1M4`, `L2M2`, `L2M4`;
- `L` is applied to the shared and both study-specific factor blocks;
- `M` is applied to every factor within every block;
- paired start seeds: 83206101--83206104.

The truth dimensions above document the later diagnostic target. Neither the
truth bundle nor any truth-derived metric is available to fitting, start
selection, holdout scoring, or dimension selection.

## Holdout splits

Every subject contributes one complete 500-variable vector to each split.

1. `split1_middle`: the existing Stage 6B interior middle-time split,
   `floor((n_i + 1)/2)`.
2. `split2_upper`: the next interior time,
   `floor((n_i + 1)/2) + 1`.

For 6--9 observed times these indices are distinct and interior, and each fit
retains at least five training times. Splits are deterministic and identical
across dimension configurations and paired start seeds.

## Endpoint accounting and reuse

The analysis has 2 splits x 4 configurations x 4 starts = 32 endpoints.

- Reused from the corrected Stage 6B v2 output for `split1_middle`:
  `L1M2`, `L2M2`, and `L1M4`, four starts each (12 endpoints).
- Newly fit for `split1_middle`: `L2M4`, four starts.
- Newly fit for `split2_upper`: all four configurations, four starts each
  (16 endpoints).

Exactly 20 new fits are permitted. Reused fits must have the same observation
bundle, split, start seed, route, and fixed-400 endpoint. Their fit objects are
hard-linked when possible and otherwise copied; their provenance and SHA256 are
recorded. No existing Stage 6B artifact is modified.

## Fit route

All new fits use the current G12 development route without modification:

- random initialization;
- Gram-unit-energy function initialization;
- one pre-score sweep;
- Jaoua default annealing `c(1, 1.9, 100)`;
- 400 ordinary T=1 sweeps after annealing;
- objective eligibility required;
- `n_cpus=1` inside every fit;
- no continuation, racing, truth access, or full-data refit.

Independent fits may run concurrently at the task layer. On the 8-vCPU,
16-GiB ECS, seven fit workers are the registered launch maximum.

## Selection order

1. For each split and fixed `(L,M)` configuration, select one of four starts by
   maximum eligible ordinary T=1 ELBO, with lexical `fit_id` tie-breaking.
2. Freeze and hash the eight start winners before reading any held-out response
   for scoring.
3. Score only those eight winners.
4. Aggregate each criterion first within subject across the two splits, then
   across the 60 independent subjects with study-stratified standard errors.
5. Apply the one-standard-error rule separately to mean squared error and the
   experimental joint negative log predictive density.

Complexity is ordered by the registered proxy

`C(L,M) = 1500 L + 141 L M`,

giving `L1M2=1782`, `L1M4=2064`, `L2M2=3564`, and `L2M4=4128`. The proxy is
used only to choose the simplest candidate inside the one-SE set; it is not a
likelihood penalty or a new fitted-model objective.

## Scores and diagnostics

The established primary interface reports all-factor held-out MSE, RMSE, and
NRMSE. The experimental joint score is a moment-matched Gaussian negative log
predictive density computed from posterior mean and covariance moments via a
diagonal-plus-low-rank representation. It is not claimed to be the exact
posterior predictive density.

Also record:

- marginal 95% interval coverage;
- mean standardized squared residual;
- marginal NLPD and mean predictive SD;
- variance floors and numerical clamps;
- per-split risks and ranks;
- Spearman rank agreement, minimum-risk winner agreement, and pairwise ordering
  concordance across the two splits.

## Truth boundary and assessment

The truth-free start-winner and dimension-selection tables must exist and their
SHA256 values must be bound in a frozen marker before the sealed truth bundle
is opened. After that boundary, evaluate the eight configuration winners with
the existing frozen evaluator.

The criterion is considered promising for broader validation only if all of
the following hold in this pilot:

- scoring has no error, non-finite value, material negative variance, or
  unexplained floor/clamp event;
- the two splits show directionally coherent configuration rankings;
- one-SE does not prefer a jointly over-specified model when a simpler model is
  statistically indistinguishable;
- the selected endpoint has scientific recovery at least comparable to the
  best of the four evaluated dimension configurations.

Failure of any item is reported, not selectively rerun. Passing this pilot does
not establish generality because it uses one data seed and only two correlated
within-subject splits.

## Prohibitions

No new data generation, model/prior/CAVI/ELBO change, truth-based start or
dimension selection, 400-to-800 continuation, full-data refit, formal v0L-V
mutation, or paper-level success claim is authorized.
