# G12 truth-free dimension-selection Stage 6B draft

Protocol ID: `G12_TRUTH_FREE_DIMENSION_SELECTION_STAGE6B_V1_DRAFT_20260901`

Draft date: 2026-09-01

Status: development draft; Stage 6B-A read-only audit completed and Stage 6B-B
experiment-local holdout interface implemented/tested; no predictive pilot is
registered, frozen, authorized, or executed

## 1. Scope and evidence boundary

This draft defines a truth-free route for selecting fitted factor counts
`(L_f,L_s[1],L_s[2])` and within-factor FPCA truncation counts `M` in
multiFSYNC. It does not modify the model, priors, CAVI updates, ordinary ELBO,
G12 initialization, Jaoua annealing, or the fixed 400-sweep development budget.

The completed Stage 6A result for data seed `83131001` is development evidence
whose truth has already been unsealed. It may be used to identify failure modes
and design this rule, but it cannot validate the rule. In particular, Stage 6A
fits used all observed time points, so their fitted-point errors cannot be
retrospectively relabelled as out-of-sample prediction errors.

No full-size predictive fit, continuation, data generation, or truth access is
authorized by this draft. The user separately authorized the Stage 6B-A
read-only audit and a nonformal Stage 6B-B software micro smoke.

## 2. Jaoua method and what can be inherited

The main method in Jaoua, Temko and Ruffieux (2026), *A scalable Bayesian
functional factor model for high-dimensional longitudinal molecular data*, is
an overfit-and-prune strategy:

In notation, Jaoua's factor count `Q` corresponds to the role-separated
multiFSYNC factor counts `L_f` and `L_s[s]`; Jaoua's within-factor FPCA count
`L^(q)` corresponds to multiFSYNC's `M_f[l]` and `M_s[[s]][l]`. It does not
correspond to multiFSYNC's factor-count symbol `L`.

1. fit a conservative factor upper bound `Qmax`;
2. compute loading inclusion probabilities `PPI_jq` and define

   `PPI_factor_q = 1 - product_j(1 - PPI_jq)`;

3. retain factors whose factor PPI exceeds a pre-specified threshold, with 0.5
   given as a median-probability-model example;
4. fit a conservative within-factor component upper bound `Lmax`, rotate the
   fitted factor trajectories into an ordered FPCA representation, and retain
   components until cumulative PVE plateaus or exceeds a high threshold such as
   99%.

Their simulations generally fit `Qmax=5` and `Lmax=5` and use annealed
variational inference. The paper reports that annealing is important because
unannealed overfitted factors can be retained and partly compensate for the mean
structure. The reported factor-count averages in a true-`Q=2` comparison are
1.52, 1.68, 1.96, 1.88 and 1.88 for sample sizes 20, 25, 30, 35 and 40,
respectively; thus the paper's rule is useful but not uniformly exact at small
sample sizes.

The separate `bayesSYNC_model_choice()` function in the local reference source
fits candidate dimensions independently and maximizes ordinary ELBO plus a
truncated-Poisson log prior on dimension. This is not the paper's principal
overfit-and-prune rule and is not treated here as independently validated
evidence for multiFSYNC.

The transferable lessons are:

- deliberate upper-bound fitting can be computationally simpler than a fully
  discrete Bayesian rank model;
- annealing and component-specific shrinkage can help remove redundant
  directions;
- factor inclusion, component variance, and prediction are different pieces of
  evidence and should not be conflated;
- post-fit PVE is a descriptive effective rank, not a posterior probability for
  the discrete rank.

The parts that cannot be copied mechanically are:

- multiFSYNC must distinguish shared and study-specific roles and detect
  duplicates or misplaced directions; a single-study factor PPI does neither;
- at `p=500`, `1-product(1-PPI_jq)` can be close to one even when many individual
  inclusion probabilities are small;
- in the Stage 6A `M=4` arm, a rigid 99% PVE rule retained all four components
  in 29 of the 36 fitted factors and never retained only two;
- Stage 6A showed that high-PPI extra or duplicate factors can have non-negligible
  complete-contribution energy.

Primary paper: <https://arxiv.org/abs/2603.20844>

## 3. Stage 6A implications

Within each Stage 6A fit configuration, the winner was selected without truth by
maximum objective-eligible ordinary ELBO. The following results motivate this
draft:

| fitted configuration | winner ELBO | observed signal NRMSE | dense signal NRMSE | feature ISE | interpretation |
|---|---:|---:|---:|---:|---|
| `truth_L1_M2` | -121854.62 | 0.170 | 0.230 | 0.094 | correct identity and good functional recovery |
| `factor_L3_M2` | -121231.87 | 0.143 | 0.894 | 0.580 | three correct directions plus duplicates/extras; poor interpolation |
| `fpca_L1_M4` | -120409.93 | 0.145 | 0.769 | 0.251 | correct factor identity but over-retained FPCA variation |

The cross-configuration ELBO values above are shown only as a failure diagnostic;
they were not eligible for Stage 6A selection. The increasingly flexible models
improve the fitted-point criterion while substantially degrading dense-grid
recovery. Within the `factor_L3_M2` arm, higher ELBO was also associated with
worse feature ISE and dense reconstruction, and the maximum-ELBO start had the
worst dense and functional recovery in that arm.

Consequently:

- training ordinary ELBO is retained as a numerical/objective diagnostic and a
  within-fixed-configuration basin-ranking quantity;
- ordinary ELBO is not used to compare fitted dimensions;
- observed-point fitted error is not used as a substitute for prediction;
- factor PPI, contribution energy, or PVE alone cannot determine dimension.

## 4. Selection architecture

The proposed rule is deliberately lexicographic rather than a weighted sum.

### Gate 1: objective and execution eligibility

A candidate fit must have a finite ordinary `T=1` ELBO, a complete endpoint,
valid dimensions, no fatal numerical error, and the registered truth-isolation
flags. Practical convergence remains a reported diagnostic rather than an
automatic scientific-quality gate.

### Gate 2: genuine held-out prediction

Dimensions are compared by prediction of observations that were excluded before
fitting. The primary loss is mean held-out negative log predictive density
whenever the posterior predictive variance implementation has passed a dedicated
calibration test. Held-out squared prediction error and absolute prediction
error are mandatory secondary quantities.

No latent signal, dense truth grid, loading identity, feature ISE, `R/P/L`, or
other truth-derived quantity may enter this step.

### Gate 3: one-standard-error simplicity rule

Let `Loss(c)` be the paired held-out loss for dimension configuration `c`. Find
the configuration with the smallest mean loss. Select the simplest nested
configuration whose mean loss is within one standard error of that minimum. The
standard error is calculated over independent subject-level prediction units,
stratified by study; in a later multi-seed confirmation, data seed is the primary
replication unit.

Complexity is ordered first by total number of fitted factors and then by total
number of within-factor FPCA components. Non-nested configurations are reported
as a partial order and are not forced into an arbitrary scalar complexity
penalty.

### Audit 4: contribution activity and cross-start stability

The selected dimension must be accompanied by an audit, not a new weighted
score. For each candidate factor report:

- factor PPI and expected active-loading count;
- loading energy and factor-process energy;
- complete-contribution energy and its share of total fitted contribution;
- cross-start recurrence of its loading direction within the same shared or
  study-specific role;
- within-role duplicate directions and cross-role leakage.

Complete-contribution energy is preferred to raw loading or process energy for
activity because it is invariant to reciprocal rescaling of loading and factor
process. Energy alone is not a selector: Stage 6A spurious factors had a median
complete-contribution RMS of about 0.086 and reached about 0.215.

Cross-start recurrence is computed after sign-invariant, one-to-one matching of
loading directions within the same role. Results are reported for all eligible
starts and for a pre-specified top-ELBO subset. Stability is not allowed to
overrule prediction automatically, because earlier development experiments show
that the scientifically correct basin can be rare among random starts.

For FPCA components, report ordered eigenvalue share, cumulative PVE,
cross-start subspace stability, and incremental held-out predictive gain. The
99%, 95%, and 90% PVE cutoffs may be shown as sensitivity summaries, but no
cutoff may be selected after examining truth.

## 5. Primary hold-out design

The first implementation should use whole-time-point holdout within subject:

- hold out exactly one registered interior observation time from every subject;
- remove the complete `p`-variable vector at that time before fitting, rather
  than holding out random cells;
- keep at least five observed times in the training data for every subject;
- use identical time-point splits and paired fit seeds for every dimension
  configuration;
- average variable-level loss within each subject-time vector, then average over
  subjects, so the `p=500` coordinates are not treated as 500 independent
  replicates;
- preserve study-stratified subject identifiers in the score table.

The package and earlier development runners contained fitted reconstruction and
posterior-variance fields, but no dedicated held-out prediction interface. The
Stage 6B-B implementation now provides an experiment-local whole-time splitter,
posterior-mean predictor, and paired point-loss scorer. It does not modify the
package API. A fully calibrated posterior predictive variance is still absent;
therefore its residual-noise-only NLPD is explicitly diagnostic and is not
eligible as the primary dimension-selection score.

An endpoint-time holdout is a separate extrapolation stress test and must not be
mixed with the primary interpolation score.

Whole-subject validation answers a different question. A completely unseen
subject has no estimated subject score, so a useful secondary design must either
score the marginal distribution of the whole subject or give the subject a
registered set of support times, infer only its scores under frozen global
parameters, and predict separate target times. This requires a verified
conditional-score interface and is deferred until the simpler time-holdout route
works.

## 6. How ordinary ELBO is used

For the first predictive pilot, ordinary ELBO has only three roles:

1. determine objective eligibility;
2. rank starts fitted to exactly the same training observations and exactly the
   same dimension configuration;
3. break a predictive tie only after the one-standard-error and simplicity rules.

It is not normalized per observation and then reused across dimensions; such a
normalization does not create a complexity penalty. It is not augmented with an
ad hoc BIC or Poisson term in the primary analysis. Stage 6A ELBO gains from
larger configurations were hundreds to more than one thousand units, whereas a
unit-rate Poisson dimension penalty would be only a few units for these changes.
More importantly, dimension-dependent variational gaps make a fixed penalty hard
to interpret.

If the held-out pilot shows that the maximum-ELBO start is repeatedly predictively
inferior within a fixed dimension, a separately registered comparison may test a
truth-free start selector: shortlist the top ELBO starts, choose by held-out
prediction, and use contribution/stability only as diagnostics. It must be
validated before replacing the current maximum-ELBO G12 winner rule.

## 7. What can be done with existing Stage 6A results

A read-only offline audit can be run without fitting to:

- cluster loading directions across the 12 starts of each configuration;
- calculate factor recurrence, duplicate/leakage rates, complete-contribution
  shares, and FPCA subspace stability;
- test the reproducibility of those summaries under pre-specified sensitivity
  grids for contribution share and cosine matching;
- verify that none of the calculations reads sealed/unsealed truth fields.

This audit can assess whether energy and stability are measurable and whether
they diagnose the known Stage 6A failure. It cannot estimate held-out predictive
loss because every Stage 6A observation was used for fitting. Naive residuals,
WAIC, or pointwise LOO computed from those conditioned fits are exploratory only,
especially because an entire `p=500` subject-time vector shares latent scores.

## 8. Minimal future predictive pilot

Before any broad dimension-selection experiment, use one development data set
and a small candidate grid:

- factor screen at fixed `M=2`: `(1,1,1)`, `(2,2,2)`, `(3,3,3)`;
- FPCA screen at the selected factor count: `M=2`, `M=3`, `M=4`;
- one registered interior-time holdout split;
- paired G12 start seeds, initially four starts per configuration for software
  and signal detection;
- no truth access until predictions, selected dimension, and audit summaries are
  frozen.

This is a development pilot, not a universal dimension-selection validation.
If it discriminates the Stage 6A-type overfit without truth, the split count,
start count, thresholds, and confirmation seeds must be frozen in a later
protocol before independent evaluation.

## 9. Decision after the pilot

The route advances only if held-out prediction prefers, or selects by the
one-standard-error rule, a simpler configuration while the energy/stability
audit explains why larger configurations are redundant. Failure modes are
separated:

- prediction favors the larger model: current evidence does not justify pruning;
- prediction is tied but stability is poor: select the simpler model and report
  structural uncertainty;
- prediction and stability both favor the larger model: retain it even if a
  development truth later shows over-selection, then revisit the model/prior in
  a separate stage;
- factor PPI alone favors the larger model: no decision, because this is already
  a known non-discriminating diagnostic at high dimension.

No threshold or success claim in this draft is frozen until the user separately
reviews and authorizes a registered pilot.

## 10. Implementation and execution record

On 2026-09-01, Stage 6B-A read exactly 36 existing Stage 6A fit objects, their
36 terminal records, and the three-row truth-free winner table. It did not read
the observation bundle, sealed/unsealed truth, or prior evaluation output, and
did not start a fit or continuation. The versioned result is stored under
`stage6ba_readonly_audit_20260901_v1/`.

The first audit attempt stopped before writing results because an integrity
assertion incorrectly required a zero loading direction to have self-cosine
one. The project contract defines every cosine involving a zero direction as
zero. Only that integrity assertion was corrected to require finite similarities
in `[0,1]`; the matching convention and all scientific summaries were left
unchanged. The clean rerun passed all integrity checks.

Stage 6B-B added only experiment-local functions for:

- removing one complete, interior `p`-variable observation vector per subject;
- preserving at least five training times;
- predicting the held-out time from fitted mean, shared/specific loadings,
  interpolated eigenfunctions, and subject FPCA scores;
- reporting subject-time-averaged RMSE, NRMSE and MAE for all factors and a
  factor-PPI `>=0.5` sensitivity mode;
- labelling residual-noise-only NLPD as uncalibrated and ineligible as the
  primary dimension score.

Targeted tests passed `9/9`, including a split-only check on the full Stage 6A
observation bundle (`60` subject-time vectors, `p=500`, no truth access). A
nonformal `S=2`, `n_s=(4,4)`, `p=20`, six-sweep micro fit passed the full
split-fit-predict-score chain on both Windows and the target Linux server, with
identical score values. No Stage 6B full-size predictive fit was started.
