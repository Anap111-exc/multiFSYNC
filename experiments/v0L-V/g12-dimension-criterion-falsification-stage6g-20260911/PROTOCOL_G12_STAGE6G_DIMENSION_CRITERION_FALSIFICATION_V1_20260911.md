# G12 Stage 6G dimension-criterion falsification protocol v1

Protocol ID: `G12_STAGE6G_DIMENSION_CRITERION_FALSIFICATION_V1_20260911`  
Date: 2026-09-11  
Status: local draft; no data generation or fit is authorized

## 1. Decision question

Can one fixed, truth-free dimension-selection rule choose both a minimal and a
non-minimal generating dimension, after optimization reachability has first
been handled within each fixed configuration?

The rule under test is the repeated paired whole-time holdout MSE/NRMSE plus
the one-standard-error simplicity rule already used in Stage 6F. This screen is
designed to falsify the rule cheaply. It is not a paper-level Monte Carlo and
does not attempt to establish dimension-selection consistency.

Stages 6D and 6F motivate the design:

- Stage 6D used three data seeds whose generating dimension was always the
  smallest candidate. The rule selected both axes correctly for two of three
  data sets and was sensitive to one high-influence subject-time vector in the
  third.
- Stage 6F used two paired holdout times on one existing data set. MSE/NRMSE
  selected the true smallest joint dimension and agreed with the best overall
  scientific recovery; the experimental joint log-score selected the largest
  model and failed.

Consequently, repeating only a true `L111/M2` setting cannot answer whether
one-SE learns dimension or merely favours the smallest candidate. Stage 6G
includes two non-minimal truths and drops joint log-score from the selection
rule.

## 2. Frozen scientific components

The following are unchanged:

- multiFSYNC model, priors, CAVI updates and ordinary ELBO;
- G12 random Gram-unit-function-energy initialization;
- one pre-score sweep and Jaoua annealing `c(1,1.9,100)`;
- `lambda_orth=0`;
- fixed endpoint of 400 ordinary T=1 sweeps;
- objective eligibility separated from practical convergence;
- maximum eligible ordinary ELBO for start selection within a fixed
  `data x split x dimension` stratum;
- loading-first truth evaluation after the truth-free selection is frozen.

No continuation, racing, alternative initialization, model change, prior
change or new scientific success metric is permitted.

## 3. Data-generating scenarios

All three scenarios use `S=2`, `n_s=(30,30)`, `p=500`, `d=0`, `K=5`, strong
loading identifiability, 10% active loadings, 6--9 irregular times,
`sigma_eps=0.3`, `mean_amp=0.6`, and score variances proportional to `1/m^2`.
The loading sparsity mechanism remains the existing fixed-support calibrated
generator; it is not silently replaced by a beta-generated sparsity scheme.

Truth functions use the validated `truth_basis_method="jaoua_splinets"` and
the current phase-shifted periodic mean. This protocol therefore changes
neither loading sparsity nor mean generation. If the supervisor later requires
a spline-generated mean, this draft must be revised before execution; an
executed v1 must never be mutated afterward.

| Scenario | Data seed | Truth-basis seed | Generating dimensions | Purpose |
|---|---:|---:|---|---|
| `s6g_minimal_L111_M2` | 84611001 | 84611901 | `L_f=1, L_s=(1,1), M=2` | detect over-selection |
| `s6g_nonminimal_shared_L211_M2` | 84611002 | 84611902 | `L_f=2, L_s=(1,1), M=2` | falsify automatic preference for the smallest factor tuple |
| `s6g_nonminimal_rank_L111_M3` | 84611003 | 84611903 | `L_f=1, L_s=(1,1), M=3` | falsify automatic preference for `M=2` |

The asymmetric `L211` truth is deliberate. The bound upstream simulator
hard-codes four factor-specific degree/periodicity/knot templates and therefore
supports at most four total truth factors; this is an implementation boundary
of that exact simulation template, not a mathematical limitation of B-spline
orthogonalisation or of the Splinets package. A uniform
`L_f=L_s1=L_s2=2` truth would have six factors and is outside the currently
validated exact-reproduction domain. `L211` has four total factors and also
tests whether prediction distinguishes a second shared factor from an extra
study-1-specific factor. Any future `Q>4` generator must be separately named,
versioned and validated rather than silently recycling these four templates.

The three data seeds and three truth-basis seeds are new relative to the
current experiment manifests. Each observation bundle and sealed truth bundle
must be written separately. Fitting and scoring receive only the observation
bundle.

## 4. Common candidate grid

Every scenario uses the same five candidates; the candidate set is not tailored
to the known generating truth.

| Candidate | `(L_f,L_s1,L_s2)` | Common M | Complexity proxy |
|---|---|---:|---:|
| `cand_L111_M2` | (1,1,1) | 2 | 1782 |
| `cand_L111_M3` | (1,1,1) | 3 | 1923 |
| `cand_L111_M4` | (1,1,1) | 4 | 2064 |
| `cand_L121_M2` | (1,2,1) | 2 | 2356 |
| `cand_L211_M2` | (2,1,1) | 2 | 2416 |

For `p=500`, `n_1=n_2=30`, and `K=5`, the registered proxy is

`C = p Q + M{(n_1+n_2)L_f + n_1 L_s1 + n_2 L_s2 + (K+2)Q}`,

where `Q=L_f+L_s1+L_s2`. It counts the main loading, subject-score and
function-coefficient dimensions sufficiently to order these candidates. It is
used only to choose the simplest member of a one-SE set; it is not inserted
into the likelihood, prior or ELBO.

The grid is intentionally compact. It includes the three generating models,
an `M=4` overfit alternative, and an extra study-specific factor that competes
with the second shared factor. It does not claim to validate exhaustive search
over every asymmetric `L_s` tuple.

## 5. Holdout splits, starts and endpoint count

Each subject supplies one complete 500-variable vector to each of two separate
fits:

1. `split1_middle`: `floor((n_i+1)/2)`;
2. `split2_upper`: `floor((n_i+1)/2)+1`.

Each split leaves at least five training times. The two splits use adjacent
times from the same subjects and are therefore paired repeated measurements,
not independent data replicates.

Four fitting seeds `84611101--84611104` are paired across all scenarios,
splits and candidates. Each fit uses `n_cpus=1`; on the 8-vCPU/16-GiB ECS,
at most seven independent fits run concurrently.

The registered count is:

`3 scenarios x 2 splits x 5 candidates x 4 starts = 120 endpoints`.

There are 30 `scenario x split x candidate` start-selection strata and hence
30 configuration winners. Four starts are a development-screening compromise,
not a replacement for the 12-start G12 evidence used for final scientific
fits.

## 6. Truth-free selection rule

### 6.1 Start selection within a fixed candidate

For each of the 30 strata, retain all terminal records. A start is eligible if
its ordinary T=1 ELBO is finite and objective eligibility is true. Select the
largest eligible ordinary ELBO, breaking an exact tie by lexical `fit_id`.
Practical convergence and warnings are reported but do not exclude an
otherwise eligible endpoint.

Ordinary ELBO values are never compared across dimension candidates.

### 6.2 Cross-dimension prediction loss

For subject `i`, split `r`, candidate `c`, let

`ell_irc = p^(-1) ||Y_ir,hold - E(Y_ir,hold | Y_ir,train, fit_c)||^2`.

This is prediction error for observed, noisy held-out responses. It does not
use the latent signal or any generated truth. First average the two splits
within subject, then average subjects with equal weight for the two studies.
The standard error is study-stratified with subject as the independent unit.

Let `c_min` have the smallest mean MSE. The frozen primary one-SE set is

`A = {c: R_hat(c) <= R_hat(c_min) + SE_hat[R_hat(c_min)]}`.

Select the member of `A` with the smallest registered complexity proxy; use
lexical `fit_config_id` only for an exact remaining tie. NRMSE uses the same
held-out response scale for all candidates and is reported for interpretation;
it cannot change the MSE ranking.

The following are diagnostics only and cannot replace or overturn the primary
selection:

- split-specific risks and candidate ranks;
- paired subject-level loss differences relative to `c_min` and their
  study-stratified SE;
- leave-one-subject-out selection stability;
- maximum subject influence on each paired risk contrast;
- practical convergence, warning counts and runtime.

The experimental joint log-score is not computed in this screen because Stage
6F already showed that it can reward extra low-rank covariance flexibility
while choosing the wrong dimension.

## 7. Truth boundary and post-unseal attribution

Before sealed truth can be opened, the following must exist and be hashed:

- all 120 terminal records, including any error records;
- the 30 start winners;
- all winner holdout predictions and subject-level losses;
- the three final truth-free dimension selections;
- the split-rank, paired-difference and influence diagnostics.

After a separate user authorization, evaluate all endpoints with the existing
loading-first evaluator. No new fit or continuation is allowed after unsealing.
For each scenario, attribute a failure in this order:

1. **Reachability:** did the true-dimension configuration contain a true-core
   or exact endpoint among its four starts?
2. **Within-configuration selection:** if such an endpoint existed, did maximum
   eligible ordinary ELBO select it?
3. **Cross-dimension selection:** if the true-configuration winner was
   scientifically valid, did holdout MSE plus one-SE select the correct full
   `(L_f,L_s1,L_s2,M)` tuple?
4. **Internal pruning:** when a larger candidate retained the true core, were
   duplicate/extra factors and excess FPCA components actually removed?

Truth-based evaluation retains the established metrics: observed and dense
signal reconstruction, loading-first identity, feature ISE and projection
floor, factor kernels and covariance operator, loadings, FPCA scores, factor
processes, complete contributions, missing/misplaced/duplicate/extra, PPI/PVE,
ELBO, runtime, warning and convergence descriptors. No new truth metric is
introduced.

## 8. Falsification decision and stop rule

This screen is **promising** only if the frozen primary rule selects all three
generating configurations exactly:

- minimal scenario -> `cand_L111_M2`;
- non-minimal shared-factor scenario -> `cand_L211_M2`;
- non-minimal rank scenario -> `cand_L111_M3`.

In particular, selecting `cand_L111_M2` in either non-minimal scenario
falsifies the hypothesis that the rule can move away from the simplest model.
Selecting an overfit or wrong-role alternative also fails the screen. Split and
leave-one-subject diagnostics describe fragility but cannot be used to edit the
selection after truth is known.

If any of the three selections is wrong, do not expand this criterion into a
larger multi-seed experiment. Record whether the cause was reachability,
within-configuration ELBO selection or cross-dimension risk, and then either
test at most one separately motivated alternative criterion or report `L/M` as
external sensitivity parameters. Do not tune thresholds on these unsealed
data.

If all three are correct, the result only authorizes drafting a separate
multi-seed confirmation; it is not itself a paper-level accuracy claim.

## 9. Operational gate and budget

Before ECS launch:

1. parse the eventual runner and all reused scoring/evaluation files;
2. run targeted tests for asymmetric `L_f/L_s`, `M=3` truth and both holdout
   rules;
3. run a generation-only micro smoke for all three truth dimensions using
   `jaoua_splinets`;
4. run one tiny no-truth end-to-end smoke;
5. bind the exact package source, Splinets version, protocol CSVs and runner
   hashes.

During execution, no failed fit is replaced by another seed. Selection may
proceed only if every one of the 30 strata has at least one objective-eligible
endpoint; otherwise the run stops incomplete and remains sealed.

Planning values from Stage 6D plus conservative asymmetric-factor estimates
give about 96.7 single-core hours. Seven outer workers imply approximately
14--18 wall-clock hours, allowing for scheduling and scoring overhead. Expected
fit storage is roughly 1.1--1.6 GiB. These are capacity estimates, not formal
runtime guarantees.

## 10. Authorization boundary

This draft does not authorize:

- generating the three registered data sets;
- launching any of the 120 fits;
- reading sealed truth or running truth evaluation;
- modifying the model, prior, CAVI, ELBO, G12 or generator;
- pushing, tagging or declaring a frozen paper experiment;
- starting the covariate-effect experiment.

After local review and runner validation, the user must separately authorize
the ECS falsification screen. Truth evaluation requires another authorization
after truth-free selections are frozen.
