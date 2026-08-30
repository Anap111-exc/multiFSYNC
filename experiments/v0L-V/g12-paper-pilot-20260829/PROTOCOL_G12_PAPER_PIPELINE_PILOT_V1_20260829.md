# G12 paper-pipeline pilot protocol v1

Date: 2026-08-29
Protocol ID: `G12_PAPER_PIPELINE_PILOT_V1_20260829`
Status: registered engineering pilot
Formal v0L-V result: `FALSE`
Formal paper Monte Carlo result: `FALSE`

## 1. Purpose and claim boundary

This pilot checks whether the next paper-level pipeline can reproducibly:

1. generate and seal full-size data;
2. run the registered G12 route and pooled bayesSYNC comparator;
3. retain complete terminal records and freeze truth-free G12 winners;
4. unseal truth only after separate authorization;
5. produce the preregistered scientific metrics; and
6. measure wall time, peak process memory, disk use, and failure modes on the
   target 4-vCPU/8-GiB Linux ECS.

Two data sets cannot estimate a scientific success probability, scenario
effect, or method superiority. Any scientific values produced by this pilot
are descriptive pipeline checks only. The pilot neither authorizes nor
replaces a formal paper Monte Carlo experiment.

## 2. Data registration

Two new seeds are registered before generation:

| data_id | seed | scenario | observation times |
|---|---:|---|---|
| `g12pp_base_01` | 82929001 | baseline strong identification | 6--9 irregular |
| `g12pp_weak_01` | 82929002 | weak loading separation, target absolute cosine 0.6 | 6--9 irregular |

Both use:

- `S=2`, `n_s=(30,30)`, `p=500`, `d=0`;
- true shared count 1 and study-specific counts `(1,1)`;
- `M=2` per factor and `K=5`;
- sparse loadings with active fraction 0.1 under the frozen information
  calibration;
- `sigma_eps=0.3`, `mean_amp=0.6`, `beta_amp=0`;
- irregular observation grids and a 201-point dense truth grid;
- the unchanged structured-data generator whose SHA256 is recorded in
  `SOURCE_BINDING.csv`.

For the weak scenario, the already reviewed Stage-4C transform sets each
shared--specific loading absolute cosine to 0.6 while preserving loading norm
and active count, then rebuilds the observations. No fit sees the scenario
label or sealed truth.

## 3. Registered fits

### G12 multiFSYNC

Each data set receives 12 preregistered random starts (24 fits total):

- `initialization="random"`;
- `function_initialization="gram_unit_energy"`;
- one bounded pre-score update;
- Jaoua annealing `c(1,1.9,100)` (99 annealed sweeps);
- ordinary `T=1` CAVI to the common fixed sweep-400 endpoint;
- public `g12_stopping_control()` with `min_t1=396`, `consecutive=5`,
  `max_t1=400`, plus diagnostic checkpoints at 380 and 400;
- `lambda_orth=0`, `n_cpus=1`, and one BLAS/OpenMP thread;
- no racing, continuation, or automatic sweep-800 extension.

Objective eligibility and practical convergence remain separate. Within each
data set, the G12 winner is the eligible endpoint with maximum ordinary `T=1`
ELBO; lexical `fit_id` is the preregistered tie break. Truth, structure labels,
NRMSE, ISE, loading recovery, factor processes, contributions, and R/P/L are
forbidden from fitting, stopping, eligibility, and selection.

### pooled bayesSYNC

The two studies are concatenated by subject and fitted once per data set,
without a study label. Each pooled fit uses `Q=3`, `L=2`, `K=5`, `n_g=51`,
Jaoua annealing `c(1,1.9,100)`, `maxit=299`, `n_cpus=1`, and the frozen
bayesSYNC reference source. There are exactly two pooled fits. Pooled and
multiFSYNC ELBO values are never compared across models.

Total registered fits: `24 + 2 = 26`.

## 4. Truth isolation and terminal policy

- Generation writes an observation-only bundle and a separate sealed-truth
  bundle for each data set.
- Fit workers and truth-free selection may read only observation bundles.
- Every fit has an independent directory containing the full fit RDS, terminal
  record, objective trace, warnings/errors, elapsed time, peak memory when
  available, and `FIT_COMPLETE.txt` or `FIT_ERROR.txt`.
- A failed fit is retained and is not selectively rerun.
- Evaluation is blocked until all 26 terminal records exist, the two G12
  winners are frozen, and a separate authorization file exactly binds the
  frozen selection hash.
- No new fit or continuation may start after truth is unsealed.

## 5. Preregistered evaluation outputs

For G12, the frozen candidate2 loading-first evaluator is reused without
changing its formulas. For pooled bayesSYNC, all three oracle loading
directions are matched one-to-one before function, process, score, and
contribution evaluation. Zero loading directions remain finite candidates with
directional similarity zero.

The output hierarchy is:

1. observed-time and dense-grid signal reconstruction;
2. loading-first shared/specific identity, missing, misplaced, duplicate, and
   extra components;
3. feature total ISE, projection floor, excess ISE, and matched/total coverage;
4. loading direction, support, and scale recovery;
5. FPCA score correlation and factor-process recovery;
6. complete-contribution and covariance-kernel/operator recovery;
7. ELBO, elapsed time, peak memory, warning/error, objective eligibility,
   practical status, and 380-to-400 fitted/RSS/PPI stability;
8. selected factor counts and retained M only as auxiliary descriptions.

Raw and canonical loading norms, factor scales, reconstruction gaps, process
error, and contribution error are retained when supplied by the evaluator.
No new success threshold is introduced in this pilot.

## 6. Engineering acceptance gates

The pilot pipeline passes only if:

- manifests contain exactly 2 data rows and 26 unique fit rows;
- both observation bundles and sealed-truth bundles are created with no truth
  field in the observation bundle;
- 26 terminal records are retained and any failure is reported without
  selective rerun;
- both G12 selections use only eligible ordinary ELBO endpoints and are frozen
  before truth access;
- each pooled baseline combines the two studies and fits exactly once without
  a study label;
- after explicit authorization, all registered evaluator tables are produced
  or a retained evaluation error identifies the blocking defect;
- wall time, per-fit elapsed time, peak process memory, and final disk usage are
  reported.

These are pipeline gates, not statistical success gates.

## 7. Compute budget

The target is the existing Alibaba Cloud Linux ECS with 4 vCPU and 8 GiB RAM.
Use three independent fit workers; every fit and BLAS/OpenMP runtime stays at
one thread. The preregistered planning envelope is:

- expected fitting wall time: 10--14 hours;
- engineering stop/review boundary: 18 hours;
- expected peak aggregate memory: below 6 GiB;
- required free disk before launch: at least 15 GiB;
- no swap-dependent completion claim.

The measured pilot values replace these planning estimates for the later
paper Monte Carlo budget. They do not change the G12 scientific route.

## 8. Authorization boundary

The current user request authorizes implementation, nonformal smoke, generation
of these two pilot data sets, and the 26 registered pilot fits. It does not by
itself authorize truth unsealing, a full paper Monte Carlo, a Git tag, a model
change, or a prior/ELBO/CAVI change. Truth evaluation requires a separate
explicit user authorization after selection freeze.
