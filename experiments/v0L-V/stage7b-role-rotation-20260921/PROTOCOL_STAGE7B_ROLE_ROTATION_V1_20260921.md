# Stage 7B: fixed-seed truth-block role rotation (developmental causal audit)

Status: candidate protocol, not a paper Monte Carlo result or a change to Stage 7A.
No model, prior, CAVI, ELBO, initialization, stopping rule, or truth-free
selection rule changes are authorized by this protocol.

## Question and competing explanations

In Stage 7A, spline-system blocks `q=1,2,3` were assigned to shared, study-1
specific, and study-2 specific factors, respectively. The six-seed no-fit audit
found that block 3 has a much larger projection residual in the fitted
O'Sullivan `K=5` space than blocks 1 and 2. It does **not** establish that the
fitted dense-domain error is caused solely by the truth block: block identity
and study role were completely confounded.

This audit asks whether the fitted error follows block 3, remains with study 2,
or depends on their interaction. It does not estimate population-level
performance or establish a new default generator.

## Paired design fixed before new fits

- Use the previously evaluated Stage 7A `s7a_03__sine` data as the baseline
  `(shared, study 1, study 2) = (q1, q2, q3)`. Do not rerun its 12 fits.
- Generate only `(q2, q3, q1)` and `(q3, q1, q2)` under the original logical
  data seed `97201003`, basis seed `97202003`, and loading seed `97203003`.
  Hold observation times, subjects, covariates, population means, FPCA scores,
  loading vectors/support, measurement errors, and beta functions fixed.
  Only the assignment of the three already-drawn truth function blocks and
  consequently the signal/response changes.
- One sine-beta arm; `S=2`, `n_s=(30,30)`, `p=500`, `d=1`, `L_f=1`,
  `L_s=(1,1)`, `M=2` per factor, `K=5`, 6--9 irregular observations per
  subject, `sigma_eps=0.3`, `mean_amp=0.6`, `beta_amp=0.3`.
- Twelve starts per new assignment use exactly Stage 7A seed-03 fit seeds
  `97210301:97210312`. Total new full-size fits: 24.
- Reuse the byte-identical Stage 7A source and fitting worker. Each fit uses
  `n_cpus=1`, G12 random function-energy initialization, one pre-score sweep,
  Jaoua annealing `c(1,1.9,100)`, and 400 ordinary `T=1` sweeps.
- Retain all COMPLETE/ERROR records. Select the maximum eligible final ordinary
  `T=1` ELBO **within each rotated dataset** without reading its truth.
  Different datasets have different responses: do not compare their raw ELBO
  values as if they were a cross-dataset model-selection criterion.
- Before any full fit, verify that the identity role mapping reconstructs the
  original Stage 7A observation and truth fields to numerical tolerance,
  verify all intended paired fields remain equal under rotations, and pass a
  separate micro smoke. Abort if either check fails.

## Post-selection scientific audit (not part of fit selection)

After both truth-free winners are frozen and truth-unsealing authorization is
recorded, compare the three assignments on: loading-first shared/specific
identity and duplicate/extra roles; by-role observed and dense signal and
complete-contribution error; matched feature ISE and projection-excess ISE;
factor-process and covariance-kernel errors; strict convergence and runtime.
Report all starts as well as selected winners, and retain failed starts.

Interpretation is conditional on one deliberately selected difficult data
seed. Error following `q3` supports a template/representation mechanism;
error remaining with study 2 supports a study-role/design mechanism; a special
problem only for `q3` in study 2 supports interaction. Mixed results remain
inconclusive. No generator or `K` change follows automatically.

Stage 7A results and source files are immutable. All Stage 7B outputs go to a
new directory outside the Git source tree. The new sealed truth is not read by
fit or winner selection.
