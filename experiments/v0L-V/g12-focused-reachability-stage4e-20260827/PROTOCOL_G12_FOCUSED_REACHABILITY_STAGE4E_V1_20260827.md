# G12 focused reachability Stage 4E protocol v1

Registration date: 2026-08-27
Status: authorized adaptive development experiment; not a formal v0L-V result

## Question and scope

Stage 4E asks a single focused question: for the three Stage-4C data sets on
which neither of the first two G12 endpoints recovered the exact shared and
study-specific identity, is the correct basin reached when G12 start coverage
is increased from two to eight starts?

The targets are fixed from the already authorized Stage-4D truth evaluation:

- `g12ss4c_base_01`: both prior endpoints selected 1/1/1 factors but confused
  shared and study-2-specific identity;
- `g12ss4c_sparse_01` and `g12ss4c_sparse_02`: both prior endpoints omitted the
  study-1-specific factor.

This targeting is adaptive and post-truth. Therefore Stage 4E is development
evidence only. It must not be pooled with Stage 4C as an unseen confirmation,
used to revise the frozen v0L-V result, or reported as a pre-registered
cross-scenario success rate.

## Frozen fit design

No new data are generated. The three exact Stage-4C observation-only bundles
are reused after SHA256 verification. Their sealed truth bundles are not copied
into the Stage-4E fitting root.

For each target data set, six new G12 starts are registered (`seed_index=3` to
`8`), giving 18 new fits and eight total endpoints after combining the two
inherited Stage-4C endpoints. New fit seeds are:

- base 01: `87266013` through `87266018`;
- sparse 01: `87266053` through `87266058`;
- sparse 02: `87266063` through `87266068`.

Every new fit retains the established route: random initialization with
Gram-unit-energy calibration, one bounded pre-score sweep, Jaoua annealing
`c(1,1.9,100)`, and the public `g12_stopping_control()` profile. Each fit has
`n_cpus=1`, all BLAS-family thread counts equal one, and a fixed 400 ordinary
`T=1` sweep endpoint (`maxit=499`). Four independent fits may run concurrently.
There is no racing, continuation, automatic 800-sweep extension, model change,
prior change, CAVI change, or ELBO change.

## Truth-free selection and later evaluation

Fitting, stopping, objective-eligibility checks, and endpoint selection cannot
read truth. After all 18 new endpoints complete, the two inherited and six new
eligible endpoints are combined within each data set. The winner is frozen by
maximum finite ordinary `T=1` ELBO, with ascending `fit_id` only as a
deterministic exact-tie rule. Structure labels, NRMSE, ISE, R/P/L, and all truth
metrics are excluded from selection.

Only after `TRUTH_FREE_SELECTION_COMPLETE.txt` exists may the truth that was
already unsealed for Stage 4D be read for Stage-4E evaluation. All 24 endpoints
are evaluated, not only winners. The primary reachability outcomes are exact
identity hits among eight starts, whether the ELBO winner is exact, and the
ELBO ordering of exact versus incorrect endpoints. Loading-first identity,
reconstruction, feature/function, factor-process, FPCA-score, covariance and
complete-contribution outputs retain the frozen evaluator semantics.

## Decision interpretation

- At least one exact endpoint and an exact ELBO winner: the basin is reachable
  under eight starts and the ordinary ELBO supports it on that data set.
- At least one exact endpoint but an incorrect ELBO winner: reachability exists,
  but the current objective/prior ordering requires investigation.
- No exact endpoint in eight starts: the failure remains a reachability or
  information/selection problem; Stage 4E alone cannot distinguish them.

Three adaptively selected data sets are insufficient for population-level
success claims. Results guide the next development decision only.

Expected ECS wall time is approximately 9–11 hours for the 18 fits with four
workers, followed by a shorter sequential truth evaluation. Expected added
disk use is roughly 1.5–2.0 GiB.
