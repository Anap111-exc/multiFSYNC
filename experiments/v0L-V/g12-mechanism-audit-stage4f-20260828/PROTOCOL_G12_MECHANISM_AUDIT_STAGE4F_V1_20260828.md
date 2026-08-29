# G12 mechanism audit Stage 4F protocol v1

Registration date: 2026-08-28

Status: authorized read-only development audit; not a formal v0L-V result

## Purpose

Stage 4F reuses the completed Stage-4E endpoints and their already authorized
truth evaluation. It does not fit or continue a model. The audit asks:

1. when the study-1-specific block loses posterior inclusion support on
   `g12ss4c_sparse_01`;
2. when loading and factor-process scale become pathological on the selected
   exact-identity endpoint of `g12ss4c_sparse_02`;
3. which component of the frozen practical-convergence rule blocks the slow
   endpoints at sweep 400; and
4. what the loading-first matched factor-process and complete-contribution
   errors are for all 24 Stage-4E endpoints and the three truth-free winners.

## Frozen inputs

- Stage 4E root:
  `/root/v0lv-g12-focused-reachability-stage4e-20260827-v1`;
- inherited Stage 4C root:
  `/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1`;
- exactly 24 objective-eligible fixed-400 endpoints, comprising six inherited
  Stage-4C endpoints and 18 new Stage-4E endpoints;
- the three winners already frozen by maximum finite ordinary T=1 ELBO before
  truth evaluation; and
- the three truth bundles already unsealed under the Stage-4E authorization.

The audit must preserve the frozen winner identities. Truth may be used only
for retrospective loading, process, contribution and mechanism evaluation.

## Outputs and interpretation

The audit exports only small CSV/text derivatives:

- checkpoint factor-PPI paths for shared, study-1-specific and
  study-2-specific blocks;
- checkpoint raw loading/function/score/process scale diagnostics;
- checkpoint global fitted/RSS/contribution scale diagnostics;
- per-endpoint final practical-gate decomposition;
- loading-first matched factor-process and complete-contribution NRMSE; and
- compact block-transition and winner summaries.

PPI transition times use diagnostic thresholds 0.5, 0.1 and 0.01. These are
not new selection rules. Same-role loading cosine and norm ratios at
checkpoints are mechanism diagnostics, not replacements for the frozen
loading-first one-to-one evaluator.

## Continuation decision boundary

No continuation belongs to the read-only audit. After the audit is frozen,
at most two separately versioned 400-to-800 continuations may be run only if
the checkpoint evidence leaves genuine late movement unresolved:

- `g12ss4c_sparse_02__G12__07`, the exact-identity but scale-pathological
  Stage-4E winner; and
- one exact-identity, slow baseline endpoint as a control.

`g12ss4c_sparse_01__G12__06` is not extended merely because its structure is
wrong: it already reached practical convergence, so additional sweeps are not
a targeted test of the observed block-collapse mechanism.

## Prohibitions

Stage 4F must not change the model, prior, CAVI, ELBO, initialization,
pre-score, annealing, stopping rule, endpoint selection or historical output.
It must not generate data, run a new start, reopen winner selection, or be
combined with formal v0L-V results.
