# G12 targeted continuation Stage 4G protocol v1

Registration date: 2026-08-28

Status: authorized adaptive development audit; not a formal v0L-V result

## Question

Stage 4F found two distinct sweep-400 states that justify a minimal 400-to-800
audit:

- `g12ss4c_sparse_02__G12__07` is the truth-free Stage-4E winner and has exact
  shared/specific identity, but its dense shared and study-1-specific process
  scale and its fitted/RSS trajectory are still moving materially;
- `g12ss4c_base_01__G12__04` is an exact-identity baseline runner-up whose only
  final gate failure is a marginal short-window RSS change. It is the stable
  control.

The experiment asks whether another 400 ordinary T=1 sweeps stabilize or
worsen the sparse endpoint, and whether the marginal baseline gate failure
resolves without scientifically important change.

Target choice is adaptive and post-truth. Results are mechanism evidence only;
they are not an unseen success rate and cannot replace Stage-4E winners.

## Fixed execution

Each source is continued from its complete sweep-400 variational state for
exactly 400 additional ordinary T=1 sweeps, ending at cumulative sweep 800.
There is no random reinitialization, pre-score or annealing. Each fit uses
`n_cpus=1`; the two independent continuations run concurrently. The public
G12 diagnostic tolerances are retained, while `consecutive=401` prevents
early termination so the fixed horizon is observed.

Checkpoints are cumulative sweeps 400, 420, ..., 800. The source-to-
continuation ELBO bridge, continuation objective eligibility, state finiteness,
source snapshot identity and complete checkpoint path must pass. Errors are
retained and are not selectively rerun.

## Isolation and endpoint semantics

The continuation runner reads observation-only bundles and cannot read truth.
It uses the same model/package source as Stage 4E. It does not select or rank
the two endpoints. After both endpoints complete, the truth already authorized
for Stage 4E may be used to compare the fixed source and continuation states.

The cumulative-800 objects are counterfactual audit endpoints. They do not
replace the cumulative-400 source objects, do not reopen the Stage-4E ELBO
winner, and do not change the current G12 budget or stopping rule.

## Outputs

- two independent continuation directories with complete RDS, ELBO,
  practical diagnostics, scale trace, warnings/errors, runtime and
  COMPLETE/ERROR marker;
- one two-row terminal table;
- source-versus-continuation truth metrics, including structure,
  reconstruction, loading, feature, process, score, covariance and complete
  contribution outputs; and
- a compact 400-to-800 mechanism trajectory and convergence-gate audit.

No further 800 extension or new fit is automatic.
