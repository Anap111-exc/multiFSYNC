# Stage 4F/4G/5A development evidence freeze

Freeze date: 2026-08-29
Scope: compact Git evidence only; full local archives remain authoritative

## Stage 4F

- Read-only audit of 24 completed fixed-400 endpoints; 13/13 QC passed.
- Sparse-01 block loss occurred before the first ordinary checkpoint at sweep
  20 and was not a late-iteration failure.
- Sparse-02 had correct loading identity but severe dense-domain and scale
  error.
- Fixed-400 slow cases were blocked by fitted/RSS movement, not PPI gates.

Git evidence is stored under
`../g12-mechanism-audit-stage4f-20260828/`. Full RDS and checkpoint tables are
not committed.

## Stage 4G

- Exactly two preregistered 400-to-800 continuations completed and passed
  evaluation QC.
- Both became numerically stable by sweep 800, but shared/specific identity was
  unchanged and scientific recovery did not improve.
- Sweep 800 is therefore retained only as a small audit budget, not a default.

Git evidence is stored under
`../g12-targeted-continuation-stage4g-20260828/`. Continuation RDS, traces,
logs, and the full archive are not committed.

## Stage 5A

- Read-only audit covered 104 standard-density G12 endpoints from 10 data sets;
  45 endpoints had exact identity and all 10 frozen maximum-ELBO winners were
  exact.
- Exact endpoints contained a small number of scale-inflated scientific
  sub-basins; the frozen winners avoided the catastrophic dense/loading/
  covariance tail.
- Strict practical convergence was 0/104 and could not serve as a structural
  success criterion. Auxiliary endpoint stability carried some conditional
  quality information but did not identify correct shared/specific structure.

Git evidence is stored under
`../g12-correct-basin-quality-audit-stage5a-20260829/`. Historical input RDS
remain outside Git.

## Frozen development decision

The next development candidate is G12 plus one bounded pre-score, Jaoua
annealing, 12 starts, fixed sweep-400 observation, and maximum finite ordinary
T=1 ELBO selection. This is a mitigation for nonconvex reachability, not a
claim that optimization is solved. The statistical model, priors, CAVI, ELBO,
and historical protocols are unchanged.
