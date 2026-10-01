# Jaoua Splinets generator transition protocol v1

Date: 2026-09-09

Status: development preparation only. This protocol does not authorize formal
data generation, fitting, continuation, truth evaluation, or replacement of
any frozen v0L-V result.

## Purpose

Add an explicit data-generation option that follows the current upstream
bayesSYNC construction of orthonormal B-spline systems through the `Splinets`
package, while preserving the historical dense-grid Gram generator byte-for-
byte at its existing default.

## Source binding

- bayesSYNC repository: `https://github.com/hruffieux/bayesSYNC`
- inspected master commit: `de326142f15c8a087f544c84f18d83511aae50f1`
- commit adding `R/simulate_data.R`: `816db2fe27bb6f596012135caaea3c8852eb8da8`
- SHA256 of bound `R/simulate_data.R`:
  `58d981720c0a639a4f941fa4a1f1353255721258b261519fce671dcdd8517639`
- Splinets package version for development verification: `1.5.1`
- Splinets upstream commit named by bayesSYNC:
  `a8f79565257da5d292acbd4fc03af62ecd3f899e`

The upstream simulation function was first committed on 2026-07-14. It is a
current author-maintained implementation consistent with the paper's stated
orthonormalised B-spline construction, but it is not treated as proof of the
exact unpublished script used for every numerical result in the paper.

## Interface and compatibility contract

- Add `truth_basis_method = c("legacy_dense_gram", "jaoua_splinets")` at the
  end of `simulate_multi_study_structured()`.
- Keep `legacy_dense_gram` as the package default so historical runners do not
  silently generate different data when re-executed.
- New simulation protocols may explicitly select `jaoua_splinets` only after
  the tests below pass.
- The new option changes only the truth-basis generator. It does not change the
  loading generator, score variances, model, prior, CAVI, ELBO, G12
  initialization, annealing, multistart selection, or evaluation formulas.

## Exact supported domain

The source implementation assumes one common FPCA count for all factors and
supports at most four factors and four components. For multiFSYNC, factors are
ordered as all shared factors followed by study-specific factors in study and
factor order. The exact option therefore requires:

- total `L_f + sum(L_s) <= 4`;
- all entries of `M_f` and `M_s` equal;
- common `M` between 2 and 4.

The lower limit records verified behaviour of Splinets 1.5.1 under the bound
knot construction. The planned 1/1/1-factor, M=2 simulations lie inside this
domain. Unsupported dimensions must fail explicitly rather than silently cycle
or duplicate factor bases.

## Required verification before any simulation study

1. Parse all package R files.
2. Run targeted structured-generator tests.
3. Verify continuous spline Gram matrices against the identity.
4. Verify observed-time values are evaluated directly by `evspline()`, not by
   interpolation from the stored dense grid.
5. Verify default and explicit `legacy_dense_gram` outputs are identical under
   the same seed.
6. Run one tiny generation-only smoke with the planned 1/1/1, M=2 dimensions.
7. Record Windows and target-Linux package versions before a new experiment is
   frozen.

## Decision boundary

Passing these checks makes the generator available for a paired sensitivity
comparison. It does not by itself justify discarding or rerunning historical
evidence. A later, separately reviewed protocol must decide whether the
Splinets generator becomes the primary paper data-generating mechanism.
