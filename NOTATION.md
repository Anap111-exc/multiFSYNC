# multiFSYNC notation convention (Scheme 1)

The paper and derivations distinguish the theoretical FPCA representation from
the representation used by CAVI.

| Role | Shared component | Study-specific component |
|---|---|---|
| Theoretical FPCA score | `zeta` (variance `lambda`) | `xi` (variance `lambda`) |
| Theoretical orthonormal eigenfunction | `phi` | `psi` |
| Computational unit-variance score | `eta` | `chi` |
| Computational unconstrained time function | `theta` | `kappa` |

For an already aligned component,

`eta = zeta / sqrt(lambda)` and `theta = sqrt(lambda) * phi`,

with the analogous relationship for `chi`, `kappa`, `xi`, and `psi`.  In a
general variational fit the working components need not already point in the
eigendirections, so the canonical theoretical quantities are recovered jointly
by the covariance eigendecomposition in `orthonormalise_multi()`.

## Backward-compatible R names

The current API keeps historical names so existing simulations do not break:

| R object or function name | Statistical meaning during CAVI |
|---|---|
| `zeta`, `mu_q_zeta`, `Sigma_q_zeta`, `update_zeta()` | computational score `eta` |
| `xi`, `mu_q_xi`, `Sigma_q_xi`, `update_xi()` | computational score `chi` |
| `nu_phi`, `mu_q_nu_phi`, `Sigma_q_nu_phi`, `update_nu_phi()` | spline coefficients of `theta` |
| `nu_psi`, `mu_q_nu_psi`, `Sigma_q_nu_psi`, `update_nu_psi()` | spline coefficients of `kappa` |
| `kappa_q_*`, `lambda_q_*` | inverse-Gamma shape/scale; denoted `alpha^q`, `beta^q` in the revised derivation |

Consequently, raw CAVI objects named `zeta`, `xi`, `Phi`, or `Psi` must not be
reported as canonical FPCA quantities.  Use the post-processed outputs from
`orthonormalise_multi()` for eigenfunctions, theoretical scores, eigenvalues,
PVE, component ordering, and sign conventions.

The post-processing covariance uses complete variational second moments. Since
coefficient uncertainty can make the posterior mean covariance have rank above
the fitted computational rank `M`, the reported point estimate is the best
rank-`M` approximation; extra directions are recorded as posterior uncertainty
and are not labelled as additional population FPCA components.

## Returned scale and indexing

When `bool_scale = TRUE`, unsuffixed raw variational objects remain on the
standardised working scale for backward compatibility. Scientific summaries are
returned on the original response scale:

- `list_mu_hat[[s]][[j]]` is the mean function for study `s`, variable `j`;
- `list_beta_hat[[j]][[r]]` is the coefficient function for variable `j`,
  covariate `r`;
- `mu_q_a_hat` and `mu_q_b_specific_hat` are canonical original-scale loadings;
- `sigsq_eps_hat` is the original-scale noise variance estimate;
- fields ending in `_original` are original-scale posterior copies.

`list_mu_hat_legacy` retains the former first-study-only layout temporarily.
Factor-level `omega` is the thesis main model; variable-factor-level `omega` is
implemented as a sensitivity analysis.

Renaming every public object is intentionally deferred because it would change
the package API and invalidate existing experiment scripts.  A future major
version can expose `eta`/`chi`/`theta`/`kappa` aliases while keeping a deprecation
period for the legacy names.
