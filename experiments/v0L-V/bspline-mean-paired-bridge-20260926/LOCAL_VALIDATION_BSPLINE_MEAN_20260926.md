# B-spline mean paired bridge: local validation record

Date: 2026-09-26

This record concerns implementation and packaging validation only. It is not a
scientific result and does not authorize truth unsealing.

## Implementation boundary

- The historical `periodic_phase` generator remains the default.
- The new `independent_splinet` option changes only the population mean truth.
- Factor functions, loadings, observation times, latent scores and observation
  noise remain paired within each logical data seed.
- The B-spline means use a periodic cubic `Splinets` system, an exact
  zero-integral coefficient subspace and continuous-L2 energy calibration.
- The fixed heat-kernel smoothness value is
  `mean_smoothness_alpha = 1.519826404906808`.
- Its calibration target is the normalized roughness of the historical
  analytic periodic mean, `64*pi^4`.

## Completed checks

- `R/generate_data_structured.R` parsed successfully.
- `run_bspline_mean_paired_bridge_20260926_v1.R` parsed successfully.
- The targeted mean-generator test completed with 64 passing expectations.
- A full `testthat` run completed with exit status 0. The 15 warnings were the
  already expected maximum-iteration warnings in legacy tests.
- The historical default generator object was byte-for-byte identical before
  and after the implementation change; its object SHA256 was
  `85191e368b5c793d9303b1611bedeb4d3320bf6454c8a2defccc039c28d8880f`.
- The runner integrity check reported two logical data seeds, two arms, 48
  fits, 12 paired starts and fitting `K = 9`.
- A from-scratch local micro smoke generated and sealed paired data, fitted
  both arms, performed truth-free selection and wrote
  `MICRO_SMOKE_COMPLETE.txt`.
- In that smoke, all pairing checks passed and the maximum residual in the
  response-difference identity was `4.440892e-16`.

## Truth boundary

The full experiment must stop after truth-free winner selection. No sealed
truth may be read and no scientific evaluation may be run without a separate,
explicit user authorization.
