# Local validation: Jaoua Splinets generator transition

Date: 2026-09-09

Environment: Windows, R 4.5.1, Splinets 1.5.1. Validation used an
ASCII-only temporary package copy because the existing Windows `pkgload/desc`
stack cannot resolve the Chinese project path. The temporary copy contained the
current working-tree package files and was not a separate model version.

## Results

- all 27 package R files parsed successfully;
- targeted structured-generator tests passed 104/104 with no warnings or
  skips;
- the full package testthat suite passed 1143/1143;
- the full suite emitted 15 expected max-iteration warnings from pre-existing
  short fitting tests and had zero failures and zero skips;
- a generation-only 1/1/1-factor, M=2 micro smoke passed;
- continuous spline Gram maximum absolute error:
  `8.8817841970012523e-16`;
- 201-point trapezoidal-grid Gram maximum absolute error:
  `5.9237999394667895e-08`;
- response = signal + noise maximum absolute error:
  `2.2204460492503131e-16`;
- generation elapsed time: `0.19` seconds;
- formal data generated: FALSE;
- fit started: FALSE.

The difference between continuous and grid Gram errors is expected: Splinets
orthogonalises and evaluates the continuous spline representation, whereas the
reported dense-grid check is a separate trapezoidal approximation.

## Remaining gate

A separate protocol decision is still required to make this option primary for
a new paper experiment; no historical experiment is replaced by this result.

## Target Linux micro smoke

Date: 2026-09-09

Environment: Alibaba Cloud Linux, R 4.5.3, Splinets 1.5.1. A 134,105-byte
development source tarball was transferred to `/tmp`, verified against SHA256
`4bd4e0a6e0c898617800f4770dcab2323105e689dece12a8e76c00ff7df8be9b`, and
installed into an isolated temporary R library. It did not replace the existing
server checkout or R library.

- generation-only micro smoke: PASS;
- factor blocks: 3 (shared 1, study-specific 1 per study);
- continuous spline Gram maximum absolute error:
  `8.8817841970012523e-16`;
- 201-point trapezoidal-grid Gram maximum absolute error:
  `5.9237999394667895e-08`;
- response = signal + noise maximum absolute error:
  `2.2204460492503131e-16`;
- generation elapsed time: `0.094` seconds;
- formal data generated: FALSE;
- fit started: FALSE.
