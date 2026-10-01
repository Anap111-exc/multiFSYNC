# Paired-confirmation micro validation

Date: 2026-09-09

## Package regression

- all package R files parsed;
- targeted structured-generator test file: 136 expectations, 0 failures,
  0 errors, 0 warnings and 0 skips;
- complete package testthat suite: 1190 expectations passed, 0 failures,
  0 errors and 0 skips;
- 15 pre-existing expected max-iteration warnings were retained.

The Windows tests used an ASCII-only temporary package copy to avoid the known
`pkgload/desc` problem with the Chinese project path. It contained the current
working-tree bytes and did not constitute a separate model implementation.

## Local end-to-end micro smoke

- source binding: 41/41;
- paired data generation/sealing: 3 logical pairs and 6 observation bundles;
- exact non-basis RNG pairing checks: PASS for all pairs;
- G12 micro fits: 2/2 complete and objective eligible;
- truth-free winners: 2/2 frozen;
- truth unsealed: FALSE;
- formal data/fits started: FALSE.

## Full development-run launch

After the Linux micro smoke passed, the two minimal operational scripts were
added and the target-server source binding was rechecked at 43/43. The new
full-size output root is `/root/jspc-output-20260909-v1`; no earlier checkout,
library or result directory was overwritten.

- launch time: 2026-09-09 20:00 China Standard Time;
- registered full-size fits: 72;
- technical-gate fits: 8, all included in the 72;
- outer workers: 7;
- per-fit `n_cpus`: 1;
- initial observed compute processes: 7 at approximately one CPU core each;
- initial error markers: 0;
- automatic `StopCharging` watcher: armed and conditional on 72/72 complete,
  zero error markers, a frozen truth-free selection and an exited runner;
- truth unseal/evaluation: not authorized and not started.

## Target Linux end-to-end micro smoke

Host: Alibaba Cloud Linux, R 4.5.3, Splinets 1.5.1. Source was installed into a
new isolated R library and did not overwrite earlier checkouts or libraries.

- uploaded archive SHA256:
  `741bf622ec78ea6ae13b7e6b6ec73c0133ee0d208f3e436d5b0ecc91f41f90b9`;
- source binding: 41/41 before the operational launch scripts were appended;
- paired data generation/sealing: PASS;
- G12 micro fits: 2/2 complete and objective eligible;
- truth-free winners: 2/2 frozen;
- elapsed wall time: 1.44 seconds;
- maximum resident set size: 104428 KiB;
- truth unsealed: FALSE;
- formal data/fits started: FALSE.
