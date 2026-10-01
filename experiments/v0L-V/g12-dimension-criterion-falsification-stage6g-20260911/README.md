# Stage 6G: compact falsification screen for the dimension criterion

Status: truth-free ECS screen authorized on 2026-09-11; truth unsealing and
continuation remain unauthorized.

This directory specifies a three-scenario screen of the truth-free paired
time-block MSE/NRMSE plus one-standard-error dimension rule. It does not modify
the multiFSYNC model, priors, CAVI, ELBO, G12 initialization, annealing, or the
400-sweep development endpoint.

The key addition relative to Stages 6D and 6F is that two scenarios have a
generating dimension that is not the simplest candidate. A common five-model
candidate set is used for all scenarios, so the grid is not tailored after
seeing which model generated a data set.

The registered design implies 3 data sets x 2 holdout splits x 5 dimension
configurations x 4 paired starts = 120 endpoints and 30 configuration/split
winners. The full runner reuses the tested Stage 6F holdout and scoring code
and adds only the registered asymmetric factor-count configurations.

Run the local, non-fitting consistency check with:

```sh
Rscript validate_stage6g_protocol_20260911_v1.R
```

The separate development-only micro smoke is:

```sh
Rscript run_stage6g_local_micro_smoke_20260911_v1.R \
  --package-root=/path/to/multiFSYNC \
  --output-root=/tmp/stage6g-micro-smoke
```

It uses reduced data and six-sweep fits only; it cannot be interpreted as a
criterion-performance experiment. The completed local result is recorded in
`LOCAL_MICRO_SMOKE_RESULT_20260911.md`.

Before server launch, targeted tests can be run with:

```sh
Rscript test_stage6g_tools_20260911_v1.R \
  --output-root=/tmp/stage6g-targeted-tests
```

The authorized Linux launcher is `launch_server_stage6g_20260911_v1.sh`. It
uses seven independent fit workers with `n_cpus=1` inside every fit. The
separate `watch_and_stop_server_stage6g_20260911_v1.sh` requests StopCharging
after the launcher exits, retaining either the complete or incomplete state.
The
authorization covers only registered generation/sealing, 120 endpoints, start
selection, winner scoring and the three truth-free dimension selections. Do
not unseal truth, run continuation, selectively rerun failures or interpret
this development screen as paper-level Monte Carlo evidence.
