# G12 paper-pipeline pilot (2026-08-29)

This directory registers a small, full-size engineering pilot for the next
paper-level simulation pipeline. It is not a paper Monte Carlo result and does
not change the model, priors, CAVI updates, ELBO, G12 route, or pooled
bayesSYNC comparator.

The registered run contains two new data sets, 24 G12 fits, and two pooled
bayesSYNC fits. Fitting and G12 winner selection are truth-free. Truth
evaluation is a separate action and requires an explicit authorization file
created only after the two G12 selections have been frozen.

Static protocol inputs live here. Full fit objects, sealed truth, logs, and
evaluation products must be written to an output directory outside Git.

Main files:

- `PROTOCOL_G12_PAPER_PIPELINE_PILOT_V1_20260829.md`
- `DATA_MANIFEST.csv`, `FIT_MANIFEST.csv`, and `METRIC_MANIFEST.csv`
- `paper_pilot_common_20260829_v1.R`
- `paper_pilot_runner_20260829_v1.R`
- `prepare_server_20260829_v1.sh`
- `launch_server_20260829_v1.sh`
- `status_server_20260829_v1.sh`

The runner stops after truth-free selection unless a valid
`TRUTH_UNSEAL_AUTHORIZATION.txt` is present in the external output root.

## Server flow

From the exact Git checkout on Linux:

```bash
export PAPER_PILOT_REPO=/root/multiFSYNC-g12-paper-pilot-20260829
export PAPER_PILOT_R_LIB=/root/multiFSYNC-g12-paper-pilot-lib-20260829
export PAPER_PILOT_OUTPUT_ROOT=/root/g12-paper-pilot-20260829-v1
bash experiments/v0L-V/g12-paper-pilot-20260829/prepare_server_20260829_v1.sh
bash experiments/v0L-V/g12-paper-pilot-20260829/launch_server_20260829_v1.sh
bash experiments/v0L-V/g12-paper-pilot-20260829/status_server_20260829_v1.sh
```

The launch performs generation/sealing, all 26 fits, truth-free G12 selection,
and resource summarization. It deliberately does not run truth evaluation.
After the user separately authorizes unsealing, bind the exact selection hash
in the external authorization file and run:

```bash
Rscript experiments/v0L-V/g12-paper-pilot-20260829/paper_pilot_runner_20260829_v1.R \
  --action=evaluate \
  --output-root=/root/g12-paper-pilot-20260829-v1
```

Full data, fit, log, and evaluation directories remain outside Git.
