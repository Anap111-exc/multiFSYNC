# Stage 6L loading-DGP robustness bridge

This directory freezes a development-only comparison between the established
controlled identifiable loading generator and a Jaoua-style blockwise
Beta/Bernoulli loading generator. It is not a formal v0L-V or paper Monte Carlo
result and must not be used to choose the generator that makes multiFSYNC look
better.

The sequence is deliberately small:

1. audit 200 paired generated datasets without fitting;
2. run one new paired dataset with four registered G12 starts per loading arm;
3. freeze truth-free winners within arms before scientific evaluation;
4. stop loading/gate tuning and move the main development effort to the
   covariate-effect simulation requested by the supervisor.

The Jaoua-style arm implements the model-intended factor-specific law
`omega_l ~ Beta(1,10)`, `gamma_jl | omega_l ~ Bernoulli(omega_l)`, and
nonzero `loading_jl ~ N(0,1)`, with an empty block forced active. It adapts that
single-study law independently to multiFSYNC shared and study-specific blocks.
Both arms then receive the same previously established v0F constant-per-active-
loading information calibration, so the bridge is not dominated by the raw
1.6/0.8 controlled-generator norms versus unscaled standard-normal norms.
It is not described as an exact byte-for-byte reproduction of the upstream R
simulator; the source-level distinction is recorded in the protocol.

`generation_audit_v1_20260915` is the retained raw-scale preflight that exposed
the scale confounding. The registered comparison uses the subsequently added
calibrated V2 audit and bridge.

## Execution

The generation-only audit is complete: 200 pairs / 400 generated datasets,
with all 8 calibrated-audit validity checks passing. A local micro chain also
passed generation/sealing, two tiny fits, within-arm truth-free selection and
post-freeze evaluation. These tiny fits are engineering checks, not recovery
evidence.

The full bridge uses eight independent fit tasks on the 8-vCPU ECS, with one
CPU/BLAS thread per fit. Upload/install into a new Stage 6L source and R-library
directory; preserve every historical directory. Run from that source tree:

```bash
Rscript experiments/v0L-V/g12-loading-dgp-robustness-stage6l-20260915/run_stage6l_loading_bridge_20260915_v1.R --action=check
bash experiments/v0L-V/g12-loading-dgp-robustness-stage6l-20260915/launch_server_stage6l_20260915_v1.sh /root/multiFSYNC-stage6l-20260915 /root/multiFSYNC-stage6l-lib-20260915 /root/multifsync-results/stage6l_20260915_v1 /root/multiFSYNC-jaoua-paired-lib-20260909
```

The launcher ends at the frozen truth-free selection, keeps truth sealed and
uses the previously configured instance-role self-stop permission only when
no other R workload is running. The separate evaluator requires explicit
truth-unseal authorization and cannot start fits or continuation. No Git push
or historical frozen-binding update is part of this step.

Its optional fourth argument reuses the already installed dependency library
(including Splinets). The new Stage 6L model library is always first in
`R_LIBS`, so this does not run or overwrite the historical multiFSYNC package.
