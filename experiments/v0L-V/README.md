# v0L-V release candidate rc1

This directory packages the self-contained `release/v0L-V-rc1` candidate. It is a release candidate, not a formally frozen or authorized experiment bundle.

## Status

- The candidate2 scientific design is unchanged.
- The 50 runtime-bound files come from `rc1_snapshot/对比实验/v0L阶段A收口与v2冻结准备_20260812/protocol_manifests_20260812_v3_candidate2/FULL_RUNTIME_SOURCE_BINDINGS.csv`.
- Windows acceptance passed 18/18 candidate2 regression checks and 12/12 end-to-end smoke stages. Only the four small successful evidence files are retained under `acceptance/windows/`.
- Native Linux smoke has not yet been run.
- The binding remains `status=candidate`, `formal_execution_authorized=FALSE`, and `frozen_bundle_created=FALSE`.
- No formal data generation, fit, continuation, or evaluation is authorized from this branch.

## Layout

```text
experiments/v0L-V/
├── README.md
├── docs/HISTORICAL_EVIDENCE_INHERITANCE.md
├── acceptance/windows/
├── server/
│   ├── SERVER_DEPENDENCIES.csv
│   └── install_server_dependencies.R
└── rc1_snapshot/
    ├── Rcode/multiFSYNC/
    ├── Rcode/bayesSYNC_ref/
    └── 对比实验/...
```

`rc1_snapshot/` is the candidate project root. Its relative paths deliberately match the original project root, and Git stores all files below it without text normalization so the manifest SHA-256 values survive Windows/Linux checkout.

## Next Linux-only acceptance step

After cloning this branch on the target Linux server:

```sh
Rscript experiments/v0L-V/server/install_server_dependencies.R
cd experiments/v0L-V/rc1_snapshot
mkdir -p _r_test_lib
R CMD INSTALL --library="$PWD/_r_test_lib" Rcode/multiFSYNC
smoke_root="$(mktemp -d /tmp/multifsync-candidate2-smoke.XXXXXX)"
Rscript "对比实验/v0L阶段A收口与v2冻结准备_20260812/run_v0lv_candidate2_e2e_smoke.R" \
  --output-dir="$smoke_root/result"
```

The native Linux smoke must finish with 12/12 `PASS`; its output must stay outside the Git checkout. The fuller instructions are in `rc1_snapshot/对比实验/v0L阶段A收口与v2冻结准备_20260812/LINUX_CANDIDATE2_SMOKE_20260812.md`.

Do not use this release candidate to generate the 10 formal datasets, start the 210 formal fits, start the 72 formal audit stages, or create a frozen binding. Those actions require later user approval after native Linux smoke and timing/memory preflight.
