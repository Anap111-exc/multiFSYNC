# v0L-V frozen1

This directory contains the self-contained v0L-V release candidate and its
separate `frozen1` binding. The scientific candidate2 design and the 50 bound
runtime source files have not changed.

## Freeze status

- Runtime source commit: `b49448f851e5401cbaf681577ecee2b4bdc72514`.
- Frozen tag: `v0.3.0-v0L-V-frozen1` (created only after all checks pass).
- The target server made a real Git checkout of the release branch at the exact
  runtime source commit; `git fsck`, 50/50 SHA-256 checks, and 50/50 byte-size
  checks passed.
- After the target server's security update, the native Linux end-to-end smoke
  passed 12/12 with exit status 0.
- The frozen binding is deliberately
  `formal_execution_authorized=FALSE`. Freezing does not authorize execution.
- Formal data generation, the 210 formal fits, the 72 audit records,
  continuation, truth unsealing, and formal evaluation have not started.

## Layout

```text
experiments/v0L-V/
|-- acceptance/
|   |-- windows/
|   `-- linux/
|-- docs/HISTORICAL_EVIDENCE_INHERITANCE.md
|-- frozen1/
|   |-- README.md
|   |-- protocol_manifests_20260813_v3_frozen1/
|   `-- verify_frozen1_binding.R
|-- rc1_snapshot/
`-- server/
```

`rc1_snapshot/` remains the immutable candidate project root. Git stores every
file below it without text normalization, preserving all 50 source hashes
across Windows and Linux checkouts.

The Linux-specific environment and installed-function signatures are in the
frozen binding. Function signature hashes are intentionally captured on Linux:
they include RDS serialization and therefore are not expected to equal the
earlier Windows candidate values even though the 50 source-file hashes match.

## Allowed next step

Only inspect or verify this freeze. A later, explicit user instruction to start
the formal experiment is required before creating a separate authorized
execution binding. Do not edit this frozen binding or move the frozen tag.
