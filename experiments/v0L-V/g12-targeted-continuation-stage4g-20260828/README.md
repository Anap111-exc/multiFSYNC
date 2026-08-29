# Stage 4G: two targeted 400-to-800 continuations

This is a minimal adaptive development audit triggered by the Stage-4F
read-only mechanism results. It continues exactly two existing Stage-4E
endpoints and never reopens truth-free winner selection.

Execution order:

```bash
Rscript build_stage4g_manifest_20260828_v1.R
Rscript build_stage4g_source_binding_20260828_v1.R
Rscript s4g_targeted_continuation_20260828_v1.R --check-only=true
Rscript s4g_targeted_continuation_20260828_v1.R
Rscript s4g_evaluate_continuations_20260828_v1.R
```

For a detached ECS run, `evaluate_stage4g_after_completion_20260828_v1.sh`
may wait for the continuation terminal marker and launch evaluation only after
successful completion. It exits without evaluation if an error is retained.

The default ECS output root is the directory containing these scripts. Inputs
remain in the completed Stage-4C/4E roots. This is not a formal v0L-V run.

The completed development report is
`REPORT_G12_TARGETED_CONTINUATION_STAGE4G_20260828_V1.md`. The `evidence/`
directory contains the compact manifest, source binding, terminal, gate,
scientific-change, process/contribution, and QC summaries. Continuation RDS,
full traces, logs, and archives remain outside Git.
