# Stage 4F: read-only mechanism audit

This directory contains a development-only audit of the completed Stage-4E
fixed-400 endpoints. Run `s4f_readonly_mechanism_audit_20260828_v1.R` from a
separate output directory on the ECS. The script reads Stage 4C/4E RDS files
but writes only compact derived CSV/text files under its own directory.

It runs no fit and no continuation, does not reopen truth-free selection, and
does not modify the frozen v0L-V result.

Default ECS command:

```bash
Rscript s4f_readonly_mechanism_audit_20260828_v1.R
```

Input roots can be overridden with `STAGE4F_STAGE4E_ROOT`,
`STAGE4F_STAGE4C_ROOT`, and `STAGE4F_SNAPSHOT_ROOT`.

The completed development report is
`REPORT_G12_MECHANISM_AUDIT_STAGE4F_20260828_V1.md`. The `evidence/`
directory contains only the compact completion, QC, transition, terminal-gate,
and process/contribution summaries used by that report. Full fit objects and
large checkpoint tables remain outside Git.
