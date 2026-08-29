# Stage 5A read-only audit

Run from the local project checkout:

```powershell
Rscript s5a_readonly_correct_basin_audit_20260829_v1.R
```

The script reads existing independent-confirmation and cross-scenario aggregate results. It recomputes process/contribution and loading-scale outputs only from the primary cohort's existing RDS files. It never calls a data generator or fitting function and refuses to overwrite its versioned output directory.

Set `STAGE5A_CHECK_ONLY=1` for input preflight without creating outputs. Input and output roots may be overridden with the `STAGE5A_*` environment variables defined near the top of the script.

The completed development report is
`REPORT_G12_CORRECT_BASIN_QUALITY_AUDIT_STAGE5A_20260829_V1.md`. The
`evidence/` directory contains only compact input bindings, completion/QC,
per-data basin, winner, scientific-quality, scale, process/contribution, ELBO,
and convergence summaries. Historical RDS inputs remain outside Git.
