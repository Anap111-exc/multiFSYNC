# G12 Stage 6E-A read-only selection attribution

This directory contains a local, read-only audit of the 180 completed Stage 6D
fits. It does not run a fit or continuation and does not use ECS.

- Protocol: `PROTOCOL_G12_STAGE6E_A_READONLY_SELECTION_ATTRIBUTION_V1_20260908.md`
- Tools: `stage6e_a_tools_20260908_v1.R`
- Runner: `run_stage6e_a_readonly_audit_20260908_v1.R`
- Targeted tests: `test_stage6e_a_readonly_audit_20260908_v1.R`
- Valid machine output: `stage6e_a_readonly_audit_20260908_v2/`
- Scientific report: `REPORT_G12_STAGE6E_A_READONLY_SELECTION_ATTRIBUTION_V1_20260908.md`

The `stage6e_a_readonly_audit_20260908_v1/` output was invalidated before
reporting because one descriptive spectral field propagated missing padded
eigenvalues. The bug did not affect risk or selection, but v2 is the only
output used in the report.

Stage 6D truth was already unsealed before this audit. Any apparent dimension
accuracy in Stage 6E-A is post-hoc exploratory evidence, not independent
validation of a selection rule.
