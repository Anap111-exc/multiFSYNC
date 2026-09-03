# Stage 6B predictive pilot v2 bugfix rerun

Date: 2026-09-03

Status: user-authorized technical rerun; truth remains sealed

The first full-size execution at
`/root/g12-truth-free-dimension-selection-stage6b-20260902-v1` is retained
read-only. All 12 registered factor-screen fits reached the post-fit code, but
eight fits with an empty warning table stopped before `fit.rds` persistence
when one phase label was assigned to a zero-row data frame. Four fits were
persisted; no start winner, factor count, FPCA cap or final truth-free selection
was frozen, and no FPCA-screen fit was started.

The technical fix:

1. assigns warning phases with `rep(phase, nrow(warnings))`, which is valid for
   both zero and nonzero rows;
2. saves every non-null fit object immediately after the fitter returns and
   before nonessential warning labelling or holdout scoring;
3. adds a targeted regression test for zero-warning tables and persistence
   ordering.

The scientific design is unchanged: the same observation-only bundle, fixed
holdout, nine registered configurations, four paired seeds, G12 initialization,
one pre-score, Jaoua annealing, 400 ordinary T=1 sweeps, one CPU per fit,
maximum-ELBO within-configuration start selection and stratified one-standard-
error dimension selection are retained exactly.

The clean rerun must use
`/root/g12-truth-free-dimension-selection-stage6b-20260903-v2`. It must rerun
all 12 factor-screen fits from scratch and must not reuse the four persisted v1
fits. Only after the factor-screen selection is frozen may the eight conditional
FPCA-screen fits run. Truth access, continuation and full-data refitting remain
forbidden until the final truth-free dimension selection is frozen.
