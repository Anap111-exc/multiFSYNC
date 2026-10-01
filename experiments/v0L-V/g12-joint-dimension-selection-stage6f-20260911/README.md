# Stage 6F: joint dimension-selection development experiment

This directory contains a small development-only experiment for assessing a
truth-free dimension-selection criterion. It does not alter the multiFSYNC
model, priors, ELBO, G12 initialization, annealing, stopping rule, or the
historical frozen v0L-V results.

The registered candidate grid is the joint 2 by 2 set
`L1M2`, `L1M4`, `L2M2`, and `L2M4`, where the common factor count is used for
the shared block and both study-specific blocks and the common FPCA cap is used
for every fitted factor. The experiment uses one existing development data set,
two deterministic paired whole-time holdout splits, and four paired G12 start
seeds per configuration.

There are 32 endpoints in the analysis:

- 12 Stage 6B endpoints are reused byte-for-byte for the first split;
- 4 new `L2M4` endpoints complete the first split;
- 16 new endpoints fit all four configurations on the second split.

Thus, exactly 20 new fits are authorized. Each new fit uses one CPU, one
pre-score sweep, Gram-unit-energy random function initialization, Jaoua's
default annealing schedule, and a fixed endpoint of 400 ordinary T=1 sweeps.
Independent fits may run concurrently.

Within each split and dimension configuration, the start winner is selected by
the largest eligible ordinary T=1 ELBO. Only after that table is frozen are the
held-out responses scored. Dimension selection compares paired holdout MSE /
NRMSE and an experimental moment-matched joint Gaussian log-score, each with a
one-standard-error rule that chooses the least complex eligible configuration.
Truth evaluation is a later, explicitly separated action.

This is a one-data-seed, two-split diagnostic. It can reject an unusable
criterion or justify a larger validation study, but it cannot establish broad
dimension-selection consistency.

See `PROTOCOL_G12_STAGE6F_JOINT_DIMENSION_SELECTION_V1_20260911.md` for the
frozen development protocol and the runner help for executable actions.
