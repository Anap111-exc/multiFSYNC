# Stage 4E: focused G12 reachability audit

This directory contains the versioned runner and manifests for an adaptive,
development-only Stage 4E audit. It reuses three Stage-4C data sets, adds six
G12 starts per data set, freezes an eight-start truth-free ELBO winner, and only
then evaluates all endpoints against truth already unsealed in Stage 4D.

The package source remains pinned to commit
`72d9a53f0ef5e9d3e5f49d9cf837958207b3c980`. No model, prior, CAVI, ELBO,
initialization route, annealing schedule, or fixed-400 stopping profile is
changed. Results are not formal v0L-V results and must not be described as an
unseen validation.

Server execution root:
`/root/v0lv-g12-focused-reachability-stage4e-20260827-v1`.

The long runner performs, in order: 18 new truth-isolated fits, truth-free
selection across 24 endpoints, then authorized truth evaluation and summary.
Any fit error is retained and stops the pipeline without selective rerunning.
