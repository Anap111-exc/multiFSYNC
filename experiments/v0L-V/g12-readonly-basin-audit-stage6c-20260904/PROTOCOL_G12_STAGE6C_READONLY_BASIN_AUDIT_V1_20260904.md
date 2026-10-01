# Stage 6C G12 read-only basin audit protocol v1

Date: 2026-09-04

## Purpose

Stage 6C asks whether independently seeded G12 endpoints form reproducible
truth-free basins, whether the maximum-ordinary-ELBO endpoint belongs to a
supported basin, and whether selecting the highest-ELBO representative of a
modal basin improves scientific recovery after the truth-free result is
frozen.

This is a retrospective development audit. It does not fit, continue, refit,
or modify any model.

## Inputs

The truth-free phase reads only `fit.rds` and `terminal_record.rds` from:

- Stage 6A: 36 endpoints, three fitted upper-bound configurations and 12
  starts per configuration;
- Stage 6B v2: 20 endpoints, five fitted upper-bound configurations and four
  starts per configuration.

Stage 6A full-data endpoints and Stage 6B holdout-training endpoints are never
pooled. Each experiment cohort and fitted upper-bound configuration is a
separate audit stratum.

The truth-free phase must not read observations, sealed truth, unsealed truth,
or any `evaluation/` file. Every terminal record must be objective eligible,
must state that truth was not used for fitting/stopping/selection, and must
state that continuation was not used.

## Basin objects

For role `r` in shared `A`, study-1-specific `B1`, and study-2-specific `B2`,
write the fitted complete contribution as

`C_r = A_r Z_r`,

where `A_r` contains fitted loadings and `Z_r` contains the fitted subject
factor paths on the common dense grid. Pairwise complete-contribution cosine
is calculated without materialising the large matrix:

`<C_1,C_2>_F = sum[(A_1' A_2) * (Z_1 Z_2')]`.

This object is invariant to matched sign, permutation, scale and within-role
rotation transformations that leave the complete fitted contribution
unchanged. Loading-subspace similarity is the mean squared canonical
correlation, with zero directions retained under the frozen convention.

Two clustering views are retained:

1. `contribution`: minimum complete-contribution cosine over `A/B1/B2`;
2. `structural`: the minimum of the three contribution cosines and the three
   loading-subspace similarities.

Complete-linkage clusters are reported at fixed similarities 0.80, 0.90 and
0.95. Similarity 0.90 is the primary descriptive threshold because it was
already used in the earlier Stage 6B-A recurrence audit; no threshold is
tuned against truth.

## Truth-free candidate rule under audit

Within each stratum:

1. retain objective-eligible endpoints;
2. find the largest cluster at a fixed clustering view and threshold;
3. if several clusters tie for largest, choose the tied cluster whose
   representative has the largest ordinary ELBO;
4. choose the maximum-ordinary-ELBO endpoint within that cluster.

This is called `modal_basin_then_elbo`. It is compared with the historical
`raw_max_elbo` rule only after the truth-free files have been written and
hashed.

## Post-freeze evaluation

The post-freeze phase may read only existing Stage 6A and Stage 6B evaluation
outputs. It reports core-direction recovery, exact factor structure, true
effective FPCA rank, dense signal NRMSE, feature ISE and covariance-operator
error. Truth is used only to evaluate the already frozen basin assignments and
candidate rules.

## Interpretation boundary

This audit contains one development data seed. Stage 6A and Stage 6B differ in
their training observations, budgets and starts per stratum. Results may
support or reject a candidate basin rule, but cannot establish a general
dimension-selection procedure.

