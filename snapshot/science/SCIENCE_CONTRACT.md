# True (2,2,2), fitted (4,4,4), FPCA 2 to 4: scientific evaluation

This is one original-CAVI arm, ten independently generated populations/datasets,
twelve registered starts per dataset. No MGP contrast is fabricated. The unit of
scientific replication is a dataset, not a start.

## Selection and data separation

All twelve registered attempts must have terminal receipts. Only complete,
finite, objective-valid endpoints passing frozen metadata, initial-state/RNG,
SHA256 and monotonicity gates enter maximum ordinary ELBO selection. Ties use
smallest start, then fit ID. Each choice and its artifacts are sealed before
reading truth or held-out responses. Prediction is computed before truth.
Scientific evaluation makes zero posterior updates and does not replace the
selected endpoint. Invalid attempts remain in the registered denominator of 120;
an unavailable winner is not silently dropped from the denominator of ten.

## Factor and FPCA counts

Formal factor counts use factor PPI >= 0.5 in each block. Truth is two shared
factors and two specific factors in each study: six distinct global true roles,
eight study-role records because each shared factor appears in both studies.
Formal NC structure retains the original one-to-one global/study loading catalog
criterion at absolute cosine 0.8, exact counts, and no extra/duplicate direction.
It is a loading-structure criterion, not a claim of perfect trajectory recovery.

FPCA diagnostics map all true factors within their own block one-to-one using
retained factors and loading cosine >= 0.9. The objective is maximum match count,
then total cosine, then the frozen lexical tie rule. Unlike the old one-factor
adapter, two truths cannot share one estimated factor. The six roles remain in
every dataset; missing roles and missing winners are explicit categories. Reported
FPCA ranks are mapped back through the factor permutation. They are postprocessed
effective ranks, not the fitted capacity of four. True function components 1 and 2
are aligned to up to four estimated components using the frozen component criterion.

The >= count rule is formal. The old > 0.5 factor count remains separately labelled
as a diagnostic; loading-support fields explicitly labelled `ppi_gt_half` preserve
their historical > threshold. No threshold is calibrated from this experiment.

## Scientific metrics

Preserve NC/Stage7M definitions, original-scale raw estimates, sign/permutation
matching, and no regression rescaling to truth. Evaluate observed times and the
201-point continuous grid; sum every fitted shared and specific factor for complete
contributions. Report shared, specific, combined factor, mean, signal errors;
per-function ISE/cosine/eigenvalues/PVE; per-component score NRMSE/correlation;
raw and FPCA-canonical loading errors, AUC/AP/support metrics; MTA, MTB1/MTB2 and
MTB using loading-subspace trace coverage; centered reconstruction, noise variance
and noise SD, pointwise variational moment intervals, and reserved-time prediction
for existing subjects. Covariate metrics are not applicable because generated and
fitted d=0. Original RDS and point-level tables preserve detailed outputs.

Unique-factor component/loading summaries give each of the six true global roles
equal weight; shared roles are not counted twice merely because there are two
studies. Complete raw contribution summaries use equal-study, then equal-subject,
then equal-variable weighting. Scientific summaries include finite N, mean, SD,
median and MCSE; count summaries include mean, bias, MAE, RMSE, under/exact/over and
missing-winner totals. Missing continuous errors remain unavailable, never zero.

## Entrypoints

* `Rscript science/evaluate_dataset.R --exec-root EX --output-root OUT --data-id ID`
* `python science/summarize.py OUT`
* `Rscript tests/science_oracle.R EX OUT/science_tests/oracle`
* `python tests/science_summary.py EX OUT/science_tests/summary`
* `Rscript science/evaluate_timing_endpoint.R --exec-root EX --output-root OUT --data-id ID --fit-path PROBE_FIT`

All production R processes require one BLAS/OpenMP thread. Frozen scientific helper
hashes are verified at load time. Package manifests bind modified adapters too.
Per-dataset completion remains `science/ID/SCIENCE_COMPLETE.json`; aggregate
completion remains `aggregate/SUMMARY.json` and
`aggregate/SCIENTIFIC_EVALUATION_COMPLETE.json`.

The small oracle is an analytic fixture, not generated experimental data or a fit.
It checks perfect recovery, all-factor observed/grid contributions, permutation
and coupled sign/scale invariance, missing/extra factors, shared/specific confusion,
PPI=0.5, M4 component matching, and selection failure handling. The summary test
checks all/one/no available winners with six-role denominators. The cloud smoke
uses the existing two-sweep resource probe solely to check full scientific schema
and finite required metrics; it has no recovery threshold and is not scientific
evidence about estimator performance. No automatic retry or endpoint overwrite.
