# G12 correct-basin quality audit Stage 5A protocol v1

Registration date: 2026-08-29

Status: authorized read-only post-truth development audit; not a formal v0L-V result

## Question

G12 reachability has already been independently checked on six new strong-identification data sets with 6--9 irregular observations and replicated on four additional 6--9-observation baseline/weak-separation data sets. Stage 5A does not repeat that experiment. It asks:

1. how heterogeneous scientific recovery is among endpoints that already have exact shared/specific identity;
2. whether the frozen truth-free maximum-ELBO winner tends to avoid scientifically poor exact-identity sub-basins;
3. whether scale, process, contribution or dense-grid pathologies recur under the original 6--9-observation design; and
4. whether strict practical convergence has any usable association with structure or scientific quality.

## Frozen inputs

Primary cohort:

- independent confirmation data `g12ic_01`--`g12ic_06`, seeds `82622001`--`82622006`;
- strong identification and 6--9 irregular observations;
- 12 G12 starts per data set, 72 endpoints total;
- the six original truth-free G12 winners remain frozen.

Secondary cohort:

- cross-scenario data `g12cs_base_01/02` and `g12cs_weak_01/02`, seeds `82623001`--`82623004`;
- 6--9 irregular observations;
- eight G12 starts per data set, 32 endpoints total;
- the four original truth-free G12 winners remain frozen.

The primary cohort has complete fit/evaluation/truth RDS and supports process, complete-contribution and loading-scale recomputation. The secondary cohort is an aggregate external replication and is not required to have local fit RDS.

## Read-only rules

- no data generation, fitting, continuation, annealing, pre-score, winner reselection or model update;
- no writes to either historical input directory;
- endpoint identity uses the frozen loading-first one-to-one matching semantics;
- factor counts or factor-PPI alone never define an exact endpoint;
- truth is used only for this post-truth audit; all winner identities were frozen before truth access;
- strict practical convergence is descriptive and never an eligibility or scientific-success filter;
- the two historical experiments remain separate cohorts and are not relabelled as one preregistered inferential experiment.

## Outputs

The audit produces:

- one 104-endpoint standard-density inventory and one 45-endpoint exact-basin table;
- frozen-winner quality, per-data reachability and within-exact metric ranks;
- within-data ELBO/scientific-metric rank associations;
- primary-cohort study-role process and complete-contribution recovery;
- primary-cohort raw/canonical loading-scale decomposition;
- convergence/structure summaries, input bindings, QC and a completion marker.

The primary analysis unit is the data set. Start-level results describe within-data endpoint heterogeneity and are not treated as independent Monte Carlo replicates.

## Decision interpretation

- Poor exact endpoints that are consistently avoided by the frozen winner are treated as managed nonconvex endpoint heterogeneity; they do not alone justify a model change.
- Recurrent poor frozen winners across independent 6--9-observation data sets, especially a joint observed/dense discrepancy, scale inflation and covariance/process/contribution degradation, would justify a separate model or regularization investigation.
- No threshold is tuned to truth in Stage 5A, and no new fit is automatically authorized by its result.
