# G12 Stage 6B truth-free predictive dimension-selection pilot

Protocol ID: `G12_STAGE6B_PREDICTIVE_PILOT_V1_20260902`

Registration date: 2026-09-02

Status: registered development pilot; execution authorized by the user on
2026-09-02; truth access, continuation and final full-data refitting are not
authorized

## 1. Question and evidence boundary

This pilot asks whether genuine held-out prediction can distinguish the
dimension-overspecified configurations that ordinary training ELBO did not
distinguish safely in Stage 6A. It does not modify the model, priors, CAVI,
ordinary ELBO, G12 initialization, pre-score, Jaoua annealing or the fixed
400-sweep development endpoint.

The pilot reuses only the observation-only bundle for Stage 6A data set
`g12s6a_01` (`data_seed=83131001`). No new data are generated. This data set's
truth was examined in the earlier Stage 6A development analysis, so this pilot
can test software operation and signal detection on a known failure mode but
cannot independently validate dimension selection. The runner must not read a
truth bundle or any Stage 6A evaluation result.

## 2. Registered split

For every subject, remove the complete `p=500` response vector at the
deterministic interior index

`floor((number_of_observed_times + 1) / 2)`.

The identical split is used for every configuration and start. Each subject
retains 5--8 irregular training times. The 60 held-out subject-time vectors are
the scoring units; the 500 variables within a vector are averaged and are not
treated as independent replicates.

## 3. Registered route and compute budget

Every fit uses:

- random G12 `gram_unit_energy` initialization;
- one pre-score sweep;
- Jaoua annealing `c(1,1.9,100)`;
- a fixed endpoint of 400 ordinary `T=1` sweeps;
- `lambda_orth=0`;
- `n_cpus=1`;
- no continuation, racing or automatic 800-sweep extension.

Four paired start seeds are fixed in `START_SEEDS_STAGE6B_PILOT_V1.csv` and are
reused across configurations. Three independent fits may run concurrently on
the four-vCPU server.

## 4. Sequential candidate grid

Phase 1 screens factor count with `M=2`:

- `(L_f,L_s[1],L_s[2])=(1,1,1)`;
- `(2,2,2)`;
- `(3,3,3)`.

This is 12 fits. Within each fixed configuration, the winner is selected before
held-out scores are inspected by maximum finite, objective-eligible ordinary
`T=1` ELBO; ties use lexicographic `fit_id`.

Phase 2 screens `M=2,3,4` only at the factor count selected in Phase 1. The
Phase-1 `M=2` fits are reused and eight conditional `M=3/4` fits are added.
Thus at most 20 fits are executed. The complete conditional grid is registered
in `FIT_CONFIGS_STAGE6B_PILOT_V1.csv`; configurations inconsistent with the
frozen Phase-1 factor count are not run.

## 5. Primary score and one-standard-error rule

The preregistered primary score is held-out signal NRMSE computed against the
actually withheld observations, not against latent signal truth. Held-out RMSE
is a mandatory companion. Because every configuration uses exactly the same
held-out responses, ranking by NRMSE, RMSE and mean subject-time MSE is
equivalent.

The one-standard-error rule is evaluated on mean subject-time MSE:

1. average squared residuals over the 500 variables within each held-out
   subject-time vector;
2. average those 60 losses across subjects;
3. estimate the standard error of that mean with studies as fixed strata,
   `sqrt(sum_s (n_s/N)^2 * var_s(loss)/n_s)`;
4. find the minimum-loss configuration;
5. choose the simplest configuration whose mean MSE is no greater than the
   minimum plus the minimum configuration's standard error.

Phase-1 simplicity is increasing common factor count `L=1,2,3`. Phase-2
simplicity is increasing common FPCA cap `M=2,3,4`. NRMSE and RMSE are reported
for interpretation, while the one-standard-error threshold remains on the MSE
scale.

All-factor posterior-mean prediction is primary. Prediction after retaining
only factors with factor PPI at least 0.5 is a sensitivity analysis. The
residual-noise-only NLPD is uncalibrated and is not eligible for selection.

## 6. Execution and failure rules

- Objective eligibility is separate from practical convergence.
- Strict/practical convergence is reported but is not a scientific-quality
  selector.
- Every registered fit has an independent directory containing the fit RDS,
  ELBO, warnings, runtime, terminal record and `FIT_COMPLETE` or `FIT_ERROR`.
- A failed or objective-ineligible registered fit is retained and stops the
  primary pilot; it is not selectively rerun.
- Within-configuration winner selection must be frozen from terminal objective
  records before held-out losses are read for cross-dimension selection.
- No latent truth, structural labels, dense-grid error, ISE or R/P/L may enter
  fit, stopping, start selection or dimension selection.

## 7. Endpoint and interpretation

The pilot stops after a truth-free factor-count selection and a truth-free FPCA
cap selection have been written and hashed. It does not unseal truth, start a
continuation or refit the selected dimension to all observations.

Support for advancing the rule requires the one-standard-error selector to
prefer a simpler configuration without a material held-out prediction penalty.
Selection of a larger configuration is retained as a negative result. Either
outcome remains single-development-data evidence and requires later independent
data seeds before a paper-level claim.

Based on Stage 6A timings, Phase 1 is expected to take about 8--11 wall-clock
hours with three workers. Total time is expected to be about 11--30 hours,
depending on the factor count selected before the `M=3/4` phase.
