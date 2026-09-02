# Stage 6B-A read-only audit and Stage 6B-B holdout interface report

Date: 2026-09-01

Status: development evidence; not a formal paper Monte Carlo result

## 1. What was executed

Stage 6B-A read the 36 existing Stage 6A `fit.rds` files, their 36
`terminal_record.rds` files and the three-row truth-free winner table. The input
allowlist excluded the observation bundle, sealed/unsealed truth and all Stage
6A evaluation outputs. No fit or continuation was started. The successful run
took 4.73 seconds and used 97,916 KiB peak resident memory on the Linux server.

The first attempt stopped before writing results because its integrity test
incorrectly expected a zero loading direction to have self-cosine one. The
frozen convention is instead that any cosine involving a zero direction equals
zero and remains in one-to-one matching. The integrity test alone was corrected
to check that similarities are finite and in `[0,1]`; no metric or Stage 6A
input was changed. The rerun completed with all seven integrity checks passing.

Stage 6B-B then implemented and tested an experiment-local whole-time-point
holdout interface. The only new fit was a nonformal `p=20`, six-sweep software
smoke; it is not scientific evidence.

## 2. Audit coverage and integrity

- fits: 36/36;
- truth-free winners: 3/3;
- candidate factor rows: 180;
- inputs: only fit endpoints, terminal records and truth-free winner table;
- terminal truth-use flags: all `FALSE`;
- zero loading directions retained in the audit: 43/180;
- per-fit contribution shares: sum to one for every fit;
- truth read: no;
- new fit or continuation: no.

The 180 candidate rows comprise 108 directions from `factor_L3_M2`, 36 from
`fpca_L1_M4`, and 36 from `truth_L1_M2`.

## 3. Factor-count over-specification (`factor_L3_M2`)

The over-specified configuration contains nine candidate factors per fit. Its
truth-free summaries show three different phenomena.

First, only a small core of loading directions is reproducible. Against the
maximum-ELBO winner, the leading shared, study-1-specific and
study-2-specific directions recur at cosine at least 0.9 in 12/12, 12/12 and
11/12 starts. The other reference directions recur much less often: shared
directions 2 and 3 recur in 6/12 and 2/12 starts; the two extra study-1
directions recur in 0/12; the two extra study-2 directions recur in 2/12 and
0/12. Correspondingly, the full three-dimensional loading subspaces are not
stable: median mean squared canonical correlations are 0.344, 0.500 and 0.503
for shared, study 1 and study 2, and median minimum canonical correlations are
0, 0 and 0.064.

Second, duplicate and misplaced directions are common. At cosine at least 0.9,
within-role duplication occurs in 6/12 shared-role fits, 2/12 study-1-role fits
and 4/12 study-2-role fits. Shared-specific leakage occurs in 6/12 fits for each
study. The maximum-ELBO winner itself has a shared within-role duplicate
(cosine 0.998) and shared-study-1 leakage (cosine 0.996).

Third, factor PPI and contribution energy do not by themselves solve the
problem. Across all 108 candidate rows, 39 are exact zero loading directions,
while 69 pass factor PPI 0.5. Six PPI-retained rows nevertheless contribute at
most 1% of fitted energy. In the winner, six of nine candidates pass PPI 0.5:
three shared candidates contribute 68.5%, 4.3% and 2.7%; the leading specific
directions contribute about 11.9% in each study; and another study-2 candidate
has PPI 1 but contributes only 0.58%. Thus PPI removes exact inactive blocks,
but can retain weak, duplicated or misplaced nonzero blocks. A contribution
cutoff alone is also unsafe because additional directions sometimes carry
non-negligible energy.

## 4. FPCA over-specification (`fpca_L1_M4`)

All 36 fitted factors pass factor PPI 0.5. The 99% cumulative-PVE rule retains
all four components for 29/36 factors and three components for 7/36; it never
returns the generating rank two. This independently reproduces the Stage 6A
finding that the current 99% PVE description is not an effective `M` selector
in this setting.

Because there is only one factor per role, there is no within-role factor
duplication, and no cross-role cosine reaches 0.9. However, loading recurrence
against the winner is only 8/12 in each role, showing that the dimension being
fixed does not remove start-to-start basin variation. Median matched FPCA
subspace mean squared canonical correlations are 0.906, 0.752 and 0.684 for the
shared, study-1 and study-2 factors. These values demonstrate measurable
subspace recurrence, but do not justify retaining four components: recurrence
is not incremental predictive utility.

## 5. Correct fitted dimension (`truth_L1_M2`) as a stability warning

Even at `(1,1,1)` factors and `M=2`, four of 36 candidate directions are zero.
Winner-referenced loading recurrence at cosine at least 0.9 is 9/12 for the
shared direction, 8/12 for study 1 and only 5/12 for study 2. There is no
within-role duplication; one non-winner start shows study-2 cross-role leakage.
The maximum-ELBO winner itself has negligible cross-role cosines.

This is the key limitation of a stability selector: a scientifically relevant
direction can be a minority basin. Therefore cross-start stability is a useful
uncertainty and redundancy audit, but a hard rule such as “retain only
directions recurring in more than half of starts” could delete a needed
study-specific direction.

## 6. Stage 6B-A conclusion

The audit successfully measures the intended truth-free diagnostics and detects
the known over-specification failure without consulting truth:

- factor over-specification produces unstable secondary directions, frequent
  duplicates and shared-specific leakage;
- factor PPI rejects exact zeros but does not reliably distinguish weak,
  duplicate and misplaced nonzero factors;
- complete-contribution energy adds useful scale-invariant activity evidence,
  but has no defensible universal cutoff here;
- 99% PVE does not reduce `M=4` to an effective rank near two;
- stability cannot be the primary selector because correct-dimension directions
  can occupy a minority basin.

Accordingly, PPI, energy, PVE and stability should remain an audit panel. A
genuine held-out prediction score is still needed to compare dimensions.

## 7. Stage 6B-B interface semantics

The new interface:

1. removes exactly one interior observation time from every subject, including
   the complete `p`-variable vector at that time;
2. requires at least five remaining time points and rejects truth-named input
   fields;
3. fits only the remaining observations;
4. predicts the held-out vector using interpolated fitted mean/eigenfunctions,
   fitted shared and study-specific loadings, and that subject's scores inferred
   from the retained times;
5. averages squared and absolute loss first over variables within a
   subject-time vector and then over subjects;
6. supports an all-factor primary reconstruction and a factor-PPI `>=0.5`
   sensitivity reconstruction.

Version 1 intentionally supports only `d=0`. It returns held-out RMSE, NRMSE
and MAE. Its NLPD uses `sigsq_eps_hat` only and omits uncertainty in means,
loadings, eigenfunctions and subject scores. It is therefore labelled
`residual_noise_only_not_calibrated_full_posterior`, and
`primary_dimension_score_ready=FALSE`. It must not be presented as calibrated
posterior predictive density.

## 8. Validation results

All five Stage 6B R files parse. Targeted tests passed 9/9, covering:

- sign-invariant one-to-one matching and the zero-loading convention;
- loading and weighted-FPCA subspace invariance;
- complete middle-time-vector removal;
- rejection of truth-contaminated observations;
- the actual Stage 6A observation-only bundle: 60 held-out subject-time vectors,
  `p=500`, and 5--8 retained times, without fitting;
- exact posterior-mean prediction for shared/specific factors under all-factor
  and PPI modes;
- point-loss aggregation and the uncalibrated-NLPD flag;
- static exclusion of truth and evaluation paths from the read-only audit.

The real-code micro smoke passed on both Windows and the target Alibaba Cloud
Linux server with `S=2`, `n_s=(4,4)`, `p=20`, `d=0`, one held-out time per
subject, one CPU, G12 Gram-unit-energy initialization, one pre-score, annealing
and six short sweeps. Both platforms produced the same eight finite
subject-time prediction records and identical aggregate scores. The Linux run
took 2.69 seconds and used 153,380 KiB peak resident memory. The single warning,
“Max iterations reached before convergence,” is expected for the deliberately
six-sweep software smoke. Its held-out NRMSE of 0.817 is not a scientific
result.

The package source under `R/`, `NAMESPACE` and `DESCRIPTION` was not changed, so
the full package testthat suite was not repeated.

## 9. Remaining boundary

Stage 6B-B is ready for point-prediction software use, not yet for a registered
dimension-selection claim. Before a predictive pilot, the protocol must choose
one of two explicit paths:

- pre-register held-out RMSE/NRMSE as the pilot's primary point-prediction loss
  and leave residual-only NLPD diagnostic; or
- first implement and calibrate the full posterior predictive variance, then
  use NLPD as primary.

No full-size Stage 6B data, fit, continuation, truth unseal or predictive
dimension selection was started in this work.
