# Stage 7C pre-unseal truth-free freeze

- Experiment: `covariate_routeb_stage7c_20260927_v1`
- Formal endpoints: 108/108 `COMPLETE`, 108/108 objective-eligible, 0 errors.
- Pools: 3 logical data seeds x 3 arms x 12 starts.
- Truth-free winners: 9/9, one per data-arm pool.
- Selection rule: maximum eligible ordinary T=1 ELBO within each data-arm pool; lexical `fit_id` tie-break.
- Cross-arm ELBO comparison: forbidden and not performed.
- Truth used for fitting, stopping, eligibility, or selection: no.
- Frozen winner file SHA256: `1e1e7791b3e975287ab195c87a801ab0fa24f97d26558a85fccc3b11af70a877`.
- Supervisor exit: 0; automatic stop request was submitted successfully.
- Scientific evaluation status at capture: not started.

The 12 machine artifacts listed in `PRE_UNSEAL_FILE_SHA256.csv` were downloaded
from the ECS result directory and hashed locally before the sealed truth was
opened. This directory is evidence of the pre-unseal state and must not be
rewritten by scientific evaluation.
