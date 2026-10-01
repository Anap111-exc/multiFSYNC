# Stage 6G local micro-smoke result (2026-09-11)

Status: **PASS**. This was a development-only interface smoke. It is not a
Stage 6G scientific result and it does not authorize the registered 120
endpoints.

## Scope

- Windows R 4.5.1; package source copied byte-for-byte to an ASCII temporary
  directory and installed into a temporary library to bypass the known
  `pkgload`/Chinese-path `DESCRIPTION` discovery failure.
- Three micro data sets (`n_s = c(4,4)`, `p = 20`, 6--7 irregular times) were
  generated with `truth_basis_method = "jaoua_splinets"` for the three Stage 6G
  truth dimensions.
- Both registered whole-time holdouts were constructed. The upper held-out
  index was exactly one position after the middle held-out index for every
  subject, and each training curve retained at least five times.
- All five candidate dimension configurations were fitted once on each split
  (10 micro fits total) with G12, one pre-score sweep, short test annealing,
  six ordinary sweeps, `lambda_orth = 0`, and `n_cpus = 1`.
- Held-out RMSE/NRMSE and the five-configuration paired one-SE code path were
  exercised without truth access.

## Checks

| Check | Result |
|---|---:|
| Protocol validator | 17/17 PASS |
| Registered full endpoints (not run) | 120 |
| Micro truth scenarios generated | 3/3 PASS |
| Continuous Gram maximum error | `2.00e-15` |
| Response identity maximum error | `2.64e-16` |
| Holdout split contracts | 2/2 PASS |
| Candidate dimension micro fits | 10/10 PASS |
| Paired five-configuration one-SE path | PASS |
| Truth used by fit or holdout scoring | FALSE |
| Formal experiment started | FALSE |

All 10 fits emitted the expected `Max iterations reached before convergence.`
warning because the smoke deliberately stops after six ordinary sweeps. The
warnings were captured rather than discarded. They are not a failure of the
full 400-sweep protocol.

The micro one-SE selection happened to be `cand_L111_M2`, while the minimum
micro loss was `cand_L211_M2`. This has **no scientific interpretation**:
there were only eight subjects, one start per configuration/split, reduced
`p`, and six sweeps. It verifies only the data flow and deterministic
simplification rule.

Elapsed time for generation plus the 10 micro fits was 4.75 seconds. No ECS,
formal data, full-size fit, continuation, or truth-unsealing action was used.
