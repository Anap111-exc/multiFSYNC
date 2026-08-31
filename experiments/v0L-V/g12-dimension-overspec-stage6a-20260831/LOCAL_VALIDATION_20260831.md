# Stage 6A local validation record

Date: 2026-08-31

Scope: nonformal Windows micro smoke; no full-size Stage 6A data or fit

Formal v0L-V or paper Monte Carlo result: `FALSE`

Validation results:

- all five Stage 6A R files parsed successfully;
- all three server shell scripts passed `bash -n` under Git Bash;
- manifest QC passed `14/14` with one data row, three fitted-dimension
  configurations, 36 fit rows, 12 starts per configuration, three registered
  capacity probes and three planned within-configuration winners;
- the targeted contract test passed dimension parsing for `(L=1,M=2)`,
  `(L=3,M=2)` and `(L=1,M=4)`, plus truth-free within-configuration ELBO
  selection;
- `SOURCE_BINDING.csv` contains 54 unique runtime/protocol files and matched
  `54/54` before the external validation copy was made;
- because the local R installation cannot reliably open the authoritative
  Chinese-character path, the unchanged worktree was copied byte-for-byte to an
  external ASCII-only temporary directory for package installation and smoke;
- a fresh install of multiFSYNC `0.3.0.9000` succeeded from that copy;
- the runner check reported `data=1`, `fits=36`, `workers=3`;
- the micro smoke generated and sealed one micro data set, fitted all three
  dimension configurations, froze three truth-free selections, authorized truth
  only inside the nonformal smoke, evaluated the three endpoints and summarized
  resources;
- fit endpoints completed and were objective eligible `3/3`; the only fit
  warnings were the expected six-sweep smoke-budget warnings;
- evaluation completed `3/3` with zero evaluation warnings/errors and each fit
  produced all 19 registered metric tables, including all-candidate activity and
  FPCA-rank activity;
- `SMOKE_COMPLETE.txt` reported `status=PASS`, `fits=3`, `g12_winners=3`,
  `pooled_endpoints=0`, and truth-free selection before smoke evaluation;
- smoke output stayed outside Git and is not a scientific result.

The multiFSYNC package source was unchanged from parent commit
`bede48a77993de85cec66604b98209cddf54e61b`; therefore the full package testthat
suite was not repeated. This validation did not generate the registered
full-size data, start any of the 36 registered fits, unseal registered truth,
start continuation, or create a formal result.
