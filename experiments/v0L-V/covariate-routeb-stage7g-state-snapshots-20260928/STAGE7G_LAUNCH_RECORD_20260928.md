# Stage 7G launch record (results pending)

- Three registered Stage 7C C1 winner replays started on ECS at
  `2026-09-28 20:59:06 +08:00`, launcher PID `7827`, using three one-thread
  workers. No other fit is registered and there is no automatic rerun.
- Isolated Stage 7G package source SHA256:
  `b3e9050f1de7154baff28b7dd986254df0b634bf85d6f50f464d00dbd4c8f17f`.
  Comparison with the frozen Stage 7C archive found exactly one changed file:
  `R/multi_core.R`, containing only the private state-copy hook. The project
  package source under `Rcode/multiFSYNC/R/` was not edited.
- Source archive transfer SHA256 matched local and ECS; isolated `R CMD
  INSTALL` completed. Runner `selftest` passed with 79 legacy driver-trace
  sweeps and 500 raw state snapshots. `prepare` verified frozen winner,
  observation, Stage 7C fit, and Stage 7F fit hashes. A two-sweep smoke saved
  states 0, 1, 2 and exactly matched final score and support arrays. The
  smoke emitted only the expected max-iterations warning at the forced
  two-sweep stop; no full diagnostic fit had started at that point.
- A passive helper PID `8260` waits for the parity marker, then runs the
  truth-assisted scientific analysis and archives the results. It does not
  launch fits. Truth is not opened by the fit runner or parity code. All
  source hashes are in ECS output `COMPLETE_SOURCE_BINDING.sha256`.
- Early status check: `FITS_STARTED.txt` present; `complete=0`, `failed=0`,
  parity pending. Stage 7F's comparable concurrent replays took 46-50
  minutes. Because Stage 7G also retains all 500 raw states, provisionally
  allow about 60-90 minutes of wall time from launch, plus archive transfer.
  This is a planning estimate, not an observed completion time.
- Decision gate: first inspect `endpoint_parity_v1/PARITY_PASS.txt`; only then
  interpret `analysis_v1/ANALYSIS_COMPLETE.txt`. If fit or parity fails, the
  passive helper stops and retains all partial outputs for diagnosis.
