# Stage 7F launch record (pending endpoint parity)

- ECS fit start marker: `2026-09-28 18:07:41 +08:00`, `/root/multifsync-results/stage7f_20260928_v1/FITS_STARTED.txt`.
- Three full replays started in parallel through `/root/multifsync-stage7f-20260928-v1/launch_stage7f_dense_replay_20260928_v1.sh`, supervisor PID `6413`, R parent `6414`, three fit workers `6431`–`6433` at the initial check. One CPU is requested per fit.
- `prepare` passed the Stage 7C source, winner, observation, and original-fit hash checks. The two-sweep trace smoke passed at `18:06:45 +08:00`.
- The first smoke attempt failed before any full fit because `result$anneal <- NULL` removed the named argument, allowing the API's default anneal. The new Stage 7F runner alone was corrected to retain the explicit `NULL` element and smoke then passed. No Stage 7C/7E file or package source was changed.
- Final uploaded runner SHA256: `eee93cb0689c2bf6ee79b118a237270b7ca94b6079d8b85020c8dcf95b160d50`. The other script/protocol hashes are in `/root/multifsync-results/stage7f_20260928_v1/COMPLETE_SOURCE_BINDING.sha256`.
- The launcher automatically runs the endpoint-parity gate after all three fits. A separate passive helper (PID `6711` at launch) runs the truth-blind trace tables only if `FITS_AND_ENDPOINT_PARITY_COMPLETE.txt` appears. It exits on fit or parity failure; no additional fit is launched.
- At this record, endpoint parity and trajectory attribution are **pending**. Check `endpoint_parity_v1/PARITY_PASS.txt` and `attribution_v1/ATTRIBUTION_TABLES_COMPLETE.txt` before drawing any mechanism conclusion. If `PARITY_FAILED.txt` appears, audit the recorded differences and stop interpretation.
- Stage 7C comparable runtimes suggest about 75–100 minutes of wall time for the three concurrent fits, plus parity and analysis. Expected fit completion is roughly `19:23–19:48 +08:00`, subject to ECS load.
