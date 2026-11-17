Independent first-endpoint resource observation; no frozen-source modification.

Deploy this amendment outside the frozen execution source tree. Pass both
the independently verified observer SHA256 and original CODE_SHA256SUMS hash:
python3 -B observe_first_endpoint.py --observer-source-sha256 OBSERVER_SHA256 --source-manifest-sha256 FROZEN_MANIFEST_SHA256

Run only after STAGE/PIPELINE_LAUNCH.json records STARTED with the 120-fit
experiment identity. The observer checks the exact pipeline PID, process group
and start ticks. It never writes the registry, scheduler ledger, fits, frozen
sources or lifecycle state, and never sends signals or invokes cloud APIs.

It checks for complete successful endpoints every30seconds. An inotify watch
on OUTPUT_ROOT wakes it immediately when FINAL_STATUS is created; hashing
also checks that marker once per1MiB chunk. Without inotify it fails closed
and does not stay resident. If FINAL_STATUS exists or the registered pipeline
exits/changes identity, it writes a termination report and exits. Event wakeup
minimizes but cannot mathematically eliminate a race with lifecycle checks;
this observer must not be mislabeled as a system service. Confirm its exit
after the first report. The report is not required by the frozen pipeline.

Report: OUT/resource_monitor/FIRST_ENDPOINT_STORAGE_REASSESSMENT.json
Status REASSESSED means a report was computed, not that resources are proven
sufficient; see forecast.status. An insufficient forecast only reports risk.
No automatic cancellation, cleanup, admission, restart or shutdown follows.

A complete endpoint requires scheduler EXIT SUCCESS, STATUS complete,400
ordinary/99annealing updates, source-manifest binding, and full/compact file
SHA256 and sizes matching producer roundtrip receipts. Python verifies file
integrity and the producer attestation, not R object identity independently.
Only verified complete directories are used for per-candidate size; partial
files are listed and excluded. Existing partial bytes already reduce free
disk. Remaining estimate uses max(probe allowance,1.25xlargest complete
candidate+worker logs), all120 minus completed candidates,512MiB science/
controller allowance,8 extra candidate slots for in-flight writes,3GiB reserve.
Failed fits counted in the conservative remainder do not authorize retries.

Worker peak RSS comes from its durable EXIT record; two-second sampling can
miss transient peaks. No hard full-runtime memory/storage bound is claimed.
Frozen source files and their manifest are rehashed before and after the
assessment. The observer manifest and test acceptance are separate artifacts.

Tests use synthetic files and fake process identities only. No ECS, fitting,
real process signals or changes to any frozen source are used by the tests.
