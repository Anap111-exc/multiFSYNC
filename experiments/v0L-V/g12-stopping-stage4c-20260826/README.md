# G12 stopping semantics stage 4C

This directory contains the registered, truth-free, independent fixed-400
development validation for the public G12 stopping profile.

Run order on the Linux ECS server:

1. build and verify `DATA_MANIFEST.csv`, `FIT_MANIFEST.csv`, and
   `SOURCE_BINDING.csv`;
2. install the exact package source and run the non-formal micro smoke;
3. set `G12SS4C_RUNNER_COMMIT` to the exact runner commit and execute
   `bash launch_stage4c_20260826_v1.sh`;
4. monitor with `bash status_stage4c_20260826_v1.sh`.

The detached runner performs six sealed data generations, 12 fixed-400 fits,
truth-free local-stability analysis, and six ordinary-ELBO selections. It does
not read sealed truth, run continuation, or start an 800-sweep extension.
