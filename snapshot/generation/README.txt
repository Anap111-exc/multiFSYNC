Registered data generation for cavi_true222_fit444_m4_20261117_v1

Scope
Ten NEW independent populations, S=2, p=500, n_s=(60,60), no covariates.
True factors=(2,2,2), true FPCA rank=2 per factor, weights=(0.8,0.2).
Fitting capacities=(4,4,4), FPCA capacity=4 are recorded in the experiment
protocol and never supplied as truth-based starting information.
Each loading column has exactly 50 nonzero entries and squared norm45;
the four factor columns within each study have disjoint nonzero support.
Noise standard deviation0.3, 6-9 irregular training times, and three new
test times within each subject's training time range. The original NC
generator and all random-generation/calibration formulas are unchanged.

Registered seeds
config/DATASET_REGISTRATION.csv contains ten unique population records.
generation/SEED_REGISTRATION.csv contains sixty new population/test seeds
and120 fitting seed records. New seeds are deterministically generated from
the experiment namespace and checked for collisions. FIT_SEED_REFERENCE.csv
copies the previous original-CAVI twelve-start seed vector by rep_id;
no previous winner or recovery result enters the seed registration.

Cloud invocation, after the parent freezes PLAN and PROTOCOL
Rscript --vanilla GENERATION_ROOT/generate_registered.R EXEC_ROOT OUTPUT_ROOT BASELINE_LIBRARY
GENERATION_ROOT is EXEC_ROOT/generation.
BASELINE_LIBRARY=/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1
OUTPUT_ROOT=/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1
The script rejects Windows and any existing data/fit_inputs/generation marker.
No replacement seeds, partial-data continuation or automatic retry are enabled.

Required parent files
EXEC_ROOT/config/PLAN.csv, exactly120 original-arm rows, matching
FIT_SEED_REFERENCE.csv on dataset_id,start,seed; twelve starts per dataset.
EXEC_ROOT/config/PROTOCOL.json and PROTOCOL_FROZEN.json, where status=FROZEN
and protocol_sha256 binds PROTOCOL.json.
EXEC_ROOT/CODE_SHA256SUMS.txt binds every packaged script/config file.

Required generation runtime files
generate_registered.R
nc_common.R
stage7m_reconstruct.R
nc_generation_core.R
nc_fixed_metadata.R
project_fit_inputs.R
GENERATION_CONTRACT.json
FIT_SEED_REFERENCE.csv
SEED_REGISTRATION.csv
config/DATASET_REGISTRATION.csv
prepare_registration.py and test_local_contract.R are local preparation/tests;
they are not part of the cloud generation call.

Output contract
OUTPUT_ROOT/data/<data_id>/observation_bundle.rds
OUTPUT_ROOT/data/<data_id>/model_metadata.rds
OUTPUT_ROOT/data/<data_id>/test_design_blind.rds
OUTPUT_ROOT/data/<data_id>/sealed_truth/truth_bundle.rds
OUTPUT_ROOT/data/<data_id>/sealed_truth/test_truth.rds
OUTPUT_ROOT/fit_inputs/<data_id>/observation_bundle.rds
OUTPUT_ROOT/fit_inputs/<data_id>/model_metadata.rds
OUTPUT_ROOT/fit_inputs/<data_id>/INPUT_SHA256SUMS.txt
DATA_BINDINGS.csv retains data_id,scenario_id,role,relative_path,sha256,bytes
with roles observation,metadata,test_design,truth,test_truth.
DGP_QC.csv, DGP_LOADING_QC.csv and FIT_INPUT_PROJECTION_AUDIT.json expose QC.
GENERATION_RUNTIME.json records R/packages/system/BLAS and source/library binding.
GENERATION_STATUS.json reports PASS,new_datasets10,reused_datasets0,
registered_fits120; DATA_READY.txt is written only after full success.
startup/SOURCE_SHA256SUMS.txt and startup/GATE_PASS.txt bind generation QC only.
The separate startup/RUNNER_ACCEPTANCE.json must also pass before fitting;
generation QC never substitutes for zero-CAVI initialization equivalence.

Fit boundary
Only fit_inputs is used by fitting. Its dimensions and fit_metadata omit
true L_f,L_s,M_f,M_s; observed Y,Z,time_obs and Linux basis/knots/transforms,
C,C_g,time_g and crossproducts are exactly preserved. True parameters and
test responses remain only in separate sealed evaluation files. Generation
reads generating truth for legality/QC, never for initialisation or selection.

Local validation
All generation R files parse. LOCAL_ACCEPTANCE.json records15 passing
contract checks, including unchanged inherited function bodies, seed budget,
loading rescaling/noise preservation and removal of true dimensions from
fit inputs. No generator was run on Windows and no local fit was started.
Linux generation QC and runtime package/source gates remain mandatory.
