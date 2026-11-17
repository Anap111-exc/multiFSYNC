"""One durable, completion-bound launch after all registered acceptance checks."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
EX=Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT=Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
STAGE=Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
EXPERIMENT=OUT.name
assert sys.platform.startswith('linux') and all(p.resolve()==p for p in (EX,OUT,STAGE))
sys.path.insert(0,str(EX))
from lifecycle import finalize as life
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
ops=Path(__file__).resolve().parent
for row in read(ops/'AMENDMENT.json')['files']:
    assert sha(ops/row['relative_path'])==row['sha256'],'Operational amendment hash differs'
target=STAGE/'PIPELINE_LAUNCH.json'
assert not target.exists() and not (OUT/'pipeline').exists() and not (OUT/'scheduler').exists(),'Prior attempt exists; no implicit restart'
baseline=STAGE/'PROCESS_BASELINE.json';review=read(baseline)
assert review['experiment_id']==EXPERIMENT and review['instance_id']=='i-bp1egmsairexm9w1rnks' and review['region']=='cn-hangzhou'
assert all(review[k] is True for k in ('all_processes_reviewed','no_other_tasks_require_instance','billing_mode_supports_stopcharging'))
assert 0<=time.time()-review['reviewed_unix']<=900
baseline_sha=sha(baseline);assert baseline_sha==(STAGE/'PROCESS_BASELINE.sha256').read_text().strip()
runtime=read(STAGE/'RUNTIME_ACCEPTANCE_V3.json')
assert runtime['status']=='PASS' and all(runtime[k] is True for k in ['stopcharging_permission_dryrun','startup_original_gate','science_tests','registration_checked'])
prepared=read(STAGE/'PREPARE_ACCEPTANCE.json')
assert prepared['status']=='PASS' and prepared['new_datasets']==10 and prepared['registered_fits']==120
assert prepared['projected_input_files_verified']==20 and prepared['formal_fits']==0
assert sha(OUT/'GENERATION_STATUS.json')==prepared['generation_status_sha256']
config=read(Path(__file__).resolve().parent/'EXECUTION_CONFIG.json');assert config['status']=='FROZEN' and config['new_fit_count']==120
assert config['arms']==['original'] and config['fitted_factor_capacities']==[4,4,4] and config['fitted_fpca_capacity']==4
resource=read(STAGE/'RESOURCE_PROBE_ACCEPTANCE_V3.json')
assert resource['status']=='PASS' and resource['registered_fit_attempts']==0 and resource['independent_of_registered_pool'] is True
assert resource['predicted_new_artifact_bytes']<=config['planned_new_artifact_bytes']
dep_path=STAGE/'RUNTIME_LIBRARY_BINDING.json';assert sha(dep_path)==runtime['runtime_library_binding_sha256']
assert sha(dep_path)==sha(EX/'RUNTIME_LIBRARY_BINDING.json'),'Reviewed and executed dependency bindings differ'
dep=read(dep_path)
for row in dep['package_files']:
    p=Path(dep['package_root'])/row['relative_path'];assert p.is_file() and sha(p)==row['sha256']
def check(name):
    v=life.process_check(life.process_inventory(),review,{'owned_process_groups':[]},os.getpid())
    life.persist(STAGE/name,v,exclusive=True);assert not v['blocking_processes'],'Other or unreviewed processes block launch'
check('PREVALIDATE_PROCESS_CHECK.json')
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',FSYNC_EXEC_ROOT=str(EX),FSYNC_OUTPUT_ROOT=str(OUT),FSYNC_PROCESS_BASELINE_SHA256=baseline_sha,R_LIBS=dep['R_LIBS'],R_LIBS_USER=dep['R_LIBS_USER'])
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
cmd=['python3','-B',str(Path(__file__).resolve().parent/'pipeline_v3.py')]
v=subprocess.run(cmd+['--validate-only'],cwd=str(EX),env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
life.persist(STAGE/'PIPELINE_VALIDATE_ONLY.json',{'exit_code':v.returncode,'stdout':v.stdout.decode(errors='replace'),'stderr':v.stderr.decode(errors='replace')},exclusive=True)
assert v.returncode==0,'Pipeline validation failed'
check('PRELAUNCH_PROCESS_CHECK.json')
life.persist(target,{'status':'PENDING','time':time.time(),'automatic_retries':0},exclusive=True)
with (STAGE/'pipeline.stdout.log').open('xb') as so,(STAGE/'pipeline.stderr.log').open('xb') as se:
    child=subprocess.Popen(cmd,cwd=str(EX),env=env,stdin=subprocess.DEVNULL,stdout=so,stderr=se,start_new_session=True)
a=Path('/proc/%d/stat'%child.pid).read_text().rsplit(')',1)[1].split()
receipt={'status':'STARTED','experiment_id':EXPERIMENT,'pid':child.pid,'pgrp':int(a[2]),'start_ticks':int(a[19]),'started_unix':time.time(),'registered_datasets':10,'registered_fits':120,'starts_per_dataset':12,'arms':['original'],'fitted_factor_capacity':[4,4,4],'fitted_fpca_capacity':4,'initial_concurrency':2,'max_concurrency':8,'automatic_retries':0,'automatic_completion_stopcharging':True,'local_return_required_before_stop':False,'process_baseline_sha256':baseline_sha,'resource_acceptance_sha256':sha(STAGE/'RESOURCE_PROBE_ACCEPTANCE_V3.json'),'external_stopped_confirmation_pending':True}
assert receipt['pid']==receipt['pgrp']
life.persist(target,receipt)
print(json.dumps(receipt,indent=2))
