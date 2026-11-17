import json,os,subprocess,sys,time,hashlib
from pathlib import Path
EX=Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT=Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
S=Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
A=S/'amendments/probe_call_v2'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
def save(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
assert sys.platform.startswith('linux') and not (OUT/'fits').exists()
old=read(S/'RUNTIME_ACCEPTANCE.json')
assert old['status']=='FAIL' and old['formal_fit_attempts']==0
failed=OUT/'resource_probe/original_L4M4_2sweep'
assert not (failed/'PROBE_FIT.rds').exists() and not (failed/'INITIAL_FIELD_OBJECT_SIZES.csv').exists()
assert 'With annealing, maxit' in read(failed/'RESOURCE_PROBE_COMPLETE.json')['error']
assert read(S/'PREPARE_ACCEPTANCE.json')['status']=='PASS'
assert read(OUT/'startup/RUNNER_ACCEPTANCE.json')['status']=='PASS'
assert read(S/'STOPCHARGING_DRYRUN.json')['status']=='PASS'
assert sha(S/'RUNTIME_LIBRARY_BINDING.json')==sha(EX/'RUNTIME_LIBRARY_BINDING.json')
dep=read(EX/'RUNTIME_LIBRARY_BINDING.json')
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',FSYNC_EXEC_ROOT=str(EX),FSYNC_OUTPUT_ROOT=str(OUT),R_LIBS=dep['R_LIBS'],R_LIBS_USER=dep['R_LIBS_USER'])
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
logs=S/'repair_probe_v2_logs';logs.mkdir()
save(logs/'START.json',{'pid':os.getpid(),'started_unix':time.time(),'previous_probe_updates':0,'new_scientific_starts':0,'new_data_generation':0})
results=[]
def run(name,cmd):
    t=time.monotonic()
    save(logs/(name+'.START.json'),{'started_unix':time.time(),'command':cmd})
    with (logs/(name+'.log')).open('x') as f:
        child=subprocess.Popen(cmd,cwd=str(EX),env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
        code=child.wait()
    row={'stage':name,'exit_code':code,'wall_seconds':time.monotonic()-t}
    save(logs/(name+'.EXIT.json'),row);results.append(row)
    assert code==0,'Correction preparation stage failed: '+name
try:
    run('corrected_resource_probe',['Rscript','--vanilla',str(A/'resource_probe.R')])
    run('science_probe',['Rscript','--vanilla',str(EX/'science/evaluate_timing_endpoint.R'),'--exec-root',str(EX),'--output-root',str(OUT),'--data-id','NC01_TRUE222_20261117__rep01','--fit-path',str(OUT/'resource_probe/original_L4M4_2sweep_v2/PROBE_FIT.rds')])
    assert read(OUT/'resource_probe/SCIENCE_SMOKE/PASS.json')['status']=='PASS'
    run('resource_assessment',['python3','-B',str(A/'assess_resource_probe.py')])
    assert read(S/'RESOURCE_PROBE_ACCEPTANCE.json')['status']=='PASS'
    runtime={'status':'PASS','stopcharging_permission_dryrun':True,'startup_original_gate':True,'science_tests':True,'registration_checked':True,'runtime_library_binding_sha256':sha(S/'RUNTIME_LIBRARY_BINDING.json'),'actual_science_smoke':True,'resource_probe_acceptance':True,'formal_fit_attempts':0,'actual_stop_requests':0,'original_failure_retained':True,'correction':'Explicit NULL list argument for independent two-update probe; original estimator and formal fit settings unchanged','original_preparation_sha256':sha(S/'RUNTIME_ACCEPTANCE.json'),'amendment_manifest_sha256':sha(A/'AMENDMENT.json'),'stages':old['stages']+results}
    save(S/'RUNTIME_ACCEPTANCE_V2.json',runtime)
    print(json.dumps({'status':'PASS','formal_fits':0,'extra_scientific_starts':0}),flush=True)
except BaseException as e:
    failure={'status':'FAIL','error_type':type(e).__name__,'error':str(e),'stages':results,'automatic_retries':0}
    save(logs/'FAILURE.json',failure);print(json.dumps(failure),flush=True);sys.exit(2)
